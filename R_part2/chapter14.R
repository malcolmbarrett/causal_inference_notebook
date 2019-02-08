# Feel free to report errors/suggestions to ehsan@stat.ubc.ca

################################################################
# PROGRAM 14.1
# Preprocessing, ranks of extreme observations, IP weights for censoring
# Data from NHEFS
################################################################

# rm(list=ls(all=TRUE))
setwd("C:\\Ehsan\\causal\\nhefs_book") # contains SAS/CSV datafile
# reading external data file
#install.packages("sas7bdat")
library(sas7bdat)
nhefs <- read.sas7bdat("nhefs_book.sas7bdat")
# or if CSV file is available:
#nhefs <- read.csv("nhefs.csv") 
nhefs <- as.data.frame(nhefs)
dim(nhefs)

# some preprocessing of the data
nhefs$id <- 1:nrow(nhefs)
nhefs$cens <- as.numeric(is.na(nhefs$wt82))
nhefs$older <- as.numeric(nhefs$age > 50 & !is.na(nhefs$age))
# install.packages("car")
library(car)
nhefs$school[is.na(nhefs$school)] <- -10
sum(nhefs$school == -10)
nhefs$education.code <- recode(nhefs$school, " 0:8 = 1; 9:11 = 2; 12 = 3; 
                         13:15 = 4; 16:hi = 5; -10 = 6 ")
nhefs$education <- recode(nhefs$education.code, " 1 = '1. 8th grade or less';
                    2 = '2. HS dropout';
                    3 = '3. HS';
                    4 = '4. College dropout';
                    5 = '5. College or more';
                    6 = NA ")

# with non-missing values in the following covariates
nhefs.original <- nhefs # Original data saved for later use
dim(nhefs)
nhefs.all <- nhefs.original[c("id", "qsmk", "sex", "race", "age", "school", 
                  "smokeintensity", "smokeyrs", "exercise", "active", 
                  "wt71", "wt82", "education", "cens", "wt82_71")]
nhefs.nonmissing <- nhefs.all[!is.na(nhefs.all$education),]
dim(nhefs.all)
dim(nhefs.nonmissing)
# [1] 1629   15
sum(is.na(nhefs.all))

# restricting data for uncensored, for comparison with observed outcome
nhefs0 <- subset(nhefs.nonmissing, cens == 0)
dim(nhefs0)
sum(is.na(nhefs0))

# ranking of extreme observations
#install.packages("Hmisc")
library(Hmisc)
describe(nhefs0$wt82_71)

# estimation of denominator of ip weights for C
cens.prob.est <- glm(as.factor(cens)~ as.factor(qsmk) + as.factor(sex) + as.factor(race) + 
                 age + I(age^2) + as.factor(education) + smokeintensity + I(smokeintensity^2) + 
                 smokeyrs + I(smokeyrs^2) + as.factor(exercise) + as.factor(active) + wt71 + 
                 I(wt71^2), data = nhefs.nonmissing, family = binomial("logit"))
summary(cens.prob.est)
nhefs0$w.cens <- 1/(1-predict(cens.prob.est, nhefs0, type = "response"))
describe(nhefs0$w.cens)
summary(nhefs0$w.cens)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.00    1.01    1.02    1.04    1.04    1.82 

##################################################################
#  PROGRAM 14.2
# G-estimation of a 1-parameter structural nested mean model
# Brute force search
# Data from NHEFS
##################################################################

nhefs.g.est <- nhefs0
nhefs.g.est$psi <- 3.446
nhefs.g.est$Hpsi <- nhefs.g.est$wt82_71 - nhefs.g.est$psi * nhefs.g.est$qsmk 

#install.packages("geepack")
require(geepack)
gee.obj <- geeglm(qsmk ~  as.factor(sex) + as.factor(race) + 
  age + I(age^2) + as.factor(education) + smokeintensity + I(smokeintensity^2) + 
  smokeyrs + I(smokeyrs^2) + as.factor(exercise) + as.factor(active) + wt71 + 
  I(wt71^2)+Hpsi, data = nhefs.g.est, weight = w.cens, 
  id=id, corstr="independence", family = binomial(logit))
summary(gee.obj)
coef(gee.obj)["Hpsi"]
#   Hpsi 
#-1.9e-06 

##################################################################
# G-estimation: Checking multiple possible values of psi*/
##################################################################
require(geepack)
data <- nhefs.g.est
grid <- seq(from = 2,to = 5, by = 0.1) # set by = 0.001 for finer estimate
j = 0
store.Hpsi.coefs <- double(length(grid))
for (i in grid){
  psi = i
  j = j+1
  data$Hpsi <- data$wt82_71 - psi * data$qsmk 
  
  gee.obj <- geeglm(qsmk ~  as.factor(sex) + as.factor(race) + 
                      age + I(age^2) + as.factor(education) + smokeintensity + I(smokeintensity^2) + 
                      smokeyrs + I(smokeyrs^2) + as.factor(exercise) + as.factor(active) + wt71 + 
                      I(wt71^2)+Hpsi, data = data, weight = w.cens, 
                    id=id, corstr="independence", family = binomial(logit))
  store.Hpsi.coefs[j] <- coef(gee.obj)["Hpsi"]
  cat("Iteration", j, "completed\n")
}
store.results <- as.data.frame(cbind(grid, abs(store.Hpsi.coefs)))
names(store.results)
names(store.results) <- c("grid", "Hpsi.est")  
store.results[store.results$Hpsi.est == min(store.results$Hpsi.est),]
# for 0.001 interval:
#      grid  Hpsi.est
#1447 3.446 1.903e-06

##################################################################
#  PROGRAM 14.3
# G-estimation for 2-parameter structural nested mean model
# Closed form estimator
# Data from NHEFS
##################################################################

##################################################################
# G-estimation: Closed form estimator linear mean models  
##################################################################
  

logit.est <- glm(as.factor(qsmk) ~ as.factor(sex) + as.factor(race) + 
                       age + I(age^2) + as.factor(education) + smokeintensity + I(smokeintensity^2) + 
                       smokeyrs + I(smokeyrs^2) + as.factor(exercise) + as.factor(active) + wt71 + 
                       I(wt71^2), data = nhefs0, weight = w.cens, family = binomial("logit"))
summary(logit.est)
nhefs0$qsmk.pred <- predict(logit.est, nhefs0, type = "response")
describe(nhefs0$qsmk.pred)
summary(nhefs0$qsmk.pred)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0514  0.1770  0.2410  0.2610  0.3230  0.7890 


# solve sum(w_c * H(psi) * (qsmk - E[qsmk | L]))  = 0
# for a single psi and H(psi) = wt82_71 - psi * qsmk
# this can be solved as psi = sum( w_c * wt82_71 * (qsmk - pqsmk)) / sum(w_c * qsmk * (qsmk - pqsmk))
  
with(nhefs0, sum( w.cens * wt82_71 * (qsmk - qsmk.pred)) / sum(w.cens * qsmk * (qsmk - qsmk.pred)))
# [1] 3.446

##################################################################
# G-estimation: Closed form estimator for 2-parameter model
##################################################################

diff = with(nhefs0, qsmk - qsmk.pred)
diff2 = with(nhefs0, w.cens * diff)

lhs = matrix(0, 2,2)
lhs[1,1] = with(nhefs0, sum( qsmk * diff2))
lhs[1,2] = with(nhefs0, sum( qsmk * smokeintensity  * diff2 ))
lhs[2,1] = with(nhefs0, sum( qsmk * smokeintensity * diff2))
lhs[2,2] = with(nhefs0, sum( qsmk * smokeintensity * smokeintensity * diff2 ))
#> lhs
#       [,1]   [,2]
#[1,]  292.1   5702
#[2,] 5701.5 153045
                                                                
rhs = matrix(0, 2,1)
rhs[1] = with(nhefs0, sum(wt82_71 * diff2 ))
rhs[2] = with(nhefs0, sum(wt82_71 * smokeintensity * diff2 ))
#> rhs
#      [,1]
#[1,]  1006
#[2,] 20901

psi = t(solve(lhs,rhs))
#> psi
#      [,1]    [,2]
#[1,] 2.859 0.03004
