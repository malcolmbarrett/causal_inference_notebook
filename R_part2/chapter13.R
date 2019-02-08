# Feel free to report errors/suggestions to ehsan@stat.ubc.ca

################################################################
# PROGRAM 13.1
# Estimating the mean outcome within levels of treatment 
# and confounders: Data from NHEFS
################################################################

# rm(list=ls(all=TRUE))
# setwd("C:\\Users\\Ehsan\\Dropbox\\causal") # contains SAS/CSV datafile
# reading external data file
# install.packages("sas7bdat")
library(sas7bdat)
nhefs <- read.sas7bdat("nhefs_book.sas7bdat")
# or if CSV file is available:
# nhefs <- read.csv("nhefs.csv") 
nhefs <- as.data.frame(nhefs)
dim(nhefs)

# some preprocessing of the data

nhefs$cens <- as.numeric(is.na(nhefs$wt82))
nhefs$older <- as.numeric(nhefs$age > 50 & !is.na(nhefs$age))
# install.packages("car")
library(car)
nhefs$education.code <- recode(nhefs$school, " 0:8 = 1; 9:11 = 2; 12 = 3; 
                         13:15 = 4; 16:hi = 5; NA = 6 ")
nhefs$education <- recode(nhefs$education.code, " 1 = '1. 8th grade or less';
                    2 = '2. HS dropout';
                    3 = '3. HS';
                    4 = '4. College dropout';
                    5 = '5. College or more';
                    6 = 'Unknown' ")

# Analysis restricted to N=1566 
# with non-missing values in the following covariates
nhefs.original <- nhefs # Original data saved for later use
nhefs$id <- 1:nrow(nhefs)
nhefs2 <- nhefs[c("id", "qsmk", "sex", "race", "age", "school", 
                  "smokeintensity", "smokeyrs", "exercise", "active", 
                  "wt71", "wt82")]
dim(nhefs2)

# restricting data for non-missing
nhefs2 <- as.data.frame(na.omit(nhefs2))
dim(nhefs2)
nhefs.id.matched <- subset(nhefs, id %in% nhefs2$id)
dim(nhefs)

# restricting data for uncensored, for comparison with observed outcome
nhefs0 <- subset(nhefs.id.matched, cens == 0)
dim(nhefs0)

table(nhefs0$qsmk) # untreated vs treated

# Estimates
glm.obj <- glm(wt82_71~ as.factor(qsmk) + as.factor(sex) + as.factor(race) + 
  age + I(age^2) + as.factor(education) + smokeintensity + I(smokeintensity^2) + 
  smokeyrs + I(smokeyrs^2) + as.factor(exercise) + as.factor(active) + wt71 + 
  I(wt71^2) + I(qsmk*smokeintensity), data = nhefs0)
summary(glm.obj)

nhefs0$meanY <- predict(glm.obj, type = "response")

# print variable values corresponding to subject with unique identifier 24770
nhefs0[nhefs0$seqn ==24770,c("meanY", "qsmk", "sex", "race", "age", "education",
                             "smokeintensity", "smokeyrs", "exercise", "active", 
                             "wt71")]

summary(nhefs0$meanY)
summary(nhefs0$wt82_71)

################################################################
# PROGRAM 13.2
# Standardizing the mean outcome to the baseline confounders
# Data from Table 2.2
################################################################

id <- c("Rheia", "Kronos", "Demeter", "Hades", "Hestia", "Poseidon", 
  "Hera", "Zeus", "Artemis", "Apollo", "Leto", "Ares", "Athena", 
  "Hephaestus", "Aphrodite", "Cyclope", "Persephone", "Hermes", 
  "Hebe", "Dionysus")
N <- length(id)
L <- c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
A <- c(0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1)
Y <- c(0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0)
interv <- rep(-1, N)
observed <- cbind(L, A, Y, interv)
untreated <- cbind(L, rep(0, N), rep(NA, N), rep(0, N))
treated <- cbind(L, rep(1, N), rep(NA, N), rep(1, N))
data <- as.data.frame(rbind(observed, untreated, treated))
data$id <- rep(id, 3)

# Estimates
glm.obj <- glm(Y~ A*L, data = data)
summary(glm.obj)

# data$meanY <- c(predict(glm.obj, data[data$interv == -1,], type = "response"),
##                predict(glm.obj, data[data$interv == 0,], type = "response"), 
#                predict(glm.obj, data[data$interv == 1,], type = "response"))

data$meanY <- c(predict(glm.obj, data, type = "response"))

#index <- sort(unique(data$interv))
#ia <- mean(data$meanY[data$interv == index[1]])
#ib <- mean(data$meanY[data$interv == index[2]])
#ic <- mean(data$meanY[data$interv == index[3]])
#predicted.values <- c(ia, ib, ic)
#cbind(index, predicted.values)

with(data, tapply(meanY, interv, mean))

################################################################
# PROGRAM 13.3
# Standardizing the mean outcome to the baseline confounders:
# Data from NHEFS
################################################################

# 1st copy: equal to original one
names(nhefs)
nhefs0$interv <- rep(-1, nrow(nhefs0))
# nhefs0$qsmk
# nhefs0$wt82_71
dim(nhefs0)
# 2nd copy: treatment set to 0, outcome to missing
nhefs.untr <- nhefs0
nhefs.untr$interv <- rep(0, nrow(nhefs.untr))
nhefs.untr$qsmk <- rep(0, nrow(nhefs.untr))
nhefs.untr$wt82_71 <- rep(NA, nrow(nhefs.untr))
dim(nhefs.untr)
# 3rd copy: treatment set to 1, outcome to missing
nhefs.tr <- nhefs0
nhefs.tr$interv <- rep(1, nrow(nhefs.tr))
nhefs.tr$qsmk <- rep(1, nrow(nhefs.tr))
nhefs.tr$wt82_71 <- rep(NA, nrow(nhefs.tr))
dim(nhefs.tr)
# create a dataset with 3 copies of each subject
onesample <- as.data.frame(rbind(nhefs0, nhefs.untr, nhefs.tr))
#onesample$qsmk
dim(onesample)

# Estimates
# linear model to estimate mean outcome conditional on treatment & confounders,
# parameters are estimated using original observations only (interv= -1),
# parameter estimates are used to predict mean outcome for observations 
# with treatment set to 0 (interv=0) and to 1 (innterv=1);

glm.obj <- glm(wt82_71~ as.factor(qsmk) + as.factor(sex) + as.factor(race) + 
  age + I(age^2) + as.factor(education) + smokeintensity + I(smokeintensity^2) + 
  smokeyrs + I(smokeyrs^2) + as.factor(exercise) + as.factor(active) + wt71 + 
  I(wt71^2), data = onesample)
summary(glm.obj)

#onesample$meanY <- c(predict(glm.obj, nhefs0, type = "response"),
#                     predict(glm.obj, nhefs.untr, type = "response"), 
#                     predict(glm.obj, nhefs.tr, type = "response"))

onesample$meanY <- predict(glm.obj, onesample, type = "response")

# estimate mean outcome in each of the groups interv=0, and interv=1;
# this mean outcome is a weighted average of the mean outcomes 
# in each combination of values of treatment and confounders, 
# that is, the standardized outcome;

#index <- sort(unique(onesample$interv))
#ia <- mean(onesample$meanY[onesample$interv == index[1]])
#ib <- mean(onesample$meanY[onesample$interv == index[2]])
#ic <- mean(onesample$meanY[onesample$interv == index[3]])
#predicted.values <- c(ia, ib, ic)
#cbind(index, predicted.values)

with(onesample, tapply(meanY, list(interv), mean))

################################################################
# PROGRAM 13.4
# Computing the 95% confidence interval of the standardized means 
# and their difference: Data from NHEFS
################################################################

# install.packages("boot")
require(boot)
# Compute basic bootstrap confidence interval.
# Below is the function to estimate mean outcome in each of the groups 
# interv=-1, interv=0 & interv=1 conditional on treatment and confounders,
# also include the mean outcome of interv = 1 vs. interv = 0
boot.fun <- function(dat, index, type){
  sampled.data <- dat[index,]
  fit <- glm(formula = wt82_71~ as.factor(qsmk) + as.factor(sex) + as.factor(race) + 
    age + I(age^2) + as.factor(education) + smokeintensity + I(smokeintensity^2) + 
    smokeyrs + I(smokeyrs^2) + as.factor(exercise) + as.factor(active) + wt71 + 
    I(wt71^2), data = sampled.data)
  sampled.data$meanY <- predict(fit, sampled.data, type = "response")
  std.mean <- with(sampled.data, tapply(meanY, list(interv), mean))
  mean.list <- as.numeric(unlist(c(std.mean, diff(std.mean[c(2,3)]))))
  return(mean.list)
    }

# Parametric g-formula: Bootstrap results using nboot samples
nboot <- 100
# number of bootstrap samples to run
# creation of bootstarp samples using fixed seed to reproduce results
set.seed(1232)
bootres <- boot(onesample, boot.fun, R = nboot)
bootres
nrow(bootres$t)
boot.stat <- matrix(NA, 4, 4)
for (i in 1:4){
boot.stat[i,] <- c(bootres$t0[i], sd(bootres$t[,i]), 
                  bootres$t0[i] - 1.96 * sd(bootres$t[,i]), 
                  bootres$t0[i] + 1.96 * sd(bootres$t[,i]))
}
dimnames(boot.stat) <- list(c("observed", "No treatment", 
                              "Treatment", "Difference"), 
                            c("mean", "sd", "LCL", "UCL"))
boot.stat

