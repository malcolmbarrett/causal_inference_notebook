# Author contact: Justin Bohn (jmb619@mail.harvard.edu)
# Last revision : 2015-08-04

#######################################################################################
#### Program 16.1                                                                  ####
#### Estimating the average causal using the standard IV estimator                 ####
#### via the calculation of sample averages                                        ####
#### Data from NHEFS                                                               ####
#######################################################################################

# load external packages (install them if necessary)
library(dplyr)      # easy, readable data manipulation
library(multcomp)   # estimation of linear contrasts from GLMs
library(ggplot2)    # graphics
library(gtools)     # easy decile creation
library(sem)        # two-stage least squares
library(geepack)    # generalized estimating equations

# clean up workspace
rm(list = ls())

# establish working directory (this will be different on your system)
setwd('~/Downloads/R_part2/')

# read in NHEFS data and preprocess
# note: this is the dplyr method, but can also be easily accomplished
#       in base R with within(). However, dplyr is faster and more readable.
nhefs <- read.csv('nhefs.csv') %>%
  mutate(
    # define censoring variable cens
    cens = ifelse(is.na(wt82_71), yes = 1, no = 0),
    # categorize the school variable
    education = cut(school, breaks = c(0, 8, 11, 12, 15, 20),
                    include.lowest = TRUE, 
                    labels = c('8th Grage or Less',
                               'HS Dropout',
                               'HS',
                               'College Dropout',
                               'College or More')),
    # establish active as a factor variable
    active = factor(active),
    # establish exercise as a factor variable
    exercise = factor(exercise),
    # create a treatment label variable
    qsmklabel = ifelse(qsmk == 1,
                       yes = 'Quit Smoking 1971-1982',
                       no = 'Did Not Quit Smoking 1971-1982'),
    # define the (proposed) instrument
    highprice = ifelse(price82 >= 1.5, yes = 1, no = 0),
    # counterfactual weight loss under no treatment, assuming an effect of 2.396
    Hpsi = wt82_71 - (qsmk * 2.396)
  ) %>%
  # provisionally ignore those with missing values on some variables
  filter(!is.na(education),
         !is.na(wt82),
         !is.na(price82))

# check categorizations to be sure
xtabs(~ school + education, data = nhefs)

# standard "Wald" estimator
# i.e.  (E[Y|Z=1] - E[Y|Z=0]) / (Pr[A=1|Z=1] - Pr[A=1|Z=0]) 

# first try it "manually"
a <- mean(nhefs$wt82_71[nhefs$highprice == 1], na.rm = TRUE) # E[Y|Z=1]
  a # print
b <- mean(nhefs$wt82_71[nhefs$highprice == 0], na.rm = TRUE) # E[Y|Z=0]
  b # print
c <- mean(nhefs$qsmk[nhefs$highprice == 1], na.rm = TRUE)    # Pr(A=1|Z=1)
  c # print
d <- mean(nhefs$qsmk[nhefs$highprice == 0], na.rm = TRUE)    # Pr(A=1|Z=0)
  d # print
(a-b)/(c-d)                                                  # the IV estimate

# then wrap it in a function
wald <- function(outcome, treatment, iv, data) {
  # the numerator, i.e., E[Y|Z=1] - E[Y|Z=0]
  num <- coef(lm(data[,outcome] ~ data[,iv]))[2]
  # the denominator, i.e., Pr[A=1|Z=1] - Pr[A=1|Z=0]
  denom <- coef(lm(data[,treatment] ~ data[,iv]))[2]
  return(as.numeric(num/denom))
}

wald('wt82_71', 'qsmk', 'highprice', nhefs)

#######################################################################################
#### Program 16.2                                                                  ####
#### Estimating the average causal using the standard IV estimator                 ####
#### via two-stage least squares regression                                        ####
#### Data from NHEFS                                                               ####
#######################################################################################

library(sem) 

model1 <- tsls(wt82_71 ~ qsmk, ~ highprice, data = nhefs)
summary(model1)
confint(model1)  # note the wide confidence intervals

#######################################################################################
#### Program 16.3                                                                  ####
#### Estimating the average causal using the standard IV estimator                 ####
#### via additive marginal structural models                                       ####
#### Data from NHEFS                                                               ####
#######################################################################################

# note: refer to programs from chapter 14 for a consideration of multiple 
#       values for psi
summary(geeglm(highprice ~ Hpsi, 
               data = subset(nhefs, !is.na(highprice) & !is.na(Hpsi)), 
               id = seqn, 
               corstr = 'independence',
               family = 'binomial'))


#######################################################################################
#### Program 16.4                                                                  ####
#### Estimating the average causal using the standard IV estimator                 ####
#### with altnerative proposed instruments                                         ####
#### Data from NHEFS                                                               ####
#######################################################################################

# note what happens to estimates....
summary(tsls(wt82_71 ~ qsmk, ~ ifelse(price82 >= 1.6, 1, 0), data = nhefs))
summary(tsls(wt82_71 ~ qsmk, ~ ifelse(price82 >= 1.7, 1, 0), data = nhefs))
summary(tsls(wt82_71 ~ qsmk, ~ ifelse(price82 >= 1.8, 1, 0), data = nhefs))
summary(tsls(wt82_71 ~ qsmk, ~ ifelse(price82 >= 1.9, 1, 0), data = nhefs))

#######################################################################################
#### Program 16.5                                                                  ####
#### Estimating the average causal using the standard IV estimator                 ####
#### Conditional on baseline covariates                                            ####
#### Data from NHEFS                                                               ####
#######################################################################################

summary(tsls(wt82_71 ~ qsmk + sex + race + age + smokeintensity + smokeyrs + 
                      exercise + active + wt71,
             ~ highprice + sex + race + age + smokeintensity + smokeyrs + exercise +
               active + wt71, data = nhefs))
