# Author contact: Justin Bohn (jmb619@mail.harvard.edu)
# Last revision : 2015-08-04

#######################################################################################
#### Program 15.1                                                                  ####
#### Estimating the average causal effect within levels of confounders             ####
#### under the assumption of effect-measure modification by smoking intensity ONLY ####
#### Data from NHEFS                                                               ####
#######################################################################################

# load external packages (install if necessary)
library(dplyr)      # easy, readable data manipulation
library(multcomp)   # estimation of linear contrasts from GLMs
library(ggplot2)    # graphics
library(gtools)     # easy decile creation

# clean up workspace
rm(list = ls())

# establish working directory (the path may be different on your system)
setwd('~/Downloads/R_part2/')

# read in NHEFS data and preprocess
nhefs <- read.csv('nhefs.csv') %>%
  mutate(
    # define censoring variable cens
    cens = ifelse(is.na(wt82_71), yes = 1, no = 0),
    # categorize the school variable
    education = cut(school, breaks = c(0, 8, 11, 12, 15, 20),
                    include.lowest = TRUE, 
                    labels = c('1. 8th Grage or Less',
                               '2. HS Dropout',
                               '3. HS',
                               '4. College Dropout',
                               '5. College or More')),
    # establish active as a factor variable
    active = factor(active),
    # establish exercise as a factor variable
    exercise = factor(exercise),
    # create a treatment label variable
    qsmklabel = ifelse(qsmk == 1,
                       yes = 'Quit Smoking 1971-1982',
                       no = 'Did Not Quit Smoking 1971-1982')
  ) %>%
  # provisionally ignore those with missing values on some variables
  filter(!is.na(education))

# check the categorizations of education to be sure
xtabs(~ school + education, data = nhefs)




# model 1: regression on covariates, allowing for some effect modification
# notes: (1) poly(x, 2) adds an orthogonal polynomial of degree 2, add the argument raw = TRUE 
#            if you want it to produce the same coefficients as would x + x^2
#        (2) x1*x2 enters the main effects of x1 and x2 and their product term
#            x1:x2 enters just the product term (necessary here for smokeintensity because
#            we want smokeintensity treated linearly in the interaction but quadratically in
#            the main effect and thus a linear term for smokeintensity is not estimable)
#        (3) observations with missing values are automatically deleted
model1 <- glm(wt82_71 ~ qsmk + sex + race + poly(age, 2, raw = TRUE) + education +
              poly(smokeintensity, 2, raw = TRUE) + poly(smokeyrs, 2, raw = TRUE) +
              + exercise + active + poly(wt71, 2, raw = TRUE) + qsmk:smokeintensity,
              data = nhefs)
summary(model1) 
confint.lm(model1)

# estimate the effect of quitting smoking at various smoking intensities
# see help(glht) for details on estimating linear contrasts

# (step 1) build the contrast matrix with all zeros
# this function builds the blank matrix 
makeContrastMatrix <- function(model, nrow, names) {
  m <- matrix(0, nrow = nrow, ncol = length(coef(model)))
  colnames(m) <- names(coef(model))
  rownames(m) <- names
  return(m)
}
K1 <- makeContrastMatrix(model1, 2, c('Effect of Quitting Smoking at Smokeintensity of 5',
                                      'Effect of Quitting Smoking at Smokeintensity of 40'))
# (step 2) fill in the relevant non-zero elements 
K1[1:2, 'qsmk'] <- 1
K1[1:2, 'qsmk:smokeintensity'] <- c(5, 40)

# (step 3) check the contrast matrix
K1 

# (step 4) estimate the contrasts, get tests and confidence intervals for them
estimates1 <- glht(model1, K1)
  summary(estimates1)
  confint(estimates1)

# model 2, regression on covariates, not allowing for effect modification
model2 <- glm(wt82_71 ~ qsmk + sex + race + poly(age, 2, raw = TRUE) + education +
                poly(smokeintensity, 2, raw = TRUE) + poly(smokeyrs, 2, raw = TRUE) +
                + exercise + active + poly(wt71, 2, raw = TRUE),
              data = nhefs)
summary(model2)
confint.lm(model2)

#######################################################################################
#### Program 15.2                                                                  ####
#### Estimating and plotting the propensity score                                  ####
#### Data from NHEFS                                                               ####
#######################################################################################

# model to estimate propensity score
model3 <- glm(qsmk ~ sex + race + poly(age, 2, raw = TRUE) + education +
              poly(smokeintensity, 2, raw = TRUE) + poly(smokeyrs, 2, raw = TRUE) +
              + exercise + active + poly(wt71, 2, raw = TRUE),
              data = nhefs, family = 'binomial')

# summary(model3) to see the results...

# predict PS values and assign to variable in dataset
nhefs$p.qsmk <- predict(model3, nhefs, 'response')

# view summary statistics for PS by qsmk (base R method)
with(nhefs, by(p.qsmk, qsmk, summary))

# view summary statistics for PS by qsmk (dplyr method)
nhefs %>%
  filter(!is.na(p.qsmk)) %>%
  group_by(qsmk) %>%
  summarize(min = min(p.qsmk),
            p25 = quantile(p.qsmk, 0.25),
            med = median(p.qsmk),
            p75 = quantile(p.qsmk, 0.75),
            max = max(p.qsmk))

# plot PS distribution by qsmk status (much more helpful)
ggplot(nhefs, aes(x = p.qsmk, fill = qsmklabel)) + geom_density(alpha = 0.2) +
  xlab('Probability of Quitting Smoking During Follow-up') +
  ggtitle('Propensity Score Distribution by Treatment Group') +
  scale_fill_discrete('') +
  theme(legend.position = 'bottom', legend.direction = 'vertical')

# alternative plot with histograms  
ggplot(nhefs, aes(x = p.qsmk, fill = qsmklabel, color = qsmklabel)) +
  geom_histogram(alpha = 0.3, position = 'identity') +
  facet_grid(qsmklabel ~ .) +
  xlab('Probability of Quitting Smoking During Follow-up') +
  ggtitle('Propensity Score Distribution by Treatment Group') +
  scale_fill_discrete('') +
  scale_color_discrete('') +
  theme(legend.position = 'bottom', legend.direction = 'vertical')

# attempt to reproduce plot from the book
nhefs %>%
  mutate(ps.grp = round(p.qsmk/0.05) * 0.05) %>%
  group_by(qsmk, qsmklabel, ps.grp) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  mutate(n2 = ifelse(qsmk == 0, yes = n, no =  -1*n)) %>%
  ggplot(aes(x = ps.grp, y = n2, fill = qsmklabel)) +
  geom_bar(stat = 'identity', position = 'identity') +
  geom_text(aes(label = n, x = ps.grp, y = n2 + ifelse(qsmk == 0, 8, -8))) +
  xlab('Probability of Quitting Smoking During Follow-up') +
  ylab('N') +
  ggtitle('Propensity Score Distribution by Treatment Group') +
  scale_fill_discrete('') +
  scale_x_continuous(breaks = seq(0, 1, 0.05)) +
  theme(legend.position = 'bottom', legend.direction = 'vertical',
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

#######################################################################################
#### Program 15.3                                                                  ####
#### Stratification on the propensity score                                        ####
#### Data from NHEFS                                                               ####
#######################################################################################

library(gtools)

# function to create deciles easily
decile <- function(x) {
  return(factor(quantcut(x, seq(0, 1, 0.1), labels = FALSE)))
}

# regression on PS deciles, not allowing for effect modification
model4 <- glm(wt82_71 ~ qsmk + decile(p.qsmk), data = nhefs)
summary(model4)
confint.lm(model4)

# regression on PS deciles, allowing for effect modification
model5 <- glm(wt82_71 ~ qsmk*decile(p.qsmk), data = nhefs)
summary(model5)

# estimate contrasts to see effect in each decile
# (step 1) set up contrast matrix
K5 <- makeContrastMatrix(model5, 10, paste('Effect of Quitting Smoking in PS Decile', 1:10))
# (step 2) fill in contrast matrix
K5[1:10, 'qsmk'] <- 1
K5[2, 'qsmk:decile(p.qsmk)2'] <- 1
K5[3, 'qsmk:decile(p.qsmk)3'] <- 1
K5[4, 'qsmk:decile(p.qsmk)4'] <- 1
K5[5, 'qsmk:decile(p.qsmk)5'] <- 1
K5[6, 'qsmk:decile(p.qsmk)6'] <- 1
K5[7, 'qsmk:decile(p.qsmk)7'] <- 1
K5[8, 'qsmk:decile(p.qsmk)8'] <- 1
K5[9, 'qsmk:decile(p.qsmk)9'] <- 1
K5[10, 'qsmk:decile(p.qsmk)10'] <- 1

# step(3) view the contrast matrix
K5

# estimate contrasts and view tests and CIs
estimates5 <- glht(model5, K5)
summary(estimates5)
confint(estimates5)

#######################################################################################
#### Program 15.4                                                                  ####
#### Standardization using the propensity score                                    ####
#### Data from NHEFS                                                               ####
#######################################################################################

# regression on the propensity score (linear term)
model6 <- glm(wt82_71 ~ qsmk + p.qsmk, data = nhefs)
summary(model6)

# standarization on the propensity score
# (step 1) create two new datasets, one with all treated and one with all untreated
treated <- nhefs
  treated$qsmk <- 1

untreated <- nhefs
  untreated$qsmk <- 0

# (step 2) predict values for everyone in each new dataset based on above model
treated$pred.y <- predict(model6, treated)
untreated$pred.y <- predict(model6, untreated)

# (step 3) compare mean weight loss had all been treated vs. that had all been untreated
mean1 <- mean(treated$pred.y, na.rm = TRUE)
mean0 <- mean(untreated$pred.y, na.rm = TRUE)
mean1
mean0
mean1 - mean0

# (step 4) bootstrap a confidence interval
# number of bootstraps
nboot <- 100
# set up a matrix to store results
boots <- data.frame(i = 1:nboot,
                    mean1 = NA,
                    mean0 = NA,
                    difference = NA)
# loop to perform the bootstrapping
nhefs <- subset(nhefs, !is.na(p.qsmk) & !is.na(wt82_71))
for(i in 1:nboot) {
  # sample with replacement
  sampl <- nhefs[sample(1:nrow(nhefs), nrow(nhefs), replace = TRUE), ]
  
  # fit the model in the bootstrap sample
  bootmod <- glm(wt82_71 ~ qsmk + p.qsmk, data = sampl)
  
  # create new datasets
  sampl.treated <- sampl %>%
    mutate(qsmk = 1)
  
  sampl.untreated <- sampl %>%
    mutate(qsmk = 0)
  
  # predict values
  sampl.treated$pred.y <- predict(bootmod, sampl.treated)
  sampl.untreated$pred.y <- predict(bootmod, sampl.untreated)
  
  # output results 
  boots[i, 'mean1'] <- mean(sampl.treated$pred.y, na.rm = TRUE)
  boots[i, 'mean0'] <- mean(sampl.untreated$pred.y, na.rm = TRUE)
  boots[i, 'difference'] <- boots[i, 'mean1'] - boots[i, 'mean0']
  
  # once loop is done, print the results
  if(i == nboot) {
    cat('95% CI for the causal mean difference\n')
    cat(mean(boots$difference) - 1.96*sd(boots$difference), 
        ',',
        mean(boots$difference) + 1.96*sd(boots$difference))
  }
}

# a more flexible and elegant way to do this is to write a function 
# to perform the model fitting, prediction, bootstrapping, and reporting all at once
# view the code contained in the file mstandardize.R to learn more

# load the code for the mstandardize() function 
# (you may need to change the filepath)
source('chapter15_mstandardize.R') 

# performt the standardization
mstandardize(formula = wt82_71 ~ qsmk + decile(p.qsmk), 
             family = 'gaussian',
             trt = 'qsmk', 
             nboot = 100, 
             data = nhefs)
