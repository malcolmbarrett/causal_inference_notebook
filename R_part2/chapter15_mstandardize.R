# Author contact: Justin Bohn (jmb619@mail.harvard.edu)
# Last revision : 2015-08-04
#
# note: This is a function provided to perform parametric standardization
#       to estimate the effect of a binary, non-time-varying treatment on
#       a binary or continuous outcome. It allows for estimation of sub-
#       population effects (i.e., effect in the treated and the untreated)
#       and obtains confidence intervals by bootstrapping from the study
#       population. For binary outcomes, the family should be specified as 
#       "binomial" which will utilize a logistic model. For continuous outcome,
#       the family should be specified as "gaussian" (the default). Other types of 
#       GLM may work, but have not been tested. Please report bugs. 

mstandardize <- function(formula,                   # the model formula 
                         family,                    # the model family
                         data,                      # dataset containing cohort
                         trt,                       # name of treatment var (coded 0,1)
                         pop = 0:1,                 # target population (see below)
                                                    #   0 = untreated 
                                                    #   1 = treated (e.g. SMR)
                                                    #   0:1 = whole population
                         nboot = 100,               # number of bootstrap samples for CIs
                         cl = 0.95,                 # confidence level for CIs
                         dig = getOption('digits'), # number of digits for printing
                         show.model = FALSE,        # logical, should the model results be printed
                         ...                        # additional arguments passed to glm()
                         ) {
  #### input checks ####
  if(!family %in% c('binomial', 'gaussian')) {
    stop('Unrecognized family supplied. Options are: binomial, gaussian')
  }
  if(any(!data[,trt] %in% 0:1)) {
    stop('Treatment indicator (', trt, ') not coded 0/1')
  }
  
  #### fit the model ###
  mod <- glm(formula, family, data, ...)
  
  #### all-treated dataset ####
  pop.idx <- data[,trt] %in% pop

  treated <- data[pop.idx, ]  # copy (and subset) original dataset
  treated[,trt] <- 1          # set all to treated
  
  #### all-untreated dataset ####
  untreated <- data[pop.idx, ] # copy (and subset) original dataset
  untreated[,trt] <- 0         # set it all to untreated

  #### obtain point estimates ####
  mean.treated <- mean(predict.glm(mod, treated, 'response'))
  mean.untreated <- mean(predict.glm(mod, untreated, 'response'))
  difference <- mean.treated - mean.untreated
  ratio <- mean.treated/mean.untreated
  
  #### bootstrapping ####
  boots <- data.frame(i = 1:nboot,
                      mean1 = NA,
                      mean0 = NA,
                      difference = NA,
                      ratio = NA)
  N <- nrow(data)
  for(i in 1:nboot) {
    # sample with replacement
    sampl <- data[sample(1:N, N, replace = TRUE), ]

    # perform standardization
    #### fit the model ###
    bmod <- glm(formula, family, data = sampl, ...)
    
    #### all-treated dataset ####
    bpop.idx <- sampl[,trt] %in% pop

    btreated <- sampl[bpop.idx, ]  # copy (and subset) original dataset
    btreated[,trt] <- 1            # set all to treated
    
    #### all-untreated dataset ####
    buntreated <- sampl[bpop.idx, ] # copy (and subset) original dataset
    buntreated[,trt] <- 0           # set all to untreated
    
    #### obtain point estimates ####
    bmean.treated <- mean(predict.glm(bmod, btreated, 'response'))
    bmean.untreated <- mean(predict.glm(bmod, buntreated, 'response'))
    bdifference <- bmean.treated - bmean.untreated
    bratio <- bmean.treated/bmean.untreated
    
    
    boots[i,'mean1'] <- bmean.treated
    boots[i,'mean0'] <- bmean.untreated
    boots[i,'difference'] <- bmean.treated - bmean.untreated
    boots[i,'ratio'] <- bmean.treated/bmean.untreated
  }
  
  #### return results ####
  Z <- qnorm(p = cl + ((1-cl)/2))
  mean1.lo <- mean.treated - (Z * sd(boots$mean1))
  mean1.hi <- mean.treated + (Z * sd(boots$mean1))
  mean0.lo <- mean.untreated - (Z * sd(boots$mean0))
  mean0.hi <- mean.untreated + (Z * sd(boots$mean0))
  difference.lo <- difference - (Z * sd(boots$difference))
  difference.hi <- difference + (Z * sd(boots$difference))
  ratio.lo <- ratio - (Z * sd(boots$ratio))
  ratio.hi <- ratio + (Z * sd(boots$ratio))
  
  out <- data.frame(
    Estimate = c(mean.treated, mean.untreated, difference, ratio),
    LCL = c(mean1.lo, mean0.lo, difference.lo, ratio.lo),
    UCL = c(mean1.hi, mean0.hi, difference.hi, ratio.hi)
  )
  
  if(family == 'binomial' & identical(pop[order(pop)],0:1)){
    rownames(out) <- c('Standardized Risk in the Treated',
                       'Standardized Risk in the Untreated',
                       'Standardized Risk Difference in Pop.',
                       'Standardized Risk Ratio in Pop.')
  } 
  if(family == 'binomial' & all(pop == 1)){
    rownames(out) <- c('Standardized Risk in the Treated',
                       'Standardized Risk in the Untreated',
                       'Standardized Risk Difference in the Treated',
                       'Standardized Risk Ratio in the Treated')
  } 
  if(family == 'binomial' & all(pop == 0)){
    rownames(out) <- c('Standardized Risk in the Treated',
                       'Standardized Risk in the Untreated',
                       'Standardized Risk Difference in the Untreated',
                       'Standardized Risk Ratio in the Untreated')
  } 
  if(family == 'gaussian' & identical(pop[order(pop)], 0:1)){
    rownames(out) <- c('Standardized Mean in the Treated',
                       'Standardized Mean in the Untreated',
                       'Standardized Mean Difference in Pop',
                       'Standardized Mean Ratio in Pop')
  } 
  if(family == 'gaussian' & all(pop == 1)){
    rownames(out) <- c('Standardized Mean in the Treated',
                       'Standardized Mean in the Untreated',
                       'Standardized Mean Difference in the Treated',
                       'Standardized Mean Ratio in the Treated')
  } 
  if(family == 'gaussian' & all(pop == 0)){
    rownames(out) <- c('Standardized Mean in the Treated',
                       'Standardized Mean in the Untreated',
                       'Standardized Mean Difference in the Untreated',
                       'Standardized Mean Ratio in the Untreated')
  }
  cat('\n')
  cat('Results of Model-Based Standardization\n')
  cat(paste0(rep('-', getOption('width')), collapse = ''), '\n')
  if(show.model) {
    print(summary(mod))
  }
  cat('  Outcome:', as.character(formula)[2], '\n')
  cat('  Model Type:', family(mod)$family, '\n')
  cat('  Model RHS:', as.character(formula)[3], '\n')
  cat('  No. of bootstrap resamples:', nboot, '\n')
  cat('  Showing', cl*100, '% confidence intervals\n')
  return(round(out, digits = dig))
} 