## This script generates simulation scenarios for the ordinal data simulations.

## These have to be slightly different from the scenarios for the MVN data for two reasons:
## 1) The residual SD from the GEE models is larger with the ordinal datasets than
##.   the normal datasets, which makes the Cohen's D lower. To get a Cohen's D of
##.   around 0.5, we have to use a larger raw treatment effect in the ordinal data.
## 2) There is a ceiling effect in the ordinal data, where the linear treatment and
##.   baseline response effects can't be so large that the expected response would
##.   be greater than the max value of the support. So we have to carefully set
##.   the effects so that the most extreme effects will be unlikely to cross this threshold.

## Varying simulation settings:
## 1) True treatment effects
##    a) Null effect 
##    b) Varying effects centered at 0
##    c) Varying with positive mean
## 2) True baseline effects
##    a) Varying with positive mean
## 3) Sample size
##    a) 40
##    b) 200

library(here)
library(readr)
source(here("Scripts", "analysis_functions.R"))

n <- c(40, 200)
n_items <- 12

target_baseline_effect <- 0.4 # baseline response effect on each item

## Scaling factors for baseline response associations. We will use these to scale quantiles
## of the standard normal distribution to create the baseline response coefficients.
baseline_gammas <- c(0, 0.2, 0.4, 0.6)

## Treatment effect coefficients
target_total_effect <- 6.85 # target treatment effect on total score
target_treatment_effect <- target_total_effect/n_items # target treatment effect per item
treatment_effect_means <- c(0, target_treatment_effect)
treatment_gammas <- c(0, 0.15, 0.3) ## Scaling factors similar to the baseline response scaling factors.

within_corr <- 0.6 # item correlation within subscale
without_corr <- 0.4 # item correlation from different subscales
subscale_sizes <- c(4, 4, 4)
subscale_list <- vector(mode="list", length=length(subscale_sizes))
for (i in 1:length(subscale_list)) {
  if (i==1) {
    subscale_list[[i]] <- c(1:subscale_sizes[i])
  } else {
    subscale_list[[i]] <- c((cumsum(subscale_sizes)[i-1]+1):cumsum(subscale_sizes)[i])
  }
}

## Coefficients for cumulative logit model
alpha0 <- c(-1.1, -0.1, 0.5, 1.2, 1.9, 2.5)
## Equivalent intercept coefficients for linear model
beta0 <- expected_response(expit(alpha0))$expected_response


baseline_effects <- lapply(1:length(baseline_gammas), function(i) {
  baseline_quantiles <- baseline_gammas[i]*(qnorm((1:n_items)/(n_items+1), target_baseline_effect, 1))
  return(baseline_quantiles - mean(baseline_quantiles) + target_baseline_effect)
})

treatment_effects <- lapply(1:length(treatment_effect_means), function(i) 
  lapply(1:length(treatment_gammas), function(j) {
    treatment_quantiles <- treatment_gammas[j]*(qnorm((1:n_items)/(n_items+1), treatment_effect_means[i], 1))
    return(treatment_quantiles - mean(treatment_quantiles) + treatment_effect_means[i])
    
  }))

simulation_scenarios <- vector(mode="list", 
                               length=length(baseline_gammas)*length(treatment_effect_means)*length(treatment_gammas)*length(n))
scenario <- 1
for (i in 1:length(n)) {
  
  for (j in 1:length(treatment_effect_means)) {
    
    for (k in 1:length(treatment_gammas)) {
      
      for (m in 1:length(baseline_gammas)) {
        simulation_scenarios[[scenario]]$n <- n[i]
        simulation_scenarios[[scenario]]$treatment_effect_mean <- treatment_effect_means[j]
        simulation_scenarios[[scenario]]$treatment_gamma <- treatment_gammas[k]
        simulation_scenarios[[scenario]]$baseline_gamma <- baseline_gammas[m]
        simulation_scenarios[[scenario]]$baseline_effect_mean <- target_baseline_effect
        simulation_scenarios[[scenario]]$treatment_effects <- treatment_effects[[j]][[k]]
        simulation_scenarios[[scenario]]$baseline_effects <- baseline_effects[[m]]
        simulation_scenarios[[scenario]]$support_list <- lapply(1:n_items, function(i) return(c(0:length(alpha0))))
        simulation_scenarios[[scenario]]$alpha0 <- alpha0
        simulation_scenarios[[scenario]]$beta0 <- beta0
        simulation_scenarios[[scenario]]$within_corr <- within_corr
        simulation_scenarios[[scenario]]$without_corr <- without_corr
        simulation_scenarios[[scenario]]$subscale_sizes <- subscale_sizes
        simulation_scenarios[[scenario]]$subscale_list <- subscale_list
        simulation_scenarios[[scenario]]$n_items <- n_items
        scenario <- scenario + 1
      }
    }
  }
}

write_rds(simulation_scenarios, file=here("Output", "Simulation Scenarios", "simulation_scenarios_ordinal.rds"))
