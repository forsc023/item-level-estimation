## Code to simulate one questionnaire dataset per simulation scenario for ordinal data simulations
## and estimate the treatment effect using gaussian-family estimation.
## Reads in the scenarios generated from "sim_scenarios_ordinal.R", generates a dataset
## from each scenario, and estimates the treatment effect on each using nine estimators:
## TS; IL-IS with independence, exchangeable, block unstructured, and block independence correlations;
## and IL-ME with independence, exchangeable, block unstructured, and block independence correlations.

## For the simulation study referenced in the main paper, the simulation iterations were
## run in parallel on the Minnesota Supercomputing Institute, with the job array IDs (1-8,000)
## used to specify the seeds.

## Change "arraynum" to run other iterations.
## We saw estimation problems (not converging, GEE getting stuck, etc) more often with
## binomial family estimation, so you might encounter these problems with other iteration numbers
## (arraynum).

arraynum <- 2

library(readr)
library(here)
library(data.table)

source(here("Scripts", "analysis_functions.R"))
#source(here("Scripts", "sim_scenarios_ordinal.R"))

simulation_scenarios <- read_rds(here("Output", "Simulation Scenarios", "simulation_scenarios_ordinal.rds"))

estimates_output <- vector(mode="list", length=length(simulation_scenarios))

## Iterate over all the simulation scenarios
  for (i in 1:length(simulation_scenarios)) {
    print(paste0("scenario = ", i))
    seed <- 100000*i + arraynum
    set.seed(seed)
    target_effects <- simulation_scenarios[[i]]$treatment_effects
    target_baseline_effects <- simulation_scenarios[[i]]$baseline_effects
    n <- simulation_scenarios[[i]]$n
    support_list <- simulation_scenarios[[i]]$support_list
    alpha0 <- simulation_scenarios[[i]]$alpha0
    within_corr <- simulation_scenarios[[i]]$within_corr
    without_corr <- simulation_scenarios[[i]]$without_corr
    subscale_list <- simulation_scenarios[[i]]$subscale_list
    subscale_sizes <- simulation_scenarios[[i]]$subscale_sizes
    n_items <- simulation_scenarios[[i]]$n_items
    
    ## Baseline effect coefficients in cumulative logit model
    alpha2 <- lapply(1:length(target_baseline_effects), function(j) {
      rep(uniroot(function(alpha0, alpha2) {
        expected_response(expit(alpha0 +  alpha2), support_list[[j]])$expected_response -
          expected_response(expit(alpha0), support_list[[j]])$expected_response - target_baseline_effects[j]
      }, c(-20, 20), alpha0=alpha0)$root,
      length(alpha0)
      )
    })

    sim_data <- simulate_dataset(n=n, n_items=n_items, alpha0=alpha0, alpha2=alpha2, 
                                 support_list=support_list, treatment_effects=target_effects,
                                 within_corr=within_corr, without_corr=without_corr, 
                                 subscale_sizes=subscale_sizes)
    sim_data_long <- sim_data$sim_data_long
    K <- max(support_list[[1]])
    sim_data_long_binomial <- sim_data_long %>% mutate(failures = K - response)
    sim_data_binomial <- sim_data$sim_data_wide %>% mutate(failures = n_items*K - score)
    
    ## TS and independence estimates first
    bin_standard_estimate <- tryCatch(get_estimates(data=sim_data_binomial, estimator="standard", 
                                                    working_cor="none", baseline="Y",
                                                    baseline_adjustment="me", covariates="baseline_score",
                                                    covariate_adjustment="me", format="N", 
                                                    dist_family=binomial(link="logit")) %>%
                                        mutate(had_warning=0, had_error=0),
                                      error=function(e) {return(data.frame(
                                        estimator="standard", baseline_name="Main Effects",
                                        estimate=NA, se=NA, lowerci=NA, upperci=NA, p=NA, 
                                        working_cor="N/A", qic=NA, qicu=NA, residual_sd=NA, had_warning=0,
                                        had_error=1
                                      ))},
                                      warning=function(w) {return(data.frame(
                                        estimator="standard", baseline_name="Main Effects",
                                        estimate=NA, se=NA, lowerci=NA, upperci=NA, p=NA, 
                                        working_cor="N/A", qic=NA, qicu=NA,
                                        residual_sd=NA, had_warning=1, had_error=0
                                      ))})
    
    bin_me_estimate_ind <- tryCatch(get_estimates(data=sim_data_long_binomial, estimator="item_maineffects", 
                                                  working_cor="independence", baseline="Y", 
                                                  baseline_adjustment="is",
                                                  covariates="baseline_response", covariate_adjustment="is", 
                                                  format="N", dist_family=binomial(link="logit")) %>%
                                      mutate(had_warning=0, had_error=0),
                                    error=function(e) {return(data.frame(
                                      estimator="item_maineffects", baseline_name="Item-Specific",
                                      estimate=NA, se=NA, lowerci=NA, upperci=NA, p=NA, 
                                      working_cor="Independence", qic=NA, qicu=NA,
                                      residual_sd=NA, had_warning=0, had_error=1
                                    ))},
                                    warning=function(w) {return(data.frame(
                                      estimator="item_maineffects", baseline_name="Item-Specific",
                                      estimate=NA, se=NA, lowerci=NA, upperci=NA, p=NA, 
                                      working_cor="Independence", qic=NA, qicu=NA,
                                      residual_sd=NA, had_warning=1, had_error=0
                                    ))})
    
    bin_item_estimate_ind <- tryCatch(get_estimates(data=sim_data_long_binomial, estimator="item_itemspecific", 
                                                    working_cor="independence", baseline="Y", 
                                                    baseline_adjustment="is",
                                                    covariates="baseline_response", covariate_adjustment="is", 
                                                    format="N", dist_family=binomial(link="logit")) %>%
                                        mutate(had_warning=0, had_error=0),
                                      error=function(e) {return(data.frame(
                                        estimator="item_itemspecific", baseline_name="Item-Specific",
                                        estimate=NA, se=NA, lowerci=NA, upperci=NA, p=NA, 
                                        working_cor="Independence", qic=NA, qicu=NA,
                                        residual_sd=NA, had_warning=0, had_error=1
                                      ))},
                                      warning=function(w) {return(data.frame(
                                        estimator="item_itemspecific", baseline_name="Item-Specific",
                                        estimate=NA, se=NA, lowerci=NA, upperci=NA, p=NA, 
                                        working_cor="Independence", qic=NA, qicu=NA,
                                        residual_sd=NA, had_warning=1, had_error=0
                                      ))})
    
    bin_item_estimate_exch <- tryCatch(get_estimates(data=sim_data_long_binomial, estimator="item_itemspecific", 
                                                       working_cor="exchangeable", baseline="Y", 
                                                       baseline_adjustment="is",
                                                       covariates="baseline_response", covariate_adjustment="is", 
                                                       format="N", dist_family=binomial(link="logit")) %>%
                                           mutate(had_warning=0, had_error=0),
                                         error=function(e) {return(data.frame(
                                           estimator="item_itemspecific", baseline_name="Item-Specific",
                                           estimate=NA, se=NA, lowerci=NA, upperci=NA, p=NA, 
                                           working_cor="Exchangeable", qic=NA, qicu=NA,
                                           residual_sd=NA, had_warning=0, had_error=1
                                         ))},
                                         warning=function(w) {return(data.frame(
                                           estimator="item_itemspecific", baseline_name="Item-Specific",
                                           estimate=NA, se=NA, lowerci=NA, upperci=NA, p=NA, 
                                           working_cor="Exchangeable", qic=NA, qicu=NA,
                                           residual_sd=NA, had_warning=1, had_error=0
                                         ))})
    bin_me_estimate_exch <- tryCatch(get_estimates(data=sim_data_long_binomial, estimator="item_maineffects", 
                                                     working_cor="exchangeable", baseline="Y", 
                                                     baseline_adjustment="is",
                                                     covariates="baseline_response", covariate_adjustment="is", 
                                                     format="N", dist_family=binomial(link="logit")) %>%
                                         mutate(had_warning=0, had_error=0),
                                       error=function(e) {return(data.frame(
                                         estimator="item_maineffects", baseline_name="Item-Specific",
                                         estimate=NA, se=NA, lowerci=NA, upperci=NA, p=NA, 
                                         working_cor="Exchangeable", qic=NA, qicu=NA,
                                         residual_sd=NA, had_warning=0, had_error=1
                                       ))},
                                       warning=function(w) {return(data.frame(
                                         estimator="item_maineffects", baseline_name="Item-Specific",
                                         estimate=NA, se=NA, lowerci=NA, upperci=NA, p=NA, 
                                         working_cor="Exchangeable", qic=NA, qicu=NA,
                                         residual_sd=NA, had_warning=1, had_error=0
                                       ))})
      
      bin_item_estimate_blockind <- tryCatch(get_estimates(data=sim_data_long_binomial, estimator="item_itemspecific", 
                                                             working_cor="userdefined",
                                                             zcor=get_zcor(sim_data$sim_data_long, "block_ind", 
                                                                           subscale_list),
                                                             working_cor_name="Block Independent", baseline="Y", 
                                                             baseline_adjustment="is",
                                                             covariates="baseline_response", covariate_adjustment="is", format="N",
                                                             dist_family=binomial(link="logit")) %>%
                                                 mutate(had_warning=0, had_error=0),
                                               error=function(e) {return(data.frame(
                                                 estimator="item_itemspecific", baseline_name="Item-Specific",
                                                 estimate=NA, se=NA, lowerci=NA, upperci=NA, p=NA, 
                                                 working_cor="Block Independent", qic=NA, qicu=NA,
                                                 residual_sd=NA, had_warning=0, had_error=1
                                               ))},
                                               warning=function(w) {return(data.frame(
                                                 estimator="item_itemspecific", baseline_name="Item-Specific",
                                                 estimate=NA, se=NA, lowerci=NA, upperci=NA, p=NA, 
                                                 working_cor="Block Independent", qic=NA, qicu=NA,
                                                 residual_sd=NA, had_warning=1, had_error=0
                                               ))})
      bin_item_estimate_blockunstr <- tryCatch(get_estimates(data=sim_data_long_binomial, estimator="item_itemspecific", 
                                                               working_cor="userdefined",
                                                               zcor=get_zcor(sim_data$sim_data_long, "block_unstr", 
                                                                             subscale_list),
                                                               working_cor_name="Block Unstructured", baseline="Y", 
                                                               baseline_adjustment="is",
                                                               covariates="baseline_response", covariate_adjustment="is", format="N",
                                                               dist_family=binomial(link="logit")) %>%
                                                   mutate(had_warning=0, had_error=0),
                                                 error=function(e) {return(data.frame(
                                                   estimator="item_itemspecific", baseline_name="Item-Specific",
                                                   estimate=NA, se=NA, lowerci=NA, upperci=NA, p=NA, 
                                                   working_cor="Block Unstructured", qic=NA, qicu=NA,
                                                   residual_sd=NA, had_warning=0, had_error=1
                                                 ))},
                                                 warning=function(w) {return(data.frame(
                                                   estimator="item_itemspecific", baseline_name="Item-Specific",
                                                   estimate=NA, se=NA, lowerci=NA, upperci=NA, p=NA, 
                                                   working_cor="Block Unstructured", qic=NA, qicu=NA,
                                                   residual_sd=NA, had_warning=1, had_error=0
                                                 ))})
        
      bin_me_estimate_blockunstr <- tryCatch(get_estimates(data=sim_data_long_binomial, estimator="item_maineffects", 
                                                             working_cor="userdefined",
                                                             zcor=get_zcor(sim_data$sim_data_long, "block_unstr", 
                                                                           subscale_list),
                                                             working_cor_name="Block Unstructured", baseline="Y", 
                                                             baseline_adjustment="is",
                                                             covariates="baseline_response", covariate_adjustment="is", format="N",
                                                             dist_family=binomial(link="logit")) %>%
                                                 mutate(had_warning=0, had_error=0),
                                               error=function(e) {return(data.frame(
                                                 estimator="item_maineffects", baseline_name="Item-Specific",
                                                 estimate=NA, se=NA, lowerci=NA, upperci=NA, p=NA, 
                                                 working_cor="Block Unstructured", qic=NA, qicu=NA,
                                                 residual_sd=NA, had_warning=0, had_error=1
                                               ))},
                                               warning=function(w) {return(data.frame(
                                                 estimator="item_maineffects", baseline_name="Item-Specific",
                                                 estimate=NA, se=NA, lowerci=NA, upperci=NA, p=NA, 
                                                 working_cor="Block Unstructured", qic=NA, qicu=NA,
                                                 residual_sd=NA, had_warning=1, had_error=0
                                               ))})
      bin_me_estimate_blockind <- tryCatch(get_estimates(data=sim_data_long_binomial, estimator="item_maineffects", 
                                                           working_cor="userdefined",
                                                           zcor=get_zcor(sim_data$sim_data_long, "block_ind", 
                                                                         subscale_list),
                                                           working_cor_name="Block Independent", baseline="Y", 
                                                           baseline_adjustment="is",
                                                           covariates="baseline_response", covariate_adjustment="is", 
                                                           format="N", dist_family=binomial(link="logit")) %>%
                                               mutate(had_warning=0, had_error=0),
                                             error=function(e) {return(data.frame(
                                               estimator="item_maineffects", baseline_name="Item-Specific",
                                               estimate=NA, se=NA, lowerci=NA, upperci=NA, p=NA, 
                                               working_cor="Block Independent", qic=NA, qicu=NA,
                                               residual_sd=NA, had_warning=0, had_error=1
                                             ))},
                                             warning=function(w) {return(data.frame(
                                               estimator="item_maineffects", baseline_name="Item-Specific",
                                               estimate=NA, se=NA, lowerci=NA, upperci=NA, p=NA, 
                                               working_cor="Block Independent", qic=NA, qicu=NA,
                                               residual_sd=NA, had_warning=1, had_error=0
                                             ))})

    total_treatment_effect = sum(sim_data$treatment_effects_used)
    estimates_output[[i]] <- rbind(bin_item_estimate_blockunstr, bin_item_estimate_blockind,
                                      bin_item_estimate_exch, bin_item_estimate_ind, bin_me_estimate_blockunstr,
                                      bin_me_estimate_blockind, bin_me_estimate_exch, bin_me_estimate_ind, bin_standard_estimate) %>%
      mutate(
        n=n,
        treatment_mean=simulation_scenarios[[i]]$treatment_effect_mean,
        treatment_gamma=simulation_scenarios[[i]]$treatment_gamma,
        baseline_mean=simulation_scenarios[[i]]$baseline_effect_mean,
        baseline_gamma=simulation_scenarios[[i]]$baseline_gamma,
        seed=seed,
        scenario=i,
        total_treatment_effect = total_treatment_effect
      )
  }

estimates <- rbindlist(estimates_output)
write_rds(estimates, file=here("Output", "Simulation Output", "Ordinal Data", "Binomial Iterations", 
                               paste0("binomial_estimates_array", arraynum, ".rds")))

