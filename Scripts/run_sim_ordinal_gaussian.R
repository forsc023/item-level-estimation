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
arraynum <- 2

library(readr)
library(here)
library(data.table)

source(here("Scripts", "analysis_functions.R"))
#source(here("Scripts", "sim_scenarios_ordinal.R"))

simulation_scenarios <- read_rds(here("Output", "Simulation Scenarios", "simulation_scenarios_ordinal.rds"))

response_summary <- vector(mode="list", length=length(simulation_scenarios))
estimates_output <- vector(mode="list", length=length(simulation_scenarios))

## Iterate over all the simulation scenarios
## Will take a few minutes to run all the scenarios
for (i in 1:length(simulation_scenarios)) {
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
  
  ## We will generate alpha1 (the coefficients in the cumulative logit model corresponding to the
  ## treatment effect) within the data generation function, by first generating the baseline data,
  ## finding the mean baseline response for each item, and generate the treatment effect coefficients
  ## to correspond to the desired linear treatment effect given the mean baseline response.
  sim_data <- simulate_dataset(n=n, n_items=n_items, alpha0=alpha0, alpha2=alpha2, 
                               support_list=support_list, treatment_effects=target_effects,
                               within_corr=within_corr, without_corr=without_corr, 
                               subscale_sizes=subscale_sizes)
  
  ## Calculate mean response by arm & by baseline response
  ## Can also recover baseline response distribution from this
  sim_data_long <- sim_data$sim_data_long
  
  baseline_var_data <- sim_data_long %>% group_by(item) %>% 
    summarize(
      baseline_mean = mean(baseline_response),
      baseline_var = var(baseline_response)
    ) %>% ungroup() %>% data.frame()
  
  response_summary[[i]] <- sim_data_long %>%
    group_by(baseline_response, item, trt) %>%
    summarize(mean_response=mean(response), n=n()) %>%
    ungroup() %>%
    data.frame() %>%
    mutate(scenario=i)
  
  ## Gaussian Family Estimates
  ## Item level - item specific estimates
  item_estimate_blockunstr_raw <- get_estimates(data=sim_data$sim_data_long, estimator="item_itemspecific", 
                                            working_cor="userdefined",
                                            zcor=get_zcor(sim_data$sim_data_long, "block_unstr", 
                                                          subscale_list),
                                            working_cor_name="Block Unstructured", baseline="Y", 
                                            baseline_adjustment="is",
                                            covariates="baseline_response", covariate_adjustment="is", 
                                            format="N", var_ratio=T)
  item_estimate_blockunstr <- item_estimate_blockunstr_raw$result %>%
    mutate(
      standardized_baseline_sd_vratio = sd(item_estimate_blockunstr_raw$residual_var_data$resid_var_ratio*target_baseline_effects),
      standardized_baseline_sd_sdratio = sd(sqrt(item_estimate_blockunstr_raw$residual_var_data$resid_var_ratio)*target_baseline_effects)
    )
  
  item_estimate_blockind_raw <- get_estimates(data=sim_data$sim_data_long, estimator="item_itemspecific", 
                                          working_cor="userdefined",
                                          zcor=get_zcor(sim_data$sim_data_long, "block_ind", 
                                                        subscale_list),
                                          working_cor_name="Block Independent", baseline="Y", 
                                          baseline_adjustment="is",
                                          covariates="baseline_response", covariate_adjustment="is", 
                                          format="N", var_ratio=T)
  item_estimate_blockind <- item_estimate_blockind_raw$result %>%
    mutate(
      standardized_baseline_sd_vratio = sd(item_estimate_blockind_raw$residual_var_data$resid_var_ratio*target_baseline_effects),
      standardized_baseline_sd_sdratio = sd(sqrt(item_estimate_blockind_raw$residual_var_data$resid_var_ratio)*target_baseline_effects)
    )
  
  item_estimate_exch_raw <- get_estimates(data=sim_data$sim_data_long, estimator="item_itemspecific", 
                                      working_cor="exchangeable", baseline="Y", 
                                      baseline_adjustment="is",
                                      covariates="baseline_response", covariate_adjustment="is", 
                                      format="N", var_ratio=T)
  item_estimate_exch <- item_estimate_exch_raw$result %>%
    mutate(
      standardized_baseline_sd_vratio = sd(item_estimate_exch_raw$residual_var_data$resid_var_ratio*target_baseline_effects),
      standardized_baseline_sd_sdratio = sd(sqrt(item_estimate_exch_raw$residual_var_data$resid_var_ratio)*target_baseline_effects)
    )
  
  item_estimate_ind_raw <- get_estimates(data=sim_data$sim_data_long, estimator="item_itemspecific", 
                                     working_cor="independence", baseline="Y", 
                                     baseline_adjustment="is",
                                     covariates="baseline_response", covariate_adjustment="is", 
                                     format="N", var_ratio=T)
  item_estimate_ind <- item_estimate_ind_raw$result %>%
    mutate(
      standardized_baseline_sd_vratio = sd(item_estimate_ind_raw$residual_var_data$resid_var_ratio*target_baseline_effects),
      standardized_baseline_sd_sdratio = sd(sqrt(item_estimate_ind_raw$residual_var_data$resid_var_ratio)*target_baseline_effects)
    )
  
  ## Item level - main effects estimate
  me_estimate_blockunstr_raw <- get_estimates(data=sim_data$sim_data_long, estimator="item_maineffects", 
                                          working_cor="userdefined",
                                          zcor=get_zcor(sim_data$sim_data_long, "block_unstr", 
                                                        subscale_list),
                                          working_cor_name="Block Unstructured", baseline="Y", 
                                          baseline_adjustment="is",
                                          covariates="baseline_response", covariate_adjustment="is", 
                                          format="N", var_ratio=T)
  me_estimate_blockunstr <- me_estimate_blockunstr_raw$result %>%
    mutate(
      standardized_baseline_sd_vratio = sd(me_estimate_blockunstr_raw$residual_var_data$resid_var_ratio*target_baseline_effects),
      standardized_baseline_sd_sdratio = sd(sqrt(me_estimate_blockunstr_raw$residual_var_data$resid_var_ratio)*target_baseline_effects)
    )
  
  
  me_estimate_blockind_raw <- get_estimates(data=sim_data$sim_data_long, estimator="item_maineffects", 
                                        working_cor="userdefined",
                                        zcor=get_zcor(sim_data$sim_data_long, "block_ind", 
                                                      subscale_list),
                                        working_cor_name="Block Independent", baseline="Y", 
                                        baseline_adjustment="is",
                                        covariates="baseline_response", covariate_adjustment="is", 
                                        format="N", var_ratio=T)
  me_estimate_blockind <- me_estimate_blockind_raw$result %>%
    mutate(
      standardized_baseline_sd_vratio = sd(me_estimate_blockind_raw$residual_var_data$resid_var_ratio*target_baseline_effects),
      standardized_baseline_sd_sdratio = sd(sqrt(me_estimate_blockind_raw$residual_var_data$resid_var_ratio)*target_baseline_effects)
    )
  
  me_estimate_exch_raw <- get_estimates(data=sim_data$sim_data_long, estimator="item_maineffects", 
                                    working_cor="exchangeable", baseline="Y", 
                                    baseline_adjustment="is",
                                    covariates="baseline_response", covariate_adjustment="is", 
                                    format="N", var_ratio=T)
  me_estimate_exch <- me_estimate_exch_raw$result %>%
    mutate(
      standardized_baseline_sd_vratio = sd(me_estimate_exch_raw$residual_var_data$resid_var_ratio*target_baseline_effects),
      standardized_baseline_sd_sdratio = sd(sqrt(me_estimate_exch_raw$residual_var_data$resid_var_ratio)*target_baseline_effects)
    )
  
  me_estimate_ind_raw <- get_estimates(data=sim_data$sim_data_long, estimator="item_maineffects", 
                                   working_cor="independence", baseline="Y", 
                                   baseline_adjustment="is",
                                   covariates="baseline_response", covariate_adjustment="is", 
                                   format="N", var_ratio=T)
  me_estimate_ind <- me_estimate_ind_raw$result %>%
    mutate(
      standardized_baseline_sd_vratio = sd(me_estimate_ind_raw$residual_var_data$resid_var_ratio*target_baseline_effects),
      standardized_baseline_sd_sdratio = sd(sqrt(me_estimate_ind_raw$residual_var_data$resid_var_ratio)*target_baseline_effects)
    )
  
  ## Total Score estimate
  standard_estimate <- get_estimates(data=sim_data$sim_data_wide, estimator="standard", working_cor="none", 
                                     baseline="Y",
                                     baseline_adjustment="me", covariates="baseline_score",
                                     covariate_adjustment="me", format="N") %>%
    mutate(standardized_baseline_sd_vratio=NA, standardized_baseline_sd_sdratio=NA)
  
  total_treatment_effect = sum(sim_data$treatment_effects_used)
  estimates <- rbind(item_estimate_blockunstr, item_estimate_blockind,
                     item_estimate_exch, item_estimate_ind, me_estimate_blockunstr,
                     me_estimate_blockind, me_estimate_exch, me_estimate_ind, standard_estimate) %>%
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
  estimates_output[[i]] <- estimates
}
result <- rbindlist(estimates_output)
write_rds(result, file=here("Output", "Simulation Output", "Ordinal Data",
                                      "Gaussian Iterations", paste0("gaussian_estimates_array", arraynum, ".rds")))
