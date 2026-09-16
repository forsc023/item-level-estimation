## Run item level GEE simulation for correlated multivariate normal data. This script
## reads in the simulation scenarios, generates a dataset for each scenario and carries
## out the treatment effect estimation for one simulation iteration. Simulations for the
## paper were run in the Minnesota Supercomputing Institute using a job array to run
## the 8,000 iterations in parallel. This script only runs one iteration.

## Change "arraynum" to run other iterations.
arraynum <- 1

library(readr)
library(here)
library(data.table)
source(here("Scripts", "analysis_functions.R"))
simulation_scenarios <- read_rds(here("Output", "Simulation Scenarios", "simulation_scenarios_mvn.rds"))

estimates_output <- vector(mode="list", length=length(simulation_scenarios))
## Iterate over all the simulation scenarios
for (i in 1:length(simulation_scenarios)) {
  seed <- 100000*i + arraynum
  set.seed(seed)
  target_effects <- simulation_scenarios[[i]]$treatment_effects
  target_baseline_effects <- simulation_scenarios[[i]]$baseline_effects
  n <- simulation_scenarios[[i]]$n
  support_list <- simulation_scenarios[[i]]$support_list
  beta0 <- simulation_scenarios[[i]]$beta0
  within_corr <- simulation_scenarios[[i]]$within_corr
  without_corr <- simulation_scenarios[[i]]$without_corr
  subscale_list <- simulation_scenarios[[i]]$subscale_list
  subscale_sizes <- simulation_scenarios[[i]]$subscale_sizes
  n_items <- simulation_scenarios[[i]]$n_items
  
  ## Variance for baseline and follow up responses
  ## To standardize the baseline coefficients, we will multiply them by the ratio of these.
  ## In general, we would want to do this at the item level; ratio of the residual variance
  ## associated with each item. In this simulation setting, however, all the items have the
  ## same variance so we can just stick with one number.
  sigma2_baseline <- 1
  sigma2_post <- 1
  standardized_baseline_effects <- target_baseline_effects*(sigma2_baseline/sigma2_post)
  standardized_baseline_sd <- sd(standardized_baseline_effects)
  
  sim_data <- simulate_dataset_mvn(n=n, n_items=n_items, beta0=beta0, 
                                   baseline_effects=target_baseline_effects,
                                   treatment_effects=target_effects, sigma2_baseline=sigma2_baseline,
                                   sigma2_post=sigma2_post,
                                   within_corr=within_corr, without_corr=without_corr,
                                   subscale_sizes=subscale_sizes)
  
  item_estimate_blockunstr <- get_estimates(data=sim_data$sim_data_long, estimator="item_itemspecific", 
                                            working_cor="userdefined",
                                            zcor=get_zcor(sim_data$sim_data_long, "block_unstr", 
                                                          subscale_list),
                                            working_cor_name="Block Unstructured", baseline="Y", 
                                            baseline_adjustment="is",
                                            covariates="baseline_response", covariate_adjustment="is", format="N")
  
  item_estimate_blockind <- get_estimates(data=sim_data$sim_data_long, estimator="item_itemspecific", 
                                          working_cor="userdefined",
                                          zcor=get_zcor(sim_data$sim_data_long, "block_ind", 
                                                        subscale_list),
                                          working_cor_name="Block Independent", baseline="Y", 
                                          baseline_adjustment="is",
                                          covariates="baseline_response", covariate_adjustment="is", 
                                          format="N")
  
  item_estimate_exch <- get_estimates(data=sim_data$sim_data_long, estimator="item_itemspecific", 
                                      working_cor="exchangeable", baseline="Y", 
                                      baseline_adjustment="is",
                                      covariates="baseline_response", covariate_adjustment="is", 
                                      format="N")
  
  item_estimate_ind <- get_estimates(data=sim_data$sim_data_long, estimator="item_itemspecific", 
                                     working_cor="independence", baseline="Y", 
                                     baseline_adjustment="is",
                                     covariates="baseline_response", covariate_adjustment="is", 
                                     format="N")
  
  ## Item level - main effects estimate
  me_estimate_blockunstr <- get_estimates(data=sim_data$sim_data_long, estimator="item_maineffects", 
                                          working_cor="userdefined",
                                          zcor=get_zcor(sim_data$sim_data_long, "block_unstr", 
                                                        subscale_list),
                                          working_cor_name="Block Unstructured", baseline="Y", 
                                          baseline_adjustment="is",
                                          covariates="baseline_response", covariate_adjustment="is", format="N")
  
  me_estimate_blockind <- get_estimates(data=sim_data$sim_data_long, estimator="item_maineffects", 
                                        working_cor="userdefined",
                                        zcor=get_zcor(sim_data$sim_data_long, "block_ind", 
                                                      subscale_list),
                                        working_cor_name="Block Independent", baseline="Y", 
                                        baseline_adjustment="is",
                                        covariates="baseline_response", covariate_adjustment="is", 
                                        format="N")
  
  me_estimate_exch <- get_estimates(data=sim_data$sim_data_long, estimator="item_maineffects", 
                                    working_cor="exchangeable", baseline="Y", 
                                    baseline_adjustment="is",
                                    covariates="baseline_response", covariate_adjustment="is", 
                                    format="N")
  
  me_estimate_ind <- get_estimates(data=sim_data$sim_data_long, estimator="item_maineffects", 
                                   working_cor="independence", baseline="Y", 
                                   baseline_adjustment="is",
                                   covariates="baseline_response", covariate_adjustment="is", 
                                   format="N")
  
  standard_estimate <- get_estimates(data=sim_data$sim_data_wide, estimator="standard", 
                                     working_cor="none", baseline="Y",
                                     baseline_adjustment="me", covariates="baseline_score",
                                     covariate_adjustment="me", format="N")
  total_treatment_effect <- sum(target_effects)
  estimates <- rbind(item_estimate_blockunstr, item_estimate_blockind,
                     item_estimate_exch, item_estimate_ind, me_estimate_blockunstr,
                     me_estimate_blockind, me_estimate_exch, me_estimate_ind, standard_estimate) %>%
    mutate(
      n=n,
      treatment_mean=simulation_scenarios[[i]]$treatment_effect_mean,
      treatment_gamma=simulation_scenarios[[i]]$treatment_gamma,
      baseline_mean=simulation_scenarios[[i]]$baseline_effect_mean,
      baseline_gamma=simulation_scenarios[[i]]$baseline_gamma,
      standardized_baseline_sd=standardized_baseline_sd,
      sigma2_baseline = sigma2_baseline,
      sigma2_post = sigma2_post,
      sigma2_ratio = sigma2_baseline/sigma2_post,
      seed=seed,
      scenario=i,
      total_treatment_effect=total_treatment_effect
    )
  estimates_output[[i]] <- estimates
}
result <- rbindlist(estimates_output)
write_rds(result, file=here("Output", "Simulation Output", "Normal Data", "Iterations", 
               paste0("estimates_array", arraynum, ".rds")))
