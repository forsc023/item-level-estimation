## Combine simulation output from different iterations into a single dataset that can then
## be used to calculate simulation metrics. If only one iteration was run, this script
## will not be helpful.

## In this script, we will read in the estimates from the ordinal data simulations with
## gaussian family estimation, but the same thing could be done for the ordinal data
## simulations with binomial family estimation or for normal data simulations

library(here)
library(readr)
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)

data_path <- here("Output", "Simulation Output", "Ordinal Data", "Gaussian Iterations")

sim_files <- list.files(path=data_path, pattern="estimates_array[0-9]+.rds",
                        full.names=T)
sim_output <- lapply(sim_files, read_rds)

estimates <- rbindlist(sim_output)
estimates <- estimates %>%
  mutate(
    estimator_type = case_when(
      estimator=="item_itemspecific" & working_cor=="Block Unstructured" ~ "is_blunstr",
      estimator=="item_itemspecific" & working_cor=="Block Independent" ~ "is_blockind",
      estimator=="item_itemspecific" & working_cor=="Exchangeable" ~ "is_exch",
      estimator=="item_itemspecific" & working_cor=="Independence" ~ "is_ind",
      estimator=="item_maineffects" & working_cor=="Block Unstructured" ~ "me_blunstr",
      estimator=="item_maineffects" & working_cor=="Block Independent" ~ "me_blockind",
      estimator=="item_maineffects" & working_cor=="Exchangeable" ~ "me_exch",
      estimator=="item_maineffects" & working_cor=="Independence" ~ "me_ind",
      estimator=="standard" ~ "std",
      .default=estimator
    ),
    bias = estimate - total_treatment_effect,
    ci_width = upperci - lowerci,
    ci_coverage = (total_treatment_effect >= lowerci & total_treatment_effect <= upperci),
    reject = ((lowerci > 0 & upperci > 0) | (lowerci < 0 & upperci < 0))
  )

#write_rds(estimates, file=here("Output", "Simulation Output", "Ordinal Data", "ordinal_gaussian_estimates.rds"))

