This repository contains simulation code used in the paper "Efficient Item-Level Treatment Effect Estimation for Questionnaire-Derived Outcomes" by Jamie Forschmiedt, Jacob Meyer, Jeni Lansing, and Thomas Murray.

Scripts include:
1. "sim_scenarios_ordinal.R", which generates the simulation scenarios used in our ordinal-data simulation study; 
2. "sim_scenarios_mvn.R", which generates the scenarios used in our normal-data simulation study;
3. "run_sim_ordinal_gaussian.R", which generates an ordinal dataset for each of the ordinal-data scenarios and calculates treatment effect estimates with the estimators described in the paper using Gaussian-family estimation;
4. "run_sim_ordinal_binomial.R", which generates an ordinal dataset for each of the ordinal-data scenarios and calculates treatment effect estimates with the estimators described in the paper using binomial-family estimation;
5. "run_sim_mvn.R", which generates a dataset with normally distributed outcomes for each of the normal-data scenarios and calculates treatment effect estimates with the estimators described in the paper;
6. "read_sim_output.R", which combines output from different iterations;
7. "analysis_functions.R", which contains functions to generate data and carry out the item-level estimation proposed in the paper, as well as other helper functions.

Some example output for individual iterations is provided in the "Output" sub-folder.

In the paper, we ran 8,000 iterations for each simulation scenario in parallel using job arrays on the Minnesota Supercomputing Institute. As such, the code provided here in "run_sim_ordinal_gaussian.R", "run_sim_ordinal_binomial.R", and "run_sim_mvn.R" only runs a single iteration rather than directly replicating the entire simulation study. 
To replicate the results described in the paper, this code can be run repeatedly, changing the "arraynum" parameter (which corresponds to the random seed) at the top of the aforementioned scripts to values 1-8000.
