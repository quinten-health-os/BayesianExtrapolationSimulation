rm(list = ls())

library(RBExT)

envs <- c("fast_cases_config", "aprepitant_approx_config", "commensurate_pp_config", "npp_config")

closeAllConnections()

outputs_config <- yaml::read_yaml(system.file("conf/outputs_config.yml", package = "RBExT"))
# Define the filenames to concatenate
file_names <- c(outputs_config$frequentist_ocs_results_filename, outputs_config$bayesian_ocs_mc_results_filename, outputs_config$bayesian_ocs_deterministic_results_filename)

results_dirs <- sapply(envs, function(env) paste0("./results/", env, "/"))

# Concatenate and save each file type
lapply(file_names, concatenate_files, results_dirs)
