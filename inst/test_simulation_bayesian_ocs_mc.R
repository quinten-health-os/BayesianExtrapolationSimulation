rm(list = ls())
library(RBExT)

source(system.file(paste0("conf/metrics_config.R"), package = "RBExT"))

env <- "pipeline_tests"
config_dir <- paste0(system.file(paste0("conf/", env), package = "RBExT"), "/")
case_studies_config_dir <- paste0(system.file(paste0("conf/case_studies"), package = "RBExT"), "/")
results_dir <- paste0("./results/", env, "/")
outputs_config <- yaml::read_yaml(system.file("conf/outputs_config.yml", package = "RBExT"))
ocs_filename <- outputs_config$bayesian_ocs_mc_results_filename

analysis_config <- yaml::read_yaml(paste0(config_dir, "analysis_config.yml"))
simulation_config <- yaml::read_yaml(paste0(config_dir, "simulation_config.yml"))
simulation_config$compute_frequentist_ocs <- FALSE
simulation_config$compute_bayesian_ocs_mc <- TRUE

case_studies <- scenarios_config$case_studies
methods <- scenarios_config$methods

if (simulation_config$delete_old_results) {
  folder_dir <- file.path(results_dir, "bayesian")
  if (file.exists(folder_dir)) {
    unlink(folder_dir, recursive = TRUE)
  }
}

# Loop over case studies, methods, and parallelization settings
for (case_study in case_studies) {
  for (method in methods) {
    # Initialize an empty list to store results
    results_list <- list()
    for (parallelization in list(TRUE)) {
      message(
        paste0(
          "Computing bayesian OCs for ",
          case_study,
          " with ",
          method,
          ", parallelization: ",
          parallelization
        )
      )
      flush.console()

      scenarios_config$case_studies <- case_study
      scenarios_config$methods <- method
      scenarios_config$parallelization <- parallelization

      success <- TRUE
      error_message <- NA
      start_time <- Sys.time()

      # Run the simulation with tryCatch
      ##################
      closeAllConnections()
      simulation_bayesian_ocs(
        env,
        simulation_config,
        analysis_config,
        config_dir,
        case_studies_config_dir
      )
      ##################
      # tryCatch(
      #   {
      #     closeAllConnections()
      #     simulation_bayesian_ocs(env, simulation_config, analysis_config, config_dir, case_studies_config_dir)
      #   },
      #   error = function(e) {
      #     success <<- FALSE
      #     error_message <<- e$message
      #   }
      # )

      end_time <- Sys.time()
      computation_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
      message(paste("Error message:", error_message, sep = " "))
      message(paste("Computation time:", computation_time, sep = " "))
      message()
      flush.console()

      # Store the results in the list
      formatted_time <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")
      results_list[[length(results_list) + 1]] <- data.frame(
        case_study = case_study,
        method = method,
        parallelization = parallelization,
        sampling_approximation = "NA",
        computation_time = computation_time,
        success = success,
        error_message = ifelse(is.na(error_message), "", error_message),
        n_replicates = scenarios_config$n_replicates,
        ndrift = scenarios_config$ndrift,
        n_samples_design_prior = simulation_config$n_samples_design_prior,
        run_time = formatted_time,
        stringsAsFactors = FALSE
      )
    }
    # Convert the list of results to a data frame
    results_df <- do.call(rbind, results_list)

    # Write the results to a CSV file
    directory <- paste0("./results/", env, "/bayesian/", case_study, "/", method)
    dir.create(directory, showWarnings = FALSE, recursive = TRUE)
    filename <- paste0(directory, "/", "logs.csv")
    readr::write_csv(results_df, filename)
  }
}

concatenate_simulation_results(results_dir, ocs_filename)
