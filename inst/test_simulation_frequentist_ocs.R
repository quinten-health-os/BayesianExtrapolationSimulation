rm(list = ls())
library(RBExT)

source(system.file(paste0("conf/metrics_config.R"), package = "RBExT"))

env <- "pipeline_tests"
config_dir <- paste0(system.file(paste0("conf/", env), package = "RBExT"), "/")
case_studies_config_dir <- paste0(system.file(paste0("conf/case_studies"), package = "RBExT"), "/")
results_dir <- paste0("./results/", env, "/")
outputs_config <- yaml::read_yaml(system.file("conf/outputs_config.yml", package = "RBExT"))
ocs_filename <- outputs_config$frequentist_ocs_results_filename

analysis_config <- yaml::read_yaml(paste0(config_dir, "analysis_config.yml"))
simulation_config <- yaml::read_yaml(paste0(config_dir, "simulation_config.yml"))
simulation_config$compute_frequentist_ocs <- TRUE
simulation_config$compute_bayesian_ocs <- FALSE

case_studies <- scenarios_config$case_studies
methods <- scenarios_config$methods

# Possible sampling approaches for each case study
sampling_approximations <- list(
  aprepitant = c(FALSE),
  belimumab = c(TRUE, FALSE),
  botox = c(FALSE),
  dapagliflozin = c(FALSE),
  mepolizumab = c(TRUE, FALSE),
  teriflunomide = c(TRUE, FALSE)
)

sampling_approximations <- list(
  aprepitant = c(FALSE),
  belimumab = c(FALSE),
  botox = c(FALSE),
  dapagliflozin = c(FALSE),
  mepolizumab = c(FALSE),
  teriflunomide = c(FALSE)
)

if (simulation_config$delete_old_results) {
  folder_dir <- file.path(results_dir, "frequentist")
  if (file.exists(folder_dir)) {
    unlink(folder_dir, recursive = TRUE)
  }
}

# Loop over case studies, methods, and parallelization settings
for (case_study in case_studies) {
  for (method in methods) {
    # Initialize an empty list to store results
    results_list <- list()
    for (parallelization in list(FALSE)) {
      scenarios_config$case_studies <- case_study
      scenarios_config$methods <- method
      scenarios_config$parallelization <- parallelization

      for (sampling_approximation in sampling_approximations[[case_study]]) {
        case_study_config <- yaml::read_yaml(paste0(case_studies_config_dir, case_study, ".yml"))
        # Update sampling_approximation in the case_study_config file
        simulation_config$sampling_approximation <- sampling_approximation

        message(
          paste0(
            "Computing frequentist OCs for ",
            case_study,
            " with ",
            method,
            ", parallelization: ",
            parallelization,
            ", sampling approximation: ",
            sampling_approximation
          )
        )
        flush.console()

        success <- TRUE
        error_message <- NA
        start_time <- Sys.time()

        # Run the simulation with tryCatch
        ##################
        closeAllConnections()
        simulation_frequentist_ocs(
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
        #     simulation_frequentist_ocs(env, simulation_config, analysis_config, config_dir, case_studies_config_dir)
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
          sampling_approximation = sampling_approximation,
          computation_time = computation_time,
          success = success,
          error_message = error_message,
          n_replicates = scenarios_config$n_replicates,
          ndrift = scenarios_config$ndrift,
          n_samples_design_prior = simulation_config$n_samples_design_prior,
          run_time = formatted_time,
          stringsAsFactors = FALSE
        )
      }
    }
    # Convert the list of results to a data frame
    results_df <- do.call(rbind, results_list)

    # Write the results to a CSV file
    directory <- paste0("./results/", env, "/frequentist/", case_study, "/", method)
    dir.create(directory, showWarnings = FALSE, recursive = TRUE)
    filename <- paste0(directory, "/", "logs.csv")
    readr::write_csv(results_df, filename)
  }
}

concatenate_simulation_results(results_dir, ocs_filename)
# simulation_analysis(env, analysis_config, config_dir, frequentist_metrics)
