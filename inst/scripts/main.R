rm(list = ls())

library(RBExT)

devtools::load_all() # FIXME

# If the envs variable is not defined as an environment variable, define it.
envs <- ifelse(Sys.getenv("envs") != "", Sys.getenv("envs"), c("fast_cases_config"))

envs = c("aprepitant_mcmc_config_light")
# envs = c("fast_cases_config")
#envs = c("tests")

options(readr.show_col_types = FALSE) #  FALSE : prevent the column specification message from appearing every time read_csv() is used.

concat_all_results <- TRUE
delete_stan_files <- FALSE
check_results_completeness <- TRUE

case_studies_config_dir <- paste0(system.file("conf/case_studies", package = "RBExT"), "/")

analysis_config <- yaml::read_yaml(system.file("conf/analysis_config.yml", package = "RBExT"))

simulation_config <- yaml::read_yaml(system.file("conf/simulation_config.yml", package = "RBExT"))

source(system.file(paste0("conf/metrics_config.R"), package = "RBExT"))

closeAllConnections()

if (delete_stan_files == TRUE) {
  # Delete all Stan files (to make sure models are recompiled)
  stan_files <- list.files(path = system.file("stan", package = "RBExT"), pattern = "\\.stan$", full.names = TRUE)
  exe_files <- list.files(path = system.file("stan", package = "RBExT"), pattern = "\\.exe$", full.names = TRUE)

  # Combine the lists of .stan and .exe files and remove them
  files_to_delete <- c(stan_files, exe_files)
  file.remove(files_to_delete)
}

outputs_config <- yaml::read_yaml(system.file("conf/outputs_config.yml", package = "RBExT"))
ocs_filename <- outputs_config$frequentist_ocs_results_filename

for (env in envs) {
  # if (!(env %in% c("aprepitant_approx_config", "aprepitant_mcmc_config", "commensurate_pp_config", "fast_cases_config", "pdccpp_config", "npp_config", "tests", "full", "pipeline_tests", "ttp"))) {
  #   stop("Invalid environment name.")
  # }

  # Set up the log file location
  dir.create(paste0("./logs/", env, "/checkpoints/"), showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0("./logs/", env, "/error_logs/"), showWarnings = FALSE, recursive = TRUE)
  # Empty the directory
  unlink(paste0("./logs/", env, "/checkpoints/*"), recursive = TRUE)
  unlink(paste0("./logs/", env, "/error_logs/*"), recursive = TRUE)

  LOGGING_FILE_PATH <- generate_log_filename(base_name = paste0("./logs/", env, "/error_logs/error_log.log"), suffix_type = "timestamp")

  # We must add trailing slashes manually as otherwise system.file will drop them.
  config_dir <- paste0(system.file(paste0("conf/", env), package = "RBExT"), "/")

  results_dir <- paste0("./results/", env, "/")
  scenarios_config <- yaml::read_yaml(paste0(config_dir, "scenarios_config.yml"))

  # Define the parallel logger
  if (scenarios_config$parallelization){
    # Parallel logger
    ParallelLogger::registerLogger(ParallelLogger::createLogger(
      name = "ParLogger",
      threshold = "INFO",
      appenders = list(
        ParallelLogger::createConsoleAppender(
          layout =  ParallelLogger::layoutSimple
        ),
        ParallelLogger::createFileAppender(
          layout =  ParallelLogger::layoutParallel,
          fileName = LOGGING_FILE_PATH
        )
      )
    ))


    ParallelLogger::registerLogger(ParallelLogger::createLogger(name = "DEFAULT_ERRORREPORT_LOGGER",
                                threshold = "FATAL",
                                appenders = list(ParallelLogger::createFileAppender(layout =  ParallelLogger::layoutErrorReport,
                                                                    fileName = paste0("./logs/", env, "/error_report.txt"),
                                                                    overwrite = TRUE,
                                                                    expirationTime = 60))))
  } else {
    # Set up the new log file
    futile.logger::flog.appender(futile.logger::appender.file(LOGGING_FILE_PATH))

    # Set the global error handler
    options(error = global_error_handler)
  }

  if (simulation_config$compute_frequentist_ocs == TRUE) {
    if (simulation_config$delete_old_results) {
      unlink(results_dir, recursive = TRUE, force = TRUE)
      dir.create(results_dir)  # This recreates the empty directory
    }

    simulation_frequentist_ocs(env = env,
    simulation_config = simulation_config,
    scenarios_config = scenarios_config,
    analysis_config = analysis_config,
    config_dir = config_dir,
    case_studies_config_dir = case_studies_config_dir,
    logging_file_path = LOGGING_FILE_PATH)

    if (check_results_completeness){
      # check_simulation_completeness(results_dir = results_dir, ocs_filename = ocs_filename, scenarios_config = scenarios_config, config_dir)
    }

    # Concatenate case_study/method results into a single file for the environment.
    concatenate_simulation_results(results_dir = results_dir, ocs_filename = ocs_filename)

    # The analysis is performed on the results concatenated at the level of the environment.
    simulation_analysis(env, analysis_config, config_dir, frequentist_metrics)
  }

  if (simulation_config$compute_bayesian_ocs_mc == TRUE) {
    if (simulation_config$delete_old_results) {
      folder_dir <- file.path(results_dir, "bayesian")
      if (file.exists(folder_dir)) {
        unlink(folder_dir, recursive = TRUE)
      }
    }

    simulation_bayesian_ocs(
      env = env,
      simulation_config = simulation_config,
      scenarios_config = scenarios_config,
      analysis_config = analysis_config,
      config_dir = config_dir,
      case_studies_config_dir = case_studies_config_dir
    )
    ocs_filename <- outputs_config$bayesian_ocs_mc_results_filename
    concatenate_simulation_results(results_dir = results_dir, ocs_filename = ocs_filename)
  }
}
#
#
# if (concat_all_results == TRUE) {
#   closeAllConnections()
#
#   envs <-  c("aprepitant_mcmc_config", "commensurate_config", "fast_cases_config", "commensurate_config", "aprepitant_approx_config", "npp_config")
#   file_names <- c(outputs_config$frequentist_ocs_results_filename, "sweet_spot.csv", outputs_config$bayesian_ocs_deterministic_results_filename)
#
#   results_dirs <- lapply(envs, function(env) paste0("./results/", env, "/"))
#
#   # Concatenate and save each file type
#   lapply(file_names, concatenate_files, results_dirs, "./results/combined/")
# }
#

