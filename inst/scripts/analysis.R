rm(list = ls())

library(RBExT)

devtools::load_all()

# If the envs variable is not defined as an environment variable, define it.
envs <- ifelse(Sys.getenv("envs") != "", Sys.getenv("envs"), c("fast_cases_config"))

envs = c("combined")

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

  case_studies = c("dapagliflozin", "belimumab", "mepolizumab", "teriflunomide", "aprepitant")
  for (case_study in case_studies){
    print(case_study)
    simulation_analysis(env = env, analysis_config = analysis_config, config_dir = config_dir, frequentist_metrics = frequentist_metrics, case_studies = c(case_study),  to_compute = c( "frequentist_power_at_equivalent_tie", "frequentist_power_at_nominal_tie"), methods = c("all"))
  }
  # simulation_analysis(env, analysis_config, config_dir, frequentist_metrics, case_studies = c("mepolizumab", "teriflunomide", "dapagliflozin", "botox", "belimumab", "aprepitant"),  to_compute = c("frequentist_power_at_equivalent_tie", "frequentist_power_at_nominal_tie", "sweet_spot", "bayesian_ocs"), methods = "conditional_power_prior", "p_value_based_PP", "pooling", "RMP", "separate", "test_then_pool_difference", "test_then_pool_equivalence", "commensurate_power_prior", "EB_PP", "NPP")
}


# results_dir <- "./results/combined"
# results_freq_df <- readr::read_csv(paste0(results_dir, "/results_frequentist.csv"))
# analyze_power_gains(results_freq_df, output_path = results_dir)
# analyze_power_loss(results_freq_df, output_path = results_dir)
# analyze_power_loss_inflated_tie(results_freq_df, output_path = results_dir)
# analyze_noninflated_tie(results_freq_df, output_path = results_dir)
