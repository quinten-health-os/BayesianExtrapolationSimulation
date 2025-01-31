#' Simulate a Scenario
#'
#' @description This function simulates a scenario based on the provided inputs and returns the results of the simulation.
#'
#' @param scenario The scenario to simulate.
#' @param simulation_config Simulation configuration
#' @param freq_filename File name for the frequentist OCs
#' @param config_dir Configuration files directory
#' @param case_studies_config_dir Case studies configurations directory
#' @return The results of the simulation.
#' @export
frequentist_ocs_scenario_simulation <- function(scenario,
                                                simulation_config,
                                                scenarios_config,
                                                freq_filename,
                                                config_dir,
                                                case_studies_config_dir) {
  set.seed(simulation_config$seed)

  target_sample_size_per_arm <- scenario$target_sample_size_per_arm[[1]]
  treatment_drift <- scenario$treatment_drift[[1]]
  control_drift <- scenario$control_drift[[1]]
  case_study <- scenario$case_study[[1]]
  method <- scenario$method[[1]]
  method_parameters <- scenario$parameters[[1]]
  source_denominator <- scenario$source_denominator[[1]]
  target_to_source_std_ratio <- scenario$target_to_source_std_ratio[[1]]

  case_study_config <- yaml::read_yaml(paste0(case_studies_config_dir, case_study, ".yml"))
  mcmc_config <- yaml::read_yaml(paste0(config_dir, "/mcmc_config.yml"))

  source_data <- SourceData$new(case_study_config, source_denominator)

  model <- Model$new()

  model <- model$create(
    case_study_config = case_study_config,
    method = method,
    method_parameters = method_parameters,
    source_data = source_data,
    mcmc_config = mcmc_config
  )

  theta_0 <- case_study_config$theta_0
  critical_value <- simulation_config$critical_value
  n_replicates <- scenarios_config$n_replicates

  summary_measure_likelihood <- case_study_config$summary_measure_likelihood
  endpoint <- case_study_config$endpoint

  sampling_approximation <- case_study_config$sampling_approximation

  # Instantiate a target data object
  target_data <- TargetDataFactory$new()

  target_data <- target_data$create(
    source_data = source_data,
    case_study_config = case_study_config,
    target_sample_size_per_arm = target_sample_size_per_arm,
    control_drift = control_drift,
    treatment_drift = treatment_drift,
    summary_measure_likelihood = summary_measure_likelihood,
    target_to_source_std_ratio = target_to_source_std_ratio
  )

  # Get the current RNG state
  computation_state <- list(rng_state = .Random.seed, computation_time = 0)

  start_time <- Sys.time()

  null_space <- case_study_config$null_space

  # Estimate the frequentist OCs for the model in the scenario considered
  sim_outputs <- model$estimate_frequentist_operating_characteristics(
    theta_0 = theta_0,
    target_data = target_data,
    n_replicates = n_replicates,
    critical_value = critical_value,
    confidence_level = simulation_config$confidence_level,
    null_space = null_space,
    n_samples_quantiles_estimation = simulation_config$n_samples_quantiles_estimation,
    case_study = case_study,
    method = method
  )

  end_time <- Sys.time()

  computation_state$computation_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  sim_outputs$posterior_parameters <- format_parameters_to_json(sim_outputs$posterior_parameters)
  scenario$parallelization <- scenarios_config$parallelization

  check_colnames(scenario, expected_colnames_scenario)

  check_colnames(data.frame(sim_outputs), expected_colnames_results)

  source_data_df <- data.frame(source_data$to_dict())

  check_colnames(source_data_df, expected_colnames_source)

  target_data_df <- data.frame(lapply(target_data$to_dict(), function(x) {
    if (is.null(x)) {
      NA
    } else {
      x
    }
  }))

  expected_colnames_target <- c(
    "summary_measure_likelihood",
    "target_sample_size_per_arm",
    "target_treatment_effect",
    "target_standard_deviation",
    "target_treatment_rate",
    "target_control_rate"
  )

  check_colnames(target_data_df, expected_colnames_target)

  additional_columns <- c("sampling_approximation", "rng_state", "computation_time")

  results_frequentist_ocs <- cbind(
    scenario,
    sim_outputs,
    source_data_df,
    target_data_df,
    sampling_approximation,
    theta_0,
    null_space,
    rng_state = jsonlite::toJSON(computation_state$rng_state, pretty = TRUE),
    computation_time = computation_state$computation_time
  )
  results_frequentist_ocs$source_denominator <- unlist(results_frequentist_ocs$source_denominator)
  results_frequentist_ocs <- results_frequentist_ocs[!duplicated(names(results_frequentist_ocs), fromLast = TRUE)]

  results_frequentist_ocs$parameters <- format_parameters_to_json(results_frequentist_ocs$parameters, escape = TRUE)


  for (col in names(results_frequentist_ocs)) {
    # Check if the column is a list
    if (is.list(results_frequentist_ocs[[col]])) {
      # Attempt to simplify to either numeric or character, depending on content
      if (is.numeric(unlist(results_frequentist_ocs[[col]][1]))) {
        results_frequentist_ocs[[col]] <- as.numeric(unlist(results_frequentist_ocs[[col]]))
      } else if (is.character(unlist(results_frequentist_ocs[[col]][1]))) {
        results_frequentist_ocs[[col]] <- as.character(unlist(results_frequentist_ocs[[col]]))
      } else {
        results_frequentist_ocs[[col]] <- unlist(results_frequentist_ocs[[col]])
      }
    }
  }

  # Reorder the columns to ease reading
  # Combine all expected column names into a single vector in the desired order
  expected_colnames_order <- c(
    expected_colnames_scenario,
    expected_colnames_target,
    expected_colnames_source,
    additional_columns,
    expected_colnames_results
  )
  expected_colnames_order <- unique(expected_colnames_order)

  results_frequentist_ocs <- results_frequentist_ocs[, expected_colnames_order, drop = FALSE]

  check_confidence_intervals(results_frequentist_ocs, c(names(inference_metrics), names(frequentist_metrics)))


  # If scenarios_config$parallelization == FALSE, the csv files containing the results are built in an iterative manner, so that in case of bugs results are not lost (not possible if parallelization is used).
  if (scenarios_config$parallelization == FALSE) {
    if (!file.exists(freq_filename)) {
      # Get the column names
      fieldnames_freq <- colnames(results_frequentist_ocs)

      # Create an empty data frame with these column names
      empty_df <- data.frame(matrix(ncol = length(fieldnames_freq), nrow = 0))
      colnames(empty_df) <- fieldnames_freq

      # Create an empty data frame with these column names
      write.table(
        empty_df,
        file = freq_filename,
        sep = ",",
        col.names = TRUE,
        row.names = FALSE,
        quote = FALSE
      )
    }

    write.table(
      results_frequentist_ocs,
      file = freq_filename,
      sep = ",",
      col.names = FALSE,
      row.names = FALSE,
      append = TRUE,
      quote = TRUE
    )
  }

  return(results_frequentist_ocs)
}

#' Run simulations based on the given environment
#'
#' @param env The environment to run the simulations in
#' @param simulation_config Simulation configuration
#' @param analysis_config Analysis configuration
#' @param config_dir Configuration directory
#' @param case_studies_config_dir Case studies configuration directory
#' @return None
#' @export
#' @importFrom foreach %dopar%
simulation_frequentist_ocs <- function(env,
                                       simulation_config,
                                       scenarios_config,
                                       analysis_config,
                                       config_dir,
                                       case_studies_config_dir,
                                       logging_file_path){
  # Create scenarios and create a summary table
  cases <- simulation_scenarios(config_dir = config_dir, scenarios_config = scenarios_config)

  scenarios_table_ranges(cases, paste0("./results/", env))

  results_dir <- paste0("./results/", env, "/", "frequentist")

  case_studies <- scenarios_config$case_studies
  methods <- scenarios_config$methods

  if (!(any("separate" %in% methods))) {
    warning("Separate analysis not included in the methods.")
    futile.logger::flog.warn("Separate analysis not included in the methods.")
  }

  if (length(case_studies) == 0 & length(methods) == 0) {
    stop("Select at least one case study and one method to run the simulation.")
  }

  for (case_study in case_studies) {
    for (method in methods) {
      scenarios_config$case_studies <- case_study
      scenarios_config$methods <- method

      cases <- simulation_scenarios(config_dir = config_dir, scenarios_config = scenarios_config)

      if (nrow(cases) == 0) {
        stop("No simulation results")
      }

      outputs_config <- yaml::read_yaml(system.file("conf/outputs_config.yml", package = "RBExT"))
      freq_filename <- paste0(
        results_dir,
        "/",
        case_study,
        "/",
        method,
        "/",
        outputs_config$frequentist_ocs_results_filename
      )

      if (nrow(cases) == 0) {
        stop("No simulation results after filtering")
      }

      # Progress bar function in R
      pb <- progress::progress_bar$new(format = "Progress : [:bar] :elapsed | eta: :eta",
                                       total = nrow(cases),
                                       width = 60)
      progress <- function(n) {
        pb$tick()
      }
      opts <- list(progress = progress)

      # Create result folder if does not exist
      dir.create(
        paste0(results_dir, "/", case_study, "/", method),
        showWarnings = FALSE,
        recursive = TRUE
      )

      # Delete existing results file. If they exist and are open, close them and delete them.
      if (file.exists(freq_filename)) {
        if (simulation_config$delete_old_results){
          freq_file <- file(freq_filename, "w")
          if (isOpen(freq_file)) {
            close(freq_file)
          }
          file.remove(freq_filename)
        } else {
          next
        }
      }

      # Temporarily disable the global error handler within the tryCatch block. Otherwise there will be interference between the tryCatch block and the global error handler.
      options(error = NULL)
      if (scenarios_config$parallelization == FALSE) {
        # Loop through cases

        for (i in 1:nrow(cases)) {
          scenario <- cases[i, ]

          tryCatch({
            results <- frequentist_ocs_scenario_simulation(
              scenario = scenario,
              simulation_config = simulation_config,
              scenarios_config = scenarios_config,
              freq_filename = freq_filename,
              config_dir = config_dir,
              case_studies_config_dir = case_studies_config_dir
            )
            progress(i)
          }, error = function(e) {
            futile.logger::flog.error(paste("Error at iteration", i, ":", e$message))
            futile.logger::flog.error("Call stack:\n%s", paste(deparse(sys.calls()), collapse = "\n"))

            # Save state
            save_state(iteration = i, scenario = scenario, env = env)
            stop(e)
          })
        }
        # At the end of the loop, close file if open
        freq_file <- file(freq_filename, "a")
        if (isOpen(freq_file)) {
          close(freq_file)
        }
      } else {
        ncores <- min(parallel::detectCores() - 1, 96) # Leave on core free

        # Set a range of ports for cluster
        options(clusterPort = c(11000, 11999))

        # Register parallel backend
        cl <- parallel::makeCluster(ncores)
        doParallel::registerDoParallel(cl)

        # Define the list of libraries to load
        required_libraries <- c("devtools")

        # Export the library paths to each worker
        paths <- .libPaths()
        parallel::clusterExport(
          cl,
          varlist = c(
            "paths",
            "required_libraries",
            "frequentist_metrics",
            "inference_metrics",
            "simulation_config",
            "freq_filename",
            "config_dir",
            "case_studies_config_dir"
          ),
          envir = environment()
        )

        # Evaluate the expression to load libraries in each worker
        parallel::clusterEvalQ(cl, {
          .libPaths(paths)
          sapply(required_libraries, library, character.only = TRUE)
          devtools::load_all()
        })

        # Parallel computation over cases
        results <- foreach::foreach(
          i = 1:nrow(cases),
          .combine = "rbind",
          .options.snow = opts
        ) %dopar% {
          worker_id <- Sys.getpid()  # Get the process ID of the worker

          scenario <- cases[i, ]

          tryCatch({
            result <- frequentist_ocs_scenario_simulation(
              scenario = scenario,
              simulation_config = simulation_config,
              scenarios_config = scenarios_config,
              freq_filename = freq_filename,
              config_dir = config_dir,
              case_studies_config_dir = case_studies_config_dir
            )
            ParallelLogger::logInfo(paste("Iteration", i, " completed by worker", worker_id))
            result  # Return the result (assuming it’s a data frame)
          }, error = function(e) {
            ParallelLogger::logError(paste("Error at iteration", i, "by worker", worker_id, ":", e$message))
            ParallelLogger::logError("Call stack:\n%s", paste(deparse(sys.calls()), collapse = "\n"))

            # Save state even if there's an error
            save_state(iteration = i, scenario = scenario, worker_id = worker_id, env = env)

            stop(e)
          })
        }
        write_csv(results, freq_filename)
        parallel::stopCluster(cl)
      }

      # Restore the global error handler
      options(error = global_error_handler)

      if (tolower(case_study) == "aprepitant" ||
          method == "commensurate_power_prior") {
        # remove remaining MCMC csv files related to the run
        package_path <- system.file("", package = "RBExT")
        stan_draws_path <- paste0(package_path, "/stan/draws/", tolower(case_study), "_", method, "/")
        file.remove(list.files(stan_draws_path, full.names = TRUE))
      }
    }
  }
}
