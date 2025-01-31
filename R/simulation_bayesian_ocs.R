library(foreach)

#' Simulate a Scenario
#'
#' @description This function simulates a scenario based on the provided inputs and returns the results of the simulation.
#'
#' @param scenario The scenario to simulate.
#' @param simulation_config Simulation configuration
#' @param bayes_filename Results were results are stored for Bayesian OCs
#' @param config_dir Configuration files directory
#' @return The results of the simulation.
#' @export
bayesian_ocs_scenario_simulation <- function(scenario,
                                             simulation_config,
                                             bayes_filename,
                                             config_dir) {
  set.seed(simulation_config$seed)

  target_sample_size_per_arm <- scenario$target_sample_size_per_arm[[1]]
  case_study <- scenario$case_study[[1]]
  method <- scenario$method[[1]]
  method_parameters <- scenario$parameters[[1]]
  source_denominator <- scenario$source_denominator[[1]]
  target_to_source_std_ratio <- scenario$target_to_source_std_ratio[[1]]

  case_study_config <- yaml::read_yaml(paste0(case_studies_config_dir, case_study, ".yml"))
  mcmc_config <- yaml::read_yaml(paste0(config_dir, "/mcmc_config.yml"))

  source_data <- SourceData$new(case_study_config = case_study_config,
                                source_denominator = source_denominator)

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

  json_parameters <- format_parameters_to_json(scenario$parameters)

  # Estimate the Bayesian OCs for the model in the scenario considered
  design_prior_types <- cbind("ui_design_prior", "analysis_prior", "source_posterior")

  # Get the current RNG state
  computation_state <- list(rng_state = .Random.seed, computation_time = 0)

  for (design_prior_type in design_prior_types) {
    results_bayesian_ocs <- estimate_bayesian_ocs(
      scenario = scenario,
      case_study_config = case_study_config,
      source_data = source_data,
      target_sample_size_per_arm = target_sample_size_per_arm,
      model = model,
      method_parameters = method_parameters,
      theta_0 = theta_0,
      n_replicates = n_replicates,
      critical_value = critical_value,
      design_prior_type = design_prior_type,
      json_parameters = json_parameters,
      simulation_config = simulation_config,
      target_to_source_std_ratio = target_to_source_std_ratio,
      mcmc_config = mcmc_config
    )

    results_bayesian_ocs$source_denominator <- unlist(results_bayesian_ocs$source_denominator)

    results_bayesian_ocs$rng_state <- jsonlite::toJSON(computation_state$rng_state, pretty = TRUE)
    results_bayesian_ocs$computation_time <- computation_state$computation_time

    if (design_prior_type == "ui_design_prior") {
      results_bayesian_ocs_ui_design_prior <- results_bayesian_ocs
    } else if (design_prior_type == "analysis_prior") {
      results_bayesian_ocs_analysis_prior <- results_bayesian_ocs
    } else if (design_prior_type == "source_posterior") {
      results_bayesian_ocs_source_posterior <- results_bayesian_ocs
    }
  }

  combined_results <- rbind(
    results_bayesian_ocs_ui_design_prior,
    results_bayesian_ocs_analysis_prior,
    results_bayesian_ocs_source_posterior
  )

  # If scenarios_config$parallelization == FALSE, the csv files containing the results are built in an iterative manner, so that in case of bugs results are not lost (not possible if parallelization is used).
  if (scenarios_config$parallelization == FALSE) {
    if (!file.exists(bayes_filename)) {
      # Get the column names
      fieldnames_bayes <- colnames(results_bayesian_ocs)

      # Create an empty data frame with these column names
      empty_df <- data.frame(matrix(ncol = length(fieldnames_bayes), nrow = 0))
      colnames(empty_df) <- fieldnames_bayes

      # Create an empty data frame with these column names
      write.table(
        empty_df,
        file = bayes_filename,
        sep = ",",
        col.names = TRUE,
        row.names = FALSE,
        quote = FALSE
      )
    }

    bayes_file <- file(bayes_filename, "a")

    write.table(
      results_bayesian_ocs,
      file = bayes_filename,
      sep = ",",
      col.names = FALSE,
      row.names = FALSE,
      append = TRUE,
      quote = TRUE
    )
  }
  return(combined_results)
}

#' Simulation Bayesian OCs
#'
#' @description This function estimates the Bayesian operating characteristics for the model in the scenario considered.
#' @param scenario Simulation scenario
#' @param case_study_config The case study configuration.
#' @param source_data The source data.
#' @param target_sample_size_per_arm The target sample size per arm.
#' @param model The model.
#' @param method_parameters The method parameters.
#' @param theta_0 The true parameter value.
#' @param n_replicates The number of replicates.
#' @param critical_value The critical value.
#' @param computation_state The computation state.
#' @param design_prior_type The design prior type.
#' @param json_parameters Method parameters in json format
#' @param simulation_config Simulation configuration
#' @param target_to_source_std_ratio Ratio between the source and target study sampling standard deviation.
#' @param mcmc_config MCMC configuration
#' @return The results of the Bayesian operating characteristics estimation.
estimate_bayesian_ocs <- function(scenario,
                                  case_study_config,
                                  source_data,
                                  target_sample_size_per_arm,
                                  model,
                                  method_parameters,
                                  theta_0,
                                  n_replicates,
                                  critical_value,
                                  computation_state,
                                  design_prior_type,
                                  json_parameters,
                                  simulation_config,
                                  target_to_source_std_ratio,
                                  mcmc_config = NULL) {
  design_prior <- DesignPrior$new()

  design_prior <- design_prior$create(
    design_prior_type = design_prior_type,
    model = model,
    source_data = source_data,
    case_study_config = case_study_config,
    simulation_config = simulation_config,
    mcmc_config = mcmc_config
  )
  start_time <- Sys.time()

  bayesian_ocs <- model$estimate_bayesian_operating_characteristics(
    design_prior = design_prior,
    theta_0 = theta_0,
    source_data = source_data,
    n_replicates = n_replicates,
    critical_value = critical_value,
    confidence_level = simulation_config$confidence_level,
    null_space = case_study_config$null_space,
    n_samples_design_prior = simulation_config$n_samples_design_prior,
    target_sample_size_per_arm = target_sample_size_per_arm,
    case_study_config = case_study_config,
    target_to_source_std_ratio = target_to_source_std_ratio,
    simulation_config = simulation_config
  )
  end_time <- Sys.time()
  computation_state$computation_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  results_bayesian_ocs <- cbind(scenario,
                                design_prior = design_prior_type,
                                bayesian_ocs,
                                source_data$to_dict(),
                                computation_state)

  results_bayesian_ocs$parameters <- json_parameters
  results_bayesian_ocs <- results_bayesian_ocs[!duplicated(names(results_bayesian_ocs), fromLast = TRUE)]
  for (col in names(results_bayesian_ocs)) {
    # Check if the column is a list
    if (is.list(results_bayesian_ocs[[col]])) {
      # Attempt to simplify to either numeric or character, depending on content
      if (col == "parameters") {
        results_bayesian_ocs[[col]] <- as.character(unlist(results_bayesian_ocs[[col]]))
      } else if (is.numeric(unlist(results_bayesian_ocs[[col]][1]))) {
        results_bayesian_ocs[[col]] <- as.numeric(unlist(results_bayesian_ocs[[col]]))
      } else if (is.character(unlist(results_bayesian_ocs[[col]][1]))) {
        results_bayesian_ocs[[col]] <- as.character(unlist(results_bayesian_ocs[[col]]))
      } else {
        results_bayesian_ocs[[col]] <- unlist(results_bayesian_ocs[[col]])
      }
    }
  }

  return(results_bayesian_ocs)
}

#' Run simulations based on the given environment
#'
#' @param env The environment to run the simulations in
#' @param simulation_config Simulation configuration
#' @param analysis_config Analysis configuration
#' @param config_dir Configuration directory
#' @param case_studies_config_dir Case studies configuration directory
#'
#' @return None
simulation_bayesian_ocs <- function(env,
                                    simulation_config,
                                    scenarios_config,
                                    analysis_config,
                                    config_dir,
                                    case_studies_config_dir) {
  results_dir <- paste0("./results/", env, "/", "bayesian")

  case_studies <- scenarios_config$case_studies
  methods <- scenarios_config$methods

  for (case_study in case_studies) {
    for (method in methods) {
      scenarios_config$case_studies <- case_study
      scenarios_config$methods <- method

      cases <- simulation_scenarios(config_dir = config_dir, scenarios_config = scenarios_config)

      if (nrow(cases) == 0) {
        stop("No simulation results")
      }

      cases <- dplyr::select(
        cases,
        -c(
          "drift",
          "control_drift",
          "treatment_drift",
          "target_treatment_effect"
        )
      )
      cases <- unique(cases)

      outputs_config <- yaml::read_yaml(system.file("conf/outputs_config.yml", package = "RBExT"))

      bayes_filename <- paste0(
        results_dir,
        "/",
        case_study,
        "/",
        method,
        "/",
        outputs_config$bayesian_ocs_mc_results_filename
      )

      if (nrow(cases) == 0) {
        stop("No simulation results after filtering")
      }

      set.seed(simulation_config$seed)

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
      if (file.exists(bayes_filename)) {
        bayes_file <<- file(bayes_filename, "w")
        if (isOpen(bayes_file)) {
          close(bayes_file)
        }
        file.remove(bayes_filename)
      }

      if (scenarios_config$parallelization == FALSE) {
        # Loop through cases
        for (i in 1:nrow(cases)) {
          scenario <- cases[i, ]
          results <- bayesian_ocs_scenario_simulation(
            scenario = scenario,
            simulation_config = simulation_config,
            bayes_filename = bayes_filename,
            config_dir = config_dir
          )
          progress(i)
        }
        # At the end of the loop, close file if open
        bayes_file <- file(bayes_filename, "a")
        if (isOpen(bayes_file)) {
          close(bayes_file)
        }
      } else {
        ncores <- min(parallel::detectCores() - 1, 8) # Leave on core free

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
            "simulation_config",
            "bayes_filename",
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
          scenario <- cases[i, ]
          bayesian_ocs_scenario_simulation(scenario,
                                           simulation_config,
                                           bayes_filename,
                                           config_dir)
        }
        write_csv(results, bayes_filename)
        parallel::stopCluster(cl)
      }
    }
  }
}
