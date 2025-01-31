#' Perform Simulation Analysis
#'
#' @description This function performs a simulation analysis for a given environment using the specified
#' analysis configuration and configuration directory. It reads the output configuration,
#' loads results, applies frequentist power analysis, and computes the sweet spot.
#'
#' @param env A character string specifying the environment for the analysis (e.g., "development", "production").
#' @param analysis_config A list containing the analysis configuration parameters.
#' @param config_dir A character string specifying the directory containing the configuration files.
#' @param frequentist_metrics List of frequentist metrics.
#'
#' @return This function does not return a value. It writes the results and sweet spot analysis
#' to CSV files in the specified results directory.
#'
#' @importFrom yaml read_yaml
#' @importFrom readr read_csv write_csv
#'
#' @export
simulation_analysis <- function(env,
                                analysis_config,
                                config_dir,
                                frequentist_metrics, case_studies = "all", to_compute = c("frequentist_power_at_equivalent_tie", "frequentist_power_at_nominal_tie", "sweet_spot", "bayesian_ocs"), methods = "all") {
  futile.logger::flog.info("Starting the analysis of results in environment %s", env)

  results_dir <- paste0("./results/", env)
  outputs_config <- yaml::read_yaml(system.file("conf/outputs_config.yml", package = "RBExT"))
  analysis_config <- yaml::read_yaml(system.file("conf/analysis_config.yml", package = "RBExT"))
  simulation_config <- yaml::read_yaml(system.file("conf/simulation_config.yml", package = "RBExT"))

  freq_filename <- paste0(results_dir,
                          "/",
                          outputs_config$frequentist_ocs_results_filename)

  # Load the results and apply frequentist_power_at_equivalent_tie
  results_freq_df <- readr::read_csv(freq_filename, col_types = frequentist_col_types, col_names = TRUE)

  if (length(case_studies) == 1 && case_studies == "all"){
    case_studies <- unique(results_freq_df$case_study)
  }

  results_freq_subset <- subset(results_freq_df, case_study %in% case_studies)

  if (length(methods) == 1 && methods == "all"){
    methods <- unique(results_freq_subset$method)
  }

  results_freq_subset <- subset(results_freq_subset, method %in% methods)

  if (nrow(results_freq_subset) == 0) {
    warning("Results dataframe is empty.")
    return()
  }

  matching_columns <- c(
    "method",
    "parameters",
    "control_drift",
    "source_denominator",
    "source_denominator_change_factor",
    "case_study",
    "target_to_source_std_ratio",
    "target_sample_size_per_arm",
    "theta_0",
    "null_space",
    "sampling_approximation",
    "summary_measure_likelihood",
    "source_sample_size_treatment",
    "source_sample_size_control",
    "endpoint",
    "source_standard_error",
    "source_treatment_effect_estimate",
    "equivalent_source_sample_size_per_arm"
  )

  if ("frequentist_power_at_equivalent_tie" %in% to_compute){
    new_results_power_df <- frequentist_power_at_equivalent_tie(results = results_freq_subset, analysis_config = analysis_config, simulation_config = simulation_config, parallelization = FALSE)

    # # Filter out the rows in results_freq_df that match new_results_power_df
    # updated_results_df <- results_freq_df %>%
    #   dplyr::filter(!dplyr::if_any(all_of(matching_columns), ~ . %in% new_results_power_df[[.col]]))
    #
    # # Combine the filtered original results with the new results
    # combined_results_df <- dplyr::bind_rows(updated_results_df, new_results_power_df)

    updated_results_df <- results_freq_df %>%
      dplyr::anti_join(new_results_power_df, by = matching_columns)

    # Append new results to existing data
    combined_results_df <- dplyr::bind_rows(updated_results_df, new_results_power_df)

    readr::write_csv(combined_results_df, freq_filename)
  }
  if ("frequentist_power_at_nominal_tie" %in% to_compute){
    new_results_power_df <- frequentist_power_at_nominal_tie(results = results_freq_subset, analysis_config = analysis_config, simulation_config = simulation_config)

    updated_results_df <- results_freq_df %>%
      dplyr::anti_join(new_results_power_df, by = matching_columns)

    # Append new results to existing data
    combined_results_df <- dplyr::bind_rows(updated_results_df, new_results_power_df)

    readr::write_csv(combined_results_df, freq_filename)
  }
  if ("sweet_spot" %in% to_compute){
    sweet_spot_df <- sweet_spot(results_freq_subset, frequentist_metrics)

    jsonlite::write_json(sweet_spot_df, paste0(results_dir, "/sweet_spot.json"))
  }
  if ("bayesian_ocs" %in% to_compute){
    # Compute Bayesian OCs (in a deterministic manner) based on the results
    results_bayesian_ocs <- compute_bayesian_ocs(results_freq_subset, env)
    bayes_filename <- paste0(results_dir,
                             "/",
                             outputs_config$bayesian_ocs_deterministic_results_filename)

    # Append new Bayesian results to existing data
    if (file.exists(bayes_filename)) {
      existing_bayesian_df <- readr::read_csv(bayes_filename, col_names = TRUE)

      # Remove rows in `existing_bayesian_df` that match the keys in `results_bayesian_ocs`
      updated_bayesian_df <- existing_bayesian_df %>%
        dplyr::anti_join(results_bayesian_ocs, by = c("case_study", "method"))

      # Append new Bayesian results
      combined_bayesian_df <- dplyr::bind_rows(updated_bayesian_df, results_bayesian_ocs)

      # Write updated results back to the CSV
      readr::write_csv(combined_bayesian_df, bayes_filename)
    } else {
      # Write new results if file does not exist
      readr::write_csv(results_bayesian_ocs, bayes_filename)
    }
  }
}

