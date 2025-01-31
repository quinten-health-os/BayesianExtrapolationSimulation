frequentist_col_types <- cols(
  method = col_character(),
  parameters = col_character(),
  target_sample_size_per_arm = col_double(),
  drift = col_double(),
  treatment_drift = col_double(),
  control_drift = col_double(),
  source_denominator = col_double(),
  source_denominator_change_factor = col_double(),
  case_study = col_character(),
  parallelization = col_logical(),
  sampling_approximation = col_logical(),
  source_treatment_effect_estimate = col_double(),
  target_treatment_effect = col_double(),
  target_to_source_std_ratio = col_double(),
  theta_0 = col_double(),
  null_space = col_character(),
  summary_measure_likelihood = col_character(),
  target_standard_deviation = col_double(),
  target_treatment_rate = col_double(),
  target_control_rate = col_double(),
  source_standard_error = col_double(),
  source_sample_size_control = col_double(),
  source_sample_size_treatment = col_double(),
  equivalent_source_sample_size_per_arm = col_double(),
  endpoint = col_character(),
  source_control_rate = col_double(),
  source_treatment_rate = col_double(),
  rng_state = col_character(),
  computation_time = col_double(),
  success_proba = col_double(),
  mcse_success_proba = col_double(),
  conf_int_success_proba_lower = col_double(),
  conf_int_success_proba_upper = col_double(),
  coverage = col_double(),
  conf_int_coverage_lower = col_double(),
  conf_int_coverage_upper = col_double(),
  mse = col_double(),
  conf_int_mse_lower = col_double(),
  conf_int_mse_upper = col_double(),
  bias = col_double(),
  conf_int_bias_lower = col_double(),
  conf_int_bias_upper = col_double(),
  posterior_mean = col_double(),
  conf_int_posterior_mean_lower = col_double(),
  conf_int_posterior_mean_upper = col_double(),
  posterior_median = col_double(),
  conf_int_posterior_median_lower = col_double(),
  conf_int_posterior_median_upper = col_double(),
  precision = col_double(),
  conf_int_precision_lower = col_double(),
  conf_int_precision_upper = col_double(),
  credible_interval_lower = col_double(),
  credible_interval_upper = col_double(),
  posterior_parameters = col_character(),
  ess_moment = col_double(),
  conf_int_ess_moment_lower = col_double(),
  conf_int_ess_moment_upper = col_double(),
  ess_precision = col_double(),
  conf_int_ess_precision_lower = col_double(),
  conf_int_ess_precision_upper = col_double(),
  ess_elir = col_double(),
  conf_int_ess_elir_lower = col_double(),
  conf_int_ess_elir_upper = col_double(),
  rhat = col_double(),
  conf_int_rhat_lower = col_double(),
  conf_int_rhat_upper = col_double(),
  mcmc_ess = col_double(),
  conf_int_mcmc_ess_lower = col_double(),
  conf_int_mcmc_ess_upper = col_double(),
  n_divergences = col_double(),
  conf_int_n_divergences_lower = col_double(),
  conf_int_n_divergences_upper = col_double(),
  warning = col_logical(),
  tie = col_double(),
  mcse_tie = col_double(),
  conf_int_tie_lower = col_double(),
  conf_int_tie_upper = col_double(),
  frequentist_power_at_equivalent_tie = col_double(),
  frequentist_power_at_equivalent_tie_lower = col_double(),
  frequentist_power_at_equivalent_tie_upper = col_double(),
  frequentist_test = col_character(),
  nominal_frequentist_power_separate = col_double(),
  nominal_frequentist_power_pooling = col_double()
)


bayesian_col_types <- cols(
  case_study = col_character(),
  method = col_character(),
  target_sample_size_per_arm = col_double(),
  parameters = col_character(),
  source_denominator_change_factor = col_double(),
  target_to_source_std_ratio = col_double(),
  design_prior_type = col_character(),
  prior_proba_success = col_double(),
  prior_proba_no_benefit = col_double(),
  prior_proba_benefit = col_double(),
  prepost_proba_FP = col_double(),
  prepost_proba_TP = col_double(),
  average_tie = col_double(),
  average_power = col_double(),
  upper_bound_proba_FP = col_double(),
  source_treatment_effect_estimate = col_double(),
  source_standard_error = col_double(),
  endpoint = col_character(),
  summary_measure_likelihood = col_character(),
  source_sample_size_control = col_double(),
  source_sample_size_treatment = col_double(),
  equivalent_source_sample_size_per_arm = col_double()
)


sweet_spot_col_types <- cols(
  method = col_character(),
  parameters = col_character(),
  target_sample_size_per_arm = col_double(),
  control_drift = col_double(),
  source_denominator = col_double(),
  source_denominator_change_factor = col_double(),
  case_study = col_character(),
  parallelization = col_logical(),
  sampling_approximation = col_logical(),
  source_treatment_effect_estimate = col_double(),
  target_to_source_std_ratio = col_double(),
  theta_0 = col_double(),
  null_space = col_character(),
  summary_measure_likelihood = col_character(),
  source_standard_error = col_double(),
  source_sample_size_control = col_double(),
  source_sample_size_treatment = col_double(),
  equivalent_source_sample_size_per_arm = col_double(),
  endpoint = col_character(),
  source_control_rate = col_double(),
  source_treatment_rate = col_double(),
  sweet_spot_lower = col_double(),
  sweet_spot_upper = col_double(),
  sweet_spot_width = col_double(),
  drift_range_upper = col_double(),
  drift_range_lower = col_double(),
  metric = col_character()
)


scenario_columns <- c("target_sample_size_per_arm", "control_drift", "source_denominator", "source_denominator_change_factor",
                      "case_study", "sampling_approximation", "source_treatment_effect_estimate",
                      "target_to_source_std_ratio", "theta_0", "null_space",
                      "summary_measure_likelihood", "source_standard_error", "source_sample_size_control",
                      "source_sample_size_treatment", "equivalent_source_sample_size_per_arm", "endpoint",
                      "source_control_rate", "source_treatment_rate")

unique_scenario_columns <- c("case_study", "target_sample_size_per_arm", "source_denominator_change_factor",
                      "target_to_source_std_ratio")


results_columns <- c("success_proba", "mcse_success_proba", "conf_int_success_proba_lower",
                     "conf_int_success_proba_upper", "coverage", "conf_int_coverage_lower",
                     "conf_int_coverage_upper", "mse", "conf_int_mse_lower", "conf_int_mse_upper",
                     "bias", "conf_int_bias_lower", "conf_int_bias_upper", "posterior_mean",
                     "conf_int_posterior_mean_lower", "conf_int_posterior_mean_upper", "posterior_median",
                     "conf_int_posterior_median_lower", "conf_int_posterior_median_upper", "precision",
                     "conf_int_precision_lower", "conf_int_precision_upper", "credible_interval_lower",
                     "credible_interval_upper", "posterior_parameters", "ess_moment",
                     "conf_int_ess_moment_lower", "conf_int_ess_moment_upper", "ess_precision",
                     "conf_int_ess_precision_lower", "conf_int_ess_precision_upper", "ess_elir",
                     "conf_int_ess_elir_lower", "conf_int_ess_elir_upper", "rhat", "conf_int_rhat_lower",
                     "conf_int_rhat_upper", "mcmc_ess", "conf_int_mcmc_ess_lower", "conf_int_mcmc_ess_upper",
                     "n_divergences", "conf_int_n_divergences_lower", "conf_int_n_divergences_upper",
                     "warning", "tie", "mcse_tie", "conf_int_tie_lower", "conf_int_tie_upper",
                     "frequentist_power_at_equivalent_tie", "frequentist_power_at_equivalent_tie_lower",
                     "frequentist_power_at_equivalent_tie_upper", "frequentist_test",
                     "nominal_frequentist_power_separate", "nominal_frequentist_power_pooling")


expected_colnames_scenario <- c(
  "method",
  "parameters",
  "target_sample_size_per_arm",
  "drift",
  "treatment_drift",
  "control_drift",
  "source_denominator",
  "source_denominator_change_factor",
  "case_study",
  "parallelization",
  "sampling_approximation",
  "source_treatment_effect_estimate",
  "target_treatment_effect",
  "target_to_source_std_ratio",
  "theta_0",
  "null_space"
)

expected_colnames_results <- c(
  "success_proba",
  "mcse_success_proba",
  "conf_int_success_proba_lower",
  "conf_int_success_proba_upper",
  "coverage",
  "conf_int_coverage_lower",
  "conf_int_coverage_upper",
  "mse",
  "conf_int_mse_lower",
  "conf_int_mse_upper",
  "bias",
  "conf_int_bias_lower",
  "conf_int_bias_upper",
  "posterior_mean",
  "conf_int_posterior_mean_lower",
  "conf_int_posterior_mean_upper",
  "posterior_median",
  "conf_int_posterior_median_lower",
  "conf_int_posterior_median_upper",
  "precision",
  "conf_int_precision_lower",
  "conf_int_precision_upper",
  "credible_interval_lower",
  "credible_interval_upper",
  "posterior_parameters",
  "ess_moment",
  "conf_int_ess_moment_lower",
  "conf_int_ess_moment_upper",
  "ess_precision",
  "conf_int_ess_precision_lower",
  "conf_int_ess_precision_upper",
  "ess_elir",
  "conf_int_ess_elir_lower",
  "conf_int_ess_elir_upper",
  "rhat",
  "conf_int_rhat_lower",
  "conf_int_rhat_upper",
  "mcmc_ess",
  "conf_int_mcmc_ess_lower",
  "conf_int_mcmc_ess_upper",
  "n_divergences",
  "conf_int_n_divergences_lower",
  "conf_int_n_divergences_upper",
  "warning"
)

expected_colnames_source <- c(
  "source_treatment_effect_estimate",
  "source_standard_error",
  "source_sample_size_control",
  "source_sample_size_treatment",
  "equivalent_source_sample_size_per_arm",
  "summary_measure_likelihood",
  "endpoint",
  "source_control_rate",
  "source_treatment_rate"
)


#' Check columns in a dataframe
#'
#' This function checks if a dataframe contains the expected columns and only the expected columns.
#'
#' @param df A dataframe to check.
#' @param expected_colnames A character vector of column names the dataframe should contain.
#'
#' @return No return value, called for side effects.
#'
#' @keywords internal
check_colnames <- function(df, expected_colnames) {
  missing_cols <- setdiff(expected_colnames, colnames(df))
  if (length(missing_cols) > 0) {
    error_message <- paste(
      "The dataframe is missing the following columns:",
      paste(missing_cols, collapse = ", ")
    )
    stop(error_message)
  }

  missing_cols <- setdiff(colnames(df), expected_colnames)
  if (length(missing_cols) > 0) {
    error_message <- paste(
      "The dataframe has unexpected columns:",
      paste(missing_cols, collapse = ", ")
    )
    stop(error_message)
  }
}

#' Remove columns from a dataframe
#'
#' This function removes specified columns from a dataframe.
#'
#' @param df A dataframe to modify.
#' @param to_remove A character vector of column names to remove.
#'
#' @return A dataframe with specified columns removed.
#'
#' @keywords internal
remove_columns_from_df <- function(df, to_remove) {
  for (key in to_remove) {
    if (key %in% colnames(df)) {
      df <- df[, -which(colnames(df) == key)]
    }
  }
  return(df)
}

#' Concatenate all simulation results to a global dataframe
#'
#' This function reads all CSV files with names starting with "results" from the specified
#' directory, concatenates them, and saves the result to a new CSV file.
#'
#' @param results_dir A character string specifying the path of the results folder.
#' @param ocs_filename A character string specifying the filename of the results dataframe.
#'
#' @return No return value, called for side effects.
#'
#' @importFrom dplyr bind_rows
#'
#' @keywords internal
concatenate_simulation_results <- function(results_dir, ocs_filename) {
  # Initialize an empty data frame to store the concatenated results
  df_results <- data.frame()

  if (grepl("frequentist", ocs_filename)) {
    ocs_type <- "frequentist"
  } else if (grepl("bayesian", ocs_filename)) {
    ocs_type <- "bayesian"
  } else {
    stop("Wrong filename defined in output config.")
  }
  folder_dir <- file.path(results_dir, ocs_type)

  # List all files in the directory tree
  files <- list.files(folder_dir, recursive = TRUE, full.names = TRUE)

  for (file_path in files) {
    if (grepl("\\.csv$", file_path)) {
      if (startsWith(basename(file_path), "results")) {
        # Read the CSV file
        df <- read.csv(file_path, stringsAsFactors = FALSE)

        # Concatenate to the global data frame
        df_results <- dplyr::bind_rows(df_results, df)
      }
    }
  }

  # Save the global data frame to a new CSV file
  readr::write_csv(df_results, file.path(results_dir, ocs_filename))
}


# Function to concatenate files
concatenate_files <- function(file_name, folders, output_path) {
  if (file_name == "results_frequentist.csv"){
    col_types <- frequentist_col_types
  } else if (file_name == "results_bayesian_simpson.csv"){
    col_types <- bayesian_col_types
  } else if (file_name == "sweet_spot.csv"){
    col_types <- sweet_spot_col_types
  } else {
    col_types <- NULL
  }

  combined_data <- bind_rows(lapply(folders, function(folder) {
    file_path <- file.path(folder, file_name)
    if (file.exists(file_path)) {
      data <- read_csv(file_path, col_types = col_types, col_names = TRUE)
      return(data)
    } else {
      NULL  # Handle cases where the file might not exist in some folders
    }
  }))


  if (!(setequal(colnames(combined_data), names(col_types$cols)))){
    warning("Invalid column names in the concatenated results. Concatenated results file not created.")
  } else {
    # Write the combined data to the root folder
    write_csv(combined_data, file.path(output_path, file_name))
  }
}


#' Concatenate all simulation logs to a global dataframe
#'
#' This function reads all CSV files with names starting with "logs" from both the "bayesian" and
#' "frequentist" subdirectories, concatenates them, adds an "OCs" column, and saves the result to a new CSV file.
#'
#' @param results_dir A character string specifying the path of the results folder.
#'
#' @return No return value, called for side effects.
#'
#' @importFrom dplyr bind_rows select everything
#'
#' @export
concatenate_simulation_logs <- function(results_dir) {
  # Initialize an empty data frame to store the concatenated results
  df_logs <- data.frame()

  ocs_types <- c("bayesian", "frequentist")

  for (ocs_type in ocs_types) {
    folder_dir <- file.path(results_dir, ocs_type)

    # List all files in the directory tree
    files <- list.files(folder_dir, recursive = TRUE, full.names = TRUE)

    for (file_path in files) {
      if (grepl("\\.csv$", file_path)) {
        if (startsWith(basename(file_path), "logs")) {
          # Read the CSV file
          df <- read.csv(file_path, stringsAsFactors = FALSE)
          # Add OCs column
          df$OCs <- ocs_type
          # Concatenate to the global data frame
          df_logs <- dplyr::bind_rows(df_logs, df)
        }
      }
    }
  }
  # Reorder columns to place OCs column at the start
  df_logs <- df_logs %>%
    dplyr::select(OCs, everything())
  # Save the global data frame to a new CSV file
  readr::write_csv(df_logs, file.path(results_dir, "logs.csv"))
}


format_simulation_output_table <- function(output) {
  # Extract the data
  data <- list(
    "Success Probability" = c(
      output$success_proba,
      output$conf_int_success_proba_lower,
      output$conf_int_success_proba_upper
    ),
    "Coverage" = c(
      output$coverage,
      output$conf_int_coverage_lower,
      output$conf_int_coverage_upper
    ),
    "MSE" = c(
      output$mse,
      output$conf_int_mse_lower,
      output$conf_int_mse_upper
    ),
    "Bias" = c(
      output$bias,
      output$conf_int_bias_lower,
      output$conf_int_bias_upper
    ),
    "Posterior Mean" = c(
      output$posterior_mean,
      output$conf_int_posterior_mean_lower,
      output$conf_int_posterior_mean_upper
    ),
    "Posterior Median" = c(
      output$posterior_median,
      output$conf_int_posterior_median_lower,
      output$conf_int_posterior_median_upper
    ),
    "Precision" = c(
      output$precision,
      output$conf_int_precision_lower,
      output$conf_int_precision_upper
    ),
    "Credible Interval" = c(
      NA,
      output$credible_interval_lower,
      output$credible_interval_upper
    ),
    "ESS Moment" = c(
      output$ess_moment,
      output$conf_int_ess_moment_lower,
      output$conf_int_ess_moment_upper
    ),
    "ESS Precision" = c(
      output$ess_precision,
      output$conf_int_ess_precision_lower,
      output$conf_int_ess_precision_upper
    ),
    "ESS ELIR" = c(
      output$ess_elir,
      output$conf_int_ess_elir_lower,
      output$conf_int_ess_elir_upper
    )
  )

  # Create a data frame
  df <- data.frame(
    Metric = names(data),
    Value = sapply(data, `[`, 1),
    CI_Lower = sapply(data, `[`, 2),
    CI_Upper = sapply(data, `[`, 3)
  )

  # Format the values
  df$Value <- format(df$Value, digits = 4, scientific = FALSE)
  df$CI_Lower <- format(df$CI_Lower, digits = 4, scientific = FALSE)
  df$CI_Upper <- format(df$CI_Upper, digits = 4, scientific = FALSE)

  # Create the confidence interval column
  df$`95% CI` <- paste0("[", df$CI_Lower, ", ", df$CI_Upper, "]")

  # Select and rename columns
  df <- df[, c("Metric", "Value", "95% CI")]

  # Print the table
  print(df, row.names = FALSE)
}

format_parameters_to_json <- function(json_parameters, escape = FALSE) {
  json_parameters <- jsonlite::toJSON(json_parameters, pretty = TRUE)
  json_parameters <- gsub("\"", "\'", json_parameters)
  json_parameters <- as.character(json_parameters)
  if (escape == TRUE) {
    json_parameters <- paste0("[\n  ", json_parameters, "\n]")
  }
  return(json_parameters)
}


format_case_study_config <- function(case_study_config) {
  if (case_study_config$endpoint == "continuous" |
      case_study_config$endpoint == "recurrent_event") {
    data <- data.frame(
      Parameter = c(
        "Name",
        "Control",
        "Summary Measure Likelihood",
        "Theta 0",
        "Endpoint",
        "Null Space",
        "Sampling Approximation",
        "Control",
        "Treatment",
        "Total",
        "Treatment Effect",
        "Standard Error",
        "Control",
        "Treatment",
        "Total",
        "Treatment Effect",
        "Standard Error"
      ),
      Value = c(
        case_study_config$name,
        case_study_config$control,
        case_study_config$summary_measure_likelihood,
        case_study_config$theta_0,
        case_study_config$endpoint,
        case_study_config$null_space,
        case_study_config$sampling_approximation,
        case_study_config$target$control,
        case_study_config$target$treatment,
        case_study_config$target$total,
        case_study_config$target$treatment_effect,
        case_study_config$target$standard_error,
        case_study_config$source$control,
        case_study_config$source$treatment,
        case_study_config$source$total,
        case_study_config$source$treatment_effect,
        case_study_config$source$standard_error
      )
    )

    # Create the table with headers and horizontal lines
    kable(data,
          col.names = c("Parameter", "Value"),
          align = "l") %>%
      kableExtra::add_header_above(c("Case Study Configuration" = 2)) %>%
      kableExtra::kable_styling("striped", full_width = F) %>%
      kableExtra::pack_rows("General", 1, 7) %>%
      kableExtra::pack_rows("Target", 8, 12) %>%
      kableExtra::pack_rows("Source", 13, 17)
  } else if (case_study_config$endpoint == "time_to_event") {
    data <- data.frame(
      Parameter = c(
        "Name",
        "Summary Measure Likelihood",
        "Theta 0",
        "Endpoint",
        "Null Space",
        "Sampling Approximation",
        "Control (Target)",
        "Treatment (Target)",
        "Total (Target)",
        "Treatment Effect (Target)",
        "Standard Error (Target)",
        "Maximum follow-up time (Target)",
        "Control (Source)",
        "Treatment (Source)",
        "Total (Source)",
        "Treatment Effect (Source)",
        "Standard Error (Source)",
        "Maximum follow-up time (Source)"
      ),
      Value = c(
        case_study_config$name,
        case_study_config$summary_measure_likelihood,
        case_study_config$theta_0,
        case_study_config$endpoint,
        case_study_config$null_space,
        case_study_config$sampling_approximation,
        case_study_config$target$control,
        case_study_config$target$treatment,
        case_study_config$target$total,
        case_study_config$target$treatment_effect,
        case_study_config$target$standard_error,
        case_study_config$target$max_follow_up_time,
        case_study_config$source$control,
        case_study_config$source$treatment,
        case_study_config$source$total,
        case_study_config$source$treatment_effect,
        case_study_config$source$standard_error,
        case_study_config$source$max_follow_up_time
      )
    )

    # Create the table with headers and horizontal lines
    kable(data,
          col.names = c("Parameter", "Value"),
          align = "l") %>%
      kableExtra::add_header_above(c("Case Study Configuration" = 2)) %>%
      kableExtra::kable_styling("striped", full_width = F) %>%
      kableExtra::pack_rows("General", 1, 6) %>%
      kableExtra::pack_rows("Target", 7, 12) %>%
      kableExtra::pack_rows("Source", 13, 18)
  } else if (case_study_config$endpoint == "binary") {
    data <- data.frame(
      Parameter = c(
        "Name",
        "Control",
        "Summary Measure Likelihood",
        "Theta 0",
        "Endpoint",
        "Null Space",
        "Sampling Approximation",
        "Target Control",
        "Target Treatment",
        "Target Total",
        "Target Responses (Control)",
        "Target Responses (Treatment)",
        "Target Treatment Effect",
        "Target Standard Error",
        "Source Control",
        "Source Treatment",
        "Source Total",
        "Source Responses (Control)",
        "Source Responses (Treatment)",
        "Source Treatment Effect",
        "Source Standard Error"
      ),
      Value = c(
        case_study_config$name,
        case_study_config$control,
        case_study_config$summary_measure_likelihood,
        case_study_config$theta_0,
        case_study_config$endpoint,
        case_study_config$null_space,
        case_study_config$sampling_approximation,
        case_study_config$target$control,
        case_study_config$target$treatment,
        case_study_config$target$total,
        case_study_config$target$responses$control,
        case_study_config$target$responses$treatment,
        case_study_config$target$treatment_effect,
        case_study_config$target$standard_error,
        case_study_config$source$control,
        case_study_config$source$treatment,
        case_study_config$source$total,
        case_study_config$source$responses$control,
        case_study_config$source$responses$treatment,
        case_study_config$source$treatment_effect,
        case_study_config$source$standard_error
      )
    )

    # Create the table with headers and horizontal lines
    kable(data,
          col.names = c("Parameter", "Value"),
          align = "l") %>%
      kableExtra::add_header_above(c("Case Study Configuration" = 2)) %>%
      kableExtra::kable_styling("striped", full_width = F) %>%
      kableExtra::pack_rows("General", 1, 7) %>%
      kableExtra::pack_rows("Target", 8, 14) %>%
      kableExtra::pack_rows("Source", 15, 21)
  }
}


compile_stan_model <- function(model_name, stan_model_code) {
  stan_directory <- paste0(system.file("stan", package = "RBExT"), "/")
  stan_model_file_path <- paste0(stan_directory, model_name, ".stan")
  writeLines(stan_model_code, con = stan_model_file_path)
  stan_exe_file_path <- paste0(stan_directory, model_name, ".exe")

  cpp_options <- list(stan_threads = TRUE)


  if (!file.exists(stan_exe_file_path)) {
    stan_model <- cmdstanr::cmdstan_model(stan_model_file_path,
                                          exe_file = stan_exe_file_path,
                                          cpp_options = cpp_options)
  } else {
    stan_model <- cmdstanr::cmdstan_model(exe_file = stan_exe_file_path, cpp_options = cpp_options)
  }
  return(stan_model)
}

load_data <- function(results_row, type, reload_data_objects = FALSE) {
  if (!(type %in% c("target", "source"))){
    stop("type must either be 'target' or 'source'.")
  }

  if (reload_data_objects) {
    case_study_config <- yaml::yaml.load_file(system.file(
      paste0("conf/case_studies/", results_row$case_study, ".yml"),
      package = "RBExT"
    ))

    source_data <- SourceData$new(
      case_study_config = case_study_config,
      source_denominator = results_row$source_denominator
    )
    if (type == "target") {
      data <- TargetDataFactory$new()
      data <- data$create(
        source_data = source_data,
        case_study_config = case_study_config,
        target_sample_size_per_arm = results_row$target_sample_size_per_arm,
        treatment_drift = results_row$treatment_drift,
        control_drift = results_row$control_drift,
        summary_measure_likelihood = source_data$summary_measure_likelihood,
        target_to_source_std_ratio = results_row$target_to_source_std_ratio
      )
    } else if (type == "source") {
      data <- source_data
    }
  } else {
    if (type == "target") {
      data <- list(
        sample_size_per_arm = results_row$target_sample_size_per_arm,
        sample_size_control = results_row$target_sample_size_per_arm,
        sample_size_treatment = results_row$target_sample_size_per_arm,
        treatment_effect = results_row$target_treatment_effect,
        standard_deviation = results_row$target_standard_deviation,
        summary_measure_likelihood = results_row$summary_measure_likelihood,
        treatment_rate = results_row$target_treatment_rate,
        control_rate = results_row$target_control_rate
      )
    } else if (type == "source") {
      data <- list(
        treatment_effect_estimate = results_row$source_treatment_effect_estimate,
        standard_error = results_row$source_standard_error,
        sample_size_control = results_row$source_sample_size_control,
        sample_size_treatment = results_row$source_sample_size_treatment,
        equivalent_source_sample_size_per_arm = results_row$equivalent_source_sample_size_per_arm,
        treatment_rate = results_row$source_treatment_rate,
        control_rate = results_row$source_control_rate
      )

      if (is.null(data$equivalent_source_sample_size_per_arm) |
          any(is.na(data$equivalent_source_sample_size_per_arm))) {
        data$equivalent_source_sample_size_per_arm <- 2 * data$sample_size_control * data$sample_size_treatment / (data$sample_size_control + source_data$sample_size_treatment)
      }
    }
  }
  return(data)
}

# Function to create boolean filter
create_boolean_filter <- function(df, conditions, exclude_key = NULL) {
  if (!is.null(exclude_key)) {
    conditions <- conditions[!names(conditions) %in% exclude_key]
  }

  # Convert all columns to numeric where possible
  df <- data.frame(lapply(df, function(x)
    as.numeric(x)))


  filter <- rep(TRUE, nrow(df))
  for (col in names(conditions)) {
    condition <- conditions[[col]]
    filter <- (df[[col]] == condition) & filter
  }
  return(filter)
}


check_confidence_intervals <- function(df, metrics) {
  # Loop through each metric and check the confidence intervals
  for (metric in metrics) {
    upper_col <- paste0("conf_int_", metric, "_upper")
    lower_col <- paste0("conf_int_", metric, "_lower")

    # Check if the metric and its CI bounds exist in the dataframe
    if (all(c(metric, lower_col, upper_col) %in% names(df))) {

      # Check if the lower bound is indeed lower than the upper bound
      invalid_bounds <- which(df[[lower_col]] > df[[upper_col]])

      if (length(invalid_bounds) > 0) {
        issue <- paste0("Lower bound of the CI is larger than the upper bound for ", metric)
        warning(issue)
      }

      # Check if the metric value lies within the confidence interval
      out_of_bounds <- which(df[[metric]] < df[[lower_col]] | df[[metric]] > df[[upper_col]])
      if (length(out_of_bounds) > 0) {
        issue <- paste0("Mean value outside the CI bounds for ", metric)
        warning(issue)
      }
    }
  }
}


# Used to extract parameters as nested list
extract_nested_parameter = function(parameters){
  # Initialize an empty list to store the nested list
  nested_list <- list()

  # Iterate over the column names of the parameters to construct the nested list
  for (colname in colnames(parameters)) {

    # Split the column name by period (".") to find the hierarchy
    split_names <- strsplit(colname, "\\.")[[1]]

    if (length(split_names) == 2) {
      if (!split_names[1] %in% names(nested_list)) {
        nested_list[[split_names[1]]] <- list()
      }

      x <- parameters[[colname]]
      x_numeric <- suppressWarnings(as.numeric(x))

      # Replace NA values with the original input
      x_clean <- ifelse(is.na(x_numeric), x, x_numeric)

      # Assign values dynamically to the second level (e.g., family or std_dev)
      nested_list[[split_names[1]]][[split_names[2]]] <- x_clean

      # For 'initial_prior', assign the value directly
    } else {
      nested_list[colname] <- parameters[[colname]]
    }
  }
  return(nested_list)
}

# Function to log errors globally
global_error_handler <- function() {
  err <- geterrmessage()  # Get the error message
  futile.logger::flog.error("Global error occurred: %s", err)

  # Optionally log other debugging info such as the call stack
  futile.logger::flog.error("Call stack:\n%s", paste(deparse(sys.calls()), collapse = "\n"))

  # Log additional variables of interest (if needed)
  # For instance, if you're in a loop or a function with certain variables
  # futile.logger::flog.error("Variable state at error - x: %s, y: %s", x, y)  # Customize as needed
}

generate_log_filename <- function(base_name = "error_log.log", suffix_type = "timestamp") {
  if (!file.exists(base_name)) {
    return(base_name)  # Return the base name if no file exists
  }

  # If the file exists, create a new filename with a suffix
  if (suffix_type == "timestamp") {
    # Add a timestamp to the filename
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    new_name <- sub(".log$", paste0("_", timestamp, ".log"), base_name)
  } else if (suffix_type == "counter") {
    # Use a counter to create unique filenames
    counter <- 1
    repeat {
      new_name <- sub(".log$", paste0("_", counter, ".log"), base_name)
      if (!file.exists(new_name)) break
      counter <- counter + 1
    }
  }

  return(new_name)
}

save_state <- function(iteration, scenario, worker_id = NULL, env) {
  if (is.null(worker_id)){
    file = paste0("./logs/", env, "/checkpoints/checkpoint_iter_", iteration, ".RData")
  } else {
    file = paste0("./logs/", env, "/checkpoints/checkpoint_iter_", worker_id, "_iter_", iteration, ".RData")
  }
  # Extract the directory path from the file path
  dir_path <- dirname(file)

  # Check if the directory exists, and if not, create it
  if (!file.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }

  save(scenario, file = file)
}

read_function_code <- function(function_obj){
  function_name <- substitute(function_obj)

  # Use deparse() to get the function's source code
  function_code <- deparse(function_obj)

  # Prepend the function name and assignment
  cat(paste0(function_name, " <- ", paste(function_code, collapse = "\n"), sep = ""))
}



check_simulation_completeness <- function(results_dir = results_dir, ocs_filename = ocs_filename, scenarios_config, config_dir) {
  for (case_study in scenarios_config$case_studies){
    case_study_config <- yaml::yaml.load_file(system.file(
      paste0("conf/case_studies/", case_study, ".yml"),
      package = "RBExT"
    ))
    for (method in scenarios_config$method){
      file_path <- paste0(results_dir, 'frequentist/', case_study, '/', method, '/', ocs_filename)

      if (file.exists(file_path)) {
        results_df <- readr::read_csv(file_path, col_types = frequentist_col_types, col_names = TRUE)
      } else {
        # File doesn't exist, continue with the rest of the code
        message("File does not exist, continuing...")
        next
      }

      null_space <- case_study_config$null_space
      theta_0 <- case_study_config$theta_0

      source_denominator_change_factors <- unique(results_df$source_denominator_change_factor)

      if (!setequal(scenarios_config$denominator_change_factor, source_denominator_change_factors) && !(case_study %in% c("botox", "dapagliflozin", "aprepitant"))){
        warning(paste0("Simulaton is incomplete, all denominator change factors are not included. Method: ", method, ", Case study : ", case_study))
      }

      for (source_denominator_change_factor in source_denominator_change_factors) {
        results_df_1 <- results_df %>%
          dplyr::filter(
            source_denominator_change_factor == !!source_denominator_change_factor  |
              is.na(source_denominator_change_factor)
          )
        target_to_source_std_ratio_range <- unique(results_df_1$target_to_source_std_ratio)

        if (!setequal(scenarios_config$target_to_source_std_ratio_range, target_to_source_std_ratio_range) && (case_study %in% c("botox", "dapagliflozin"))){
          warning(paste0("Simulaton is incomplete, all target to source std ratios are not included. Method: ", method, ", Case study : ", case_study))
        }

        for (target_to_source_std_ratio in target_to_source_std_ratio_range) {
          results_df_2 <- results_df_1 %>%
            dplyr::filter(target_to_source_std_ratio == !!target_to_source_std_ratio |
                            is.na(target_to_source_std_ratio))


          if (nrow(results_df_2) == 0) {
            stop("Dataframe is empty")
          }

          theta_0 <- case_study_config$theta_0

          target_sample_sizes <- unique(results_df_2$target_sample_size_per_arm)

          if (!(length(scenarios_config$sample_size_factors) == length(target_sample_sizes))){
            warning(paste0("Simulaton is incomplete, all sample size factors are not included. Method: ", method, ", Case study : ", case_study))
          }

          for (target_sample_size_per_arm in target_sample_sizes) {
            results_df_3 <- results_df_2 %>%
              dplyr::filter(target_sample_size_per_arm == !!target_sample_size_per_arm)

            mandatory_drift_values <- important_drift_values(unique(results_df_3$source_treatment_effect_estimate), case_study_config)

            # Make sure the mandatory drift values are in the drift values
            drift_values <- unique(results_df_3$drift)

            if (sum(!(mandatory_drift_values %in% drift_values))>1){
              stop(paste0("Simulaton is incomplete, some mandatory drift values are not included. Method: ", method, ", Case study : ", case_study))
            }

            drift_values <- setdiff(drift_values, mandatory_drift_values)
            if (length(drift_values) < (scenarios_config$ndrift - length(mandatory_drift_values))){
              stop(paste0("Simulaton is incomplete, all drift values are not included. Method: ", method, ", Case study : ", case_study))
            }
          }
        }
      }
    }
  }
}


is_approx_equal <- function(x, y, tolerance = .Machine$double.eps^0.5) {
  abs(x - y) < tolerance
}

# Function to compare rows, ignoring NAs
compare_ignore_na <- function(row, ref_row) {
  # Only compare non-NA elements in the reference row
  non_na_indices <- !is.na(ref_row)
  all(row[non_na_indices] == ref_row[non_na_indices])
}
