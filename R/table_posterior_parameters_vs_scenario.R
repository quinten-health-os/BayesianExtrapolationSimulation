#' Generate a table of posterior parameters vs drift
#'
#' This function creates tables comparing posterior parameters to drift for a given method, case study, and sample size.
#'
#' @param results_metrics_df A data frame containing the results metrics.
#' @param method The method to filter the results by.
#' @param case_study The case study to filter the results by.
#' @param target_sample_size_per_arm The target sample size per arm to filter the results by.
#' @param control_drift A logical value indicating whether to control for drift.
#' @param xvars A list containing the x-variable configurations.
#'
#' @return This function does not return a value. It creates and saves tables in HTML, PDF, and LaTeX formats.
#'
#' @import dplyr
#' @import knitr
#' @export
table_posterior_parameters_vs_drift <- function(results_metrics_df,
                                                method,
                                                case_study,
                                                target_sample_size_per_arm,
                                                control_drift,
                                                xvars,
                                                source_denominator_change_factor,
                                                target_to_source_std_ratio) {
  # Filter the data frame based on method, case_study, and target_sample_size_per_arm
  results_metrics_df <- results_metrics_df %>%
    dplyr::filter(
      method == !!method,
      case_study == !!case_study,
      target_sample_size_per_arm == !!target_sample_size_per_arm,
      source_denominator_change_factor == !!source_denominator_change_factor  |
        is.na(source_denominator_change_factor),
      target_to_source_std_ratio == !!target_to_source_std_ratio |
        is.na(target_to_source_std_ratio)
    )

  if (control_drift) {
    xvar <- xvars[["control_drift"]]
    xvar$name <- "control_drift"
    # Remove non-zero treatment drift
    results_metrics_df <- results_metrics_df %>%
      dplyr::filter(treatment_drift == 0)
  } else {
    xvar <- xvars[["drift"]]
    results_metrics_df <- results_metrics_df %>%
      dplyr::filter(control_drift == 0)
  }

  source_treatment_effect_estimate <- unique(results_metrics_df$source_treatment_effect_estimate)[1]

  case_study_config <- yaml::read_yaml(paste0(case_studies_config_dir, case_study, ".yml"))

  # mandatory_drift_values <- important_drift_values(source_treatment_effect_estimate, case_study_config)

  posterior_parameters_df <- get_parameters(results_metrics_df[, "posterior_parameters"])


  posterior_parameters_df[] <- lapply(posterior_parameters_df, format_num, digits = 4, scientific = TRUE)
  # Identify columns with unique values and drop them
  cols_to_keep <- sapply(posterior_parameters_df, function(col) length(unique(col)) > 1)
  posterior_parameters_df <- posterior_parameters_df[, cols_to_keep, drop = FALSE]
  posterior_parameters_df <- convert_CI_columns(posterior_parameters_df)

  if (nrow(posterior_parameters_df) == 0 || ncol(posterior_parameters_df) == 0){
    return()
  }

  prior_parameters_df <- get_parameters(results_metrics_df[, "parameters"])
  prior_parameters_df[] <- lapply(prior_parameters_df, format_num)
  # Identify columns with unique values and drop them
  cols_to_keep <- sapply(prior_parameters_df, function(col) length(unique(col)) > 1)

  prior_params <- colnames(prior_parameters_df)
  unique_prior_params <- prior_params[sapply(prior_parameters_df, function(col) length(unique(col)) == 1)]
  updated_prior_parameters_df <- prior_parameters_df[, cols_to_keep, drop = FALSE]

  keys <- names(methods_dict[[method]])

  # Function to rename columns using dplyr
  rename_parameters_dplyr <- function(df, prefix) {
    df %>%
      rename_with(
        .cols = any_of(keys),
        .fn = ~ sprintf("%s %s", prefix, sapply(.x, function(key) methods_dict[[method]][[key]]$parameter_notation))
      )
  }

  # Rename columns in prior_parameters_df
  updated_prior_parameters_df <- rename_parameters_dplyr( updated_prior_parameters_df, "Prior")

  keys <- names(empirical_bayes_hyperparameters[[method]])

  print(method)

  rename_posterior_parameters_dplyr <- function(df) {
    df %>%
      rename_with(
        .cols = everything(), # Apply to all columns, or replace `everything()` with specific columns
        .fn = ~ sapply(.x, function(key) {
          parameter_notation <- empirical_bayes_hyperparameters[[method]][[key]][["parameter_notation"]]
          suffix <- empirical_bayes_hyperparameters[[method]][[key]][["suffix"]]
          sprintf("%s %s", parameter_notation, suffix)
        })
      )
  }

  # Rename columns in posterior_parameters_df
  posterior_parameters_df <- rename_posterior_parameters_dplyr(posterior_parameters_df)

  data_table <- cbind(results_metrics_df$drift, updated_prior_parameters_df, posterior_parameters_df)

  colnames(data_table)[1] <- "Drift"

  data_table <- data_table[order(data_table[['Drift']]),]
  # Drop the row names (reset index)
  row.names(data_table) <- NULL

  directory <- file.path(tables_dir, case_study)
  if (!dir.exists(directory)) {
    dir.create(directory,
               showWarnings = FALSE,
               recursive = TRUE)
  }

  filename <- paste0(
    case_study,
    "_",
    method,
    "_prior_vs_posterior_parameters_sample_size_",
    target_sample_size_per_arm
  )

  filename <- format_filename(filename = filename, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)

  file_path <- file.path(directory, filename)

  prior_params = ""
  for (key in unique_prior_params){
    if (key == "initial_prior"){
      next
    }
    prior_params = paste0(prior_params, sprintf("%s %s = %s", "Prior", methods_dict[[method]][[key]]$parameter_notation, unique(prior_parameters_df[key])))
  }


  title <- sprintf(
    "%s, %s, $N_T/2 =$ %s",
    str_to_title(case_study),
    methods_labels[[method]]$label,
    target_sample_size_per_arm
  )


  title <- format_title(title = title, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor, as_latex = TRUE)

  if (method == "commensurate_power_prior"){ # TODO: replace this hard-coded part
    colnames(data_table) <- c(
      "Drift",
      "Heterogeneity prior",
      "Prior $\\alpha$",
      "Prior $\\beta$",
      "Prior $\\sigma_\\tau$",
      "Posterior $\\mu_\\tau$",
      "Posterior $\\sigma_\\tau$",
      "Posterior $\\mu_\\gamma$",
      "Posterior $\\sigma_\\gamma$"
    )
  }

  # Replace all "_" with " " in string values
  data_table[] <- lapply(data_table, function(x) {
    if (is.character(x)) {
      gsub("_", " ", x)
    } else {
      x
    }
  })

  # Capitalize each beginning of words, ignoring NA and "NA" strings
  data_table[] <- lapply(data_table, function(x) {
    if (is.character(x)) {
      sapply(x, function(y) {
        if (!is.na(y) && !grepl("^\\s*NA\\s*$", y)) { # Ignore actual NA and "NA" strings
          str_to_title(trimws(y)) # Capitalize words and remove any leading/trailing spaces
        } else {
          y # Leave NA and "NA" strings as they are
        }
      })
    } else {
      x
    }
  })

  export_table(data_table = data_table, ncollapses = 1, title = title, file_path = file_path)
}


#' Generate multiple posterior parameter tables
#'
#' This function generates multiple tables for posterior parameters across different case studies, methods, and sample sizes.
#'
#' @param results_metrics_df A data frame containing the results metrics.
#' @param metrics A list or vector of metrics to include in the analysis.
#'
#' @return This function does not return a value. It calls other functions to create and save multiple tables.
#'
#' @import dplyr
#' @import yaml
#'
#' @export
posterior_parameters_tables <- function(results_metrics_df, metrics) {
  # Get the list of case studies
  case_studies <- unique(results_metrics_df$case_study)

  # Get the list of methods
  methods <- unique(results_metrics_df$method)

  treatment_effects <- c("partially_consistent", "consistent")

  for (case_study in case_studies) {
    results_metrics_df1 <- results_metrics_df[results_metrics_df$case_study == case_study,]
    for (method in methods) {
      results_metrics_df2 <- results_metrics_df1[results_metrics_df1$method == method,]
      target_sample_sizes <- unique(results_metrics_df2$target_sample_size_per_arm)
      for (target_sample_size_per_arm in target_sample_sizes) {

        results_metrics_df3 <- results_metrics_df2[results_metrics_df2$target_sample_size_per_arm == target_sample_size_per_arm,]

        target_to_source_std_ratios <- unique(results_metrics_df3$target_to_source_std_ratio)

        for (target_to_source_std_ratio in target_to_source_std_ratios) {
          results_metrics_df4 <- results_metrics_df3 %>%
            dplyr::filter(
              target_to_source_std_ratio == !!target_to_source_std_ratio |
                is.na(target_to_source_std_ratio)
            )

          source_denominator_change_factors <- unique(results_metrics_df4$source_denominator_change_factor)

          for (source_denominator_change_factor in source_denominator_change_factors){
            results_metrics_df5 <- results_metrics_df4 %>%
              dplyr::filter(
                source_denominator_change_factor == !!source_denominator_change_factor |
                  is.na(source_denominator_change_factor)
              )

            # Call the posterior_parameters_vs_drift function
            table_posterior_parameters_vs_drift(
              results_metrics_df,
              method,
              case_study,
              target_sample_size_per_arm,
              control_drift = FALSE,
              xvars = xvars,
              source_denominator_change_factor = source_denominator_change_factor,
              target_to_source_std_ratio = target_to_source_std_ratio
            )
          }
        }
      }
    }
  }
}
