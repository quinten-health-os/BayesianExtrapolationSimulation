#' Generate a table of hyperparameters estimated using Empirical Bayes vs drift
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
table_empirical_bayes_hyperparameters_vs_drift <- function(results_metrics_df,
                                                method,
                                                case_study,
                                                control_drift,
                                                xvars,
                                                source_denominator_change_factor,
                                                target_to_source_std_ratio,
                                                parameters_combinations) {

  if (nrow(parameters_combinations) > 1){
    stop("The function was designed for a single parameters combination.")
  }

  # Filter the data frame based on method, case_study, and target_sample_size_per_arm
  results_metrics_df <- results_metrics_df %>%
    dplyr::filter(
      method == !!method,
      case_study == !!case_study,
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

  results_df <- merge(results_metrics_df, parameters_combinations[, ])

  prior_parameters_df <- get_parameters(results_df[, "parameters", drop = FALSE])

  param_values <- unlist(get_parameters(unique(results_df[, "parameters", drop = FALSE])))

  posterior_parameters_names <- colnames(get_parameters(unique(results_metrics_df[, "posterior_parameters", drop = FALSE])))

  if (length(posterior_parameters_names) == 1 && posterior_parameters_names == "parameters"){
    return()
  }

  # Remove all elements in posterior_parameters_names that start with "conf_int"
  posterior_parameters_names <- posterior_parameters_names[!grepl("^conf_int", posterior_parameters_names)]

  method_parameters_names <- names(methods_dict[[method]])

  # Select hyperparameters by taking all elements in posterior_parameters_names that are not in method_parameters_names
  selected_parameters <- setdiff(posterior_parameters_names, method_parameters_names)

  parameters_label <- make_labels_from_parameters(get_parameters(parameters_combinations), method)
  parameters_str <- convert_params_to_str(methods_dict[[method]], parameters_combinations)

  if (parameters_label == "") {
    parameters_label_title <- ""
  } else {
    parameters_label_title <- paste0(", ", parameters_label)
  }

  for (key in  selected_parameters) {
    results_metrics_subdf <- results_df

    results_metrics_subdf <- results_df

    posterior_parameters_df <- get_parameters(results_metrics_subdf[, "posterior_parameters", drop = FALSE])

    filter_results_metrics_df <- results_metrics_subdf

    drift <- unlist(filter_results_metrics_df[, xvar[["name"]]])
    y <- unlist(posterior_parameters_df[, key])

    conf_int_lower <- unlist(posterior_parameters_df[, paste0("conf_int_lower_", key)])
    conf_int_upper <- unlist(posterior_parameters_df[, paste0("conf_int_upper_", key)])

    n_digits = 4
    table_data <- data.frame(
      drift = format_num(drift, digits = n_digits),
      y = format_num(y, digits = n_digits),
      ymin = format_num(conf_int_lower, digits = n_digits),
      ymax = format_num(conf_int_upper, digits = n_digits),
      parameters = format_num(filter_results_metrics_df$parameters),
      method = method,
      target_sample_size_per_arm = factor(filter_results_metrics_df$target_sample_size_per_arm)
    )

    # Combine all table data
    table_data_df <- dplyr::bind_rows(table_data)

    data_table <- table_data_df %>%
      dplyr::mutate(
        drift = drift,
        hyperparameter = y,
        CI_low = ymin,
        CI_upper = ymax,
        CI = paste0("[", CI_low, ", ", CI_upper, "]")
      ) %>%
      dplyr::select(drift, target_sample_size_per_arm, hyperparameter, CI) %>%
      dplyr::arrange(drift)


    hyperparameter_colname <- sprintf("%s%s", empirical_bayes_hyperparameters[[method]][[key]][["parameter_notation"]], empirical_bayes_hyperparameters[[method]][[key]][["suffix"]])


    data_table <- data_table %>%
      dplyr::rename(
        "Drift" = drift,
        !!hyperparameter_colname := hyperparameter,
        "95% CI" = CI,
        "$N_T/2$" = target_sample_size_per_arm
      )

    # Merge Mean and CI columns
    data_table <- data_table %>%
      unite(!!hyperparameter_colname, '95% CI', sep = " ")

    # Sort the table
    data_table <- data_table[order(data_table$Drift, data_table[, "$N_T/2$"]), ]
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
      "_",
      empirical_bayes_hyperparameters[[method]][[key]]$parameter_label,
      "_vs_",
      xvar$name,
      "_cat_target_sample_size_per_arm",
      parameters_str
    )


    filename <- format_filename(filename = filename, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)

    file_path <- file.path(directory, filename)

    if (parameters_label_title == ""){
      title <-
        sprintf(
          "%s, %s",
          str_to_title(case_study),
          methods_labels[[method]]$label
        )
    } else {
      title <-
        sprintf(
          "%s, %s%s",
          str_to_title(case_study),
          methods_labels[[method]]$label,
          parameters_label_title
        )
    }
    title <- format_title(title = title, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor, as_latex = TRUE)

    export_table(data_table = data_table, ncollapses = 1, title = title, file_path = file_path)
  }
}



#' This function tables for hyperparameters updated using Empirical Bayes across different case studies, methods, and sample sizes.
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
table_empirical_bayes_hyperparameters_vs_scenario <- function(results_metrics_df, metrics) {
  # Get the list of case studies
  case_studies <- unique(results_metrics_df$case_study)

  # Get the list of methods
  methods <- unique(results_metrics_df$method)

  for (case_study in case_studies) {
    results_metrics_df1 <- results_metrics_df[results_metrics_df$case_study == case_study,]
    for (method in methods) {
      results_metrics_df2 <- results_metrics_df1[results_metrics_df1$method == method,]

      target_to_source_std_ratios <- unique(results_metrics_df2$target_to_source_std_ratio)

      for (target_to_source_std_ratio in target_to_source_std_ratios) {
        results_metrics_df4 <- results_metrics_df2 %>%
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

          # Get the different parameters combinations studies for this method
          parameters_combinations <- unique(results_metrics_df5[results_metrics_df5$method == method, "parameters"])

          for (i in 1:nrow(parameters_combinations)) {
            table_empirical_bayes_hyperparameters_vs_drift(
              results_metrics_df,
              method,
              case_study,
              control_drift = FALSE,
              xvars = xvars,
              source_denominator_change_factor = source_denominator_change_factor,
              target_to_source_std_ratio = target_to_source_std_ratio,
              parameters_combinations = parameters_combinations[i, , drop = FALSE]
            )
          }
        }
      }
    }
  }
}
