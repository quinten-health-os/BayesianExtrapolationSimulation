#' Generate tables for Bayesian metrics
#'
#'
#' @param results_metrics_df A data frame containing the results and metrics.
#' @param case_study Character string specifying the case study.
#' @return This function doesn't return a value but saves the generated tables as HTML, PDF, and LaTeX files.
#'
#' @import dplyr
#' @import knitr
#' @export
table_bayesian_metrics <- function(results_metrics_df,
                                                  case_study,
                                                  target_sample_size_per_arm,
                                                  source_denominator_change_factor = 1,
                                                  target_to_source_std_ratio = 1) {

  results_df <- results_metrics_df %>%
    dplyr::filter(
      case_study == !!case_study,
      target_sample_size_per_arm == !!target_sample_size_per_arm,
      source_denominator_change_factor == !!source_denominator_change_factor  |
        is.na(source_denominator_change_factor),
      target_to_source_std_ratio == !!target_to_source_std_ratio |
        is.na(target_to_source_std_ratio)
    )

  results_df <- results_df %>%
    dplyr::mutate(
      design_prior_type = case_when(
        design_prior_type == "ui_design_prior" ~ "UI design prior",
        design_prior_type == "analysis_prior" ~ "Analysis prior",
        design_prior_type == "source_posterior" ~ "Source posterior",
        TRUE ~ design_prior_type
      )
    )

  results_df$Method <- format_results_df_parameters(results_df)

  results_df <- results_df %>% select(-c("case_study", "method", "parameters", "target_sample_size_per_arm", "source_denominator_change_factor", "target_to_source_std_ratio", "source_treatment_effect_estimate", "source_standard_error", "endpoint", "summary_measure_likelihood", "source_sample_size_control", "source_sample_size_treatment", "equivalent_source_sample_size_per_arm"))

  for (i in seq_along(bayesian_metrics)) {
    metric_name <- bayesian_metrics[[i]]$name

    if (metric_name == "prior_proba_benefit" || metric_name == "prior_proba_no_benefit"){
      # Select the metric of interest
      results_df <- results_df[results_df$design_prior_type == "Analysis prior",]
      data_table <- results_df %>% select(c("Method", metric_name))
      data_wide <- data_table[!is.na(data_table[,metric_name]), ]

      if (metric_name == "prior_proba_benefit"){
        colnames(data_wide)[colnames(data_wide) == metric_name] <- "Prior probability of benefit"
      } else if (metric_name == "prior_proba_no_benefit"){
        colnames(data_wide)[colnames(data_wide) == metric_name] <- "Prior probability of no benefit"
      }
    } else {
      # Select the metric of interest
      data_table <- results_df %>% select(c("Method", "design_prior_type", metric_name))

      # Pivot the table to get one column per design prior type
      data_wide <- data_table %>%
        tidyr::pivot_wider(
          names_from = design_prior_type,
          values_from = metric_name
        )
    }

    directory <- file.path(tables_dir, case_study)
    if (!dir.exists(directory)) {
      dir.create(directory, showWarnings = FALSE, recursive = TRUE)
    }

    filename <- paste0(
      metric_name, "_", case_study,
      "_target_sample_size_per_arm_",
      target_sample_size_per_arm
    )


    title <- sprintf(
        "%s, %s, $N_T/2 = %s$",
        bayesian_metrics[[i]]$label,
        str_to_title(case_study),
        target_sample_size_per_arm
      )

    title <- format_title(title = title, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor, as_latex = TRUE)
    filename <- format_filename(filename = filename, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)

    file_path <- file.path(directory, filename)

    export_table(data_table = data_wide, ncollapses = 1, title = title, file_path = file_path)
  }
}


table_bayesian_metrics_sample_sizes_comparison <- function(results_metrics_df,
                                   case_study,
                                   source_denominator_change_factor = 1,
                                   target_to_source_std_ratio = 1) {

  results_df <- results_metrics_df %>%
    dplyr::filter(
      case_study == !!case_study,
      source_denominator_change_factor == !!source_denominator_change_factor  |
        is.na(source_denominator_change_factor),
      target_to_source_std_ratio == !!target_to_source_std_ratio |
        is.na(target_to_source_std_ratio)
    )

  results_df <- results_df %>%
    dplyr::mutate(
      design_prior_type = case_when(
        design_prior_type == "ui_design_prior" ~ "UI design prior",
        design_prior_type == "analysis_prior" ~ "Analysis prior",
        design_prior_type == "source_posterior" ~ "Source posterior",
        TRUE ~ design_prior_type
      )
    )

  results_df$Method <- format_results_df_parameters(results_df)

  results_df <- results_df %>% select(-c("case_study", "method", "parameters", "source_denominator_change_factor", "target_to_source_std_ratio", "source_treatment_effect_estimate", "source_standard_error", "endpoint", "summary_measure_likelihood", "source_sample_size_control", "source_sample_size_treatment", "equivalent_source_sample_size_per_arm"))

  for (i in seq_along(bayesian_metrics)) {
    metric_name <- bayesian_metrics[[i]]$name

    # Select the metric of interest
    results_df <- results_df[results_df$design_prior_type == "Analysis prior",]
    data_table <- results_df %>% select(c("Method", metric_name, "target_sample_size_per_arm"))
    data_table <- data_table[!is.na(data_table[,metric_name]), ]

    # Convert to wide format
    wide_data_table <- data_table %>%
      pivot_wider(
        names_from = target_sample_size_per_arm, # Columns become values of `target_sample_size_per_arm`
        values_from = metric_name               # Values filled from `metric_name`
      )

    wide_data_table <- wide_data_table %>%
      rename_with(
        ~ ifelse(grepl("^[0-9]+$", .), paste0("$N_T/2 = ", ., "$"), .)
      )

    wide_data_table_filtered <- wide_data_table %>%
      select(where(~ all(!is.na(.))))


    directory <- file.path(tables_dir, case_study)
    if (!dir.exists(directory)) {
      dir.create(directory, showWarnings = FALSE, recursive = TRUE)
    }

    filename <- paste0(
      case_study, "_", metric_name, "_analysis_prior"
    )

    title <- sprintf(
      "%s, %s",
      bayesian_metrics[[i]]$label,
      str_to_title(case_study)
    )

    title <- format_title(title = title, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor, as_latex = TRUE)
    filename <- format_filename(filename = filename, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)

    file_path <- file.path(directory, filename)

    export_table(data_table = wide_data_table_filtered, ncollapses = 1, title = title, file_path = file_path)
  }
}


table_bayesian_metrics_across_case_studies <- function(results_metrics_df) {

  source_denominator_change_factor <- 1
  target_to_source_std_ratio <- 1

  results_df <- results_metrics_df %>%
    dplyr::filter(
      source_denominator_change_factor == !!source_denominator_change_factor  |
        is.na(source_denominator_change_factor),
      target_to_source_std_ratio == !!target_to_source_std_ratio |
        is.na(target_to_source_std_ratio)
    )


  results_df <- results_df %>%
    dplyr::mutate(
      design_prior_type = case_when(
        design_prior_type == "ui_design_prior" ~ "UI design prior",
        design_prior_type == "analysis_prior" ~ "Analysis prior",
        design_prior_type == "source_posterior" ~ "Source posterior",
        TRUE ~ design_prior_type
      )
    )

  results_df$Method <- format_results_df_parameters(results_df)

  results_df <- results_df %>% select(-c("method", "parameters", "source_denominator_change_factor", "target_to_source_std_ratio", "source_treatment_effect_estimate", "source_standard_error", "endpoint", "summary_measure_likelihood", "source_sample_size_control", "source_sample_size_treatment", "equivalent_source_sample_size_per_arm"))



  for (i in seq_along(bayesian_metrics)) {
    metric_name <- bayesian_metrics[[i]]$name


    if (metric_name != "prior_proba_benefit" && metric_name != "prior_proba_no_benefit"){
      next
    }

    # Select the metric of interest
    results_df <- results_df[results_df$design_prior_type == "Analysis prior",]
    data_table <- results_df %>% select(c("Method", metric_name, "case_study"))
    data_table <- data_table[!is.na(data_table[,metric_name]), ]

    # Remove rows where metric_name is NA
    data_table <- data_table %>%
      filter(!is.na(metric_name))

    # The prior proba of (no) benefit does not depend on the sample size
    data_table <- data_table %>%
      dplyr::distinct(case_study, Method, .keep_all = TRUE)


    # Convert to wide format
    wide_data_table <- data_table %>%
      pivot_wider(
        names_from = c('case_study'),
        values_from = metric_name
      )

    wide_data_table_filtered <- wide_data_table %>%
      select(where(~ all(!is.na(.))))

    wide_data_table_filtered  <- wide_data_table_filtered  %>%
      rename_with(~ paste0(toupper(substr(., 1, 1)), substr(., 2, nchar(.))))

    directory <- file.path(tables_dir)
    if (!dir.exists(directory)) {
      dir.create(directory, showWarnings = FALSE, recursive = TRUE)
    }

    filename <- paste0(metric_name, "_analysis_prior"
    )

    title <- bayesian_metrics[[i]]$label

    # title <- format_title(title = title, case_study = NA, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor, as_latex = TRUE)
    filename <- format_filename(filename = filename, case_study = NA, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)

    file_path <- file.path(directory, filename)

    export_table(data_table = wide_data_table_filtered, ncollapses = 1, title = title, file_path = file_path)
  }
}

#' Generate tables for methods operating characteristics
#'
#' This function creates tables for various Bayesian operating characteristics across different case studies and methods.
#'
#' @param results_metrics_df A data frame containing the results and metrics.
#' @param metrics A vector of metric names to include in the tables.
#'
#' @return This function doesn't return a value but calls other functions to generate and save tables.
#'
#' @import dplyr
#' @import yaml
#'
#' @export
table_bayesian_ocs <- function(results_metrics_df, metrics, source_denominator_change_factor = 1,
                               target_to_source_std_ratio = 1) {
  case_studies <- unique(results_metrics_df$case_study)
  methods <- unique(results_metrics_df$method)

  results_metrics_df <- results_metrics_df %>%
    dplyr::filter(
      source_denominator_change_factor == !!source_denominator_change_factor  |
        is.na(source_denominator_change_factor),
      target_to_source_std_ratio == !!target_to_source_std_ratio |
        is.na(target_to_source_std_ratio)
    )

  table_bayesian_metrics_across_case_studies(
    results_metrics_df =  results_metrics_df
  )

  for (case_study in case_studies) {
    results_df <- results_metrics_df %>%
      dplyr::filter(case_study == !!case_study)

    target_sample_sizes_per_arm <- unique(results_df$target_sample_size_per_arm[results_df$case_study == case_study])

    for (source_denominator_change_factor in unique(results_df$source_denominator_change_factor)) {
      for (target_to_source_std_ratio in unique(results_df$target_to_source_std_ratio)) {
        table_bayesian_metrics_sample_sizes_comparison(
          results_metrics_df = results_df,
          case_study = case_study,
          source_denominator_change_factor = source_denominator_change_factor,
          target_to_source_std_ratio = target_to_source_std_ratio
        )

        for (sample_size_per_arm in target_sample_sizes_per_arm) {

         table_bayesian_metrics(
                  results_metrics_df = results_df,
                  case_study = case_study,
                  target_sample_size_per_arm = sample_size_per_arm,
                  source_denominator_change_factor = source_denominator_change_factor,
                  target_to_source_std_ratio = target_to_source_std_ratio
                )
        }
      }
    }
  }
}
