#' Generate a power vs tie table
#'
#' This function creates a table comparing power (or power difference) and type I error (TIE)
#' for different methods under specific conditions.
#'
#' @param results_metrics_df A data frame containing the results and metrics.
#' @param case_study The specific case study to filter the results.
#' @param target_sample_size_per_arm The target sample size per arm to filter the results.
#' @param treatment_effect Character string specifying the treatment effect: "consistent", "partially_consistent", or "no_effect".
#' @param power_difference Logical; if TRUE, calculate power difference instead of power.
#' @param metrics A vector of metric names (not used in the current implementation).
#'
#' @return This function doesn't return a value but generates and saves tables in various formats (HTML, PDF, LaTeX).
#'
#' @import dplyr
#' @import knitr
#' @import kableExtra
#'
#' @export
table_power_vs_tie <- function(results_metrics_df,
                               case_study,
                               target_sample_size_per_arm,
                               treatment_effect,
                               power_difference,
                               metrics,
                               source_denominator_change_factor,
                               target_to_source_std_ratio) {
  results_df <- results_metrics_df %>%
    dplyr::filter(
      target_sample_size_per_arm == !!target_sample_size_per_arm,
      case_study == !!case_study,
      source_denominator_change_factor == !!source_denominator_change_factor  |
        is.na(source_denominator_change_factor),
      target_to_source_std_ratio == !!target_to_source_std_ratio |
        is.na(target_to_source_std_ratio)
    )

  if (nrow(results_df) == 0){
    return()
  }

  results_df$Method <- format_results_df_parameters(results_df)


  theta_0 <- unique(results_df$theta_0)
  results_df_tie <- results_df %>% dplyr::filter(target_treatment_effect == theta_0)

  if (treatment_effect == "consistent") {
    results_df <- results_df %>%
      dplyr::filter(target_treatment_effect == source_treatment_effect_estimate)
    treatment_effect_label <- "consistent treatment effect"
  } else if (treatment_effect == "no_effect") {
    stop("The treatment effect should not be non-positive")
  } else if (treatment_effect == "partially_consistent") {
    treatment_effect_label <- "partially consistent treatment effect"
    results_df <- results_df %>%
      dplyr::filter(abs(
        target_treatment_effect - source_treatment_effect_estimate / 2
      ) < 1e-4)
  } else {
    stop('treatment_effect must be "consistent", "no_effect" or "partially_consistent"')
  }

  if (power_difference) {
    results_df <- results_df %>% dplyr::mutate(
      success_proba = success_proba - frequentist_power_at_equivalent_tie,
      conf_int_success_proba_lower = conf_int_success_proba_lower - frequentist_power_at_equivalent_tie,
      conf_int_success_proba_upper = conf_int_success_proba_upper - frequentist_power_at_equivalent_tie
    )
    power_label <- "Power difference"
  } else {
    power_label <- "Power"
  }


  data_table <- results_df %>%
    dplyr::select(
      Method,
      tie,
      conf_int_tie_lower,
      conf_int_tie_upper,
      success_proba,
      conf_int_success_proba_lower,
      conf_int_success_proba_upper
    ) %>%
    dplyr::mutate(
      tie = format_num(tie),
      error_low_tie = format_num(conf_int_tie_lower),
      error_upper_tie = format_num(conf_int_tie_upper),
      CI_tie = paste0("[", error_low_tie, ", ", error_upper_tie, "]"),
      success_proba = format_num(success_proba),
      error_low_power = format_num(conf_int_success_proba_lower),
      error_upper_power = format_num(conf_int_success_proba_upper),
      CI_power = paste0("[", error_low_power, ", ", error_upper_power, "]")
    ) %>%
    dplyr::select(Method, tie, CI_tie, success_proba, CI_power) %>%
    dplyr::arrange(Method, tie)


  data_table <- data_table %>%
    # dplyr::group_by(method) %>%
    dplyr::rename(
      "Method" = Method,
      "TIE" = tie,
      "Confidence Interval TIE" = CI_tie,
      !!power_label := success_proba,
      "Confidence Interval Power" = CI_power
    )

  # Merge Power and CI columns
  data_table <- data_table %>%
    unite(!!power_label, power_label, 'Confidence Interval Power', sep = " ")

  data_table <- data_table %>%
    unite("TIE", "TIE", 'Confidence Interval TIE', sep = " ")

  directory <- file.path(tables_dir, case_study)
  if (!dir.exists(directory)) {
    dir.create(directory, showWarnings = FALSE, recursive = TRUE)
  }

  if (power_difference) {
    power_label <- "_power_difference_vs_tie_"
  } else {
    power_label <-  "_power_vs_tie_"
  }

  filename <- paste0(
        case_study,
        "_",
        power_label,
        "target_sample_size_per_arm_",
        target_sample_size_per_arm,
        "_",
        treatment_effect
      )

  title <-  sprintf(
      "%s, $N_T/2 = %s$, %s",
      str_to_title(case_study),
      target_sample_size_per_arm,
      treatment_effect_label
    )

  filename <- format_filename(filename = filename, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)

  title <- format_title(title = title, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor, as_latex = TRUE)

  file_path <- file.path(directory, filename)

  export_table(data_table = data_table, ncollapses = 1, title = title, file_path = file_path)
}

#' Generate power vs tie tables for multiple conditions
#'
#' This function generates multiple power vs tie tables for various combinations of
#' case studies, sample sizes, and treatment effects.
#'
#' @param results_metrics_df A data frame containing the results and metrics.
#' @param metrics A vector of metric names (not used in the current implementation).
#'
#' @return This function doesn't return a value but generates multiple tables using table_power_vs_tie().
#'
#' @import dplyr
#'
#' @export
power_vs_tie_tables <- function(results_metrics_df, metrics) {
  # Get the list of case studies
  case_studies <- unique(results_metrics_df$case_study)

  treatment_effects <- c("partially_consistent", "consistent")

  for (case_study in case_studies) {
    # Get the list of target sample sizes
    target_sample_sizes <- unique(results_metrics_df$target_sample_size_per_arm[results_metrics_df$case_study == case_study])

    for (target_sample_size in target_sample_sizes) {
      for (treatment_effect in treatment_effects) {
        for (source_denominator_change_factor in unique(results_metrics_df$source_denominator_change_factor)){
          for (target_to_source_std_ratio in unique(results_metrics_df$target_to_source_std_ratio)){
            table_power_vs_tie(
              results_metrics_df,
              case_study,
              target_sample_size,
              treatment_effect = treatment_effect,
              power_difference = FALSE,
              metrics = metrics,
              source_denominator_change_factor = source_denominator_change_factor,
              target_to_source_std_ratio = target_to_source_std_ratio
            )

            table_power_vs_tie(
              results_metrics_df,
              case_study,
              target_sample_size,
              treatment_effect = treatment_effect,
              power_difference = TRUE,
              metrics = metrics,
              source_denominator_change_factor = source_denominator_change_factor,
              target_to_source_std_ratio = target_to_source_std_ratio
            )
          }
        }
      }
    }
  }
}
