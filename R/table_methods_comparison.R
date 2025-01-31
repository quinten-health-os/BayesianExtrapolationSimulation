#' Generate a comparison table
#'
#' This function generates a comparison table for different methods and treatment effects.
#'
#' @param results_df A dataframe containing the results to be compared.
#' @param case_study The case study to filter the results.
#' @param target_sample_size_per_arm The target sample size per arm to filter the results.
#' @param x_metric The metric to be used for comparison.
#'
#' @return A kable object representing the comparison table.
#'
#' @import dplyr
#' @import knitr
#' @export
generate_comparison_table <- function(results_df,
                                      x_metric, source_denominator_change_factor, target_to_source_std_ratio) {
  if (x_metric %in% names(frequentist_metrics)) {
    metric <- frequentist_metrics[[x_metric]]
  } else if (x_metric %in% names(inference_metrics)) {
    metric <- inference_metrics[[x_metric]]
  } else if (x_metric %in% names(bayesian_metrics)) {
    metric <- bayesian_metrics[[x_metric]]
  } else {
    stop("Metric is not supported")
  }


  x_metric_name <- metric$name
  x_metric_uncertainty_lower <- paste0(metric$metric_uncertainty, "_lower")
  x_metric_uncertainty_upper <- paste0(metric$metric_uncertainty, "_upper")
  x_metric_label <- metric$label


  # Calculate important target treatment effects
  important_target_treatment_effects <- round(results_df %>%
                                                pull(source_treatment_effect_estimate) %>%
                                                unique() %>%
                                                first() %>%
                                                `*`(c(0, 1 / 2, 1)),
                                              4)

  df <- results_df %>%
    dplyr::mutate(
      means = !!rlang::sym(x_metric_name),
      error_low = !!rlang::sym(x_metric_uncertainty_lower),
      error_upper = !!rlang::sym(x_metric_uncertainty_upper)
    )

  # Find the closest values in df$target_treatment_effect to important_target_treatment_effects
  closest_values <- sapply(important_target_treatment_effects, function(x) {
    df$target_treatment_effect[which.min(abs(df$target_treatment_effect - x))]
  })
  closest_values <- round(closest_values, digits = 4)


  df <- df %>%
    mutate(target_treatment_effect = round(target_treatment_effect, digits = 4))

  # Filter the dataframe to keep only the rows with the closest values
  df <- df %>%
    dplyr::filter(target_treatment_effect %in% closest_values) %>%
    dplyr::arrange(target_treatment_effect)

  if (length(unique(df$target_treatment_effect)) != length(important_target_treatment_effects)) {
    stop(
      "The filtered dataframe cannot contain more unique target treatment effect values than there are in important_target_treatment_effects"
    )
  }

  # Factorize target treatment effect
  df$target_treatment_effect <- factor(
    df$target_treatment_effect,
    levels = unique(df$target_treatment_effect),
    labels = c(
      "No treatment effect",
      "Partially consistent treatment effect",
      "Consistent treatment effect"
    )
  )


  df$Method <- format_results_df_parameters(df)

  # Create separate data frames for each effect
  no_effect_data <- df[df[, "target_treatment_effect"] == "No treatment effect", ]
  partially_consistent_effect_data <- df[df[, "target_treatment_effect"] == "Partially consistent treatment effect", ]
  consistent_effect_data <- df[df[, "target_treatment_effect"] == "Consistent treatment effect", ]

  # Combine the data into a single table
  combined_data <- dplyr::bind_rows(
    no_effect_data %>% dplyr::mutate(effect = "No treatment effect"),
    partially_consistent_effect_data %>% dplyr::mutate(effect = "Partially consistent"),
    consistent_effect_data %>% dplyr::mutate(effect = "Consistent")
  )

  # Select relevant columns and create the table
  data_table <- combined_data %>%
    dplyr::select(Method, effect, means, error_low, error_upper) %>%
    dplyr::mutate(
      means = format_num(means),
      error_low = format_num(error_low),
      error_upper = format_num(error_upper),
      CI = paste0("[", error_low, ", ", error_upper, "]")
    ) %>%
    dplyr::select(Method, effect, means, CI) %>%
    dplyr::arrange(Method, effect) # Order by labels and effect

  data_table <- data_table %>%
    dplyr::rename(
      "Method" = Method,
      "Treatment Effect" = effect,
      !!x_metric_label := means,
      !!metric$uncertainty_label := CI
    )

  case_study <- unique(df$case_study)
  target_sample_size_per_arm <- unique(df$target_sample_size_per_arm)
  target_to_source_std_ratio <- unique(df$target_to_source_std_ratio)
  source_denominator_change_factor <- unique(df$source_denominator_change_factor)

  title <- sprintf(
      "%s, $N_T/2 = $ %s",
      str_to_title(case_study),
      target_sample_size_per_arm
    )
  title <- format_title(title = title, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor, as_latex = TRUE)

  directory <- file.path(tables_dir, case_study)
  if (!dir.exists(directory)) {
    dir.create(directory, showWarnings = FALSE, recursive = TRUE)
  }

  filename <- paste0(
    case_study,
    "_",
    x_metric,
    "_table_target_sample_size_per_arm_",
    target_sample_size_per_arm
  )

  filename <- format_filename(filename = filename, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)

  file_path <- file.path(directory, filename)

  # Wide table version
  wide_data_table <- data_table %>%
    pivot_wider(
      names_from = c("Treatment Effect"),
      values_from = c(x_metric_label, metric$uncertainty_label),
      names_sep = " "
    ) %>%
    unite("No treatment effect", c(paste(x_metric_label, "No treatment effect"), paste(metric$uncertainty_label, "No treatment effect")), sep = " ", remove = TRUE) %>%
    unite("Partially consistent", c(paste(x_metric_label, "Partially consistent"), paste(metric$uncertainty_label, "Partially consistent")), sep = " ", remove = TRUE) %>%
    unite("Consistent", c(paste(x_metric_label, "Consistent"), paste(metric$uncertainty_label, "Consistent")), sep = " ", remove = TRUE)

  export_table(data_table = wide_data_table, ncollapses = 1, column_names = c("Method", "Consistent", "No treatment effect", "Partially consistent"), title = paste(metric$label,  title, sep =", "), file_path = file_path)
}


#' Methods comparison in a table
#'
#' This function creates tables comparing methods with respect to different metrics for various case studies and sample sizes.
#'
#' @param results_metrics_df The results metrics dataframe.
#' @param metrics List of frequentist metrics and their properties.
#'
#' @return None
#'
#' @import dplyr
#'
#' @export
table_methods_comparison <- function(results_metrics_df, metrics) {
  if (nrow(results_metrics_df) == 0) {
    stop("The dataframe is empty")
  }

  # Get the list of case studies
  case_studies <- unique(results_metrics_df$case_study)


  for (metric in names(metrics)) {
    for (case_study in case_studies) {
      results_df1 <- results_metrics_df %>% dplyr::filter(
        case_study == !!case_study
      )
      target_sample_sizes <- unique(results_df1$target_sample_size_per_arm)
      for (target_sample_size_per_arm in target_sample_sizes) {
        results_df2 <- results_df1 %>% dplyr::filter(
          target_sample_size_per_arm == !!target_sample_size_per_arm
        )
        source_denominator_change_factors <- unique(results_df2$source_denominator_change_factor)
        for (source_denominator_change_factor in source_denominator_change_factors){
          results_df3 <- results_df2 %>% dplyr::filter(
            source_denominator_change_factor == !!source_denominator_change_factor
          )
          target_to_source_std_ratios <- unique(results_df3$target_to_source_std_ratio)
          for (target_to_source_std_ratio in target_to_source_std_ratios){
            results_df4 <- results_df3 %>% dplyr::filter(
              target_to_source_std_ratio == !!target_to_source_std_ratio |
                is.na(target_to_source_std_ratio)
            )
            generate_comparison_table(results_df4, metric, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)
          }
        }
      }
    }
  }
}
