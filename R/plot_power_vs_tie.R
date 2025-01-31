#' Function to generate a power vs tie plot
#'
#' @param results_metrics_df The dataframe containing the results and metrics
#' @param case_study The case study name
#' @param target_sample_size_per_arm The target sample size per arm
#' @param treatment_effect The treatment effect type ("consistent", "no_effect", "partially_consistent")
#' @param power_difference Logical indicating whether to calculate power difference
#' @param metrics The metrics to include in the plot
#'
#' @return None
#'
#' @export
power_vs_tie <- function(results_metrics_df,
                         case_study,
                         target_sample_size_per_arm,
                         treatment_effect,
                         power_difference,
                         metrics,
                         source_denominator_change_factor,
                         target_to_source_std_ratio) {

  directory <- paste0(figures_dir, case_study)
  if (!dir.exists(directory)) {
    dir.create(directory, showWarnings = FALSE, recursive = TRUE)
  }

  results_df <- results_metrics_df %>%
    dplyr::filter(
      target_sample_size_per_arm == !!target_sample_size_per_arm,
      case_study == !!case_study,
      source_denominator_change_factor == !!source_denominator_change_factor  |
        is.na(source_denominator_change_factor),
      target_to_source_std_ratio == !!target_to_source_std_ratio |
        is.na(target_to_source_std_ratio)
    )

  labels <- format_results_df_parameters(results_df)
  methods <- format_results_df_methods(results_df)

  results_df$label <- labels
  results_df$method <- methods

  theta_0 <- unique(results_df$theta_0)
  results_df_tie <- results_df %>% dplyr::filter(target_treatment_effect == theta_0)

  if (treatment_effect == "consistent") {
    results_df <- results_df %>%
      dplyr::filter(target_treatment_effect == source_treatment_effect_estimate)
  } else if (treatment_effect == "no_effect") {
    stop("The treatment effect should not be non-positive")
  } else if (treatment_effect == "partially_consistent") {
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
    power_label <- "power_difference"
  } else {
    power_label <- "power"
  }

  filename <- paste0(
    case_study,
    "_",
    power_label,
    "_sample_size=",
    target_sample_size_per_arm,
    "_treatment_effect=",
    treatment_effect
  )

  filename <- format_filename(filename = filename, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)

  file_path <- file.path(directory, filename)

  if (file.exists(paste0(file_path, ".pdf")) && file.exists(paste0(file_path, ".png")) && remake_figures == FALSE){
    return()
  }


  if (power_difference == FALSE) {
    ymin <- 0
    ymax <- 1
  } else if (power_difference == TRUE) {
    y_range <- sort(range(results_df[, "success_proba"]))
    ymin <- y_range[1]
    ymax <- y_range[2]
  }
  # xmin <- 0
  # xmax <- 1
  x_range <- sort(range(results_df[, "tie"]))
  xmin <- x_range[1]
  xmax <- x_range[2]

  x_width <- xmax - xmin
  y_width <- ymax - ymin

  y_cap_size <- x_width * relative_error_cap_width / 2 # Adjust relative_error_cap_width to control the relative cap size
  x_cap_size <- y_width * relative_error_cap_width / 2 # Adjust relative_error_cap_width to control the relative cap size

  cap_size <- min(x_cap_size, y_cap_size)
  x_cap_size <- cap_size
  y_cap_size <- cap_size

  treatment_effects_labels <- list(
    no_effect = "No effect",
    partially_consistent = "Partially consistent effect",
    consistent = "Consistent effect"
  )

  plot_title <- sprintf(
      "%s, %s %s, %s",
      str_to_title(case_study),
      "$N_T/2 = $",
      target_sample_size_per_arm,
      treatment_effects_labels[[treatment_effect]]
    )

  plot_title <- format_title(title = plot_title, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)

  plt <- ggplot(results_df, aes(x = tie)) +
    # coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
    geom_point(aes(
      y = success_proba,
      color = as.factor(method),
      shape = as.factor(method)
    )) +
    geom_errorbar(
      aes(
        ymin = conf_int_success_proba_lower,
        ymax = conf_int_success_proba_upper,
        color =  as.factor(method)
      ),
      width = y_cap_size
    ) +
    geom_errorbarh(
      aes(
        xmin = conf_int_tie_lower,
        xmax = conf_int_tie_upper,
        y = success_proba,
        color = as.factor(method)
      ),
      height = x_cap_size
    ) +
    ggplot2::labs(
      title = plot_title,
      x = "TIE",
      y = ifelse(
        power_difference,
        "Power difference compared to a frequentist test",
        "Power"
      ),
      color = "Methods",
      # Title for the color legend
      shape = "Methods" # Title for the shape legend
    )

  if (power_difference == FALSE) {
    # Add the nominal frequentist power and power at equivalent TIE
    plt <- plt +
      geom_point(
      aes(
        y = !!sym("nominal_frequentist_power_separate"),
        color = "Nominal frequentist power",
        shape = "Nominal frequentist power"
      )
    ) +
      geom_point(
        aes(
          y = !!sym("frequentist_power_at_equivalent_tie"),
          color = "Frequentist power at equivalent TIE",
          shape = "Frequentist power at equivalent TIE"
        )
      )

    plt <- plt +
      ggplot2::geom_vline(xintercept = analysis_config$nominal_tie, color = "black", linetype = "dashed") + # Add horizontal line at x = 0.025
      ggplot2::scale_x_continuous(breaks = function(x) unique(c(pretty(x), analysis_config$nominal_tie)))  # Add x-tick at x = 0.025
  }

  plt <- plt +
    scale_color_viridis_d() # Apply viridis color scale

  # Save plot

  plot.size <- set_size(textwidth)
  fig_width_in <- plot.size[1]
  fig_height_in <- plot.size[2]

  export_plots(plt, file_path, fig_width_in, fig_height_in, type = "pdf")
  export_plots(plt, file_path, fig_width_in, fig_height_in, type = "png")
}

#' Plot frequentist methods operating characteristics
#'
#' @description This function generates plots for the operating characteristics of frequentist methods.
#'
#' @param results_metrics_df The data frame containing the results and metrics.
#' @param metrics A list of metrics to be plotted.
#'
#' @return None
#'
#' @examples NA
power_vs_tie_plots <- function(results_metrics_df, metrics) {
  # Get the list of case studies
  case_studies <- unique(results_metrics_df$case_study)
  treatment_effects <- c("partially_consistent", "consistent")

  for (case_study in case_studies) {
    results_metrics_df1 <- results_metrics_df[results_metrics_df$case_study == case_study,]

    # Get the list of target sample sizes
    target_sample_sizes <- unique(results_metrics_df1$target_sample_size_per_arm[results_metrics_df1$case_study == case_study])

    for (target_sample_size in target_sample_sizes) {
      results_metrics_df2 <- results_metrics_df1[results_metrics_df1$target_sample_size_per_arm == target_sample_size,]

      for (treatment_effect in treatment_effects) {

        target_to_source_std_ratios <- unique(results_metrics_df2$target_to_source_std_ratio)

        for (target_to_source_std_ratio in target_to_source_std_ratios) {
          results_metrics_df3 <- results_metrics_df2 %>%
            dplyr::filter(
              target_to_source_std_ratio == !!target_to_source_std_ratio |
                is.na(target_to_source_std_ratio)
            )

          source_denominator_change_factors <- unique(results_metrics_df3$source_denominator_change_factor)

          for (source_denominator_change_factor in source_denominator_change_factors){
            results_metrics_df4 <- results_metrics_df3 %>%
              dplyr::filter(
                source_denominator_change_factor == !!source_denominator_change_factor |
                  is.na(source_denominator_change_factor)
              )

            power_vs_tie(
              results_metrics_df4,
              case_study,
              target_sample_size,
              treatment_effect = treatment_effect,
              power_difference = FALSE,
              metrics = metrics,
              source_denominator_change_factor = source_denominator_change_factor,
              target_to_source_std_ratio = target_to_source_std_ratio
            )

            power_vs_tie(
              results_metrics_df4,
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
