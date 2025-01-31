plot_metric_vs_ess <- function(results_metrics_df,
                         case_study,
                         target_sample_size_per_arm,
                         treatment_effect,
                         ess_method,
                         metric,
                         source_denominator_change_factor = 1,
                         target_to_source_std_ratio = 1) {

  if (!(metric$name %in% c("mse", "bias", "coverage", "precision"))){
   return()
  }

  directory <- paste0(figures_dir, case_study)
  if (!dir.exists(directory)) {
    dir.create(directory, showWarnings = FALSE, recursive = TRUE)
  }

  filename <- paste0(
    case_study,
    "_",
    metric$name,
    "_vs_",
    ess_method$name,
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

  if (treatment_effect == "consistent") {
    results_df <- results_df %>%
      dplyr::filter(target_treatment_effect == source_treatment_effect_estimate)
  } else if (treatment_effect == "no_effect") {
    theta_0 <- unique(results_df$theta_0)
    results_df <- results_df %>%
      dplyr::filter(target_treatment_effect == theta_0)
  } else if (treatment_effect == "partially_consistent") {
    results_df <- results_df %>%
      dplyr::filter(abs(
        target_treatment_effect - source_treatment_effect_estimate / 2
      ) < 1e-4)
  } else {
    stop('treatment_effect must be "consistent", "no_effect" or "partially_consistent"')
  }

  y_range <- sort(range(results_df[, metric$name]))
  ymin <- y_range[1]
  ymax <- y_range[2]

  x_range <- sort(range(results_df[,  ess_method$name]))
  xmin <- x_range[1]
  xmax <- x_range[2]

  x_width <- xmax - xmin
  y_width <- ymax - ymin

  y_cap_size <- x_width * relative_error_cap_width / 2 # Adjust relative_error_cap_width to control the relative cap size
  x_cap_size <- y_width * relative_error_cap_width / 2 # Adjust relative_error_cap_width to control the relative cap size

  # cap_size <- min(x_cap_size, y_cap_size)
  # x_cap_size <- cap_size
  # y_cap_size <- cap_size

  treatment_effects_labels <- list(
    no_effect = "No effect",
    partially_consistent = "Partially consistent effect",
    consistent = "Consistent effect"
  )

  plot_title <- sprintf(
    "%s, $N_T/2 = $ %s, %s",
    str_to_title(case_study),
    target_sample_size_per_arm,
    treatment_effects_labels[[treatment_effect]]
  )

  plot_title <- format_title(title = plot_title, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)


  metric_name <- rlang::sym(metric$name)
  metric_uncertainty_lower <- rlang::sym(paste0(metric$metric_uncertainty, "_lower"))
  metric_uncertainty_upper <- rlang::sym(paste0(metric$metric_uncertainty, "_upper"))
  metric_label <- metric$label

  ess_name <- rlang::sym(ess_method$name)
  ess_uncertainty_lower <- rlang::sym(paste0(ess_method$metric_uncertainty, "_lower"))
  ess_uncertainty_upper <- rlang::sym(paste0(ess_method$metric_uncertainty, "_upper"))
  ess_label <- ess_method$label


  plt <- ggplot(results_df, aes(x = !!ess_name)) +
    geom_point(aes(
      y = !!metric_name,
      color = as.factor(method),
      shape = as.factor(method)
    )) +
    geom_errorbar(
      aes(
        ymin = !!metric_uncertainty_lower,
        ymax = !!metric_uncertainty_upper,
        color =  as.factor(method)
      ),
      width = y_cap_size
    ) +
    geom_errorbarh(
      aes(
        xmin = !!ess_uncertainty_lower,
        xmax = !!ess_uncertainty_upper,
        y = !!metric_name,
        color = as.factor(method)
      ),
      height = x_cap_size
    ) +
    ggplot2::labs(
      title = plot_title,
      x = ess_label,
      y = metric_label,
      color = "Methods",
      # Title for the color legend
      shape = "Methods" # Title for the shape legend
    )

  plt <- plt +
    scale_color_viridis_d()   # Apply viridis color scale

  # Save plot
  plot.size <- set_size(textwidth)
  fig_width_in <- plot.size[1]
  fig_height_in <- plot.size[2]

  export_plots(plt, file_path, fig_width_in, fig_height_in, type = "pdf")
  export_plots(plt, file_path, fig_width_in, fig_height_in, type = "png")
}



#' Plot methods operating characteristics
#'
#' @description This function generates plots for the operating characteristics of different methods.
#'
#' @param results_metrics_df The data frame containing the results and metrics.
#' @param metrics A list of metrics to be plotted.
#'
#' @return None
#'
#' @examples NA
plot_metrics_vs_ess <- function(results_metrics_df, frequentist_metrics, inference_metrics) {
  # Get the list of case studies
  case_studies <- unique(results_metrics_df$case_study)

  ess_methods <- c("ess_elir", "ess_moment", "ess_precision")

  treatment_effects <-  treatment_effects <- c("no_effect", "partially_consistent", "consistent")

  for (case_study in case_studies) {
    case_study_config <- yaml::yaml.load_file(system.file(
      paste0("conf/case_studies/", case_study, ".yml"),
      package = "RBExT"
    ))

    filtered_results_metrics_df <- results_metrics_df %>%
      dplyr::filter(case_study == !!case_study)

    if (nrow(filtered_results_metrics_df) == 0) {
      warning("Dataframe is empty")
      next
    }


    target_sample_sizes <- unique(filtered_results_metrics_df$target_sample_size_per_arm)

    for (metric in frequentist_metrics) {
      for (sample_size in target_sample_sizes) {
        for (ess_method in  ess_methods){
          for (treatment_effect in treatment_effects){
            plot_metric_vs_ess(
              results_metrics_df = filtered_results_metrics_df,
              case_study = case_study,
              ess_method = inference_metrics[[ess_method]],
              target_sample_size_per_arm = sample_size,
              treatment_effect = treatment_effect,
              metric = metric,
            )
          }
        }
      }
    }
  }
}
