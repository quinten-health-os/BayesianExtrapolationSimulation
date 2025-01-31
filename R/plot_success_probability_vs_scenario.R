plot_baseline_success_proba_vs_xvar <- function(plt,
                                                data,
                                                xvar_name,
                                                selected_metric_name,
                                                cap_size,
                                                markersize,
                                                join_points,
                                                label = NULL) {
  if (nrow(data) == 0) {
    stop("The dataframe is empty.")
  }

  plt <- plt +
    geom_point(
      data = data,
      ggplot2::aes(
        x = !!sym(xvar_name),
        y = !!sym(selected_metric_name),
        color = label
      ),
      position = position_dodge(width = cap_size),
      size = markersize,
    )

  if (join_points == TRUE) {
    plt <- plt + geom_line(
      #  Add lines between points
      data = data,
      ggplot2::aes(
        x = !!sym(xvar_name),
        y = !!sym(selected_metric_name),
        color = label
      ),
      position = position_dodge(width = cap_size)
    )
  }
  return(plt)
}

plot_success_proba_vs_xvar <- function(plt,
                                       data,
                                       xvar_name,
                                       selected_metric_name,
                                       selected_metric_uncertainty_lower,
                                       selected_metric_uncertainty_upper,
                                       cap_size,
                                       markersize,
                                       join_points,
                                       label = NULL) {
  if (nrow(data) == 0) {
    warning("The dataframe is empty.")
    return()
  }

  plt <- plt +
    geom_errorbar(
      data = data,
      ggplot2::aes(
        x = !!sym(xvar_name),
        y = !!sym(selected_metric_name),
        ymin = !!sym(selected_metric_uncertainty_lower),
        ymax = !!sym(selected_metric_uncertainty_upper),
        color = label
      ),
      position = position_dodge(width = cap_size),
      width = cap_size
    ) +
    geom_point(
      data = data,
      ggplot2::aes(
        x = !!sym(xvar_name),
        y = !!sym(selected_metric_name),
        color = label
      ),
      position = position_dodge(width = cap_size),
      size = markersize,
    )
  if (join_points == TRUE) {
    plt <- plt + geom_line(
      #  Add lines between points
      data = data,
      ggplot2::aes(
        x = !!sym(xvar_name),
        y = !!sym(selected_metric_name),
        color = label
      ),
      position = position_dodge(width = cap_size)
    )
  }

  plt <- plt +
    ylim(0, 1)
  return(plt)
}


#' Function to generate a metric vs drift plot
#'
#' @param metric The metric to plot
#' @param results_metrics_df The dataframe containing the results and metrics
#' @param theta_0 The true treatment effect
#' @param case_study The case study name
#' @param method The method name
#' @param target_sample_size_per_arm The target sample size per arm
#' @param parameters_combinations The combinations of parameters to filter the results by
#' @param xvars A list of x-variables for control drift and treatment drift
#' @param target_to_source_std_ratio Ratio between the target and source studies sampling standard deviation.
#' @param join_points Whether to join the point with a line or not.
#'
#' @return None
#'
#' @export
plot_success_proba_vs_drift <- function(metric,
                                        results_metrics_df,
                                        theta_0,
                                        case_study,
                                        method,
                                        target_sample_size_per_arm = NULL,
                                        target_to_source_std_ratio,
                                        source_denominator_change_factor = NULL,
                                        parameters_combinations = NULL,
                                        xvars,
                                        join_points = FALSE,
                                        baseline_success_proba) {


  xvar <- xvars$drift

  results_df <- results_metrics_df %>%
    dplyr::filter(
      method == !!method,
      target_sample_size_per_arm == !!target_sample_size_per_arm,
      source_denominator_change_factor == !!source_denominator_change_factor,
      target_to_source_std_ratio == !!target_to_source_std_ratio |
        is.na(target_to_source_std_ratio)
    )
  results_df <- cbind(results_df, get_parameters(results_df[, "parameters"]))

  results_df <- merge(results_df, parameters_combinations[, ])


  if (nrow(results_df) == 0) {
    warning("Dataframe is empty")
    return()
  }

  xvar_name <- rlang::sym(xvar$name)
  xvar_label <- rlang::sym(xvar$label)

  if (metric %in% names(frequentist_metrics)) {
    selected_metric <- frequentist_metrics[[metric]]
  } else if (metric %in% names(inference_metrics)) {
    selected_metric <- inference_metrics[[metric]]
  } else if (metric %in% names(bayesian_metrics)) {
    selected_metric <- bayesian_metrics[[metric]]
  } else {
    stop("Metric is not supported")
  }

  if (is.null(selected_metric$name)) {
    stop("Metric name is NULL")
  }
  selected_metric_name <- rlang::sym(selected_metric$name)
  selected_metric_uncertainty_lower <- rlang::sym(paste0(selected_metric$metric_uncertainty, "_lower"))
  selected_metric_uncertainty_upper <- rlang::sym(paste0(selected_metric$metric_uncertainty, "_upper"))
  selected_label <- selected_metric$label

  parameters_str <- convert_params_to_str(methods_dict[[method]], parameters_combinations)

  figure_name <- paste0(
    case_study,
    "_",
    method,
    "_",
    selected_metric_name,
    "_vs_",
    xvar_name,
    "_target_sample_size_per_arm=",
    target_sample_size_per_arm,
    "_",
    parameters_str, "_comparisons")

  figure_name <- format_filename(filename = figure_name, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)

  directory <- paste0(figures_dir, case_study)
  if (!dir.exists(directory)) {
    dir.create(directory, showWarnings = FALSE, recursive = TRUE)
  }

  file_path <- file.path(directory, figure_name)

  if (file.exists(paste0(file_path, ".pdf")) && file.exists(paste0(file_path, ".png")) && remake_figures == FALSE){
    return()
  }

  # Calculate the range of the x-axis
  x_range <- range(results_df[, xvar$name])
  x_width <- diff(x_range)
  # Set a relative cap size (necessary, otherwise the cap size will depend on the x-axis width)

  cap_size <- x_width * relative_error_cap_width / 2 # Adjust relative_error_cap_width to control the relative cap size

  plt <- ggplot2::ggplot()
  plt <- plot_success_proba_vs_xvar(
    plt,
    results_df,
    xvar_name,
    selected_metric_name,
    selected_metric_uncertainty_lower,
    selected_metric_uncertainty_upper,
    cap_size,
    markersize,
    join_points,
    label = "Pr(Success) with borrowing"
  )

  color_groups <- c(
    "Pr(Success) with borrowing" = "black",
    "Nominal Pr(Success) (pooling)" = "red",
    "Freq. Pr(Success) at eq. TIE" = "green",
    "Nominal Pr(Success)" = "blue"
  )

  if ("at_equivalent_TIE" %in% baseline_success_proba) {
    plt <- plot_success_proba_vs_xvar(
      plt = plt,
      data = results_df,
      xvar_name = xvar_name,
      selected_metric_name = "frequentist_power_at_equivalent_tie",
      selected_metric_uncertainty_lower = "frequentist_power_at_equivalent_tie_lower",
      selected_metric_uncertainty_upper = "frequentist_power_at_equivalent_tie_upper",
      cap_size = cap_size,
      markersize = markersize,
      join_points = join_points,
      label = "Freq. Pr(Success) at eq. TIE"
    )
  }

  if ("at_nominal_TIE" %in% baseline_success_proba) {
    plt <- plot_baseline_success_proba_vs_xvar(
      plt,
      data = results_df,
      xvar_name = xvar_name,
      selected_metric_name = "nominal_frequentist_power_separate",
      cap_size = cap_size,
      markersize = markersize,
      join_points = join_points,
      label = "Nominal Pr(Success)"
    )
    # plt <- plot_baseline_success_proba_vs_xvar(plt,
    #                                            data = results_df,
    #                                            xvar_name = xvar_name,
    #                                            selected_metric_name = "nominal_frequentist_power_pooling",
    #                                            cap_size = cap_size,
    #                                            markersize= markersize,
    #                                            join_points = join_points,
    #                                            label = "Nominal Pr(Success) (pooling)")
  }

  drift_theta_0 <- theta_0 - unique(results_df$source_treatment_effect_estimate)
  plt <- plt +
    geom_vline(xintercept = drift_theta_0,
               linetype = "dashed",
               color = "black") #+
  # annotate("text", x = drift_theta_0, y = 0, label = TeX("$\\theta_0 - \\hat{\\theta}_0$"), vjust = -0.1 )

  param_label <- make_labels_from_parameters(parameters_combinations, method)

  if (param_label == "") {
    param_label <- paste0("")
  } else {
    param_label <- paste0(", ", param_label)
  }

  plot_title <- sprintf(
    "%s, %s%s, $N_T/2 = $ %s",
    str_to_title(case_study),
    methods_labels[[method]]$label,
    param_label,
    target_sample_size_per_arm
  )

  plot_title <- format_title(title = plot_title, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)

  # Add colors, labels and theme
  plt <- plt +
    ggplot2::labs(x = xvar_label, y = selected_label, title = plot_title) +
    theme_bw() +
    scale_color_manual(values = color_groups) +
    labs(color = "")+ # Legend title
    theme(legend.position = "bottom")

  plot.size <- set_size(textwidth)
  fig_width_in <- plot.size[1]
  fig_height_in <- plot.size[2]

  export_plots(
    plt = plt,
    file_path = file_path,
    fig_width_in,
    fig_height_in,
    type = "pdf"
  )
  export_plots(
    plt = plt,
    file_path = file_path,
    fig_width_in,
    fig_height_in,
    type = "png"
  )
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
plot_success_proba_vs_scenario <- function(results_metrics_df, metrics) {
  # Get the list of case studies
  case_studies <- unique(results_metrics_df$case_study)

  # Get the list of methods
  methods <- unique(results_metrics_df$method)

  # Remove non-zero control drift
  results_metrics_df <- results_metrics_df %>% dplyr::filter(control_drift == 0)


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

    theta_0 <- case_study_config$theta_0

    target_sample_sizes <- unique(filtered_results_metrics_df$target_sample_size_per_arm)

    for (method in methods) {
      if (sum(filtered_results_metrics_df$method == method) == 0) {
        next
      }
      # Get the different parameters combinations studies for this method
      parameters_combinations <- unique(get_parameters(filtered_results_metrics_df[filtered_results_metrics_df$method == method, "parameters"]))

      metric <- "success_proba"

      for (sample_size in target_sample_sizes) {
        for (i in 1:nrow(parameters_combinations)) {
          for (source_denominator_change_factor in unique(filtered_results_metrics_df$source_denominator_change_factor)) {
            for (target_to_source_std_ratio in unique(filtered_results_metrics_df$target_to_source_std_ratio)) {
              plot_success_proba_vs_drift(
                metric = metric,
                results_metrics_df = filtered_results_metrics_df,
                case_study = case_study,
                method = method,
                source_denominator_change_factor = source_denominator_change_factor,
                target_to_source_std_ratio = target_to_source_std_ratio,
                theta_0 = theta_0,
                target_sample_size_per_arm = sample_size,
                parameters_combinations = parameters_combinations[i, , drop = FALSE],
                xvars = xvars,
                baseline_success_proba = c("at_equivalent_TIE", "at_nominal_TIE"),
                join_points = TRUE
              )
            }
          }
        }
      }
    }
  }
}
