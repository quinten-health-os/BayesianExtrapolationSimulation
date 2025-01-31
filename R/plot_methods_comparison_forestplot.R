#' Create a forest plot
#'
#' @description This function creates a forest plot using the provided data.
#'
#' @param data The data for creating the forest plot.
#' @param title The title of the forest plot.
#' @param ylabel A logical value indicating whether to display the y-axis label.
#' @param x_metric_name Name of the metric on the x-axis
#' @param x_metric_uncertainty_lower Lower limit of the metric on the x-axis
#' @param x_metric_uncertainty_upper Upper limit of the metric on the x-axis
#' @param x_metric_label Label of the metric on the x-axis
#'
#' @return A ggplot2::ggplot( object representing the forest plot.
#' @keywords internal
forest_subplot <- function(data,
                           title,
                           ylabel,
                           x_metric_name,
                           x_metric_uncertainty_lower,
                           x_metric_uncertainty_upper,
                           x_metric_label,
                           methods_labels,
                           legend = FALSE,
                           sort_by = FALSE) {
  if (sort_by == "methods_parameters"){
    # Make sure that the data are sorted according to the parameters values
    for (method in unique(data$method)) {
      data_subset <- data[data$method == method, ]

      params <- get_parameters(data_subset[, "parameters"])

      column_names <- colnames(params)
      sorted_index <- do.call(order, lapply(params[column_names], as.factor))

      # Reorder the data :
      data[data$method == method, ] <- data_subset[sorted_index, ]
    }
  } else if (sort_by == "value"){
    x_metric_name_str <- as.character(x_metric_name)
    data <- data[order(data[[x_metric_name_str]]), ]
  }


  # Process the row of a results dataframe to create a Method + Parameters label
  labels <- sapply(1:nrow(data), function(i) {
    process_method_parameters_label(data[i, ], methods_labels)
  })

  data$rows <- seq(1, nrow(data))

  plt <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      y = rows,
      x = !!x_metric_name,
      xmin = !!x_metric_uncertainty_lower,
      xmax = !!x_metric_uncertainty_upper
    )
  ) +
    geom_pointrange(size = 0.001) +
    ggplot2::labs(title = title, x = x_metric_label) +
    theme(
      axis.title.y = element_blank(),
      text = element_text(family = font, size = small_text_size/2),
      axis.text.y = element_text(size = small_text_size/2),
      strip.text.x = element_text(size = small_text_size/2),
      legend.key = element_blank(),
      # strip.background = element_blank(),
      panel.grid.major = element_line(color = "gray80", linetype = "solid",  size = 0.15),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = small_text_size/2),
      plot.title = element_text(hjust = 0.5),
      plot.background = element_blank()
    ) +
    scale_y_continuous(breaks = data$rows, labels = labels) # LaTeX labels
    # scale_x_continuous(
    #   breaks = pretty(data[[x_metric_name]], n = 5), # Set 10 evenly spaced breaks
    #   labels = function(x) format(x, nsmall = 2)     # Optional: Format to 2 decimal places
    # )
  if (!ylabel) {
    plt <- plt + theme(axis.text.y = element_blank(), # Removes y-axis tick labels
                       axis.ticks.y = element_blank()) # Removes y-axis ticks
  }
  if (x_metric_name == "posterior_mean" ||
      x_metric_name == "posterior_median" || x_metric_name == "point_estimate_mean") {
    theta_0 <- unique(data$theta_0)
    if (legend) {
      plt <- plt + geom_vline(aes(xintercept = theta_0, color = "theta_0"), linetype = "dotted") +
        scale_color_manual(
          name = "",
          values = c("theta_0" = "black"),
          labels = expression(theta[0])
        )
    } else {
      plt <- plt + geom_vline(xintercept = theta_0, linetype = "dotted")
    }
  } else if (x_metric_name == "ess_moment" ||
             x_metric_name == "ess_precision" ||
             x_metric_name == "ess_elir") {
    if (is.null(data$equivalent_source_sample_size_per_arm) |
        any(is.na(data$equivalent_source_sample_size_per_arm))) {
      equivalent_source_sample_size_per_arm <- 2 * data$source_sample_size_control * data$source_sample_size_treatment / (data$source_sample_size_control + data$source_sample_size_treatment)
    } else {
      equivalent_source_sample_size_per_arm <- data$equivalent_source_sample_size_per_arm
    }
    target_sample_size_per_arm <- unique(data$target_sample_size_per_arm)
    source_sample_size_per_arm <- unique(round(equivalent_source_sample_size_per_arm))

    if (legend) {
      plt <- plt +
        geom_vline(
          aes(xintercept = target_sample_size_per_arm, color = "Target sample size per arm"),
          linetype = "dotted"
        ) +
        geom_vline(
          aes(xintercept = source_sample_size_per_arm, color = "Source sample size per arm"),
          linetype = "dotted"
        ) +
        scale_color_manual(
          name = "",
          values = c(
            "Target sample size per arm" = "black",
            "Source sample size per arm" = "red"
          ),
          labels = c(
            "Target sample size per arm" = expression(N[T]/2),
            "Source sample size per arm" = expression(N[S]/2)
          )
        )
    } else {
      plt <- plt +
        geom_vline(
          xintercept = target_sample_size_per_arm,
          linetype = "dotted",
          "color" = "black"
        ) +
        geom_vline(
          xintercept = source_sample_size_per_arm,
          linetype = "dotted",
          "color" = "red"
        )
    }
  }
  # } else if (x_metric_name == "success_proba" |
  #            x_metric_name == "average_tie" |
  #            x_metric_name == "average_power" |
  #            x_metric_name == "prior_proba_no_benefit" |
  #            x_metric_name == "prior_proba_benefit" |
  #            x_metric_name == "prepost_proba_FP" |
  #            x_metric_name == "prepost_proba_TP" |
  #            x_metric_name == "upper_bound_proba_FP" |
  #            x_metric_name == "prior_proba_success") {
  #   # # Make sure the x range are within (0,1)
  #   # plot_build <- ggplot_build(p)
  #   # default_x_limits <- plot_build$layout$panel_params[[1]]$x.range
  #   #
  #   # # Default x-axis limits
  #   # xl <- default_x_limits[1]
  #   # xu <- default_x_limits[2]
  #   #
  #   # # Compute new limits
  #   # new_xl <- max(xl, 0)
  #   # new_xu <- min(xu, 1)
  #
  #
  #   # plt <- plt + scale_x_continuous(limits = c(new_xl, new_xu))
  # }
  return(plt)
}


#' Create a forest plot
#'
#' @description This function creates a forest plot using the provided data.
#'
#' @param data The data for creating the forest plot.
#' @param title The title of the forest plot.
#' @param ylabel A logical value indicating whether to display the y-axis label.
#' @param x_metric_name Name of the metric on the x-axis
#' @param x_metric_uncertainty_lower Lower limit of the metric on the x-axis
#' @param x_metric_uncertainty_upper Upper limit of the metric on the x-axis
#' @param x_metric_label Label of the metric on the x-axis
#'
#' @return A ggplot2::ggplot( object representing the forest plot.
#' @keywords internal
forest_subplot_no_uncertainty <- function(data,
                                          title,
                                          ylabel,
                                          x_metric_name,
                                          x_metric_label,
                                          methods_labels,
                                          legend = FALSE,
                                          sort_by = FALSE) {

  if (sort_by == "methods_parameters"){
    # Make sure that the data are sorted according to the parameters values
    for (method in unique(data$method)) {
      data_subset <- data[data$method == method, ]

      params <- get_parameters(data_subset[, "parameters"])

      column_names <- colnames(params)
      sorted_index <- do.call(order, lapply(params[column_names], as.factor))

      # Reorder the data :
      data[data$method == method, ] <- data_subset[sorted_index, ]
    }
  } else if (sort_by == "value"){
    x_metric_name_str <- as.character(x_metric_name)
    data <- data[order(data[[x_metric_name_str]]), ]
  }

  labels <- sapply(1:nrow(data), function(i) {
    process_method_parameters_label(data[i, ], methods_labels)
  })

  data$rows <- seq(1, nrow(data))

  plt <- ggplot2::ggplot(data, ggplot2::aes(y = rows, x = !!x_metric_name)) +
    ggplot2::geom_point(size = 0.001) + # Scatter plot
    ggplot2::labs(title = title, x = x_metric_label) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      text = element_text(family = font, size = text_size/2),
      axis.text.y = element_text(size = text_size/2),
      strip.text.x = element_text(size = text_size/2),
      legend.key = element_blank(),
      # strip.background = element_blank(),
      panel.grid.major = element_line(color = "gray80", linetype = "solid",  size = 0.15),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = text_size/2),
      plot.title = element_text(hjust = 0.5)
      # Center the title
      # plot.background = element_blank()
    ) +
    scale_y_continuous(breaks = data$rows, labels = labels) # LaTeX labels

  if (!ylabel) {
    plt <- plt + theme(axis.text.y = element_blank(), # Removes y-axis tick labels
                       axis.ticks.y = element_blank()) # Removes y-axis ticks
  }

  return(plt)
}

forest_combined_plot <- function(data,
                                 x_metric_name,
                                 x_metric_uncertainty_lower,
                                 x_metric_uncertainty_upper,
                                 x_metric_label,
                                 methods_labels,
                                 legend = TRUE,
                                 sorting_by = "value") {

  # Sort data based on the specified sorting method
  if (sorting_by == "methods_parameters") {
    for (method in unique(data$method)) {
      data_subset <- data[data$method == method, ]
      params <- get_parameters(data_subset[, "parameters"])
      column_names <- colnames(params)
      sorted_index <- do.call(order, lapply(params[column_names], as.factor))
      data[data$method == method, ] <- data_subset[sorted_index, ]
    }
  } else if (sorting_by == "value") {
    x_metric_name_str <- as.character(x_metric_name)
    data <- data[order(data[[x_metric_name_str]]), ]
  }

  # Generate labels for the y-axis
  labels <- sapply(1:nrow(data), function(i) {
    process_method_parameters_label(data[i, ], methods_labels)
  })
  data$rows <- seq(1, nrow(data))

  # Create the plot
  plt <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      y = rows,
      x = !!x_metric_name,
      xmin = !!x_metric_uncertainty_lower,
      xmax = !!x_metric_uncertainty_upper,
      color = effect_type
    )
  ) +
    geom_pointrange(size = 0.001) +
    ggplot2::labs(
      title = "Add title",
      x = x_metric_label,
      color = ""
    ) +
    theme(
      axis.title.y = element_blank(),
      text = element_text(family = font, size = small_text_size/2),
      axis.text.y = element_text(size = small_text_size/2),
      strip.text.x = element_text(size = small_text_size/2),
      legend.key = element_blank(),
      # strip.background = element_blank(),
      panel.grid.major = element_line(color = "gray80", linetype = "solid",  size = 0.15),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = small_text_size/2),
      plot.title = element_text(hjust = 0.5)
      # plot.background = element_blank()
    ) +
    scale_y_continuous(breaks = data$rows, labels = labels) +
    scale_color_manual(values = c(
      "No effect" = "blue",
      "Partially consistent" = "orange",
      "Consistent (no drift)" = "green"
    ))

  # Add vertical lines for specific cases
  if (x_metric_name == "posterior_mean" || x_metric_name == "posterior_median" || x_metric_name == "point_estimate_mean") {
    theta_0 <- unique(data$theta_0)
    plt <- plt + geom_vline(xintercept = theta_0, linetype = "dotted", color = "black")
  } else if (x_metric_name %in% c("ess_moment", "ess_precision", "ess_elir")) {
    target_sample_size_per_arm <- unique(data$target_sample_size_per_arm)
    equivalent_source_sample_size_per_arm <- ifelse(
      is.null(data$equivalent_source_sample_size_per_arm) | any(is.na(data$equivalent_source_sample_size_per_arm)),
      2 * data$source_sample_size_control * data$source_sample_size_treatment / (data$source_sample_size_control + data$source_sample_size_treatment),
      data$equivalent_source_sample_size_per_arm
    )
    source_sample_size_per_arm <- unique(round(equivalent_source_sample_size_per_arm))

    plt <- plt +
      geom_vline(xintercept = target_sample_size_per_arm, linetype = "dotted", color = "black") +
      geom_vline(xintercept = source_sample_size_per_arm, linetype = "dotted", color = "red")
  }

  return(plt)
}


#' Generate a forest plot
#'
#' @description This function generates a forest plot based on the provided results dataframe, selected case study, and selected target sample size per arm.
#'
#' @param results_freq_df The results dataframe.
#' @param selected_case_study The selected case study.
#' @param selected_target_sample_size_per_arm The selected target sample size per arm.
#' @param x_metric Metric on the x-axis
#'
#' @return None
forest_plot <- function(results_freq_df, x_metric, panels = TRUE) {
  selected_case_study <- unique(results_freq_df$case_study)
  selected_target_sample_size_per_arm <- unique(results_freq_df$target_sample_size_per_arm)

  if (x_metric %in% names(frequentist_metrics)) {
    metric <- frequentist_metrics[[x_metric]]
  } else if (x_metric %in% names(inference_metrics)) {
    metric <- inference_metrics[[x_metric]]
  } else if (x_metric %in% names(bayesian_metrics)) {
    metric <- bayesian_metrics[[x_metric]]
  } else {
    stop("Metric is not supported")
  }

  if (!(metric$name %in% colnames(results_freq_df))) {
    warning("Metric not in input dataframe.")
    return()
  }

  x_metric_name <- rlang::sym(metric$name)
  x_metric_uncertainty_lower <- rlang::sym(paste0(metric$metric_uncertainty, "_lower"))
  x_metric_uncertainty_upper <- rlang::sym(paste0(metric$metric_uncertainty, "_upper"))
  x_metric_label <- metric$label


  directory <- file.path(figures_dir, selected_case_study)

  if (!dir.exists(directory)) {
    dir.create(directory, showWarnings = FALSE, recursive = TRUE)
  }

  filename <- paste0(
    selected_case_study,
    "_",
    x_metric,
    "_forest_plot_target_sample_size_per_arm_",
    selected_target_sample_size_per_arm
  )

  case_study = unique(results_freq_df$case_study)
  target_to_source_std_ratio = unique(results_freq_df$target_to_source_std_ratio)
  source_denominator_change_factor = unique(results_freq_df$source_denominator_change_factor)

  filename <- format_filename(
    filename = filename,
    case_study = case_study,
    target_to_source_std_ratio = target_to_source_std_ratio,
    source_denominator_change_factor = source_denominator_change_factor
  )

  file_path <- file.path(directory, filename)

  if (file.exists(paste0(file_path, ".pdf")) && file.exists(paste0(file_path, ".png")) && remake_figures == FALSE){
    return()
  }

  if (x_metric %in% names(bayesian_metrics)) {
    df <- results_freq_df %>%
      dplyr::mutate(means = !!x_metric_name)
  } else {
    df <- results_freq_df %>%
      dplyr::mutate(
        means = !!x_metric_name,
        error_low = !!x_metric_uncertainty_lower,
        error_upper = !!x_metric_uncertainty_upper
      )
  }

  round_decimal <- 4
  # Calculate important target treatment effects
  important_target_treatment_effects <- round(
    df %>%
      pull(source_treatment_effect_estimate) %>%
      unique() %>%
      first() %>%
      `*`(c(0, 1 / 2, 1)),
    round_decimal
  )

  # The following is to handle problems with the way R handles floating point precision.
  # Find the closest values in df$target_treatment_effect to important_target_treatment_effects
  closest_values <- sapply(important_target_treatment_effects, function(x) {
    df$target_treatment_effect[which.min(abs(df$target_treatment_effect - x))]
  })

  # Filter the dataframe based on approximate equality
  df <- df %>%
    rowwise() %>%
    filter(any(sapply(closest_values, function(val) is_approx_equal(target_treatment_effect, val)))) %>%
    ungroup() %>%
    dplyr::arrange(target_treatment_effect)

  # Custom function to find the closest value
  closest_value <- function(x, values) {
    values[which.min(abs(values - x))]
  }

  # Assign closest values to target_treatment_effect
  df <- df %>%
    mutate(target_treatment_effect= sapply(target_treatment_effect, closest_value, values = closest_values))


  # # Filter the dataframe to keep only the rows with the closest values
  # df <- df %>%
  #   dplyr::filter(target_treatment_effect %in% closest_values) %>%
  #   dplyr::arrange(target_treatment_effect)

  if (length(unique(df$target_treatment_effect)) != 3) {
    stop("Only three main treatment effect values should be selected")
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

  no_effect_data <- df[df[, "target_treatment_effect"] == "No treatment effect", ]
  partially_consistent_effect_data <- df[df[, "target_treatment_effect"] == "Partially consistent treatment effect", ]
  consistent_effect_data <- df[df[, "target_treatment_effect"] == "Consistent treatment effect", ]

  if (nrow(no_effect_data) != nrow(partially_consistent_effect_data) || nrow(consistent_effect_data) != nrow(partially_consistent_effect_data) || nrow(no_effect_data) != nrow(consistent_effect_data)){
    warning("Inconsistency in the number of rows in the partially consistent, consistent, and no-effect data. Skipping this plot.")
    return()
  }

  # Ensure that the methods and parameters are correctly ordered
  no_effect_data <- no_effect_data %>%
    arrange(method, parameters)

  partially_consistent_effect_data <- partially_consistent_effect_data %>%
    arrange(method, parameters)

  consistent_effect_data <- consistent_effect_data %>%
    arrange(method, parameters)

  # Sort the dataframes by the value of interest
  x_metric_name_str <- as.character(x_metric_name)
  new_order <- order(no_effect_data[[x_metric_name_str]])
  no_effect_data <- no_effect_data[new_order, ]
  partially_consistent_effect_data <- partially_consistent_effect_data[new_order, ]
  consistent_effect_data <- consistent_effect_data[new_order, ]


  if (x_metric_name == "success_proba"){
    x_metric_label_no_effect <- "TIE"
    x_metric_label_partially_consistent <- "Power"
    x_metric_label_consistent <- "Power"
  } else {
    x_metric_label_no_effect <- x_metric_label
    x_metric_label_partially_consistent <- x_metric_label
    x_metric_label_consistent <- x_metric_label
  }


  if (panels == TRUE){
    # Create forest plots for each effect
    plot_no_effect <- forest_subplot(
      no_effect_data,
      "No effect",
      ylabel = TRUE,
      x_metric_name,
      x_metric_uncertainty_lower,
      x_metric_uncertainty_upper,
      x_metric_label_no_effect,
      methods_labels = methods_labels,
      legend = FALSE,
      sort_by = FALSE
    )
    plot_partially_consistent_effect <- forest_subplot(
      partially_consistent_effect_data,
      "Partially consistent",
      ylabel = FALSE,
      x_metric_name,
      x_metric_uncertainty_lower,
      x_metric_uncertainty_upper,
      x_metric_label_partially_consistent,
      methods_labels = methods_labels,
      legend = FALSE,
      sort_by = FALSE
    )
    plot_consistent_effect <- forest_subplot(
      consistent_effect_data,
      "Consistent (no drift)",
      ylabel = FALSE,
      x_metric_name,
      x_metric_uncertainty_lower,
      x_metric_uncertainty_upper,
      x_metric_label_consistent,
      methods_labels = methods_labels,
      legend = TRUE,
      sort_by = FALSE
    )

    # Convert to grobs
    grob_no_effect <- ggplotGrob(plot_no_effect)
    grob_partially_consistent_effect <- ggplotGrob(plot_partially_consistent_effect)
    grob_consistent_effect <- ggplotGrob(plot_consistent_effect)

    grid.newpage()

    # Note that adding the Tikz export to the export_plots function does not work.
    if (plots_to_latex == TRUE) {
      tikzDevice::tikz(file = paste0(file_path, ".tex"), width = 5.92)
    }

    plt <- gridExtra::grid.arrange(
      grob_no_effect,
      grob_partially_consistent_effect,
      grob_consistent_effect,
      ncol = 3,
      widths = c(0.44, 0.25, 0.31)
    )
    if (plots_to_latex == TRUE) {
      dev.off()
    }
  } else {
    # Combine data from all three cases into one data frame
    combined_data <- rbind(
      cbind(no_effect_data, effect_type = "No effect"),
      cbind(partially_consistent_effect_data, effect_type = "Partially consistent"),
      cbind(consistent_effect_data, effect_type = "Consistent (no drift)")
    )

    plt <- forest_combined_plot(
      combined_data,
      x_metric_name = "posterior_mean",
      x_metric_uncertainty_lower = "lower_bound",
      x_metric_uncertainty_upper = "upper_bound",
      x_metric_label = "Posterior Mean",
      methods_labels = methods_labels,
      legend = TRUE
    )
  }
  plot_size <- set_size(textwidth)
  fig_width_in <- plot_size[1]
  fig_height_in <- plot_size[2]

  export_plots(plt, file_path, fig_width_in, fig_height_in, type = "pdf", forest_plot = TRUE)
  export_plots(plt, file_path, fig_width_in, fig_height_in, type = "png", forest_plot = TRUE)
}


#' Generate a forest plot for bayesian metrics
#'
#' @description This function generates a forest plot based on the provided results dataframe, selected case study, and selected target sample size per arm.
#'
#' @param results_bayes_df The results dataframe.
#' @param selected_case_study The selected case study.
#' @param selected_target_sample_size_per_arm The selected target sample size per arm.
#' @param x_metric Metric on the x-axis
#'
#' @return None
forest_plot_bayesian <- function(results_bayes_df, x_metric) {
  selected_case_study <- unique(results_bayes_df$case_study)
  selected_target_sample_size_per_arm <- unique(results_bayes_df$target_sample_size_per_arm)


  directory <- file.path(figures_dir, selected_case_study)
  if (!dir.exists(directory)) {
    dir.create(directory, showWarnings = FALSE, recursive = TRUE)
  }

  if (x_metric %in% names(bayesian_metrics)) {
    metric <- bayesian_metrics[[x_metric]]
  } else {
    stop("Metric is not supported")
  }

  if (!(x_metric %in% colnames(results_bayes_df))) {
    warning("Metric not in input dataframe.")
    return()
  }

  x_metric_name <- rlang::sym(metric$name)
  x_metric_label <- metric$label


  filename <- paste0(
    selected_case_study,
    "_",
    x_metric,
    "_forest_plot_target_sample_size_per_arm_",
    selected_target_sample_size_per_arm
  )

  file_path <- file.path(directory, filename)

  if (file.exists(paste0(file_path, ".pdf")) && file.exists(paste0(file_path, ".png")) && remake_figures == FALSE){
    return()
  }

  round_decimal <- 4

  df <- results_bayes_df %>%
    dplyr::mutate(means = !!x_metric_name)

  ui_design_prior_data <- df[df[, "design_prior_type"] == "ui_design_prior", ]
  analysis_design_prior_data <- df[df[, "design_prior_type"] == "analysis_prior", ]
  source_posterior_design_prior_data <- df[df[, "design_prior_type"] == "source_posterior", ]


  # Sort the dataframes by the value of interest
  x_metric_name_str <- as.character(x_metric_name)
  new_order <- order(ui_design_prior_data[[x_metric_name_str]])
  ui_design_prior_data <- ui_design_prior_data[new_order, ]
  analysis_design_prior_data <- analysis_design_prior_data[new_order, ]
  source_posterior_design_prior_data <- source_posterior_design_prior_data[new_order, ]

  # Create forest plots for each effect
  plot_ui_design_prior <- forest_subplot_no_uncertainty(
    ui_design_prior_data,
    "UI design prior",
    ylabel = TRUE,
    x_metric_name,
    x_metric_label,
    methods_labels = methods_labels,
    legend = FALSE
  )
  plot_analysis_design_prior <- forest_subplot_no_uncertainty(
    analysis_design_prior_data,
    "Analysis prior",
    ylabel = FALSE,
    x_metric_name,
    x_metric_label,
    methods_labels = methods_labels,
    legend = FALSE
  )
  plot_source_posterior_design_prior <- forest_subplot_no_uncertainty(
    source_posterior_design_prior_data,
    "Source posterior as design prior",
    ylabel = FALSE,
    x_metric_name,
    x_metric_label,
    methods_labels = methods_labels,
    legend = FALSE
  )

  # Convert to grobs
  grob_ui_design_prior <- ggplotGrob(plot_ui_design_prior)
  grob_analysis_design_prior <- ggplotGrob(plot_analysis_design_prior)
  grob_source_posterior_design_prior <- ggplotGrob(plot_source_posterior_design_prior)

  grid.newpage()
  plt <- gridExtra::grid.arrange(
    grob_ui_design_prior,
    grob_analysis_design_prior,
    grob_source_posterior_design_prior,
    ncol = 3,
    widths = c(0.46, 0.23, 0.31)
  )

  plot.size <- set_size(textwidth)
  fig_width_in <- plot.size[1]
  fig_height_in <- plot.size[2]

  export_plots(plt, file_path, fig_width_in, fig_height_in, type = "pdf", forest_plot = TRUE)
  export_plots(plt, file_path, fig_width_in, fig_height_in, type = "png", forest_plot = TRUE)
}


#' Plot methods comparison
#'
#' @description This function plots the comparison of methods based on the provided results metrics dataframe.
#'
#' @param results_metrics_df The results metrics dataframe.
#' @param metrics Metrics for which to create forest plots.
#'
#' @return None
forest_plot_methods_comparison <- function(results_metrics_df, metrics) {
  case_studies <- unique(results_metrics_df$case_study)

  results_metrics_df <- results_metrics_df[(results_metrics_df[, "source_denominator_change_factor"] == 1 &
                                              !is.na(results_metrics_df[, "source_denominator_change_factor"])), ]

  for (metric in names(metrics)) {
    if (metric %in% c('tie', 'power')){
      next
    }

    for (case_study in case_studies) {
      results_metrics_df1 <- results_metrics_df[results_metrics_df$case_study == case_study,]

      target_sample_sizes <- unique(results_metrics_df1$target_sample_size_per_arm)
      target_sample_sizes <- target_sample_sizes[!target_sample_sizes %in% NA]
      for (target_sample_size_per_arm in target_sample_sizes) {
        results_metrics_df2 <- results_metrics_df1[results_metrics_df1$target_sample_size_per_arm == target_sample_size_per_arm,]

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

            forest_plot(results_metrics_df4, metric)
          }
        }
      }
    }
  }
}

#' Plot methods comparison
#'
#' @description This function plots the comparison of methods based on the provided results metrics dataframe.
#'
#' @param results_metrics_df The results metrics dataframe.
#' @param metrics Metrics for which to create forest plots.
#'
#' @return None
forest_plot_methods_comparison_bayesian <- function(results_metrics_df, metrics) {
  case_studies <- unique(results_metrics_df$case_study)

  results_metrics_df <- results_metrics_df[(results_metrics_df[, "source_denominator_change_factor"] == 1 &
                                              !is.na(results_metrics_df[, "source_denominator_change_factor"])), ]

  for (metric in names(metrics)) {
    for (case_study in case_studies) {
      results_metrics_df1 <- results_metrics_df[results_metrics_df$case_study == case_study,]

      target_sample_sizes <- unique(results_metrics_df1$target_sample_size_per_arm)
      target_sample_sizes <- target_sample_sizes[!target_sample_sizes %in% NA]
      for (target_sample_size_per_arm in target_sample_sizes) {
        results_metrics_df2 <- results_metrics_df1[results_metrics_df1$target_sample_size_per_arm == target_sample_size_per_arm,]

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
            forest_plot_bayesian(results_metrics_df4, metric)
          }
        }
      }
    }
  }
}
