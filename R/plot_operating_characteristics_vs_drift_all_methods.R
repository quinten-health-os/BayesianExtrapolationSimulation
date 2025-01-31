plot_element_metric_vs_xvar_methods <- function(plt,
                                        data,
                                        xvar_name,
                                        selected_metric_name,
                                        selected_metric_uncertainty_lower,
                                        selected_metric_uncertainty_upper,
                                        cap_size,
                                        markersize,
                                        labels_title,
                                        join_points,
                                        color = NULL,
                                        is_baseline = FALSE,
                                        dodging = TRUE,
                                        parameters_labels = NULL) {
  if (nrow(data) == 0) {
    stop("The dataframe is empty.")
  }


  parameters_labels <- setNames(object = parameters_labels, nm = data$parameters_labels)

  if (dodging){
    if (!is.null(selected_metric_uncertainty_lower)){
      plt <- plt +
        geom_errorbar(
          data = data,
          ggplot2::aes(
            x = !!xvar_name,
            y = !!selected_metric_name,
            ymin = !!selected_metric_uncertainty_lower,
            ymax = !!selected_metric_uncertainty_upper,
            color = factor(parameters_labels, levels = unique(data$parameters_labels)),
            # Used to keep the parameters ordered
            group = factor(parameters_labels, levels = unique(data$parameters_labels)),
            shape = factor(parameters_labels, levels = unique(data$parameters_labels))
          ),
          position = position_dodge(width = cap_size),
          width = cap_size
        )
    }
    plt <- plt +
      geom_point(
        data = data,
        ggplot2::aes(
          x = !!xvar_name,
          y = !!selected_metric_name,
          color = factor(parameters_labels, levels = unique(data$parameters_labels)),
          # Used to keep the parameters ordered
          group = factor(parameters_labels, levels = unique(data$parameters_labels)),
          shape = factor(parameters_labels, levels = unique(data$parameters_labels))
        ),
        position = position_dodge(width = cap_size),
        size = markersize,
      )
    } else {
      if (!is.null(selected_metric_uncertainty_lower)){
        plt <- plt +
          geom_errorbar(
            data = data,
            ggplot2::aes(
              x = !!xvar_name,
              y = !!selected_metric_name,
              ymin = !!selected_metric_uncertainty_lower,
              ymax = !!selected_metric_uncertainty_upper,
              color = factor(parameters_labels, levels = unique(data$parameters_labels)),
              # Used to keep the parameters ordered
              group = factor(parameters_labels, levels = unique(data$parameters_labels)),
              shape = factor(parameters_labels, levels = unique(data$parameters_labels))
            ),
            width = cap_size
          )
      }
      plt <- plt +
        geom_point(
          data = data,
          ggplot2::aes(
            x = !!xvar_name,
            y = !!selected_metric_name,
            color = factor(parameters_labels, levels = unique(data$parameters_labels)),
            # Used to keep the parameters ordered
            group = factor(parameters_labels, levels = unique(data$parameters_labels)),
            shape = factor(parameters_labels, levels = unique(data$parameters_labels))
          ),
          size = markersize
        )
  }

  plt <- plt + labs(color = labels_title) + scale_color_viridis_d(labels = parameters_labels, guide = "legend") + scale_fill_viridis_d(labels = parameters_labels, guide = "legend") +
    scale_shape_manual(values = 1:length(parameters_labels), labels = parameters_labels, guide = "legend") +
    guides(color = guide_legend(title = "Methods"), shape = guide_legend(title = "Methods"))


  if (join_points == TRUE) {
    if (dodging){
      plt <- plt + geom_line(
        #  Add lines between points
        data = data,
        ggplot2::aes(
          x = !!xvar_name,
          y = !!selected_metric_name,
          color = factor(parameters_labels, levels = unique(data$parameters_labels)),
          # Used to keep the parameters ordered
          group = factor(parameters_labels, levels = unique(data$parameters_labels)),
        ),
        position = position_dodge(width = cap_size)
      )
    } else {
      plt <- plt + geom_line(
        #  Add lines between points
        data = data,
        ggplot2::aes(
          x = !!xvar_name,
          y = !!selected_metric_name,
          color = factor(parameters_labels, levels = unique(data$parameters_labels)),
          group = factor(parameters_labels, levels = unique(data$parameters_labels)),
        )
      )
    }
  }

  return(plt)
}


#' Function to generate a metric vs drift plot
#'
#' @param metric The metric to plot
#' @param results_metrics_df The dataframe containing the results and metrics
#' @param theta_0 The true treatment effect
#' @param case_study The case study name
#' @param method The method name
#' @param category The category to group the results by (either "parameters" or "target_sample_size_per_arm")
#' @param target_sample_size_per_arm The target sample size per arm
#' @param parameters_combinations The combinations of parameters to filter the results by
#' @param xvars A list of x-variables for control drift and treatment drift
#' @param target_to_source_std_ratio Ratio between the target and source studies sampling standard deviation.
#' @param join_points Whether to join the point with a line or not.
#'
#' @return None
#'
#' @export
plot_metric_vs_drift_methods <- function(metric,
                                 results_metrics_df,
                                 theta_0,
                                 case_study,
                                 target_sample_size_per_arm,
                                 xvars,
                                 target_to_source_std_ratio = 1,
                                 join_points = FALSE,
                                 source_denominator_change_factor = 1) {

  if (!(metric %in% colnames(results_metrics_df))) {
    warning("Metric not in input dataframe.")
    return()
  }

  # Remove non-zero control drift
  results_df <- results_metrics_df %>% dplyr::filter(control_drift == 0)

  results_df <- results_df %>%
    dplyr::arrange(method, drift)

  xvar <- xvars$drift

  boolean_filter <- c()
  # Sort by method before filtering on parameters
  results_df <- results_df %>% arrange(method)
  # Filter on the important parameters values for each method
  for (method in unique(results_df$method)){

    parameters <- get_parameters(results_df[results_df$method == method, "parameters", drop = FALSE])
    filter_per_method = rep(TRUE, nrow(parameters))
    for (key in names(methods_dict[[method]])) {
      # Loop over the method's parameters

      if (method == "commensurate_power_prior"){
        filter <- (parameters$heterogeneity_prior.family == "half_normal" &  parameters$heterogeneity_prior.std_dev == 1) | (parameters$heterogeneity_prior.family == "inverse_gamma" &  parameters$heterogeneity_prior.alpha == 0.3333)
      } else {
        # Get the important values for which we will make plots
        if (is.null(methods_dict[[method]][[key]][["important_values"]])) {
          important_parameters_values <- methods_dict[[method]][[key]][["range"]]
        } else {
          important_parameters_values <- methods_dict[[method]][[key]][["important_values"]]
        }

        filter <- parameters[,key] %in% important_parameters_values
      }
      filter_per_method <- filter & filter_per_method
    }
    boolean_filter <- c(boolean_filter, filter_per_method)
  }

  results_df <- results_df[boolean_filter,]

  color_map <- viridis::viridis(n = 256)

  results_df$rows <- seq(1, nrow(results_df))

  # Process the row of a results dataframe to create a Method + Parameters label
  parameters_labels <- sapply(1:nrow(results_df), function(i) {
    process_method_parameters_label(results_df[i, ], methods_labels)
  })

  results_df$parameters_labels <- lapply(1:nrow(results_df), function(i) {
    process_method_parameters_label(results_df[i, ], methods_labels)
  })

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

  # Calculate the range of the x-axis
  x_range <- range(results_df[, xvar$name])
  x_width <- diff(x_range)
  # Set a relative cap size (necessary, otherwise the cap size will depend on the x-axis width)

  cap_size <- x_width * relative_error_cap_width # Adjust relative_error_cap_width to control the relative cap size

  labels_title <- "Methods"
  plt <- ggplot2::ggplot()
  plt <- plot_element_metric_vs_xvar_methods(
    plt,
    results_df,
    xvar_name,
    selected_metric_name,
    selected_metric_uncertainty_lower,
    selected_metric_uncertainty_upper,
    cap_size,
    markersize,
    labels_title,
    join_points,
    color,
    dodging = FALSE,
    parameters_labels =  parameters_labels
  )

  if(metric %in% c("tie", "success_proba")){
    plt <- plt +
      ggplot2::geom_hline(yintercept = analysis_config$nominal_tie, color = "black", linetype = "dashed") + # Add horizontal line at y = 0.05
      ggplot2::scale_y_continuous(breaks = function(x) unique(c(pretty(x), analysis_config$nominal_tie)))  # Add y-tick at y = 0.05
  }

  case_study_config <- yaml::read_yaml(paste0(case_studies_config_dir, case_study, ".yml"))
  source_treatment_effect_estimate <- unique(results_df$source_treatment_effect_estimate)[1]
  mandatory_drift_values <- sort(important_drift_values(source_treatment_effect_estimate, case_study_config))

  mandatory_drift_values <- data.frame(
    drift_value = mandatory_drift_values,
    effect_label = c("No effect", "Partially consistent effect", "Consistent effect")
  )

  plt <- plt +
    geom_vline(data = mandatory_drift_values,
               aes(xintercept = drift_value, linetype = effect_label), color = "black") +
    scale_linetype_manual(values = c("No effect" = "dotted",
                                   "Partially consistent effect" = "dashed",
                                   "Consistent effect" = "solid")) +
    labs(linetype = "Treatment effect")

  title <- sprintf(
          "%s, $N_T/2 = %s$",
          str_to_title(case_study),
          target_sample_size_per_arm
        )

  figure_name <- paste0(
    case_study,
    "_",
    selected_metric_name,
    "_vs_",
    xvar_name,
    "_all_methods_target_sample_size_",
    target_sample_size_per_arm
  )

  title <- format_title(title = title, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)
  figure_name <- format_filename(filename = figure_name, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)


  # Add colors, labels and theme
  plt <- plt +
    ggplot2::labs(x = xvar_label, y = selected_label, title = title)

  # Save plot
  directory <- paste0(figures_dir, case_study)
  if (!dir.exists(directory)) {
    dir.create(directory, showWarnings = FALSE, recursive = TRUE)
  }

  plot.size <- set_size(426)
  fig_width_in <- plot.size[1]
  fig_height_in <- plot.size[2]*1.15

  file_path <- file.path(directory, figure_name)

  if (file.exists(paste0(file_path, ".pdf")) && file.exists(paste0(file_path, ".png")) && remake_figures == FALSE){
    return()
  }
  export_plots(plt, file_path, fig_width_in, fig_height_in, type = "pdf")
  export_plots(plt, file_path, fig_width_in, fig_height_in, type = "png")


}

#' Function to plot metric vs sample size
#'
#' @description This function plots a metric against the sample size for a given method, case study, and parameters combinations.
#'
#' @param metric The metric to plot.
#' @param results_metrics_df The data frame containing the results and metrics.
#' @param case_study The case study to plot the metric for.
#' @param theta_0 Boundary of the null hypothesis space.
#'
#' @return None
#'
#' @export
plot_metric_vs_sample_size_methods <- function(metric,
                                       results_metrics_df,
                                       case_study,
                                       theta_0 = 0,
                                       target_to_source_std_ratio = 1,
                                       source_denominator_change_factor = 1,
                                       dodging = FALSE,
                                       join_points = FALSE,
                                       target_treatment_effect = target_treatment_effect) {


  if (!(metric %in% colnames(results_metrics_df))) {
    warning("Metric not in input dataframe.")
    return()
  }
  # Remove non-zero control drift
  results_df <- results_metrics_df %>% dplyr::filter(control_drift == 0)

  # Filter results_metrics_df on the selected case study and method
  results_df <- results_metrics_df %>%
    dplyr::filter(
      case_study == !!case_study,
      source_denominator_change_factor == !!source_denominator_change_factor |
        is.na(source_denominator_change_factor),
      target_to_source_std_ratio == !!target_to_source_std_ratio |
        is.na(target_to_source_std_ratio),
    )

  if (nrow(results_df) == 0) {
    warning("Dataframe is empty")
    return()
  }



  case_study_config <- yaml::read_yaml(paste0(case_studies_config_dir, case_study, ".yml"))

  # Filter on important drift value
  source_treatment_effect_estimate <- unique(results_df$source_treatment_effect_estimate)[1]
  mandatory_drift_values <- important_drift_values(source_treatment_effect_estimate, case_study_config)

  if ("drift" %in% colnames(results_df)) {
    # Find the closest values in df$target_treatment_effect to important_target_treatment_effects
    closest_values <- sapply(mandatory_drift_values, function(x) {
      results_df$drift[which.min(abs(results_df$drift - x))]
    })

    # Filter the dataframe to keep only the rows with the closest values
    results_df <- results_df %>%
      dplyr::filter(drift %in% closest_values) %>%
      dplyr::arrange(drift)

    if (length(unique(results_df$drift)) != 3) {
      stop("Only three main treatment effect values should be selected")
    }


    closest_values <- sort(closest_values)
    results_df <- results_df %>%
      mutate(effect_name = case_when(
        drift == closest_values[1] ~ "No effect",
        drift == closest_values[2] ~ "Partially consistent effect",
        drift == closest_values[3] ~ "Consistent effect"
      ))
    # Filter on the target_treatmnent_effect
    results_df <- results_df %>%
      dplyr::filter(effect_name == !! target_treatment_effect)
  }

  if (nrow(results_df) == 0) {
    warning("Dataframe is empty")
    return()
  }


  boolean_filter <- c()
  # Sort by method before filtering on parameters
  results_df <- results_df %>% arrange(method)
  # Filter on the important parameters values for each method
  for (method in unique(results_df$method)){

    parameters <- get_parameters(results_df[results_df$method == method, "parameters", drop = FALSE])
    filter_per_method = rep(TRUE, nrow(parameters))
    for (key in names(methods_dict[[method]])) {
      # Loop over the method's parameters

      if (method == "commensurate_power_prior"){
        filter <- (parameters$heterogeneity_prior.family == "half_normal" &  parameters$heterogeneity_prior.std_dev == 1) | (parameters$heterogeneity_prior.family == "inverse_gamma" &  parameters$heterogeneity_prior.alpha == 0.3333)
      } else {
        # Get the important values for which we will make plots
        if (is.null(methods_dict[[method]][[key]][["important_values"]])) {
          important_parameters_values <- methods_dict[[method]][[key]][["range"]]
        } else {
          important_parameters_values <- methods_dict[[method]][[key]][["important_values"]]
        }

        filter <- parameters[,key] %in% important_parameters_values
      }
      filter_per_method <- filter & filter_per_method
    }
    boolean_filter <- c(boolean_filter, filter_per_method)
  }

  results_df <- results_df[boolean_filter,]


  # Process the row of a results dataframe to create a Method + Parameters label
  parameters_labels <- sapply(1:nrow(results_df), function(i) {
    process_method_parameters_label(results_df[i, ], methods_labels, method_name = FALSE)
  })

  results_df$parameters_labels <- lapply(1:nrow(results_df), function(i) {
    process_method_parameters_label(results_df[i, ], methods_labels, method_name = FALSE)
  })

  results_df$methods_parameters_labels <- lapply(1:nrow(results_df), function(i) {
    process_method_parameters_label(results_df[i, ], methods_labels, method_name = TRUE, as_latex = FALSE)
  })

  parameters_labels <- setNames(object = parameters_labels, nm = results_df$parameters_labels)

  methods <- format_results_df_methods(results_df)
  results_df$method <- methods

  # Extract the metric of interest
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

  # Calculate the range of the x-axis
  x_range <- range(results_df$target_sample_size_per_arm)
  x_width <- diff(x_range)
  # Set a relative cap size (necessary, otherwise the cap size will depend on the x-axis width)
  cap_size <- x_width * relative_error_cap_width # Adjust relative_error_cap_width to control the relative cap size

  # Calculate minimum difference between x-values
  x_values <- sort(unique(results_df$target_sample_size_per_arm))
  if (length(x_values) > 1) {
    x_diffs <- diff(x_values)
    min_x_diff <- min(x_diffs)
  } else {
    min_x_diff <- 1
  }

  # Determine number of groups
  n_groups <- length(unique(results_df$parameters_labels))

  # Set fractions for cap width and dodge width
  cap_width_fraction <- 0.5
  dodge_width_fraction <- 0.8

  # Calculate cap_width and dodge_width
  cap_width <- min_x_diff * cap_width_fraction
  total_dodge_width <- min_x_diff * dodge_width_fraction
  dodge_width <- total_dodge_width / (n_groups - 1)

  results_df <- results_df %>%
    group_by(parameters_labels) %>%
    mutate(
      jitter_offset = jitter(0, amount = 4)[1],  # Generate a single jitter offset per group
      jittered_x = target_sample_size_per_arm + jitter_offset  # Add the same offset to all rows in the group
    ) %>%
    ungroup()

  # Split data by 'method'
  results_df_split <- split(results_df, results_df$method)

  # Unique parameters and methods
  unique_parameters <- unique(results_df$parameters_labels)
  unique_methods <- sort(unique(results_df$method))

  # Shapes for methods
  shape_codes <- c(16, 2, 15, 18, 1, 8, 4, 3, 7, 9)

  shapes <- setNames(shape_codes[1:length(unique_methods)], unique_methods)


  # Initialize the plot
  plt <- ggplot(mapping = aes(x = target_sample_size_per_arm, y = !!selected_metric_name))

  # Iterate over each method and add layers
  plt <- plt + purrr::imap(results_df_split, function(df_subset, method_name) {
    df_subset <- df_subset[!is.na(df_subset[[selected_metric_name]]), ]


    # Get unique parameters for this method
    method_parameters <- unique(df_subset$parameters_labels)

    if (length(method_parameters) == 0){
      return()
    }

    # Generate a large color palette
    full_palette <- scales::hue_pal()(50)  # Generate a larger palette
    # Randomly sample 'length(method_parameters)' colors from the full palette
    method_colors <- setNames(
      sample(full_palette, length(method_parameters), replace = FALSE),
      method_parameters
    )

    # Get the shape for this method
    method_shape <- shapes[method_name]

    # Create the plot layers for this method
    layers <- list(
      # Add points for this method
      geom_point(
        data = df_subset,
        aes(x = jittered_x, color = factor(parameters_labels, levels = method_parameters)),
        shape = method_shape,
        size = 2
      ),
      geom_line(
        data = df_subset,
        aes(x = jittered_x, color = factor(parameters_labels, levels = method_parameters))
      ),
      # Add vertical error bars
      geom_errorbar(
        data = df_subset,
        aes(
          x = jittered_x,
          ymin = !!sym(paste0("conf_int_", selected_metric_name, "_lower")),
          ymax = !!sym(paste0("conf_int_", selected_metric_name, "_upper")),
          color = factor(parameters_labels, levels = method_parameters)
        ),
        width = cap_width
      ),
      # Define the color scale for parameters within this method
      scale_color_manual(
        values = method_colors,
        labels = method_parameters,
        name = method_name,
        drop = TRUE,
        guide = guide_legend(override.aes = list(shape = method_shape), ncol = 1)
      )
    )


    # Reset the color scale for the next method
    layers <- c(layers, list(new_scale_color()))

    layers
  })

  if(metric %in% c("tie", "success_proba")){
    plt <- plt +
      ggplot2::geom_hline(yintercept = analysis_config$nominal_tie, color = "black", linetype = "dashed") + # Add horizontal line at y = 0.05
      ggplot2::scale_y_continuous(breaks = function(x) unique(c(pretty(x), analysis_config$nominal_tie)))  # Add y-tick at y = 0.05
  }

  plot_title = sprintf(
    "%s, %s",
    str_to_title(case_study),
    target_treatment_effect
  )
  plot_title <- format_title(title = plot_title, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)

  # Add shape scale for methods and labels
  plt <- plt +
    scale_shape_manual(
      values = shapes,
      name = "Methods",
      guide = guide_legend(order = 1, nrow = 2)
    ) +
    labs(
      title = plot_title,
      x = "Target study sample size per arm",
      y = selected_label
    )

  # Apply theme settings
  plt <- plt + theme_bw() + theme(
    axis.text = element_text(family = font, size = text_size / 2),
    axis.text.y = element_text(family = font, size = small_text_size / 2),
    axis.text.x = element_text(family = font, size = small_text_size / 2),
    axis.title = element_text(family = font, size = text_size / 2),
    plot.title = element_text(family = font, size = text_size / 2),
    legend.text = element_text(family = font, size = small_text_size / 2),
    legend.title = element_text(family = font, size = text_size / 2),
    legend.key.size = unit(0.1, "cm"),  # Key size
    legend.spacing.x = unit(0.03, "cm"),
    legend.spacing.y = unit(0.01, "cm"),  # Narrow vertical spacing
    legend.position = "bottom",
    legend.direction = "vertical",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.box = "horizontal",             # Stack legends vertically
    legend.box.just = "center"             # Align legends to the left
  )

  plot.size <- set_size(textwidth)
  fig_width_in <- plot.size[1]*1.35
  fig_height_in <- plot.size[2]*1.5

  # plt <- plt +
  #   ggplot2::labs(title = title, x = "Target study sample size per arm",
  #                 y = selected_label)

  directory <- paste0(figures_dir, case_study)
  if (!dir.exists(directory)) {
    dir.create(directory, showWarnings = FALSE, recursive = TRUE)
  }

  figure_name <- paste0(
    case_study,
    "_", target_treatment_effect,
    "_",
    selected_metric_name,
    "_vs_target_sample_size_per_arm_all_methods"
  )

  figure_name <- format_filename(filename = figure_name, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)


  # plot.size <- set_size(426)
  # fig_width_in <- plot.size[1]
  # fig_height_in <- plot.size[2]

  file_path <- file.path(directory, figure_name)

  if (file.exists(paste0(file_path, ".pdf")) && file.exists(paste0(file_path, ".png")) && remake_figures == FALSE){
    return()
  }
  export_plots(plt, file_path, fig_width_in, fig_height_in, type = "pdf", adjust_theme = FALSE)
  export_plots(plt, file_path, fig_width_in, fig_height_in, type = "png", adjust_theme = FALSE)

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
plot_metric_vs_scenario_methods <- function(results_metrics_df, metrics) {
  # Get the list of case studies
  case_studies <- unique(results_metrics_df$case_study)

  treatment_effects_labels <- c("No effect", "Partially consistent effect", "Consistent effect")


  for (case_study in case_studies) {
    case_study_config <- yaml::yaml.load_file(system.file(
      paste0("conf/case_studies/", case_study, ".yml"),
      package = "RBExT"
    ))

    results_metrics_df1 <- results_metrics_df %>%
      dplyr::filter(case_study == !!case_study)

    if (nrow(results_metrics_df1) == 0) {
      warning("Dataframe is empty")
      next
    }

    theta_0 <- case_study_config$theta_0

    target_sample_sizes <- unique(results_metrics_df1$target_sample_size_per_arm)

    for (metric in names(metrics)) {
      if (!(metric %in% c("tie"))){
        for (source_denominator_change_factor in unique(results_metrics_df1$source_denominator_change_factor)) {
          results_metrics_df2 <- results_metrics_df1 %>%
            dplyr::filter(source_denominator_change_factor == !!source_denominator_change_factor | is.na(source_denominator_change_factor))

          for (target_to_source_std_ratio in unique(results_metrics_df2$target_to_source_std_ratio)) {
            results_metrics_df3 <- results_metrics_df2 %>%
              dplyr::filter(
                target_to_source_std_ratio == !!target_to_source_std_ratio |
                  is.na(target_to_source_std_ratio),
              )

            for (target_sample_size_per_arm in target_sample_sizes) {
              results_metrics_df4 <- results_metrics_df3 %>%
                dplyr::filter(
                  target_sample_size_per_arm == !!target_sample_size_per_arm
                )


              if (length(unique(results_metrics_df4$method))<2){
                next
              }
              plot_metric_vs_drift_methods(
                metric = metric,
                results_metrics_df = results_metrics_df4,
                case_study = case_study,
                theta_0 = theta_0,
                source_denominator_change_factor = source_denominator_change_factor,
                target_to_source_std_ratio = target_to_source_std_ratio,
                target_sample_size_per_arm = target_sample_size_per_arm,
                xvars = xvars,
                join_points= FALSE
              )
          }

          for (target_treatment_effect in treatment_effects_labels){
            plot_metric_vs_sample_size_methods(
              metric = metric,
              results_metrics_df = results_metrics_df3,
              case_study = case_study,
              theta_0 = theta_0,
              target_to_source_std_ratio = target_to_source_std_ratio,
              target_treatment_effect = target_treatment_effect,
              join_points = TRUE,
              dodging = TRUE
            )
            }
          }
        }
      }
    }
  }
}
