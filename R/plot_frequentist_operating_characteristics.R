plot_element_metric_vs_xvar <- function(plt,
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
                                        parameters_labels = NA) {
  if (nrow(data) == 0) {
    stop("The dataframe is empty.")
  }

  if (is_baseline == FALSE) {
    if (!(labels_title == "Parameters")) {
      if (is.null(data$label)) {
        stop("The dataframe does not contain a 'label' column")
      }

      plt <- plt +
        geom_errorbar(
          data = data,
          ggplot2::aes(
            x = !!xvar_name,
            y = !!selected_metric_name,
            ymin = !!selected_metric_uncertainty_lower,
            ymax = !!selected_metric_uncertainty_upper,
            color = as.factor(label),
            group = as.factor(label),
          ),
          position = position_dodge(width = cap_size),
          width = cap_size
        ) +
        geom_point(
          data = data,
          ggplot2::aes(
            x = !!xvar_name,
            y = !!selected_metric_name,
            group = as.factor(label),
            color = as.factor(label),
          ),
          position = position_dodge(width = cap_size),
          size = markersize,
        ) +
        scale_color_viridis_d() + # Apply the viridis color map for discrete data
        ggplot2::labs(color = labels_title)

      if (join_points == TRUE) {
        plt <- plt + geom_line(
          #  Add lines between points
          data = data,
          ggplot2::aes(
            x = !!xvar_name,
            y = !!selected_metric_name,
            group = as.factor(label),
            color = as.factor(label),
          ),
          position = position_dodge(width = cap_size)
        )
      }
    } else {
      parameters_labels <- setNames(object = parameters_labels, nm = data$parameters)

      if (length(unique(data$parameters)) > 1) {
        plt <- plt +
          geom_errorbar(
            data = data,
            ggplot2::aes(
              x = !!xvar_name,
              y = !!selected_metric_name,
              ymin = !!selected_metric_uncertainty_lower,
              ymax = !!selected_metric_uncertainty_upper,
              color = factor(parameters, levels = unique(data$parameters)),
              # Used to keep the parameters ordered
              group = factor(parameters, levels = unique(data$parameters)),
            ),
            position = position_dodge(width = cap_size),
            width = cap_size
          ) +
          geom_point(
            data = data,
            ggplot2::aes(
              x = !!xvar_name,
              y = !!selected_metric_name,
              color = factor(parameters, levels = unique(data$parameters)),
              # Used to keep the parameters ordered
              group = factor(parameters, levels = unique(data$parameters)),
            ),
            position = position_dodge(width = cap_size),
            size = markersize,
          ) +
          scale_color_viridis_d(labels = parameters_labels) + labs(color = labels_title)

        if (join_points == TRUE) {
          plt <- plt + geom_line(
            #  Add lines between points
            data = data,
            ggplot2::aes(
              x = !!xvar_name,
              y = !!selected_metric_name,
              color = factor(parameters, levels = unique(data$parameters)),
              # Used to keep the parameters ordered
              group = factor(parameters, levels = unique(data$parameters)),
            ),
            position = position_dodge(width = cap_size)
          )
        }
      } else {
        plt <- plt +
          geom_errorbar(
            data = data,
            ggplot2::aes(
              x = !!xvar_name,
              y = !!selected_metric_name,
              ymin = !!selected_metric_uncertainty_lower,
              ymax = !!selected_metric_uncertainty_upper
            ),
            position = position_dodge(width = cap_size),
            width = cap_size
          ) +
          geom_point(
            data = data,
            ggplot2::aes(
              x = !!xvar_name,
              y = !!selected_metric_name
            ),
            position = position_dodge(width = cap_size),
            size = markersize,
          )
        if (join_points == TRUE) {
          plt <- plt + geom_line(
            #  Add lines between points
            data = data,
            ggplot2::aes(
              x = !!xvar_name,
              y = !!selected_metric_name
            ),
            position = position_dodge(width = cap_size)
          )
        }
      }
    }
  } else {
    plt <- plt +
      geom_errorbar(
        data = data,
        ggplot2::aes(
          x = !!xvar_name,
          y = !!selected_metric_name,
          ymin = !!selected_metric_uncertainty_lower,
          ymax = !!selected_metric_uncertainty_upper,
        ),
        position = position_dodge(width = cap_size),
        width = cap_size,
        color = color
      ) +
      geom_point(
        data = data,
        ggplot2::aes(
          x = !!xvar_name,
          y = !!selected_metric_name,
        ),
        position = position_dodge(width = cap_size),
        size = markersize,
        color = color
      )
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
#' @param control_drift Logical indicating whether to filter results by control drift (TRUE) or treatment drift (FALSE)
#' @param target_sample_size_per_arm The target sample size per arm
#' @param parameters_combinations The combinations of parameters to filter the results by
#' @param xvars A list of x-variables for control drift and treatment drift
#' @param add_baselines Logical indicating whether to add baselines to the plot
#' @param target_to_source_std_ratio Ratio between the target and source studies sampling standard deviation.
#' @param join_points Whether to join the point with a line or not.
#'
#' @return None
#'
#' @export
plot_metric_vs_drift <- function(metric,
                                 results_metrics_df,
                                 theta_0,
                                 case_study,
                                 method,
                                 category,
                                 control_drift,
                                 target_sample_size_per_arm,
                                 parameters_combinations,
                                 xvars,
                                 add_baselines = FALSE,
                                 target_to_source_std_ratio = 1,
                                 join_points = FALSE,
                                 source_denominator_change_factor = 1,
                                 analysis_config) {

  if (!(metric %in% colnames(results_metrics_df))) {
    warning("Metric not in input dataframe.")
    return()
  }

  # Filter results by control drift
  if (control_drift) {
    # Remove non-zero treatment drift
    results_metrics_df <- results_metrics_df %>% dplyr::filter(treatment_drift == 0)
    xvar <- xvars$control_drift
  } else {
    # Remove non-zero control drift
    results_metrics_df <- results_metrics_df %>% dplyr::filter(control_drift == 0)
    xvar <- xvars$drift
  }

  # Add baselines if required
  if (add_baselines) {
    separate_df <- results_metrics_df %>% dplyr::filter(method == "separate")
    pooling_df <- results_metrics_df %>% dplyr::filter(method == "pooling")
  }

  # Filter results by method
  results_df <- results_metrics_df %>% dplyr::filter(method == !!method)

  categories <- unique(results_df[[category]])
  color_map <- viridis::viridis(n = 256)

  # Filter on the scenarios of main interest
  if (category == "parameters") {
    results_df <- results_df %>%
      dplyr::filter(
        target_sample_size_per_arm == !!target_sample_size_per_arm,
        source_denominator_change_factor == !!source_denominator_change_factor |
          is.na(source_denominator_change_factor),
        target_to_source_std_ratio == !!target_to_source_std_ratio |
          is.na(target_to_source_std_ratio),
      )

    params <- get_parameters(results_df[, "parameters"])

    # Filter on the main parameter values of interest
    # TODO: remove this hard-coded part, use important_parameter_values instead
    if (method == "RMP"){
      main_values_filter <- params$prior_weight %in% c(0, 0.5, 1)
      results_df <- results_df[main_values_filter,]
      params <- params[main_values_filter,]
    }

    #  Make sure that the data are sorted according to the parameters values
    column_names <- colnames(params)
    sorted_index <- do.call(order, lapply(params[column_names], as.factor))
    # Reorder the data :
    results_df[results_df$method == method, ] <- results_df[sorted_index, ]


  } else if (category == "source_denominator_change_factor") {
    if (!(case_study %in% c("teriflunomide", "belimumab", "mepolizumab"))){
      return()
    } #TODO : remove this hard-coded piece
    results_df <- results_df %>%
      dplyr::filter(
        target_sample_size_per_arm == !!target_sample_size_per_arm,
        target_to_source_std_ratio == !!target_to_source_std_ratio |
          is.na(target_to_source_std_ratio),
      )
  } else if (category == "target_to_source_std_ratio") {
    if (!(case_study %in% c("botox", "dapagliflozin"))){
      return()
    } #TODO : remove this hard-coded piece
    results_df <- results_df %>%
      dplyr::filter(
        target_sample_size_per_arm == !!target_sample_size_per_arm,
        source_denominator_change_factor == !!source_denominator_change_factor |
          is.na(source_denominator_change_factor)
      )
  } else if (category == "target_sample_size_per_arm") {
    results_df <- results_df %>%
      dplyr::filter(
        target_to_source_std_ratio == !!target_to_source_std_ratio | is.na(target_to_source_std_ratio),
        source_denominator_change_factor == !!source_denominator_change_factor |
          is.na(source_denominator_change_factor)
      )
  } else {
    stop("Category is not supported")
  }

  if (nrow(results_df) == 0){
    return()
  }

  parameters_labels <- NA
  if (category == "target_sample_size_per_arm") {
    results_df$label <- results_df$target_sample_size_per_arm
    labels_title <- "Target sample size per arm"
  } else if (category == "parameters") {
    parameters_labels <- unname(sapply(1:nrow(results_df), function(i) {
      process_method_parameters_label(results_df[i, ], methods_labels, method_name = FALSE)
    }))

    labels_title <- "Parameters"
  } else if (category == "source_denominator_change_factor") {
    results_df$label <- results_df$source_denominator_change_factor
    labels_title <- "Source denom. change factor"
  } else if (category == "target_to_source_std_ratio") {
    results_df$label <- results_df$target_to_source_std_ratio
    labels_title <- latex2exp::TeX("$\\sigma_T/\\sigma_S$") # Ratio between target and source standard deviation
  } else {
    stop("Not implemented for this category")
  }

  results_df <- cbind(results_df, get_parameters(results_df[, "parameters"]))


  if (category != "parameters") {
    # Filter rows that match parameters_combinations
    results_df <- merge(results_df, parameters_combinations[, ])
    param_values <- unlist(unique(results_df[, category]))
  } else {
    param_values <- unlist(get_parameters(unique(results_df[, category, drop = FALSE])))
  }

  # Do not make the plot if all values for the parameters of interest are NA
  if (all(sapply(param_values, is.na))) {
    return()
  }

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


  colors <- list()
  labels <- list()
  for (i in seq_along(categories)) {
    combination <- categories[i]

    norm <- i / length(categories)
    # Get the color
    colors <- cbind(colors, color_map[as.integer(norm * 255) + 1])

    if (category == "source_denominator_change_factor") {
      if (is.na(combination)) {
        next
      } else {
        results_df[results_df[, "source_denominator_change_factor"] == combination, "label"] <- combination
      }
    }
  }

  if (nrow(results_df) == 0) {
    warning("Dataframe is empty")
    next
  }

  # Calculate the range of the x-axis
  x_range <- range(results_df[, xvar$name])
  x_width <- diff(x_range)
  # Set a relative cap size (necessary, otherwise the cap size will depend on the x-axis width)

  cap_size <- x_width * relative_error_cap_width # Adjust relative_error_cap_width to control the relative cap size

  plt <- ggplot2::ggplot()
  plt <- plot_element_metric_vs_xvar(
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
    parameters_labels = parameters_labels
  )

  if (control_drift) {
    plt <- plt +
      geom_vline(xintercept = 0,
                 linetype = "dashed",
                 color = "gray")
  } else {
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

  }

  if(metric %in% c("tie", "success_proba")){
    plt <- plt +
      ggplot2::geom_hline(yintercept = analysis_config$nominal_tie, color = "black", linetype = "dashed") + # Add horizontal line at y = 0.05
      ggplot2::scale_y_continuous(breaks = function(x) unique(c(pretty(x), analysis_config$nominal_tie)))  # Add y-tick at y = 0.05
  }


  if (add_baselines &&
      nrow(separate_df) != 0 && nrow(pooling_df) != 0) {
    if (category != "target_sample_size_per_arm") {
      separate_df <- separate_df %>%
        dplyr::filter(target_sample_size_per_arm == !!target_sample_size_per_arm)

      pooling_df <- pooling_df %>%
        dplyr::filter(target_sample_size_per_arm == !!target_sample_size_per_arm)
    }

    if (category == "target_sample_size_per_arm") {
      separate_df$label <- separate_df$target_sample_size_per_arm
      pooling_df$label <- pooling_df$target_sample_size_per_arm
    } else if (category == "parameters") {
      separate_df$label <- separate_df$parameters
      pooling_df$label <- pooling_df$parameters
    } else if (category == "source_denominator_change_factor") {
      separate_df$label <- separate_df$source_denominator_change_factor
      pooling_df$label <- pooling_df$source_denominator_change_factor
    } else if (category == "target_to_source_std_ratio") {
      separate_df$label <- separate_df$target_to_source_std_ratio
      pooling_df$label <- pooling_df$target_to_source_std_ratio
    }

    plt <- plot_element_metric_vs_xvar(
      plt,
      separate_df,
      xvar_name,
      selected_metric_name,
      selected_metric_uncertainty_lower,
      selected_metric_uncertainty_upper,
      cap_size,
      markersize,
      labels_title,
      join_points,
      color = "red",
      is_baseline = TRUE,
      parameters_labels = parameters_labels
    )

    plt <- plot_element_metric_vs_xvar(
      plt,
      pooling_df,
      xvar_name,
      selected_metric_name,
      selected_metric_uncertainty_lower,
      selected_metric_uncertainty_upper,
      cap_size,
      markersize,
      labels_title,
      join_points,
      color = "black",
      is_baseline = TRUE,
      parameters_labels = parameters_labels
    )
  }

  if (!(category == "parameters")){
    parameters_label <- make_labels_from_parameters(parameters_combinations, method)
    parameters_str <- convert_params_to_str(methods_dict[[method]], parameters_combinations)

    if (parameters_label == "") {
      parameters_label_title <- ""
    } else {
      parameters_label_title <- paste0(", ", parameters_label)
    }
  }


  if (category == "target_sample_size_per_arm") {
    title <- sprintf(
        "%s, %s %s",
        str_to_title(case_study),
        methods_labels[[method]]$label,
        parameters_label_title
      )

    figure_name <- paste0(
      case_study,
      "_",
      method,
      "_",
      selected_metric_name,
      "_vs_",
      xvar_name,
      "_cat_target_sample_size_parameters_",
      parameters_str
    )
  } else if (category == "source_denominator_change_factor") {
    title <- sprintf(
        "%s, %s %s, %s %s",
        str_to_title(case_study),
        methods_labels[[method]]$label,
        parameters_label_title,
        "$N_T/2 = $",
        target_sample_size_per_arm
      )
    figure_name <- paste0(
      case_study,
      "_",
      method,
      "_",
      selected_metric_name,
      "_vs_",
      xvar_name,
      "_cat_source_denominator_sample_size=",
      target_sample_size_per_arm,
      parameters_str
    )

    source_denominator_change_factor <- NA
  } else if (category == "parameters") {
    title <- sprintf(
        "%s, %s, $N_T/2 = $ %s",
        str_to_title(case_study),
        methods_labels[[method]]$label,
        target_sample_size_per_arm
      )

    figure_name <- paste0(
      case_study,
      "_",
      method,
      "_",
      selected_metric_name,
      "_vs_",
      xvar_name,
      "_cat_parameters_sample_size=",
      target_sample_size_per_arm
    )
  } else if (category == "target_to_source_std_ratio") {
    title <- sprintf(
        "%s, %s%s, $N_T/2 = $%s",
        str_to_title(case_study),
        methods_labels[[method]]$label,
        parameters_label_title,
        target_sample_size_per_arm
      )

    figure_name <- paste0(
      case_study,
      "_",
      method,
      "_",
      selected_metric_name,
      "_vs_",
      xvar_name,
      "_cat_target_to_source_std_ratio",
      "_target_sample_size_per_arm=",
      target_sample_size_per_arm,
      parameters_str
    )

    target_to_source_std_ratio <- NA
  } else {
    stop("Not implemented for this category")
  }

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
  fig_height_in <- plot.size[2]

  file_path <- file.path(directory, figure_name)

  if (file.exists(paste0(file_path, ".pdf")) && file.exists(paste0(file_path, ".png")) && remake_figures == FALSE){
    return()
  }
  export_plots(plt, file_path, fig_width_in, fig_height_in, type = "pdf")
  export_plots(plt, file_path, fig_width_in, fig_height_in, type = "png")
}

#' Function to plot metric vs parameters
#'
#' @description This function plots a metric against the parameters for a given method, case study, and sample size.
#'
#' @param results_metrics_df The data frame containing the results and metrics.
#' @param metric The metric to plot.
#' @param case_study The case study to plot the metric for.
#' @param method The method to plot the metric for.
#' @param target_sample_size_per_arm The target sample size per arm to plot the metric for.
#' @param theta_0 Boundary of the null hypothesis space.
#' @param add_baselines A logical value indicating whether to add baselines to the plot.
#'
#' @return None
#'
#' @export
plot_metric_vs_parameters <- function(results_metrics_df,
                                      metric,
                                      case_study = "belimumab",
                                      method = "RMP",
                                      target_sample_size_per_arm = 93,
                                      theta_0 = NULL,
                                      add_baselines = FALSE,
                                      target_to_source_std_ratio = 1,
                                      source_denominator_change_factor = 1,
                                      analysis_config) {
  if (!(metric %in% colnames(results_metrics_df))) {
    warning("Metric not in input dataframe.")
    return()
  }

  # Filter results_metrics_df on the selected case study, target sample arm and method
  results_df <- results_metrics_df %>%
    dplyr::filter(
      case_study == !!case_study,
      target_sample_size_per_arm == !!target_sample_size_per_arm,
      method == !!method,
      source_denominator_change_factor == !!source_denominator_change_factor |
        is.na(source_denominator_change_factor),
      target_to_source_std_ratio == !!target_to_source_std_ratio |
        is.na(target_to_source_std_ratio)
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
  }

  parameters_df <- get_parameters(results_df[, "parameters"])

  methods_parameters <- methods_dict[[method]]


  color_map <- viridis::viridis(length(unique(results_df$target_treatment_effect)))

  # Effects label
  treatment_effects_names <- c("No effect",
                               "Partially consistent effect",
                               "Consistent effect")

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


  # Loop over the different parameters that take different values in the simulation study
  for (i in 1:ncol(parameters_df)) {
    # Select the other parameters
    other_parameters <- unique(parameters_df[,-i])

    if (is.null(nrow(other_parameters)) || nrow(other_parameters) == 0){
      n_params_loop = 1
    } else {
      n_params_loop = nrow(other_parameters)
    }
    for (j in seq(n_params_loop)){
      if (is.null(nrow(other_parameters)) || nrow(other_parameters) == 0){
        other_params_label <- ""
        other_params_str <- ""
        parameters_subdf <- parameters_df
        results_subdf <- results_df
      } else {
        other_params_label <- make_labels_from_parameters(parameter = other_parameters[j,], method = method)
        other_params_label <- paste0(other_params_label, ", ")
        other_params_str <- convert_params_to_str(methods_dict[[method]], other_parameters[j,])

        # Filter the dataframe on the value of these other parameters
        matching_filter <- apply(parameters_df[,-i], 1, function(row) all(row == other_parameters[j,]))
        parameters_subdf <- parameters_df[matching_filter,]

        results_subdf <- results_df[matching_filter,]
      }

      if (length(unique(parameters_subdf[, i])) < 2) {
        next
      } else {
        parameter_values <- as.numeric(parameters_subdf[, i])
      }

      # Calculate the range of the x-axis
      x_range <- range(parameter_values)
      x_width <- diff(x_range)
      # Set a relative cap size (necessary, otherwise the cap size will depend on the x-axis width)
      cap_size <- x_width * relative_error_cap_width # Adjust relative_error_cap_width to control the relative cap size

      results_subdf$parameter_values <- parameter_values

      plot_title =  sprintf(
          "%s, %s, %s $N_T/2 = $ %s",
          str_to_title(case_study),
          methods_labels[[method]]$label,
          other_params_label,
          target_sample_size_per_arm
        )

      plot_title <- format_title(title = plot_title, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)

      # Produce the plot
      if(!(metric %in% c("tie"))){
        plt <- ggplot2::ggplot(
          results_subdf,
          ggplot2::aes(
            x = parameter_values,
            y = !!rlang::sym(selected_metric_name),
            group = as.factor(target_treatment_effect),
            color = as.factor(target_treatment_effect)
          )
        )
      } else {
        plt <- ggplot2::ggplot(
          results_subdf,
          ggplot2::aes(
            x = parameter_values,
            y = !!rlang::sym(selected_metric_name)
          )
        )
      }

      if(metric %in% c("tie", "success_proba")){
        plt <- plt +
          ggplot2::geom_hline(yintercept = analysis_config$nominal_tie, color = "black", linetype = "dashed") + # Add horizontal line at y = 0.05
          ggplot2::scale_y_continuous(breaks = function(x) unique(c(pretty(x), analysis_config$nominal_tie)))  # Add y-tick at y = 0.05
      }

      plt <- plt + geom_point(position = position_dodge(width = cap_size)) +
      geom_errorbar(
        ggplot2::aes(
          x = parameter_values,
          ymin = !!rlang::sym(selected_metric_uncertainty_lower),
          ymax = !!rlang::sym(selected_metric_uncertainty_upper)
        ),
        position = position_dodge(width = cap_size),
        width = cap_size
      )

      if(!(metric %in% c("tie"))){
        plt <- plt + scale_color_manual(values = color_map,
                                        name = NULL,
                                        labels = treatment_effects_names)
        plt <- plt +
          ggplot2::labs(
            color = "Target treatment effect",
            group = "Target treatment effect"
          )
      }

      plt <- plt +
      ggplot2::labs(
        title = plot_title,
        x = latex2exp::TeX(methods_parameters[[i]]$parameter_notation),
        y = selected_metric$label
      )

      # Add a y-line at theta = 0
      if (add_baselines) {
        plt <- plt + geom_hline(yintercept = theta_0, linetype = "dashed")
      }


      # Save plot
      directory <- paste0(figures_dir, case_study)
      if (!dir.exists(directory)) {
        dir.create(directory,
                   showWarnings = FALSE,
                   recursive = TRUE)
      }

      figure_name <- paste0(
        case_study,
        "_",
        method,
        "_",
        selected_metric_name,
        "_vs_parameters_target_sample_size_per_arm_",
        target_sample_size_per_arm
      )

      if (other_params_str != ""){
        figure_name <- paste0(figure_name, "_", other_params_str)
      }


      figure_name <- format_filename(filename = figure_name, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)


      plot.size <- set_size(426)
      fig_width_in <- plot.size[1]
      fig_height_in <- plot.size[2]

      file_path <- file.path(directory, figure_name)

      if (file.exists(paste0(file_path, ".pdf")) && file.exists(paste0(file_path, ".png")) && remake_figures == FALSE){
        next
      }
      export_plots(plt, file_path, fig_width_in, fig_height_in, type = "pdf")
      export_plots(plt, file_path, fig_width_in, fig_height_in, type = "png")

    }
  }
}

#' Function to plot metric vs sample size
#'
#' @description This function plots a metric against the sample size for a given method, case study, and parameters combinations.
#'
#' @param metric The metric to plot.
#' @param results_metrics_df The data frame containing the results and metrics.
#' @param case_study The case study to plot the metric for.
#' @param method The method to plot the metric for.
#' @param parameters_combinations The parameter combinations to filter the data on.
#' @param theta_0 Boundary of the null hypothesis space.
#' @param add_baselines A logical value indicating whether to add baselines to the plot.
#'
#' @return None
#'
#' @export
plot_metric_vs_sample_size <- function(metric,
                                       results_metrics_df,
                                       case_study,
                                       method,
                                       parameters_combinations,
                                       theta_0 = 0,
                                       add_baselines = FALSE,
                                       target_to_source_std_ratio = 1,
                                       source_denominator_change_factor = 1,
                                       analysis_config) {


  if (!(metric %in% colnames(results_metrics_df))) {
    warning("Metric not in input dataframe.")
    return()
  }

  # Filter results_metrics_df on the selected case study and method
  results_df <- results_metrics_df %>%
    dplyr::filter(
      case_study == !!case_study,
      method == !!method,
      source_denominator_change_factor == !!source_denominator_change_factor |
        is.na(source_denominator_change_factor),
      target_to_source_std_ratio == !!target_to_source_std_ratio |
        is.na(target_to_source_std_ratio),
    )

  if (nrow(results_df) == 0) {
    warning("Dataframe is empty")
    next
  }
  # Format the parameters (from the json string)
  results_df <- cbind(results_df, get_parameters(results_df[, "parameters"]))

  # Filter rows that match parameters_combinations
  results_df <- merge(results_df, parameters_combinations[, ])

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
  }

  # Plot setup
  color_map <- viridis::viridis(length(unique(results_df$target_treatment_effect)))

  # Effects label
  effects <- c("No effect",
               "Partially consistent effect",
               "Consistent effect")

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

  # Produce the plot
  if(!(metric %in% c("tie"))){
    plt <- ggplot2::ggplot(
      results_df,
      ggplot2::aes(
        x = target_sample_size_per_arm,
        y = !!selected_metric_name,
        group = as.factor(target_treatment_effect),
        color = as.factor(target_treatment_effect)
      )
    ) +
      geom_point(position = position_dodge(width = cap_size)) +
      geom_errorbar(
        ggplot2::aes(
          x = target_sample_size_per_arm,
          ymin = !!selected_metric_uncertainty_lower,
          ymax = !!selected_metric_uncertainty_upper,
          group = as.factor(target_treatment_effect),
          color = as.factor(target_treatment_effect)
        ),
        position = position_dodge(width = cap_size),
        width = cap_size
      ) +
      scale_color_manual(values = color_map,
                         name = NULL,
                         labels = effects)
  } else {
    plt <- ggplot2::ggplot(
      results_df,
      ggplot2::aes(
        x = target_sample_size_per_arm,
        y = !!selected_metric_name
      )
    ) +
      geom_point(position = position_dodge(width = cap_size)) +
      geom_errorbar(
        ggplot2::aes(
          x = target_sample_size_per_arm,
          ymin = !!selected_metric_uncertainty_lower,
          ymax = !!selected_metric_uncertainty_upper
        ),
        position = position_dodge(width = cap_size),
        width = cap_size
      )
  }

  if(metric %in% c("tie", "success_proba")){
    plt <- plt +
      ggplot2::geom_hline(yintercept = analysis_config$nominal_tie, color = "black", linetype = "dashed") + # Add horizontal line at y = 0.05
      ggplot2::scale_y_continuous(breaks = function(x) unique(c(pretty(x), analysis_config$nominal_tie)))  # Add y-tick at y = 0.05
  }

  param_label <- make_labels_from_parameters(parameters_combinations, method)

  if (param_label == "") {
    param_label <- paste0("")
  } else {
    param_label <- paste0(", ", param_label)
  }

  params_str <- convert_params_to_str(methods_dict[[method]], parameters_combinations)


  title = sprintf(
      "%s, %s%s",
      str_to_title(case_study),
      methods_labels[[method]]$label,
      param_label
    )

  title <- format_title(title = title, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)

  plt <- plt +
    ggplot2::labs(title = title, x = "Target study sample size per arm",
                  y = selected_label)


  # Add a y-line at theta = 0
  if (add_baselines) {
    plt <- plt + geom_hline(yintercept = theta_0, linetype = "dashed")
  }

  directory <- paste0(figures_dir, case_study)
  if (!dir.exists(directory)) {
    dir.create(directory, showWarnings = FALSE, recursive = TRUE)
  }

  figure_name <- paste0(
    case_study,
    "_",
    method,
    "_",
    selected_metric_name,
    "_vs_target_sample_size_per_arm_",
    params_str
  )

  figure_name <- format_filename(filename = figure_name, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)


  plot.size <- set_size(426)
  fig_width_in <- plot.size[1]
  fig_height_in <- plot.size[2]

  file_path <- file.path(directory, figure_name)

  if (file.exists(paste0(file_path, ".pdf")) && file.exists(paste0(file_path, ".png")) && remake_figures == FALSE){
    return()
  }
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
plot_metric_vs_scenario <- function(results_metrics_df, metrics) {

  analysis_config <- yaml::read_yaml(system.file("conf/analysis_config.yml", package = "RBExT"))

  # Get the list of case studies
  case_studies <- unique(results_metrics_df$case_study)

  # Get the list of methods
  methods <- unique(results_metrics_df$method)

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

      if (method %in% c("separate", "pooling")) {
        add_baselines <- FALSE
      } else {
        add_baselines <- FALSE
      }

      for (metric in names(metrics)) {
        for (sample_size in target_sample_sizes) {
          plot_metric_vs_parameters(
            filtered_results_metrics_df,
            case_study = case_study,
            method = method,
            target_sample_size_per_arm = sample_size,
            metric = metric,
            theta_0 = theta_0,
            analysis_config = analysis_config
          )
        }

        if (!(metric %in% c("tie"))){
          for (sample_size in target_sample_sizes) {
            for (i in 1:nrow(parameters_combinations)) {
              # Run plot_metric_vs_drift function with control_drift equal to FALSE
              categories <- c("source_denominator_change_factor",
                              "target_to_source_std_ratio")

              for (category in categories) {
                if (length(unique(filtered_results_metrics_df[[category]])) == 1) {
                  next
                } else {
                  plot_metric_vs_drift(
                    metric = metric,
                    results_metrics_df = filtered_results_metrics_df,
                    case_study = case_study,
                    method = method,
                    theta_0 = theta_0,
                    category = category,
                    control_drift = FALSE,
                    target_sample_size_per_arm = sample_size,
                    parameters_combinations = parameters_combinations[i, , drop = FALSE],
                    add_baselines = FALSE,
                    xvars = xvars,
                    analysis_config = analysis_config
                  )
                }
              }

              plot_metric_vs_drift(
                metric = metric,
                results_metrics_df = filtered_results_metrics_df,
                case_study = case_study,
                method = method,
                theta_0 = theta_0,
                category = "parameters",
                control_drift = FALSE,
                target_sample_size_per_arm = sample_size,
                parameters_combinations = NULL,
                add_baselines = add_baselines,
                xvars = xvars,
                analysis_config = analysis_config
              )
            }
          }
        }

        for (i in 1:nrow(parameters_combinations)) {
          if (!(metric %in% c("tie"))){
            plot_metric_vs_drift(
              metric = metric,
              results_metrics_df = filtered_results_metrics_df,
              case_study = case_study,
              method = method,
              theta_0 = theta_0,
              category = "target_sample_size_per_arm",
              control_drift = FALSE,
              target_sample_size_per_arm = NULL,
              parameters_combinations = parameters_combinations[i, , drop = FALSE],
              add_baselines = add_baselines,
              xvars = xvars,
              analysis_config = analysis_config
            )
          }

          plot_metric_vs_sample_size(
            metric = metric,
            results_metrics_df = filtered_results_metrics_df,
            case_study = case_study,
            method = method,
            parameters_combinations = parameters_combinations[i, , drop = FALSE],
            theta_0 = theta_0,
            add_baselines = add_baselines,
            analysis_config = analysis_config
          )
        }
      }
    }
  }
}
