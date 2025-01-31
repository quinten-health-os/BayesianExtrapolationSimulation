#' @title Plot Empirical Bayes Hyperparameters vs Drift
#' @description This function plots the hyperparameters updated using empirical Bayes against the drift for a given method, case study sample size, and control drift option.
#'
#' @param results_metrics_df The data frame containing the results and metrics.
#' @param method The method to plot the parameters for.
#' @param case_study The case study to plot the parameters for.
#' @param control_drift A logical value indicating whether to filter the data based on control drift.
#' @param xvars Variable on the x axis
#'
#' @return None
#' @export
plot_empirical_bayes_hyperparameters_vs_drift <- function(results_metrics_df,
                                               method,
                                               case_study,
                                               control_drift,
                                               xvars,
                                               source_denominator_change_factor,
                                               target_to_source_std_ratio,
                                               parameters_combinations) {

  if (nrow(parameters_combinations) >1){
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

  for (key in selected_parameters) {
    results_metrics_subdf <- results_df

    posterior_parameters_df <- get_parameters(results_metrics_subdf[, "posterior_parameters", drop = FALSE])

    filter_results_metrics_df <- results_metrics_subdf

    drift <- unlist(filter_results_metrics_df[, xvar[["name"]]])
    y <- unlist(posterior_parameters_df[, key])

    conf_int_lower <- unlist(posterior_parameters_df[, paste0("conf_int_lower_", key)])
    conf_int_upper <- unlist(posterior_parameters_df[, paste0("conf_int_upper_", key)])

    plot_data <- data.frame(
      drift = drift,
      y = y,
      ymin = conf_int_lower,
      ymax = conf_int_upper,
      parameters = filter_results_metrics_df$parameters,
      method = method,
      target_sample_size_per_arm = factor(filter_results_metrics_df$target_sample_size_per_arm)
    )

    # Combine all plot data
    plot_data_df <- dplyr::bind_rows(plot_data)

    # Calculate the range of the x-axis
    x_range <- range(drift)
    x_width <- diff(x_range)
    # Set a relative cap size (necessary, otherwise the cap size will depend on the x-axis width)
    cap_size <- x_width * relative_error_cap_width # Adjust relative_error_cap_width to control the relative cap size


    if (parameters_label_title == ""){
      plot_title <-
        sprintf(
          "%s, %s",
          str_to_title(case_study),
          methods_labels[[method]]$label
        )
    } else {
      plot_title <-
        sprintf(
          "%s, %s%s",
          str_to_title(case_study),
          methods_labels[[method]]$label,
          parameters_label_title
        )
    }

    plot_title <- format_title(title = plot_title, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)

    plt <- ggplot2::ggplot()
    plt <- plt + geom_errorbar(
      data = plot_data_df,
      ggplot2::aes(
        x = drift,
        y = y,
        ymin = ymin,
        ymax = ymax,
        color = target_sample_size_per_arm,
        group = target_sample_size_per_arm
      ),
      position = position_dodge(width = cap_size / 2),
      width = cap_size
    ) +
      geom_point(
        data = plot_data_df,
        ggplot2::aes(
          x = drift,
          y = y,
          color = target_sample_size_per_arm,
          group = target_sample_size_per_arm
        ),
        position = position_dodge(width = cap_size / 2),
        size = markersize,
      ) +
      ggplot2::labs(color = TeX("$N_T/2$"), group = TeX("$N_T/2$")) +
      scale_color_viridis_d() + # Apply the viridis color map for discrete data
      ggtitle(plot_title)

    plt <- plt + ggplot2::labs(x = xvar[["label"]], y = TeX(sprintf("%s%s", empirical_bayes_hyperparameters[[method]][[key]][["parameter_notation"]], empirical_bayes_hyperparameters[[method]][[key]][["suffix"]])))

    hyperparameter_range <- empirical_bayes_hyperparameters[[method]][[key]][["range"]]

    if (length(hyperparameter_range) == 2 && !is.na(hyperparameter_range[1]) && !is.na(hyperparameter_range[2])){
      plt <- plt + ylim(hyperparameter_range[1], hyperparameter_range[2])
    }

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
      empirical_bayes_hyperparameters[[method]][[key]]$parameter_label,
      "_vs_",
      xvar$name,
      "_cat_target_sample_size_per_arm",
      parameters_str
    )

    figure_name <- format_filename(filename = figure_name, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)

    plot.size <- set_size(textwidth)
    fig_width_in <- plot.size[1]
    fig_height_in <- plot.size[2]

    file_path <- file.path(directory, figure_name)
    export_plots(plt, file_path, fig_width_in, fig_height_in, type = "pdf")
    export_plots(plt, file_path, fig_width_in, fig_height_in, type = "png")
  }
}


#' @title Plot the posterior parameters
#' @description  Plot the hyperparameters updated using empirical Bayes as a function of drift, for each scenario.
#'
#' @param results_metrics_df The data frame containing the results and metrics.
#'
#' @return None
#' @export
empirical_bayes_parameters_plots <- function(results_metrics_df) {
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
            plot_empirical_bayes_hyperparameters_vs_drift(
              results_metrics_df5,
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
