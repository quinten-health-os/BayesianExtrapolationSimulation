plot_posterior_vs_prior_parameters <- function(input_df,
                                               method,
                                               case_study,
                                               target_sample_size_per_arm,
                                               source_denominator_change_factor,
                                               target_to_source_std_ratio) {

  if (method == "commensurate_power_prior"){
    warning("Not implemented for the commensurate power prior, refer to the corresponding table instead.")
    return()
  }
  results_metrics_df <- input_df %>%
    dplyr::filter(
      case_study == !!case_study,
      target_sample_size_per_arm == !!target_sample_size_per_arm,
      method == !!method,
      source_denominator_change_factor == !!source_denominator_change_factor  |
        is.na(source_denominator_change_factor),
      target_to_source_std_ratio == !!target_to_source_std_ratio |
        is.na(target_to_source_std_ratio)
    )

  if (nrow(results_metrics_df) == 0) {
    warnings("Dataframe is empty.")
    return()
  }

  source_treatment_effect_estimate <- unique(results_metrics_df$source_treatment_effect_estimate)[1]

  case_study_config <- yaml::read_yaml(paste0(case_studies_config_dir, case_study, ".yml"))

  mandatory_drift_values <- important_drift_values(source_treatment_effect_estimate, case_study_config)

  treatment_effects_names <- c("No effect",
                               "Partially consistent effect",
                               "Consistent effect")

  color_map <- viridis::viridis(length(mandatory_drift_values))

  posterior_parameters_df <- get_parameters(results_metrics_df[, "posterior_parameters"])
  prior_parameters_df <- get_parameters(results_metrics_df[, "parameters"])

  # Loop over the prior parameters
  k <- 0
  for (key in names(methods_dict[[method]])) {
    k <- k+1
    # Check if the parameter is updated
    if (method != "commensurate_power_prior" & !key %in% colnames(posterior_parameters_df)) {
      next
    }

    if (length(methods_dict[[method]][[key]][['range']]) < 2){
      next
    }

    # Select the other parameters
    other_parameters <- unique(prior_parameters_df[,-k])

    if (is.null(nrow(other_parameters)) || nrow(other_parameters) == 0){
      n_params_loop = 1
    } else {
      n_params_loop = nrow(other_parameters)
    }

    # Loop over other prior parameters values
    for (j in seq(n_params_loop)){
      if ((is.null(nrow(other_parameters)) || nrow(other_parameters) == 0) || method == "commensurate_power_prior"){
        other_params_label <- ""
        other_params_str <- ""
        prior_parameters_subdf <- prior_parameters_df
        posterior_parameters_subdf <- posterior_parameters_df
        results_metrics_subdf <- results_metrics_df
      } else {
        other_params_label <- make_labels_from_parameters(parameter = other_parameters[j,], method = method)
        other_params_label <- paste0(other_params_label, ", ")
        other_params_str <- convert_params_to_str(methods_dict[[method]], other_parameters[j,])

        # Filter the dataframe on the value of these other parameters

        # matching_filter <- apply(prior_parameters_df[,-k], 1, function(row) compare_ignore_na(row, other_parameters[j,]))

        matching_filter <- apply(prior_parameters_df[,-k], 1, function(row) all(row == other_parameters[j,]))
        prior_parameters_subdf <- prior_parameters_df[matching_filter,]

        posterior_parameters_subdf <- posterior_parameters_df[matching_filter,]

        results_metrics_subdf <- results_metrics_df[matching_filter,]
      }


      # Loop over the posterior parameters
      posterior_colnames <- colnames(posterior_parameters_subdf)
      posterior_keys <- posterior_colnames[!grepl("^conf_int", posterior_colnames)]

      for (posterior_key in posterior_keys){

        # Combine all necessary data based on the drift values into a single data frame
        combined_data <- do.call(rbind, lapply(mandatory_drift_values, function(drift) {
          filter_df <- abs(results_metrics_subdf$drift - drift) < .Machine$double.eps ^ 0.5
          posterior_parameters <- posterior_parameters_subdf[filter_df, ]
          prior_parameters <- prior_parameters_subdf[filter_df, ]
          filtered_results_metrics_df <- results_metrics_df[filter_df, ]

          parameters_labels <- unname(sapply(1:nrow(filtered_results_metrics_df), function(i) {
            process_method_parameters_label(filtered_results_metrics_df[i, ],
                                            methods_labels,
                                            method_name = FALSE)
          }))
          parameters_labels <- setNames(object = parameters_labels, nm = filtered_results_metrics_df$parameters)

          x <- as.numeric(prior_parameters[[key]])
          y <- as.numeric(posterior_parameters[[posterior_key]])
          conf_int_lower <- unlist(posterior_parameters[[paste0("conf_int_lower_", posterior_key)]])
          conf_int_upper <- unlist(posterior_parameters[[paste0("conf_int_upper_", posterior_key)]])

          if (is.null(x)) {
            stop("x is null")
          }

          data.frame(
            x = x,
            y = y,
            conf_int_lower = conf_int_lower,
            conf_int_upper = conf_int_upper,
            drift = drift
          )
        }))

        plot_title <- sprintf(
          "%s, %s, %s $N_T/2 = $ %s",
          str_to_title(case_study),
          methods_labels[[method]]$label,
          other_params_label,
          target_sample_size_per_arm
        )

        plot_title <- format_title(title = plot_title, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)

        # Calculate the range of the x-axis and the cap size
        x_range <- range(combined_data$x)
        x_width <- diff(x_range)
        cap_size <- x_width * relative_error_cap_width # Adjust relative_error_cap_width to control the relative cap size

        # Create the plot
        plt <- ggplot(combined_data, aes(x = x, y = y, ymin = conf_int_lower, ymax = conf_int_upper, color = as.factor(drift))) +
          geom_errorbar(width = cap_size, position = position_dodge(width = cap_size / 2)) +
          geom_point(size = markersize, position = position_dodge(width = cap_size / 2)) +
          scale_color_manual(values = color_map, name = NULL, labels = treatment_effects_names) +
          ggplot2::labs(
            title = plot_title,
            x = TeX(sprintf("%s %s", "Prior", methods_dict[[method]][[key]]$parameter_notation)),
            y = TeX(sprintf("%s %s", "Posterior", methods_dict[[method]][[posterior_key]]$parameter_notation)),
            color = "Target treatment effect",
            group = "Target treatment effect"
          )

        # Save the plot as PDF and PNG
        directory <- paste0(figures_dir, case_study)
        if (!dir.exists(directory)) {
          dir.create(directory, showWarnings = FALSE, recursive = TRUE)
        }

        figure_name <- paste0(
          case_study,
          "_posterior_",
          methods_dict[[method]][[posterior_key]]$parameter_label,
          "_vs_prior_",
          methods_dict[[method]][[key]]$parameter_label,
          "_",
          method,
          "_sample_size_",
          target_sample_size_per_arm
        )

        if (other_params_str != ""){
          figure_name <- paste0(figure_name, "_", other_params_str)
        }

        figure_name <- format_filename(filename = figure_name, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)


        plot.size <- set_size(textwidth)
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
}


#' Function to plot posterior parameters vs drift
#'
#' @description This function plots the posterior parameters against the drift for a given method, case study, sample size, and control drift option.
#'
#' @param results_metrics_df The data frame containing the results and metrics.
#' @param method The method to plot the parameters for.
#' @param case_study The case study to plot the parameters for.
#' @param target_sample_size_per_arm The sample size to plot the parameters for.
#' @param control_drift A logical value indicating whether to filter the data based on control drift.
#' @param xvars Variable on the x axis
#'
#' @return None
#'
#' @export
plot_posterior_parameters_vs_drift <- function(results_metrics_df,
                                               method,
                                               case_study,
                                               target_sample_size_per_arm,
                                               control_drift,
                                               xvars,
                                               source_denominator_change_factor,
                                               target_to_source_std_ratio) {
  # Filter the data frame based on method, case_study, and target_sample_size_per_arm
  results_metrics_df <- results_metrics_df %>%
    dplyr::filter(
      method == !!method,
      case_study == !!case_study,
      target_sample_size_per_arm == !!target_sample_size_per_arm,
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

  posterior_parameters_names <- colnames(get_parameters(results_metrics_df[, "posterior_parameters"]))


  reference_values_all_parameters <- list()
  for (key in names(methods_dict[[method]])) {
    reference_values_all_parameters[[key]] <- methods_dict[[method]][[key]][["important_values"]]
  }

  # Get important parameters values combinations across parameters
  #Initialize an empty list to store important values for each key
  all_important_parameters_values <- list()

  # Iterate over all keys
  for (key in names(methods_dict[[method]])) {

    # Get the important values or the range if important values are not defined
    if (is.null(methods_dict[[method]][[key]][["important_values"]])) {
      important_parameters_values <- methods_dict[[method]][[key]][["range"]]
    } else {
      important_parameters_values <- methods_dict[[method]][[key]][["important_values"]]
    }

    # Append the important values for the current key to the list
    all_important_parameters_values[[key]] <- important_parameters_values
  }

  # Generate all combinations of important values for all keys
  parameters_combinations <- expand.grid(all_important_parameters_values)

  parameters_names <- names(methods_dict[[method]])

  if (method == "commensurate_power_prior"){
    parameters_combinations <- get_parameters(results_metrics_df[, "parameters"])
  }

  for (key in parameters_names) {
    # Loop over the method's parameters
    if (method != "commensurate_power_prior" & !key %in% posterior_parameters_names) {
      next
    }

    # Get the important prior values for which we will make plots (one per important value)
    if (is.null(methods_dict[[method]][[key]][["important_values"]])) {
      important_parameters_values <- methods_dict[[method]][[key]][["range"]]
    } else {
      important_parameters_values <- methods_dict[[method]][[key]][["important_values"]]
    }

    if (!(method == "commensurate_power_prior")){

      # Other prior parameters can take different values, so we filter based on their reference value
      # Filter based on reference values of other parameters
      filter_conditions <- reference_values_all_parameters
      filter_conditions <- filter_conditions[!names(filter_conditions) %in% key]

      # Filter the results dataframe based on the filter conditions
      # Create a boolean filter for df1
      prior_parameters_df <- get_parameters(results_metrics_df[, "parameters"])
      boolean_filter <- create_boolean_filter(prior_parameters_df, filter_conditions, exclude_key = key)


      # Apply the filter to the dataframe
      results_metrics_subdf <- results_metrics_df[boolean_filter, ]
    } else {
      results_metrics_subdf <- results_metrics_df
    }

    posterior_parameters_df <- get_parameters(results_metrics_subdf[, "posterior_parameters"])

    prior_parameters_df <- get_parameters(results_metrics_subdf[, "parameters"])

    # Get the labels for the prior parameters
    parameters_labels <- unname(sapply(1:nrow(results_metrics_subdf), function(i) {
      process_method_parameters_label(results_metrics_subdf[i, ], methods_labels, method_name = FALSE)
    }))
    parameters_labels <- setNames(object = parameters_labels, nm = results_metrics_subdf$parameters)

    plot_data <- list()
    prior_param_labels <- list()
    for (i in 1:nrow(parameters_combinations)) {
      prior_param <- parameters_combinations[i,]

      if (method == "commensurate_power_prior"){
        # Compare each row of prior_parameters_df with the single row of prior_param
        filter_df <- apply(prior_parameters_df, 1, function(row) compare_ignore_na(row, as.vector(prior_param)))
      } else {
        filter_df <- prior_parameters_df[key] == prior_param[[key]]
      }

      filter_results_metrics_df <- results_metrics_subdf[filter_df, ]

      drift <- unlist(filter_results_metrics_df[, xvar[["name"]]])

      if (method == "commensurate_power_prior"){
        keys <- c("heterogeneity_parameter_mean", "heterogeneity_parameter_std", "power_parameter_mean", "power_parameter_std")
        y <- unlist(posterior_parameters_df[filter_df, keys])
        conf_int_lower <- unlist(posterior_parameters_df[filter_df, paste0("conf_int_lower_", keys)])
        conf_int_upper <- unlist(posterior_parameters_df[filter_df, paste0("conf_int_upper_", keys)])
      } else {
        y <- unlist(posterior_parameters_df[filter_df, key])
        conf_int_lower <- unlist(posterior_parameters_df[filter_df, paste0("conf_int_lower_", key)])
        conf_int_upper <- unlist(posterior_parameters_df[filter_df, paste0("conf_int_upper_", key)])
      }

      plot_data[[i]] <- data.frame(
        drift = drift,
        y = y,
        ymin = conf_int_lower,
        ymax = conf_int_upper,
        parameters = filter_results_metrics_df$parameters,
        method = method
      )

      prior_param_labels[[i]] <- parameters_labels[filter_df]
    }

    # Combine all plot data
    plot_data_df <- dplyr::bind_rows(plot_data)

    parameters_labels <- unname(sapply(1:nrow(plot_data_df), function(i) {
      process_method_parameters_label(
        plot_data_df[i, ],
        methods_labels,
        method_name = FALSE,
        parameters_colname = "parameters"
      )
    }))

    parameters_labels <- parameters_labels[!duplicated(plot_data_df$parameters)]

    # Calculate the range of the x-axis
    x_range <- range(drift)
    x_width <- diff(x_range)
    # Set a relative cap size (necessary, otherwise the cap size will depend on the x-axis width)
    cap_size <- x_width * relative_error_cap_width # Adjust relative_error_cap_width to control the relative cap size

    plot_title <-
      sprintf(
        "%s, %s, $N_T/2 =$ %s",
        str_to_title(case_study),
        methods_labels[[method]]$label,
        target_sample_size_per_arm
      )

    plot_title <- format_title(title = plot_title, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)

    plt <- ggplot2::ggplot()
    plt <- plt + geom_errorbar(
      data = plot_data_df,
      ggplot2::aes(
        x = drift,
        y = y,
        ymin = ymin,
        ymax = ymax,
        color = parameters,
        group = parameters
      ),
      position = position_dodge(width = cap_size / 2),
      width = cap_size
    ) +
      geom_point(
        data = plot_data_df,
        ggplot2::aes(
          x = drift,
          y = y,
          color = parameters,
          group = parameters
        ),
        position = position_dodge(width = cap_size / 2),
        size = markersize,
      ) +
      scale_color_viridis_d(labels = parameters_labels) + # Apply the viridis color map for discrete data
      ggplot2::labs(color = "Prior parameters value") +
      ggtitle(plot_title)

    plt <- plt + ggplot2::labs(x = xvar[["label"]], y = TeX(sprintf("Posterior %s", methods_dict[[method]][[key]][["parameter_notation"]])))


    case_study_config <- yaml::read_yaml(paste0(case_studies_config_dir, case_study, ".yml"))
    source_treatment_effect_estimate <- unique(results_metrics_df$source_treatment_effect_estimate)[1]
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
      "_posterior_",
      methods_dict[[method]][[key]]$parameter_label,
      "_vs_",
      xvar$name,
      "_sample_size_",
      target_sample_size_per_arm
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


#' Plot frequentist methods operating characteristics
#'
#' @description This function generates plots for the operating characteristics of frequentist methods.
#'
#' @param results_metrics_df The data frame containing the results and metrics.
#'
#' @return None
#'
#' @examples NA
posterior_parameters_plots <- function(results_metrics_df) {
  # Get the list of case studies
  case_studies <- unique(results_metrics_df$case_study)

  # Get the list of methods
  methods <- unique(results_metrics_df$method)

  for (case_study in case_studies) {
    results_metrics_df1 <- results_metrics_df[results_metrics_df$case_study == case_study,]
    for (method in methods) {
      results_metrics_df2 <- results_metrics_df1[results_metrics_df1$method == method,]
      target_sample_sizes <- unique(results_metrics_df2$target_sample_size_per_arm)
      for (target_sample_size_per_arm in target_sample_sizes) {

        results_metrics_df3 <- results_metrics_df2[results_metrics_df2$target_sample_size_per_arm == target_sample_size_per_arm,]

        target_to_source_std_ratios <- unique(results_metrics_df3$target_to_source_std_ratio)

        for (target_to_source_std_ratio in target_to_source_std_ratios) {
          results_metrics_df4 <- results_metrics_df3 %>%
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

            # Call the posterior_parameters_vs_drift function
            plot_posterior_parameters_vs_drift(
              results_metrics_df5,
              method,
              case_study,
              target_sample_size_per_arm,
              control_drift = FALSE,
              xvars = xvars,
              source_denominator_change_factor = source_denominator_change_factor,
              target_to_source_std_ratio = target_to_source_std_ratio
            )
            # Call the posterior_vs_prior_parameters function
            plot_posterior_vs_prior_parameters(
              results_metrics_df5,
              method,
              case_study,
              target_sample_size_per_arm = target_sample_size_per_arm,
              source_denominator_change_factor = source_denominator_change_factor,
              target_to_source_std_ratio = target_to_source_std_ratio
            )
          }
        }
      }
    }
  }
}
