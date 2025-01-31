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
bayesian_metric_vs_parameters <- function(results_metrics_df,
                                          metric,
                                          case_study = "belimumab",
                                          method = "RMP",
                                          target_sample_size_per_arm = 93,
                                          theta_0 = NULL,
                                          add_baselines = FALSE,
                                          source_denominator_change_factor = 1,
                                          target_to_source_std_ratio = 1) {

  # Filter results_metrics_df on the selected case study, target sample arm and method
  results_df <- results_metrics_df %>%
    dplyr::filter(
      case_study == !!case_study,
      target_sample_size_per_arm == !!target_sample_size_per_arm,
      method == !!method,
      source_denominator_change_factor == !!source_denominator_change_factor  |
        is.na(source_denominator_change_factor),
      target_to_source_std_ratio == !!target_to_source_std_ratio |
        is.na(target_to_source_std_ratio)
    )

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

  if (nrow(results_df) == 0) {
    warning("Dataframe is empty")
    return()
  }

  results_df[results_df[, "design_prior_type"] == "analysis_prior", "label"] <- "Analysis prior"
  results_df[results_df[, "design_prior_type"] == "source_posterior", "label"] <- "Source posterior"
  results_df[results_df[, "design_prior_type"] == "ui_design_prior", "label"] <- "UI prior"

  # results_df <- results_df %>%
  #   filter(!is.na(prepost_proba_TP) & !is.na(average_tie) & !is.na(average_power) & !is.na(prior_proba_success) & !is.na(prior_proba_no_benefit) & !is.na(prior_proba_benefit) & !is.na(prepost_proba_FP) & !is.na(upper_bound_proba_FP))

  design_priors <- unlist(unique(results_df[, "label"]))


  parameters_df <- get_parameters(results_df[, "parameters"])

  methods_parameters <- methods_dict[[method]]

  if (metric %in% names(bayesian_metrics)) {
    selected_metric <- bayesian_metrics[[metric]]
  } else {
    stop("Metric is not supported")
  }

  if (is.null(selected_metric$name)) {
    stop("Metric name is NULL")
  }

  selected_metric_name <- rlang::sym(selected_metric$name)
  selected_label <- selected_metric$label

  labels_title <- "Design prior"

  color_map <- viridis::viridis(length(design_priors))


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

      plot_title =
        sprintf(
          "%s, %s, %s $N_T/2 = $ %s",
          str_to_title(case_study),
          methods_labels[[method]]$label,
          other_params_label,
          target_sample_size_per_arm
        )

      plot_title <- format_title(title = plot_title, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)

      if (metric %in% c("prior_proba_benefit", "prior_proba_no_benefit")){
        # Only plot the prior probability of (no) benefit for the analysis prior
        results_subdf <- results_subdf %>%
          dplyr::filter(
            label == "Analysis prior"
          )

      }

      if (metric %in% c("prior_proba_benefit", "prior_proba_no_benefit")){

        # Produce the plot
        plt <- ggplot2::ggplot(
          results_subdf,
          ggplot2::aes(
            x = parameter_values,
            y = !!selected_metric_name
          )
        ) +
          geom_point() +
          geom_line()

        plt <- plt +
          ggplot2::labs(
            title = plot_title,
            x = latex2exp::TeX(methods_parameters[[i]]$parameter_notation),
            y = selected_metric$label
          )
      } else {
        # Produce the plot
        plt <- ggplot2::ggplot(
          results_subdf,
          ggplot2::aes(
            x = parameter_values,
            y = !!selected_metric_name,
            group = as.factor(label),
            color = as.factor(label)
          )
        ) +
          geom_point() +
          geom_line()
        plt <- plt + scale_color_manual(values = color_map,
                           name = NULL,
                           labels = design_priors) +
        ggplot2::labs(
          title = plot_title,
          x = latex2exp::TeX(methods_parameters[[i]]$parameter_notation),
          y = selected_metric$label,
          color = labels_title,
          group = labels_title
        )
      }

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

      plot.size <- set_size(textwidth)
      fig_width_in <- plot.size[1]
      fig_height_in <- plot.size[2]

      file_path <- file.path(directory, figure_name)

      if (file.exists(paste0(file_path, ".pdf")) && file.exists(paste0(file_path, ".png")) && remake_figures == FALSE){
        return()
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
#' @param target_to_source_std_ratio Ratio between the sampling standard deviation in the source study and in the target study.
#'
#' @return None
#'
#' @export
bayesian_metric_vs_sample_size <- function(metric,
                                           input_df,
                                           case_study,
                                           method,
                                           parameters_combinations,
                                           theta_0 = 0,
                                           add_baselines = FALSE,
                                           target_to_source_std_ratio = 1,
                                           source_denominator_change_factor = 1,
                                           design_prior = NULL) {
  # Filter on the selected case study, method, etc
  results_df <- input_df %>%
    dplyr::filter(
      case_study == !!case_study,
      source_denominator_change_factor == !!source_denominator_change_factor  |
        is.na(source_denominator_change_factor),
      target_to_source_std_ratio == !!target_to_source_std_ratio |
        is.na(target_to_source_std_ratio)
    )

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



  # results_df <- results_df %>%
  #   filter(!is.na(prepost_proba_TP) & !is.na(average_tie) & !is.na(average_power) & !is.na(prior_proba_success) & !is.na(prior_proba_no_benefit) & !is.na(prior_proba_benefit) & !is.na(prepost_proba_FP) & !is.na(upper_bound_proba_FP))

  # columns_to_check <- c("prepost_proba_TP", "average_tie", "average_power",
  #                       "prior_proba_success", "prior_proba_no_benefit",
  #                       "prior_proba_benefit", "prepost_proba_FP", "upper_bound_proba_FP")
  #
  # # Filter columns that exist in the data frame
  # existing_columns <- columns_to_check[columns_to_check %in% colnames(results_df)]
  #
  # # Dynamically filter the data frame using only the existing columns
  # results_df <- results_df %>%
  #   filter(if_all(all_of(existing_columns), ~ !is.na(.)))

  if (!(metric %in% colnames(results_df))){
    return()
  }
  results_df <- results_df[!is.na(results_df[,metric]), ]

  if (nrow(results_df) == 0) {
    warning("Dataframe is empty")
    return()
  }

  if (!(is.null(method))){
    results_df <- results_df %>%
      dplyr::filter(method == !!method)

    # Format the parameters (from the json string)
    results_df <- cbind(results_df, get_parameters(results_df[, "parameters"]))
    # Filter rows that match parameters_combinations
    results_df <- merge(results_df, parameters_combinations[, ])

    results_df[results_df[, "design_prior_type"] == "analysis_prior", "label"] <- "Analysis prior"
    results_df[results_df[, "design_prior_type"] == "source_posterior", "label"] <- "Source posterior"
    results_df[results_df[, "design_prior_type"] == "ui_design_prior", "label"] <- "UI prior"

    labels_title <- "Design prior"

    design_priors_names <- unique(results_df[, "label"])

    color_map <- viridis::viridis(length(design_priors_names))

    param_label <- make_labels_from_parameters(parameters_combinations, method)

    if (param_label == "") {
      param_label <- paste0("")
    } else {
      param_label <- paste0(", ", param_label)
    }

    figure_name <- paste0(
      case_study,
      "_",
      method,
      "_",
      selected_metric_name,
      "_vs_target_sample_size_per_arm_",
      convert_params_to_str(methods_dict[[method]], parameters_combinations)
    )

    plot_title <-
      sprintf(
        "%s, %s %s",
        str_to_title(case_study),
        methods_labels[[method]]$label,
        param_label
      )

  } else if (!(is.null(design_prior))) {

    labels_title <- "Methods"

    results_df <- results_df %>%
      dplyr::filter(design_prior_type == !!design_prior)

    if (design_prior == "analysis_prior") {
      design_prior_label <- "Analysis prior"
    } else if (design_prior == "source_posterior"){
      design_prior_label <- "Source posterior"
    } else if (design_prior == "ui_design_prior"){
      design_prior_label <- "UI design prior"
    } else {
      stop("Unknown design prior type.")
    }


    results_df[, "label"] <- design_prior_label

    color_map <- viridis::viridis(length(unique(results_df$method)))

    figure_name <- paste0(
      case_study,
      "_",
      selected_metric_name,
      "_vs_target_sample_size_per_arm_all_methods_",
      design_prior
    )

    plot_title <-
      sprintf(
        "%s, Design prior : %s",
        str_to_title(case_study),
        design_prior_label
      )

    boolean_filter <- c()
    # Sort by method before filtering on parameters
    results_df <- results_df %>% arrange(method)
    # Filter on the important parameters values for each method
    for (selected_method in unique(results_df$method)){
      parameters <- get_parameters(results_df[results_df$method == selected_method, "parameters", drop = FALSE])
      filter_per_method = rep(TRUE, nrow(parameters))
      for (key in names(methods_dict[[selected_method]])) {
        # Loop over the method's parameters
        # Get the important values for which we will make plots
        if (is.null(methods_dict[[selected_method]][[key]][["important_values"]])) {
          important_parameters_values <- methods_dict[[selected_method]][[key]][["range"]]
        } else {
          important_parameters_values <- methods_dict[[selected_method]][[key]][["important_values"]]
        }
        if (selected_method == "commensurate_power_prior" && key == "heterogeneity_prior"){
          key <- paste0(key, ".family")
        }
        filter <- parameters[,key] %in% important_parameters_values
        filter_per_method <- filter & filter_per_method
      }
      boolean_filter <- c(boolean_filter, filter_per_method)
    }

    results_df <- results_df[boolean_filter,]

    color_map <- viridis::viridis(n = 256)

    if (nrow(results_df) == 0){
      return()
    }
    results_df$rows <- seq(1, nrow(results_df))

    # Process the row of a results dataframe to create a Method + Parameters label
    parameters_labels <- sapply(1:nrow(results_df), function(i) {
      process_method_parameters_label(results_df[i, ], methods_labels)
    })

    results_df$parameters_labels <- lapply(1:nrow(results_df), function(i) {
      process_method_parameters_label(results_df[i, ], methods_labels)
    })

  } else {
    stop("Either method or design_prior must be specified.")
  }

  if (nrow(results_df) == 0) {
    warning("Dataframe is empty")
    return()
  }

  case_study_config <- yaml::read_yaml(paste0(case_studies_config_dir, case_study, ".yml"))

  # # Filter on important drift value
  # source_treatment_effect_estimate <- unique(results_df$source_treatment_effect_estimate)[1]
  # mandatory_drift_values <- important_drift_values(source_treatment_effect_estimate, case_study_config)
#
#   if ("drift" %in% colnames(results_df)) {
#     # Find the closest values in df$target_treatment_effect to important_target_treatment_effects
#     closest_values <- sapply(mandatory_drift_values, function(x) {
#       results_df$drift[which.min(abs(results_df$drift - x))]
#     })
#
#     # Filter the dataframe to keep only the rows with the closest values
#     results_df <- results_df %>%
#       dplyr::filter(drift %in% closest_values) %>%
#       dplyr::arrange(drift)
#
#     if (length(unique(results_df$drift)) != 3) {
#       stop("Only three main treatment effect values should be selected")
#     }
#   }

  # # Plot setup
  # color_map <- viridis::viridis(length(mandatory_drift_values))



  selected_label <- selected_metric$label

  # Calculate the range of the x-axis
  x_range <- range(results_df$target_sample_size_per_arm)
  x_width <- diff(x_range)
  # Set a relative cap size (necessary, otherwise the cap size will depend on the x-axis width)
  cap_size <- x_width * relative_error_cap_width # Adjust relative_error_cap_width to control the relative cap size


  plot_title <- format_title(title = plot_title, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)
  figure_name <- format_filename(filename = figure_name, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)


  if (!(is.null(method))){
    if (metric %in% c("prior_proba_benefit", "prior_proba_no_benefit")){
       # Only plot the prior probability of (no) benefit for the analysis design prior.
      results_df <- results_df %>%
        dplyr::filter(
          label == "Analysis prior"
        )

      # Produce the plot
      plt <- ggplot2::ggplot(
        results_df,
        ggplot2::aes(
          x = target_sample_size_per_arm,
          y = !!selected_metric_name
        )
      ) +
        geom_point() +
        # geom_point(position = position_dodge(width = cap_size / 2)) +
        geom_line() +
        scale_x_continuous(name = "Target study sample size per arm") +
        scale_y_continuous(name = selected_label)
    } else {
      # Produce the plot
      design_priors_names <- c("Analysis prior",  "Source posterior", "UI design prior")

      plt <- ggplot2::ggplot(
        results_df,
        ggplot2::aes(
          x = target_sample_size_per_arm,
          y = !!selected_metric_name,
          group = as.factor(label),
          color = as.factor(label)
        )
      ) +
        geom_point() +
        # geom_point(position = position_dodge(width = cap_size / 2)) +
        geom_line() +
        scale_x_continuous(name = "Target study sample size per arm") +
        scale_y_continuous(name = selected_label) +
        ggplot2::labs(color = labels_title) +
        scale_color_manual(values = color_map,
                           name = NULL,
                           labels = design_priors_names)
    }
  } else {
    if (metric %in% c("prior_proba_benefit", "prior_proba_no_benefit") & (design_prior != "analysis_prior")){
      return() # The prior proba of (no) benefit is only relevant for an analysis prior.
    }
    # Produce the plot
    design_priors_names <- c("Analysis prior",  "Source posterior", "UI design prior")

    parameters_labels <- setNames(object = parameters_labels, nm = results_df$parameters_labels)

    plt <- ggplot2::ggplot(
      results_df,
      ggplot2::aes(
        x = target_sample_size_per_arm,
        y = !!selected_metric_name,
        group = factor(parameters_labels, levels = unique(results_df$parameters_labels)),
        color = factor(parameters_labels, levels = unique(results_df$parameters_labels)),
        shape = factor(parameters_labels, levels = unique(results_df$parameters_labels))
      )
    )

    plt <- plt +
      geom_point(size = 1) + # Adjust size if needed
      geom_line() +
      scale_x_continuous(name = "Target study sample size per arm") +
      scale_y_continuous(name = selected_label) +
      ggplot2::labs(color = "Methods") +
      scale_color_viridis_d(labels = parameters_labels, guide = "legend") +
      scale_fill_viridis_d(labels = parameters_labels, guide = "legend") +
      scale_shape_manual(values = 1:length(parameters_labels), labels = parameters_labels, guide = "legend") +
      guides(color = guide_legend(title = "Methods"), shape = guide_legend(title = "Methods"))
  }

  plt <- plt +
    ggplot2::labs(title = plot_title)

  # Add a y-line at theta = 0
  if (add_baselines) {
    plt <- plt + geom_hline(yintercept = theta_0, linetype = "dashed")
  }

  # Save plot
  directory <- paste0(figures_dir, case_study)
  if (!dir.exists(directory)) {
    dir.create(directory, showWarnings = FALSE, recursive = TRUE)
  }

  plot.size <- set_size(textwidth)
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
bayesian_ocs_plots <- function(results_metrics_df, metrics) {
  # Get the list of case studies
  case_studies <- unique(results_metrics_df$case_study)

  design_prior_types <- c("analysis_prior", "source_posterior", "ui_design_prior")

  for (case_study in case_studies) {
    results_df0 <- results_metrics_df %>%
      dplyr::filter(
        case_study == !!case_study
      )

    conf <- yaml::yaml.load_file(system.file(paste0(
      "conf/case_studies/", case_study, ".yml"
    )))
    theta_0 <- conf$theta_0

    for (source_denominator_change_factor in unique(results_df0$source_denominator_change_factor)){
      results_df1 <- results_df0 %>%
        dplyr::filter(
          source_denominator_change_factor == !!source_denominator_change_factor  |
            is.na(source_denominator_change_factor)
        )

      for (target_to_source_std_ratio in unique(results_df1$target_to_source_std_ratio)){

        results_df2 <- results_df1 %>%
          dplyr::filter(
            target_to_source_std_ratio == !!target_to_source_std_ratio |
              is.na(target_to_source_std_ratio)
          )


        # Get the list of methods
        methods <- unique(results_df2$method)

        for (method in methods) {
          # Get the different parameters combinations studies for this method

          results_df3 <- results_df2 %>%
            dplyr::filter(method == !!method, )

          parameters_combinations <- unique(get_parameters(results_df3[, "parameters"]))

          if (method %in% c("separate", "pooling")) {
            add_baselines <- FALSE
          } else {
            add_baselines <- TRUE
          }

          target_sample_sizes_per_arm <- unique(results_df3$target_sample_size_per_arm)


          for (metric in names(metrics)) {
            for (sample_size_per_arm in target_sample_sizes_per_arm) {
              bayesian_metric_vs_parameters(
                results_metrics_df = results_df3,
                case_study = case_study,
                method = method,
                target_sample_size_per_arm = sample_size_per_arm,
                metric = metric,
                theta_0 = theta_0,
                source_denominator_change_factor = source_denominator_change_factor,
                target_to_source_std_ratio = target_to_source_std_ratio
              )
            }

            for (i in 1:nrow(parameters_combinations)) {
              #Plot the results with different design priors on the same plot, for a single method
              bayesian_metric_vs_sample_size(
                metric = metric,
                input_df = results_df3,
                case_study = case_study,
                method = method,
                parameters_combinations = parameters_combinations[i, , drop = FALSE],
                theta_0 = theta_0,
                add_baselines = add_baselines,
                source_denominator_change_factor = source_denominator_change_factor,
                target_to_source_std_ratio = target_to_source_std_ratio
              )
            }
          # }
        # }
          }
        }
        for (metric in names(metrics)) {
          for (design_prior in design_prior_types){
            # Plot the results for a single design priors for multiple methods on the same plot
            bayesian_metric_vs_sample_size(
              metric = metric,
              input_df = results_df2,
              case_study = case_study,
              method = NULL,
              parameters_combinations = NULL,
              theta_0 = theta_0,
              add_baselines = add_baselines,
              source_denominator_change_factor = source_denominator_change_factor,
              target_to_source_std_ratio = target_to_source_std_ratio,
              design_prior = design_prior
            )
          }
        }
      }
    }
  }
}
