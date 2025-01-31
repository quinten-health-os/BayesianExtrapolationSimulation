library(dplyr)
library(ggplot2)

plot_effect_denominator_change = function(results_df, metric, figures_dir, treatment_effect = NULL, design_prior = NULL){

  if (nrow(results_df)<2){
    return()
  }

  if (length(unique(results_df$source_denominator_change_factor)) < 2){
    return()
  }

  results_df <- results_df %>%
    dplyr::group_by(parameters) %>%
    dplyr::filter(n_distinct(source_denominator_change_factor) > 1) %>%
    dplyr::ungroup()

  if (nrow(results_df) == 0){
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

  if (any(metric %in% names(bayesian_metrics))) {

    if (design_prior == "analysis_prior") {
      design_prior_label <- "Analysis prior"
    } else if (design_prior == "source_posterior"){
      design_prior_label <- "Source posterior"
    } else if (design_prior == "ui_design_prior"){
      design_prior_label <- "UI design prior"
    } else {
      stop("Unknown design prior type.")
    }

    results_df <- results_df %>%
      dplyr::filter(design_prior_type == design_prior)
  }

  if (nrow(results_df) == 0){
    return()
  }

  # Process the row of a results dataframe to create a Method + Parameters label
  parameters_labels <- sapply(1:nrow(results_df), function(i) {
    process_method_parameters_label(results_df[i, ], methods_labels)
  })
  # Process the rows to create parameter labels
  results_df$parameters_labels <- lapply(1:nrow(results_df), function(i) {
    process_method_parameters_label(results_df[i, ], methods_labels, method_name = TRUE)
  })

  results_df$methods_parameters_labels <- lapply(1:nrow(results_df), function(i) {
    process_method_parameters_label(results_df[i, ], methods_labels, method_name = TRUE, as_latex = FALSE)
  })

  parameters_labels <- setNames(object = results_df$parameters_labels, nm = results_df$parameters_labels)
  methods <- format_results_df_methods(results_df)
  results_df$method <- methods
  methods <- unique(methods)

  target_sample_size_per_arm <- unique(results_df$target_sample_size_per_arm)

  if (any(metric %in% names(frequentist_metrics)) || any(metric %in% names(inference_metrics))) {
    treatment_effects_labels <- list(
      no_effect = "No effect",
      partially_consistent = "Partially consistent effect",
      consistent = "Consistent effect"
    )

    if (treatment_effect == "consistent") {
      results_df <- results_df %>%
        dplyr::filter(target_treatment_effect == source_treatment_effect_estimate)
    } else if (treatment_effect == "no_effect") {
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

    if (nrow(results_df) == 0){
      return()
    }

    plot_title <- sprintf(
      "%s, $N_T/2 = $%s, %s",
      str_to_title(unique(results_df$case_study)),
      target_sample_size_per_arm,
      treatment_effects_labels[[treatment_effect]]
    )

    filename <- paste0(
      unique(results_df$case_study),
      "_",
      metric$name,
      "_vs_source_denominator_change_factor",
      "_treatment_effect=",
      treatment_effect,
      "_target_sample_size_per_arm=",target_sample_size_per_arm
    )

  } else if (any(metric %in% names(bayesian_metrics))) {


    if (design_prior == "analysis_prior") {
      design_prior_label <- "Analysis prior"
    } else if (design_prior == "source_posterior"){
      design_prior_label <- "Source posterior"
    } else if (design_prior == "ui_design_prior"){
      design_prior_label <- "UI design prior"
    } else {
      stop("Unknown design prior type.")
    }


    plot_title <- sprintf(
      "%s, $N_T/2 = $%s, %s",
      str_to_title(unique(results_df$case_study)),
      target_sample_size_per_arm,
      design_prior_label
    )

    filename <- paste0(
      unique(results_df$case_study),
      "_",
      metric$name,
      "_vs_source_denominator_change_factor",
      "_target_sample_size_per_arm=",target_sample_size_per_arm,
      design_prior_label
    )
  }


  plot_title <- format_title(
    title = plot_title,
    case_study = unique(results_df$case_study))

  filename <- format_filename(
    filename = filename,
    case_study = unique(results_df$case_study))

  directory <- paste0(figures_dir, unique(results_df$case_study))
  file_path <- file.path(directory, filename)

  if (file.exists(paste0(file_path, ".pdf")) && file.exists(paste0(file_path, ".png")) && remake_figures == FALSE){
    return()
  }

  # Calculate the range of the x-axis
  x_range <- range(results_df$source_denominator_change_factor)
  x_width <- diff(x_range)
  cap_size <- x_width * relative_error_cap_width # Adjust relative_error_cap_width to control the relative cap size

  labels_title <- "Methods"


  selected_metric_uncertainty_lower <- rlang::sym(paste0(metric$metric_uncertainty, "_lower"))
  selected_metric_uncertainty_upper <- rlang::sym(paste0(metric$metric_uncertainty, "_upper"))


  if (c(selected_metric_uncertainty_lower) %in% colnames(results_df) && c(selected_metric_uncertainty_upper) %in% colnames(results_df)){
    plt <- ggplot2::ggplot()
    plt <- plot_element_metric_vs_xvar_methods(
      plt,
      results_df,
      sym("source_denominator_change_factor"),
      sym(metric$name),
      selected_metric_uncertainty_lower,
      selected_metric_uncertainty_upper,
      cap_size,
      markersize,
      labels_title,
      join_points = TRUE,
      color,
      dodging = FALSE,
      parameters_labels =  parameters_labels
    )
  } else {
    plt <- ggplot2::ggplot()
    plt <- plot_element_metric_vs_xvar_methods(
      plt,
      results_df,
      sym("source_denominator_change_factor"),
      sym(metric$name),
      NULL,
      NULL,
      cap_size,
      markersize,
      labels_title,
      join_points = TRUE,
      color,
      dodging = FALSE,
      parameters_labels =  parameters_labels
    )
  }

  plt <- plt +
    labs(
      title = plot_title,
      x = "Source Denominator Change Factor",
      y = metric$label
    ) +
    theme_minimal()


  plot.size <- set_size(textwidth)
  fig_width_in <- plot.size[1]
  fig_height_in <- plot.size[2]

  # Export the plots
  export_plots(plt, file_path, fig_width_in, fig_height_in, type = "pdf")
  export_plots(plt, file_path, fig_width_in, fig_height_in, type = "png")

  return()
}

loop_plot_effect_denominator_change = function(results_df, metrics){
  # Define the columns that uniquely define a scenario
  scenario_vars <- c("target_sample_size_per_arm",
                     "drift", "case_study") #"method", "parameters",
  existing_vars <- intersect(scenario_vars, colnames(results_df))

  # Get a data frame of unique scenarios
  unique_scenarios <- results_df %>%
    distinct(across(all_of(existing_vars)))

  unique_scenarios  <- subset(results_df,  case_study %in% c("teriflunomide", "belimumab", "mepolizumab"))

  main_treatment_effects <-  list(
    no_effect = "no_effect",
    partially_consistent = "partially_consistent",
    consistent = "consistent"
  )

  # Loop over each scenario
  for (i in seq_len(nrow(unique_scenarios))) {
    # Extract the scenario as a single row tibble
    scenario <- unique_scenarios[i, ]

    # Filter the main data frame for this scenario

    scenario_df <- results_df %>%
      filter(
        #method == scenario$method,
        #parameters == scenario$parameters,
        target_sample_size_per_arm == scenario$target_sample_size_per_arm,
        if ("drift" %in% colnames(results_df)) drift == scenario$drift else TRUE,
        case_study == scenario$case_study
      )

    for (metric in metrics){
      if (metric %in% names(bayesian_metrics)){
        for (design_prior in c('analysis_prior', 'source_posterior', 'ui_design_prior')){
          plot_effect_target_denominator_change(scenario_df, metric, figures_dir, design_prior = design_prior)
        }
      } else {
        for (treatment_effect in main_treatment_effects){
          plot_effect_denominator_change(scenario_df, metric, figures_dir, treatment_effect = treatment_effect)
        }
      }
    }
  }
}


