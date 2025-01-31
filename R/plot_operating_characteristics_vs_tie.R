#' Function to generate an operating characteristic vs tie plot
#'
#' @param results_metrics_df The dataframe containing the results and metrics
#' @param case_study The case study name
#' @param target_sample_size_per_arm The target sample size per arm
#' @param treatment_effect The treatment effect type ("consistent", "no_effect", "partially_consistent")
#' @param operating_characteristic The operating characteristic to plot (e.g., "power", "type_1_error")
#' @param power_difference Logical indicating whether to calculate difference for the operating characteristic
#'
#' @return None
#'
#' @export
operating_characteristic_vs_tie <- function(
    results_metrics_df,
    case_study,
    target_sample_size_per_arm,
    treatment_effect,
    operating_characteristic,
    source_denominator_change_factor,
    target_to_source_std_ratio
) {

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

  theta_0 <- unique(results_df$theta_0)
  results_df_tie <- results_df %>% dplyr::filter(target_treatment_effect == theta_0)

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

  # Process the row of a results dataframe to create a Method + Parameters label
  results_df$parameters_labels <- lapply(1:nrow(results_df), function(i) {
    process_method_parameters_label(results_df[i, ], methods_labels, method_name = FALSE)
  })

  results_df$parameters_labels_not_latex <- lapply(1:nrow(results_df), function(i) {
    process_method_parameters_label(results_df[i, ], methods_labels, method_name = FALSE, as_latex = FALSE)
  })

  results_df$methods_parameters_labels <- lapply(1:nrow(results_df), function(i) {
    process_method_parameters_label(results_df[i, ], methods_labels, method_name = TRUE, as_latex = FALSE)
  })

  parameters_labels <- setNames(object = results_df$parameters_labels, nm = results_df$parameters_labels)
  methods <- format_results_df_methods(results_df)
  results_df$method <- methods
  methods <- unique(methods)

  filename <- paste0(
    case_study,
    "_",
    operating_characteristic$name,
    "_vs_tie_sample_size=",
    target_sample_size_per_arm,
    "_treatment_effect=",
    treatment_effect
  )

  filename <- format_filename(
    filename = filename,
    case_study = case_study,
    target_to_source_std_ratio = target_to_source_std_ratio,
    source_denominator_change_factor = source_denominator_change_factor
  )
  file_path <- file.path(directory, filename)

  if (file.exists(paste0(file_path, ".pdf")) && file.exists(paste0(file_path, ".png")) && remake_figures == FALSE) {
    return()
  }

  y_range <- sort(range(results_df[, operating_characteristic$name]))
  ymin <- y_range[1]
  ymax <- y_range[2]

  x_range <- sort(range(results_df[, "tie"]))
  xmin <- x_range[1]
  xmax <- x_range[2]

  x_width <- xmax - xmin
  y_width <- ymax - ymin

  y_cap_size <- x_width * relative_error_cap_width / 2
  x_cap_size <- y_width * relative_error_cap_width / 2
  cap_size <- min(x_cap_size, y_cap_size)
  x_cap_size <- cap_size
  y_cap_size <- cap_size

  treatment_effects_labels <- list(
    no_effect = "No effect",
    partially_consistent = "Partially consistent effect",
    consistent = "Consistent effect"
  )

  plot_title <- sprintf(
    "%s, $N_T/2 = $%s, %s",
    str_to_title(case_study),
    target_sample_size_per_arm,
    treatment_effects_labels[[treatment_effect]]
  )

  plot_title <- format_title(
    title = plot_title,
    case_study = case_study,
    target_to_source_std_ratio = target_to_source_std_ratio,
    source_denominator_change_factor = source_denominator_change_factor
  )


  # Split data by 'method'
  results_df_split <- split(results_df, results_df$method)

  # Unique parameters and methods
  unique_parameters <- unique(results_df$parameters_labels)
  unique_methods <- sort(unique(results_df$method))

  # Shapes for methods
  shape_codes <- c(16, 2, 15, 18, 1, 8, 4, 3, 7, 9)

  shapes <- setNames(shape_codes[1:length(unique_methods)], unique_methods)

  # Initialize the plot
  plt <- ggplot(mapping = aes(x = tie, y = !!sym(operating_characteristic$name)))

  # Iterate over each method and add layers
  plt <- plt + purrr::imap(results_df_split, function(df_subset, method_name) {
    df_subset <- df_subset[!is.na(df_subset[[operating_characteristic$name]]), ]

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
        aes(
          x = tie,
          y = !!sym(operating_characteristic$name),
          color = factor(parameters_labels, levels = method_parameters)
        ),
        shape = method_shape,
        size = 2
      ),
      # Add vertical error bars
      geom_errorbar(
        data = df_subset,
        aes(
          x = tie,
          ymin = !!sym(paste0("conf_int_", operating_characteristic$name, "_lower")),
          ymax = !!sym(paste0("conf_int_", operating_characteristic$name, "_upper")),
          color = factor(parameters_labels, levels = method_parameters)
        ),
        width = y_cap_size
      ),
      # Add horizontal error bars
      geom_errorbarh(
        data = df_subset,
        aes(
          y = !!sym(operating_characteristic$name),
          xmin = conf_int_tie_lower,
          xmax = conf_int_tie_upper,
          color = factor(parameters_labels, levels = method_parameters)
        ),
        height = x_cap_size
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

  # Add vertical line
  plt <- plt +
    ggplot2::geom_vline(
      xintercept = analysis_config$nominal_tie,
      color = "black",
      linetype = "dashed"
    ) +
    ggplot2::scale_x_continuous(
      breaks = function(x) unique(c(pretty(x), analysis_config$nominal_tie))
    )

  # Add shape scale for methods
  plt <- plt +
    scale_shape_manual(
      values = shapes,
      name = "Methods",
      guide = guide_legend(order = 1, nrow = 2)
    ) +
    labs(
      title = plot_title,
      x = "TIE",
      y = operating_characteristic$label
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
    legend.key.size = unit(0.1, "cm"),
    legend.spacing.x = unit(0.03, "cm"),
    legend.spacing.y = unit(0.01, "cm"),
    legend.position = "bottom",
    legend.direction = "vertical",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.box = "horizontal",
    legend.box.just = "center"
  )

  plot.size <- set_size(textwidth)
  fig_width_in <- plot.size[1] * 1.35
  fig_height_in <- plot.size[2] * 1.5

  # Export the plots
  export_plots(plt, file_path, fig_width_in, fig_height_in, type = "pdf", adjust_theme = FALSE)
  export_plots(plt, file_path, fig_width_in, fig_height_in, type = "png", adjust_theme = FALSE)
}


#' Plot methods' operating characteristics
#'
#' @description OCs vs TIE
#'
#' @param results_metrics_df The data frame containing the results and metrics.
#' @param metrics A list of metrics to be plotted.
#'
#' @return None
#'
#' @examples NA
operating_characteristics_vs_tie_plots <- function(results_metrics_df, metrics) {
  # Get the list of case studies
  case_studies <- unique(results_metrics_df$case_study)

  treatment_effects <- c("partially_consistent", "consistent" ) #, "no_effect"


  for (case_study in case_studies) {
    results_metrics_df1 <- results_metrics_df[results_metrics_df$case_study == case_study,]

    # Get the list of target sample sizes
    target_sample_sizes <- unique(results_metrics_df1$target_sample_size_per_arm)

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

          for (source_denominator_change_factor in source_denominator_change_factors) {
            results_metrics_df4 <- results_metrics_df3 %>%
              dplyr::filter(
                source_denominator_change_factor == !!source_denominator_change_factor |
                  is.na(source_denominator_change_factor)
              )

            for (operating_characteristic in metrics){
              operating_characteristic_vs_tie(
                results_metrics_df4,
                case_study,
                target_sample_size,
                treatment_effect = treatment_effect,
                operating_characteristic = operating_characteristic,
                source_denominator_change_factor = source_denominator_change_factor,
                target_to_source_std_ratio = target_to_source_std_ratio
              )
            }
          }
        }
      }
    }
  }
}

bayesian_operating_characteristic_vs_tie <- function(results_metrics_df,
                                                     case_study,
                                                     target_sample_size_per_arm,
                                                     operating_characteristic,
                                                     source_denominator_change_factor,
                                                     target_to_source_std_ratio,
                                                     tie_type,
                                                     design_prior_type,
                                                     design_prior_type_tie) {

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

  if (tie_type == "average_tie") {
    results_df_tie <- results_metrics_df %>%
      dplyr::filter(
        design_prior_type == !!design_prior_type_tie
      )

    results_df <- results_df %>% select(-average_tie)
    results_df <- results_df %>%
      dplyr::filter(
        design_prior_type == !!design_prior_type
      )

    results_df$average_tie <- results_df_tie$average_tie

    if (design_prior_type_tie == "ui_design_prior") {
      design_prior_tie_label <- "UI design prior"
    } else if (design_prior_type_tie == "analysis_prior") {
      design_prior_tie_label <- "Analysis prior"
    } else if (design_prior_type_tie == "source_posterior") {
      design_prior_tie_label <- "Source posterior as design prior"
    } else {
      stop("Design prior type not supported.")
    }
  }

  # Process the rows to create parameter labels
  results_df$parameters_labels <- lapply(1:nrow(results_df), function(i) {
    process_method_parameters_label(results_df[i, ], methods_labels, method_name = FALSE)
  })

  results_df$parameters_labels_not_latex <- lapply(1:nrow(results_df), function(i) {
    process_method_parameters_label(results_df[i, ], methods_labels, method_name = FALSE, as_latex = FALSE)
  })

  results_df$methods_parameters_labels <- lapply(1:nrow(results_df), function(i) {
    process_method_parameters_label(results_df[i, ], methods_labels, method_name = TRUE, as_latex = FALSE)
  })

  parameters_labels <- setNames(object = results_df$parameters_labels, nm = results_df$parameters_labels)
  methods <- format_results_df_methods(results_df)
  results_df$method <- methods
  methods <- unique(methods)

  if (tie_type == "average_tie") {
    filename <- paste0(
      case_study,
      "_", design_prior_type, "_",
      operating_characteristic$name,
      "_vs_", design_prior_type_tie, "_", tie_type, "_sample_size=",
      target_sample_size_per_arm
    )
  } else {
    filename <- paste0(
      case_study,
      "_",
      operating_characteristic$name,
      "_vs_", tie_type, "_sample_size=",
      target_sample_size_per_arm
    )
  }

  filename <- format_filename(
    filename = filename,
    case_study = case_study,
    target_to_source_std_ratio = target_to_source_std_ratio,
    source_denominator_change_factor = source_denominator_change_factor
  )
  file_path <- file.path(directory, filename)

  if (file.exists(paste0(file_path, ".pdf")) && file.exists(paste0(file_path, ".png")) && remake_figures == FALSE) {
    return()
  }

  if (!(operating_characteristic$name %in% colnames(results_df))) {
    return()
  }

  y_range <- sort(range(results_df[, operating_characteristic$name]))
  ymin <- y_range[1]
  ymax <- y_range[2]

  x_range <- sort(range(results_df[, tie_type]))
  xmin <- x_range[1]
  xmax <- x_range[2]

  x_width <- xmax - xmin
  y_width <- ymax - ymin

  y_cap_size <- x_width * relative_error_cap_width / 2
  x_cap_size <- y_width * relative_error_cap_width / 2
  cap_size <- min(x_cap_size, y_cap_size)
  x_cap_size <- cap_size
  y_cap_size <- cap_size

  if (design_prior_type == "ui_design_prior") {
    design_prior_label <- "UI design prior"
  } else if (design_prior_type == "analysis_prior") {
    design_prior_label <- "Analysis prior"
  } else if (design_prior_type == "source_posterior") {
    design_prior_label <- "Source posterior as design prior"
  } else {
    stop("Design prior type not supported.")
  }

  plot_title <- sprintf(
    "%s, $N_T/2 = $%s, %s",
    str_to_title(case_study),
    target_sample_size_per_arm,
    design_prior_label
  )

  plot_title <- format_title(
    title = plot_title,
    case_study = case_study,
    target_to_source_std_ratio = target_to_source_std_ratio,
    source_denominator_change_factor = source_denominator_change_factor
  )

  # Load necessary libraries for plotting
  library(ggnewscale)
  library(purrr)

  # Split data by 'method'
  results_df_split <- split(results_df, results_df$method)

  # Unique parameters and methods
  unique_parameters <- unique(results_df$parameters_labels)
  unique_methods <- sort(unique(results_df$method))

  # Shapes for methods
  shape_codes <- c(16, 2, 15, 18, 1, 8, 4, 3, 7, 9)

  shapes <- setNames(shape_codes[1:length(unique_methods)], unique_methods)

  # Determine x-axis variable and label
  if (tie_type == "tie") {
    x_var <- "tie"
    x_label <- "TIE"
  } else {
    x_var <- "average_tie"
    x_label <- paste0("Average TIE with ", design_prior_tie_label)
  }

  if (all(is.na(results_df[[x_var]])) || all(is.na(results_df[[operating_characteristic$name]]))){
    return()
  }

  # Initialize the plot
  plt <- ggplot(mapping = aes_string(x = x_var, y = operating_characteristic$name))

  # Iterate over each method and add layers
  plt <- plt + purrr::imap(results_df_split, function(df_subset, method_name) {
    df_subset <- df_subset[!is.na(df_subset[[operating_characteristic$name]]), ]

    # Get unique parameters for this method
    method_parameters <- unique(df_subset$parameters_labels)

    if (length(method_parameters) == 0){
      return()
    }

    # # Generate a large color palette
    # full_palette <- scales::hue_pal()(50)  # Generate a larger palette
    # # Randomly sample 'length(method_parameters)' colors from the full palette
    # method_colors <- setNames(
    #   sample(full_palette, length(method_parameters), replace = FALSE),
    #   method_parameters
    # )
    #
    method_colors <- setNames(
        scales::hue_pal()(length(method_parameters)),
        method_parameters
      )

    # Get the shape for this method
    method_shape <- shapes[method_name]

    # Create the plot layers for this method
    layers <- list(
      # Add points for this method
      geom_point(
        data = df_subset,
        aes(color = factor(parameters_labels, levels = method_parameters)),
        shape = method_shape,
        size = 2
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


    # Add horizontal error bars if tie_type is "tie"
    if (tie_type == "tie") {
      layers <- c(layers, list(
        geom_errorbarh(
          data = df_subset,
          aes(
            xmin = conf_int_tie_lower,
            xmax = conf_int_tie_upper,
            y = !!sym(operating_characteristic$name),
            color = factor(parameters_labels, levels = method_parameters)
          ),
          height = x_cap_size
        )
      ))
    }

    # Reset the color scale for the next method
    layers <- c(layers, list(new_scale_color()))

    layers
  })

  # Add shape scale for methods and labels
  plt <- plt +
    scale_shape_manual(
      values = shapes,
      name = "Methods",
      guide = guide_legend(order = 1, nrow = 2)
    ) +
    labs(
      title = plot_title,
      x = x_label,
      y = operating_characteristic$label
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


  # Export the plots
  export_plots(plt, file_path, fig_width_in, fig_height_in, type = "pdf", adjust_theme = FALSE)
  export_plots(plt, file_path, fig_width_in, fig_height_in, type = "png", adjust_theme = FALSE)

}


#' Plot methods operating characteristics
#'
#' @description OCs vs TIE.
#'
#' @param results_metrics_df The data frame containing the results and metrics.
#' @param metrics A list of metrics to be plotted.
#'
#' @return None
#'
#' @examples NA
bayesian_operating_characteristics_vs_tie_plots <- function(results_bayes_df, results_freq_df, bayesian_metrics) {
 common_columns <- intersect(names(results_bayes_df), names(results_freq_df))

  # Select the results that correspond to TIE computation.
  results_freq_df_tie <- results_freq_df[results_freq_df$target_treatment_effect == results_freq_df$theta_0, c(
    common_columns,
    c(
      "tie",
      "conf_int_tie_lower",
      "conf_int_tie_upper"
    )
  )]

  common_columns <- intersect(names(results_freq_df_tie),common_columns)

  # Now merge results_bayes_df with the filtered results_freq_selected on the common columns specified in 'by'
  results_bayes_df_merged <- results_bayes_df %>%
  left_join(results_freq_df_tie, by =  common_columns)

  # Get the list of case studies
  case_studies <- unique(results_bayes_df_merged$case_study)

  for (case_study in case_studies) {
    results_metrics_df1 <- results_bayes_df_merged[results_bayes_df_merged$case_study == case_study,]

    # Get the list of target sample sizes
    target_sample_sizes <- unique(results_metrics_df1$target_sample_size_per_arm)

    for (target_sample_size in target_sample_sizes) {
      results_metrics_df2 <- results_metrics_df1[results_metrics_df1$target_sample_size_per_arm == target_sample_size,]

      target_to_source_std_ratios <- unique(results_metrics_df2$target_to_source_std_ratio)

      for (target_to_source_std_ratio in target_to_source_std_ratios) {
        results_metrics_df3 <- results_metrics_df2 %>%
          dplyr::filter(
            target_to_source_std_ratio == !!target_to_source_std_ratio |
              is.na(target_to_source_std_ratio)
          )

        source_denominator_change_factors <- unique(results_metrics_df3$source_denominator_change_factor)

        for (source_denominator_change_factor in source_denominator_change_factors) {
          results_metrics_df4 <- results_metrics_df3 %>%
            dplyr::filter(
              source_denominator_change_factor == !!source_denominator_change_factor |
                is.na(source_denominator_change_factor)
            )

          design_prior_types <- unique(results_metrics_df4$design_prior_type)

          #design_prior_types <- "analysis_prior" # FIXME
          for (operating_characteristic in bayesian_metrics){
            for (design_prior_type in design_prior_types){
              bayesian_operating_characteristic_vs_tie(
                results_metrics_df4,
                case_study,
                target_sample_size,
                operating_characteristic = operating_characteristic,
                source_denominator_change_factor = source_denominator_change_factor,
                target_to_source_std_ratio = target_to_source_std_ratio,
                tie_type = "tie",
                design_prior_type = design_prior_type
              )

              for (design_prior_type_tie in design_prior_types){
                bayesian_operating_characteristic_vs_tie(
                  results_metrics_df4,
                  case_study,
                  target_sample_size,
                  operating_characteristic = operating_characteristic,
                  source_denominator_change_factor = source_denominator_change_factor,
                  target_to_source_std_ratio = target_to_source_std_ratio,
                  tie_type = "average_tie",
                  design_prior_type = design_prior_type,
                  design_prior_type_tie = design_prior_type_tie
                )
              }
            }
          }
        }
      }
    }
  }
}
