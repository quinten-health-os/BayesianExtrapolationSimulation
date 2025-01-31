#' Generate a metric vs drift table
#'
#' @param metric The metric to be plotted
#' @param results_metrics_df A dataframe containing the results metrics
#' @param theta_0 The null hypothesis value
#' @param case_study The case study being analyzed
#' @param method The method being used
#' @param category The category of analysis ('parameters' or 'target_sample_size_per_arm')
#' @param control_drift Boolean indicating whether to control for drift
#' @param target_sample_size_per_arm The target sample size per arm
#' @param parameters_combinations The combinations of parameters
#' @param xvars A list containing x-axis variable information
#'
#' @return Generates and saves tables in HTML, PDF, and LaTeX formats
#' @export
#'
#' @importFrom dplyr filter mutate select arrange rename group_by ungroup
#' @importFrom rlang sym
#' @importFrom knitr kable
#' @importFrom kableExtra kable_styling row_spec collapse_rows
table_metric_vs_drift <- function(metric,
                                  results_metrics_df,
                                  theta_0,
                                  case_study,
                                  method,
                                  category,
                                  control_drift,
                                  target_sample_size_per_arm,
                                  parameters_combinations,
                                  xvars,
                                  source_denominator_change_factor = 1,
                                  target_to_source_std_ratio = 1,
                                  wide_table = TRUE) {
  # Filter results by case study and target sample size per arm
  results_metrics_df <- results_metrics_df %>%
    dplyr::filter(case_study == !!case_study) %>%
    dplyr::filter(method == !!method)


  if (control_drift) {
    results_metrics_df <- results_metrics_df %>% dplyr::filter(treatment_drift == 0)
    xvar <- xvars$control_drift
  } else {
    results_metrics_df <- results_metrics_df %>% dplyr::filter(control_drift == 0)
    xvar <- xvars$drift
  }

  results_df <- results_metrics_df %>% dplyr::filter(method == !!method)

  categories <- unique(results_df[[category]])

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

    # Make sure that the data are sorted according to the parameters values
    params <- get_parameters(results_df[, "parameters"])
    column_names <- colnames(params)
    sorted_index <- do.call(order, lapply(params[column_names], as.factor))
    # Reorder the data :
    results_df[results_df$method == method, ] <- results_df[sorted_index, ]
  } else if (category == "source_denominator") {
    results_df <- results_df %>%
      dplyr::filter(
        target_sample_size_per_arm == !!target_sample_size_per_arm,
        target_to_source_std_ratio == !!target_to_source_std_ratio |
          is.na(target_to_source_std_ratio),
      )
  } else if (category == "target_to_source_std_ratio") {
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

  results_df$rows <- seq(1, nrow(results_df))

  parameters_labels <- NA
  if (category == "target_sample_size_per_arm") {
    results_df$label <- results_df$target_sample_size_per_arm
    labels_title <- "Target sample size per arm"
  } else if (category == "parameters") {
    parameters_labels <- unname(sapply(1:nrow(results_df), function(i) {
      process_method_parameters_label(results_df[i, ], methods_labels, method_name = FALSE)
    }))

    labels_title <- "Parameters"
  } else if (category == "source_denominator") {
    results_df$label <- results_df$source_denominator
    labels_title <- "Source denominator change factor"
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
    next
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

  labels <- list()
  for (i in seq_along(categories)) {
    combination <- categories[i]
    if (category == "source_denominator") {
      if (is.na(combination)) {
        next
      } else {
        results_df[results_df[, "source_denominator"] == combination, "label"] <- combination
      }
    }
  }


  if (!(category == "parameters")){
    parameters_label <- make_labels_from_parameters(parameters_combinations, method)
    parameters_str <- convert_params_to_str(methods_dict[[method]], parameters_combinations)

    if (parameters_label == "") {
      parameters_label_title <- ""
    } else {
      parameters_label_title <- paste0(", ", parameters_label)
    }

    results_df$labels <- parameters_labels
  }


  if (category == "target_sample_size_per_arm") {
    title <- sprintf(
      "%s, %s, %s %s",
      frequentist_metrics[[metric]]$label,
      str_to_title(case_study),
      methods_labels[[method]]$label,
      parameters_label_title
    )

    filename <- paste0(
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
  } else if (category == "source_denominator") {
    title <- sprintf(
      "%s, %s, %s %s, %s %s",
      frequentist_metrics[[metric]]$label,
      str_to_title(case_study),
      methods_labels[[method]]$label,
      parameters_label_title,
      "$N_T/2 = $",
      target_sample_size_per_arm
    )
    filename <- paste0(
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
      "%s, %s, %s, $N_T/2 = $ %s",
      frequentist_metrics[[metric]]$label,
      str_to_title(case_study),
      methods_labels[[method]]$label,
      target_sample_size_per_arm
    )
    filename <- paste0(
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
      "%s, %s, %s %s, %s %s",
      frequentist_metrics[[metric]]$label,
      str_to_title(case_study),
      methods_labels[[method]]$label,
      parameters_label_title,
      "$N_T/2 = $",
      target_sample_size_per_arm
    )
    filename <- paste0(
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

  title <- format_title(title = title, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor, as_latex = TRUE)
  filename <- format_filename(filename = filename, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)


  if (category == "target_sample_size_per_arm") {
    results_df$label <- results_df$target_sample_size_per_arm
    labels_title <- "Target sample size per arm"
  } else if (category == "parameters") {
    results_df$label <- format_results_df_parameters(results_df, include_method_name = FALSE)
    labels_title <- "Parameters"
  } else if (category == "source_denominator") {
    results_df$label <- results_df$source_denominator
    labels_title <- "Source denominator change factor"
  } else if (category == "target_to_source_std_ratio") {
    results_df$label <- results_df$target_to_source_std_ratio
    labels_title <- expression(Tex("$\\sigma_T/\\sigma_S")) # Ratio between target and source standard deviation
  } else {
    stop("Not implemented for this category")
  }

  data_table <- results_df %>%
    dplyr::select(
      label,
      !!xvar_name,
      !!selected_metric_name,
      !!selected_metric_uncertainty_lower,
      !!selected_metric_uncertainty_upper
    ) %>%
    dplyr::mutate(
      !!xvar_label := format_num(!!xvar_name),
      !!selected_metric$label := format_num(!!selected_metric_name),
      error_low = format_num(!!selected_metric_uncertainty_lower),
      error_upper = format_num(!!selected_metric_uncertainty_upper),
      !!selected_metric$uncertainty_label := paste0("[", error_low, ", ", error_upper, "]")
    ) %>%
    dplyr::select(
      label,
      !!xvar_label,
      !!selected_metric$label,
      !!selected_metric$uncertainty_label
    ) %>%
    dplyr::rename(
      "$\\delta$" = "Drift in treatment effect",
    )


  # Merge Mean and CI columns
  data_table <- data_table %>%
    unite(!!selected_metric$label, selected_metric$label, selected_metric$uncertainty_label, sep = " ")
  if (wide_table && !(all(is.na(data_table$label) | data_table$label == ""))){
    # Pivot the table so that Parameters become columns and Sample Size per Arm is in rows
    data_table_formatted <- data_table %>%
      tidyr::pivot_wider(
        names_from = label,  # The column to pivot (label)
        values_from = !!rlang::sym(selected_metric$label)  # The values for the new columns
      )
    longtable = FALSE
  } else {
    data_table_formatted <- data_table
    longtable = TRUE
  }


  directory <- file.path(tables_dir, case_study)
  if (!dir.exists(directory)) {
    dir.create(directory, showWarnings = FALSE, recursive = TRUE)
  }

  file_path <- file.path(directory, filename)


  export_table(data_table = data_table_formatted, ncollapses = 1, title = title, file_path = file_path)
}

#' Plot metric vs parameters
#'
#' @param results_metrics_df A dataframe containing the results metrics
#' @param metric The metric to be plotted
#' @param case_study The case study being analyzed
#' @param method The method being used
#' @param target_sample_size_per_arm The target sample size per arm
#' @param theta_0 The null hypothesis value
#'
#' @return Generates and saves tables in HTML, PDF, and LaTeX formats
#' @export
#'
#' @importFrom dplyr filter mutate select arrange rename group_by ungroup
#' @importFrom rlang sym
#' @importFrom knitr kable
#' @importFrom kableExtra kable_styling row_spec collapse_rows
table_metric_vs_parameters <- function(results_metrics_df,
                                       metric,
                                       case_study = "belimumab",
                                       method = "RMP",
                                       target_sample_size_per_arm = 93,
                                       theta_0 = NULL,
                                       source_denominator_change_factor = 1,
                                       target_to_source_std_ratio = 1) {
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

  if (nrow(results_df) == 0) {
    warning("Dataframe is empty")
    next
  }

  case_study_config <- yaml::read_yaml(paste0(case_studies_config_dir, case_study, ".yml"))

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

      results_subdf$parameter_values <- parameter_values

      data_table <- results_subdf %>%
        dplyr::mutate(
          means = format_num(!!rlang::sym(selected_metric$name)),
          error_low = format_num(!!rlang::sym(
            paste0(selected_metric$metric_uncertainty, "_lower")
          )),
          error_upper = format_num(!!rlang::sym(
            paste0(selected_metric$metric_uncertainty, "_upper")
          )),
          CI = paste0("[", error_low, ", ", error_upper, "]")
        ) %>%
        dplyr::select(parameter_values, means, CI) %>%
        dplyr::arrange(parameter_values)

      data_table <- data_table %>%
        # dplyr::group_by(parameter_values) %>%
        dplyr::rename(
          "Parameter" = parameter_values,
          !!selected_metric$label := means,
          !!selected_metric$uncertainty_label := CI
        )

      data_table <- data_table %>%
        unite(!!selected_metric$label, selected_metric$label, selected_metric$uncertainty_label, sep = " ")


      directory <- file.path(tables_dir, case_study)
      if (!dir.exists(directory)) {
        dir.create(directory,
                   showWarnings = FALSE,
                   recursive = TRUE)
      }

      filename <- paste0(
        case_study,
        "_",
        method,
        "_metrics_vs_parameters_target_sample_size_per_arm_",
        target_sample_size_per_arm
      )

      if (other_params_str != ""){
        filename <- paste0(filename, "_", other_params_str)
      }

      title <- sprintf(
        "%s, $N_T/2 = $ %s, Source denominator change factor = %s",
        str_to_title(case_study),
        methods_labels[[method]]$label,
        other_params_label,
        target_sample_size_per_arm
      )
      title <- format_title(title = title, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)
      filename <- format_filename(filename = filename, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)

      file_path <- file.path(directory, filename)

      export_table(data_table = data_table, ncollapses = 1, title = title, file_path = file_path)

    }
  }
}

#' Plot metric vs sample size
#'
#' @param metric The metric to be plotted
#' @param results_metrics_df A dataframe containing the results metrics
#' @param case_study The case study being analyzed
#' @param method The method being used
#' @param parameters_combinations The combinations of parameters
#' @param theta_0 The null hypothesis value
#'
#' @return Generates and saves tables in HTML, PDF, and LaTeX formats
#' @export
#'
#' @importFrom dplyr filter mutate select arrange rename group_by ungroup
#' @importFrom rlang sym
#' @importFrom knitr kable
#' @importFrom kableExtra kable_styling row_spec collapse_rows
table_metric_vs_sample_size <- function(metric,
                                        results_metrics_df,
                                        case_study,
                                        method,
                                        theta_0 = 0,
                                        source_denominator_change_factor = 1,
                                        target_to_source_std_ratio = 1,
                                        wide_table = TRUE) {



  results_df <- results_metrics_df %>%
    dplyr::filter(
      case_study == !!case_study,
      method == !!method,
      source_denominator_change_factor == !!source_denominator_change_factor  |
        is.na(source_denominator_change_factor),
      target_to_source_std_ratio == !!target_to_source_std_ratio |
        is.na(target_to_source_std_ratio)
    )

  results_df <- cbind(results_df, get_parameters(results_df[, "parameters"]))

  case_study_config <- yaml::read_yaml(paste0(case_studies_config_dir, case_study, ".yml"))

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

  # Effects label
  effects <- c("No effect",
               "Partially consistent effect",
               "Consistent effect")

  # Factorize target treatment effect
  results_df$target_treatment_effect <- factor(
    results_df$target_treatment_effect,
    levels = unique(results_df$target_treatment_effect),
    labels = effects
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

  results_df$Parameters <- format_results_df_parameters(results_df, include_method_name = FALSE)

  for (treatment_effect in effects){
    # Filter on the treatment effect
    results_df_filtered <- results_df %>%
      dplyr::filter(target_treatment_effect %in% treatment_effect)

    data_table <- results_df_filtered %>%
      dplyr::mutate(
        means = format_num(!!rlang::sym(selected_metric$name)),
        error_low = format_num(!!rlang::sym(
          paste0(selected_metric$metric_uncertainty, "_lower")
        )),
        error_upper = format_num(!!rlang::sym(
          paste0(selected_metric$metric_uncertainty, "_upper")
        )),
        CI = paste0("[", error_low, ", ", error_upper, "]"),
        combined_mean_CI = paste0(means, " ", CI)  # Combine mean and CI into one column
      ) %>%
      dplyr::select(Parameters, target_sample_size_per_arm, combined_mean_CI) %>%
      dplyr::arrange(target_sample_size_per_arm)

    data_table <- data_table %>%
      dplyr::rename(
        "$N_T/2$" = target_sample_size_per_arm,
        !!selected_metric$label := combined_mean_CI  # Rename the combined column
      )

    if (wide_table && !(all(is.na(data_table$Parameters) | data_table$Parameters == ""))){
      # Pivot the table so that Parameters become columns and Sample Size per Arm is in rows
      data_table_formatted <- data_table %>%
        tidyr::pivot_wider(
          names_from = Parameters,  # The column to pivot (Parameters)
          values_from = !!rlang::sym(selected_metric$label)  # The values for the new columns
        )
      longtable = FALSE
    } else {
      data_table_formatted <- data_table
      longtable = TRUE
    }


    directory <- file.path(tables_dir, case_study)
    if (!dir.exists(directory)) {
      dir.create(directory, showWarnings = FALSE, recursive = TRUE)
    }
    filename <- paste0(
      case_study,
      "_",
      method,
      "_",
      frequentist_metrics[[metric]]$name,
      "_vs_target_sample_size_per_arm_",
      gsub(" ", "_", treatment_effect)
    )

    title <-  TeX(
      sprintf(
        "%s, %s, %s, %s",
        frequentist_metrics[[metric]]$label,
        str_to_title(case_study),
        methods_labels[[method]]$label,
        treatment_effect
      )
    )

    title <- format_title(title = title, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor, as_latex = TRUE)
    filename <- format_filename(filename = filename, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)


    file_path <- file.path(directory, filename)

    export_table(data_table = data_table_formatted, ncollapses = 1, title = title, file_path = file_path, longtable = longtable)
  }
}


table_metric_vs_drift_scenario_cat <- function(metric,
                                  results_metrics_df,
                                  theta_0,
                                  case_study,
                                  method,
                                  category,
                                  control_drift,
                                  target_sample_size_per_arm,
                                  parameters_combinations,
                                  xvars,
                                  wide_table = TRUE) {
  # Filter results by case study and target sample size per arm
  results_metrics_df <- results_metrics_df %>%
    dplyr::filter(case_study == !!case_study) %>%
    dplyr::filter(method == !!method)


  if (control_drift) {
    results_metrics_df <- results_metrics_df %>% dplyr::filter(treatment_drift == 0)
    xvar <- xvars$control_drift
  } else {
    results_metrics_df <- results_metrics_df %>% dplyr::filter(control_drift == 0)
    xvar <- xvars$drift
  }

  results_df <- results_metrics_df %>% dplyr::filter(method == !!method)

  #categories <- unique(results_df[["parameters"]])

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

    # Make sure that the data are sorted according to the parameters values
    params <- get_parameters(results_df[, "parameters"])
    column_names <- colnames(params)
    sorted_index <- do.call(order, lapply(params[column_names], as.factor))
    # Reorder the data :
    results_df[results_df$method == method, ] <- results_df[sorted_index, ]
  } else if (category == "source_denominator") {
    results_df <- results_df %>%
      dplyr::filter(
        target_sample_size_per_arm == !!target_sample_size_per_arm,
        target_to_source_std_ratio == !!target_to_source_std_ratio |
          is.na(target_to_source_std_ratio),
      )
  } else if (category == "target_to_source_std_ratio") {
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

  results_df$rows <- seq(1, nrow(results_df))

  parameters_labels <- NA
  if (category == "target_sample_size_per_arm") {
    results_df$label <- results_df$target_sample_size_per_arm
    labels_title <- "Target sample size per arm"
  } else if (category == "parameters") {
    parameters_labels <- unname(sapply(1:nrow(results_df), function(i) {
      process_method_parameters_label(results_df[i, ], methods_labels, method_name = FALSE)
    }))

    labels_title <- "Parameters"
  } else if (category == "source_denominator") {
    results_df$label <- results_df$source_denominator
    labels_title <- "Source denominator change factor"
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
    next
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

  labels <- list()
  for (i in seq_along(categories)) {
    combination <- categories[i]
    if (category == "source_denominator") {
      if (is.na(combination)) {
        next
      } else {
        results_df[results_df[, "source_denominator"] == combination, "label"] <- combination
      }
    }
  }


  if (!(category == "parameters")){
    parameters_label <- make_labels_from_parameters(parameters_combinations, method)
    parameters_str <- convert_params_to_str(methods_dict[[method]], parameters_combinations)

    if (parameters_label == "") {
      parameters_label_title <- ""
    } else {
      parameters_label_title <- paste0(", ", parameters_label)
    }

    results_df$labels <- parameters_labels
  }


  if (category == "target_sample_size_per_arm") {
    title <- sprintf(
      "%s, %s, %s %s",
      frequentist_metrics[[metric]]$label,
      str_to_title(case_study),
      methods_labels[[method]]$label,
      parameters_label_title
    )

    filename <- paste0(
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
  } else if (category == "source_denominator") {
    title <- sprintf(
      "%s, %s, %s %s, %s %s",
      frequentist_metrics[[metric]]$label,
      str_to_title(case_study),
      methods_labels[[method]]$label,
      parameters_label_title,
      "$N_T/2 = $",
      target_sample_size_per_arm
    )
    filename <- paste0(
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
      "%s, %s, %s, $N_T/2 = $ %s",
      frequentist_metrics[[metric]]$label,
      str_to_title(case_study),
      methods_labels[[method]]$label,
      target_sample_size_per_arm
    )
    filename <- paste0(
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
      "%s, %s, %s %s, %s %s",
      frequentist_metrics[[metric]]$label,
      str_to_title(case_study),
      methods_labels[[method]]$label,
      parameters_label_title,
      "$N_T/2 = $",
      target_sample_size_per_arm
    )
    filename <- paste0(
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

  title <- format_title(title = title, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor, as_latex = TRUE)
  filename <- format_filename(filename = filename, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)


  if (category == "target_sample_size_per_arm") {
    results_df$label <- results_df$target_sample_size_per_arm
    labels_title <- "Target sample size per arm"
  } else if (category == "parameters") {
    results_df$label <- format_results_df_parameters(results_df, include_method_name = FALSE)
    labels_title <- "Parameters"
  } else if (category == "source_denominator") {
    results_df$label <- results_df$source_denominator
    labels_title <- "Source denominator change factor"
  } else if (category == "target_to_source_std_ratio") {
    results_df$label <- results_df$target_to_source_std_ratio
    labels_title <- expression(Tex("$\\sigma_T/\\sigma_S")) # Ratio between target and source standard deviation
  } else {
    stop("Not implemented for this category")
  }

  data_table <- results_df %>%
    dplyr::select(
      label,
      !!xvar_name,
      !!selected_metric_name,
      !!selected_metric_uncertainty_lower,
      !!selected_metric_uncertainty_upper
    ) %>%
    dplyr::mutate(
      !!xvar_label := format_num(!!xvar_name),
      !!selected_metric$label := format_num(!!selected_metric_name),
      error_low = format_num(!!selected_metric_uncertainty_lower),
      error_upper = format_num(!!selected_metric_uncertainty_upper),
      !!selected_metric$uncertainty_label := paste0("[", error_low, ", ", error_upper, "]")
    ) %>%
    dplyr::select(
      label,
      !!xvar_label,
      !!selected_metric$label,
      !!selected_metric$uncertainty_label
    ) %>%
    dplyr::rename(
      "$\\delta$" = "Drift in treatment effect",
    )


  # Merge Mean and CI columns
  data_table <- data_table %>%
    unite(!!selected_metric$label, selected_metric$label, selected_metric$uncertainty_label, sep = " ")
  if (wide_table && !(all(is.na(data_table$label) | data_table$label == ""))){
    # Pivot the table so that Parameters become columns and Sample Size per Arm is in rows
    data_table_formatted <- data_table %>%
      tidyr::pivot_wider(
        names_from = label,  # The column to pivot (label)
        values_from = !!rlang::sym(selected_metric$label)  # The values for the new columns
      )
    longtable = FALSE
  } else {
    data_table_formatted <- data_table
    longtable = TRUE
  }


  directory <- file.path(tables_dir, case_study)
  if (!dir.exists(directory)) {
    dir.create(directory, showWarnings = FALSE, recursive = TRUE)
  }

  file_path <- file.path(directory, filename)


  export_table(data_table = data_table_formatted, ncollapses = 1, title = title, file_path = file_path)

}



#' Generate tables for methods operating characteristics
#'
#' @param results_metrics_df A dataframe containing the results metrics
#' @param metrics A list of metrics to be analyzed
#'
#' @return Generates and saves multiple tables for different metrics and methods
#' @export
#'
#' @importFrom dplyr filter
tables_metric_vs_scenario <- function(results_metrics_df, metrics) {
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
      parameters_combinations <- unique(get_parameters(filtered_results_metrics_df[filtered_results_metrics_df$method == method, "parameters"]))


      for (metric in names(metrics)) {
        for (sample_size in target_sample_sizes) {

          # table_metric_vs_drift(
          #   metric = metric,
          #   results_metrics_df = filtered_results_metrics_df,
          #   case_study = case_study,
          #   method = method,
          #   theta_0 = theta_0,
          #   category = "target_to_source_std_ratio",
          #   control_drift = FALSE,
          #   target_sample_size_per_arm = sample_size,
          #   parameters_combinations = NULL,
          #   xvars = xvars
          # )
          # table_metric_vs_drift_scenario_cat(
          #   metric = metric,
          #   results_metrics_df = filtered_results_metrics_df,
          #   case_study = case_study,
          #   method = method,
          #   theta_0 = theta_0,
          #   control_drift = FALSE,
          #   target_sample_size_per_arm = sample_size,
          #   parameters_combinations = NULL,
          #   xvars = xvars,
          #   category = "target_to_sourc_std_ratio"
          # )
          #
          # table_metric_vs_drift_scenario_cat(
          #   metric = metric,
          #   results_metrics_df = filtered_results_metrics_df,
          #   case_study = case_study,
          #   method = method,
          #   theta_0 = theta_0,
          #   control_drift = FALSE,
          #   target_sample_size_per_arm = sample_size,
          #   parameters_combinations = NULL,
          #   xvars = xvars,
          #   category = "source_denominator_change_factor"
          # )
        }

        for (source_denominator_change_factor in unique(filtered_results_metrics_df$source_denominator_change_factor)) {
          for (target_to_source_std_ratio in unique(filtered_results_metrics_df$target_to_source_std_ratio)) {
            for (sample_size in target_sample_sizes) {
              for (i in 1:nrow(parameters_combinations)) {
                table_metric_vs_drift(
                  metric = metric,
                  results_metrics_df = filtered_results_metrics_df,
                  case_study = case_study,
                  method = method,
                  theta_0 = theta_0,
                  category = "parameters",
                  control_drift = FALSE,
                  target_sample_size_per_arm = sample_size,
                  parameters_combinations = NULL,
                  xvars = xvars
                )
              }
            }

            table_metric_vs_sample_size(
              metric = metric,
              results_metrics_df = filtered_results_metrics_df,
              case_study = case_study,
              method = method,
              theta_0 = theta_0,
              wide_table = TRUE
            )
          }
        }
      }
    }
  }
}
