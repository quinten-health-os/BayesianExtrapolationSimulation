#' Compute the Hellinger distance between two normal distributions.
#'
#' @param mu1 Mean of the first normal distribution.
#' @param sigma1 Standard deviation of the first normal distribution.
#' @param mu2 Mean of the second normal distribution.
#' @param sigma2 Standard deviation of the second normal distribution.
#'
#' @return The Hellinger distance between the two normal distributions.
#'
hellinger_distance <- function(mu1, sigma1, mu2, sigma2) {
  sqrt(1 - sqrt(2 * sigma1 * sigma2 / (sigma1 ^ 2 + sigma2 ^ 2)) * exp(-0.25 * (mu1 - mu2) ^
                                                                         2 / (sigma1 ^ 2 + sigma2 ^ 2)))
}

#' Compute the drift range for a given simulation and case study configuration.
#'
#' @param scenarios_config Configuration for the simulation.
#' @param case_study_config Configuration for the case study.
#'
#' @return A vector representing the drift range.
compute_drift_range <- function(scenarios_config, case_study_config) {
  # theta_0 <- case_study_config$theta_0
  source_treatment_effect <- case_study_config$source$treatment_effect

  if (case_study_config$summary_measure_likelihood == "normal") {
    sigma1 <- case_study_config$source$standard_error
    sigma2 <- case_study_config$target$standard_error


    max_drift <- abs(source_treatment_effect) * 10
    target_hellinger <- function(mu_diff) {
      hellinger_distance(0, sigma1, mu_diff, sigma2) - 0.9
    }
    solver_result <- stats::uniroot(target_hellinger, c(0, max_drift))
    upper_bound <- solver_result$root
    lower_bound <- -upper_bound

    # if (lower_bound > theta_0 - source_treatment_effect) {
    #   lower_bound <- theta_0 - source_treatment_effect
    # }
    # upper_bound <- abs(lower_bound) / 2
  } else if (case_study_config$summary_measure_likelihood == "binomial") {
    source_control_rate <- case_study_config$source$responses$control / case_study_config$source$control
    source_treatment_rate <- case_study_config$source$responses$treatment / case_study_config$source$treatment
    source_treatment_effect <- source_treatment_rate - source_control_rate
    lower_bound <- max(-1 - source_treatment_effect, -source_treatment_rate)
    upper_bound <- min(1 - source_treatment_effect, 1 - source_treatment_rate)
  } else {
    stop("Not implemented for other treatment effect distributions")
  }
  drift_range <- seq(lower_bound, upper_bound, length.out = scenarios_config$ndrift)
  return(drift_range)
}

#' Get important drift values for a given source treatment effect and case study configuration.
#'
#' @param source_treatment_effect The source treatment effect.
#' @param case_study_config Configuration for the case study.
#'
#' @return A vector of important drift values to consider.
important_drift_values <- function(source_treatment_effect,
                                   case_study_config) {
  vec <- c(-source_treatment_effect, -source_treatment_effect / 2, 0)
  if (case_study_config$summary_measure_likelihood == "normal") {
    return(vec)
  } else if (case_study_config$summary_measure_likelihood == "binomial") {
    lower_bound <- -1 - source_treatment_effect
    upper_bound <- 1 - source_treatment_effect
    vec <- vec[vec <= upper_bound & vec >= lower_bound]
    return(vec)
  }
}

#' Compute the control drift range for a given source treatment effect.
#'
#' @param source_treatment_effect The source treatment effect.
#'
#' @return A vector representing the control drift range.
compute_control_drift_range <- function(source_treatment_effect) {
  return(c(0))
}

#' Compute the source denominator range for a given source data, simulation configuration, and case study configuration.
#'
#' @param source_data The source data.
#' @param scenarios_config Configuration for the simulation.
#' @param case_study_config Configuration for the case study.
#'
#' @return A vector representing the source denominator range.
compute_source_denominator_range <- function(source_data,
                                             scenarios_config,
                                             case_study_config) {
  # Compute, for each value of the source denominator, the corresponding numerator so that the treatment effec stays the same
  # Default behavior, used when the endpoint is continuous
  source_denominator_range <- c(NA)
  source_treatment_rate_range <- NA

  denominator_change_factor <- unlist(scenarios_config$denominator_change_factor)
  if (source_data$endpoint == "binary" &&
      case_study_config$summary_measure_likelihood == "normal") {
    odds_control <- source_data$control_rate / (1 - source_data$control_rate)
    source_denominator_range <- denominator_change_factor * odds_control

    source_treatment_rate_range <- 1 / (1 + exp(-(
      source_data$treatment_effect_estimate + log(source_denominator_range)
    )))
  } else if (source_data$endpoint == "time_to_event" ||
             source_data$endpoint == "recurrent_event") {
    source_denominator_range <- source_data$control_rate * denominator_change_factor

    source_treatment_rate_range <- source_denominator_range * exp(source_data$treatment_effect_estimate)
  } else {
    denominator_change_factor <- c(1)
  }

  return(
    list(
      source_denominator = source_denominator_range,
      source_denominator_change_factor = denominator_change_factor
    )
  )
}

#' Combine parameters for a given method.
#'
#' @param method The method name.
#'
#' @return A list of combinations of parameters for the method.
combine_parameters <- function(method) {
  if (is.null(methods_dict[[method]])) {
    return(NULL)
  } else {
    param_info <- methods_dict[[method]]
    ranges <- lapply(param_info, function(x) {
      x$range
    })
    combinations <- expand.grid(ranges)
    return(combinations)
  }
}

#' Generate simulation scenarios based on the configuration directory.
#'
#' @param config_dir The directory containing the configuration files.
#' @param scenarios_config The simulation configuration file.
#'
#' @return A list of simulation scenarios.
simulation_scenarios <- function(config_dir, scenarios_config) {
  source(paste0(config_dir, "methods_config.R"))
  case_studies <- scenarios_config$case_studies

  methods_configurations_df <- data.frame()

  # We store the parameters in a dataframe with two columns : method and parameters
  for (method in names(methods_dict)) {
    if (method %in% scenarios_config$methods) {
      parameters_combinations <- combine_parameters(method)

      converted_parameters_df <- data.frame(method = method, parameters = I(apply(parameters_combinations, 1, function(row) {
        params_list <- as.list(row)
        names(params_list) <- colnames(parameters_combinations)
        return(params_list)
      })))

      methods_configurations_df <- rbind(methods_configurations_df, converted_parameters_df)
    }
  }

  if (any(!(scenarios_config$methods %in% names(methods_dict)))) {
    stop("Some methods listed in the simulation configuration need configuration.")
  }

  scenarios_configurations <- list()

  for (case_study in case_studies) {
    case_study_config <- yaml::read_yaml(paste0(case_studies_config_dir, case_study, ".yml"))
    source_data <- ObservedSourceData$new(case_study_config)

    sampling_approximation <- case_study_config$sampling_approximation
    theta_0 <- case_study_config$theta_0
    null_space <- case_study_config$null_space

    source_denominator <- compute_source_denominator_range(source_data, scenarios_config, case_study_config)
    source_denominator <- data.frame(source_denominator)

    drift_range <- compute_drift_range(scenarios_config, case_study_config)

    source_treatment_effect_estimate <- source_data$treatment_effect_estimate

    case_study_config <- yaml::read_yaml(paste0(case_studies_config_dir, case_study, ".yml"))

    mandatory_drift_values <- important_drift_values(source_treatment_effect_estimate, case_study_config)
    drift_range <- sort(c(drift_range, mandatory_drift_values))

    control_drift_range <- compute_control_drift_range(source_treatment_effect_estimate)

    if (case_study_config$endpoint == "continuous") {
      target_to_source_std_ratio_range <- scenarios_config$target_to_source_std_ratio_range
    } else {
      target_to_source_std_ratio_range <- NA
    }

    n <- nrow(source_denominator)

    # Create drift/control drift pairs
    drift_combinations <- expand.grid(
      drift = drift_range,
      control_drift = control_drift_range,
      target_to_source_std_ratio = target_to_source_std_ratio_range
    )

    source_denominator_repeated <- source_denominator[rep(seq_len(n), each = nrow(drift_combinations)), ]
    drift_combinations_repeated <- drift_combinations[rep(seq_len(nrow(drift_combinations)), times = n), ]

    drift_combinations <- cbind(source_denominator_repeated,
                                drift_combinations_repeated)

    drift_combinations <- unique(drift_combinations)

    # Exclude combinations that would result in rates that are outside [0,1]
    if (case_study_config$summary_measure_likelihood == "binomial") {
      lower_bound <- pmax(
        -source_data$treatment_rate - drift_combinations$control_drift,
        -1 - source_data$treatment_effect_estimate
      )
      upper_bound <- pmin(
        1 - source_data$treatment_rate - drift_combinations$control_drift,
        1 - source_data$treatment_effect_estimate
      )

      condition <- drift_combinations$drift >= lower_bound &
        drift_combinations$drift <= upper_bound
      drift_combinations <- subset(drift_combinations, condition)
    }

    total_target_sample_sizes <- (source_data$sample_size_control + source_data$sample_size_treatment) / scenarios_config$sample_size_factors
    current_scenario_configurations <- list()

    for (i in 1:nrow(drift_combinations)) {
      drift <- drift_combinations[i, "drift"]
      control_drift <- drift_combinations[i, "control_drift"]
      source_denominator <- drift_combinations[i, "source_denominator"]
      source_denominator_change_factor <- drift_combinations[i, "source_denominator_change_factor"]
      treatment_drift <- drift + control_drift
      target_treatment_effect <- source_treatment_effect_estimate + drift
      for (total_target_sample_size in total_target_sample_sizes) {
        target_sample_size_per_arm <- floor(total_target_sample_size / 2)

        current_scenario_configurations[[length(current_scenario_configurations) + 1]] <- list(
          target_sample_size_per_arm = target_sample_size_per_arm,
          drift = drift,
          treatment_drift = treatment_drift,
          control_drift = control_drift,
          source_denominator = source_denominator,
          source_denominator_change_factor = source_denominator_change_factor,
          case_study = case_study,
          sampling_approximation = sampling_approximation,
          source_treatment_effect_estimate = source_treatment_effect_estimate,
          target_treatment_effect = target_treatment_effect,
          target_to_source_std_ratio = drift_combinations[i, "target_to_source_std_ratio"],
          theta_0 = theta_0,
          null_space = null_space
        )
      }
    }
    scenarios_configurations <- c(scenarios_configurations,
                                  current_scenario_configurations)
  }
  # Cartesian product between methods configurations and scenarios
  scenarios_configurations_df <- data.frame(do.call(rbind, scenarios_configurations))

  cases_df <- merge(methods_configurations_df,
                    scenarios_configurations_df,
                    by = NULL)

  return(cases_df)
}

scenarios_table_ranges <- function(cases_df, results_path) {
  dir.create(results_path, showWarnings = FALSE, recursive = TRUE)
  file_path <- paste0(results_path, "/drift_ranges.tex")

  # Initialize an empty dataframe
  df <- data.frame()


  for (case_study in unique(cases_df$case_study)) {
    subdf <- cases_df[cases_df$case_study == case_study, ]
    subdf$theta_0 <- as.numeric(subdf$theta_0)
    subdf$source_treatment_effect_estimate <- as.numeric(subdf$source_treatment_effect_estimate)

    # Compute the ranges
    drift_range <- range(subdf$drift)
    treatment_range <- range(subdf$target_treatment_effect)

    # Bind rows to df
    df <- rbind(
      df,
      data.frame(
        Drug = paste0(toupper(substring(case_study, 1, 1)), substring(case_study, 2)),
        Drift_range_lower = signif(drift_range[1], 3),
        Drift_range_upper = signif(drift_range[2], 3),
        Drift_with_theta_T_equals_theta_0 = signif(
          unique(subdf$theta_0 - subdf$source_treatment_effect_estimate),
          3
        ),
        Treatment_effect_range_lower = signif(treatment_range[1], 3),
        Treatment_effect_range_upper = signif(treatment_range[2], 3)
      )
    )
  }

  # Function to generate the LaTeX table
  generate_latex_table <- function(df, file_path) {
    sink(file_path)

    cat(
      "\\begin{table}[h]
\\centering
\\begin{tabular}{|l|c|c|c|}
\\textbf{Case study} & \\textbf{Drift range} & \\textbf{Drift with $\\theta_T^{(true)} = \\theta_0$} & \\textbf{Treatment effect range} \\\\
\\midrule
"
    )

    for (i in 1:nrow(df)) {
      # Concatenate ranges
      drift_range <- paste0("[",
                            df$Drift_range_lower[i],
                            ",",
                            df$Drift_range_upper[i],
                            "]")
      treatment_range <- paste0(
        "[",
        df$Treatment_effect_range_lower[i],
        ",",
        df$Treatment_effect_range_upper[i],
        "]"
      )

      cat(
        paste0(
          "\\hline \\text{",
          df$Drug[i],
          "} & ",
          drift_range,
          " & ",
          df$Drift_with_theta_T_equals_theta_0[i],
          " & ",
          treatment_range,
          " \\\\
"
        )
      )
    }

    cat(
      "\\hline
    \\end{tabular}
\\caption{Drift ranges considered for each case study.}
\\label{tab:drift_ranges}
\\end{table}"
    )

    # Close the file connection
    sink()
  }

  # Generate the LaTeX table
  generate_latex_table(df, file_path)
}
