check_probability_value <- function(p) {
  if (length(p) == 0) {
    return(NA)
  }

  if (is.na(p)) {
    return(NA)
  }

  if (p < 0) {
    warnings("Invalid value for the probability, should be between 0 and 1. Rounding the value to 0.")
    futile.logger::flog.warn("Invalid value for the probability, should be between 0 and 1. Rounding the value to 0.")
    p <- 0
  }

  if (p > 1) {
    warnings("Invalid value for the probability, should be between 0 and 1. Rounding the value to 1.")
    futile.logger::flog.warn("Invalid value for the probability, should be between 0 and 1. Rounding the value to 0.")
    p <- 1
  }
  return(p)
}

#' Calculate the probability of false positive
#'
#' @description This function calculates the probability of false positive based on the conditional probability of success,
#' design prior samples, theta_0, and null space.
#'
#' @param conditional_proba_success The conditional probability of success.
#' @param design_prior_samples The design prior samples.
#' @param theta_0 The value of theta_0.
#' @param null_space The null space (either "left" or "right").
#'
#' @return The probability of false positive.
#'
#' @export
#'
#' @examples NA
preposterior_proba_FP_MC <- function(conditional_proba_success,
                                     design_prior_samples,
                                     theta_0,
                                     null_space) {
  if (null_space == "left") {
    samples_in_null_space <- design_prior_samples <= theta_0
  } else if (null_space == "right") {
    samples_in_null_space <- design_prior_samples >= theta_0
  }

  p <- mean(conditional_proba_success * samples_in_null_space)

  p <- sapply(p, check_probability_value)

  assertions::assert_number(p)
  return(p)
}

#' Calculate the probability of true positive
#'
#' @description This function calculates the probability of true positive based on the conditional probability of success,
#' design prior samples, theta_0, and null space.
#'
#' @param conditional_proba_success The conditional probability of success.
#' @param design_prior_samples The design prior samples.
#' @param theta_0 The value of theta_0.
#' @param null_space The null space (either "left" or "right").
#'
#' @return The probability of true positive.
#'
#' @export
#'
#' @examples NA
preposterior_proba_TP_MC <- function(conditional_proba_success,
                                     design_prior_samples,
                                     theta_0,
                                     null_space) {
  if (null_space == "left") {
    samples_in_alt_space <- design_prior_samples > theta_0
  } else if (null_space == "right") {
    samples_in_alt_space <- design_prior_samples < theta_0
  }

  p <- mean(conditional_proba_success * samples_in_alt_space)

  p <- sapply(p, check_probability_value)

  assertions::assert_number(p)
  return(p)
}

#' Calculate the average type 1 error
#'
#' @description This function calculates the average type 1 error based on the preposterior probability of false positive
#' and the prior probability of no benefit.
#'
#' @param prepost_proba_FP The preposterior probability of false positive.
#' @param prior_proba_no_benefit The prior probability of no benefit.
#'
#' @return The average type 1 error.
#'
#' @export
#'
#' @examples NA
average_tie <- function(prepost_proba_FP,
                        prior_proba_no_benefit) {
  if (is.na(prepost_proba_FP) | is.na(prior_proba_no_benefit)) {
    return(NA)
  }

  assertions::assert_number(prepost_proba_FP)
  assertions::assert_number(prior_proba_no_benefit)

  p <- prepost_proba_FP / prior_proba_no_benefit
  p <- sapply(p, check_probability_value)

  assertions::assert_number(p)
  return(p)
}

#' Calculate the average power
#'
#' @description This function calculates the average power based on the preposterior probability of true positive
#' and the prior probability of no benefit.
#'
#' @param prepost_proba_TP The preposterior probability of true positive.
#' @param prior_proba_no_benefit The prior probability of no benefit.
#'
#' @return The average power.
#'
#' @export
#'
#' @examples NA
average_power <- function(prepost_proba_TP,
                          prior_proba_no_benefit) {
  if (is.na(prepost_proba_TP) | is.na(prior_proba_no_benefit)) {
    return(NA)
  }

  assertions::assert_number(prepost_proba_TP)
  assertions::assert_number(prior_proba_no_benefit)

  p <- prepost_proba_TP / (1 - prior_proba_no_benefit)
  p <- sapply(p, check_probability_value)

  assertions::assert_number(p)
  return(p)
}

#' Calculate the upper bound probability of false positive
#'
#' @description This function calculates the upper bound probability of false positive based on the model, prior probability of no benefit,
#' source data, theta_0, target sample size per arm, case study configuration, number of replicates, confidence level,
#' null space, and critical value.
#'
#' @param model The model.
#' @param prior_proba_no_benefit The prior probability of no benefit.
#' @param source_data The source data.
#' @param theta_0 The value of theta_0.
#' @param target_sample_size_per_arm The target sample size per arm.
#' @param case_study_config The case study configuration.
#' @param n_replicates The number of replicates.
#' @param confidence_level The confidence level.
#' @param null_space The null space (either "left" or "right").
#' @param critical_value The critical value.
#'
#' @return The upper bound probability of false positive, defined as \eqn{Pr(Study success|\theta_T = \theta_0) \times Pr(\theta_T \leq \theta_0)}
#'
#' @examples NA
#'
#' @export
upper_bound_proba_FP_MC <- function(model,
                                    prior_proba_no_benefit,
                                    source_data,
                                    theta_0,
                                    target_sample_size_per_arm,
                                    case_study_config,
                                    n_replicates,
                                    confidence_level,
                                    null_space,
                                    critical_value) {
  # We need to estimate Pr(Study success|\theta_T = \theta_0)
  treatment_drift <- theta_0 - source_data$treatment_effect_estimate

  target_data <- TargetDataFactory$new()
  target_data <- target_data$create(
    source_data = source_data,
    case_study_config = case_study_config,
    target_sample_size_per_arm = target_sample_size_per_arm,
    treatment_drift = treatment_drift,
    control_drift = 0,
    summary_measure_likelihood = source_data$summary_measure_likelihood
  )

  results <- model$simulation_for_given_treatment_effect(
    target_data = target_data,
    n_replicates = n_replicates,
    critical_value = critical_value,
    theta_0 = theta_0,
    confidence_level = confidence_level,
    null_space = null_space,
    to_return = c("test_decision")
  )

  p <- prior_proba_no_benefit * mean(results$test_decisions)

  p <- sapply(p, check_probability_value)
  assertions::assert_number(p)
  return(p)
}


# The goal here is to perform numerical integration based on a predefined range of values for the treatment effect. This range of values should be the same for all design priors (but not all methods) methods. Since we integrate the conditional power over different design priors and different bounds, the range should


upper_bound_proba_FP <- function(prior_proba_no_benefit,
                                 treatment_effect_values,
                                 conditional_proba_success, theta_0) {
  if (is.na(prior_proba_no_benefit)) {
    return(NA)
  }

  assertions::assert_number(prior_proba_no_benefit)

  p <- prior_proba_no_benefit * conditional_proba_success[treatment_effect_values == theta_0]

  p <- sapply(p, check_probability_value)

  if (is.list(p)) {
    return(NA)
  }

  if (length(p)>1){
    stop("Length of p is >1.")
  }

  if (is.na(p)){
    return(p)
    #stop("p is NA.")
  }

  assertions::assert_number(p)
  return(p)
}

# Preposterior probability of a False Positive
preposterior_proba_FP <- function(conditional_proba_success,
                                  treatment_effect_values,
                                  theta_0,
                                  null_space,
                                  design_prior_pdf) {
  if (null_space == "left") {
    values_in_null_space <- treatment_effect_values <= theta_0
  } else if (null_space == "right") {
    values_in_null_space <- treatment_effect_values >= theta_0
  }

  design_prior_pdf <- design_prior_pdf[values_in_null_space]

  if (sum(values_in_null_space) < 2) {
    p <- NA
    warning(
      "The number of values in the null space is less than 2, the preposterior probability of FP cannot be computed."
    )
  } else {
    p <- Bolstad2::sintegral(
      x = treatment_effect_values[values_in_null_space],
      fx = design_prior_pdf * conditional_proba_success[values_in_null_space],
      n.pts = sum(values_in_null_space)
    )$int

    p <- sapply(p, check_probability_value)
    assertions::assert_number(p)
  }
  return(p)
}

preposterior_proba_TP <- function(conditional_proba_success,
                                  treatment_effect_values,
                                  theta_0,
                                  null_space,
                                  design_prior_pdf) {
  if (null_space == "left") {
    values_in_alt_space <- treatment_effect_values > theta_0
  } else if (null_space == "right") {
    values_in_alt_space <- treatment_effect_values < theta_0
  }

  if (sum(values_in_alt_space) < 2) {
    p <- NA
  } else {
    p <- Bolstad2::sintegral(
      x = treatment_effect_values[values_in_alt_space],
      fx = design_prior_pdf[values_in_alt_space] * conditional_proba_success[values_in_alt_space],
      n.pts = sum(values_in_alt_space)
    )$int

    p <- sapply(p, check_probability_value)
    assertions::assert_number(p)
  }
  return(p)
}


prior_proba_success <- function(conditional_proba_success,
                                treatment_effect_values,
                                design_prior_pdf) {
  if (length(conditional_proba_success) < 2) {
    return(NA)
  }

  p <- Bolstad2::sintegral(
    x = treatment_effect_values,
    fx = design_prior_pdf * conditional_proba_success,
    n.pts = length(treatment_effect_values)
  )$int

  p <- sapply(p, check_probability_value)
  assertions::assert_number(p)
  return(p)
}

compute_bayesian_ocs <- function(results_freq_df, env) {
  # Initialize an empty dataframe to store the results
  results_bayesian_ocs <- data.frame(
    case_study = character(),
    method = character(),
    target_sample_size_per_arm = numeric(),
    parameters_combination = character(),
    design_prior_type = character(),
    prior_proba_success = numeric(),
    prior_proba_no_benefit = numeric(),
    prior_proba_benefit = numeric(),
    prepost_proba_FP_SI = numeric(),
    average_tie_SI = numeric(),
    prepost_proba_TP_SI = numeric(),
    average_power_SI = numeric(),
    upper_bound_proba_FP_SI = numeric()
  )

  config_dir <- paste0(system.file(paste0("conf/", env), package = "RBExT"), "/")
  scenarios_config <- yaml::yaml.load_file(paste0(config_dir, "scenarios_config.yml"))
  simulaton_config <- yaml::yaml.load_file(system.file("conf/simulation_config.yml", package = "RBExT"))

  if (is.null(simulation_config)) {
    stop("Simulation config is NULL")
  }

  design_prior_types <- cbind("ui_design_prior", "analysis_prior", "source_posterior")

  # Get the list of case studies
  case_studies <- unique(results_freq_df$case_study)

  mcmc_config <- yaml::read_yaml(paste0(config_dir, "/mcmc_config.yml"))

  for (case_study in case_studies) {
    case_study_config <- yaml::yaml.load_file(system.file(
      paste0("conf/case_studies/", case_study, ".yml"),
      package = "RBExT"
    ))
    results_df_0 <- results_freq_df %>%
      dplyr::filter(case_study == !!case_study)
    null_space <- case_study_config$null_space
    theta_0 <- case_study_config$theta_0


    source_denominator_change_factors <- unique(results_df_0$source_denominator_change_factor)

    for (source_denominator_change_factor in source_denominator_change_factors) {
      results_df_1 <- results_df_0 %>%
        dplyr::filter(
          source_denominator_change_factor == !!source_denominator_change_factor  |
            is.na(source_denominator_change_factor)
        )

      source_data <- load_data(results_df_1[1, ],
                               type = "source",
                               reload_data_objects = TRUE)

      source_data_df <- data.frame(source_data$to_dict())
      to_remove <- c("source_control_rate", "source_treatment_rate")

      source_data_df <- remove_columns_from_df(source_data_df, to_remove)


      target_to_source_std_ratio_range <- unique(results_df_1$target_to_source_std_ratio)
      for (target_to_source_std_ratio in target_to_source_std_ratio_range) {
        results_df_2 <- results_df_1 %>%
          dplyr::filter(target_to_source_std_ratio == !!target_to_source_std_ratio |
                          is.na(target_to_source_std_ratio))


        if (nrow(results_df_2) == 0) {
          warning("Dataframe is empty")
          next
        }

        theta_0 <- case_study_config$theta_0

        # Get the list of methods
        methods <- unique(results_df_2$method)

        for (method in methods) {
          results_df_3 <- results_df_2 %>%
            dplyr::filter(method == !!method)

          target_sample_sizes <- unique(results_df_3$target_sample_size_per_arm)

          for (target_sample_size_per_arm in target_sample_sizes) {
            results_df_4 <- results_df_3 %>%
              dplyr::filter(target_sample_size_per_arm == !!target_sample_size_per_arm)

            # Get the different parameters combinations studies for this method
            parameters_combinations <- data.frame(parameters = unique(results_df_4[, "parameters"]))

            for (i in 1:nrow(parameters_combinations)) {
              # We unpack the parameter inside this loop (and not inside the previous one), because for some methods such as the commensurate power prior, there is a nested parameters structure which implies that they cannot all be stored in a single dataframe.
              method_parameters <- as.list(get_parameters(parameters_combinations[i, , drop = FALSE]))
              # convert strings to numeric or boolean if possible
              method_parameters <- data.frame(lapply(method_parameters, convert_if_possible))

              # The following applies to a single row dataframe, to recover the nested list structure
              method_params <- extract_nested_parameter(method_parameters)

              results_df <- results_df_4 %>%
                dplyr::filter(parameters == !!unlist(parameters_combinations[i, ]))

              conditional_proba_success <- results_df$success_proba
              treatment_effect_values <- results_df$target_treatment_effect

              model <- Model$new()
              model <- model$create(
                case_study_config = case_study_config,
                method = method,
                method_parameters = method_params,
                source_data = source_data,
                mcmc_config = mcmc_config
              )

              for (design_prior_type in design_prior_types) {
                if (design_prior_type == "analysis_prior" &&
                    model$empirical_bayes) {
                  # We cannot compute the Bayesian OCs with an analysis design prior for empirical Bayes methods.
                  results_to_add <- data.frame(
                    case_study = case_study,
                    method = method,
                    target_sample_size_per_arm = target_sample_size_per_arm,
                    parameters = parameters_combinations[i, ],
                    source_denominator_change_factor = source_denominator_change_factor,
                    target_to_source_std_ratio = target_to_source_std_ratio,
                    design_prior_type = design_prior_type,
                    prior_proba_success = NA,
                    prior_proba_no_benefit = NA,
                    prior_proba_benefit = NA,
                    prepost_proba_FP = NA,
                    prepost_proba_TP = NA,
                    average_tie = NA,
                    average_power = NA,
                    upper_bound_proba_FP = NA
                  )
                } else {
                  design_prior <- DesignPrior$new()
                  design_prior <- design_prior$create(
                    design_prior_type = design_prior_type,
                    model = model,
                    source_data = source_data,
                    case_study_config = case_study_config,
                    simulation_config = simulation_config,
                    mcmc_config = mcmc_config,
                    case_study = case_study
                  )

                  design_prior_pdf <- design_prior$pdf(treatment_effect_values)

                  if (null_space == "left") {
                    prior_proba_no_benefit <- design_prior$cdf(theta_0)
                  } else if (null_space == "right") {
                    prior_proba_no_benefit <- 1 - design_prior$cdf(theta_0)
                  }

                  prior_proba_benefit <- 1 - prior_proba_no_benefit

                  prior_proba_success_SI <- prior_proba_success(
                    conditional_proba_success = conditional_proba_success,
                    treatment_effect_values = treatment_effect_values,
                    design_prior_pdf = design_prior_pdf
                  )

                  prepost_proba_FP_SI <- preposterior_proba_FP(
                    conditional_proba_success = conditional_proba_success,
                    treatment_effect_values = treatment_effect_values,
                    theta_0 = theta_0,
                    null_space = null_space,
                    design_prior_pdf = design_prior_pdf
                  )

                  prepost_proba_TP_SI <- preposterior_proba_TP(
                    conditional_proba_success = conditional_proba_success,
                    treatment_effect_values = treatment_effect_values,
                    theta_0 = theta_0,
                    null_space = null_space,
                    design_prior_pdf = design_prior_pdf
                  )

                  average_tie_SI <- average_tie(
                    prepost_proba_FP = prepost_proba_FP_SI,
                    prior_proba_no_benefit = prior_proba_no_benefit
                  )

                  average_power_SI <- average_power(
                    prepost_proba_TP = prepost_proba_TP_SI,
                    prior_proba_no_benefit = prior_proba_no_benefit
                  )

                  if(length(conditional_proba_success[treatment_effect_values == theta_0]) > 1){
                    stop("Too many values corresponding to TIE.")
                  }
                  upper_bound_proba_FP_SI <- upper_bound_proba_FP(
                    prior_proba_no_benefit = prior_proba_no_benefit,
                    treatment_effect_values = treatment_effect_values,
                    conditional_proba_success = conditional_proba_success,
                    theta_0 = theta_0
                  )

                  results_to_add <- data.frame(
                    case_study = case_study,
                    method = method,
                    target_sample_size_per_arm = target_sample_size_per_arm,
                    parameters = parameters_combinations[i, ],
                    source_denominator_change_factor = source_denominator_change_factor,
                    target_to_source_std_ratio = target_to_source_std_ratio,
                    design_prior_type = design_prior_type,
                    prior_proba_success = prior_proba_success_SI,
                    prior_proba_no_benefit = prior_proba_no_benefit,
                    prior_proba_benefit = prior_proba_benefit,
                    prepost_proba_FP = prepost_proba_FP_SI,
                    prepost_proba_TP = prepost_proba_TP_SI,
                    average_tie = average_tie_SI,
                    average_power = average_power_SI,
                    upper_bound_proba_FP = upper_bound_proba_FP_SI
                  )
                }
                results_to_add <- cbind(results_to_add, source_data_df)
                results_bayesian_ocs <- rbind(results_bayesian_ocs, results_to_add)
              }
            }
          }
        }
      }
    }
  }
  return(results_bayesian_ocs)
}
