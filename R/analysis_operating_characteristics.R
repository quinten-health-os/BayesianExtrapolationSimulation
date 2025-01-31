#' Compute the frequentist power
#'
#' @description This function computes the power of a test for a given significance level
#' @description This function computes the power of a test for a given significance level
#'
#' @param alpha The significance level.
#' @param target_data The target data containing sample size, treatment effect, standard deviation, and summary measure distribution.
#' @param frequentist_test The type of frequentist test ("t-test" or "z-test").
#' @param theta_0 The null hypothesis value.
#' @param target_data Target data object
#' @param frequentist_test Type of frequentist test to apply, either z-test or t-test
#' @param theta_0 Boundary of the null hypothesis space
#' @param null_space Side of the null space, either left or right.
#'
#' @return The power of the test.
#'
#' @export
compute_freq_power <- function(alpha,
                               target_data,
                               frequentist_test,
                               theta_0,
                               null_space,
                               simulation_config,
                               case_study = NULL,
                               n_replicates = 1000) {
  if (null_space == "left") {
    alternative <- "greater"
  } else if (null_space == "right") {
    alternative <- "less"
  } else {
    stop("Null space must be either 'left' or 'right'")
  }

  if (is.na(alpha)){
    return(NA)
  }

  assertions::assert_number(target_data$treatment_effect)
  assertions::assert_number(target_data$standard_deviation)
  assertions::assert_number(target_data$sample_size_per_arm)

  power <- NA # Default value in case of an unsupported distribution

  case_study_config <- yaml::read_yaml(paste0(case_studies_config_dir, case_study, ".yml"))

  if (target_data$summary_measure_likelihood == "normal") {
    if (target_data$endpoint == "normal"  || target_data$endpoint == "continuous"  || case_study == "mepolizumab"){
      # In this case, we use an analytical computation of power
      effect_size <- (target_data$treatment_effect - theta_0) / target_data$standard_deviation

      if (frequentist_test == "t-test") {
        # Use pwr::pwr.t.test for a t-test power calculation
        power <- pwr::pwr.t.test(
          d = effect_size,
          n = target_data$sample_size_per_arm,
          sig.level = alpha,
          type = "one.sample",
          alternative = alternative
        )$power
      } else if (frequentist_test == "z-test") {
        # Calculate the power of the z-test
        power <- pwr::pwr.norm.test(
          d = effect_size,
          n = target_data$sample_size_per_arm,
          sig.level = alpha,
          alternative = alternative
        )$power
      } else {
        stop("Unsupported test type.")
      }

      conf_int_power <- c(power, power)
    } else {
      # In this case, we cannot use an analytical computation of power
      set.seed(simulation_config$seed)

      theta_0 <- case_study_config$theta_0

      # Estimate the frequentist OCs for the model in the scenario considered
      # Generate data for n_replicates clinical trials
      target_data_samples <- target_data$generate(n_replicates)

      test_decisions = numeric(n_replicates)
      for (r in 1:nrow(target_data_samples)) {
        target_data$sample <- target_data_samples[r, ]
        if (frequentist_test == 't-test'){
          test <- BSDA::tsum.test(
            mean.x = target_data$sample$treatment_effect_estimate,
            mu = theta_0,
            alternative = alternative,
            s.x = target_data$sample$standard_deviation,
            n.x = target_data$sample$sample_size_per_arm
          )
        } else {
          stop("Only implemented for a t-test.")
        }
        test_decisions[r] <- test$p.value < alpha
      }
      power <- mean(test_decisions)
      conf_int_power <- binom.test(sum(test_decisions), length(test_decisions), conf.level = 0.95)$conf.int
    }
  } else if (target_data$summary_measure_likelihood == "binomial") {
    # Compute Cohen's h
    h <- pwr::ES.h(target_data$treatment_rate, target_data$control_rate)

    power <- pwr::pwr.2p2n.test(
      h = h,
      n1 = target_data$sample_size_per_arm,
      n2 = target_data$sample_size_per_arm,
      sig.level = alpha,
      alternative = alternative
    )$power

    conf_int_power <- c(power, power)
  } else {
    stop("Unsupported likelihood type.")
  }

  return(list(power = power, conf_int_power = conf_int_power))
}

#' Compute the frequentist power
#'
#' @description This function computes the power of a test for a pooled analysis at a given significance level
#'
#' @param alpha The significance level.
#' @param target_data The target data containing sample size, treatment effect, standard deviation, and summary measure distribution.
#' @param source_data Source study data
#' @param frequentist_test The type of frequentist test ("t-test" or "z-test").
#' @param theta_0 The null hypothesis value.
#' @param target_data Target data object
#' @param frequentist_test Type of frequentist test to apply, either z-test or t-test
#' @param theta_0 Boundary of the null hypothesis space
#' @param null_space Side of the null space, either left or right.
#'
#' @return The power of the test.
#'
#' @export
compute_freq_power_pooling <- function(alpha,
                                       target_data,
                                       source_data,
                                       frequentist_test,
                                       theta_0,
                                       null_space,
                                       simulation_config,
                                       case_study =  NULL,
                                       n_replicates = 1000) {

  if (null_space == "left") {
    alternative <- "greater"
  } else if (null_space == "right") {
    alternative <- "less"
  } else {
    stop("Null space must be either 'left' or 'right'")
  }

  assertions::assert_number(target_data$treatment_effect)
  assertions::assert_number(target_data$standard_deviation)
  assertions::assert_number(target_data$sample_size_per_arm)

  power <- NA # Default value in case of an unsupported distribution

  if (target_data$summary_measure_likelihood == "normal") {
    if (target_data$endpoint == "normal"  || target_data$endpoint == "continuous"  || case_study == "mepolizumab"){
      target_treatment_effect_standard_error <- target_data$standard_deviation / sqrt(target_data$sample_size_per_arm)

      pooled_treatment_effect <- (
        source_data$treatment_effect_estimate / (
          source_data$standard_error ^ 2 / target_treatment_effect_standard_error ^
            2 + 1
        )
      ) + (
        target_data$treatment_effect / (
          1 + target_treatment_effect_standard_error ^ 2 / source_data$standard_error ^
            2
        )
      )


      pooled_standard_error_2 <- 1 / (1 / source_data$standard_error ^ 2 + 1 / target_treatment_effect_standard_error ^
                                        2)

      pooled_variance <- pooled_standard_error_2 * (
        target_data$sample_size_per_arm + source_data$equivalent_source_sample_size_per_arm
      )

      effect_size <- (pooled_treatment_effect - theta_0) / sqrt(pooled_variance)


      if (frequentist_test == "t-test") {
        # Use pwr::pwr.t.test for a t-test power calculation
        power <- pwr::pwr.t.test(
          d = effect_size,
          n = target_data$sample_size_per_arm + source_data$equivalent_source_sample_size_per_arm,
          sig.level = alpha,
          type = "one.sample",
          alternative = alternative
        )$power
      } else if (frequentist_test == "z-test") {
        # Calculate the power of the z-test
        power <- pwr::pwr.norm.test(
          d = effect_size,
          n = target_data$sample_size_per_arm + source_data$equivalent_source_sample_size_per_arm,
          sig.level = alpha,
          alternative = alternative
        )$power
      } else {
        stop("Only implemented for a t-test or z-test.")
      }
      conf_int_power <- c(power, power)
    } else {
      # In this case, we cannot use an analytical computation of power
      set.seed(simulation_config$seed)
      case_study_config <- yaml::read_yaml(paste0(case_studies_config_dir, case_study, ".yml"))


      theta_0 <- case_study_config$theta_0
      null_space <- case_study_config$null_space

      # Estimate the frequentist OCs for the model in the scenario considered
      # Generate data for n_replicates clinical trials
      target_data_samples <- target_data$generate(n_replicates)

      test_decisions = numeric(n_replicates)
      for (r in 1:nrow(target_data_samples)) {
        target_data$sample <- target_data_samples[r, ]

        target_treatment_effect_standard_error <- target_data$sample$standard_deviation / sqrt(target_data$sample$sample_size_per_arm)

        pooled_treatment_effect <- (
          source_data$treatment_effect_estimate / (
            source_data$standard_error ^ 2 / target_treatment_effect_standard_error ^
              2 + 1
          )
        ) + (
          target_data$sample$treatment_effect_estimate / (
            1 + target_treatment_effect_standard_error ^ 2 / source_data$standard_error ^
              2
          )
        )


        pooled_standard_error_2 <- 1 / (1 / source_data$standard_error ^ 2 + 1 / target_treatment_effect_standard_error ^
                                          2)

        pooled_variance <- pooled_standard_error_2 * (
          target_data$sample_size_per_arm + source_data$equivalent_source_sample_size_per_arm
        )

        effect_size <- (pooled_treatment_effect - theta_0) / sqrt(pooled_variance)

        if (frequentist_test == 't-test'){
          if (null_space == "left"){
            alternative = "greater"
          } else if (null_space == "right"){
            alternative = "less"
          }

          test <- BSDA::tsum.test(
            mean.x = pooled_treatment_effect,
            mu = theta_0,
            alternative = alternative,
            s.x = target_data$sample$standard_deviation,
            n.x = target_data$sample$sample_size_per_arm,
          )
        } else {
          stop("Only implemented for a t-test.")
        }
        test_decisions[r] <- test$p.value < alpha
      }

      power <- mean(test_decisions)
      conf_int_power <- binom.test(sum(test_decisions), length(test_decisions), conf.level = 0.95)$conf.int
    }
  } else if (target_data$summary_measure_likelihood == "binomial") {
    pooled_sample_size_treatment <- target_data$sample_size_per_arm + source_data$sample_size_treatment
    treatment_rate_pooled <- (
      target_data$treatment_rate * target_data$sample_size_per_arm + source_data$treatment_rate * source_data$sample_size_treatment
    ) / (pooled_sample_size_treatment)

    pooled_sample_size_control <- target_data$sample_size_per_arm + source_data$sample_size_control
    control_rate_pooled <- (
      target_data$control_rate * target_data$sample_size_per_arm + source_data$control_rate * source_data$sample_size_control
    ) / (pooled_sample_size_control)

    # Compute Cohen's h
    h <- pwr::ES.h(treatment_rate_pooled, control_rate_pooled)
    power <- pwr::pwr.2p2n.test(
      h = h,
      n1 = pooled_sample_size_treatment,
      n2 = pooled_sample_size_control,
      sig.level = alpha,
      alternative = alternative
    )$power

    conf_int_power = c(power, power)
  } else {
    stop("This likelihood is not supported.")
  }

  assertions::assert_number(power)
  return(list(power = power, conf_int_power = conf_int_power))
}


compute_power_with_tie_ci <- function(alpha,
                                      target_data,
                                      frequentist_test,
                                      theta_0,
                                      null_space,
                                      simulation_config,
                                      case_study = NULL,
                                      n_replicates = 100,
                                      n_samples = 100) {
  if (is.na(alpha$conf_int_upper) || is.na(alpha$conf_int_lower)) {
    return(NA)
  }

  # Generate samples of alpha (TIE) based on the confidence interval
  alpha_samples <- rnorm(n_samples, mean = alpha$mean, sd = ((alpha$conf_int_upper - alpha$conf_int_lower) / (2 * 1.96)))

  # Ensure alpha values stay within valid bounds (0, 1)
  alpha_samples <- alpha_samples[alpha_samples > 0 & alpha_samples < 1]

  # Compute power for each sampled alpha
  power_samples <- sapply(alpha_samples, function(alpha) {
    compute_freq_power(
      alpha = alpha,
      target_data = target_data,
      frequentist_test = frequentist_test,
      theta_0 = theta_0,
      null_space = null_space,
      case_study = case_study,
      simulation_config = simulation_config,
      n_replicates = n_replicates
    )$power
  })

  # Compute mean and confidence interval of power
  mean_power <- mean(power_samples, na.rm = TRUE)
  power_ci <- quantile(power_samples, probs = c(0.025, 0.975), na.rm = TRUE)

  # Takes into account the uncertainty on alpha
  return(list(power = mean_power, conf_int_power = power_ci))
}


#' Compute the frequentist power at equivalent tie
#'
#' @description This function computes the frequentist power at equivalent tie for a given set of results and analysis configuration.
#'
#' @param results The results data frame.
#' @param analysis_config The analysis configuration.
#'
#' @return The final results data frame with power and frequentist test columns added.
#'
#' @export
frequentist_power_at_equivalent_tie <- function(results, analysis_config, simulation_config, parallelization = FALSE) {
  if (nrow(results) == 0) {
    stop("The results dataframe is empty.")
  }

  # Remove the tie, mcse_tie, conf_int_tie_lower and conf_int_tie_upper if they exist
  results <- results[, !(
    names(results) %in% c(
      "tie",
      "mcse_tie",
      "conf_int_tie_lower",
      "conf_int_tie_upper",
      "frequentist_power_at_equivalent_tie",
      "frequentist_power_at_equivalent_tie_lower",
      "frequentist_power_at_equivalent_tie_upper",
      "frequentist_test"
    )
  )]

  matching_columns <- c(
    "method",
    "parameters",
    "control_drift",
    "source_denominator",
    "source_denominator_change_factor",
    "case_study",
    "target_to_source_std_ratio",
    "target_sample_size_per_arm",
    "theta_0",
    "null_space",
    "sampling_approximation",
    "summary_measure_likelihood",
    "source_sample_size_treatment",
    "source_sample_size_control",
    "endpoint",
    "source_standard_error",
    "source_treatment_effect_estimate",
    "equivalent_source_sample_size_per_arm"
  )

  for (case_study in unique(results$case_study)){
    if (sum(results[results$case_study == case_study, ]['target_treatment_effect'] == results$theta_0) == 0){
      stop("theta_0 not included among the target study treatment effects considered in the simulation study!")
    }
  }

  # Select the results that correspond to TIE computation.
  results_freq_df_tie <- results[results$target_treatment_effect == results$theta_0, c(
    matching_columns,
    c(
      "success_proba",
      "mcse_success_proba",
      "conf_int_success_proba_lower",
      "conf_int_success_proba_upper"
    )
  )]

  # For these results the probability of success corresponds to TIE, so we rename the columns accordingly.
  results_freq_df_tie <- results_freq_df_tie %>%
    dplyr::rename(
      tie = success_proba,
      mcse_tie = mcse_success_proba,
      conf_int_tie_lower = conf_int_success_proba_lower,
      conf_int_tie_upper = conf_int_success_proba_upper
    )

  # We merge the TIE results onto the corresponding scenarios of the results dataframe.
  results <- dplyr::left_join(results, results_freq_df_tie, by = matching_columns)

  frequentist_test <- analysis_config[["frequentist_test"]]

  if (parallelization == TRUE){
    # Set up parallel backend
    n_cores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)

    required_libraries <- c(
      "pwr", "dplyr", "yaml", "BSDA"
    )

    # Export necessary functions and objects to the cluster
    paths <- .libPaths()
    parallel::clusterExport(cl,
                  varlist = c("paths", "required_libraries", "config_dir", "case_studies_config_dir"),
                  envir = environment())

    # Load required libraries in workers
    parallel::clusterEvalQ(cl, {
      .libPaths(paths)
      sapply(required_libraries, library, character.only = TRUE)
      devtools::load_all()
    })

    # Use foreach for parallel computation
    results_list <- foreach(i = seq_len(nrow(results)), .packages = c("dplyr", "yaml", "pwr", "BSDA")) %dopar% {
      if (is.na(results$tie[i])) {
        warning("TIE is NA")
        return(NULL)
      }

      target_data <- load_data(results[i, ], type = "target", reload_data_objects = TRUE)
      source_data <- load_data(results[i, ], type = "source", reload_data_objects = TRUE)

      alpha <- list(
        mean = results$tie[i],
        conf_int_lower = results$conf_int_tie_lower[i],
        conf_int_upper = results$conf_int_tie_upper[i]
      )

      power_estimation <- compute_power_with_tie_ci(
        alpha = alpha,
        target_data = target_data,
        frequentist_test = frequentist_test,
        theta_0 = results$theta_0[i],
        null_space = results$null_space[i],
        case_study = results[i, ]$case_study,
        simulation_config = simulation_config
      )

      list(
        frequentist_power_at_equivalent_tie = power_estimation$power,
        frequentist_power_at_equivalent_tie_lower = power_estimation$conf_int_power[1],
        frequentist_power_at_equivalent_tie_upper = power_estimation$conf_int_power[2],
        frequentist_test = frequentist_test
      )
    }

    # Stop parallel backend
    parallel::stopCluster(cl)

    # Combine results into the dataframe
    results <- cbind(results, do.call(rbind, results_list))
  } else {
    # Progress bar function in R
    progress_bar <- function(n) {
      pb <- txtProgressBar(min = 0,
                           max = n,
                           style = 3)
      return(function(i) {
        setTxtProgressBar(pb, i)
      })
    }

    frequentist_test <- analysis_config[["frequentist_test"]]
    # Iterate through rows and compute power at equivalent TIE.
    for (i in seq_len(nrow(results))) {
      target_data <- load_data(results[i, ],
                               type = "target",
                               reload_data_objects = TRUE)

      source_data <- load_data(results[i, ],
                               type = "source",
                               reload_data_objects = TRUE)

      if (is.na(results$tie[i])){
        warning("TIE is NA")
        next
      }

      alpha = list(mean = results$tie[i], conf_int_lower = results$conf_int_tie_lower[i], conf_int_upper = results$conf_int_tie_upper[i])

      power_estimation <- compute_power_with_tie_ci (
        alpha = alpha,
        target_data = target_data,
        frequentist_test = frequentist_test,
        theta_0 = results$theta_0[i],
        null_space = results$null_space[i],
        case_study =  results[i, ]$case_study,
        simulation_config = simulation_config
      )

      results$frequentist_power_at_equivalent_tie[i] <- power_estimation$power
      results$frequentist_power_at_equivalent_tie_lower[i] <- power_estimation$conf_int_power[1]
      results$frequentist_power_at_equivalent_tie_upper[i] <- power_estimation$conf_int_power[2]

    results$frequentist_test[i] <- frequentist_test

    # Update progress bar
    pb <- progress_bar(nrow(results))(i)
    }
  }

  return(results)
}


frequentist_power_at_nominal_tie <- function(results, analysis_config, simulation_config) {
  if (nrow(results) == 0) {
    stop("The results dataframe is empty.")
  }

  results <- results[, !(
    names(results) %in% c(
     "nominal_frequentist_power_separate",
     "nominal_frequentist_power_pooling"
    )
  )]

  nominal_tie <- analysis_config[["nominal_tie"]]
  frequentist_test <- analysis_config[["frequentist_test"]]

  # Initialize the new columns
  results$nominal_frequentist_power_separate <- NA
  results$nominal_frequentist_power_pooling <- NA

  # Initialize the new columns
  results$nominal_frequentist_power_separate <- NA
  results$nominal_frequentist_power_pooling <- NA

  # Progress bar function in R
  progress_bar <- function(n) {
    pb <- txtProgressBar(min = 0,
                         max = n,
                         style = 3)
    return(function(i) {
      setTxtProgressBar(pb, i)
    })
  }

  # Iterate through rows and compute power
  for (i in seq_len(nrow(results))) {
    target_data <- load_data(results[i, ],
                             type = "target",
                             reload_data_objects = TRUE)
    source_data <- load_data(results[i, ],
                             type = "source",
                             reload_data_objects = TRUE)



    nominal_frequentist_power_separate <- compute_freq_power(
      alpha = nominal_tie,
      target_data = target_data,
      frequentist_test = frequentist_test,
      theta_0 = results$theta_0[i],
      null_space = results$null_space[i],
      case_study =  results[i, ]$case_study,
      simulation_config = simulation_config
    )

    results$nominal_frequentist_power_separate[i] <- nominal_frequentist_power_separate$power

    nominal_frequentist_power_pooling <- compute_freq_power_pooling(
      alpha = nominal_tie,
      target_data = target_data,
      source_data = source_data,
      frequentist_test = frequentist_test,
      theta_0 = results$theta_0[i],
      null_space = results$null_space[i],
      case_study =  results[i, ]$case_study,
      simulation_config = simulation_config
    )

    results$nominal_frequentist_power_pooling[i] <- nominal_frequentist_power_pooling$power

    # Update progress bar
    pb <- progress_bar(nrow(results))(i)
  }
  return(results)
}
