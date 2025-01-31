#' @title TestThenPool class
#' @description This class implements the Test-Then-Pool framework for Bayesian borrowing in clinical trials.
#' @details The TestThenPool class provides methods for testing and inference with the Test-then-Pool method.
#' @field pooling Pooling model
#' @field separate Separate model
#' @field pool Pooling indicator
#' @field inference_method Inference method
#' @field summary_measure_likelihood Summary measure distribution
#' @field source_treatment_effect_estimate Source treatment effect estimate
#' @field source_standard_error Source standard error
#' @field method Method name
#' @field empirical_bayes Indicator that the method uses empirical Bayes
#' @export
TestThenPool <- R6::R6Class(
  "TestThenPool",
  inherit = Model,
  list(
    pooling = NULL,
    separate = NULL,
    pool = NULL,
    inference_method = NULL,
    summary_measure_likelihood = NULL,
    source_treatment_effect_estimate = NULL,
    source_standard_error = NULL,
    method = "test_then_pool",
    empirical_bayes = TRUE,

    #' @description Initializes a TestThenPool object.
    #' @param prior The prior distribution for the treatment effect.
    #' @param mcmc_config Configuration for the MCMC sampling.
    initialize = function(prior, mcmc_config = NULL) {
      super$initialize()
      self$pool <- NULL
      self$summary_measure_likelihood <- prior$source$summary_measure_likelihood

      if (self$summary_measure_likelihood == "binomial") {
        if (is.null(mcmc_config)) {
          stop("mcmc_config is null.")
        }
        self$separate <- BinomialSeparate$new(prior = prior, mcmc_config = mcmc_config)
        self$pooling <- BinomialPooling$new(prior = prior, mcmc_config = mcmc_config)
      } else if (self$summary_measure_likelihood == "normal") {
        self$separate <- SeparateGaussian$new(prior = prior)
        self$pooling <- PoolGaussian$new(prior = prior)
      } else {
        stop("Distribution not supported")
      }
      self$pooling$prior <- prior
      self$separate$prior <- prior

      self$source_treatment_effect_estimate <- prior$source$treatment_effect_estimate
      self$source_standard_error <- prior$source$standard_error
    },

    #' @description Performs the test with the Test-then-Pool method.
    #' @param target_data The data from the target study.
    test = function(target_data) {
      stop("Subclasses must implement the 'test' method.")
    },

    #' @description p-value of the test
    #' @param target_data The data from the target study.
    test_pvalue = function(target_data) {
      stop("Subclasses must implement the 'test' method.")
    },

    #' @description Performs inference with the Test-then-Pool method.
    #' @param target_data The data from the target study.
    inference = function(target_data) {
      self$test(target_data) # This will set the value of self$pool, and determine whether we should use pooling or separate
      assertions::assert_logical(self$pool)
      if (length(self$pool) == 0) {
        stop("Test result is invalid.")
      }

      if (self$pool) {
        self$prior_mean <- self$pooling$prior_mean
        self$prior_var <- self$pooling$prior_var

        # Perform a pooled analysis
        self$pooling$inference(target_data)
        self$post_mean <- self$pooling$post_mean
        self$post_var <- self$pooling$post_var
      } else {
        self$prior_mean <- self$separate$prior_mean
        self$prior_var <- self$separate$prior_var

        # Perform a separate analysis
        self$separate$inference(target_data)
        self$post_mean <- self$separate$post_mean
        self$post_var <- self$separate$post_var
      }

      if (is.numeric((self$post_mean)) &&
          is.numeric((self$post_var))) {
        return("Success")
      } else {
        stop("Non-numeric moments")
      }
    },
    #' @description Performs the empirical Bayes update with the Test-then-Pool method.
    #' @param target_data The data from the target study.
    empirical_bayes_update = function(target_data) {
      if (self$pool) {
        self$pooling$empirical_bayes_update(target_data)
      } else {
        self$separate$empirical_bayes_update(target_data)
      }
    },

    #' @description Calculates the credible interval with the Test-then-Pool method.
    #' @param level The confidence level for the credible interval (default is 0.95).
    credible_interval = function(level = 0.95) {
      if (self$pool) {
        return(self$pooling$credible_interval(level))
      } else {
        return(self$separate$credible_interval(level))
      }
    },
    #' @description Return the median of the posterior distribution.
    #' @param ... Optional argument
    posterior_median = function(...) {
      if (self$pool) {
        return(self$pooling$posterior_median(...))
      } else {
        return(self$separate$posterior_median(...))
      }
    },

    #' @description Calculates the posterior probability density function with the Test-then-Pool method.
    #' @param target_treatment_effect The treatment effect of interest.
    posterior_pdf = function(target_treatment_effect) {
      if (self$pool) {
        return(self$pooling$posterior_pdf(target_treatment_effect))
      } else {
        return(self$separate$posterior_pdf(target_treatment_effect))
      }
    },

    #' @description Calculates the prior probability density function with the Test-then-Pool method.
    #' @param target_treatment_effect The treatment effect of interest.
    prior_pdf = function(target_treatment_effect) {
      if (self$pool) {
        return(self$pooling$prior_pdf(target_treatment_effect))
      } else {
        return(self$separate$prior_pdf(target_treatment_effect))
      }
    },

    #' @description Calculates the posterior cumulative distribution function with the Test-then-Pool method.
    #' @param target_treatment_effect The treatment effect of interest.
    posterior_cdf = function(target_treatment_effect) {
      if (self$pool) {
        return(self$pooling$posterior_cdf(target_treatment_effect))
      } else {
        return(self$separate$posterior_cdf(target_treatment_effect))
      }
    },
    #' @description CDF of the prior distribution
    #' @param target_treatment_effect Treatment effect value in the target study.
    prior_cdf = function(target_treatment_effect) {
      if (self$pool) {
        return(self$pooling$prior_cdf(target_treatment_effect))
      } else {
        return(self$separate$prior_cdf(target_treatment_effect))
      }
    },

    #' @description Samples from the prior distribution with the Test-then-Pool method.
    #' @param n_samples The number of samples to generate.
    sample_prior = function(n_samples) {
      if (self$pool) {
        return(self$pooling$sample_prior(n_samples))
      } else {
        return(self$separate$sample_prior(n_samples))
      }
    },

    #' @description Samples from the posterior distribution with the Test-then-Pool method.
    #' @param n_samples The number of samples to generate.
    sample_posterior = function(n_samples) {
      if (self$pool) {
        return(self$pooling$sample_posterior(n_samples))
      } else {
        return(self$separate$sample_posterior(n_samples))
      }
    },

    #' @description Converts the prior distribution to the RBesT format.
    #' @param ... Optional argument
    #' @return None
    prior_to_RBesT = function(...) {
      # test-then-pool can be seen as a form of empirical Bayes, therefore the prior is only fully defined once the data is observed.
      if (self$pool) {
        self$RBesT_prior <- self$pooling$prior_to_RBesT(...)
      } else {
        self$RBesT_prior <- self$separate$prior_to_RBesT(...)
      }
      self$RBesT_prior_normix <- self$RBesT_prior
    },

    #' @description Convert the posterior distribution to RBesT format
    #' @param target_data Target study data
    #' @param simulation_config Simulation configuration
    posterior_to_RBesT = function(target_data, simulation_config) {
      if (self$pool) {
        self$pooling$posterior_to_RBesT(target_data, simulation_config)
        self$RBesT_posterior <- self$pooling$RBesT_posterior
        self$RBesT_posterior_normix <- self$pooling$RBesT_posterior_normix
      } else {
        self$separate$posterior_to_RBesT(target_data, simulation_config)
        self$RBesT_posterior <- self$separate$RBesT_posterior
        self$RBesT_posterior_normix <- self$separate$RBesT_posterior_normix
      }
    },

    #' @description Return the test decision based on the posterior distribution. The decision rule is: \eqn{P(\theta_T > \theta_0 \mid \mathbf{D}_S, \mathbf{D}_T) > ) \eta} (if null_space is right).
    #' @param critical_value Critical value
    #' @param theta_0 Boundary of the null hypothesis space
    #' @param null_space Side of the null hypothesis space
    #' @param confidence_level Confidence level of the test
    test_decision = function(critical_value,
                             theta_0,
                             null_space,
                             confidence_level) {
      if (self$pool) {
        return(
          self$pooling$test_decision(
            critical_value = critical_value,
            theta_0 = theta_0,
            null_space = null_space,
            confidence_level = confidence_level
          )
        )
      } else {
        return(
          self$separate$test_decision(
            critical_value = critical_value,
            theta_0 = theta_0,
            null_space = null_space,
            confidence_level = confidence_level
          )
        )
      }
    },

    #' @description
    #' Plot the pooling test decision as a function of drift in treatment effect
    #' @param source_treatment_effect_estimate Treatment effect estimate in the source study
    #' @param target_data Target study data
    #' @param min_drift Minimum drift value
    #' @param max_drift Maximum drift value
    #' @param resolution Number of points on the drift grid.
    #' @return The prior PDF.
    plot_test_vs_drift = function(source_treatment_effect_estimate,
                                  target_data,
                                  min_drift,
                                  max_drift,
                                  resolution) {
      drifts <- seq(min_drift, max_drift, length.out = resolution)

      test_results <- sapply(drifts, function(drift) {
        target_data_drift <- duplicate(target_data)
        target_data_drift$sample$treatment_effect_estimate <- source_treatment_effect_estimate + drift
        self$test(target_data_drift)
      })

      plot_data <- data.frame(Drift = drifts, Test = test_results * 1)

      plt <- ggplot(plot_data, aes(x = Drift, y = Test)) +
        geom_line(color = "blue") +
        labs(title = "Pooling Test a Function of Drift", x = "Drift", y = "Pooling or not") +
        theme_minimal()

      return(plt)
    },
    #' @description
    #' Plot the pooling test p-value as a function of drift in treatment effect
    #' @param source_treatment_effect_estimate Treatment effect estimate in the source study
    #' @param target_data Target study data
    #' @param min_drift Minimum drift value
    #' @param max_drift Maximum drift value
    #' @param resolution Number of points on the drift grid.
    #' @return The prior PDF.
    plot_test_pvalue_vs_drift = function(source_treatment_effect_estimate,
                                         target_data,
                                         min_drift,
                                         max_drift,
                                         resolution) {
      drifts <- seq(min_drift, max_drift, length.out = resolution)

      p_values <- sapply(drifts, function(drift) {
        target_data_drift <- duplicate(target_data)
        target_data_drift$sample$treatment_effect_estimate <- source_treatment_effect_estimate + drift
        self$test_pvalue(target_data_drift)
      })

      plot_data <- data.frame(Drift = drifts, Pvalue = p_values)

      plt <- ggplot(plot_data, aes(x = Drift, y = Pvalue)) +
        geom_line(color = "blue") +
        labs(title = "P-value a Function of Drift", x = "Drift", y = "p-value") +
        theme_minimal()

      return(plt)
    }
  )
)

#' @title TestThenPoolEquivalence class
#' @description This class extends the TestThenPool class to implement the Test-Then-Pool framework for equivalence trials.
#' @details The TestThenPoolEquivalence class provides methods for testing and inference in equivalence trials using the Test-Then-Pool framework.
#' @field significance_level Significance level for the test
#' @field equivalence_margin Equivalence margin for the test
#' @field method Method name
#' @export
TestThenPoolEquivalence <- R6::R6Class(
  classname = "TestThenPoolEquivalence",
  inherit = TestThenPool,
  public = list(
    significance_level = NULL,
    equivalence_margin = NULL,
    method = "test_then_pool_equivalence",
    #' @description Initializes a TestThenPoolEquivalence object.
    #' @param prior The prior distribution for the treatment effect.
    #' @param mcmc_config Configuration for the MCMC sampling.
    initialize = function(prior, mcmc_config = NULL) {
      super$initialize(prior = prior, mcmc_config = mcmc_config)
      self$significance_level <- unlist(prior$method_parameters$significance_level)
      self$equivalence_margin <- unlist(prior$method_parameters$equivalence_margin)

      if (is.null(self$significance_level)) {
        stop("Significance level is not defined.")
      }

      if (is.null(self$equivalence_margin)) {
        stop("Equivalence margin is not defined.")
      }
    },

    #' @description Performs the test in the Test-Then-Pool framework for equivalence.
    #' @param target_data The data from the target study.
    #' @param test_type Test type
    test_pvalue = function(target_data, test_type = "t-test") {
      assertions::assert_number(target_data$sample$treatment_effect_estimate)
      assertions::assert_number(target_data$sample$treatment_effect_standard_error)
      assertions::assert_number(target_data$sample_size_per_arm)

      if (!(
        self$summary_measure_likelihood == "normal" ||
        self$summary_measure_likelihood == "binomial"
      )) {
        stop("Distribution not supported for this method.")
      }

      if (test_type == "z-test") {
        # H0a : θS - θT > λ
        left <- BSDA::zsum.test(
          mean.x = self$prior$source$treatment_effect_estimate,
          mean.y = target_data$sample$treatment_effect_estimate,
          mu = self$equivalence_margin,
          alternative = "less",
          sigma.x = self$prior$source$standard_error * sqrt(
            self$prior$source$equivalent_source_sample_size_per_arm
          ),
          sigma.y = target_data$sample$standard_deviation,
          n.x = self$prior$source$equivalent_source_sample_size_per_arm,
          n.y = target_data$sample_size_per_arm
        )

        # Equivalent to :
        # z_left <- (
        #   self$source_treatment_effect_estimate - target_data$sample$treatment_effect_estimate - self$equivalence_margin
        # ) / std
        # p_left <- pnorm(z_left, lower.tail = TRUE)
        #
        # std <- sqrt(
        #   target_data$sample$treatment_effect_standard_error ^ 2 + self$source_standard_error ^
        #     2
        # )

        # H0b : θS - θT < -λ
        right <- BSDA::zsum.test(
          mean.x = self$prior$source$treatment_effect_estimate,
          mean.y = target_data$sample$treatment_effect_estimate,
          mu = -self$equivalence_margin,
          alternative = "greater",
          sigma.x = self$prior$source$standard_error * sqrt(
            self$prior$source$equivalent_source_sample_size_per_arm
          ),
          sigma.y = target_data$sample$standard_deviation,
          n.x = self$prior$source$equivalent_source_sample_size_per_arm,
          n.y = target_data$sample_size_per_arm
        )

        # Equivalent to :
        # z_right <- (
        #   self$source_treatment_effect_estimate - target_data$sample$treatment_effect_estimate + self$equivalence_margin
        # ) / std
        # p_right <- pnorm(z_right, lower.tail = FALSE)
      } else if (test_type == "t-test") {
        # H0a : θS - θT > λ
        left <- BSDA::tsum.test(
          mean.x = self$prior$source$treatment_effect_estimate,
          mean.y = target_data$sample$treatment_effect_estimate,
          mu = self$equivalence_margin,
          alternative = "less",
          s.x = self$prior$source$standard_error * sqrt(
            self$prior$source$equivalent_source_sample_size_per_arm
          ),
          s.y = target_data$sample$standard_deviation,
          n.x = self$prior$source$equivalent_source_sample_size_per_arm,
          n.y = target_data$sample_size_per_arm
        )

        # H0b : θS - θT < -λ
        right <- BSDA::tsum.test(
          mean.x = self$prior$source$treatment_effect_estimate,
          mean.y = target_data$sample$treatment_effect_estimate,
          mu = -self$equivalence_margin,
          alternative = "greater",
          s.x = self$prior$source$standard_error * sqrt(
            self$prior$source$equivalent_source_sample_size_per_arm
          ),
          s.y = target_data$sample$standard_deviation,
          n.x = self$prior$source$equivalent_source_sample_size_per_arm,
          n.y = target_data$sample_size_per_arm
        )
      } else {
        stop("Not implemented for this test. Must be either 't-test' or 'z-test'.")
      }

      # Take the maximum p_value of the two tests
      p_value <- pmax(left$p.value, right$p.value)

      return(p_value)
    },
    #' @param target_data Target study data
    #' @param test_type Frequentist test
    test = function(target_data, test_type = "t-test") {
      p_value <- self$test_pvalue(target_data, test_type = test_type)

      reject_null <- p_value < self$significance_level

      # If the null is rejected, this means the data can be pooled
      self$pool <- reject_null

      self$posterior_parameters$pool <- self$pool

      assertions::assert_logical(self$pool)
      if (length(self$pool) == 0) {
        stop("Test result is invalid.")
      }
      return(self$pool)
    }
  )
)

#' @title TestThenPoolDifference class
#' @description This class extends the TestThenPool class to implement the Test-Then-Pool framework for difference trials.
#' @details The TestThenPoolDifference class provides methods for testing and inference in difference trials using the Test-Then-Pool framework.
#' @field significance_level Significance level for the test
#' @field method Method name
#' @export
TestThenPoolDifference <- R6::R6Class(
  "TestThenPoolDifference",
  inherit = TestThenPool,
  public = list(
    significance_level = NULL,
    method = "test_then_pool_difference",
    #' @description
    #' Instantiate model
    #' @param prior Prior
    #' @param mcmc_config MCMC configuration
    initialize = function(prior, mcmc_config = NULL) {
      super$initialize(prior = prior, mcmc_config = mcmc_config)
      self$significance_level <- prior$method_parameters$significance_level[[1]]

      if (is.null(self$significance_level)) {
        stop("Significance level is not defined.")
      }
    },
    #' @param target_data Target study data
    #' @param test_type Frequentist test
    test_pvalue = function(target_data, test_type = "t-test") {
      assertions::assert_number(target_data$sample$treatment_effect_estimate)
      assertions::assert_number(target_data$sample$treatment_effect_standard_error)
      assertions::assert_number(target_data$sample_size_per_arm)

      if (!(
        self$summary_measure_likelihood == "normal" ||
        self$summary_measure_likelihood == "binomial"
      )) {
        stop("Distribution not supported for this method.")
      }

      # Null Hypothesis (H0): There is no difference between the means of the two studies
      # H0 : θS - θT = 0
      if (test_type == "z-test") {
        test_result <- BSDA::zsum.test(
          mean.x = self$prior$source$treatment_effect_estimate,
          mean.y = target_data$sample$treatment_effect_estimate,
          mu = 0,
          alternative = "two.sided",
          sigma.x = self$prior$source$standard_error * sqrt(
            self$prior$source$equivalent_source_sample_size_per_arm
          ),
          sigma.y = target_data$sample$standard_deviation,
          n.x = self$prior$source$equivalent_source_sample_size_per_arm,
          n.y = target_data$sample_size_per_arm
        )
      } else if (test_type == "t-test") {
        test_result <- BSDA::tsum.test(
          mean.x = self$prior$source$treatment_effect_estimate,
          mean.y = target_data$sample$treatment_effect_estimate,
          mu = 0,
          alternative = "two.sided",
          s.x = self$prior$source$standard_error * sqrt(
            self$prior$source$equivalent_source_sample_size_per_arm
          ),
          s.y = target_data$sample$standard_deviation,
          n.x = self$prior$source$equivalent_source_sample_size_per_arm,
          n.y = target_data$sample_size_per_arm
        )
      } else {
        stop("Not implemented for this test. Must be either 't-test' or 'z-test'.")
      }

      return(test_result$p.value)
    },
    #' @param target_data Target study data
    #' @param test_type Frequentist test
    test = function(target_data, test_type = "t-test") {
      p_value <- self$test_pvalue(target_data, test_type = test_type)

      reject_null <- p_value < self$significance_level

      # If the null is rejected, this means the data can be pooled
      self$pool <- !reject_null

      self$posterior_parameters$pool <- self$pool

      assertions::assert_logical(self$pool)
      if (length(self$pool) == 0) {
        stop("Test result is invalid.")
      }
      return(self$pool)
    }
  )
)
