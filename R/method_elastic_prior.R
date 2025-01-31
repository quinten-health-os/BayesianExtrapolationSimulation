#' GaussianElasticPriorLogistic class
#'
#' @description This class represents a Gaussian Elastic Prior model for logistic regression.
#' It inherits from the ConjugateGaussian class.
#'
#' @field tuning_param_a Numeric. Tuning parameter a for the elastic function.
#' @field tuning_param_b Numeric. Tuning parameter b for the elastic function.
#' @field target_standard_error Numeric. Standard error of the target study.
#' @field target_treatment_effect_estimate Numeric. Treatment effect estimate of the target study.
#' @field target_sample_size_per_arm Numeric. Sample size per arm in the target study.
#' @field source_sample_size_per_arm Numeric. Sample size per arm in the source study.
#' @field source_data_variance Numeric. Variance of the source data.
#' @field method Character. Method used for the analysis (default: "elastic_prior").
#' @field shape_param Numeric. Shape parameter for the elastic function.
#' @field prior List. Prior information for the analysis.
#'
#' @export
GaussianElasticPriorLogistic <- R6::R6Class(
  "GaussianElasticPriorLogistic",
  inherit = ConjugateGaussian,
  public = list(
    tuning_param_a = NULL,
    tuning_param_b = NULL,
    target_standard_error = NULL,
    target_treatment_effect_estimate = NULL,
    target_sample_size_per_arm = NULL,
    source_sample_size_per_arm = NULL,
    source_data_variance = NULL,
    method = "elastic_prior",
    shape_param = NULL,
    prior = NULL,

    #' @description Create an object from the GaussianElasticPriorLogistic class
    #'
    #' @param prior List. The prior information for the analysis.
    #' @param summary_measure_likelihood List. Summary measure likelihood.
    #'
    #' @return A new `GaussianElasticPriorLogistic` object.
    initialize = function(prior, summary_measure_likelihood) {
      if (!(prior$method_parameters$initial_prior[[1]] == "noninformative")) {
        stop("Only implemented for a noninformative initial prior")
      }
      super$initialize(prior, summary_measure_likelihood)

      self$source_sample_size_per_arm <- prior$source$equivalent_source_sample_size_per_arm

      self$source_data_variance <- prior$source$standard_error ^ 2 * self$source_sample_size_per_arm

      self$shape_param <- prior$method_parameters$shape_param[[1]]
    },

    #' @description Updates the vague prior variance based on empirical Bayes approach.
    #'
    #' @param target_data List. The target data for the analysis.
    #'
    #' @return None. Updates object fields internally.
    empirical_bayes_update = function(target_data) {
      target_sample_size_per_arm <- target_data$sample_size_per_arm
      target_standard_error <- target_data$sample$treatment_effect_standard_error
      target_treatment_effect_estimate <- target_data$sample$treatment_effect_estimate

      if (is.null(self$target_sample_size_per_arm) ||
          target_sample_size_per_arm != self$target_sample_size_per_arm ||
          is.null(self$target_standard_error) ||
          target_standard_error != self$target_standard_error) {
        self$method_calibration(target_data = target_data)
      }

      self$target_sample_size_per_arm <- target_sample_size_per_arm
      self$target_standard_error <- target_standard_error
      self$target_treatment_effect_estimate <- target_treatment_effect_estimate

      self$setup_prior(target_data)
    },

    #' @description Calibrates the method by determining tuning parameters.
    #'
    #' @param target_data List. The target data for the analysis.
    #'
    #' @return None. Updates tuning parameters internally.
    method_calibration = function(target_data) {
      ## obtain calibrated a and b in smooth elastic function
      tuning_params <- self$determine_tuning_parameters(target_data = target_data)
      self$tuning_param_a <- tuning_params$a
      self$tuning_param_b <- tuning_params$b
    },

    #' @description Determines tuning parameters a and b for the logistic elastic function.
    #'
    #' @param target_data List. The target data for the analysis.
    #'
    #' @return List with components:
    #' \itemize{
    #'   \item{a}{Numeric. Tuning parameter a in elastic function.}
    #'   \item{b}{Numeric. Tuning parameter b in elastic function.}
    #' }
    determine_tuning_parameters = function(target_data) {
      shape_param <- self$shape_param
      clinically_meaningful_diff <- self$prior$method_parameters$clinically_meaningful_diff[[1]]
      percentile_homogeneous <- self$prior$method_parameters$percentile_homogeneous[[1]]
      percentile_heterogeneous <- self$prior$method_parameters$percentile_heterogeneous[[1]]
      value_elastic_heterogeneous <- self$prior$method_parameters$value_elastic_heterogeneous[[1]]
      value_elastic_homogeneous <- self$prior$method_parameters$value_elastic_homogeneous[[1]]
      num_simulations <- self$prior$method_parameters$num_simulations[[1]]

      source_treatment_effect_estimate <- self$prior$source$treatment_effect_estimate

      target_treatment_effect_values <- c(
        source_treatment_effect_estimate,
        source_treatment_effect_estimate + clinically_meaningful_diff,
        source_treatment_effect_estimate - clinically_meaningful_diff
      )
      congruence_measures <- matrix(NA,
                                    num_simulations,
                                    length(target_treatment_effect_values))

      target_sample_size_per_arm <- target_data$sample$sample_size_per_arm

      target_data_variance <- target_data$sample$treatment_effect_standard_error ^
        2 * target_sample_size_per_arm
      target_standard_error <- target_data$sample$treatment_effect_standard_error

      for (i in 1:num_simulations) {
        for (j in 1:length(target_treatment_effect_values)) {
          target_treatment_effect_estimate <- rnorm(1,
                                                    target_treatment_effect_values[j],
                                                    target_standard_error)

          pooled_variance <- ((self$source_sample_size_per_arm - 1) * self$source_data_variance + (target_sample_size_per_arm - 1) * target_data_variance
          ) / (self$source_sample_size_per_arm + target_sample_size_per_arm - 2)

          congruence_measures[i, j] <- max(self$source_sample_size_per_arm,
                                           target_sample_size_per_arm) ^ (-1 / 4) * abs(source_treatment_effect_estimate - target_treatment_effect_estimate) / (
                                             sqrt(
                                               pooled_variance / self$source_sample_size_per_arm + pooled_variance / target_sample_size_per_arm
                                             )
                                           )
        }
      }
      quantile_homogeneous <- quantile(congruence_measures[, 1], percentile_homogeneous)
      quantile_heterogeneous1 <- quantile(congruence_measures[, 2], percentile_heterogeneous)
      quantile_heterogeneous2 <- quantile(congruence_measures[, 3], percentile_heterogeneous)
      KS_homogeneous <- quantile_homogeneous
      KS_heterogeneous <- min(quantile_heterogeneous1, quantile_heterogeneous2)
      b <- log((1 - value_elastic_homogeneous) * value_elastic_heterogeneous / ((1 - value_elastic_heterogeneous) * value_elastic_homogeneous
      )
      ) / ((log(KS_homogeneous)) ^ shape_param - (log(KS_heterogeneous)) ^ shape_param)
      a <- log((1 - value_elastic_homogeneous) / value_elastic_homogeneous) - b * (log(KS_homogeneous)) ^
        shape_param
      return(list(a = a, b = b))
    },

    #' @description Sets up the prior based on target data.
    #'
    #' @param target_data List. The target data for the analysis.
    #'
    #' @return None. Updates prior variance internally.
    setup_prior = function(target_data) {
      target_treatment_effect_estimate <- target_data$sample$treatment_effect_estimate
      target_standard_error <- target_data$sample$treatment_effect_standard_error
      target_sample_size_per_arm <- target_data$sample_size_per_arm

      target_data_variance <- target_standard_error * target_sample_size_per_arm

      # calculate statistic and g(t)
      pooled_variance <- ((self$source_sample_size_per_arm - 1) * self$source_data_variance + (target_sample_size_per_arm - 1) * target_data_variance
      ) / (self$source_sample_size_per_arm + target_sample_size_per_arm - 2)

      congruence_measure <- max(self$source_sample_size_per_arm,
                                target_sample_size_per_arm) ^ (-1 / 4) * abs(self$prior$source$treatment_effect_estimate - target_treatment_effect_estimate) / (
                                  sqrt(
                                    pooled_variance / self$source_sample_size_per_arm + pooled_variance / target_sample_size_per_arm
                                  )
                                )

      elastic_function_value <- 1 / (1 + exp(
        self$tuning_param_a + self$tuning_param_b * (log(congruence_measure)) ^
          self$shape_param
      ))

      if (elastic_function_value == 0) {
        elastic_function_value <- 0.00001
      }

      self$prior_var <- self$prior_var / elastic_function_value
    }
  )
)

#' GaussianElasticPriorStep class
#'
#' @description This class represents an Elastic Prior model.
#' It inherits from the Model class.
#'
#' @field congruence_threshold Congruence threshold
#' @field target_standard_error Standard error on the treatment effect in the target study
#' @field target_treatment_effect_estimate Estimate of the treatment effect in the target study
#' @field target_sample_size_per_arm Sample size per arm in the target study
#' @field source_sample_size_per_arm Sample size per arm in the source study
#' @field source_data_variance Estimated sampling variance in the source study
#' @field method Name of the method
#'
#' @return An R6 class object representing a GaussianElasticPriorStep model
#'
#' @export
#'
GaussianElasticPriorStep <- R6::R6Class(
  "GaussianElasticPriorStep",
  inherit = ConjugateGaussian,
  public = list(
    congruence_threshold = NULL,
    target_standard_error = NULL,
    target_treatment_effect_estimate = NULL,
    target_sample_size_per_arm = NULL,
    source_sample_size_per_arm = NULL,
    source_data_variance = NULL,
    method = "elastic_prior",
    #' @description Initializes the GaussianElasticPriorStep object.
    #'
    #' @param prior The prior information for the analysis.
    initialize = function(prior) {
      if (!(prior$method_parameters$initial_prior[[1]] == "noninformative")) {
        stop("Only implemented for a noninformative initial prior")
      }
      super$initialize(prior)
      self$method <- prior$method_parameters$variant[[1]]

      # The implementation does not handle cases with different sample sizes in each arm, therefore we compute an equivalent sample size per arm
      self$source_sample_size_per_arm <- self$prior$source$equivalent_source_sample_size_per_arm

      self$source_data_variance <- self$prior$source$standard_error ^ 2 * self$source_sample_size_per_arm
    },

    #' @description Updates the vague prior variance based on empirical Bayes approach.
    #'
    #' @param target_data The target data for the analysis.
    empirical_bayes_update = function(target_data) {
      target_sample_size_per_arm <- target_data$sample_size_per_arm
      target_standard_error <- target_data$sample$treatment_effect_standard_error
      target_treatment_effect_estimate <- target_data$sample$treatment_effect_estimate

      if (is.null(self$target_sample_size_per_arm) ||
          target_sample_size_per_arm != self$target_sample_size_per_arm ||
          is.null(self$target_standard_error) ||
          target_standard_error != self$target_standard_error ||
          is.null(self$target_treatment_effect_estimate) ||
          target_treatment_effect_estimate != self$target_treatment_effect_estimate) {
        self$method_calibration(
          target_sample_size_per_arm,
          target_standard_error,
          target_treatment_effect_estimate
        )
      }
      self$target_sample_size_per_arm <- target_sample_size_per_arm
      self$target_standard_error <- target_standard_error
      self$target_treatment_effect_estimate <- target_treatment_effect_estimate

      self$setup_prior(target_data)
    },

    #' @description Calibrates the method based on target data.
    #'
    #' @param target_sample_size_per_arm Sample size per arm in the target study
    #' @param target_standard_error Standard error of the treatment effect in the target study
    #' @param target_treatment_effect_estimate Estimate of the treatment effect in the target study
    method_calibration = function(target_sample_size_per_arm,
                                  target_standard_error,
                                  target_treatment_effect_estimate) {
      target_data_variance <- target_standard_error ^ 2 * target_sample_size_per_arm # sample variance (estimate of the true sampling variance)
      congruence_measures <- self$generate_congruence_measure_distribution(
        target_sample_size_per_arm = target_sample_size_per_arm,
        target_treatment_effect_estimate = target_treatment_effect_estimate,
        target_data_variance = target_data_variance
      )

      self$congruence_threshold <- quantile(congruence_measures,
                                            self$prior$method_parameters$congruence_quantile[[1]])
    },

    #' @description Generate distribution of congruence measure T in the homogeneous case
    #'
    #' @param target_sample_size_per_arm Sample size per arm in the target study
    #' @param target_treatment_effect_estimate Treatment effect estimate in the target study
    #' @param target_data_variance Estimated sampling variance of the target data
    #' @return A numeric vector of congruence measures
    generate_congruence_measure_distribution = function(target_sample_size_per_arm,
                                                        target_treatment_effect_estimate,
                                                        target_data_variance) {
      source_sample_size_per_arm <- self$source_sample_size_per_arm
      source_treatment_effect_estimate <- self$prior$source$treatment_effect_estimate
      source_data_variance <- self$source_data_variance
      num_simulations <- self$prior$method_parameters$num_simulations[[1]]

      congruence_measures <- numeric(num_simulations)
      for (i in 1:num_simulations) {
        pooled_variance <- ((source_sample_size_per_arm - 1) * source_data_variance + (target_sample_size_per_arm - 1) * target_data_variance
        ) / (source_sample_size_per_arm + target_sample_size_per_arm - 2)

        congruence_measures[i] <- max(source_sample_size_per_arm, target_sample_size_per_arm) ^
          (-1 / 4) * abs(source_treatment_effect_estimate - target_treatment_effect_estimate) / (
            sqrt(
              pooled_variance / source_sample_size_per_arm + pooled_variance / target_sample_size_per_arm
            )
          )
      }
      return(congruence_measures)
    },

    #' @description Set up the prior based on target data
    #'
    #' @param target_data The target data for the analysis
    setup_prior = function(target_data) {
      target_treatment_effect_estimate <- target_data$sample$treatment_effect_estimate
      target_standard_error <- target_data$sample$treatment_effect_standard_error
      target_sample_size_per_arm <- target_data$sample_size_per_arm

      target_data_variance <- target_standard_error ^ 2 * target_sample_size_per_arm # sample variance (estimate of the true sampling variance)

      # calculate statistic and g(t)
      pooled_variance <- ((self$source_sample_size_per_arm - 1) * self$source_data_variance + (target_sample_size_per_arm - 1) * target_data_variance
      ) / (self$source_sample_size_per_arm + target_sample_size_per_arm - 2)

      congruence_measure <- max(self$source_sample_size_per_arm,
                                target_sample_size_per_arm) ^ (-1 / 4) * abs(self$prior$source$treatment_effect_estimate - target_treatment_effect_estimate) / (
                                  sqrt(
                                    pooled_variance / self$source_sample_size_per_arm + pooled_variance / target_sample_size_per_arm
                                  )
                                )

      assertions::assert_number(self$congruence_threshold)
      if (congruence_measure <= self$congruence_threshold) {
        elastic_function_value <- 1
      } else {
        elastic_function_value <- 0.00001 # 0
      }
      self$prior_var <- self$prior_var / elastic_function_value
    }
  )
)
