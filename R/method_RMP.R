#' Calculate the variance of a mixture distribution
#'
#' This function computes the variance of a mixture of two distributions.
#'
#' @param w Numeric. The weight of the first distribution in the mixture.
#' @param sigma2A Numeric. The variance of the first distribution.
#' @param muA Numeric. The mean of the first distribution.
#' @param sigma2B Numeric. The variance of the second distribution.
#' @param muB Numeric. The mean of the second distribution.
#'
#' @return Numeric. The total variance of the mixture distribution.
#'
#' @keywords internal
mixture_variance <- function(w, sigma2A, muA, sigma2B, muB) {
  # Weighted sum of the variances
  weighted_var_sum <- w * sigma2A + (1 - w) * sigma2B

  # Weighted sum of the square of the means
  weighted_mean_sq_sum <- w * muA ^ 2 + (1 - w) * muB ^ 2

  # Square of the weighted sum of the means
  mean_weighted_sum_sq <- (w * muA + (1 - w) * muB) ^ 2

  # Total variance of the mixture
  variance <- weighted_var_sum + weighted_mean_sq_sum - mean_weighted_sum_sq

  return(variance)
}


#' GaussianRMP_RBesT class
#'
#' This class represents a Gaussian Robust Mixture Prior model.
#' It inherits from the Model class. Inference is performed using RBesT.
#'
#' @description
#' A class for Gaussian Robust Mixture Prior models using RBesT for inference.
#'
#' @details
#' This class extends the Model_RBesT class and provides methods for
#' initializing the model, updating priors, calculating posterior moments,
#' and converting distributions to RBesT format.
#'
#' @field w The weight of the prior distribution.
#' @field vague_prior_mean The mean of the vague prior distribution.
#' @field vague_prior_variance The variance of the vague prior distribution.
#' @field info_prior_mean The mean of the informative prior distribution.
#' @field info_prior_variance The variance of the informative prior distribution.
#' @field wpost The weight of the posterior distribution.
#' @field vague_posterior_mean The mean of the vague posterior distribution.
#' @field info_posterior_mean The mean of the informative posterior distribution.
#' @field vague_posterior_variance The variance of the vague posterior distribution.
#' @field info_posterior_variance The variance of the informative posterior distribution.
#' @field method Method name
#'
#' @export
GaussianRMP_RBesT <- R6::R6Class(
  "GaussianRMP_RBesT",
  inherit = Model_RBesT,
  public = list(
    w = NULL,
    vague_prior_mean = NULL,
    vague_prior_variance = NULL,
    info_prior_mean = NULL,
    info_prior_variance = NULL,
    wpost = NULL,
    vague_posterior_mean = NULL,
    info_posterior_mean = NULL,
    vague_posterior_variance = NULL,
    info_posterior_variance = NULL,
    method = "RMP",

    #' @description
    #' Initialize a new GaussianRMP_RBesT object.
    #' @param prior A list containing prior information for the analysis.
    initialize = function(prior) {
      if (!(prior$method_parameters$initial_prior[[1]] == "noninformative")) {
        stop("Only implemented for a noninformative initial prior")
      }
      super$initialize(prior)
      self$w <- unlist(prior$method_parameters$prior_weight)
      self$vague_prior_mean <- prior$vague_mean

      self$info_prior_mean <- prior$source$treatment_effect_estimate
      self$info_prior_variance <- prior$source$standard_error ^ 2

      if (prior$method_parameters$empirical_bayes[[1]]) {
        self$empirical_bayes <- TRUE
        self$vague_prior_variance <- NULL
      } else {
        self$empirical_bayes <- FALSE
        self$vague_prior_variance <- prior$source$standard_error ^ 2 *
          prior$source$equivalent_source_sample_size_per_arm # information provided by a single subject per arm in the source study, corresponds to the approach used in Best et al, 2023 (Belimumab)

        self$prior_to_RBesT()

        prior_summary <- summary(self$RBesT_prior)

        self$prior_mean <- prior_summary['mean']
        self$prior_var <- prior_summary['sd']^2
      }

      self$wpost <- NULL
      self$vague_posterior_mean <- NULL
      self$vague_posterior_variance <- NULL
      self$info_posterior_mean <- NULL
      self$info_posterior_variance <- NULL

      self$posterior_parameters <- list(prior_weight = NA)
    },

    #' @description
    #' Update the vague prior variance based on empirical Bayes approach.
    #' @param target_data A list containing the target data for the analysis.
    empirical_bayes_update = function(target_data) {
      if (self$empirical_bayes) {
        self$vague_prior_variance <- (target_data$sample$treatment_effect_standard_error ^
                                        2) * target_data$sample_size_per_arm # sample standard deviation in the target study

        # information provided by a single subject per arm in the target study (as explained in the protocol)
        # Note that in Best et al, 2021 (Mepolizumab), this is the information provided by a single subject in the target study

        info = c(
          self$w,
          self$info_prior_mean,
          sqrt(self$info_prior_variance)
        )

        vague = c(
          1 - self$w,
          self$vague_prior_mean,
          sqrt(self$vague_prior_variance)
        )

        # We use the following conditions to avoid bug in RBesT ELIR computation if self$wpost == 0
        if (self$w == 1){
          self$RBesT_prior <- RBesT::mixnorm(
            info = info
          )
        } else if (self$w == 0){
          self$RBesT_prior <- RBesT::mixnorm(
            vague = vague
          )
        } else {
          self$RBesT_prior <- RBesT::mixnorm(
            info = info,
            vague = vague
          )
        }

        prior_summary <- summary(self$RBesT_prior)
        self$prior_mean <- prior_summary['mean']
        self$prior_var <- prior_summary['sd']^2

        RBesT::sigma(self$RBesT_prior) <- target_data$sample$standard_deviation

        self$RBesT_prior_normix <- self$RBesT_prior
      } else {
        if (is.null(self$vague_prior_variance)) {
          stop("The vague prior variance must be defined if empirical_bayes is False") # information provided by a single subject per arm in the source study : Best et al (2023, case study on Belimumab)
        }
      }
    },

    #' @description
    #' Calculate the posterior moments based on the target data.
    #' @param target_data A list containing the target data for the analysis.
    posterior_moments = function(target_data) {
      self$RBesT_posterior <- RBesT::postmix(
        self$RBesT_prior,
        m = target_data$sample$treatment_effect_estimate,
        se = target_data$sample$treatment_effect_standard_error
      )

      inference_results <- summary(self$RBesT_posterior)
      names(inference_results) <- c("mean", "standard_deviation", "cri95L", "median", "cri95U")
      self$posterior_summary <- inference_results

      # The following if statement is used because of a bug in the computation of ELIR with RBesT that requires to define a single component of the RMP if w = 0 or w = 1
      if (self$w == 0){
        self$wpost <- 0
        self$vague_posterior_mean <-  self$RBesT_posterior[2]
        self$vague_posterior_variance <- self$RBesT_posterior[3]^ 2
      } else if (self$w == 1){
        self$wpost <- 1
        self$info_posterior_mean <-  self$RBesT_posterior[2]
        self$info_posterior_variance <- self$RBesT_posterior[3]^ 2
      } else {
        self$wpost <- self$RBesT_posterior[1]
        self$info_posterior_mean <- self$RBesT_posterior[2]
        self$info_posterior_variance <- self$RBesT_posterior[3] ^ 2
        self$vague_posterior_mean <- self$RBesT_posterior[5]
        self$vague_posterior_variance <- self$RBesT_posterior[6]^ 2
      }

      self$posterior_parameters$prior_weight <- self$wpost

      self$post_mean <- self$posterior_summary["mean"]

      self$post_var <- unname(self$posterior_summary["standard_deviation"] ^ 2)
      self$post_median <- self$posterior_summary["median"]
      assertions::assert_number(self$post_mean)
    },


    #' @description
    #' Convert the prior distribution to the RBesT format.
    #' @param ... Additional arguments.
    prior_to_RBesT = function(...) {
      if (is.null(self$vague_prior_variance)) {
        stop(
          "The vague prior variance must be defined. If empirical_bayes == TRUE, the prior can only be fully specified after the data has been observed."
        )
      }

      info = c(
        self$w,
        self$info_prior_mean,
        sqrt(self$info_prior_variance)
      )

      vague = c(
        1 - self$w,
        self$vague_prior_mean,
        sqrt(self$vague_prior_variance)
      )

      # We use the following conditions to avoid bug in RBesT ELIR computation if self$wpost == 0
      if (self$w == 1){
        self$RBesT_prior <- RBesT::mixnorm(
          info = info
        )
      } else if (self$w == 0){
        self$RBesT_prior <- RBesT::mixnorm(
          vague = vague
        )
      } else {
        self$RBesT_prior <- RBesT::mixnorm(
          info = info,
          vague = vague
        )
      }

      self$RBesT_prior_normix <- self$RBesT_prior
      return(self$RBesT_prior)
    },

    #' @description
    #' Convert the posterior distribution to the RBesT format.
    #' @param target_data Target data for the analysis.
    #' @param ... Additional arguments.
    posterior_to_RBesT = function(target_data, ...) {
      RBesT::sigma(self$RBesT_posterior) = target_data$sample$standard_deviation

      self$RBesT_posterior_normix <- self$RBesT_posterior

      return(self$RBesT_posterior)
    },


    #' @description
    #' Print a summary of the model attributes
    #'
    #' This method creates and prints a formatted table of key model attributes.
    #'
    #' @return A printed data frame displaying the following model attributes:
    #'   \itemize{
    #'     \item Posterior Weight
    #'     \item Vague Posterior Mean
    #'     \item Vague Posterior Variance
    #'     \item Informative Posterior Mean
    #'     \item Informative Posterior Variance
    #'   }
    #'   All numeric values are formatted to 6 decimal places.
    print_model_summary = function() {
      # Create a data frame with the model attributes
      df <- data.frame(
        Attribute = c(
          "Posterior Weight",
          "Vague Posterior Mean",
          "Vague Posterior Variance",
          "Informative Posterior Mean",
          "Informative Posterior Variance"
        ),
        Value = c(
          self$wpost,
          self$vague_posterior_mean,
          self$vague_posterior_variance,
          self$info_posterior_mean,
          self$info_posterior_variance
        )
      )

      # Format the values to 6 decimal places
      df$Value <- format(df$Value, digits = 6, nsmall = 6)

      # Print the table
      print(df, row.names = FALSE)
    }
  )
)


#' GaussianRMP class
#'
#' This class represents a Gaussian Robust Mixture Prior model.
#' It inherits from the Model class.
#'
#' @description
#' A class for Gaussian Robust Mixture Prior models.
#'
#' @details
#' This class extends the Model class and provides methods for
#' initializing the model, updating priors, calculating posterior moments,
#' and sampling from prior and posterior distributions.
#'
#' @field w The weight of the prior distribution.
#' @field vague_prior_mean The mean of the vague prior distribution.
#' @field vague_prior_variance The variance of the vague prior distribution.
#' @field info_prior_mean The mean of the informative prior distribution.
#' @field info_prior_variance The variance of the informative prior distribution.
#' @field wpost The weight of the posterior distribution.
#' @field vague_posterior_mean The mean of the vague posterior distribution.
#' @field info_posterior_mean The mean of the informative posterior distribution.
#' @field vague_posterior_variance The variance of the vague posterior distribution.
#' @field info_posterior_variance The variance of the informative posterior distribution.
#' @field method Name of the method
#'
#' @keywords internal
GaussianRMP <- R6::R6Class(
  "GaussianRMP",
  inherit = Model,
  public = list(
    w = NULL,
    vague_prior_mean = NULL,
    vague_prior_variance = NULL,
    info_prior_mean = NULL,
    info_prior_variance = NULL,
    wpost = NULL,
    vague_posterior_mean = NULL,
    info_posterior_mean = NULL,
    vague_posterior_variance = NULL,
    info_posterior_variance = NULL,
    method = "RMP",

    #' @description
    #' Initialize a new GaussianRMP object.
    #' @param prior A list containing prior information for the analysis.
    initialize = function(prior) {
      if (!(prior$method_parameters$initial_prior[[1]] == "noninformative")) {
        stop("Only implemented for a noninformative initial prior")
      }
      super$initialize()
      self$w <- unlist(prior$method_parameters$prior_weight)

      if (prior$method_parameters$empirical_bayes[[1]]) {
        self$empirical_bayes <- TRUE
        self$vague_prior_variance <- NULL
      } else {
        self$empirical_bayes <- FALSE
        self$vague_prior_variance <- prior$source$standard_error ^ 2 *
          prior$source$equivalent_source_sample_size_per_arm # information provided by a single subject per arm in the source study, corresponds to the approach used in Best et al, 2023 (Belimumab)
      }
      self$vague_prior_mean <- prior$vague_mean

      self$info_prior_mean <- prior$source$treatment_effect_estimate
      self$info_prior_variance <- prior$source$standard_error ^ 2


      self$wpost <- NULL
      self$vague_posterior_mean <- NULL
      self$vague_posterior_variance <- NULL
      self$info_posterior_mean <- NULL
      self$info_posterior_variance <- NULL

      self$posterior_parameters <- list(prior_weight = NA)
    },

    #' @description
    #' Update the vague prior variance based on empirical Bayes approach.
    #' @param target_data The target data for the analysis.
    empirical_bayes_update = function(target_data) {
      if (self$empirical_bayes) {
        self$vague_prior_variance <- target_data$sample$treatment_effect_standard_error ^
          2 * target_data$sample_size_per_arm # sample standard deviation in the target study

        # information provided by a single subject per arm in the target study (as explained in the protocol)
        # Note that in Best et al, 2021 (Mepolizumab), this is the information provided by a single subject in the target study
      } else {
        if (is.null(self$vague_prior_variance)) {
          stop("The vague prior variance must be defined if empirical_bayes is False") # information provided by a single subject per arm in the source study : Best et al (2023, case study on Belimumab)
        }
      }
    },

    #' @description
    #' Calculate the posterior weight based on the target data.
    #' @param target_data The target data for the analysis.
    #' @return The posterior weight.
    prior_weight = function(target_data) {
      sample_mean <- target_data$sample$treatment_effect_estimate
      sample_size_per_arm <- target_data$sample_size_per_arm
      sample_standard_error <- target_data$sample$treatment_effect_standard_error

      if (is.null(sample_mean) || is.null(sample_size_per_arm)) {
        stop("One or more of the required variables (sample_mean, sample_size_per_arm) is NULL.")
      }



      log_likelihood_vague <- (
        -0.5
        * (sample_mean - self$vague_prior_mean) ^ 2
        / (self$vague_prior_variance + sample_standard_error ^ 2)
      ) - 0.5 * log(sample_standard_error ^ 2 + self$vague_prior_variance)

      log_likelihood_source <- (
        -0.5
        * (sample_mean - self$info_prior_mean) ^ 2
        / (self$info_prior_variance + sample_standard_error ^ 2)
      ) - 0.5 * log(
        target_data$sample$treatment_effect_standard_error ^ 2 + self$info_prior_variance
      )

      log_unnormalized_post_weights <- c(log(self$w) + log_likelihood_source,
                                         log(1 - self$w) + log_likelihood_vague)

      post_weights <- exp(
        log_unnormalized_post_weights - matrixStats::logSumExp(log_unnormalized_post_weights)
      )

      wpost <- post_weights[1]

      assertions::assert_number(wpost)
      if (is.na(wpost)) {
        stop(
          "Posterior w is NaN, most likely due to both the vague and informative likelihood equal to 0."
        )
      }
      self$wpost <- wpost
      self$posterior_parameters$prior_weight <- self$wpost
      return(self$wpost)
    },

    #' @description  Calculate the prior probability density function (PDF) for a given target treatment effect.
    #' @param target_treatment_effect The target treatment effect.
    #' @return The prior PDF.
    prior_pdf = function(target_treatment_effect) {
      prior_dist <- (1 - self$w) * dnorm(
        target_treatment_effect,
        mean = self$vague_prior_mean,
        sd = sqrt(self$vague_prior_variance)
      ) + self$w * dnorm(
        target_treatment_effect,
        mean = self$info_prior_mean,
        sd = sqrt(self$info_prior_variance)
      )
      return(prior_dist)
    },

    #' @description
    #' Calculate the prior cumulative distribution function (CDF) for a given target treatment effect.
    #' @param target_treatment_effect The target treatment effect.
    #' @return The prior CDF.
    prior_cdf = function(target_treatment_effect) {
      prior_proba <- (1 - self$w) * pnorm(
        target_treatment_effect,
        mean = self$vague_prior_mean,
        sd = sqrt(self$vague_prior_variance)
      ) + self$w * pnorm(
        target_treatment_effect,
        mean = self$info_prior_mean,
        sd = sqrt(self$info_prior_variance)
      )
      return(prior_proba)
    },

    #' @description
    #' Calculate the posterior moments based on the target data.
    #' @param target_data A list containing the target data for the analysis.
    posterior_moments = function(target_data) {
      sample_mean <- target_data$sample$treatment_effect_estimate
      standard_error <- target_data$sample$treatment_effect_standard_error

      self$prior_weight(target_data)


      vague_posterior_mean <- self$vague_prior_mean / (self$vague_prior_variance / standard_error ^
                                                         2 + 1) + sample_mean / (1
                                                                                 + standard_error ^ 2
                                                                                 / self$vague_prior_variance)
      vague_posterior_variance <- 1 / (1 / self$vague_prior_variance + 1 / standard_error ^
                                         2)

      info_posterior_mean <- self$info_prior_mean / (self$info_prior_variance / standard_error ^
                                                       2 + 1) + sample_mean / (1 + standard_error ^ 2 / self$info_prior_variance)
      info_posterior_variance <- 1 / (1 / self$info_prior_variance + 1 / standard_error ^
                                        2)

      self$vague_posterior_mean <- vague_posterior_mean
      self$vague_posterior_variance <- vague_posterior_variance
      self$info_posterior_mean <- info_posterior_mean
      self$info_posterior_variance <- info_posterior_variance

      self$post_mean <- self$posterior_mean()
      self$post_var <- self$posterior_variance()

    },

    #' @description
    #' Calculate the posterior probability density function (PDF) for a given target treatment effect.
    #' @param target_treatment_effect The target treatment effect.
    #' @return The posterior PDF.
    posterior_pdf = function(target_treatment_effect) {
      posterior_dist <- (1 - self$wpost) * dnorm(
        target_treatment_effect,
        mean = self$vague_posterior_mean,
        sd = sqrt(self$vague_posterior_variance)
      ) + self$wpost * dnorm(
        target_treatment_effect,
        mean = self$info_posterior_mean,
        sd = sqrt(self$info_posterior_variance)
      )
      return(posterior_dist)
    },

    #' @description
    #' Calculate the posterior cumulative distribution function (CDF) for a given target treatment effect.
    #' @param target_treatment_effect The target treatment effect.
    #' @return The posterior CDF.
    posterior_cdf = function(target_treatment_effect) {
      posterior_proba <- (1 - self$wpost) * pnorm(
        target_treatment_effect,
        mean = self$vague_posterior_mean,
        sd = sqrt(self$vague_posterior_variance)
      ) + self$wpost * pnorm(
        target_treatment_effect,
        mean = self$info_posterior_mean,
        sd = sqrt(self$info_posterior_variance)
      )
      return(posterior_proba)
    },

    #' @description
    #' Sample from the prior distribution.
    #' @param n_samples The number of samples to generate.
    #' @return The samples from the prior distribution.
    sample_prior = function(n_samples) {
      component_choice <- sample(
        x = c(0, 1),
        size = n_samples,
        replace = TRUE,
        prob = c(1 - self$w, self$w)
      )

      samples <- rep(0, n_samples)

      samples[component_choice == 0] <- rnorm(
        n = sum(component_choice == 0),
        mean = self$vague_prior_mean,
        sd = sqrt(self$vague_prior_variance)
      )

      samples[component_choice == 1] <- rnorm(
        n = sum(component_choice == 1),
        mean = self$info_prior_mean,
        sd = sqrt(self$info_prior_variance)
      )

      return(samples)
    },

    #' @description
    #' Sample from the posterior distribution.
    #' @param n_samples The number of samples to generate.
    #' @return The samples from the posterior distribution.
    sample_posterior = function(n_samples) {
      component_choice <- sample(
        x = c(0, 1),
        size = n_samples,
        replace = TRUE,
        prob = c(1 - self$wpost, self$wpost)
      )

      samples <- rep(0, n_samples)

      samples[component_choice == 0] <- rnorm(
        n = sum(component_choice == 0),
        mean = self$vague_posterior_mean,
        sd = sqrt(self$vague_posterior_variance)
      )

      samples[component_choice == 1] <- rnorm(
        n = sum(component_choice == 1),
        mean = self$info_posterior_mean,
        sd = sqrt(self$info_posterior_variance)
      )

      return(samples)
    },

    #' @description
    #' Calculate the posterior mean.
    #' @return The posterior mean.
    posterior_mean = function() {
      return((1 - self$wpost) * self$vague_posterior_mean + self$wpost * self$info_posterior_mean)
    },

    #' @description
    #' Calculate the posterior variance.
    #' @return The posterior variance.
    posterior_variance = function() {
      post_var <- mixture_variance(
        self$wpost,
        self$info_posterior_variance,
        self$info_posterior_mean,
        self$vague_posterior_variance,
        self$vague_posterior_mean
      )
      return(post_var)
    },

    #' @description
    #' Convert the prior distribution to the RBesT format.
    #' @param ... Additional arguments.
    prior_to_RBesT = function(...) {
      if (is.null(self$vague_prior_variance)) {
        stop(
          "The vague prior variance must be defined. If empirical_bayes == TRUE, the prior can only be fully specified after the data has been observed."
        )
      }

      info = c(
        self$w,
        self$info_prior_mean,
        sqrt(self$info_prior_variance)
      )

      vague = c(
        1 - self$w,
        self$vague_prior_mean,
        sqrt(self$vague_prior_variance)
      )

      # We use the following conditions to avoid bug in RBesT ELIR computation if self$wpost == 0
      if (self$w == 1){
        self$RBesT_prior <- RBesT::mixnorm(
          info = info
        )
      } else if (self$w == 0){
        self$RBesT_prior <- RBesT::mixnorm(
          vague = vague
        )
      } else {
        self$RBesT_prior <- RBesT::mixnorm(
          info = info,
          vague = vague
        )
      }

      self$RBesT_prior_normix <- self$RBesT_prior

      return(self$RBesT_prior)
    },

    #' @description
    #' Convert the posterior distribution to the RBesT format.
    #' @param target_data Target data for the analysis.
    #' @param ... Additional arguments.
    posterior_to_RBesT = function(target_data, ...) {
      info = c(
        self$wpost,
        self$info_posterior_mean,
        sqrt(self$info_posterior_variance)
      )

      vague = c(
        1 - self$wpost,
        self$vague_posterior_mean,
        sqrt(self$vague_posterior_variance)
      )

      # We use the following conditions to avoid bug in RBesT ELIR computation if self$wpost == 0
      if (self$wpost == 0){
        self$RBesT_posterior <- RBesT::mixnorm(
          vague = vague
        )
      } else if (self$wpost == 1){
        self$RBesT_posterior <- RBesT::mixnorm(
          info = info
        )
      } else {
        self$RBesT_posterior <- RBesT::mixnorm(
          info = info,
          vague = vague
        )
      }

      RBesT::sigma(self$RBesT_posterior) = target_data$sample$standard_deviation

      self$RBesT_posterior_normix <- self$RBesT_posterior

      return(self$RBesT_posterior)
    },

    #' @description
    #' Print a summary of the model attributes
    #'
    #' This method creates and prints a formatted table of key model attributes.
    #'
    #' @return A printed data frame displaying the following model attributes:
    #'   \itemize{
    #'     \item Posterior Weight
    #'     \item Vague Posterior Mean
    #'     \item Vague Posterior Variance
    #'     \item Informative Posterior Mean
    #'     \item Informative Posterior Variance
    #'   }
    #'   All numeric values are formatted to 6 decimal places.
    print_model_summary = function() {
      # Create a data frame with the model attributes
      df <- data.frame(
        Attribute = c(
          "Posterior Weight",
          "Vague Posterior Mean",
          "Vague Posterior Variance",
          "Informative Posterior Mean",
          "Informative Posterior Variance"
        ),
        Value = c(
          self$wpost,
          self$vague_posterior_mean,
          self$vague_posterior_variance,
          self$info_posterior_mean,
          self$info_posterior_variance
        )
      )

      # Format the values to 6 decimal places
      df$Value <- format(df$Value, digits = 6, nsmall = 6)

      # Print the table
      print(df, row.names = FALSE)
    }
  )
)

#' TruncatedGaussianRMP class
#'
#' @description This class represents a Bayesian borrowing model with a truncated Gaussian prior. It inherits from the MCMCModel class.
#'
#' @field w The weight parameter for the mixture prior
#' @field vague_prior_mean The mean of the vague prior distribution
#' @field vague_posterior_mean The mean of the vague posterior distribution
#' @field vague_prior_variance The variance of the vague prior distribution
#' @field vague_posterior_variance The variance of the vague posterior distribution
#' @field info_prior_mean The mean of the informative prior distribution
#' @field info_posterior_mean The mean of the informative posterior distribution
#' @field info_prior_variance The variance of the informative prior distribution
#' @field info_posterior_variance The variance of the informative posterior distribution
#' @field wpost The posterior weight parameter for the mixture prior
#' @field method Method name
#'
#' @return An R6 class object representing a TruncatedGaussianRMP model
#' @export
TruncatedGaussianRMP <- R6::R6Class(
  "TruncatedGaussianRMP",
  inherit = MCMCModel,
  public = list(
    w = NULL,
    vague_prior_mean = NULL,
    vague_posterior_mean = NULL,
    vague_prior_variance = NULL,
    vague_posterior_variance = NULL,
    info_prior_mean = NULL,
    info_posterior_mean = NULL,
    info_prior_variance = NULL,
    info_posterior_variance = NULL,
    wpost = NULL,
    method = "RMP",

    #' @description Initialize the TruncatedGaussianRMP object
    #'
    #' @param prior The prior object
    #' @param mcmc_config The MCMC configuration parameters
    #'
    initialize = function(prior, mcmc_config) {
      if (!(prior$method_parameters$initial_prior[[1]] == "noninformative")) {
        stop("Only implemented for a noninformative initial prior")
      }
      super$initialize(prior = prior, mcmc_config = mcmc_config)

      stan_model_code <- "
        data {
      int<lower = 0> n_treatment;
      int<lower = 0> n_control;

      int<lower = 0, upper = n_treatment> n_successes_treatment;
      int<lower = 0, upper = n_control> n_successes_control;

      real<lower = 0, upper = 1> w;

      real<lower = -1, upper = 1> vague_mean;
      real<lower = 0> vague_sd;
      real<lower = -1, upper = 1> info_mean;
      real<lower = 0> info_sd;
    }
    transformed data {
      int<lower = 0, upper = 1> debug = 0;
    }
    parameters {
      real<lower = 0, upper = 1> control_rate;
      real<lower = - control_rate, upper = 1 - control_rate> target_treatment_effect;
    }
    transformed parameters {
      real treatment_rate = control_rate + target_treatment_effect;

      real vague_normalizing_constant = normal_cdf(1-control_rate | vague_mean, vague_sd) - normal_cdf(-control_rate | vague_mean, vague_sd);
      real info_normalizing_constant = normal_cdf(1-control_rate | info_mean, info_sd) - normal_cdf(-control_rate | info_mean, info_sd);
      real log_vague_normalizing_constant = log(vague_normalizing_constant);
      real log_info_normalizing_constant = log(info_normalizing_constant);
    }
    model {
      // robust mixture prior

      control_rate ~ uniform(0, 1);

      target += log_mix(w,
                normal_lpdf(target_treatment_effect | vague_mean, vague_sd) - log_vague_normalizing_constant,
                normal_lpdf(target_treatment_effect | info_mean, info_sd) - log_info_normalizing_constant);

      n_successes_control ~ binomial(n_control, control_rate);
      n_successes_treatment ~ binomial(n_treatment, treatment_rate);
    }
    "

      model_name <- "truncated_gaussian_rmp"
      self$stan_model <- compile_stan_model(model_name, stan_model_code)

      self$w <- unlist(prior$method_parameters$prior_weight)
      self$vague_prior_mean <- prior$vague_mean

      self$mcmc <- TRUE

      self$mcmc_config <- mcmc_config

      if (prior$method_parameters$empirical_bayes[[1]]) {
        self$empirical_bayes <- TRUE
        self$vague_prior_variance <- NULL
      } else {
        self$empirical_bayes <- FALSE

        self$vague_prior_variance <- prior$source$standard_error ^ 2 *
          prior$source$equivalent_source_sample_size_per_arm # information provided by a single subject per arm in the source study, corresponds to the approach used in Best et al, 2023 (Belimumab)
      }

      self$info_prior_mean <- prior$source$treatment_effect_estimate
      self$info_prior_variance <- prior$source$standard_error ^ 2

      self$vague_posterior_mean <- NULL
      self$vague_posterior_variance <- NULL
      self$info_posterior_mean <- NULL
      self$info_posterior_variance <- NULL

      self$method <- "RMP"
      self$posterior_parameters <- list(prior_weight = NA)

      self$summary_measure_likelihood <- "binomial"
    },

    #' @description Update the empirical Bayes parameters
    #'
    #' @param target_data The target data for inference
    #'
    empirical_bayes_update = function(target_data) {
      if (self$empirical_bayes) {
        self$vague_prior_variance <- target_data$sample$treatment_effect_standard_error ^
          2 * target_data$sample_size_per_arm # sample standard deviation in the target study

        # information provided by a single subject per arm in the target study (as explained in the protocol)
        # Note that in Best et al, 2021 (Mepolizumab), this is the information provided by a single subject in the target study
      } else {
        if (is.null(self$vague_prior_variance)) {
          stop("The vague prior variance must be defined if empirical_bayes is False.") # information provided by a single subject per arm in the source study : Best et al (2023, case study on Belimumab)
        }
      }
    },
    #' @description Prepare the data for use with Stan
    #'
    #' @param target_data The target data for inference
    #'
    #' @return A list of data prepared for Stan
    prepare_data = function(target_data) {
      data_list <- list(
        w = self$w,
        n_treatment = as.integer(target_data$sample_size_treatment),
        n_control = as.integer(target_data$sample_size_control),
        n_successes_treatment = as.integer(
          target_data$sample_size_treatment * target_data$sample$sample_treatment_rate
        ),
        n_successes_control = as.integer(
          target_data$sample_size_control * target_data$sample$sample_control_rate
        ),
        vague_mean = self$vague_prior_mean,
        vague_sd = sqrt(self$vague_prior_variance),
        info_mean = self$info_prior_mean,
        info_sd = sqrt(self$info_prior_variance)
      )
      return(data_list)
    },
    #' @description Sample from the prior distribution
    #'
    #' @param n_samples Number of samples to generate
    #' @return A vector of samples from the prior distribution
    sample_prior = function(n_samples) {

      assert_number(n_samples)
      component_choice <- sample(
        x = c(0, 1),
        size = n_samples,
        replace = TRUE,
        prob = c(1 - self$w, self$w)
      )

      samples <- rep(0, n_samples)

      lower <- -1 # Lower truncation point
      upper <- 1 # Upper truncation point

      samples[component_choice == 0] <- truncnorm::rtruncnorm(
        n = sum(component_choice == 0),
        mean = self$vague_prior_mean,
        sd = sqrt(self$vague_prior_variance),
        a = lower,
        b = upper
      )

      samples[component_choice == 1] <- truncnorm::rtruncnorm(
        n = sum(component_choice == 1),
        mean = self$info_prior_mean,
        sd = sqrt(self$info_prior_variance),
        a = lower,
        b = upper
      )
      return(samples)
    },

    #' @description Inference
    #'
    #' @param target_data Target data object
    #' @return Indicator whether inference succeeded or not
    inference = function(target_data) {
      fit_success <- super$inference(target_data)

      lower <- -1 # Lower truncation point
      upper <- 1 # Upper truncation point

      prior_predictive_proba_vague <- truncnorm::dtruncnorm(
        x = target_data$sample$treatment_effect_estimate,
        a = lower,
        b = upper,
        mean = self$vague_prior_mean,
        sd = sqrt(self$vague_prior_variance)
      )

      prior_predictive_proba_info <- truncnorm::dtruncnorm(
        x = target_data$sample$treatment_effect_estimate,
        a = lower,
        b = upper,
        mean = self$info_prior_mean,
        sd = sqrt(self$info_prior_variance)
      )

      self$posterior_parameters$prior_weight <- self$w * prior_predictive_proba_vague / (self$w * prior_predictive_proba_vague + (1 - self$w) * prior_predictive_proba_info)

      return(fit_success)
    },

    #' @description
    #' Print a summary of the model attributes
    #'
    #' This method creates and prints a formatted table of key model attributes.
    #'
    #' @return A printed data frame displaying the following model attributes:
    #'   \itemize{
    #'     \item Posterior Weight
    #'     \item Vague Posterior Mean
    #'     \item Vague Posterior Variance
    #'     \item Informative Posterior Mean
    #'     \item Informative Posterior Variance
    #'   }
    #'   All numeric values are formatted to 6 decimal places.
    print_model_summary = function() {
      # Create a data frame with the model attributes
      df <- data.frame(
        Attribute = c(
          "Posterior Weight",
          "Vague Posterior Mean",
          "Vague Posterior Variance",
          "Informative Posterior Mean",
          "Informative Posterior Variance"
        ),
        Value = c(
          self$wpost,
          self$vague_posterior_mean,
          self$vague_posterior_variance,
          self$info_posterior_mean,
          self$info_posterior_variance
        )
      )

      # Format the values to 6 decimal places
      df$Value <- format(df$Value, digits = 6, nsmall = 6)

      # Print the table
      print(df, row.names = FALSE)
    }
  )
)
