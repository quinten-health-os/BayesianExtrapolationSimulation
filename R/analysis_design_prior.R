#' @title DesignPrior Class
#' @description A class representing the design prior for Bayesian borrowing.
#' @field design_prior_type Type of design prior.
#' @field RBesT_model RBesT version of the model.
#' @export
DesignPrior <- R6::R6Class(
  "DesignPrior",
  public = list(
    design_prior_type = NULL,
    RBesT_model = NULL,

    #' @description Creates a new instance of DesignPrior.
    #' @param design_prior_type The type of design prior.
    #' @param model The model object.
    #' @param source_data The source data object.
    #' @param case_study_config Case study configuration
    #' @param simulation_config Simulation configuration
    #' @param mcmc_config MCMC configuration
    #' @param case_study Case study name
    #' @return A new instance of DesignPrior.
    create = function(design_prior_type,
                      model,
                      source_data,
                      case_study_config,
                      simulation_config,
                      mcmc_config,
                      case_study) {
      if (design_prior_type == "ui_design_prior") {
        design_prior <- UnitInformationDesignPrior$new(
          source_data = source_data,
          case_study_config = case_study_config,
          case_study = case_study,
          simulation_config = simulation_config,
          mcmc_config = mcmc_config
        )
      } else if (design_prior_type == "analysis_prior") {
        if (model$empirical_bayes) {
          stop("Analysis prior cannot be used as a design prior for methods that rely on empirical Bayes.")
        }

        design_prior <- AnalysisPriorDesignPrior$new(
          model = model,
          case_study_config = case_study_config,
          simulation_config = simulation_config,
          mcmc_config = mcmc_config
        )
      } else if (design_prior_type == "source_posterior") {
        design_prior <- SourcePosteriorDesignPrior$new(
          source_data = source_data,
          case_study_config = case_study_config,
          simulation_config = simulation_config,
          mcmc_config = mcmc_config
        )
      } else {
        stop("Not implemented for other types of design prior.")
      }
      return(design_prior)
    },

    #' @description Samples from the design prior.
    #' @param n_samples The number of samples to generate.
    sample = function(n_samples) {
      stop("The subclass must implement a sample() method")
    },

    #' @description Computes the cumulative distribution function (CDF) of the design prior.
    #' @param x The value at which to evaluate the CDF.
    cdf = function(x) {
      stop("The subclass must implement a cdf() method")
    },

    #' @description Computes the PDF of the design prior.
    #' @param x The value at which to evaluate the PDF.
    pdf = function(x) {
      stop("The subclass must implement a pdf() method")
    }
  )
)


binomial_product_integral <- function(x,
                                      alpha_control,
                                      beta_control,
                                      alpha_treatment,
                                      beta_treatment) {
  epsilon <- 1e-8
  lower_limit <- max(epsilon, -x + epsilon)
  upper_limit <- min(1 - epsilon, 1 - x - epsilon)

  # Handle cases where lower_limit > upper_limit
  if (lower_limit >= upper_limit) {
    return(0)
  }

  # Perform integration with error handling
  result <- tryCatch({
    integrate(
      function(y) {
        # Ensure y is within [0, 1] to avoid non-finite values
        valid_range <- (y >= 0) &
          (y <= 1) & ((y + x) >= 0) & ((y + x) <= 1)
        # Return 0 for out-of-bound values
        ifelse(
          valid_range,
          dbeta(y, alpha_control, beta_control) * dbeta(y + x, alpha_treatment, beta_treatment),
          0
        )
      },
      lower = lower_limit,
      upper = upper_limit,
      subdivisions = 1000,
      rel.tol = .Machine$double.eps ^ 0.25
    )$value
  }, error = function(e) {
    warning(paste("Integration failed at x =", x, ":", e$message))
    futile.logger::flog.error(paste("Integration failed at x =", x, ":", e$message))
    return(NA)
  })
  return(result)
}

#' @title UnitInformationDesignPrior class
#' @description A class representing the unit information design prior for Bayesian borrowing.
#' @field parameters Parameters for the prior
#' @field summary_measure_likelihood Treatment effect distribution
#' @export
UnitInformationDesignPrior <- R6::R6Class(
  "UnitInformationDesignPrior",
  inherit = DesignPrior,
  public = list(
    parameters = NULL,
    summary_measure_likelihood = NULL,

    #' @description Initializes a new instance of UnitInformationDesignPrior.
    #' @param source_data The source data object.
    #' @param case_study_config Case study configuration.
    #' @param case_study Case study name
    #' @param simulation_config Simulation config
    #' @param mcmc_config MCMC configuration
    initialize = function(source_data,
                          case_study_config,
                          case_study,
                          simulation_config,
                          mcmc_config = NULL) {
      summary_measure_likelihood <- source_data$summary_measure_likelihood
      self$summary_measure_likelihood <- summary_measure_likelihood

      if (summary_measure_likelihood == "normal") {
        # Ideally the variance should depend on the target study variance, but we don't know it, and assyle that it is the same as the source study variance
        variance <- (
          source_data$standard_error ^ 2 * source_data$equivalent_source_sample_size_per_arm
        )
        self$parameters <- list(
          mean = source_data$treatment_effect_estimate,
          variance = variance,
          sd = sqrt(variance)
        ) # information provided by a single subject per arm in the source study
      } else if (summary_measure_likelihood == "binomial") {
        # Start by fitting a separate analysis model on the source data
        method_parameters <- list(initial_prior = "noninformative", empirical_bayes = FALSE)

        if (is.null(mcmc_config)) {
          stop("MCMC config is not provided.")
        }

        separate_model <- Model$new()
        separate_model <- separate_model$create(
          case_study_config = case_study_config,
          method = "separate",
          method_parameters = method_parameters,
          source_data = source_data,
          mcmc_config = mcmc_config
        )


        # The rationale is the following:
        # 1. Define a target data object with the same properties as the source study
        # 2. Perform a separate analysis
        # 3. Approximate the resulting posterior distribution with a mixture distribution
        # 4. Scale the parameters of the distribution to get a Unit Information distribution
        target_data <- TargetDataFactory$new()
        target_data <- target_data$create(
          source_data = source_data,
          case_study_config = case_study_config,
          target_sample_size_per_arm = source_data$equivalent_source_sample_size_per_arm,
          control_drift = 0, # The actual value does matter here as we set the treatment effect estimate after
          treatment_drift = 0, # The actual value does matter here as we set the treatment effect estimate after
          summary_measure_likelihood = case_study_config$summary_measure_likelihood,
          target_to_source_std_ratio = 1
        )

        # The target data object used for the analysis has the same properties as the source data, but centered in theta0.
        target_data$sample_size_treatment <- source_data$sample_size_treatment
        target_data$sample_size_control <- source_data$sample_size_control
        target_data$sample$sample_treatment_rate <- source_data$treatment_rate
        target_data$sample$sample_control_rate <- source_data$control_rate
        target_data$sample$sample_size_per_arm <- source_data$equivalent_source_sample_size_per_arm
        target_data$sample$treatment_effect_estimate <- case_study_config$theta_0
        target_data$sample$treatment_effect_standard_error <- source_data$standard_error

        # Perform inference using a separate analysis of the source study
        separate_model$inference(target_data)

        #  Convert the posterior to a mixture distribution and compute the prior ESS
        separate_model$posterior_to_RBesT(target_data, simulation_config)

        mixture_approximation <-  separate_model$RBesT_posterior
        ESS <- RBesT::ess(mixture_approximation, method = "moment")

        # Extract the weights, a and b parameters
        weights <- mixture_approximation[1, ]
        a_params <- mixture_approximation[2, ]
        b_params <- mixture_approximation[3, ]

        # Adjust the a and b parameters by dividing by lambda
        a_params_new <- a_params / ESS
        b_params_new <- b_params / ESS


        # Create the new mixture of Beta distributions with adjusted parameters
        components_new <- lapply(1:length(weights), function(i) {
          c(weights[i], a_params_new[i], b_params_new[i])
        })

        # Construct the new mixture
        ui_mixture <- do.call(mixbeta, components_new)

        self$RBesT_model <- ui_mixture

        # Check that the resulting distribution has an ESS of 1
        #RBesT::ess(ui_mixture, method = "moment")
      } else {
        stop("Other distributions not supported")
      }

      self$design_prior_type <- "ui_design_prior"
    },

    #' @description Samples from the unit information design prior.
    #' @param n_samples The number of samples to generate.
    sample = function(n_samples) {
      if (self$summary_measure_likelihood == "normal") {
        return(rnorm(n_samples, self$parameters$mean, self$parameters$sd))
      } else if (self$summary_measure_likelihood == "binomial") {
        return(rmix(self$RBesT_model, n_samples))
      } else {
        stop("Not implemented for other distributions.")
      }
    },

    #' @description Computes the cumulative distribution function (CDF) of the unit information design prior.
    #' @param x The value at which to evaluate the CDF.
    cdf = function(x) {
      if (self$summary_measure_likelihood == "normal") {
        return(pnorm(x, self$parameters$mean, self$parameters$sd))
      } else if (self$summary_measure_likelihood == "binomial") {
        return(integrate(self$pdf, lower = -1, upper = x)$value)
      } else {
        stop("Not implemented for other distributions.")
      }
    },

    #' @description Computes the cumulative distribution function (CDF) of the unit information design prior.
    #' @param x The value at which to evaluate the CDF.
    pdf = function(x) {
      if (self$summary_measure_likelihood == "normal") {
        return(dnorm(x, self$parameters$mean, self$parameters$sd))
      } else if (self$summary_measure_likelihood == "binomial") {
        result <- sapply(x, function(x) {
          binomial_product_integral(
            x,
            alpha_control = 0.5,
            beta_control = 0.5,
            alpha_treatment = 0.5,
            beta_treatment = 0.5
          )
        })
        return(result)
      } else {
        stop("Not implemented for other distributions.")
      }
    }
  )
)

#' @title SourcePosteriorDesignPrior class
#' @description A class representing the source posterior design prior for Bayesian borrowing.
#' @field parameters Parameters
#' @field summary_measure_likelihood Treatment effect distribution
#' @field n_successes_control Number of successes in the control arm
#' @field n_successes_treatment Number of successes in the treatment arm
#' @field n_control Number of participants in the control arm
#' @field n_treatment Number of participants in the treatment arm
#' @export
SourcePosteriorDesignPrior <- R6::R6Class(
  "SourcePosteriorDesignPrior",
  inherit = DesignPrior,
  public = list(
    parameters = NULL,
    summary_measure_likelihood = NULL,
    n_successes_control = NULL,
    n_successes_treatment = NULL,
    n_control = NULL,
    n_treatment = NULL,

    #' @description Initializes a new instance of SourcePosteriorDesignPrior.
    #' @param source_data The source data object.
    #' @param case_study_config Case study configuration
    #' @param simulation_config Simulation configuration
    #' @param mcmc_config MCMC configuration
    initialize = function(source_data,
                          case_study_config,
                          simulation_config,
                          mcmc_config = NULL) {
      summary_measure_likelihood <- source_data$summary_measure_likelihood
      self$summary_measure_likelihood <- summary_measure_likelihood
      self$parameters <- list(
        mean = source_data$treatment_effect_estimate,
        variance = source_data$standard_error ^ 2,
        sd = source_data$standard_error
      )

      if (self$summary_measure_likelihood == "binomial") {
        self$n_successes_control <- int(source_data$control_rate * source_data$sample_size_control)
        self$n_successes_treatment <- int(source_data$treatment_rate * source_data$sample_size_treatment)
        self$n_control <- source_data$sample_size_control
        self$n_treatment <- source_data$sample_size_treatment
      }

      self$design_prior_type <- "source_posterior"
    },

    #' @description Samples from the source posterior design prior.
    #' @param n_samples The number of samples to generate.
    sample = function(n_samples) {
      if (self$summary_measure_likelihood == "normal") {
        return(rnorm(n_samples, self$parameters$mean, self$parameters$sd))
      } else if (self$summary_measure_likelihood == "binomial") {
        return(rmix(self$RBesT_model, n_samples))
      } else {
        stop("Not implemented for other distributions.")
      }
    },

    #' @description Computes the cumulative distribution function (CDF) of the source posterior design prior.
    #' @param x The value at which to evaluate the CDF.
    cdf = function(x) {
      if (self$summary_measure_likelihood == "normal") {
        return(pnorm(x, self$parameters$mean, self$parameters$sd))
      } else if (self$summary_measure_likelihood == "binomial") {
        return(integrate(self$pdf, lower = -1, upper = x)$value)
      } else {
        stop("Not implemented for other distributions.")
      }
    },

    #' @description Computes the cumulative distribution function (PDF) of the source posterior design prior.
    #' @param x The value at which to evaluate the PDF.
    pdf = function(x) {
      if (self$summary_measure_likelihood == "normal") {
        return(dnorm(x, self$parameters$mean, self$parameters$sd))
      } else if (self$summary_measure_likelihood == "binomial") {
        # PDF for each arm
        result <- sapply(x, function(x) {
          binomial_product_integral(
            x,
            alpha_control = 1 + self$n_successes_control,
            beta_control = 1 + self$n_control - self$n_successes_control,
            alpha_treatment = 1 + self$n_successes_treatment,
            beta_treatment = 1 + self$n_treatment - self$n_successes_treatment
          )
        })
      } else {
        stop("Not implemented for other distributions.")
      }
    }
  )
)

#' @title AnalysisPriorDesignPrior class
#' @description A class representing the analysis prior design prior for Bayesian borrowing.
#' @field model Model used to define the analysis prior
#' @export
AnalysisPriorDesignPrior <- R6::R6Class(
  "AnalysisPriorDesignPrior",
  inherit = DesignPrior,
  public = list(
    model = NULL,

    #' @description Initializes a new instance of AnalysisPriorDesignPrior.
    #' @param model The model object.
    #' @param case_study_config Case study configuration
    #' @param simulation_config Simulation configuration
    #' @param mcmc_config MCMC configuration
    initialize = function(model,
                          case_study_config,
                          simulation_config,
                          mcmc_config) {
      if (model$empirical_bayes == TRUE) {
        stop(
          "It is not possible to define an analysis design prior for a method that uses empirical Bayes."
        )
      }

      model$prior_to_RBesT(simulation_config$n_samples_mixture_approx)
      self$model <- model
      self$design_prior_type <- "analysis_prior"
    },

    #' @description Samples from the analysis prior design prior.
    #' @param n_samples The number of samples to generate.
    sample = function(n_samples) {
      return(self$model$sample_prior(n_samples))
    },

    #' @description Computes the cumulative distribution function (CDF) of the analysis prior design prior.
    #' @param x The value at which to evaluate the CDF.
    cdf = function(x) {
      return(self$model$prior_cdf(x))
    },

    #' @description Computes the PDF of the analysis prior design prior.
    #' @param x The value at which to evaluate the PDF.
    pdf = function(x) {
      return(self$model$prior_pdf(x))
    }
  )
)
