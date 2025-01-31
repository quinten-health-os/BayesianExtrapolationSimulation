#' BinomialCPP Class
#'
#' @description A class representing a Binomial Conditional Power Prior model.
#' @field power_parameter The power parameter for the model
#' @field method Method name
#' @field stan_prior Stan prior model
#' @field stan_prior_code Stan code for the prior
#' @export
BinomialCPP <- R6::R6Class(
  "BinomialCPP",
  inherit = MCMCModel,
  public = list(
    power_parameter = NULL,
    method = "CPP",
    stan_prior = NULL,
    stan_prior_code = "
        data {
          int<lower = 0> n_treatment_source; // number of patients in the treatment group of the source study
          int<lower = 0> n_control_source; // number of patients in the control group of the source study

          int<lower = 0, upper = n_treatment_source> successes_treatment_source;  // number of successes in the treatment group of the source study
          int<lower = 0, upper = n_control_source> successes_control_source; // number of successes in the control group of the source study

          real<lower = 0, upper = 1> power_parameter; // power parameter
        }
        parameters {
          real<lower = 0, upper = 1> control_rate_source;
          real<lower = 0, upper = 1> control_rate_target;
          real<lower = -fmin(control_rate_target, control_rate_source), upper = 1 - fmax(control_rate_target, control_rate_source)> target_treatment_effect; // The power prior approach assumes that the treatment effect is the same in the source and target studies
        }
        transformed parameters {
          real treatment_rate_source = control_rate_source + target_treatment_effect;
          real treatment_rate_target = control_rate_target + target_treatment_effect;
        }
        model {
          control_rate_target ~ uniform(0, 1);
          control_rate_source ~ uniform(0, 1);

          // Conditioning on the control rates constrains the distribution of the treatment effect
          target_treatment_effect ~ uniform(-fmin(control_rate_target, control_rate_source), 1 - fmax(control_rate_target, control_rate_source)); // Initial prior on the treatment effect


          target += power_parameter * (binomial_lpmf(successes_control_source |n_control_source, control_rate_source) +
          binomial_lpmf(successes_treatment_source |n_treatment_source, treatment_rate_source)); // Likelihood of the source study p(D_S|target_treatment_effect)
          // Note that here, we implicitly borrow from the two arms
        }
        ",

    #' @description Initialize an instance of the BinomialCPP class.
    #' @param prior The prior object
    #' @param mcmc_config The MCMC configuration parameters
    initialize = function(prior, mcmc_config) {
      if (!(prior$method_parameters$initial_prior[[1]] == "noninformative")) {
        stop("Only implemented for a noninformative initial prior")
      }

      super$initialize(prior = prior, mcmc_config = mcmc_config)

      self$summary_measure_likelihood <- "binomial"

      self$stan_model_code <- "
        data {
          int<lower = 0> n_treatment_source; // number of patients in the treatment group of the source study
          int<lower = 0> n_control_source; // number of patients in the control group of the source study

          int<lower = 0, upper = n_treatment_source> successes_treatment_source;  // number of successes in the treatment group of the source study
          int<lower = 0, upper = n_control_source> successes_control_source; // number of successes in the control group of the source study

          int<lower = 0> n_treatment_target; // number of patients in the treatment group of the target study
          int<lower = 0> n_control_target; // number of patients in the control group of the target study

          int<lower = 0, upper = n_treatment_target> successes_treatment_target;  // number of successes in the treatment group of the target study
          int<lower = 0, upper = n_control_target> successes_control_target; // number of successes in the control group of the target study

          real<lower = 0, upper = 1> power_parameter; // power parameter
        }
        parameters {
          real<lower = 0, upper = 1> control_rate_source;
          real<lower = 0, upper = 1> control_rate_target;
          real<lower = -fmin(control_rate_target, control_rate_source), upper = 1 - fmax(control_rate_target, control_rate_source)> target_treatment_effect; // The power prior approach assumes that the treatment effect is the same in the source and target studies
        }
        transformed parameters {
          real treatment_rate_source = control_rate_source + target_treatment_effect;
          real treatment_rate_target = control_rate_target + target_treatment_effect;
        }
        model {
          control_rate_target ~ uniform(0, 1);
          control_rate_source ~ uniform(0, 1);

          // Conditioning on the control rates constrains the distribution of the treatment effect
          target_treatment_effect ~ uniform(-fmin(control_rate_target, control_rate_source), 1 - fmax(control_rate_target, control_rate_source)); // Initial prior on the treatment effect


          target += power_parameter * (binomial_lpmf(successes_control_source |n_control_source, control_rate_source) +
          binomial_lpmf(successes_treatment_source |n_treatment_source, treatment_rate_source)); // Likelihood of the source study p(D_S|target_treatment_effect)
          // Note that here, we implicitly borrow from the two arms

          successes_control_target ~ binomial(n_control_target, control_rate_target);
          successes_treatment_target ~ binomial(n_treatment_target, treatment_rate_target);
        }
        "

      model_name <- "binomial_cpp"
      self$stan_model <- compile_stan_model(model_name, self$stan_model_code)

      self$power_parameter <- prior$method_parameters$power_parameter[[1]]
    },

    #' @description Prepare data for the model
    #' @param target_data The target study data
    #' @return A list of prepared data
    prepare_data = function(target_data) {
      # Specify the parameters
      data_list <- list(
        power_parameter = self$power_parameter,
        n_treatment_source = as.integer(self$prior$source$sample_size_treatment),
        n_control_source = as.integer(self$prior$source$sample_size_control),
        successes_treatment_source = as.integer(
          self$prior$source$sample_size_treatment * self$prior$source$treatment_rate
        ),
        successes_control_source = as.integer(
          self$prior$source$sample_size_control * self$prior$source$control_rate
        ),
        n_treatment_target = as.integer(target_data$sample_size_treatment),
        n_control_target = as.integer(target_data$sample_size_control),
        successes_treatment_target = as.integer(
          target_data$sample_size_treatment * target_data$sample$sample_treatment_rate
        ),
        successes_control_target = as.integer(
          target_data$sample_size_control * target_data$sample$sample_control_rate
        )
      )
      return(data_list)
    },

    #' @description Sample from the prior using Stan
    draw_mcmc_prior = function() {
      model_name <- "binomial_cpp_prior"
      self$stan_prior <- compile_stan_model(model_name, self$stan_prior_code)
      data_list <- list(
        power_parameter = self$power_parameter,
        n_treatment_source = as.integer(self$prior$source$sample_size_treatment),
        n_control_source = as.integer(self$prior$source$sample_size_control),
        successes_treatment_source = as.integer(
          self$prior$source$sample_size_treatment * self$prior$source$treatment_rate
        ),
        successes_control_source = as.integer(
          self$prior$source$sample_size_control * self$prior$source$control_rate
        )
      )

      # Sample from the posterior
      self$prior_draws <- self$stan_prior$sample(
        data = data_list,
        chains = self$mcmc_config$num_chains,
        parallel_chains = self$mcmc_config$parallel_chains,
        iter_sampling = self$mcmc_config$chain_length,
        iter_warmup = self$mcmc_config$tune,
        threads_per_chain = self$mcmc_config$threads_per_chain,
        output_dir = self$draws_dir
      )
    }
  )
)

#' Gaussian_NPP class
#'
#' @description This class represents a Gaussian model using the NPP (Noninformative Power Prior) approach.
#' It inherits from the Model class.
#'
#' @field power_parameter_mean The mean of the prior on the power parameter.
#' @field power_parameter_std The standard deviation of the prior on the power parameter.
#' @field target_treatment_effect_estimate Estimate of the treatment effect in the target study
#' @field target_treatment_effect_standard_error Standard error on the estimate of the treatment effect in the target study
#' @field p Parameter of the initial Beta prior on the power parameter
#' @field q Parameter of the initial Beta prior on the power parameter
#' @field prior_normconst Normalization constant of the prior
#' @field posterior_normconst Normalization constant of the posterior
#' @field max_posterior_pdf Maximum value of the posterior p.d.f., used for rejection sampling of the posterior p.d.f.
#' @field method Name of the method
#' @field summary_measure_likelihood Summary measure likelihood
#'
#' @export
Gaussian_NPP <- R6::R6Class(
  "Gaussian_NPP",
  inherit = Model,
  public = list(
    power_parameter_mean = NULL,
    power_parameter_std= NULL,
    target_treatment_effect_estimate = NULL,
    target_treatment_effect_standard_error = NULL,
    p = NULL,
    q = NULL,
    prior_normconst = NULL,
    posterior_normconst = NULL,
    max_posterior_pdf = NULL,
    # Max of the pdf
    method = "NPP",
    summary_measure_likelihood = "normal",

    #' @description Initialize a new Gaussian_NPP object.
    #' @param prior Prior object containing method parameters.
    #' @return A new Gaussian_NPP object.
    initialize = function(prior) {
      super$initialize()

      if (!(prior$method_parameters$initial_prior[[1]] == "noninformative")) {
        stop("Only implemented for a noninformative initial prior")
      }

      self$power_parameter_mean <- prior$method_parameters$power_parameter_mean[[1]]
      self$power_parameter_std <- prior$method_parameters$power_parameter_std[[1]]

      if (self$power_parameter_mean < 0 ||
          self$power_parameter_mean > 1) {
        stop("Invalid parameter for the mean of the prior on the power parameter")
      }

      if (self$power_parameter_std < 0) {
        stop("Invalid parameter for the std of the prior on the power parameter")
      }


      if (self$power_parameter_std > 0.5 | self$power_parameter_std <= 0){
        stop("The standard deviation of a Beta distribution is between 0 and 0.5.")
      }

      if (self$power_parameter_mean > 1 | self$power_parameter_mean < 0){
        stop("The mean of a Beta distribution is between 0 and 1.")
      }
      xi <- self$power_parameter_mean
      omega <- self$power_parameter_std ^ 2 / (xi * (1 - xi) - self$power_parameter_std ^ 2)
      self$p <- xi / omega
      self$q <- (1 - xi) / omega

      # Because the behavior of the beta distribution changes a lot at 0 and 1 for p and q above or below 1, we do the following to avoid issues due to limited numerical precision:
      tolerance <- .Machine$double.eps^0.5
      if (isTRUE(all.equal(self$p, 1, tolerance = tolerance))){
        self$p <- 1
      }
      if (isTRUE(all.equal(self$q, 1, tolerance = tolerance))){
        self$q <- 1
      }

      if (self$p <=0 | self$q <= 0){
        stop("Parameters p and q of the Beta distribution must be >0.")
      }
      self$n_components_mixture_approx <- seq(1, 2)
    },

    #' @description Calculate the unnormalized posterior power parameter PDF.
    #' @param power_parameter Power parameter value.
    #' @param target_data Target study data.
    #' @return Unnormalized posterior power parameter PDF value.
    unnormalized_posterior_power_parameter_pdf = function(power_parameter, target_data) {
      # Numerator of p(\power_parameter|D_T, D_S)

      pdf <- numeric(length(power_parameter))

      for (i in seq(1, length(power_parameter))){
        pdf[i] <- dnorm(
          x = target_data$sample$treatment_effect_estimate,
          mean = self$prior$source$treatment_effect_estimate,
          sd = sqrt(
            target_data$sample$treatment_effect_standard_error ^ 2 + self$prior$source$standard_error ^
              2 / power_parameter[i]
          )
        ) *
          dbeta(
            x = power_parameter[i],
            shape1 = self$p,
            shape2 = self$q
          )
      }
      if (any(is.na(pdf))){
        stop("Unnormalized posterior PDF of the power parameter is evaluated to NA.")
      }
      return(pdf)
    },

    #' @description Return the posterior distribution of the power parameter
    #' @param power_parameter Power parameter value.
    #' @param target_data Target study data.
    #' @return Posterior power parameter PDF value.
    power_parameter_posterior_pdf = function(power_parameter, target_data) {
      normconst <- self$normalizing_constant_power_parameter(target_data)
      pdf <- self$unnormalized_posterior_power_parameter_pdf(power_parameter, target_data)


      norm_pdf <- pdf / normconst

      if (any(is.na(norm_pdf))){
        stop("Normalized posterior PDF of the power parameter is evaluated to NA.")
      }

      return(norm_pdf)
    },

    #' @description Calculate the normalizing constant for the power parameter.
    #' @param target_data Target study data.
    #' @return Normalizing constant value.
    normalizing_constant_power_parameter = function(target_data) {
      # Numerator of p(\power_parameter|D_T, D_S)
      intFun <- function(power_parameter) {
        return(
          self$unnormalized_posterior_power_parameter_pdf(power_parameter, target_data)
        )
      }

      # Normalizing constant of p(\power_parameter|D_T, D_S)
      normconst <- integrate(
        f = intFun,
        lower = 0,
        upper = 1,
        rel.tol = 1e-8,
        stop.on.error = FALSE
      )$value

      if (is.na(normconst)){
        stop("Normalizing constant of the power parameter posterior PDF is NA.")
      }

      if (normconst == 0) {
        stop("Normalization constant of the power parameter posterior PDF is zero")
      }
      return(normconst)
    },

    #' @description Perform inference on the target data.
    #' @param target_data Target study data.
    #' @return A string indicating the success status of the inference.
    inference = function(target_data) {
      # Return the mean of the posterior distribution of power_parameter:
      self$prior_normconst <- self$normalizing_constant_power_parameter(target_data)

      self$target_treatment_effect_estimate <- target_data$sample$treatment_effect_estimate
      self$target_treatment_effect_standard_error <- target_data$sample$treatment_effect_standard_error

      self$posterior_parameters$power_parameter_mean <- integrate(function(power_parameter) {
        power_parameter * self$unnormalized_posterior_power_parameter_pdf(power_parameter, target_data = target_data)
      },
      lower = 0,
      upper = 1,
      rel.tol = 1e-8)$value / self$prior_normconst

      power_parameter_var <- integrate(function(power_parameter) {
        power_parameter ^ 2 * self$unnormalized_posterior_power_parameter_pdf(power_parameter, target_data = target_data)
      },
      lower = 0,
      upper = 1,
      rel.tol = 1e-8)$value/self$prior_normconst - self$posterior_parameters$power_parameter_mean ^
        2

      self$posterior_parameters$power_parameter_std <- sqrt(power_parameter_var)

      # Numerical integration for mean

      bound <- Inf # 30 Using -Inf/Inf may result in error because of non-finite function value
      self$post_mean <- integrate(
        function(theta) {
          theta * self$posterior_pdf(theta)
        },
        lower = -bound,
        upper = bound,
        rel.tol = 1e-8,
        stop.on.error = FALSE
      )$value


      # Numerical integration for variance
      self$post_var <- integrate(
        function(theta) {
          theta ^ 2 * self$posterior_pdf(theta)
        },
        lower = -bound, # Using -Inf/Inf would result in error because of non-finite function value
        upper = bound,
        rel.tol = 1e-8,
        stop.on.error = FALSE
      )$value - self$post_mean ^ 2

      if (self$post_var == 0){
        stop("Posterior variance is 0. Most likely due to numerical errors when computing the integral. One solution is to adapt the integration range.")
      }

      self$post_median <- stats::uniroot(
        function(theta) {
          self$posterior_cdf(x = theta) - 0.5
        },
        interval = c(
          self$post_mean - 3 * sqrt(self$post_var),
          self$post_mean + 3 * sqrt(self$post_var)
        ),
        tol = 0.001
      )$root

      self$max_posterior_pdf <- -optimize(
        function(x)
          - self$posterior_pdf(x),
        interval = c(
          self$post_mean - 3 * sqrt(self$post_var),
          self$post_mean + 3 * sqrt(self$post_var)
        )
      )$objective

      fit_success <- "Success"
      return(fit_success)
    },

    #' @description Calculate the credible interval.
    #' @param level Credible interval level (default: 0.95).
    #' @return A vector containing the lower and upper bounds of the credible interval.
    credible_interval = function(level = 0.95) {
      hdi <- HDInterval::inverseCDF(c((1-level)/2, (1+level)/2), function(x) {
        self$posterior_cdf(x)
      })
      return(hdi)
    },

    #' @description Return the median of the posterior distribution.
    #' @param ... Optional arguments
    #' @return Median of the posterior distribution.
    posterior_median = function(...) {
      return(self$post_median)
    },

    #' @description Sample from the posterior distribution.
    #' @param n_samples Number of samples to draw from the posterior distribution.
    #' @return A vector of samples from the posterior distribution.
    sample_posterior = function(n_samples) {
      # Sample from the posterior using rejection sampling

      # Define a proposal distribution
      q <- function(x) {
        dnorm(x, self$post_mean, sqrt(self$post_var))
      }

      samples <- numeric(n_samples)
      count <- 0

      # Constant c such that c * q(x) >= p(x)
      c <- self$max_posterior_pdf / q(self$post_mean)

      while (count < n_samples) {
        x_proposal <- rnorm(1, self$post_mean, sqrt(self$post_var)) # Sampling from the proposal distribution q(x)
        u <- runif(1) # Uniform random number [0, 1]
        # Acceptance check
        if (u < self$posterior_pdf(x_proposal) / (c * q(x_proposal))) {
          count <- count + 1
          samples[count] <- x_proposal
        }
      }
      return(samples)
    },

    #' @description Sample from the prior distribution.
    #' @param n_samples Number of samples to draw from the prior distribution.
    #' @return A vector of samples from the prior distribution.
    sample_prior = function(n_samples) {
      # Sample from the beta prior
      power_parameter_samples <- rbeta(n_samples, shape1 = self$p, shape2 = self$q)

      samples <- numeric(n_samples)

      # Sample theta from the normal distribution for each sampled power_parameter
      for (i in 1:n_samples) {
        sigma <- sqrt(self$prior$source$standard_error ^ 2 / power_parameter_samples[i])
        samples[i] <- rnorm(1,
                            mean = self$prior$source$treatment_effect_estimate,
                            sd = sigma)
      }
      return(samples)
    },

    #' @description Calculate the posterior cumulative distribution function (CDF).
    #' @param x Vector of points to evaluate the CDF at.
    #' @return Vector of CDF values corresponding to the input points.
    posterior_cdf = function(x) {
      integral <- numeric(length(x))

      for (i in seq(1, length(x))) {
        xi <- x[i]

        # The integrate function has bad behavior for large upper bound, with a non-monotonic behavior despite a positive integrand. The following decomposition solves the issue:
        integral_part_1 <- integrate(
          self$posterior_pdf,
          lower = -Inf,
          upper = 0,
          rel.tol = 1e-8,
          stop.on.error = FALSE
        )$value

        integral_part_2 <- integrate(
          self$posterior_pdf,
          lower = 0,
          upper = xi,
          rel.tol = 1e-8,
          stop.on.error = FALSE
        )$value

        integral[i] <- integral_part_1 + integral_part_2
      }

      # The value of the CDF may be inaccurate for large x values due to double integration.
      return(integral)
    },

    #' @description Calculate the posterior probability density function (PDF).
    #' @param x Vector of points to evaluate the PDF at.
    #' @return Vector of PDF values corresponding to the input points.
    posterior_pdf = function(x) {
      # marginal posterior of theta p(\theta| D_T, D_S): integration of p(\theta, power_parameter | D_T, D_S) with respect to power_parameter
      if (is.null(self$prior_normconst)) {
        stop("Call self$inference(target_data) before accessing posterior_pdf().")
      }

      integral <- numeric(length(x))

      for (i in seq(1, length(x))) {
        xi <- x[i]
        integral[i] <- integrate(
          f = function(power_parameter) {
            dnorm(
              x = self$target_treatment_effect_estimate,
              mean = xi,
              sd = self$target_treatment_effect_standard_error
            ) *
              dnorm(
                x = xi,
                mean = self$prior$source$treatment_effect_estimate,
                sd = sqrt(self$prior$source$standard_error ^ 2 / power_parameter)
              ) *
              dbeta(x = power_parameter,
                    shape1 = self$p,
                    shape2 = self$q)
          },
          lower = 0,
          upper = 1,
          rel.tol = 1e-8,
          stop.on.error = FALSE
        )$value
      }
      return(integral / self$prior_normconst)
    },

    #' @description Calculate the prior probability density function (PDF).
    #' @param x Vector of points to evaluate the PDF at.
    #' @return Vector of PDF values corresponding to the input points.
    prior_pdf = function(x) {
      # p(\theta|D_S): integration of p(\theta, power_parameter | D_S) with respect to power_parameter

      integral <- numeric(length(x))

      for (i in seq(1, length(x))) {
        xi <- x[i]
        integral[i] <- integrate(
          f = function(power_parameter) {
            dnorm(
              x = xi,
              mean = self$prior$source$treatment_effect_estimate,
              sd = sqrt(self$prior$source$standard_error ^ 2 / power_parameter)
            ) *
              dbeta(x = power_parameter,
                    shape1 = self$p,
                    shape2 = self$q)
          },
          lower = 0,
          upper = 1,
          rel.tol = 1e-8,
          stop.on.error = FALSE
        )$value
      }
      return(integral)
    },

    #' @description Plot posterior probability density function (PDF) of the power parameter
    #' @param target_data Target study data
    #' @return A plot
    plot_power_parameter_posterior_pdf = function(target_data) {
      x_values <- seq(0, 1, length.out = 100)

      posterior_pdf <- model$power_parameter_posterior_pdf(x_values, target_data)

      df <- data.frame(x = x_values, posterior_pdf = posterior_pdf)

      plt <- ggplot2::ggplot(df, ggplot2::aes(x = x_values)) +
        geom_line(ggplot2::aes(y = posterior_pdf),
                  color = "red",
                  size = 1) +
        ggplot2::labs(
          title = "Posterior distribution of the power parameter",
          x = "Power parameter",
          y = "Density",
          labels = c("Posterior pdf")
        ) +
        theme_minimal() +
        scale_color_manual(values = c("red")) +
        guides(color = guide_legend(title = NULL)) +
        scale_fill_manual(
          name = "PDF",
          labels = c("Posterior pdf"),
          values = c("red")
        )
      return(plt)
    },

    #' @description
    #' Plot the power parameter as a function of drift in treatment effect
    #' @param source_treatment_effect_estimate Treatment effect estimate in the source study
    #' @param target_data Target study data
    #' @param min_drift Minimum drift value
    #' @param max_drift Maximum drift value
    #' @param resolution Number of points on the drift grid.
    #' @return The prior PDF.
    plot_power_parameter_vs_drift = function(source_treatment_effect_estimate,
                                             target_data,
                                             min_drift,
                                             max_drift,
                                             resolution) {
      # Generate a sequence of differences
      drifts <- seq(min_drift, max_drift, length.out = resolution)

      # Calculate power parameters for each difference
      power_parameters <- sapply(drifts, function(drift) {
        target_data_drift <- duplicate(target_data)
        target_data_drift$sample$treatment_effect_estimate <- source_treatment_effect_estimate + drift
        self$power_parameter_estimation(target_data)
      })


      # Create a data frame for plotting
      plot_data <- data.frame(Drift = drifts, PowerParameter = power_parameters)

      # Plot using ggplot2
      ggplot(plot_data, aes(x = Drift, y = PowerParameter)) +
        geom_line(color = "blue") +
        labs(title = "Power Parameter as a Function of Drift", x = "Drift", y = "Power Parameter") +
        theme_minimal()
    }
  )
)

g_function <- function(log_tau) {
  return(pmax(log_tau, 1))
}

#' GaussianCommensuratePowerPrior class
#'
#' @description This class represents a Gaussian Commensurate Power Prior model.
#'
#' @field method Method name
#' @field heterogeneity_prior_family Heterogeneity prior family (half_normal, inverse_gamma)
#'
#' @export
GaussianCommensuratePowerPrior <- R6::R6Class(
  "GaussianCommensuratePowerPrior",
  inherit = MCMCModel,
  public = list(
    method = "commensurate_power_prior",
    heterogeneity_prior_family = NULL,

    #' @description Initialize the GaussianCommensuratePowerPrior object
    #'
    #' @param prior The prior object
    #' @param mcmc_config The MCMC configuration parameters
    #'
    initialize = function(prior, mcmc_config) {
      if (!(prior$method_parameters$initial_prior[[1]] == "noninformative")) {
        stop("Only implemented for a noninformative initial prior")
      }
      super$initialize(prior, mcmc_config)
      self$summary_measure_likelihood <- "normal"
      self$stan_model_code <- "
                                      functions {
                                        real g_function(real log_tau) {
                                          return fmax(log_tau, 1); // Define g(τ)
                                        }
                                      }

                                      data {
                                        int<lower=0,upper=2> prior_type; // Type of prior on the commensurability parameter; 0 for inverse_gamma, 1 for half normal, 2 for Cauchy
                                        int<lower=0> NS;       // Number of samples in dataset S
                                        int<lower=0> NT;       // Number of samples in dataset T
                                        real<lower=0> target_sampling_variance;  // Known variance of the target data
                                        real<lower=0> prior_variance; // Prior variance component
                                        real source_treatment_effect_estimate;      // Estimate from dataset S
                                        real target_treatment_effect_estimate;      // Estimate from dataset T
                                        real std_dev; // Parameter for the half normal prior
                                        real location; // Location parameter for the Cauchy prior
                                        real scale; // Scale parameter for the Cauchy prior
                                        real alpha; // Parameter for the inverse_gamma prior
                                        real beta; // Pocation parameter for the inverse_gamma prior
                                      }

                                      parameters {
                                        real<lower=0, upper=1> power_parameter;   // Power parameter
                                        real<lower=0> tau;     // Commensurability parameter
                                        real target_treatment_effect;          // Target treatment effect
                                      }

                                      transformed parameters {
                                        real<lower=0>  u = power_parameter * NS + tau * prior_variance;
                                        real<lower=0> tau2 = tau^2;
                                        real log_tau = log(tau);

                                      }

                                      model {
                                        // Priors
                                        if (prior_type == 0) {
                                          tau2 ~ inv_gamma(alpha, beta);
                                        } else if (prior_type == 1){
                                          tau ~ normal(0, std_dev) T[0, ];
                                        } else if (prior_type == 2){
                                          log_tau ~ cauchy(location, scale);
                                        }


                                        power_parameter ~ beta(g_function(log_tau), 1);

                                        target_treatment_effect ~ normal(
    (power_parameter * NS * tau * target_sampling_variance * source_treatment_effect_estimate + NT * u * target_treatment_effect_estimate) /
    (power_parameter * NS * tau * target_sampling_variance + NT * u),
    sqrt((u * target_sampling_variance) / (power_parameter * NS * tau * target_sampling_variance + NT * u)));
                                      }
                                      "

      self$heterogeneity_prior_family <- prior$method_parameters$heterogeneity_prior$family

      model_name <- paste0("gaussian_commensurate_pp_", self$heterogeneity_prior_family)
      self$stan_model <- compile_stan_model(model_name, self$stan_model_code)

      self$posterior_parameters <- list(
        heterogeneity_parameter_mean = NA,
        heterogeneity_parameter_std = NA,
        power_parameter_mean = NA,
        power_parameter_std = NA
      )
    },
    #' @description Prepare the data to be used by CmdStanR
    #' @param target_data Target study data
    prepare_data = function(target_data) {
      # Specify the parameters
      data <- list(
        NS = as.integer(self$prior$source$equivalent_source_sample_size_per_arm),
        # Number of samples in dataset S
        NT = as.integer(target_data$sample_size_per_arm),
        # Number of samples in dataset T
        target_sampling_variance = target_data$sample$treatment_effect_standard_error ^
          2 * target_data$sample_size_per_arm,
        # Known variance
        prior_variance = self$prior$source$standard_error ^ 2 * self$prior$source$equivalent_source_sample_size_per_arm,
        # Source variance
        source_treatment_effect_estimate = self$prior$source$treatment_effect_estimate,
        # Estimate from source dataset
        target_treatment_effect_estimate = target_data$sample$treatment_effect_estimate # Estimate from target dataset
      )

      heterogeneity_prior <- self$prior$method_parameters$heterogeneity_prior

      # Initialize parameters list based on the prior family
      parameters <- switch(
        self$heterogeneity_prior_family,
        half_normal = list(
          prior_type = 1,
          std_dev = heterogeneity_prior$std_dev,
          location = 0,
          scale = 0,
          alpha = 0,
          beta = 0
        ),
        cauchy = list(
          prior_type = 2,
          location = heterogeneity_prior$location,
          scale = heterogeneity_prior$scale,
          std_dev = 0,
          alpha = 0,
          beta = 0
        ),
        inverse_gamma = list(
          prior_type = 0,
          alpha = heterogeneity_prior$alpha,
          beta = heterogeneity_prior$beta,
          location = 0,
          scale = 0,
          std_dev = 0
        ),
        stop("Prior type not implemented")
      )
      data <- c(data, parameters)

      return(data)
    },
    #' @description Compute posterior parameters
    compute_posterior_parameters = function() {
      self$posterior_parameters <- list(
        heterogeneity_parameter_mean = self$fit_summary[self$fit_summary$variable == "tau", "mean"] %>% dplyr::pull(mean),
        heterogeneity_parameter_std = self$fit_summary[self$fit_summary$variable == "tau", "sd"] %>% dplyr::pull(sd),
        power_parameter_mean = self$fit_summary[self$fit_summary$variable == "power_parameter", "mean"] %>% dplyr::pull(mean),
        power_parameter_std = self$fit_summary[self$fit_summary$variable == "power_parameter", "sd"] %>% dplyr::pull(sd)
      )
    },
    #' @description Draw samples from the prior distribution. Based on equation 8 in Hobbs et al (2011).
    #' @param n_samples Number of samples
    sample_prior = function(n_samples) {
      heterogeneity_prior <- self$prior$method_parameters$heterogeneity_prior
      if (self$heterogeneity_prior_family == "cauchy"){
        log_tau_samples <- rcauchy(n = n_samples, location = heterogeneity_prior$location, scale = heterogeneity_prior$scale)
        tau_samples <- exp(log_tau_samples)
      } else if (self$heterogeneity_prior_family == "half_normal"){
        tau_samples <- extraDistr::rhnorm(n = n_samples, sigma =  heterogeneity_prior$std_dev)
        log_tau_samples <- log(tau_samples)
      } else if (self$heterogeneity_prior_family == "inverse_gamma"){
        # Sample the power_parameter from the Beta distribution
        tau2_samples <- extraDistr::rinvgamma(n = n_samples, alpha = heterogeneity_prior$alpha, beta = heterogeneity_prior$beta)
        tau_samples <- sqrt(tau2_samples)
        log_tau_samples <- log(sqrt(tau2_samples))
      } else {
        stop(paste0("Heterogeneity prior not implemented : ", self$heterogeneity_prior_family))
      }

      power_parameter_samples <- rbeta(n_samples, g_function(log_tau_samples), 1) # power_parameter ~ beta(g_function(log_tau), 1)

      prior_variance <- self$prior$source$standard_error ^ 2 * self$prior$source$equivalent_source_sample_size_per_arm

      NS <- as.integer(self$prior$source$equivalent_source_sample_size_per_arm)

      std <- sqrt(1 / tau_samples + prior_variance /
                    (power_parameter_samples * NS))

      # Filter out Inf elements
      std <- sapply(std, function(x)
        if (!is.infinite(x))
          x
        else
          NULL)
      # Removing NULL elements from the list (which replaced the Inf elements)
      std <- unlist(std[!sapply(std, is.null)])

      # Compute the target treatment effect for each sample
      treatment_effect_samples <- rnorm(
        n_samples,
        mean = self$prior$source$treatment_effect_estimate,
        sd = std
      )

      if (any(is.na(treatment_effect_samples))) {
        stop("Some samples are NA.")
      }

      return(treatment_effect_samples)
    },
    #' @description Joint prior p.d.f. Based on equation (8) in Hobbs et al (2011).
    #' @param treatment_effect Treatment effect
    #' @param gamma Power parameter
    #' @param tau Heterogeneity parameter
    joint_prior_pdf = function(treatment_effect, gamma, tau) {
      source_variance <- self$prior$source$standard_error ^ 2 * self$prior$source$equivalent_source_sample_size_per_arm

      NS <- self$prior$source$equivalent_source_sample_size_per_arm

      # Compute the normal distribution component
      sigma_T_squared <- 1 / tau + (source_variance / (gamma * NS))
      normal_component <- dnorm(
        treatment_effect,
        mean = self$prior$source$treatment_effect_estimate,
        sd = sqrt(sigma_T_squared)
      )

      # Compute the Cauchy prior on log(tau)
      log_tau <- log(tau)

      # Compute the beta distribution component
      g_tau <- g_function(log_tau)
      beta_component <- dbeta(gamma, shape1 = g_tau, shape2 = 1)

      heterogeneity_prior <- self$prior$method_parameters$heterogeneity_prior

      if (self$heterogeneity_prior_family == "cauchy"){
        prior_tau <- dcauchy(log_tau,
                             location = heterogeneity_prior$location,
                             scale = heterogeneity_prior$scale) / tau # Adjust for the Jacobian of the transformation

      } else if (self$heterogeneity_prior_family ==  "inverse_gamma"){
        tau2 <- tau ^ 2

        prior_tau <- 2 * tau * extraDistr::dinvgamma(tau2,
                                         alpha = heterogeneity_prior$alpha,
                                         beta = 1 / heterogeneity_prior$beta) # Adjust for the Jacobian of the transformation
      } else if (self$heterogeneity_prior_family =="half_normal"){
        prior_tau <- extraDistr::dhnorm(tau, sigma = heterogeneity_prior$std_dev)
      } else {
        stop(paste0("Heterogeneity prior not implemented : ", self$heterogeneity_prior_family))
      }

      # Return the product of components
      return(normal_component * beta_component * prior_tau)
    },

    #' @description Integrate out gamma and tau to get the marginal PDF for treatment_effect
    #' @param treatment_effect Treatment effect
    unnormalized_prior_pdf = function(treatment_effect) {
      return(
        pracma::integral2(
          function(gamma, tau) {
            self$joint_prior_pdf(treatment_effect, gamma, tau)
          },
          xmin = 0,
          xmax = 1,
          ymin = 0.001,
          ymax = 100,
          reltol = 1e-6
        )$Q
      ) # Double integration
    },
    #' @description Prior p.d.f. of the treatment effect
    #' @param treatment_effect Treatment effect
    prior_pdf = function(treatment_effect) {
      # Vectorize the marginal PDF function
      unnormalized_prior_pdf_vec <- Vectorize(self$unnormalized_prior_pdf)

      # Normalize the marginal PDF
      normalizing_constant <- integrate(unnormalized_prior_pdf_vec,
                                        lower = -Inf,
                                        upper = Inf)$value
      normalized_pdf_treatment_effect <- unnormalized_prior_pdf_vec(treatment_effect) / normalizing_constant
      return(normalized_pdf_treatment_effect)
    }
  )
)
