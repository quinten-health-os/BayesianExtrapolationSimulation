% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/method_model.R
\name{MCMCModel}
\alias{MCMCModel}
\title{MCMCModel class}
\description{
This class represents a Bayesian borrowing model using MCMC sampling.
It inherits from the Model class.
}
\section{Super class}{
\code{\link[RBExT:Model]{RBExT::Model}} -> \code{MCMCModel}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{stan_model_code}}{Code of the Stan model}

\item{\code{stan_model}}{The compiled Stan model}

\item{\code{fit}}{The MCMC fit object}

\item{\code{fit_summary}}{Summary of the Stan fit}

\item{\code{treatment_effect_summary}}{Summary statistics of the treatment effect posterior distribution}

\item{\code{credible_interval_97.5}}{The upper bound of the credible interval}

\item{\code{credible_interval_2.5}}{The lower bound of the credible interval}

\item{\code{mcmc_config}}{The MCMC configuration parameters}

\item{\code{mcmc_ess}}{MCMC ESS}

\item{\code{n_divergences}}{Number of divergences in MCMC inference}

\item{\code{rhat}}{r-hat statistics}

\item{\code{draws_dir}}{Directory where to store MCMC draws (used by Stan)}

\item{\code{prior_draws}}{Draws from the prior}

\item{\code{prior_pdf_approx}}{Approximation to the prior probability density function}

\item{\code{prior_cdf_approx}}{Approximation to the prior cumulative density function}

\item{\code{posterior_pdf_approx}}{Approximation to the posterior probability density function}

\item{\code{posterior_cdf_approx}}{Approximation to the posterior cumulative density function}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-MCMCModel-new}{\code{MCMCModel$new()}}
\item \href{#method-MCMCModel-check_mcmc_config}{\code{MCMCModel$check_mcmc_config()}}
\item \href{#method-MCMCModel-prepare_data}{\code{MCMCModel$prepare_data()}}
\item \href{#method-MCMCModel-inference}{\code{MCMCModel$inference()}}
\item \href{#method-MCMCModel-credible_interval}{\code{MCMCModel$credible_interval()}}
\item \href{#method-MCMCModel-posterior_median}{\code{MCMCModel$posterior_median()}}
\item \href{#method-MCMCModel-sample_posterior}{\code{MCMCModel$sample_posterior()}}
\item \href{#method-MCMCModel-check_data}{\code{MCMCModel$check_data()}}
\item \href{#method-MCMCModel-compute_posterior_parameters}{\code{MCMCModel$compute_posterior_parameters()}}
\item \href{#method-MCMCModel-draw_mcmc_prior}{\code{MCMCModel$draw_mcmc_prior()}}
\item \href{#method-MCMCModel-posterior_pdf}{\code{MCMCModel$posterior_pdf()}}
\item \href{#method-MCMCModel-posterior_cdf}{\code{MCMCModel$posterior_cdf()}}
\item \href{#method-MCMCModel-prior_pdf}{\code{MCMCModel$prior_pdf()}}
\item \href{#method-MCMCModel-prior_cdf}{\code{MCMCModel$prior_cdf()}}
\item \href{#method-MCMCModel-sample_prior}{\code{MCMCModel$sample_prior()}}
\item \href{#method-MCMCModel-clone}{\code{MCMCModel$clone()}}
}
}
\if{html}{\out{
<details><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="RBExT" data-topic="Model" data-id="create"><a href='../../RBExT/html/Model.html#method-Model-create'><code>RBExT::Model$create()</code></a></span></li>
<li><span class="pkg-link" data-pkg="RBExT" data-topic="Model" data-id="empirical_bayes_update"><a href='../../RBExT/html/Model.html#method-Model-empirical_bayes_update'><code>RBExT::Model$empirical_bayes_update()</code></a></span></li>
<li><span class="pkg-link" data-pkg="RBExT" data-topic="Model" data-id="estimate_bayesian_operating_characteristics"><a href='../../RBExT/html/Model.html#method-Model-estimate_bayesian_operating_characteristics'><code>RBExT::Model$estimate_bayesian_operating_characteristics()</code></a></span></li>
<li><span class="pkg-link" data-pkg="RBExT" data-topic="Model" data-id="estimate_frequentist_operating_characteristics"><a href='../../RBExT/html/Model.html#method-Model-estimate_frequentist_operating_characteristics'><code>RBExT::Model$estimate_frequentist_operating_characteristics()</code></a></span></li>
<li><span class="pkg-link" data-pkg="RBExT" data-topic="Model" data-id="plot_pdfs"><a href='../../RBExT/html/Model.html#method-Model-plot_pdfs'><code>RBExT::Model$plot_pdfs()</code></a></span></li>
<li><span class="pkg-link" data-pkg="RBExT" data-topic="Model" data-id="plot_posterior_pdf"><a href='../../RBExT/html/Model.html#method-Model-plot_posterior_pdf'><code>RBExT::Model$plot_posterior_pdf()</code></a></span></li>
<li><span class="pkg-link" data-pkg="RBExT" data-topic="Model" data-id="plot_prior_pdf"><a href='../../RBExT/html/Model.html#method-Model-plot_prior_pdf'><code>RBExT::Model$plot_prior_pdf()</code></a></span></li>
<li><span class="pkg-link" data-pkg="RBExT" data-topic="Model" data-id="posterior_mean"><a href='../../RBExT/html/Model.html#method-Model-posterior_mean'><code>RBExT::Model$posterior_mean()</code></a></span></li>
<li><span class="pkg-link" data-pkg="RBExT" data-topic="Model" data-id="posterior_moments"><a href='../../RBExT/html/Model.html#method-Model-posterior_moments'><code>RBExT::Model$posterior_moments()</code></a></span></li>
<li><span class="pkg-link" data-pkg="RBExT" data-topic="Model" data-id="posterior_quantile"><a href='../../RBExT/html/Model.html#method-Model-posterior_quantile'><code>RBExT::Model$posterior_quantile()</code></a></span></li>
<li><span class="pkg-link" data-pkg="RBExT" data-topic="Model" data-id="posterior_to_RBesT"><a href='../../RBExT/html/Model.html#method-Model-posterior_to_RBesT'><code>RBExT::Model$posterior_to_RBesT()</code></a></span></li>
<li><span class="pkg-link" data-pkg="RBExT" data-topic="Model" data-id="prior_ESS"><a href='../../RBExT/html/Model.html#method-Model-prior_ESS'><code>RBExT::Model$prior_ESS()</code></a></span></li>
<li><span class="pkg-link" data-pkg="RBExT" data-topic="Model" data-id="prior_to_RBesT"><a href='../../RBExT/html/Model.html#method-Model-prior_to_RBesT'><code>RBExT::Model$prior_to_RBesT()</code></a></span></li>
<li><span class="pkg-link" data-pkg="RBExT" data-topic="Model" data-id="prior_treatment_benefit"><a href='../../RBExT/html/Model.html#method-Model-prior_treatment_benefit'><code>RBExT::Model$prior_treatment_benefit()</code></a></span></li>
<li><span class="pkg-link" data-pkg="RBExT" data-topic="Model" data-id="simulation_for_given_treatment_effect"><a href='../../RBExT/html/Model.html#method-Model-simulation_for_given_treatment_effect'><code>RBExT::Model$simulation_for_given_treatment_effect()</code></a></span></li>
<li><span class="pkg-link" data-pkg="RBExT" data-topic="Model" data-id="test_decision"><a href='../../RBExT/html/Model.html#method-Model-test_decision'><code>RBExT::Model$test_decision()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MCMCModel-new"></a>}}
\if{latex}{\out{\hypertarget{method-MCMCModel-new}{}}}
\subsection{Method \code{new()}}{
Initialize the MCMCModel object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MCMCModel$new(prior, mcmc_config)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{prior}}{The prior object}

\item{\code{mcmc_config}}{The MCMC configuration parameters}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MCMCModel-check_mcmc_config"></a>}}
\if{latex}{\out{\hypertarget{method-MCMCModel-check_mcmc_config}{}}}
\subsection{Method \code{check_mcmc_config()}}{
Check validity of the MCMC configuration
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MCMCModel$check_mcmc_config(mcmc_config)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{mcmc_config}}{MCMC configuration}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MCMCModel-prepare_data"></a>}}
\if{latex}{\out{\hypertarget{method-MCMCModel-prepare_data}{}}}
\subsection{Method \code{prepare_data()}}{
Prepare the data for inference.
Subclasses must implement the 'prepare_data' method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MCMCModel$prepare_data(target_data)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{target_data}}{The target data for inference}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MCMCModel-inference"></a>}}
\if{latex}{\out{\hypertarget{method-MCMCModel-inference}{}}}
\subsection{Method \code{inference()}}{
Perform inference using MCMC sampling
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MCMCModel$inference(target_data)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{target_data}}{The target data for inference}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MCMCModel-credible_interval"></a>}}
\if{latex}{\out{\hypertarget{method-MCMCModel-credible_interval}{}}}
\subsection{Method \code{credible_interval()}}{
Calculate the credible interval
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MCMCModel$credible_interval(level = 0.95)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{level}}{The confidence level for the credible interval (default is 0.95)}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The credible interval as a numeric vector
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MCMCModel-posterior_median"></a>}}
\if{latex}{\out{\hypertarget{method-MCMCModel-posterior_median}{}}}
\subsection{Method \code{posterior_median()}}{
Get the posterior median
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MCMCModel$posterior_median(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{Optional argument}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The posterior median as a numeric value
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MCMCModel-sample_posterior"></a>}}
\if{latex}{\out{\hypertarget{method-MCMCModel-sample_posterior}{}}}
\subsection{Method \code{sample_posterior()}}{
Sample from the posterior distribution
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MCMCModel$sample_posterior(n_samples)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n_samples}}{The number of samples to draw from the posterior distribution}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The sampled treatment effect values as a numeric vector
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MCMCModel-check_data"></a>}}
\if{latex}{\out{\hypertarget{method-MCMCModel-check_data}{}}}
\subsection{Method \code{check_data()}}{
Check the validity of the data
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MCMCModel$check_data(data_list)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{data_list}}{The list of data elements}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MCMCModel-compute_posterior_parameters"></a>}}
\if{latex}{\out{\hypertarget{method-MCMCModel-compute_posterior_parameters}{}}}
\subsection{Method \code{compute_posterior_parameters()}}{
Compute the posterior parameters. If there are posterior borrowing parameters,
the following method must be overriden in the subclass.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MCMCModel$compute_posterior_parameters()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MCMCModel-draw_mcmc_prior"></a>}}
\if{latex}{\out{\hypertarget{method-MCMCModel-draw_mcmc_prior}{}}}
\subsection{Method \code{draw_mcmc_prior()}}{
Draw samples from the prior using MCMC
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MCMCModel$draw_mcmc_prior()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MCMCModel-posterior_pdf"></a>}}
\if{latex}{\out{\hypertarget{method-MCMCModel-posterior_pdf}{}}}
\subsection{Method \code{posterior_pdf()}}{
Posterior PDF
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MCMCModel$posterior_pdf(target_treatment_effect)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{target_treatment_effect}}{Point at which to evaluate the posterior PDF}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MCMCModel-posterior_cdf"></a>}}
\if{latex}{\out{\hypertarget{method-MCMCModel-posterior_cdf}{}}}
\subsection{Method \code{posterior_cdf()}}{
Calculates the posterior cumulative distribution function (CDF) for a given target treatment effect.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MCMCModel$posterior_cdf(target_treatment_effect)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{target_treatment_effect}}{The target treatment effect.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The posterior CDF.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MCMCModel-prior_pdf"></a>}}
\if{latex}{\out{\hypertarget{method-MCMCModel-prior_pdf}{}}}
\subsection{Method \code{prior_pdf()}}{
Prior PDF
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MCMCModel$prior_pdf(
  target_treatment_effect,
  n_samples_quantile_estimation = 10000
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{target_treatment_effect}}{Point at which to evaluate the prior PDF}

\item{\code{n_samples_quantile_estimation}}{Number of samples used to estimate the quantiles of the distribution}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MCMCModel-prior_cdf"></a>}}
\if{latex}{\out{\hypertarget{method-MCMCModel-prior_cdf}{}}}
\subsection{Method \code{prior_cdf()}}{
Prior CDF
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MCMCModel$prior_cdf(
  target_treatment_effect,
  n_samples_quantile_estimation = 10000
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{target_treatment_effect}}{Point at which to evaluate the prior CDF}

\item{\code{n_samples_quantile_estimation}}{Number of samples used to estimate the quantiles of the distribution}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MCMCModel-sample_prior"></a>}}
\if{latex}{\out{\hypertarget{method-MCMCModel-sample_prior}{}}}
\subsection{Method \code{sample_prior()}}{
Sample from the prior distribution
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MCMCModel$sample_prior(n_samples)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n_samples}}{Number of samples to draw}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A vector of samples
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MCMCModel-clone"></a>}}
\if{latex}{\out{\hypertarget{method-MCMCModel-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MCMCModel$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
