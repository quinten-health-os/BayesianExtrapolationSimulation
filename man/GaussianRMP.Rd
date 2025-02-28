% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/method_RMP.R
\name{GaussianRMP}
\alias{GaussianRMP}
\title{GaussianRMP class}
\description{
A class for Gaussian Robust Mixture Prior models.
}
\details{
This class represents a Gaussian Robust Mixture Prior model.
It inherits from the Model class.

This class extends the Model class and provides methods for
initializing the model, updating priors, calculating posterior moments,
and sampling from prior and posterior distributions.
}
\keyword{internal}
\section{Super class}{
\code{\link[RBExT:Model]{RBExT::Model}} -> \code{GaussianRMP}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{w}}{The weight of the prior distribution.}

\item{\code{vague_prior_mean}}{The mean of the vague prior distribution.}

\item{\code{vague_prior_variance}}{The variance of the vague prior distribution.}

\item{\code{info_prior_mean}}{The mean of the informative prior distribution.}

\item{\code{info_prior_variance}}{The variance of the informative prior distribution.}

\item{\code{wpost}}{The weight of the posterior distribution.}

\item{\code{vague_posterior_mean}}{The mean of the vague posterior distribution.}

\item{\code{info_posterior_mean}}{The mean of the informative posterior distribution.}

\item{\code{vague_posterior_variance}}{The variance of the vague posterior distribution.}

\item{\code{info_posterior_variance}}{The variance of the informative posterior distribution.}

\item{\code{method}}{Name of the method}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-GaussianRMP-new}{\code{GaussianRMP$new()}}
\item \href{#method-GaussianRMP-empirical_bayes_update}{\code{GaussianRMP$empirical_bayes_update()}}
\item \href{#method-GaussianRMP-prior_weight}{\code{GaussianRMP$prior_weight()}}
\item \href{#method-GaussianRMP-prior_pdf}{\code{GaussianRMP$prior_pdf()}}
\item \href{#method-GaussianRMP-prior_cdf}{\code{GaussianRMP$prior_cdf()}}
\item \href{#method-GaussianRMP-posterior_moments}{\code{GaussianRMP$posterior_moments()}}
\item \href{#method-GaussianRMP-posterior_pdf}{\code{GaussianRMP$posterior_pdf()}}
\item \href{#method-GaussianRMP-posterior_cdf}{\code{GaussianRMP$posterior_cdf()}}
\item \href{#method-GaussianRMP-sample_prior}{\code{GaussianRMP$sample_prior()}}
\item \href{#method-GaussianRMP-sample_posterior}{\code{GaussianRMP$sample_posterior()}}
\item \href{#method-GaussianRMP-posterior_mean}{\code{GaussianRMP$posterior_mean()}}
\item \href{#method-GaussianRMP-posterior_variance}{\code{GaussianRMP$posterior_variance()}}
\item \href{#method-GaussianRMP-prior_to_RBesT}{\code{GaussianRMP$prior_to_RBesT()}}
\item \href{#method-GaussianRMP-posterior_to_RBesT}{\code{GaussianRMP$posterior_to_RBesT()}}
\item \href{#method-GaussianRMP-print_model_summary}{\code{GaussianRMP$print_model_summary()}}
\item \href{#method-GaussianRMP-clone}{\code{GaussianRMP$clone()}}
}
}
\if{html}{\out{
<details><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="RBExT" data-topic="Model" data-id="create"><a href='../../RBExT/html/Model.html#method-Model-create'><code>RBExT::Model$create()</code></a></span></li>
<li><span class="pkg-link" data-pkg="RBExT" data-topic="Model" data-id="credible_interval"><a href='../../RBExT/html/Model.html#method-Model-credible_interval'><code>RBExT::Model$credible_interval()</code></a></span></li>
<li><span class="pkg-link" data-pkg="RBExT" data-topic="Model" data-id="estimate_bayesian_operating_characteristics"><a href='../../RBExT/html/Model.html#method-Model-estimate_bayesian_operating_characteristics'><code>RBExT::Model$estimate_bayesian_operating_characteristics()</code></a></span></li>
<li><span class="pkg-link" data-pkg="RBExT" data-topic="Model" data-id="estimate_frequentist_operating_characteristics"><a href='../../RBExT/html/Model.html#method-Model-estimate_frequentist_operating_characteristics'><code>RBExT::Model$estimate_frequentist_operating_characteristics()</code></a></span></li>
<li><span class="pkg-link" data-pkg="RBExT" data-topic="Model" data-id="inference"><a href='../../RBExT/html/Model.html#method-Model-inference'><code>RBExT::Model$inference()</code></a></span></li>
<li><span class="pkg-link" data-pkg="RBExT" data-topic="Model" data-id="plot_pdfs"><a href='../../RBExT/html/Model.html#method-Model-plot_pdfs'><code>RBExT::Model$plot_pdfs()</code></a></span></li>
<li><span class="pkg-link" data-pkg="RBExT" data-topic="Model" data-id="plot_posterior_pdf"><a href='../../RBExT/html/Model.html#method-Model-plot_posterior_pdf'><code>RBExT::Model$plot_posterior_pdf()</code></a></span></li>
<li><span class="pkg-link" data-pkg="RBExT" data-topic="Model" data-id="plot_prior_pdf"><a href='../../RBExT/html/Model.html#method-Model-plot_prior_pdf'><code>RBExT::Model$plot_prior_pdf()</code></a></span></li>
<li><span class="pkg-link" data-pkg="RBExT" data-topic="Model" data-id="posterior_median"><a href='../../RBExT/html/Model.html#method-Model-posterior_median'><code>RBExT::Model$posterior_median()</code></a></span></li>
<li><span class="pkg-link" data-pkg="RBExT" data-topic="Model" data-id="posterior_quantile"><a href='../../RBExT/html/Model.html#method-Model-posterior_quantile'><code>RBExT::Model$posterior_quantile()</code></a></span></li>
<li><span class="pkg-link" data-pkg="RBExT" data-topic="Model" data-id="prior_ESS"><a href='../../RBExT/html/Model.html#method-Model-prior_ESS'><code>RBExT::Model$prior_ESS()</code></a></span></li>
<li><span class="pkg-link" data-pkg="RBExT" data-topic="Model" data-id="prior_treatment_benefit"><a href='../../RBExT/html/Model.html#method-Model-prior_treatment_benefit'><code>RBExT::Model$prior_treatment_benefit()</code></a></span></li>
<li><span class="pkg-link" data-pkg="RBExT" data-topic="Model" data-id="simulation_for_given_treatment_effect"><a href='../../RBExT/html/Model.html#method-Model-simulation_for_given_treatment_effect'><code>RBExT::Model$simulation_for_given_treatment_effect()</code></a></span></li>
<li><span class="pkg-link" data-pkg="RBExT" data-topic="Model" data-id="test_decision"><a href='../../RBExT/html/Model.html#method-Model-test_decision'><code>RBExT::Model$test_decision()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GaussianRMP-new"></a>}}
\if{latex}{\out{\hypertarget{method-GaussianRMP-new}{}}}
\subsection{Method \code{new()}}{
Initialize a new GaussianRMP object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GaussianRMP$new(prior)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{prior}}{A list containing prior information for the analysis.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GaussianRMP-empirical_bayes_update"></a>}}
\if{latex}{\out{\hypertarget{method-GaussianRMP-empirical_bayes_update}{}}}
\subsection{Method \code{empirical_bayes_update()}}{
Update the vague prior variance based on empirical Bayes approach.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GaussianRMP$empirical_bayes_update(target_data)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{target_data}}{The target data for the analysis.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GaussianRMP-prior_weight"></a>}}
\if{latex}{\out{\hypertarget{method-GaussianRMP-prior_weight}{}}}
\subsection{Method \code{prior_weight()}}{
Calculate the posterior weight based on the target data.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GaussianRMP$prior_weight(target_data)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{target_data}}{The target data for the analysis.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The posterior weight.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GaussianRMP-prior_pdf"></a>}}
\if{latex}{\out{\hypertarget{method-GaussianRMP-prior_pdf}{}}}
\subsection{Method \code{prior_pdf()}}{
Calculate the prior probability density function (PDF) for a given target treatment effect.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GaussianRMP$prior_pdf(target_treatment_effect)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{target_treatment_effect}}{The target treatment effect.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The prior PDF.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GaussianRMP-prior_cdf"></a>}}
\if{latex}{\out{\hypertarget{method-GaussianRMP-prior_cdf}{}}}
\subsection{Method \code{prior_cdf()}}{
Calculate the prior cumulative distribution function (CDF) for a given target treatment effect.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GaussianRMP$prior_cdf(target_treatment_effect)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{target_treatment_effect}}{The target treatment effect.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The prior CDF.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GaussianRMP-posterior_moments"></a>}}
\if{latex}{\out{\hypertarget{method-GaussianRMP-posterior_moments}{}}}
\subsection{Method \code{posterior_moments()}}{
Calculate the posterior moments based on the target data.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GaussianRMP$posterior_moments(target_data)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{target_data}}{A list containing the target data for the analysis.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GaussianRMP-posterior_pdf"></a>}}
\if{latex}{\out{\hypertarget{method-GaussianRMP-posterior_pdf}{}}}
\subsection{Method \code{posterior_pdf()}}{
Calculate the posterior probability density function (PDF) for a given target treatment effect.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GaussianRMP$posterior_pdf(target_treatment_effect)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{target_treatment_effect}}{The target treatment effect.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The posterior PDF.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GaussianRMP-posterior_cdf"></a>}}
\if{latex}{\out{\hypertarget{method-GaussianRMP-posterior_cdf}{}}}
\subsection{Method \code{posterior_cdf()}}{
Calculate the posterior cumulative distribution function (CDF) for a given target treatment effect.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GaussianRMP$posterior_cdf(target_treatment_effect)}\if{html}{\out{</div>}}
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
\if{html}{\out{<a id="method-GaussianRMP-sample_prior"></a>}}
\if{latex}{\out{\hypertarget{method-GaussianRMP-sample_prior}{}}}
\subsection{Method \code{sample_prior()}}{
Sample from the prior distribution.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GaussianRMP$sample_prior(n_samples)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n_samples}}{The number of samples to generate.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The samples from the prior distribution.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GaussianRMP-sample_posterior"></a>}}
\if{latex}{\out{\hypertarget{method-GaussianRMP-sample_posterior}{}}}
\subsection{Method \code{sample_posterior()}}{
Sample from the posterior distribution.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GaussianRMP$sample_posterior(n_samples)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n_samples}}{The number of samples to generate.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The samples from the posterior distribution.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GaussianRMP-posterior_mean"></a>}}
\if{latex}{\out{\hypertarget{method-GaussianRMP-posterior_mean}{}}}
\subsection{Method \code{posterior_mean()}}{
Calculate the posterior mean.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GaussianRMP$posterior_mean()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
The posterior mean.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GaussianRMP-posterior_variance"></a>}}
\if{latex}{\out{\hypertarget{method-GaussianRMP-posterior_variance}{}}}
\subsection{Method \code{posterior_variance()}}{
Calculate the posterior variance.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GaussianRMP$posterior_variance()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
The posterior variance.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GaussianRMP-prior_to_RBesT"></a>}}
\if{latex}{\out{\hypertarget{method-GaussianRMP-prior_to_RBesT}{}}}
\subsection{Method \code{prior_to_RBesT()}}{
Convert the prior distribution to the RBesT format.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GaussianRMP$prior_to_RBesT(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{Additional arguments.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GaussianRMP-posterior_to_RBesT"></a>}}
\if{latex}{\out{\hypertarget{method-GaussianRMP-posterior_to_RBesT}{}}}
\subsection{Method \code{posterior_to_RBesT()}}{
Convert the posterior distribution to the RBesT format.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GaussianRMP$posterior_to_RBesT(target_data, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{target_data}}{Target data for the analysis.}

\item{\code{...}}{Additional arguments.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GaussianRMP-print_model_summary"></a>}}
\if{latex}{\out{\hypertarget{method-GaussianRMP-print_model_summary}{}}}
\subsection{Method \code{print_model_summary()}}{
Print a summary of the model attributes

This method creates and prints a formatted table of key model attributes.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GaussianRMP$print_model_summary()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A printed data frame displaying the following model attributes:
\itemize{
\item Posterior Weight
\item Vague Posterior Mean
\item Vague Posterior Variance
\item Informative Posterior Mean
\item Informative Posterior Variance
}
All numeric values are formatted to 6 decimal places.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GaussianRMP-clone"></a>}}
\if{latex}{\out{\hypertarget{method-GaussianRMP-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GaussianRMP$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
