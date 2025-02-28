% Generated by roxygen2: do not edit by hand
% Please edit documentation in
%   R/table_empirical_bayes_hyperparameters_vs_scenario.R
\name{table_empirical_bayes_hyperparameters_vs_drift}
\alias{table_empirical_bayes_hyperparameters_vs_drift}
\title{Generate a table of hyperparameters estimated using Empirical Bayes vs drift}
\usage{
table_empirical_bayes_hyperparameters_vs_drift(
  results_metrics_df,
  method,
  case_study,
  control_drift,
  xvars,
  source_denominator_change_factor,
  target_to_source_std_ratio,
  parameters_combinations
)
}
\arguments{
\item{results_metrics_df}{A data frame containing the results metrics.}

\item{method}{The method to filter the results by.}

\item{case_study}{The case study to filter the results by.}

\item{control_drift}{A logical value indicating whether to control for drift.}

\item{xvars}{A list containing the x-variable configurations.}

\item{target_sample_size_per_arm}{The target sample size per arm to filter the results by.}
}
\value{
This function does not return a value. It creates and saves tables in HTML, PDF, and LaTeX formats.
}
\description{
Generate a table of hyperparameters estimated using Empirical Bayes vs drift
}
