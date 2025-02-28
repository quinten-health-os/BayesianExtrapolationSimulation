% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_bayesian_ocs.R
\name{estimate_bayesian_ocs}
\alias{estimate_bayesian_ocs}
\title{Simulation Bayesian OCs}
\usage{
estimate_bayesian_ocs(
  scenario,
  case_study_config,
  source_data,
  target_sample_size_per_arm,
  model,
  method_parameters,
  theta_0,
  n_replicates,
  critical_value,
  computation_state,
  design_prior_type,
  json_parameters,
  simulation_config,
  target_to_source_std_ratio,
  mcmc_config = NULL
)
}
\arguments{
\item{scenario}{Simulation scenario}

\item{case_study_config}{The case study configuration.}

\item{source_data}{The source data.}

\item{target_sample_size_per_arm}{The target sample size per arm.}

\item{model}{The model.}

\item{method_parameters}{The method parameters.}

\item{theta_0}{The true parameter value.}

\item{n_replicates}{The number of replicates.}

\item{critical_value}{The critical value.}

\item{computation_state}{The computation state.}

\item{design_prior_type}{The design prior type.}

\item{json_parameters}{Method parameters in json format}

\item{simulation_config}{Simulation configuration}

\item{target_to_source_std_ratio}{Ratio between the source and target study sampling standard deviation.}

\item{mcmc_config}{MCMC configuration}
}
\value{
The results of the Bayesian operating characteristics estimation.
}
\description{
This function estimates the Bayesian operating characteristics for the model in the scenario considered.
}
