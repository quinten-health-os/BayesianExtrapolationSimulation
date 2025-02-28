% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_data.R
\name{sample_aggregate_binary_data}
\alias{sample_aggregate_binary_data}
\title{Sample aggregate binary data based on the rate and sample size}
\usage{
sample_aggregate_binary_data(rate, n, n_replicates)
}
\arguments{
\item{rate}{The rate of success (1s) in the generated data}

\item{n}{The total sample size}

\item{n_replicates}{Number of replicates to sample}
}
\value{
A vector of sampled rates
}
\description{
This function samples aggregate binary data based on the rate and sample size.
}
