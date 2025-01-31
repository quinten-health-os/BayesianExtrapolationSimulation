# R Bayesian Extrapolation Tool

A collection of R tools to :

- Study frequentist and Bayesian operating characteristics of clinical trial designs leveraging Bayesian partial extrapolation (also known as borrowing).
- Analyze clinical trial data using Bayesian partial extrapolation methods.

Copyright 2024 Quinten Health, under exclusive licence to the European Medicines Agency. Sharing and distribution are prohibited.

- Core contributor: Tristan Fauvel
- Contributors: Pascal Godbillot
- Project management: Marie Génin
- Leads: Julien Tanniou, Billy Amzal
- Sponsors: European Medicines Agency, Quinten Health

## Quick start

### Installation

Download the latest release.

You can then install the RBExT package by running :

```
install.packages("/path/to/RBExT_0.0.1.tar.gz", repos = NULL, type = "source")
```

and load it using:

```
library(RBExT)
```

### Running and analyzing simulations

- To launch a simulation study, run main.R

- To plot the results of the simulation study, run plots.R

- To create the results tables, run tables.R

## Design logic

There are three main objects in the simulation framework: source data, target data, and models.
Source data are defined based on the case study configuration. Target data depend on the source data, and the specific scenario (drift, sample size). Models correspond to specific methods, and depend on the source data. They implement methods to estimate their operating characteristics. The target data objects implement a sample method, allowing to sample many replicates of the target study.

## Access documentation

To access documentation :
browseURL("docs/index.html")

## Methods and configurations requiring MCMC inference

MCMC inference with Stan is used in the following scenarios :

- RMP x aprepitant
- separate x aprepitant
- pooling x aprepitant
- conditional_power_prior x aprepitant

When using MCMC inference, the method inherits from the MCMCModel class. The initialize() method contains a string with the Stan code. This code is saved in a temporary Stan file which is used to compile the model (using `self$stan_model <- cmdstanr::cmdstan_model(stan_file_path)`). At inference

### Comparison with existing packages

- psborrow2 implements Bayesian dynamic borrowing in scenarios where a two-arm RCT is supplemented with external data on the control arm.

- RBEsT focuses on meta-analytic and mixture models

## Contributing

New contributors are always welcome. Please have a look at the [contribution guidelines](CONTRIBUTING.md) on how to get started and to make sure your code complies with our guidelines.

## Rationale of the design

The core of the code implements the inference logic: a source data class represents source data, which are observed. Target data are represented in a different class. They different as we can sample target data replicates. When initializing a model, source data are provided to it, as well as methods parameters and MCMC configuration (if applicable), so as to define a prior. Inference is performed by using the inference method on target data, which updates properties of the model object with moments of the treatment effect posterior distribution and, if applicable, computes the posterior distribution of borrowing parameters. Note that inference proceeds in two steps: first, if the method uses empirical Bayes, prior parameters are updated based on the target data sample, second Bayes" rule is applied.

When launching a simulation study, all scenarios are generated based on the provided configuration files. These scenarios are then sequentially treated: source data and target data are initiated, and a model is defined. Target data replicates are then generated depending on the target data characteristics. The model is then fitted on each data replicate, and the corresponding inference metrics and frequentist metrics are computed.
The data from the simulation are then analyzed: the sweet spot is computed for each metric of interest, and the Bayesian operating characteristics are estimated.
