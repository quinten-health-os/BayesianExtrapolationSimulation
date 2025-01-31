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



\subsection{Rationale of the design}

The core of the code implements the inference logic: a source data class represents source data, which are observed. Target data are represented in a different class. They different as we can sample target data replicates. When initializing a model, source data are provided to it, as well as methods parameters and MCMC configuration (if applicable), so as to define a prior. Inference is performed by using the inference method on target data, which updates properties of the model object with moments of the treatment effect posterior distribution and, if applicable, computes the posterior distribution of borrowing parameters. Note that inference proceeds in two steps: first, if the method uses empirical Bayes, prior parameters are updated based on the target data sample, second Bayes" rule is applied.

When launching a simulation study, all scenarios are generated based on the provided configuration files. These scenarios are then sequentially treated: source data and target data are initiated, and a model is defined. Target data replicates are then generated depending on the target data characteristics. The model is then fitted on each data replicate, and the corresponding inference metrics and frequentist metrics are computed. 
The data from the simulation are then analyzed: the sweet spot is computed for each metric of interest, and the Bayesian operating characteristics are estimated. 

\subsection{Use of existing code and packages}

We tried, as far as possible, to reuse existing methods implementations. Not only did we need to perform inference with the method of interest, but, to compute the prior ESS using the different approaches we included (ELIR, difference between the moment-based or precision-based ESS of the posterior and the target study sample size per arm), we needed to sample from the prior and the posterior distribution.  
We identified several packages and code repositories that could potentially be used. 
The \textit{RBesT} package (\url{https://opensource.nibr.com/RBesT/}) implements inference with conjugate mixture priors and Bayesian meta-analysis. It also allows for approximating distributions using mixture distributions, based on samples, and to compute ELIR and moment-based ESS. We therefore used this package for ESS computation.
For Gaussian endpoints, we used  \textit{RBesT} for Separate Bayesian analysis and Pooling, as well as for the use of Gaussian Robust Mixture Priors.

The \textit{NPP} package (\url{https://cran.r-project.org/web/packages/NPP/index.html}) contains an implementation of the Normalized Power Prior for normally distributed endpoints. However, it uses a custom MCMC implementation, whereas an analytical posterior is available in this case (see Section \ref{subsec:NPP}). Moreover,  the \textit{NPP} package requires individual-level data as input. These two elements prompted us to write a custom implementation for computational efficiency.

The \textit{historicalborrow }package (\url{https://wlandau.github.io/historicalborrow/index.html}) is focused on control group borrowing.  It includes hierarchical models such as the MAC model, as well as pooled and separate analysis. It also implements simulation routines, but with limited flexibility: it would not allow us to simulate time-to-event or recurrent event data. Therefore, we did not use this package.

The \textit{PowerPriorVari }repository (\url{https://github.com/lxt3/PowerPriorVari}), which implements variations of the power prior to borrow from a single source study \parencite{thompson_dynamic_2021} does not contain reusable code.

The \textit{ESS }repository  (\url{https://github.com/DKFZ-biostats/ESS}), contains implementations of different methods to estimate different prior ESS measures for a wide variety of models and is well documented. However, to use a limited number of packages, we only relied on \textit{RBesT} for ESS computation.

\textit{psborrow2 } (\url{https://genentech.github.io/psborrow2/index.html}) is an R package for conducting Bayesian dynamic borrowing analyses and simulation studies. However, it focuses on the borrowing of an external control arm and only implements the hierarchical commensurate prior. It was therefore of limited use for our study.

The \textit{StudyPrior} mostly focuses on the case of binomial likelihood. For Gaussian likelihood, the package implements the Empirical Bayes Power Prior and the Normalized Power Prior. However, for the Empirical Bayes Power Prior, we used the implementation provided in \textcite{nikolakopoulos_dynamic_2018}, and for the Normalized Power Prior, we used a custom implementation of the analytical posterior to avoid relying on computationally expensive approximations.

The \textit{hdbayes} package was released in April 2024 (\url{https://github.com/ethan-alt/hdbayes}), well after the start of the implementation phase of the project. It implements a variety of Bayesian borrowing methods for generalized linear models. However, it always uses Stan for inference, whereas we were able to use analytical posteriors in several cases, which provided significant computational gains. 

