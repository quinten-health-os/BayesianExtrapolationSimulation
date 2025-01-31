# Create a nested list with the plot settings
frequentist_metrics <<- list(
  success_proba = list(
    name = "success_proba",
    metric_uncertainty = "conf_int_success_proba",
    label = "Pr(Study success)",
    uncertainty_label = "95% CI", # "MCSE",
    baseline = FALSE,
    larger_is_better = TRUE
  ),
  mse = list(
    name = "mse",
    metric_uncertainty = "conf_int_mse",
    label = "MSE",
    uncertainty_label = "95% CI", # "95% CI",
    baseline = 0,
    larger_is_better = FALSE
  ),
  bias = list(
    name = "bias",
    metric_uncertainty = "conf_int_bias",
    label = "Bias",
    uncertainty_label = "95% CI",
    baseline = 0,
    larger_is_better = FALSE
  ),
  coverage = list(
    name = "coverage",
    metric_uncertainty = "conf_int_coverage",
    label = "Coverage of the 95% CrI",
    uncertainty_label = "95% CI",
    baseline = FALSE,
    larger_is_better = TRUE
  ),
  precision = list(
    name = "precision",
    metric_uncertainty = "conf_int_precision",
    label = "Half width of the 95% CrI",
    uncertainty_label = "95% CI",
    baseline = 0,
    larger_is_better = FALSE
  ),
  tie = list(
    name = "tie",
    metric_uncertainty = "conf_int_tie",
    label = "TIE",
    uncertainty_label = "95% CI",
    baseline = 0.05,
    larger_is_better = FALSE
  )
)

inference_metrics <<- list(
  posterior_mean = list(
    name = "posterior_mean",
    metric_uncertainty = "conf_int_posterior_mean",
    label = "Posterior mean",
    uncertainty_label = "95% CI",
    baseline = "theta_0"
  ),
  posterior_median = list(
    name = "posterior_median",
    metric_uncertainty = "conf_int_posterior_median",
    label = "Posterior median",
    uncertainty_label = "95% CI", # "95% CI",
    baseline = "theta_0"
  ),
  point_estimate_mean = list(
    name = "posterior_mean",
    metric_uncertainty = "credible_interval",
    label = "Posterior mean",
    uncertainty_label = "95% CrI",
    baseline = "theta_0"
  ),
  ess_precision = list(
    name = "ess_precision",
    metric_uncertainty = "conf_int_ess_precision",
    label = "Precision based ESS",
    uncertainty_label = "95% CI",
    baseline = FALSE
  ),
  ess_moment = list(
    name = "ess_moment",
    metric_uncertainty = "conf_int_ess_moment",
    label = "Moment based ESS",
    uncertainty_label = "95% CI",
    baseline = FALSE
  ),
  ess_elir = list(
    name = "ess_elir",
    metric_uncertainty = "conf_int_ess_elir",
    label = "ELIR",
    uncertainty_label = "95% CI",
    baseline = FALSE
  )
)

bayesian_metrics <<- list(
  average_tie = list(
    name = "average_tie",
    metric_uncertainty = "conf_int_average_tie",
    label = "Average TIE",
    uncertainty_label = "95% CI",
    baseline = FALSE
  ),
  average_power = list(
    name = "average_power",
    metric_uncertainty = "conf_int_average_power",
    label = "Average power",
    uncertainty_label = "95% CI",
    baseline = FALSE
  ),
  prior_proba_no_benefit = list(
    name = "prior_proba_no_benefit",
    metric_uncertainty = "credible_interval_prior_proba_no_benefit",
    label = "Pr(No benefit)",
    uncertainty_label = "95% CI",
    baseline = FALSE
  ),
  prior_proba_benefit = list(
    name = "prior_proba_benefit",
    metric_uncertainty = "credible_interval_prior_proba_benefit",
    label = "Pr(Benefit)",
    uncertainty_label = "95% CI",
    baseline = FALSE
  ),
  prepost_proba_FP = list(
    name = "prepost_proba_FP",
    metric_uncertainty = "conf_int_prepost_proba_FP",
    label = "Pre-posterior proba. of FP",
    uncertainty_label = "95% CI",
    baseline = FALSE
  ),
  prepost_proba_TP = list(
    name = "prepost_proba_TP",
    metric_uncertainty = "conf_int_prepost_proba_TP",
    label = "Pre-posterior proba. of TP",
    uncertainty_label = "95% CI",
    baseline = FALSE
  ),
  upper_bound_proba_FP = list(
    name = "upper_bound_proba_FP",
    metric_uncertainty = "conf_int_upper_bound_proba_FP",
    label = "Upper bound proba. of FP",
    uncertainty_label = "95% CI",
    baseline = FALSE
  ),
  prior_proba_success = list(
    name = "prior_proba_success",
    metric_uncertainty = "credible_interval_prior_proba_success",
    label = "Prior Pr(Study success)",
    uncertainty_label = "95% CI",
    baseline = FALSE
  )
)


# Create a nested list with the plot settings
sweet_spots_metrics <<- list(
  success_proba = list(
    name = "success_proba",
    label = "Sweet spot width - Pr(Study success)"
  ),
  success_proba_smaller_than_nominal_TIE = list(
    name = "success_proba_smaller_than_nominal_TIE",
    label = "Sweet spot width - Pr(Study success) < nominal TIE"
  ),
  power_larger_than_nominal = list(
    name = "power_larger_than_nominal",
    label = "Sweet spot width - Pr(Study success) > nominal"
  ),
  mse = list(
    name = "mse",
    label = "Sweet spot width - MSE"
  ),
  bias = list(
    name = "bias",
    label = "Sweet spot width - Bias"
  ),
  coverage = list(
    name = "coverage",
    label = "Sweet spot width - coverage of the 95% CrI"
  ),
  precision = list(
    name = "precision",
    label = "Sweet spot width - half width of the 95% CrI"
  ),
  tie = list(
    name = "tie",
    label = "Sweet spot width - TIE"
  )
)
