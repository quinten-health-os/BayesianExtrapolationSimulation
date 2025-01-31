library(testthat)

env <- "pipeline_tests"
config_dir <- paste0("../../inst/conf/", env, "/")
case_studies_config_dir <- "../../inst/conf/case_studies/"
results_dir <- paste0("../../results/", env, "/")
outputs_config <- yaml::read_yaml(system.file("conf/outputs_config.yml", package = "RBExT"))
ocs_filename_frequentist <- outputs_config$frequentist_ocs_results_filename
ocs_filename_bayesian <- outputs_config$bayesian_ocs_deterministic_results_filename

case_studies <- c(
  "belimumab",
  "botox",
  "dapagliflozin",
  "mepolizumab",
  "teriflunomide",
  "aprepitant"
)
methods <- c(
  "separate",
  "pooling",
  "conditional_power_prior",
  "RMP",
  "NPP",
  "PDCCPP",
  "predictive_test_then_pool",
  "EB_PP",
  "test_then_pool_equivalence",
  "test_then_pool_difference",
  "p_value_based_PP",
  "ElasticPrior",
  "commensurate_power_prior"
)

test_results_format <- function(df, config) {
  for (col_name in names(config)) {
    col_config <- config[[col_name]]

    test_that(paste("Testing column", col_name), {
      # Test if column in dataframe
      expect_true(col_name %in% colnames(df))

      non_na_values <- df[[col_name]][!is.na(df[[col_name]])]
      if (length(non_na_values) > 0) {
        # Test type
        if (!is.null(col_config$type)) {
          expect_true(
            all(sapply(non_na_values, class) == col_config$type),
            info = paste(
              "Column",
              col_name,
              "should be of type",
              col_config$type
            )
          )
        }

        # Test accepted values for character type
        if (!is.null(col_config$type) &&
          col_config$type == "character" &&
          !is.null(col_config$accepted_values)) {
          expect_true(
            all(non_na_values %in% col_config$accepted_values),
            info = paste(
              "Column",
              col_name,
              "contains invalid values. Accepted values are:",
              paste(col_config$accepted_values, collapse = ", ")
            )
          )
        }

        # Test range
        if (!is.null(col_config$range)) {
          expect_true(
            all(
              non_na_values >= col_config$range[1] &
                non_na_values <= col_config$range[2]
            ),
            info = paste(
              "Values in column",
              col_name,
              "should be in range",
              col_config$range
            )
          )
        }
      }

      # Test absence of missing values
      if (!is.null(col_config$not_missing) &&
        col_config$not_missing) {
        expect_true(all(!is.na(df[[col_name]])),
          info = paste("Column", col_name, "should not have missing values")
        )
      }

      # Test value comparison between columns
      if (!is.null(col_config$compare)) {
        compare_col <- col_config$compare$column
        operator <- col_config$compare$operator
        df_compare <- eval(parse(
          text = paste("df[[col_name]]", operator, "df[[compare_col]]")
        ))
        expect_true(
          all(df_compare),
          info = paste(
            "Values in column",
            col_name,
            "should be",
            operator,
            "values in column",
            compare_col
          )
        )
      }
    })
  }
}

df_frequentist <- readr::read_csv(paste0(results_dir, "/", ocs_filename_frequentist))
config_frequentist <- list(
  case_study = list(
    type = "character",
    accepted_values = case_studies,
    not_missing = TRUE
  ),
  method = list(
    type = "character",
    accepted_values = methods,
    not_missing = TRUE
  ),
  parameters = list(type = "character", not_missing = TRUE),
  drift = list(
    type = "numeric",
    range = c(-Inf, Inf),
    not_missing = TRUE
  ),
  treatment_drift = list(
    type = "numeric",
    range = c(-Inf, Inf),
    not_missing = TRUE
  ),
  control_drift = list(
    type = "numeric",
    range = c(-Inf, Inf),
    not_missing = TRUE
  ),
  source_denominator = list(
    type = "numeric",
    range = c(0, Inf),
    not_missing = FALSE
  ),
  source_denominator_change_factor = list(
    type = "numeric",
    range = c(0, Inf),
    not_missing = FALSE
  ),
  parallelization = list(type = "logical", not_missing = TRUE),
  success_proba = list(
    type = "numeric",
    range = c(0, 1),
    not_missing = TRUE
  ),
  mcse_success_proba = list(
    type = "numeric",
    range = c(0, 1),
    not_missing = TRUE
  ),
  conf_int_success_proba_lower = list(
    type = "numeric",
    range = c(0, 1),
    not_missing = TRUE,
    compare = list(column = "conf_int_success_proba_upper", operator = "<=")
  ),
  conf_int_success_proba_upper = list(
    type = "numeric",
    range = c(0, 1),
    not_missing = TRUE,
    compare = list(column = "conf_int_success_proba_lower", operator = ">=")
  ),
  ess_moment = list(
    type = "numeric",
    range = c(-Inf, Inf),
    not_missing = TRUE
  ),
  ess_precision = list(
    type = "numeric",
    range = c(-Inf, Inf),
    not_missing = TRUE
  ),
  source_treatment_effect_estimate = list(
    type = "numeric",
    range = c(-Inf, Inf),
    not_missing = TRUE
  ),
  source_total_sample_size = list(
    type = "numeric",
    range = c(0, Inf),
    not_missing = TRUE
  ),
  source_standard_error = list(
    type = "numeric",
    range = c(0, Inf),
    not_missing = TRUE
  ),
  endpoint = list(
    type = "character",
    accepted_values = c("continuous", "binary", "recurrent_event", "time_to_event"),
    not_missing = TRUE
  ),
  source_sample_size_control = list(
    type = "numeric",
    range = c(0, Inf),
    not_missing = TRUE
  ),
  source_sample_size_treatment = list(
    type = "numeric",
    range = c(0, Inf),
    not_missing = TRUE
  ),
  summary_measure_likelihood = list(
    type = "character",
    accepted_values = c("Normal", "binomial"),
    not_missing = TRUE
  ),
  target_sample_size_per_arm = list(
    type = "numeric",
    range = c(0, Inf),
    not_missing = TRUE
  ),
  target_treatment_effect = list(
    type = "numeric",
    range = c(-Inf, Inf),
    not_missing = TRUE
  ),
  #   target_standard_deviation = list(type = "numeric", range = c(0, Inf), not_missing = TRUE),
  sampling_approximation = list(type = "logical", not_missing = TRUE),
  theta_0 = list(
    type = "numeric",
    range = c(0, 1),
    not_missing = TRUE
  ),
  null_space = list(
    type = "character",
    accepted_values = c("left", "right"),
    not_missing = TRUE
  ),
  coverage = list(
    type = "numeric",
    range = c(0, 1),
    not_missing = TRUE
  ),
  conf_int_coverage_lower = list(
    type = "numeric",
    range = c(0, 1),
    not_missing = TRUE,
    compare = list(column = "conf_int_coverage_upper", operator = "<=")
  ),
  conf_int_coverage_upper = list(
    type = "numeric",
    range = c(0, 1),
    not_missing = TRUE,
    compare = list(column = "conf_int_coverage_lower", operator = ">=")
  ),
  mse = list(
    type = "numeric",
    range = c(0, Inf),
    not_missing = TRUE
  ),
  conf_int_mse_lower = list(
    type = "numeric",
    range = c(0, Inf),
    not_missing = TRUE,
    compare = list(column = "conf_int_mse_upper", operator = "<=")
  ),
  conf_int_mse_upper = list(
    type = "numeric",
    range = c(0, Inf),
    not_missing = TRUE,
    compare = list(column = "conf_int_mse_lower", operator = ">=")
  ),
  bias = list(
    type = "numeric",
    range = c(-Inf, Inf),
    not_missing = TRUE
  ),
  conf_int_bias_lower = list(
    type = "numeric",
    range = c(-Inf, Inf),
    not_missing = TRUE,
    compare = list(column = "conf_int_bias_upper", operator = "<=")
  ),
  conf_int_bias_upper = list(
    type = "numeric",
    range = c(-Inf, Inf),
    not_missing = TRUE,
    compare = list(column = "conf_int_bias_lower", operator = ">=")
  ),
  posterior_mean = list(
    type = "numeric",
    range = c(-Inf, Inf),
    not_missing = TRUE
  ),
  conf_int_posterior_mean_lower = list(
    type = "numeric",
    range = c(-Inf, Inf),
    not_missing = TRUE,
    compare = list(column = "conf_int_posterior_mean_upper", operator = "<=")
  ),
  conf_int_posterior_mean_upper = list(
    type = "numeric",
    range = c(-Inf, Inf),
    not_missing = TRUE,
    compare = list(column = "conf_int_posterior_mean_lower", operator = ">=")
  ),
  posterior_median = list(
    type = "numeric",
    range = c(-Inf, Inf),
    not_missing = TRUE
  ),
  conf_int_posterior_median_lower = list(
    type = "numeric",
    range = c(-Inf, Inf),
    not_missing = TRUE,
    compare = list(column = "conf_int_posterior_median_upper", operator = "<=")
  ),
  conf_int_posterior_median_upper = list(
    type = "numeric",
    range = c(-Inf, Inf),
    not_missing = TRUE,
    compare = list(column = "conf_int_posterior_median_lower", operator = ">=")
  ),
  precision = list(
    type = "numeric",
    range = c(0, Inf),
    not_missing = TRUE
  ),
  conf_int_precision_lower = list(
    type = "numeric",
    range = c(0, Inf),
    not_missing = TRUE,
    compare = list(column = "conf_int_precision_upper", operator = "<=")
  ),
  conf_int_precision_upper = list(
    type = "numeric",
    range = c(0, Inf),
    not_missing = TRUE,
    compare = list(column = "conf_int_precision_lower", operator = ">=")
  ),
  credible_interval_lower = list(
    type = "numeric",
    range = c(-Inf, Inf),
    not_missing = TRUE,
    compare = list(column = "credible_interval_upper", operator = "<=")
  ),
  credible_interval_upper = list(
    type = "numeric",
    range = c(-Inf, Inf),
    not_missing = TRUE,
    compare = list(column = "credible_interval_lower", operator = ">=")
  ),
  posterior_parameters = list(type = "character", not_missing = TRUE),
  rng_state = list(type = "numeric", not_missing = TRUE),
  computation_time = list(
    type = "numeric",
    range = c(0, Inf),
    not_missing = TRUE
  ),
  #   target_to_source_std_ratio = list(type = "numeric", range = c(0, Inf), not_missing = TRUE),
  tie = list(
    type = "numeric",
    range = c(0, 1),
    not_missing = TRUE
  ),
  mcse_tie = list(
    type = "numeric",
    range = c(0, Inf),
    not_missing = TRUE
  ),
  conf_int_tie_lower = list(
    type = "numeric",
    range = c(0, 1),
    not_missing = TRUE,
    compare = list(column = "conf_int_tie_upper", operator = "<=")
  ),
  conf_int_tie_upper = list(
    type = "numeric",
    range = c(0, 1),
    not_missing = TRUE,
    compare = list(column = "conf_int_tie_lower", operator = ">=")
  ),
  #   power = list(type = "numeric", range = c(0, 1), not_missing = TRUE),
  frequentist_test = list(
    type = "character",
    accepted_values = c(analysis_config[["frequentist_test"]]),
    not_missing = TRUE
  )
)

df_bayesian <- readr::read_csv(paste0(results_dir, "/", ocs_filename_bayesian))
config_bayesian <- list(
  case_study = list(
    type = "character",
    accepted_values = case_studies,
    not_missing = TRUE
  ),
  method = list(
    type = "character",
    accepted_values = methods,
    not_missing = TRUE
  ),
  parameters = list(type = "character", not_missing = TRUE),
  target_sample_size_per_arm = list(
    type = "numeric",
    range = c(0, Inf),
    not_missing = TRUE
  ),
  source_denominator = list(
    type = "numeric",
    range = c(0, Inf),
    not_missing = FALSE
  ),
  source_denominator_change_factor = list(
    type = "numeric",
    range = c(0, Inf),
    not_missing = FALSE
  ),
  parallelization = list(type = "logical", not_missing = TRUE),
  sampling_approximation = list(type = "logical", not_missing = TRUE),
  design_prior = list(
    type = "character",
    accepted_values = c("ui_design_prior", "analysis_prior", "source_posterior"),
    not_missing = TRUE
  ),
  #   average_tie = list(type = "numeric", range = c(0, 1), not_missing = FALSE),
  average_power = list(
    type = "numeric",
    range = c(0, Inf),
    not_missing = FALSE
  ),
  prior_proba_no_benefit = list(
    type = "numeric",
    range = c(0, 1),
    not_missing = FALSE
  ),
  prepost_proba_FP = list(
    type = "numeric",
    range = c(0, 1),
    not_missing = FALSE
  ),
  prepost_proba_TP = list(
    type = "numeric",
    range = c(0, 1),
    not_missing = FALSE
  ),
  upper_bound_proba_FP = list(
    type = "numeric",
    range = c(0, 1),
    not_missing = FALSE
  ),
  prior_proba_success = list(
    type = "numeric",
    range = c(0, 1),
    not_missing = FALSE
  ),
  source_treatment_effect_estimate = list(
    type = "numeric",
    range = c(-Inf, Inf),
    not_missing = TRUE
  ),
  source_standard_error = list(
    type = "numeric",
    range = c(0, Inf),
    not_missing = TRUE
  ),
  endpoint = list(
    type = "character",
    accepted_values = c("continuous", "binary", "recurrent_event", "time_to_event"),
    not_missing = TRUE
  ),
  summary_measure_likelihood = list(
    type = "character",
    accepted_values = c("Normal", "binomial"),
    not_missing = TRUE
  ),
  source_sample_size_control = list(
    type = "numeric",
    range = c(0, Inf),
    not_missing = TRUE
  ),
  source_sample_size_treatment = list(
    type = "numeric",
    range = c(0, Inf),
    not_missing = TRUE
  ),
  rng_state = list(type = "numeric", not_missing = TRUE),
  computation_time = list(
    type = "numeric",
    range = c(0, Inf),
    not_missing = TRUE
  )
)

test_results_format(df_frequentist, config_frequentist)
test_results_format(df_bayesian, config_bayesian)
