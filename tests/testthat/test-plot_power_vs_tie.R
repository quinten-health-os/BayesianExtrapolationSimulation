# # Load the CSV file for testing

# env <- "pipeline_tests"
# config_dir <- paste0("../../inst/conf/", env, "/")
# case_studies_config_dir <- "../../inst/conf/case_studies/"
# results_dir <- paste0("../../results/", env, "/")
# outputs_config <- yaml::read_yaml(system.file("conf/outputs_config.yml", package = "RBExT"))
# ocs_filename_frequentist <- outputs_config$csv_filenames[1]
# test_results_metrics_df <- read_csv(paste0(results_dir, "/", ocs_filename_frequentist))

# # Unit tests for power_vs_tie function
# test_that("power_vs_tie function works correctly", {
#   case_study <- "belimumab"
#   target_sample_size_per_arm <- 50
#   treatment_effect <- "strong"
#   power_difference <- FALSE
#   metrics <- c("success_proba", "conf_int_success_proba_lower", "conf_int_success_proba_upper")

#   # Test for strong treatment effect
#   expect_silent(power_vs_tie(test_results_metrics_df, case_study, target_sample_size_per_arm, treatment_effect, power_difference, metrics))

#   # Test for moderate treatment effect
#   treatment_effect <- "moderate"
#   expect_silent(power_vs_tie(test_results_metrics_df, case_study, target_sample_size_per_arm, treatment_effect, power_difference, metrics))

#   # Test for null treatment effect (should stop with error)
#   treatment_effect <- "null"
#   expect_error(power_vs_tie(test_results_metrics_df, case_study, target_sample_size_per_arm, treatment_effect, power_difference, metrics), "The treatment effect should not be non-positive")

#   # Test for invalid treatment effect
#   treatment_effect <- "invalid"
#   expect_error(power_vs_tie(test_results_metrics_df, case_study, target_sample_size_per_arm, treatment_effect, power_difference, metrics), 'treatment_effect must be "strong", "null" or "moderate"')

#   # Test with power_difference = TRUE
#   treatment_effect <- "strong"
#   power_difference <- TRUE
#   expect_silent(power_vs_tie(test_results_metrics_df, case_study, target_sample_size_per_arm, treatment_effect, power_difference, metrics))
# })

# # Unit tests for power_vs_tie_plots function
# test_that("power_vs_tie_plots function works correctly", {
#   metrics <- c("success_proba", "conf_int_success_proba_lower", "conf_int_success_proba_upper")

#   # Mock configuration file reading
#   mock_read_yaml <- function(file) {
#     list(env = "test_env", config = list())
#   }
#   assignInNamespace("yaml.load_file", mock_read_yaml, ns = "yaml")

#   expect_silent(power_vs_tie_plots(test_results_metrics_df, metrics))
# })
