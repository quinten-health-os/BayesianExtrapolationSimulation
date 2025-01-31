# Sample target data for testing
sample_results <- data.frame(
  target_sample_size_per_arm = c(30, 30, 30),
  target_treatment_effect = c(0, 0.5, 1),
  target_standard_deviation = c(1, 1, 1),
  summary_measure_likelihood = c("normal", "normal", "normal"),
  success_proba = c(0.8, 0.6, 0.7),
  mcse_success_proba = c(0.02, 0.03, 0.025),
  conf_int_success_proba_lower = c(0.75, 0.55, 0.65),
  conf_int_success_proba_upper = c(0.85, 0.65, 0.75),
  theta_0 = c(0, 0, 0),
  null_space = c("left", "right", "left")
)

sample_analysis_config <- list(
  frequentist_test = "t-test"
)

test_that("frequentist_power_at_equivalent_tie works with valid input", {
  results <- sample_results
  analysis_config <- sample_analysis_config

  # Mock compute_freq_power function for testing
  mock_compute_freq_power <- function(alpha, target_data, frequentist_test, theta_0, null_space) {
    return(0.8) # Return a dummy power value
  }

  with_mock(
    `compute_freq_power` = mock_compute_freq_power,
    {
      final_results <- frequentist_power_at_equivalent_tie(results, analysis_config)

      expect_true("frequentist_power" %in% names(final_results))
      expect_true("frequentist_test" %in% names(final_results))
      expect_equal(final_results$frequentist_power[1], 0.8)
      expect_equal(final_results$frequentist_test[1], "t-test")
    }
  )
})

test_that("frequentist_power_at_equivalent_tie throws error if treatment effect 0 is not included", {
  results <- sample_results
  results <- results[results$target_treatment_effect != 0, ]
  analysis_config <- sample_analysis_config

  expect_error(
    frequentist_power_at_equivalent_tie(results, analysis_config),
    "Treatment effect = 0 not included."
  )
})

test_that("frequentist_power_at_equivalent_tie works with different configurations", {
  results <- sample_results
  analysis_config <- list(frequentist_test = "z-test")

  # Mock compute_freq_power function for testing
  mock_compute_freq_power <- function(alpha, target_data, frequentist_test, theta_0, null_space) {
    return(0.7) # Return a dummy power value
  }

  with_mock(
    `compute_freq_power` = mock_compute_freq_power,
    {
      final_results <- frequentist_power_at_equivalent_tie(results, analysis_config)

      expect_true("frequentist_power" %in% names(final_results))
      expect_true("frequentist_test" %in% names(final_results))
      expect_equal(final_results$frequentist_power[1], 0.7)
      expect_equal(final_results$frequentist_test[1], "z-test")
    }
  )
})
