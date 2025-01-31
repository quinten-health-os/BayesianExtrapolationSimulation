# Test standard_error_log_odds_ratio function
test_that("standard_error_log_odds_ratio calculates the standard error correctly", {
  se <- standard_error_log_odds_ratio(10, 20, 30, 40)
  expect_equal(se, 0.449, tolerance = 1e-3)
})

# Test generate_binary_data_from_rate function
test_that("generate_binary_data_from_rate generates binary data correctly", {
  data <- generate_binary_data_from_rate(0.3, 100)
  expect_equal(length(data), 100)
  expect_equal(sum(data), 30)
})

# Test rate_from_drift_logOR function
test_that("rate_from_drift_logOR calculates the success rate correctly", {
  rate <- rate_from_drift_logOR(0.5, 0.2)
  expect_equal(rate, 0.292, tolerance = 1e-3)
})

# Test rate_from_drift_logRR function
test_that("rate_from_drift_logRR calculates the rate correctly", {
  rate <- rate_from_drift_logRR(0.5, 0.2)
  expect_equal(rate, 0.33, tolerance = 1e-3)
})

# Test compute_ORs function
test_that("compute_ORs computes odds ratios correctly", {
  result <- compute_ORs(100, 200, 30, 40)
  expect_equal(result$log_odds_ratio, -0.539, tolerance = 1e-3)
  expect_equal(result$treatment_rate, 0.2)
  expect_equal(result$control_rate, 0.3)
  expect_equal(result$std_err_log_odds_ratio, 0.279, tolerance = 1e-3)
})

# Test compute_log_odds_ratio_from_rates function
test_that("compute_log_odds_ratio_from_rates computes log odds ratio correctly", {
  log_or <- compute_log_odds_ratio_from_rates(0.2, 0.3)
  expect_equal(log_or, -0.539, tolerance = 1e-3)
})

# Test compute_log_odds_ratio_from_counts function
test_that("compute_log_odds_ratio_from_counts computes log odds ratio correctly", {
  log_or <- compute_log_odds_ratio_from_counts(30, 40, 70, 160, TRUE)
  expect_equal(log_or, -0.539, tolerance = 1e-3)
})

# Test sample_log_odds_ratios function
test_that("sample_log_odds_ratios samples log odds ratios correctly", {
  result <- sample_log_odds_ratios(100, 200, 0.3, 0.2, 10)
  expect_equal(length(result$log_odds_ratio), 10)
  expect_equal(length(result$treatment_rate), 10)
  expect_equal(length(result$control_rate), 10)
  expect_equal(length(result$std_err_log_odds_ratio), 10)
})

# Test sample_rate_ratios function
test_that("sample_rate_ratios samples rate ratios correctly", {
  result <- sample_rate_ratios(0.3, 0.2, 10, 100, 200)
  expect_equal(length(result), 10)
})

# Test sample_aggregate_normal_data function
test_that("sample_aggregate_normal_data samples aggregate normal data correctly", {
  result <- sample_aggregate_normal_data(10, 5, 10, 100)
  expect_equal(length(result$sample_mean), 10)
  expect_equal(length(result$sample_standard_error), 10)
})

# Test sample_aggregate_binary_data function
test_that("sample_aggregate_binary_data samples aggregate binary data correctly", {
  result <- sample_aggregate_binary_data(0.3, 100, 10)
  expect_equal(length(result), 10)
})
