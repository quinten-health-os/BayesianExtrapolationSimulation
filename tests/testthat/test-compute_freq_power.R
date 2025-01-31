# Sample target data for testing
target_data <- list(
  sample_size_per_arm = 30,
  treatment_effect = 1,
  standard_deviation = 1
)

test_that("compute_freq_power works for t-test with left null space", {
  alpha <- 0.05
  frequentist_test <- "t-test"
  theta_0 <- 0
  null_space <- "left"
  result <- compute_freq_power(alpha, target_data, frequentist_test, theta_0, null_space)

  effect_size <- (target_data$treatment_effect - theta_0) / target_data$standard_deviation
  expected_power <- pwr::pwr.t.test(d = effect_size, n = target_data$sample_size_per_arm, sig.level = alpha, type = "two.sample", alternative = "greater")$power

  expect_equal(result, expected_power, tolerance = 1e-6)
})

test_that("compute_freq_power works for t-test with right null space", {
  alpha <- 0.05
  frequentist_test <- "t-test"
  theta_0 <- 0
  null_space <- "right"
  result <- compute_freq_power(alpha, target_data, frequentist_test, theta_0, null_space)

  effect_size <- (target_data$treatment_effect - theta_0) / target_data$standard_deviation
  expected_power <- pwr::pwr.t.test(d = effect_size, n = target_data$sample_size_per_arm, sig.level = alpha, type = "two.sample", alternative = "less")$power

  expect_equal(result, expected_power, tolerance = 1e-6)
})

test_that("compute_freq_power works for z-test with left null space", {
  alpha <- 0.05
  frequentist_test <- "z-test"
  theta_0 <- 0
  null_space <- "left"
  result <- compute_freq_power(alpha, target_data, frequentist_test, theta_0, null_space)

  effect_size <- (target_data$treatment_effect - theta_0) / target_data$standard_deviation
  expected_power <- pwr::pwr.norm.test(d = effect_size, n = target_data$sample_size_per_arm, sig.level = alpha, alternative = "greater")$power

  expect_equal(result, expected_power, tolerance = 1e-6)
})

test_that("compute_freq_power works for z-test with right null space", {
  alpha <- 0.05
  frequentist_test <- "z-test"
  theta_0 <- 0
  null_space <- "right"
  result <- compute_freq_power(alpha, target_data, frequentist_test, theta_0, null_space)

  effect_size <- (target_data$treatment_effect - theta_0) / target_data$standard_deviation
  expected_power <- pwr::pwr.norm.test(d = effect_size, n = target_data$sample_size_per_arm, sig.level = alpha, alternative = "less")$power

  expect_equal(result, expected_power, tolerance = 1e-6)
})

test_that("compute_freq_power throws error for invalid null space", {
  alpha <- 0.05
  frequentist_test <- "t-test"
  theta_0 <- 0
  null_space <- "middle"

  expect_error(
    compute_freq_power(alpha, target_data, frequentist_test, theta_0, null_space),
    "Null space must be either 'left' or 'right'"
  )
})
