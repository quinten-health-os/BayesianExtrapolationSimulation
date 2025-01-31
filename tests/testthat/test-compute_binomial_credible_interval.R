test_that("compute_binomial_credible_interval works with vector samples", {
  samples <- c(1, 0, 1, 1, 0, 1, 0, 1, 1, 0)
  confidence_level <- 0.95
  result <- compute_binomial_credible_interval(samples, confidence_level)

  expect_type(result, "list")
  expect_true("mean" %in% names(result))
  expect_true("conf_int" %in% names(result))
  expect_equal(result$mean, mean(samples))
  expect_equal(nrow(result$conf_int), 1)
})

test_that("compute_binomial_credible_interval works with matrix samples", {
  samples <- matrix(c(1, 0, 1, 1, 0, 1, 0, 1, 1, 0), nrow = 2)
  confidence_level <- 0.95
  result <- compute_binomial_credible_interval(samples, confidence_level)

  expect_type(result, "list")
  expect_true("mean" %in% names(result))
  expect_true("conf_int" %in% names(result))
  expect_equal(result$mean, mean(samples))
  expect_equal(nrow(result$conf_int), nrow(samples))
})

test_that("compute_binomial_credible_interval handles different confidence levels", {
  samples <- c(1, 0, 1, 1, 0, 1, 0, 1, 1, 0)
  confidence_level <- 0.99
  result <- compute_binomial_credible_interval(samples, confidence_level)

  expect_type(result, "list")
  expect_true("mean" %in% names(result))
  expect_true("conf_int" %in% names(result))
  expect_equal(result$mean, mean(samples))
  expect_equal(nrow(result$conf_int), 1)
})

test_that("compute_binomial_credible_interval handles edge cases", {
  samples <- c(1, 1, 1, 1, 1) # All successes
  confidence_level <- 0.95
  result <- compute_binomial_credible_interval(samples, confidence_level)

  expect_equal(result$mean, mean(samples))
  expect_true(all(result$conf_int$lower <= result$mean & result$conf_int$upper >= result$mean))

  samples <- c(0, 0, 0, 0, 0) # All failures
  result <- compute_binomial_credible_interval(samples, confidence_level)

  expect_equal(result$mean, mean(samples))
  expect_true(all(result$conf_int$lower <= result$mean & result$conf_int$upper >= result$mean))
})
