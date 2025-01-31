# Mock model objects for testing
model_mcmc_true <- list(
  mcmc = TRUE,
  credible_interval_97.5 = 1.5,
  credible_interval_2.5 = 0.5
)

model_mcmc_false <- list(
  mcmc = FALSE,
  posterior_cdf = function(x) {
    ifelse(x < 1, 0.3, 0.7)
  }
)

test_that("test_decision works for mcmc model with left null_space", {
  result <- model_mcmc_true$test_decision(critical_value = 0.975, theta_0 = 0.3, null_space = "left", confidence_level = 0.95)
  expect_true(result)
})

test_that("test_decision works for mcmc model with right null_space", {
  result <- model_mcmc_true$test_decision(critical_value = 0.975, theta_0 = 2, null_space = "right", confidence_level = 0.95)
  expect_true(result)
})

test_that("test_decision works for non-mcmc model with left null_space", {
  result <- model_mcmc_false$test_decision(critical_value = 0.975, theta_0 = 0.5, null_space = "left", confidence_level = 0.95)
  expect_false(result)
})

test_that("test_decision works for non-mcmc model with right null_space", {
  result <- model_mcmc_false$test_decision(critical_value = 0.975, theta_0 = 1.5, null_space = "right", confidence_level = 0.95)
  expect_false(result)
})

test_that("test_decision throws error for invalid critical_value with mcmc model", {
  expect_error(
    model_mcmc_true$test_decision(critical_value = 0.95, theta_0 = 0.3, null_space = "left", confidence_level = 0.95),
    "Only the 97.5 and 2.5 percentiles are computed, so you can only consider "
  )
})
