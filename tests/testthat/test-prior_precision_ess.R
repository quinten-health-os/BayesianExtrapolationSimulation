# Mock the RBesT::ess function to return a fixed value for testing
mocked_gaussian_mix_precision_ess <- function(model, method = "moment") {
  return(200)
}

test_that("gaussian_mix_moment_ess calculates ESS correctly with provided sigma", {
  with_mock(`gaussian_mix_precision_ess` = mocked_gaussian_mix_precision_ess, {
    rbest_model <- list()
    target_data <- list(
      sample = list(treatment_effect_standard_error = 0.2),
      summary_measure_likelihood = "normal",
      target_sample_size_per_arm = 100
    )
    expected_result <- 200 - target_data$sample_size_per_arm
    result <- prior_precision_ess(rbest_model, target_data)
    expect_equal(result, expected_result, tolerance = 1e-4)
  })
})
