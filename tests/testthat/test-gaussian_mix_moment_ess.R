# Mock a mixture object for testing
mock_mixture <- function() {
  mix <- mixnorm(rob = c(1, 0, 1), inf = c(1, 0, 1), sigma = 1) # Simplified example of a mixture object
  return(mix)
}

# Mock the RBesT::sigma function to return a fixed value for testing
mock_sigma <- function(mix) {
  return(2)
}

test_that("gaussian_mix_moment_ess calculates ESS correctly with provided sigma", {
  with_mock(
    `RBesT::sigma` = mock_sigma,
    {
      mix <- mock_mixture()
      expected_result <- 2^2 / 1^2
      result <- gaussian_mix_moment_ess(mix)
      expect_equal(result, expected_result, tolerance = 1e-4)
    }
  )
})
