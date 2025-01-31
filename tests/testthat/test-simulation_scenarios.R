# Test for the hellinger_distance function
test_that("hellinger_distance computes the correct Hellinger distance", {
  expect_equal(hellinger_distance(0, 1, 0, 1), 0)
  expect_equal(hellinger_distance(0, 1, 1, 1), 0.343, tolerance = 1e-3)
  expect_equal(hellinger_distance(0, 1, 0, 2), 0.325, tolerance = 1e-3)
})

# Test for the compute_drift_range function
test_that("compute_drift_range computes the correct drift range", {
  simulation_config <- list(ndrift = 10)
  case_study_config <- list(
    theta_0 = 0,
    source = list(treatment_effect = 1, standard_error = 0.5),
    target = list(treatment_effect = 2, standard_error = 1),
    summary_measure_likelihood = "normal"
  )

  expected_drift_range <- c(
    -2.78,
    -2.32,
    -1.85,
    -1.39,
    -9.27e-01,
    -4.63e-01,
    4.44e-16,
    4.63e-01,
    9.27e-01,
    1.390
  )

  expect_equal(
    compute_drift_range(simulation_config, case_study_config),
    expected_drift_range,
    tolerance = 1e-3
  )
})

# Test for the important_drift_values function
test_that("important_drift_values computes the correct important drift values", {
  source_treatment_effect <- 1
  case_study_config <- list(summary_measure_likelihood = "normal")

  expected_important_drift_values <- c(-1, -0.5, 0)

  expect_equal(
    important_drift_values(source_treatment_effect, case_study_config),
    expected_important_drift_values
  )
})

# Test for the compute_control_drift_range function
test_that("compute_control_drift_range computes the correct control drift range", {
  source_treatment_effect <- 1

  expected_control_drift_range <- c(0)

  expect_equal(
    compute_control_drift_range(source_treatment_effect),
    expected_control_drift_range
  )
})

# Test for the compute_source_denominator_range function
test_that("compute_source_denominator_range computes the correct source denominator range", {
  source_data <- list(
    endpoint = "binary",
    control_rate = 0.2,
    treatment_effect_estimate = 0.5
  )
  simulation_config <- list(denominator_change_factor = c(0.8, 1, 1.2))
  case_study_config <- list(summary_measure_likelihood = "normal")

  expected_source_denominator_range <- list(
    source_denominator = c(0.20, 0.25, 0.30),
    source_denominator_change_factor = c(0.8, 1, 1.2)
  )
  expect_equal(
    compute_source_denominator_range(source_data, simulation_config, case_study_config),
    expected_source_denominator_range
  )
})
