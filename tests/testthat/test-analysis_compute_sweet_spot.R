test_that("Test Case 1: Normal case with a two crossing points", {
  x <- 1:10
  y <- c(3, 5, 7, 9, 11, 8, 6, 4, 2, 0)
  ref_value <- 5
  expected_result <- list(sweet_spot_bounds = c(2, 7.5), sweet_spot_width = 5.5)
  result <- sweet_spot_determination(x, y, ref_value)
  expect_equal(result, expected_result)
})

test_that("Test Case 2: Normal case with three crossing points", {
  x <- 1:10
  y <- c(3, 5, 7, 9, 11, 8, 6, 4, 6, 8)
  ref_value <- 6
  expected_result <- list(sweet_spot_bounds = c(2.5, 7, 9), sweet_spot_width = 4.5)
  result <- sweet_spot_determination(x, y, ref_value)
  expect_equal(result, expected_result)
})

test_that("Test Case 3: No crossing points", {
  x <- 1:10
  y <- c(3, 5, 7, 9, 11, 8, 6, 4, 2, 0)
  ref_value <- 12
  expected_result <- list(sweet_spot_bounds = c(NA, NA), sweet_spot_width = NA)
  result <- sweet_spot_determination(x, y, ref_value)
  expect_equal(result, expected_result)
})

test_that("Test Case 4: Empty y_values vector", {
  x <- 1:10
  y <- numeric(0)
  ref_value <- 5
  expected_result <- list(sweet_spot_bounds = c(NA, NA), sweet_spot_width = NA)
  result <- sweet_spot_determination(x, y, ref_value)
  expect_equal(result, expected_result)
})

test_that("Test Case 5: All y_values are equal to reference_value", {
  x <- 1:10
  y <- rep(5, 10)
  ref_value <- 5
  expected_result <- list(sweet_spot_bounds = c(NA, NA), sweet_spot_width = NA)
  result <- sweet_spot_determination(x, y, ref_value)
  expect_equal(result, expected_result)
})
