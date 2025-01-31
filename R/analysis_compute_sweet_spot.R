#' Determine Sweet Spot Bounds and Width
#'
#' @description This function calculates the bounds and width of the "sweet spot" where a given `y_values` vector
#' crosses a specified `reference_value`. It uses linear interpolation for precise crossing points.
#'
#' @param x_values A numeric vector of x values.
#' @param y_values A numeric vector of y values.
#' @param reference_value A numeric reference value to determine the crossing points.
#'
#' @return A list containing:
#' \itemize{
#'     \item sweet_spot_bounds A numeric vector of length 2 with the lower and upper bounds of the sweet spot.
#'     \item sweet_spot_widthA numeric value indicating the width of the sweet spot.
#'   }
sweet_spot_determination <- function(x_values, y_values, reference_value = 0, larger_is_better, return_NA_at_boundaries = FALSE) {
  # Check for non-zero length input vectors
  if (length(x_values) == 0 || length(y_values) == 0) {
    warning("Input vectors must have non-zero length.")
    sweet_spot <- data.frame(start = NA, end = NA, width = NA)
    return(sweet_spot)
  }

  xrange <- range(x_values)

  # Adjust y_values and reference_value if larger_is_better is FALSE
  if (!larger_is_better) {
    y_values <- -y_values
    reference_value <- -reference_value
  }

  # Compute the sign of y_values relative to reference_value
  sign_y <- sign(y_values - reference_value)
  sign_y[sign_y == 0] <- 1  # Treat zero as positive

  # Find indices where the sign changes (crossings)
  crossings <- which(diff(sign_y) != 0)

  # Handle the case with no crossings
  if (length(crossings) == 0) {
    if (return_NA_at_boundaries) {
      sweet_spot <- data.frame(start = NA, end = NA, width = NA, total_width = NA)
    } else {
      if (all(sign_y == 1)) {
        # Sweet spot covers the entire range
        sweet_spot <- data.frame(
          start = xrange[1],
          end = xrange[2],
          width = diff(xrange),
          total_width = diff(xrange)
        )
      } else {
        # No sweet spot within the range
        sweet_spot <- data.frame(start = NA, end = NA, width = 0, total_width = 0)
      }
    }
    return(sweet_spot)
  }

  # Calculate precise crossing points using linear interpolation
  precise_crossings <- numeric(length(crossings))
  for (i in seq_along(crossings)) {
    idx <- crossings[i]
    y_pair <- y_values[c(idx, idx + 1)]
    x_pair <- x_values[c(idx, idx + 1)]
    precise_crossings[i] <- approx(y = x_pair, x = y_pair, xout = reference_value)$y
  }

  # Initialize intervals between crossings
  intervals_start <- c(xrange[1], precise_crossings)
  intervals_end <- c(precise_crossings, xrange[2])
  n_intervals <- length(intervals_start)

  # Determine the initial sign of the function
  initial_sign <- sign_y[1]

  # Compute interval signs (alternating signs)
  interval_signs <- initial_sign * (-1)^(seq_len(n_intervals) - 1)

  # Identify sweet spots where the sign is positive
  sweet_spot_indices <- which(interval_signs == 1)
  sweet_spot <- data.frame(
    start = intervals_start[sweet_spot_indices],
    end = intervals_end[sweet_spot_indices]
  )

  # Adjust bounds if return_NA_at_boundaries is TRUE
  if (return_NA_at_boundaries) {
    sweet_spot$start[sweet_spot$start == xrange[1]] <- NA
    sweet_spot$end[sweet_spot$end == xrange[2]] <- NA
  }

  # Calculate the width of each sweet spot
  sweet_spot$width <- sweet_spot$end - sweet_spot$start

  # Warning message for multiple sweet spots
  if (nrow(sweet_spot) > 1) {
    warning("More than one sweet spot.")
  }

  sweet_spot$total_width <- sum(sweet_spot$width)
  # Return the results
  return(sweet_spot)
}

compare_metrics <- function(merged_df, metric, based_on_CI = based_on_CI, suffix = ""){
  if (based_on_CI){
    if (metric$larger_is_better){
     difference <- merged_df[[paste0("conf_int_", metric$name, "_lower")]] - merged_df[[paste0("conf_int_", metric$name, "_upper", suffix)]]
    } else {
      difference <- merged_df[[paste0("conf_int_", metric$name, "_upper")]] - merged_df[[paste0("conf_int_", metric$name, "_lower", suffix)]]
    }
  } else {
    difference <- merged_df[[metric$name]] - merged_df[[paste0(metric$name, suffix)]]
  }

  return(difference)
}


#' Calculate Sweet Spots for Multiple Metrics
#'
#' @description This function calculates the sweet spots for various metrics over multiple case studies and methods.
#' It determines where the metric differences between two methods cross a reference value.
#'
#' @param results_freq_df A dataframe containing the results frequency data with columns
#' "case_study", "method", "source_denominator_change_factor", "control_drift", and "target_sample_size_per_arm".
#' @param metrics A list of metrics to evaluate. Each metric should have a `name` attribute.
#'
#' @return A dataframe containing the sweet spots for each metric, case study, method, and combination
#' of parameters.
#'
#' @examples
#' \dontrun{
#' results <- sweet_spot(results_freq_df, metrics)
#' }
sweet_spot <- function(results_freq_df, metrics, based_on_CI = TRUE) {

  # Initialize an empty dataframe to store the results
  final_df <- data.frame(
    metric = character(),
    method = character(),
    parameters = character(),
    target_sample_size_per_arm = numeric(),
    control_drift = numeric(),
    source_denominator = numeric(),
    source_denominator_change_factor = numeric(),
    case_study = character(),
    parallelization = logical(),
    sampling_approximation = logical(),
    source_treatment_effect_estimate = numeric(),
    target_to_source_std_ratio = numeric(),
    theta_0 = numeric(),
    null_space = character(),
    summary_measure_likelihood = character(),
    source_standard_error = numeric(),
    source_sample_size_control = numeric(),
    source_sample_size_treatment = numeric(),
    equivalent_source_sample_size_per_arm = numeric(),
    endpoint = character(),
    source_control_rate = numeric(),
    source_treatment_rate = numeric(),
    sweet_spot_lower = numeric(),
    sweet_spot_upper = numeric(),
    sweet_spot_width = numeric(),
    drift_range_upper = numeric(),
    drift_range_lower = numeric(),
    stringsAsFactors = FALSE
  )

  # Loop through each unique combination to calculate sweet spots
  for (case_study in unique(results_freq_df$case_study)) {
      results_freq_df1 = results_freq_df[results_freq_df$case_study == case_study, ]
      for (target_sample_size_per_arm in unique(results_freq_df1$target_sample_size_per_arm)) {
        results_freq_df2 = results_freq_df1[results_freq_df1$target_sample_size_per_arm == target_sample_size_per_arm, ]
        for (source_denominator_change_factor in unique(results_freq_df2$source_denominator_change_factor)) {
          results_freq_df3 = results_freq_df2[results_freq_df2$source_denominator_change_factor == source_denominator_change_factor, ]
          for (target_to_source_std_ratio in unique(results_freq_df3$target_to_source_std_ratio)){
            results_freq_df4 <- results_freq_df3 %>%
              dplyr::filter(
                target_to_source_std_ratio == !!target_to_source_std_ratio |
                  is.na(target_to_source_std_ratio)
              )

            df <- results_freq_df4[results_freq_df4$control_drift == 0, ]

            separate_df <- df[df$method == "separate", ]

            if (nrow(separate_df) == 0) {
              warning(
                "A separate analysis but be performed for comparison and sweet spot computation."
              )
            }

            if (nrow(df) > 0) {
              # Round the values in the dataframe as R does not necessarily save the data with the same precison across different environmnents.
              separate_df <- separate_df %>%
                mutate_if(is.numeric, round, digits = 12)
              df <- df %>%
                mutate_if(is.numeric, round, digits = 12)

              merged_df <- dplyr::left_join(df, separate_df, by = c("drift", "target_treatment_effect", scenario_columns), suffix = c("", "_separate"))

              if (nrow(merged_df) != nrow(df)){
                stop("The number of rows in the merged df does not match the number of rows in the results df.")
              }

              # Sort the dataframe by drift in ascending order
              merged_df <- merged_df %>%
                dplyr::arrange(drift)

              for (metric in metrics) {
                # Depending on whether the sweet spot is determined based on the mean or on the bounds of the CI, compute the difference based on which it will be computed.
                merged_df$metric_diff <- compare_metrics(merged_df, metric, based_on_CI = based_on_CI, suffix = "_separate")

                # Determine the sweet spot location and width
                for (method in unique(merged_df$method)){
                  parameters_combinations <- unique(merged_df[merged_df$method == method,]$parameters)
                  for (parameters in parameters_combinations){
                    method_df <- merged_df[merged_df$parameters == parameters & merged_df$method == method,]

                    # Drop the rows where metric_diff is NA.
                    method_df <- method_df[!is.na(method_df$metric_diff), ]
                    sweet_spot <- sweet_spot_determination(method_df$drift,
                                                           method_df$metric_diff,
                                                           larger_is_better = metric$larger_is_better)

                    if (nrow(method_df) == 0){
                      warning("Dataframe is empty")
                      next
                    }

                    if (nrow(sweet_spot) > 1){
                      sweet_spot <- data.frame(
                        start = I(list(sweet_spot$start)),
                        end = I(list(sweet_spot$end)),
                        width = I(list(sweet_spot$width)),
                        total_width = I(list(sweet_spot$total_width))
                      )
                    }
                    drift_range <- range(method_df$drift)


                    df_to_add <- cbind(unique(method_df[, c("method", "parameters", scenario_columns)]), data.frame(
                      sweet_spot_lower = sweet_spot$start,
                      sweet_spot_upper = sweet_spot$end,
                      sweet_spot_width = sweet_spot$width,
                      sweet_spot_total_width = sweet_spot$total_width,
                      drift_range_upper = drift_range[2],
                      drift_range_lower = drift_range[1]
                    ))

                    df_to_add$metric <- metric$name

                    # Append the results to final_df
                    final_df <- rbind(
                      final_df,
                        df_to_add
                      )

                    if (metric$name == "success_proba") {
                      method_df$metric_diff <- method_df[["success_proba"]] - analysis_config$nominal_tie

                      # Determine the range of drift values for which the success proba is smaller than the nominal TIE
                      sweet_spot <- sweet_spot_determination(method_df$drift,
                                                             method_df$metric_diff,
                                                             larger_is_better = FALSE)

                      df_to_add <- cbind(unique(method_df[, c("method", "parameters", scenario_columns)]), data.frame(
                        sweet_spot_lower = sweet_spot$start,
                        sweet_spot_upper = sweet_spot$end,
                        sweet_spot_width = sweet_spot$width,
                        sweet_spot_total_width = sweet_spot$total_width,
                        drift_range_upper = drift_range[2],
                        drift_range_lower = drift_range[1]
                      )
                      )

                      df_to_add$metric <- "success_proba_smaller_than_nominal_TIE"

                      final_df <- rbind(
                        final_df,
                        df_to_add
                      )

                      method_df$metric_diff <- compare_metrics(method_df, metric, based_on_CI = based_on_CI, suffix = "_separate")

                      # Determine the range of drift values for which the success proba is higher than the nominal power (the power of a frequentist analysis at nominal TIE)
                      sweet_spot <- sweet_spot_determination(method_df$drift,
                                                             method_df$metric_diff,
                                                             larger_is_better =  TRUE)


                      df_to_add <- cbind(unique(method_df[, c("method", "parameters", scenario_columns)]), data.frame(
                        sweet_spot_lower = sweet_spot$start,
                        sweet_spot_upper = sweet_spot$end,
                        sweet_spot_width = sweet_spot$width,
                        sweet_spot_total_width = sweet_spot$total_width,
                        drift_range_upper = drift_range[2],
                        drift_range_lower = drift_range[1]
                      )
                      )

                      df_to_add$metric <- "power_larger_than_nominal"

                      final_df <- rbind(
                        final_df,
                        df_to_add
                      )
                    }
                  }
                }
              }
            }
          }
        }
      }
  }

  return(final_df)
}
