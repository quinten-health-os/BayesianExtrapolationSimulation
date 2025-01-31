analyze_power_gains <- function(results_freq_df, output_path){
  # Check whether there are cases where the success proba is higher in the alternative space than the frequentist power at equivalent TIE
  subdf <- subset(results_freq_df,
                  (null_space == "left" & target_treatment_effect > theta_0) |
                    (null_space == "right" & target_treatment_effect < theta_0))

  # Determine whether the power of the method of interest is higher than the power of a separate analysis at equivalent TIE
  # power_gain_cases <- subset(subdf, conf_int_success_proba_lower > frequentist_power_at_equivalent_tie_upper)

  subdf$se_success_proba <- (subdf$conf_int_success_proba_upper - subdf$conf_int_success_proba_lower) / (2 * 1.96)
  subdf$se_power_at_equivalent_tie <- (subdf$frequentist_power_at_equivalent_tie_upper - subdf$frequentist_power_at_equivalent_tie_lower) / (2 * 1.96)

  n <- 12
  # Apply the hypothesis test for each row in the dataframe
  subdf$p_value <- mapply(function(x, y, se_x, se_y) {
    # Perform a two-sample z-test using BSDA's z.test function
    # Determine whether the mean of the first sample (mean.x) is significantly greater than the mean of the second sample (mean.y).
    z_test <- BSDA::zsum.test(
      mean.x = x,             # First sample mean
      mean.y = y,             # Second sample mean
      sigma.x = se_x*sqrt(n),    # Standard error of x
      n.x = n,
      sigma.y = se_y*sqrt(n),     # Standard error of y
      n.y = n,
      alternative = "greater"
    )

    # Extract the p-value
    z_test$p.value
  },
  x = subdf$success_proba,
  y = subdf$frequentist_power_at_equivalent_tie,
  se_x = subdf$se_success_proba,
  se_y = subdf$se_power_at_equivalent_tie
  )


  # Subset based on a significance level (e.g., alpha = 0.05)
  power_gain_cases <- subset(subdf, p_value < 0.05)


  # Determine whether the power of the method of interest is higher than the power of the frequentist method at the nominal TIE (otherwise there is no point in borrowing)
  power_gain_cases <- subset(power_gain_cases, conf_int_success_proba_lower > nominal_frequentist_power_separate)


  # Detect anomalies where a separate analysis leads to power gains
  separate_analysis_gain <- subset(power_gain_cases, method == "separate")

  write.csv(file = paste0(output_path, "/", "power_gain_cases.csv"), power_gain_cases)#, power_gain_cases[,expected_colnames_scenario])

  write.csv(file = paste0(output_path, "/", "power_gain_cases_separate_analysis.csv"), separate_analysis_gain)

}

analyze_power_loss <- function(results_freq_df, output_path){
  # Check whether there are cases where the success proba is lower in the alternative space than the frequentist power at equivalent TIE while TIE is inflated compared to the nominal TIE.
  subdf <- subset(results_freq_df,
                  (null_space == "left" & target_treatment_effect > theta_0) |
                    (null_space == "right" & target_treatment_effect < theta_0))

  power_loss_cases <- subset(subdf, conf_int_success_proba_upper < frequentist_power_at_equivalent_tie_lower)


  write.csv(file = paste0(output_path, "/","power_loss_cases.csv"), power_loss_cases[,expected_colnames_scenario])
}

analyze_power_loss_inflated_tie <- function(results_freq_df, output_path){
  # Check whether there are cases where the success proba is lower in the alternative space than the frequentist power at equivalent TIE while TIE is inflated compared to the nominal TIE.
  subdf <- subset(results_freq_df,
                  (null_space == "left" & target_treatment_effect > theta_0) |
                    (null_space == "right" & target_treatment_effect < theta_0))

  subdf <- subset(subdf, conf_int_tie_lower > 0.025)

  subdf$se_success_proba <- (subdf$conf_int_success_proba_upper - subdf$conf_int_success_proba_lower) / (2 * 1.96)
  subdf$se_power_at_equivalent_tie <- (subdf$frequentist_power_at_equivalent_tie_upper - subdf$frequentist_power_at_equivalent_tie_lower) / (2 * 1.96)

  #power_loss_cases <- subset(subdf, conf_int_success_proba_upper < frequentist_power_at_equivalent_tie_lower)
  # Perform the one-sided hypothesis test
  n <- 12
  subdf$p_value <- mapply(function(x, y, se_x, se_y) {
    # Perform a one-sided z-test for success_proba < power_at_equivalent_tie
    #  determine whether the mean of the first sample (mean.x) is significantly less than the mean of the second sample (mean.y).
    z_test <- BSDA::zsum.test(
      mean.x = x,             # First sample mean (success_proba)
      mean.y = y,             # Second sample mean (power_at_equivalent_tie)
      sigma.x = se_x*sqrt(n),    # Standard error of x
      sigma.y = se_y*sqrt(n),    # Standard error of y
      n.x = n,
      n.y = n,
      alternative = "less"  # Specify the one-sided alternative hypothesis
    )

    # Extract the p-value
    z_test$p.value
  },
  x = subdf$success_proba,
  y = subdf$frequentist_power_at_equivalent_tie,
  se_x = subdf$se_success_proba,
  se_y = subdf$se_power_at_equivalent_tie
  )

  # Subset based on a significance level (e.g., alpha = 0.05)
  power_loss_cases <- subset(subdf, p_value < 0.05)


  write.csv(file = paste0(output_path, "/", "power_loss_inflated_tie_cases.csv"), power_loss_cases[,expected_colnames_scenario])

}

analyze_noninflated_tie <- function(results_freq_df, output_path){
  # Check whether there are cases where the TIE is not inflated

  non_inflated_tie_cases <- subset(results_freq_df, conf_int_tie_upper < 0.025)
  non_inflated_tie_cases <- subset(non_inflated_tie_cases, drift == - source_treatment_effect_estimate)

  # Remove the cases where there was no borrowing
  # filter <- non_inflated_tie_cases$conf_int_ess_moment_upper >= 0 & non_inflated_tie_cases$conf_int_ess_moment_lower >= 0
  # non_inflated_tie_cases <- non_inflated_tie_cases[!filter,]


  write.csv(file = paste0(output_path, "/", "noninflated_tie_cases.csv"), non_inflated_tie_cases)

}
