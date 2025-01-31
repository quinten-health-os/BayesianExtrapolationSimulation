# Required Libraries
library(ggplot2)
library(RBExT)

# Parameters for the simulation
historical_control_success_rate <- 0.65  # Historical control success rate
effect_size <- 0.12  # 12% effect size for treatment arm
n_trials <- 10000  # Number of simulations
sample_size <- 200  # Sample size per arm
control_proportions <- seq(0.50, 0.80, length.out = 30)  # Range of true control proportions

# Containers for results
mse_separate <- c()
mse_pooled <- c()
mse_single <- c()

psuccess_separate <- c()
psuccess_pooled <- c()
psuccess_single <- c()

# Pooled analysis (borrow control data from historical trial)
historical_control_data <- c(rep(0,35), rep(1, 65))

true_treatment <- 0.7

confidence_level <- 0.95

conf_int_proba_success_separate <- matrix(NA, nrow = length(control_proportions), ncol = 2)

conf_int_proba_success_pooled <- matrix(NA, nrow = length(control_proportions), ncol = 2)

conf_int_proba_success_single <- matrix(NA, nrow = length(control_proportions), ncol = 2)

# Simulation loop
for (j in seq_along(control_proportions)) {

  true_control <- control_proportions[j]
  estimates_separate <- c()
  estimates_pooled <- c()
  estimates_single <- c()

  test_result_separate <- c()
  test_result_pooled <- c()
  test_result_single <- c()

  for (i in 1:n_trials) {

    # Generate control and treatment arm data
    control_data <- rbinom(sample_size, 1, true_control)
    treatment_data <- rbinom(sample_size, 1, true_treatment)

    p_treatment <- mean(treatment_data)

    # Separate analysis
    p_control_separate <- mean(control_data)
    estimates_separate <- c(estimates_separate, p_control_separate)


    separate <- BSDA::tsum.test(
      mean.x = p_treatment,
      mean.y = p_control_separate,
      mu = 0,
      alternative = "greater",
      s.x = sqrt(p_treatment * (1 - p_treatment)),
      s.y = sqrt(p_control_separate * (1 - p_control_separate)),
      n.x = sample_size,
      n.y = sample_size
    )

    test_result_separate <- c(test_result_separate, separate$p.value < 0.05)


    pooled_control_data <- c(control_data, historical_control_data)
    p_control_pooled <- mean(pooled_control_data)
    estimates_pooled <- c(estimates_pooled, p_control_pooled)

    pooled <- BSDA::tsum.test(
      mean.x = p_treatment,
      mean.y = p_control_pooled,
      mu = 0,
      alternative = "greater",
      s.x = sqrt(p_treatment * (1 - p_treatment)),
      s.y = sqrt(p_control_pooled * (1 - p_control_pooled)),
      n.x = sample_size,
      n.y = length(pooled_control_data)
    )

    test_result_pooled <- c(test_result_pooled, pooled$p.value < 0.05)

    # Single-arm analysis (compare treatment arm to historical control)
    p_control_single <- mean(historical_control_data)
    estimates_single <- c(estimates_single, p_control_single)

    single <- BSDA::tsum.test(
      mean.x = p_treatment,
      alternative = "greater",
      mu = p_control_single,
      s.x = sqrt(p_treatment * (1 - p_treatment)),
      n.x = sample_size
    )

    test_result_single <- c(test_result_single, single$p.value < 0.05)
  }

  # Calculate MSE for each design
  mse_separate <- c(mse_separate, mean((estimates_separate - true_control)^2))
  mse_pooled <- c(mse_pooled, mean((estimates_pooled - true_control)^2))
  mse_single <- c(mse_single, mean((estimates_single - true_control)^2))

  psuccess_separate <- c(psuccess_separate, mean(test_result_separate))
  psuccess_pooled <- c(psuccess_pooled, mean(test_result_pooled))
  psuccess_single <- c(psuccess_single, mean(test_result_single))

  conf_int_proba_success_separate[j,] <- binom.test(sum(test_result_separate), length(test_result_separate), conf.level = confidence_level)$conf.int

  conf_int_proba_success_pooled[j,] <- binom.test(sum(test_result_pooled), length(test_result_pooled), conf.level = confidence_level)$conf.int

  conf_int_proba_success_single[j,] <- binom.test(sum(test_result_single), length(test_result_single), conf.level = confidence_level)$conf.int

}

ci_low_separate <- conf_int_proba_success_separate[,1]
ci_low_pooled <- conf_int_proba_success_pooled[,1]
ci_low_single <- conf_int_proba_success_single[,1]

ci_high_separate <- conf_int_proba_success_separate[,2]
ci_high_pooled <- conf_int_proba_success_pooled[,2]
ci_high_single <- conf_int_proba_success_single[,2]


# Prepare data for plotting
  data_ci <- data.frame(
    control_proportion = rep(control_proportions, times=3),
    Power = c(psuccess_separate, psuccess_pooled, psuccess_single),
    ci_low = c(ci_low_separate, ci_low_pooled, ci_low_single),
    ci_high = c(ci_high_separate, ci_high_pooled, ci_high_single),
    Design = rep(c("Separate", "Pooled", "Single Arm"), each=length(control_proportions))
)


# Plot MSE
ggplot(data, aes(x=control_proportion, y=MSE, color=Design)) +
  geom_line(size= LINEWIDTH) +
  # geom_vline(xintercept = historical_control_success_rate, linetype = "dashed", color = "red", size = LINEWIDTH) +
  # Add vline with a legend
  geom_vline(aes(xintercept = historical_control_success_rate, linetype = "Historical Control Success Rate"), color = "black", size = LINEWIDTH) +

  # + geom_vline(xintercept = true_treatment, linetype = "True treatment proportion", color = "blue", size = LINEWIDTH)
  labs(title="Mean Squared Error (MSE)", x="True Control Proportion", y="MSE of Control Estimate") +
  theme_minimal()
plt <- ggplot(data_ci, aes(x = control_proportion, y = Power, color = Design, fill = Design)) +
  geom_line(size = LINEWIDTH) +

  # Add confidence intervals as shaded regions, suppressing the legend
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.2, color = NA, show.legend = FALSE) +

  # Add hline with a legend
  geom_hline(aes(yintercept = 0.05, linetype = "Significance threshold (0.05)"), color = "black", size = LINEWIDTH) +

  # Add title and axis labels
  labs(title = "Probability of success as a function of the true success proportion in the control arm",
       x = "True Control Proportion", y = "Probability of success", linetype = "", color = "") +

  # Set theme
  scale_color_viridis_d() +  # Apply viridis color scale
  scale_fill_viridis_d() +  # Use same scale for fill
  theme_bw() +
  theme(
    text = element_text(family = font, size = text_size),
    axis.text.y = element_text(size = text_size),
    axis.text.x = element_text(size = text_size),
    legend.key = element_blank(),
    # strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +

  # Add custom y-ticks including 0.05
  scale_y_continuous(breaks = c(0.05, seq(0, 1, by = 0.1))) +
  scale_x_continuous(breaks = c(c(0.5, 0.6, 0.8), historical_control_success_rate, true_treatment), labels = c(c(0.5, 0.6, 0.8), TeX("$p_S^{c}$"), TeX("$p_T^{t}= 0.7$")), expand = c(0, 0))
  # scale_x_continuous(breaks = c(seq(0, 1, by = 0.1), historical_control_success_rate, true_treatment), labels = c(seq(0, 1, by = 0.1), TeX("$p_S^{(c)}$"), TeX("$p_T^{(t)}$")), expand = c(0, 0))

print(plt)


# Export the plot
export_plots(
  plt = plt,
  file_path = "./figures/Viele_et_al_variant_with_CI",
  fig_width_in,
  fig_height_in,
  type = "pdf"
)
export_plots(
  plt = plt,
  file_path = "./figures/Viele_et_al_variant_with_CI",
  fig_width_in,
  fig_height_in,
  type = "png"
)
