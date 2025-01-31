# Required Libraries
library(RBExT)
library(ggplot2)
library(gridExtra)

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

type1_error_separate <- c()
type1_error_pooled <- c()
type1_error_single <- c()

power_separate <- c()
power_pooled <- c()
power_single <- c()

# Pooled analysis (borrow control data from historical trial)
historical_control_data <- c(rep(0,35), rep(1, 65))

# Simulation loop
for (true_control in control_proportions) {

  estimates_separate <- c()
  estimates_pooled <- c()
  estimates_single <- c()

  reject_null_separate <- c()
  reject_null_pooled <- c()
  reject_null_single <- c()

  power_detected_separate <- c()
  power_detected_pooled <- c()
  power_detected_single <- c()

  for (i in 1:n_trials) {

    # Generate control and treatment arm data
    control_data <- rbinom(sample_size, 1, true_control)
    null_treatment_data <- rbinom(sample_size, 1, true_control)
    alt_treatment_data <- rbinom(sample_size, 1, true_control + effect_size)

    p_treatment_null <- mean(null_treatment_data)
    p_treatment_alt <- mean(alt_treatment_data)

    # Separate analysis
    p_control_separate <- mean(control_data)
    estimates_separate <- c(estimates_separate, p_control_separate)

    # Create contingency table for Fisher's Exact Test
    table_separate_null <- matrix(c(
      sum(null_treatment_data == 1), sum(null_treatment_data == 0),
      sum(control_data == 1), sum(control_data == 0)
    ), nrow = 2)

    # Fisher's Exact Test for null hypothesis
    fisher_null_separate <- fisher.test(table_separate_null, alternative = "greater")
    reject_null_separate <- c(reject_null_separate, fisher_null_separate$p.value < 0.05)

    table_separate_alt <- matrix(c(
      sum(alt_treatment_data == 1), sum(alt_treatment_data == 0),
      sum(control_data == 1), sum(control_data == 0)
    ), nrow = 2)

    # Fisher's Exact Test for alternative hypothesis
    fisher_alt_separate <- fisher.test(table_separate_alt, alternative = "greater")
    power_detected_separate <- c(power_detected_separate, fisher_alt_separate$p.value < 0.05)


    # Pooled analysis
    pooled_control_data <- c(control_data, historical_control_data)
    p_control_pooled <- mean(pooled_control_data)
    estimates_pooled <- c(estimates_pooled, p_control_pooled)

    table_pooled <- matrix(c(
      sum(null_treatment_data == 1), sum(null_treatment_data == 0),
      sum(pooled_control_data == 1), sum(pooled_control_data == 0)
    ), nrow = 2)

    # Fisher's Exact Test for pooled data
    fisher_pooled <- fisher.test(table_pooled, alternative = "greater")
    reject_null_pooled <- c(reject_null_pooled, fisher_pooled$p.value < 0.05)

    table_pooled_alt <- matrix(c(
      sum(alt_treatment_data == 1), sum(alt_treatment_data == 0),
      sum(pooled_control_data == 1), sum(pooled_control_data == 0)
    ), nrow = 2)

    # Fisher's Exact Test for alternative hypothesis with pooled data
    fisher_alt_pooled <- fisher.test(table_pooled_alt, alternative = "greater")
    power_detected_pooled <- c(power_detected_pooled, fisher_alt_pooled$p.value < 0.05)

    # Single-arm analysis (compare treatment arm to historical control)
    p_control_single <- mean(historical_control_data)
    estimates_single <- c(estimates_single, p_control_single)

    table_single <- matrix(c(
      sum(null_treatment_data == 1), sum(null_treatment_data == 0),
      sum(historical_control_data == 1), sum(historical_control_data == 0)
    ), nrow = 2)

    # Fisher's Exact Test for single-arm analysis
    fisher_single <- fisher.test(table_single, alternative = "greater")
    reject_null_single <- c(reject_null_single, fisher_single$p.value < 0.05)

    table_single_alt <- matrix(c(
      sum(alt_treatment_data == 1), sum(alt_treatment_data == 0),
      sum(historical_control_data == 1), sum(historical_control_data == 0)
    ), nrow = 2)

    # Fisher's Exact Test for alternative hypothesis with single-arm data
    fisher_alt_single <- fisher.test(table_single_alt, alternative = "greater")
    power_detected_single <- c(power_detected_single, fisher_alt_single$p.value < 0.05)
  }

  # Calculate MSE for each design
  mse_separate <- c(mse_separate, mean((estimates_separate - true_control)^2))
  mse_pooled <- c(mse_pooled, mean((estimates_pooled - true_control)^2))
  mse_single <- c(mse_single, mean((estimates_single - true_control)^2))

  # Calculate Type I error rates
  type1_error_separate <- c(type1_error_separate, mean(reject_null_separate))
  type1_error_pooled <- c(type1_error_pooled, mean(reject_null_pooled))
  type1_error_single <- c(type1_error_single, mean(reject_null_single))

  # Calculate Power
  power_separate <- c(power_separate, mean(power_detected_separate))
  power_pooled <- c(power_pooled, mean(power_detected_pooled))
  power_single <- c(power_single, mean(power_detected_single))
}

# Prepare data for plotting
data <- data.frame(
  control_proportion = rep(control_proportions, times=3),
  MSE = c(mse_separate, mse_pooled, mse_single),
  Type1Error = c(type1_error_separate, type1_error_pooled, type1_error_single),
  Power = c(power_separate, power_pooled, power_single),
  Design = rep(c("Separate", "Pooled", "Single Arm"), each=length(control_proportions))
)

# Plot MSE
ggplot(data, aes(x=control_proportion, y=MSE, color=Design)) +
  geom_line(size=1) +
  labs(title="Mean Squared Error (MSE)", x="True Control Proportion", y="MSE of Control Estimate") +
  theme_minimal()

# Plot Type I Error
plot1 <- ggplot(data, aes(x=control_proportion, y=Type1Error, color=Design)) +
  geom_line(size=1) +
  # Add hline with a legend
  geom_hline(aes(yintercept = 0.05,
                 linetype = "Significance Threshold (0.05)"), color = "black", size = LINEWIDTH) +
  geom_vline(aes(xintercept = historical_control_success_rate,
                 linetype = "Historical Control"), color = "black", size = LINEWIDTH) +
  labs(title="Type I Error", x="True Control Proportion", y="Type I Error") +
  # Set theme
  scale_color_viridis_d()  # Apply viridis color scale


# Plot Power
plot2 <-  ggplot(data, aes(x=control_proportion, y=Power, color=Design)) +
  geom_line(size=1) +
  geom_vline(aes(xintercept = historical_control_success_rate,
                 linetype = "Historical Control"), color = "black", size = LINEWIDTH) +
  labs(title="Power for 12 percent difference", x="True Control Proportion", y="Power") +
  # Set theme
  scale_color_viridis_d() +
  scale_x_continuous(breaks = c(seq(0, 1, by = 0.1), historical_control_success_rate), labels = c(seq(0, 1, by = 0.1), TeX("$p_S^{(c)}$")), expand = c(0, 0))

# Arrange the plots side by side using grid.arrange()
plt <- grid.arrange(plot1, plot2, ncol=2, widths=c(1, 1.4))


plot.size <- set_size(textwidth)
fig_width_in <- plot.size[1]
fig_height_in <- plot.size[2]

export_plots(
  plt = plt,
  file_path = "./figures/Viele_et_al_fisher",
  fig_width_in,
  fig_height_in,
  type = "pdf"
)
export_plots(
  plt = plt,
  file_path = "./figures/Viele_et_al_fisher",
  fig_width_in,
  fig_height_in,
  type = "png"
)
