# Required Libraries
library(ggplot2)
library(BSDA)  # for tsum.test

# Parameters for the simulation
historical_control_mean <- 0.65  # Historical control mean
historical_control_sd <- 1  # Historical control standard deviation
effect_size <- 0.2
n_trials <- 10000  # Number of simulations
sample_size <- 200  # Sample size per arm
control_means <- seq(0.3, 1.00, length.out = 40)  # Range of true control means

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

# Historical control data (Normally distributed)
historical_control_data <- rnorm(100, mean = historical_control_mean, sd = historical_control_sd)

# Simulation loop
for (true_control_mean in control_means) {

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

    # Generate control and treatment arm data (Normally distributed)
    control_data <- rnorm(sample_size, mean = true_control_mean, sd = historical_control_sd)
    null_treatment_data <- rnorm(sample_size, mean = true_control_mean, sd = historical_control_sd)
    alt_treatment_data <- rnorm(sample_size, mean = true_control_mean + effect_size, sd = historical_control_sd)

    p_treatment_null <- mean(null_treatment_data)
    p_treatment_alt <- mean(alt_treatment_data)

    # Separate analysis
    p_control_separate <- mean(control_data)
    estimates_separate <- c(estimates_separate, p_control_separate)

    separate_null <- BSDA::tsum.test(
      mean.x = p_treatment_null,
      mean.y = p_control_separate,
      mu = 0,
      alternative = "greater",
      s.x = sd(null_treatment_data),
      s.y = sd(control_data),
      n.x = sample_size,
      n.y = sample_size
    )

    reject_null_separate <- c(reject_null_separate, separate_null$p.value < 0.05)

    separate_alt <- BSDA::tsum.test(
      mean.x = p_treatment_alt,
      mean.y = p_control_separate,
      mu = 0,
      alternative = "greater",
      s.x = sd(alt_treatment_data),
      s.y = sd(control_data),
      n.x = sample_size,
      n.y = sample_size
    )

    power_detected_separate <- c(power_detected_separate, separate_alt$p.value < 0.05)

    # Pooled analysis (include historical control data)
    pooled_control_data <- c(control_data, historical_control_data)
    p_control_pooled <- mean(pooled_control_data)
    estimates_pooled <- c(estimates_pooled, p_control_pooled)

    pooling_null <- BSDA::tsum.test(
      mean.x = p_treatment_null,
      mean.y = p_control_pooled,
      mu = 0,
      alternative = "greater",
      s.x = sd(null_treatment_data),
      s.y = sd(pooled_control_data),
      n.x = sample_size,
      n.y = length(pooled_control_data)
    )

    reject_null_pooled <- c(reject_null_pooled, pooling_null$p.value < 0.05)

    pooling_alt <- BSDA::tsum.test(
      mean.x = p_treatment_alt,
      mean.y = p_control_pooled,
      mu = 0,
      alternative = "greater",
      s.x = sd(alt_treatment_data),
      s.y = sd(pooled_control_data),
      n.x = sample_size,
      n.y = length(pooled_control_data)
    )

    power_detected_pooled <- c(power_detected_pooled, pooling_alt$p.value < 0.05)

    # Single-arm analysis (compare treatment arm to historical control)
    p_control_single <- mean(historical_control_data)
    estimates_single <- c(estimates_single, p_control_single)

    # Type I error and power
    single_null <- BSDA::tsum.test(
      mean.x = p_treatment_null,
      alternative = "greater",
      mu = p_control_single,
      s.x = sd(null_treatment_data),
      n.x = sample_size
    )

    single_alt <- BSDA::tsum.test(
      mean.x = p_treatment_alt,
      alternative = "greater",
      mu = p_control_single,
      s.x = sd(alt_treatment_data),
      n.x = sample_size
    )

    reject_null_single <- c(reject_null_single, single_null$p.value < 0.05)
    power_detected_single <- c(power_detected_single, single_alt$p.value < 0.05)
  }

  # Calculate MSE for each design
  mse_separate <- c(mse_separate, mean((estimates_separate - true_control_mean)^2))
  mse_pooled <- c(mse_pooled, mean((estimates_pooled - true_control_mean)^2))
  mse_single <- c(mse_single, mean((estimates_single - true_control_mean)^2))

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
  control_mean = rep(control_means, times=3),
  MSE = c(mse_separate, mse_pooled, mse_single),
  Type1Error = c(type1_error_separate, type1_error_pooled, type1_error_single),
  Power = c(power_separate, power_pooled, power_single),
  Design = rep(c("Separate", "Pooled", "Single Arm"), each=length(control_means))
)

# Plot MSE
ggplot(data, aes(x=control_mean, y=MSE, color=Design)) +
  geom_line(size=1) +
  labs(title="Mean Squared Error (MSE)", x="True Control mean", y="MSE of Control Estimate", color = "") +
  theme_minimal()

# Plot Type I Error
plot1 <- ggplot(data, aes(x=control_mean, y=Type1Error, color=Design)) +
  geom_line(size=1) +
  # Add hline with a legend
  geom_hline(aes(yintercept = 0.05,
                 linetype = "Significance Threshold (0.05)"), color = "black", size = LINEWIDTH) +
  geom_vline(aes(xintercept = historical_control_success_rate,
                 linetype = "Historical Control"), color = "black", size = LINEWIDTH) +
  labs(title="Type I Error", x="True Control mean", y="Type I Error", linetype = "") +
  # Set theme
  scale_color_viridis_d() + # Apply viridis color scale
  theme_bw() +
  theme(
    text = element_text(family = font, size = text_size),
    axis.text.y = element_text(size = text_size),
    axis.text.x = element_text(size = text_size),
    legend.key = element_blank(),
    # strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none" # Remove the legend
  ) +
  scale_x_continuous(breaks = c(seq(0, 1, by = 0.1), historical_control_success_rate), labels = c(seq(0, 1, by = 0.1), TeX("$p_S^{(c)}$")), expand = c(0, 0))


# Plot Power
plot2 <-  ggplot(data, aes(x=control_mean, y=Power, color=Design)) +
  geom_line(size=LINEWIDTH) +
  geom_vline(aes(xintercept = historical_control_success_rate,
                 linetype = "Historical Control"), color = "black", size = LINEWIDTH) +
  labs(title="Power for an effect size of 0.2", x="True Control mean", y="Power", linetype = "", color = "") +
  # Set theme
  scale_color_viridis_d() + # Apply viridis color scale
  theme_bw() +
  theme(
    text = element_text(family = font, size = text_size),
    axis.text.y = element_text(size = text_size),
    axis.text.x = element_text(size = text_size),
    legend.key = element_blank(),
    # strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Arrange the plots side by side using grid.arrange()
plt <- grid.arrange(plot1, plot2, ncol=2, widths=c(1, 1.4))


plot.size <- set_size(textwidth)
fig_width_in <- plot.size[1]
fig_height_in <- plot.size[2]

export_plots(
  plt = plt,
  file_path = "./figures/Viele_et_al_ttest_gaussian",
  fig_width_in,
  fig_height_in,
  type = "pdf"
)
export_plots(
  plt = plt,
  file_path = "./figures/Viele_et_al_ttest_gaussian",
  fig_width_in,
  fig_height_in,
  type = "png"
)
