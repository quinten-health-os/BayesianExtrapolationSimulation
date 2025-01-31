# Required Libraries
library(ggplot2)
library(reshape2)

library(RBExT)

# Parameters for the simulation
historical_control_success_rate <- 0.65  # Historical control success rate
effect_size <- 0.12  # 12% effect size for treatment arm
n_trials <- 5000  # Number of simulations
sample_size <- 200  # Sample size per arm
control_proportions <- seq(0.50, 0.80, length.out = 30)  # Range of true control proportions
treatment_proportions <- seq(0.50, 0.80, length.out = 30)  # Range of true treatment proportions

# Containers for results
psuccess_diff <- matrix(NA, nrow = length(control_proportions), ncol = length(treatment_proportions))

# historical trial control data
historical_control_data <- c(rep(0,35), rep(1, 65))

# Simulation loop
for (i in seq_along(control_proportions)) {
  for (j in seq_along(treatment_proportions)) {

    true_control <- control_proportions[i]
    true_treatment <- treatment_proportions[j]

    test_result_separate <- c()
    test_result_pooled <- c()

    for (k in 1:n_trials) {

      # Generate control and treatment arm data
      control_data <- rbinom(sample_size, 1, true_control)
      treatment_data <- rbinom(sample_size, 1, true_treatment)

      p_treatment <- mean(treatment_data)

      # Separate analysis
      p_control_separate <- mean(control_data)

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

      # Pooled analysis
      pooled_control_data <- c(control_data, historical_control_data)
      p_control_pooled <- mean(pooled_control_data)

      pooling <- BSDA::tsum.test(
        mean.x = p_treatment,
        mean.y = p_control_pooled,
        mu = 0,
        alternative = "greater",
        s.x = sqrt(p_treatment * (1 - p_treatment)),
        s.y = sqrt(p_control_pooled * (1 - p_control_pooled)),
        n.x = sample_size,
        n.y = length(pooled_control_data)
      )

      test_result_pooled <- c(test_result_pooled, pooling$p.value < 0.05)
    }

    # Calculate Pr(Success) for each method
    psuccess_separate <- mean(test_result_separate)
    psuccess_pooled <- mean(test_result_pooled)

    # Store the difference in success probability between pooled and separate methods
    psuccess_diff[i, j] <- psuccess_pooled - psuccess_separate
  }
}

# Convert to data frame for plotting
diff_data <- expand.grid(control_proportion = control_proportions, treatment_proportion = treatment_proportions)
diff_data$diff_psuccess <- as.vector(psuccess_diff)


# Plot the heatmap with diagonal highlight
heatmap_plot <- ggplot(diff_data, aes(x = control_proportion, y = treatment_proportion, fill = diff_psuccess)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  # Color gradient
  geom_vline(xintercept = historical_control_success_rate, linetype = "dashed", color = "black", size = 0.5) +  # Add vline for historical control
  labs(title = "Difference in Success Probability (Pooled vs Separate)",
       x = "True Control Proportion",
       y = "True Treatment Proportion",
       fill = "Difference in\nPr(Success)") +  # Label for color bar
  theme_minimal() +
  theme(
    text = element_text(family = font, size = text_size),
    axis.text.y = element_text(size = text_size, margin = margin(r = 0)),
    axis.text.x = element_text(size = text_size, margin = margin(t = 0)),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.ticks = element_line(color = "black"),  # Customize tick appearance
    axis.ticks.length = unit(0.3, "cm"),  # Control the length of the ticks
     axis.line = element_line(color = "black"),  # Optional: Add axis lines for clarity

    legend.position = "right") +  # Add colorbar on the right side
  scale_x_continuous(breaks = c(seq(0, 1, by = 0.1), 0.65), labels = c(seq(0, 1, by = 0.1), TeX("$p_S^{c}= 0.65$")), expand = c(0, 0)) +

  scale_y_continuous(
    expand = c(0, 0)  # Remove padding between y-axis and heatmap
  )  +
  # Highlight the diagonal line (TIE)
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray", size = 0.5) +

  # Add label for the diagonal
  annotate("text", x = 0.8, y = 0.8, label = "TIE", color = "black", size = 3, angle = 0)

# Display the plot
print(heatmap_plot)


export_plots(
  plt = heatmap_plot,
  file_path = "./figures/Viele_et_al_variant_heatmap",
  fig_width_in,
  fig_height_in,
  type = "pdf"
)
export_plots(
  plt = heatmap_plot,
  file_path = "./figures/Viele_et_al_variant_heatmap",
  fig_width_in,
  fig_height_in,
  type = "png"
)
