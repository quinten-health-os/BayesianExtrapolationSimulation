library(RBExT)

devtools::load_all()
results_dir <- "./results/combined"
results_freq_df <- readr::read_csv(paste0(results_dir, "/results_frequentist.csv"))

results_freq_df <- results_freq_df %>%
  mutate_if(is.numeric, round, digits = 4)

#(results_freq_df, output_path = results_dir)
#analyze_power_loss(results_freq_df, output_path = results_dir)
analyze_power_loss_inflated_tie(results_freq_df, output_path = results_dir)
analyze_noninflated_tie(results_freq_df, output_path = results_dir)
