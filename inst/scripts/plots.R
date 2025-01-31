library(RBExT)
library(tidyverse)
library(readr)
library(ggplot2)
library(extrafont)
library(tikzDevice)
library(grid)
library(ggnewscale)
library(dplyr)
library(purrr)
library(tibble)

devtools::load_all()

env <- "combined"

results_dir <- paste0("./results/", env)
config_dir <- paste0(system.file(paste0("conf/", env), package = "RBExT"), "/")
case_studies_config_dir <- paste0(system.file(paste0("conf/case_studies"), package = "RBExT"), "/")

figures_dir <<- paste0("./figures/", env, "_figures/")

source(system.file(paste0("conf/plots_config.R"), package = "RBExT"))
source(system.file(paste0("conf/methods_plots_config.R"), package = "RBExT"))
source(system.file(paste0("conf/metrics_config.R"), package = "RBExT"))

source(paste0(config_dir, "methods_config.R"))

results_freq_df <- readr::read_csv(paste0(results_dir, "/results_frequentist.csv"))
results_bayes_df <- readr::read_csv(paste0(results_dir, "/results_bayesian_simpson.csv"))
sweet_spot_df <-  read_json(paste0(results_dir, "/sweet_spot.json"), simplifyVector = TRUE)

# Filter the results with the methods we are interest in:
# results_freq_df <- subset(results_freq_df, method %in% c("conditional_power_prior", "EB_PP", "p_value_based_PP", "PDCCPP", "pooling", "RMP", "separate", "test_then_pool_difference", "test_then_pool_equivalence", "NPP", "commensurate_power_prior"))

# results_bayes_df <- subset(results_bayes_df, method %in% c("NPP"))
#
# results_bayes_df <- subset(results_bayes_df, case_study %in% c("aprepitant"))
# #
# #
#results_freq_df <- subset(results_freq_df, case_study %in% c("belimumab"))
# results_freq_df  <- subset(results_freq_df,  case_study %in% c("teriflunomide"))
#results_freq_df <- subset(results_freq_df, method %in% c("separate", "pooling", "commensurate_power_prior"))

#results_freq_df <- subset(results_freq_df, case_study %in% c("botox"))

# results_bayes_df <- subset(results_bayes_df, case_study %in% c("botox"))
# results_freq_df <- subset(results_freq_df, case_study %in% c("belimumab"))
#sweet_spot_df <- subset(sweet_spot_df, case_study %in% c("botox"))

#sweet_spot_df <- subset(sweet_spot_df, method %in% c("conditional_power_prior"))

remake_figures <<- FALSE


# loop_plot_effect_target_standard_deviation(results_freq_df, frequentist_metrics)
# loop_plot_effect_denominator_change(results_freq_df, frequentist_metrics)
#
# loop_plot_effect_target_standard_deviation(results_freq_df, inference_metrics)
# loop_plot_effect_denominator_change(results_freq_df, inference_metrics)

# loop_plot_effect_target_standard_deviation(results_bayes_df, bayesian_metrics)
# loop_plot_effect_denominator_change(results_bayes_df, bayesian_metrics)



# plot_metric_vs_scenario_methods(results_freq_df, inference_metrics) # DONE
# plot_metric_vs_scenario_methods(results_freq_df, frequentist_metrics)# DONE
#
# forest_plot_methods_comparison(results_freq_df, inference_metrics)# DONE
# forest_plot_methods_comparison(results_freq_df, frequentist_metrics)# DONE
#
# plot_success_proba_vs_scenario(results_freq_df, frequentist_metrics)# DONE

# plot_metric_vs_scenario(results_freq_df, inference_metrics)
# plot_metric_vs_scenario(results_freq_df, frequentist_metrics)
plot_metrics_vs_ess(results_freq_df, frequentist_metrics, inference_metrics)

empirical_bayes_parameters_plots(results_freq_df)

posterior_parameters_plots(results_freq_df)

power_vs_tie_plots(results_freq_df, frequentist_metrics)

operating_characteristics_vs_tie_plots(results_freq_df, frequentist_metrics)
operating_characteristics_vs_tie_plots(results_freq_df, inference_metrics)
#
#
#
# bayesian_operating_characteristics_vs_tie_plots(results_bayes_df, results_freq_df, bayesian_metrics) # DONE
#
# forest_plot_methods_comparison_bayesian(results_bayes_df, bayesian_metrics)# DONE
# forest_plot_sweet_spots_comparison(sweet_spot_df, frequentist_metrics)# DONE
# forest_plot_sweet_spots_comparison(sweet_spot_df, frequentist_metrics)# DONE
# plot_sweet_spot_width_vs_scenario(sweet_spot_df, sweet_spots_metrics)# DONE
# bayesian_ocs_plots(results_bayes_df, bayesian_metrics)# DONE
