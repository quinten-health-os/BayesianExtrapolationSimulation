
library(RBExT)
devtools::load_all()

source(system.file(paste0("conf/plots_config.R"), package = "RBExT"))

source(system.file(paste0("conf/metrics_config.R"), package = "RBExT"))

env <- "combined"

config_dir <- paste0(system.file(paste0("conf/", env), package = "RBExT"), "/")
source(paste0(config_dir, "methods_config.R"))

results_dir <- paste0("./results/", env)

tables_dir <<- paste0("./tables/", env, "_tables")

case_studies_config_dir <- paste0(system.file(paste0("conf/case_studies"), package = "RBExT"), "/")
source(system.file(paste0("conf/plots_config.R"), package = "RBExT"))
source(system.file(paste0("conf/methods_plots_config.R"), package = "RBExT"))
source(system.file(paste0("conf/metrics_config.R"), package = "RBExT"))

source(paste0(config_dir, "methods_config.R"))

results_freq_df <- readr::read_csv(paste0(results_dir, "/results_frequentist.csv"))
results_bayes_df <- readr::read_csv(paste0(results_dir, "/results_bayesian_simpson.csv"))

# results_freq_df <- subset(results_freq_df, case_study %in% c("botox"))
# results_bayes_df <- subset(results_bayes_df, case_study %in% c("belimumab"))

#results_freq_df <- subset(results_freq_df, method %in% c("commensurate_power_prior"))

results_freq_df <- results_freq_df %>%
  mutate_if(is.numeric, round, digits = 4)

results_bayes_df <- results_bayes_df  %>%
  mutate_if(is.numeric, round, digits = 4)


# noninflated_cases_df <- readr::read_csv(paste0(results_dir, "/noninflated_tie_cases.csv"))
# table_noninflated_tie(noninflated_cases_df)
#
# power_gain_cases_df <- readr::read_csv(paste0(results_dir, "/power_gain_cases.csv"))
# table_power_gain(power_gain_cases_df)
#
# power_loss_cases_df <- readr::read_csv(paste0(results_dir, "/power_loss_cases.csv"))
# table_power_loss(power_loss_cases_df)

power_loss_cases_inflated_tie_df <- readr::read_csv(paste0(results_dir, "/power_loss_inflated_tie_cases.csv"))
# table_power_loss_cases_inflated_tie(power_loss_cases_inflated_tie_df) #Done

power_vs_tie_tables(results_freq_df, frequentist_metrics)

# Make tables with columns method-parameters, consistent, partially consistent and no effect
table_methods_comparison(results_freq_df, frequentist_metrics)

table_bayesian_ocs(results_bayes_df, bayesian_metrics)

tables_metric_vs_scenario(results_freq_df, frequentist_metrics)

table_empirical_bayes_hyperparameters_vs_scenario(results_freq_df, frequentist_metrics)

posterior_parameters_tables(results_freq_df, frequentist_metrics)
