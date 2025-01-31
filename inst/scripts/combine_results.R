library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(ggnewscale)
# Define the folder path (adjust based on your environment)
folder_path <- paste0(system.file("results/simulation_runs", package = "RBExT"), "/")

# Function to extract date from folder name
extract_date <- function(folder_name) {
  str_extract(folder_name, "\\d{2}_\\d{2}")
}

# List all folders containing "results_frequentist.csv"
folders <- list.dirs(folder_path, recursive = FALSE)

# Read and combine all CSV files
all_data <- folders %>%
  map_df(~ {
    # Construct file path
    csv_file <- file.path(.x, "results_frequentist.csv")

    # Read CSV if it exists
    if (file.exists(csv_file)) {
      # Extract date from the folder name
      folder_date <- extract_date(basename(.x))

      # Read the CSV file and add the extracted date
      read_csv(csv_file, col_types = cols()) %>%
        mutate(date = as.Date(folder_date, format = "%d_%m"))
    } else {
      NULL
    }
  })

# Function to round numeric columns to 6 decimal places
round_numeric_columns <- function(df) {
  df %>%
    mutate(across(where(is.numeric), ~ round(.x, 8)))
}

# Apply rounding to all numeric columns
all_data <- round_numeric_columns(all_data)

# Define the key columns for identifying duplicates
key_columns <- c(
  "method",
  "parameters",
  "control_drift",
  "source_denominator",
  "source_denominator_change_factor",
  "case_study",
  "target_to_source_std_ratio",
  "target_sample_size_per_arm",
  "theta_0",
  "null_space",
  "sampling_approximation",
  "summary_measure_likelihood",
  "source_sample_size_treatment",
  "source_sample_size_control",
  "endpoint",
  "source_standard_error",
  "source_treatment_effect_estimate",
  "equivalent_source_sample_size_per_arm",
  "target_treatment_effect",
  "target_sample_size_per_arm"
)

# Deduplicate rows, keeping only the most recent entry for each combination of key columns
deduplicated_data <- all_data %>%
  arrange(desc(date)) %>%
  distinct(across(all_of(key_columns)), .keep_all = TRUE)

# View the final result
print(deduplicated_data)

# Optionally, save the concatenated and deduplicated data to a new CSV
output_file <- paste0(folder_path, "/combined_results.csv")
write_csv(deduplicated_data, output_file)

cat("Deduplication complete. Results saved to:", output_file, "\n")



# Check for duplicates based on key columns
duplicates_count <- deduplicated_data %>%
  group_by(across(all_of(key_columns))) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  nrow()

# Display the total number of duplicates
cat("Number of duplicate rows based on scenarios columns:", duplicates_count, "\n")



#
# # List all folders containing "results_frequentist.csv"
# folders <- list.dirs(folder_path, recursive = FALSE)
#
# # Initialize an empty vector to store folder names
# folders_with_column <- c()
#
# # Function to check if a CSV file contains the specified column
# check_column_presence <- function(folder) {
#   csv_file <- file.path(folder, "results_frequentist.csv")
#   if (file.exists(csv_file)) {
#     data <- read_csv(csv_file, col_types = cols())
#     if ("conf_int_post_mean_upper" %in% colnames(data)) {
#       return(TRUE)
#     }
#   }
#   return(FALSE)
# }
#
# # Iterate over each folder and check for the column
# for (folder in folders) {
#   if (check_column_presence(folder)) {
#     folders_with_column <- c(folders_with_column, basename(folder))
#   }
# }
#
# # Display the list of folders with the specified column
# print("Folders containing the column 'conf_int_post_mean_upper':")
# print(folders_with_column)
