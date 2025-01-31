# Function to format numbers in scientific notation with 2 significant digits
format_num <- function(x, digits = 2, scientific = FALSE) {
  format(x, digits = digits, scientific = scientific)
}

#' Function to return markers from a list
#'
#' @param i The index of the marker to return.
#' @param markers_list A list of markers.
#' @return The marker at the specified index.
#' @keywords internal
markers <- function(i, markers_list) {
  return(markers_list[i %% length(markers_list)])
}

#' Function to set figure dimensions
#'
#' @param width The desired width of the figure (in pts).
#' @param fraction The fraction of the width to use (default is 1).
#' @param aspect_ratio The desired aspect ratio of the figure (default is golden ratio).
#' @return A vector containing the width and height of the figure.
#' @keywords internal
set_size <- function(width,
                     fraction = 1,
                     aspect_ratio = NULL) {
  # Width of figure (in pts)
  fig_width_pt <- width * fraction

  # Convert from pt to inches
  inches_per_pt <- 1 / 72.27

  # Golden ratio to set aesthetic figure height
  # https://disq.us/p/2940ij3
  golden_ratio <- (5 ^ 0.5 - 1) / 2

  if (is.null(aspect_ratio)) {
    aspect_ratio <- golden_ratio
  }

  # Figure width in inches
  fig_width_in <- fig_width_pt * inches_per_pt
  # Figure height in inches
  fig_height_in <- fig_width_in * aspect_ratio

  fig_dim <- c(fig_width_in, fig_height_in)

  return(fig_dim)
}


#' Function to format uncertainty
#'
#' @param yerr_input The input uncertainty values.
#' @param y The central values.
#' @param metric The metric information.
#' @return The formatted uncertainty.
#' @keywords internal
format_uncertainty <- function(yerr_input, y, metric) {
  # Check if metric uncertainty is a 95% CI or a MCSE
  if (grepl("conf_int|credible_interval", metric$metric_uncertainty)) {
    stopifnot(ncol(yerr_input) == 2)

    # Convert CI into error for error bars
    yerr <- yerr_input - y

    yerr[, 1] <- -yerr[, 1]

    # Check for negative errors and replace with 0
    close_to_zero <- abs(yerr) < 1e-15
    yerr[close_to_zero] <- 0

    # Check for negative errors
    if (any(yerr[, 1] < 0) || any(yerr[, 2] < 0)) {
      stop("yerr bounds are negative")
    }
  }

  # Return formatted uncertainty
  return(yerr)
}


#' Function to convert parameters in dataframe to string
#'
#' @param method The method name.
#' @param parameters The parameters dataframe.
#' @return The string representation of the parameters.
#' @keywords internal
convert_params_to_str <- function(method, parameters) {
  # Round floating point parameters to 2 decimal places
  parameters <- map_if(parameters, is.numeric, round, 2)

  if (is.null(names(parameters))) {
    stop("parameters must be a dataframe")
  }
  # Convert parameters to string representation
  labels <- paste(Filter(Negate(is.null), sapply(names(parameters), function(key) {
    if (length(method[[key]]$range) > 1) {
      paste0(method[[key]]$parameter_label, "=", parameters[[key]])
    }
  })), collapse = "_")
  return(labels)
}

#' Function to make labels from parameters. Return a label formatted in Tex, for example "$\\xi_\\gamma$ = 0.5, $\\sigma_\\gamma$ = 0.1"
#'
#' @param parameter The parameter dataframe.
#' @param method The method name.
#' @return The labels generated from the parameters.
#' @keywords internal
make_labels_from_parameters <- function(parameter, method) {
  # Round floating point parameters to 2 decimal places
  parameter[sapply(parameter, is.numeric)] <- lapply(parameter[sapply(parameter, is.numeric)], round, 2)

  # Select parameters with non-singleton ranges or display keyword = TRUE
  selection <- list()

  colnames(parameter) <- gsub("heterogeneity_prior\\.", "", colnames(parameter))
  if (method == "commensurate_power_prior"){
    if (is.null(parameter$family)){
      if (!is.null(parameter$alpha) && !is.na(parameter$alpha)){
        parameter$family = "inverse_gamma"
      } else if (!is.null(parameter$std_dev) && !is.na(parameter$std_dev)){
        parameter$family = "half_normal"
      } else if (!is.null(parameter$location) && !is.na(parameter$location)){
        parameter$family = "cauchy"
      }
    }
    if (parameter$family == "half_normal"){
      parameters_list <- list(heterogeneity_prior_family = parameter$family, std_dev = parameter$std_dev)
      label <- paste0("$\\tau \\sim HN(", round(as.numeric(parameters_list$std_dev),2), ")$")
    } else if (parameter$family == "inverse_gamma"){
      parameters_list <- list(heterogeneity_prior_family = parameter$family, alpha = parameter$alpha, beta = parameter$beta)
      label <- paste0("$\\tau \\sim IG(\\alpha = ",  round( as.numeric(parameters_list$alpha),2),", \\beta = ",  round( as.numeric(parameters_list$beta),2), ")$")
    } else if (parameter$family == "cauchy"){
      parameters_list <- list(heterogeneity_prior_family = parameter$family, location =  parameter$location, scale =  parameter$scale)
      label <- paste0("$\\tau \\sim Cauchy(x_0 = ",  round(as.numeric(parameter$location),2),"\\gamma = ", round( as.numeric(parameter$scale),2), ")$")
    } else {
      stop("Heterogeneity prior family not implemented.")
    }
  } else {
    for (parameter_name in names(parameter)) {
      if (is.null(methods_dict[[method]][[parameter_name]])){
        stop("Parameter not listed in the method configuration.")
      }

      if (!is.null(methods_dict[[method]][[parameter_name]]$display)) {
        selection <- c(selection, methods_dict[[method]][[parameter_name]]$display)
      } else {
        selection <- c(selection, length(methods_dict[[method]][[parameter_name]]$range) > 1)
      }
    }

    selection <- unlist(selection)

    label <- ""
    if (!is.null(selection)) {
      i <- 0
      parameters_to_display <- parameter[, selection, drop = FALSE]
      for (parameter_name in colnames(parameters_to_display)) {
        # Combine column names and values and convert selected parameters to label

        new_label <- sprintf(
          "%s = %s",
          methods_dict[[method]][[parameter_name]]$parameter_notation,
          as.character(round(as.numeric(parameter[[parameter_name]]),2))
        )

        if (i > 0) {
          label <- paste(label, new_label, sep = ", ")
        } else {
          label <- paste0(label, new_label)
        }

        i <- i + 1
      }
    }
  }

  return(label)
}



process_method_parameters_label <- function(row,
                                            methods_labels,
                                            method_name = TRUE,
                                            parameters_colname = "parameters",
                                            as_latex = TRUE){
  # Process the row of a results dataframe to create a Method + Parameters label
  method <- unlist(row["method"])

  if (is.null(methods_labels[[method]])) {
    stop(
      paste0(
        "Method label is not defined. Add a label in methods_config.R for method ",
        method
      )
    )
  }

  parameters_df <- get_parameters(row[parameters_colname])

  # Return a label formatted in Tex, for example "$\\xi_\\gamma$ = 0.5, $\\sigma_\\gamma$ = 0.1"
  param_label <- make_labels_from_parameters(parameter = parameters_df, method = method)


  if (method_name == TRUE) {
    # Add the method name to the label
    if (param_label == "") {
      full_label <- methods_labels[[method]]$label
    } else {
      full_label <- paste(methods_labels[[method]]$label, param_label, sep = ", ")
    }
  } else {
    if (param_label == "") {
      full_label <- ""
    } else {
      full_label <- param_label
    }
  }

  if (as_latex){
    full_label <- latex2exp::TeX(full_label)
  }

  if (is.null(full_label)) {
    stop("Label is null.")
  }

  return(full_label)
}


#' Function to return a dataframe of parameters from json strings
#'
#' @description This function returns, for a dataframe containing parameters (in json strings) for a given method, an R dataframe.
#'
#' @param parameters_df The dataframe containing parameters in JSON strings.
#' @return An R dataframe with parsed parameters.
#' @keywords internal
get_parameters <- function(parameters_df) {
  if ("posterior_parameters" %in% colnames(parameters_df)) {
    # Process the posterior parameters
    parameters_df$parameters <- parameters_df$posterior_parameters
    parameters_df <- parameters_df[, !(names(parameters_df) %in% c("posterior_parameters")), drop = FALSE]
  }

  if (is.null(parameters_df$parameters)){
    stop("The input parameters dataframe does not contain parameters.")
  }

  if (all(unique(parameters_df$parameters) == "[]")) {
    return(parameters_df)
  }

  N <- nrow(parameters_df)

  # Check that the input is a dataframe
  if (!is.data.frame(parameters_df)) {
    if (length(parameters_df) == 1) {
      parameters_df <- gsub("'", "\"", parameters_df)
      json_parameters_df <- jsonlite::fromJSON(parameters_df, simplifyDataFrame = TRUE)
      if (length(names(json_parameters_df)) == 1 &&
          (json_parameters_df) == "parameters") {
        json_parameters_df <- unlist(unlist(json_parameters_df$parameters))
      } else {
        json_parameters_df <- unlist(json_parameters_df)
      }
      return(json_parameters_df)
    } else {
      stop("Input is not a dataframe")
    }
  } else {
    # Replace single quotes with double quotes
    parameters_df$parameters <- gsub("'", "\"", parameters_df$parameters)
  }

  results_list <- list()

  if (nrow(parameters_df) > 0) {
    # Loop through each row in the dataframe
    for (i in 1:nrow(parameters_df)) {
      json_parameters_df <- data.frame(jsonlite::fromJSON(parameters_df$parameters[i], simplifyDataFrame = TRUE))

      if (is.null(names(json_parameters_df))) {
        parameter_values <- NA
      } else if (length(names(json_parameters_df)) == 1 &&
                 names(json_parameters_df) == "parameters") {
        parameter_values <- unlist(unlist(json_parameters_df$parameters))
      } else {
        parameter_values <- unlist(json_parameters_df)
      }

      results_list[[i]] <- parameter_values
    }
  }

  # Combine the list of data frames into a single data frame, removing any NULL entries
  results_list <- results_list[!sapply(results_list, is.null)]

  # Get the number of columns for each element in the list
  num_cols <- sapply(results_list, length)

  # Check if all elements have the same number of columns
  if (length(results_list) != 0 && length(unique(num_cols)) != 1) {

    # Step 1: Extract all unique column names across elements in the list
    all_columns <- unique(unlist(lapply(results_list, names)))

    # Step 2: Ensure each element has the same columns by adding missing ones as NA
    results_list_aligned <- lapply(results_list, function(x) {
      missing_cols <- setdiff(all_columns, names(x))
      x[missing_cols] <- NA
      x <- x[all_columns]  # Reorder to keep column order consistent
      return(x)
    })

    # Step 3: Bind rows without error
    results_list <- results_list_aligned

    # Now, results_df will have consistent columns

    # stop(
    #   "Error: The elements in the list have different numbers of columns. The function should only be applied to a dataframe containing parameters for a single method."
    # )
  }

  results_df <- do.call(rbind, results_list)

  if (!is.data.frame(results_df)) {
    results_df <- data.frame(results_df)
  }

  if (nrow(results_df) != N) {
    stop("Output does not have the same number of rows as the input df.")
  }

  return(results_df)
}

# Define a function to convert strings to numeric or boolean if possible
convert_if_possible <- function(x) {
  output <- vector("list", length(x)) # Initialize output as a list
  for (i in seq(length(x))) {
    if (is.character(x[i])) {
      # Try to convert to numeric
      numeric_value <- suppressWarnings(as.numeric(x[i]))
      if (!is.na(numeric_value)) {
        output[[i]] <- numeric_value
        next
      }

      # Try to convert to boolean
      if (tolower(x[i]) == "true") {
        output[[i]] <- TRUE
      } else if (tolower(x[i]) == "false") {
        output[[i]] <- FALSE
      } else {
        output[[i]] <- x[i] # Keep as character if no conversion occurs
      }
    }
  }
  # Return the original value if no conversion was possible
  return(unlist(output))
}

export_plots <- function(plt,
                         file_path,
                         fig_width_in,
                         fig_height_in,
                         type = "pdf", forest_plot = FALSE, adjust_theme = TRUE) {

  if (!(forest_plot) && adjust_theme == TRUE){
    plt <- plt + theme_bw() + theme(
      axis.text = element_text(family = font, size = text_size/2),
      axis.text.y = element_text(family = font, size = small_text_size/2),
      axis.text.x = element_text(family = font, size = small_text_size/2),
      axis.title = element_text(family = font, size = text_size/2),
      plot.title = element_text(family = font, size = text_size/2),
      legend.text = element_text(family = font, size = small_text_size/2),
      legend.title = element_text(family = font, size = text_size/2),
      legend.key = element_blank(),
      # strip.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  }

  if (type == "pdf") {
    ggplot2::ggsave(
      filename = paste0(file_path, ".pdf"),
      plot = plt,
      device = cairo_pdf,
      width = fig_width_in,
      height = fig_height_in,
      dpi = dpi
    )
  } else {
    ggplot2::ggsave(
      filename = paste0(file_path, ".png"),
      plot = plt,
      width = fig_width_in,
      height = fig_height_in,
      dpi = dpi
    )
  }

}

format_title <- function(title, case_study, target_to_source_std_ratio = NA, source_denominator_change_factor = NA, as_latex = FALSE){
  if (case_study %in% c("botox", "dapagliflozin") & !is.na(target_to_source_std_ratio)){
    title <- paste0(title, ", $\\sigma_T/\\sigma_S = $", target_to_source_std_ratio)
  }

  if (!is.na(source_denominator_change_factor)){
    if (!(case_study %in% c("botox", "dapagliflozin", "aprepitant"))){
      title = paste0(title, ", Source denominator change factor = ",source_denominator_change_factor)
    }
  }

  if (as_latex == FALSE){
    title <- latex2exp::TeX(title)
  }

  return(title)
}

format_filename <- function(filename, case_study, target_to_source_std_ratio = NA, source_denominator_change_factor = NA){
  if (case_study %in% c("botox", "dapagliflozin")){
    filename <- paste0(filename,
                       "_target_to_source_std_ratio=",
                       target_to_source_std_ratio
    )
  }

  if (!is.na(source_denominator_change_factor)){
    if (!(case_study %in% c("botox", "dapagliflozin", "aprepitant"))){
      filename <- paste0(filename, "_source_denominator_change_factor=",  source_denominator_change_factor)
    }
  }

  # Convert to lowercase and replace spaces with underscores
  filename <- tolower(filename)
  filename <- gsub(" ", "_", filename)

  return(filename)
}

format_results_df_parameters <- function(results_df, include_method_name = TRUE){
  labels <- c()
  for (i in seq_len(nrow(results_df))) {
    method <- unlist(results_df[i, "method"])

    parameters_df <- get_parameters(results_df[i, "parameters", drop = FALSE])

    param_label <- make_labels_from_parameters(parameter = parameters_df, method = method)

    if (include_method_name){
      if (param_label == "") {
        full_label <- methods_labels[[method]]$label
      } else {
        full_label <- paste(methods_labels[[method]]$label, param_label, sep = ", ")
      }
    } else {
      full_label <- param_label
    }
    # Append the label to the labels vector
    labels <- c(labels, full_label)
  }
  return(labels)
}

format_results_df_methods <- function(results_df){
  methods <- c()
  for (i in seq_len(nrow(results_df))) {
    method <- unlist(results_df[i, "method"])
    methods <- c(methods, methods_labels[[method]]$label)
  }
  return(methods)
}


convert_CI_columns <- function(table_data_df) {
  # Get the names of all columns that start with "conf_int"
  conf_int_cols <- grep("^conf_int", colnames(table_data_df), value = TRUE)

  # Extract unique parameter names (e.g., "power_parameter_mean")
  unique_params <- unique(sub("conf_int_(lower|upper)_", "", conf_int_cols))

  # Dynamically combine parameter values with their confidence intervals
  for (param in unique_params) {
    # Find the lower and upper confidence interval columns
    lower_col <- paste0("conf_int_lower_", param)
    upper_col <- paste0("conf_int_upper_", param)

    # Check if the parameter column exists
    if (param %in% colnames(table_data_df)) {
      # Combine the value and confidence interval into a single string
      table_data_df[[param]] <- paste0(
        table_data_df[[param]], " [", table_data_df[[lower_col]], ", ", table_data_df[[upper_col]], "]"
      )
    } else {
      # If the parameter column doesn't exist, create a new one with just the CI
      table_data_df[[param]] <- paste0(
        "[", table_data_df[[lower_col]], ", ", table_data_df[[upper_col]], "]"
      )
    }
  }
  # Drop the original confidence interval columns
  cols_to_keep <- !grepl("^conf_int", colnames(table_data_df))
  table_data_df <- table_data_df[, cols_to_keep, drop = FALSE]
  return(table_data_df)
}


com_pp_params_filtering = function(key, parameter){ # TODO : remove ?
  if (is.null(parameter$family)){
    if ((!is.null(parameter$alpha) && !is.na(parameter$alpha))){
      parameter$family = "inverse_gamma"
    } else if (!is.null(parameter$std_dev) && !is.na(parameter$std_dev)){
      parameter$family = "half_normal"
    } else if (!is.null(parameter$location) && !is.na(parameter$location)){
      parameter$family = "cauchy"
    }

    if ((!is.null(parameter$heterogeneity_prior.alpha) && !is.na(parameter$heterogeneity_prior.alpha))){
      parameter$family = "inverse_gamma"
    } else if (!is.null(parameter$heterogeneity_prior.std_dev) && !is.na(parameter$heterogeneity_prior.std_dev)){
      parameter$family = "half_normal"
    } else if (!is.null(parameter$heterogeneity_prior.location) && !is.na(parameter$heterogeneity_prior.location)){
      parameter$family = "cauchy"
    }
  }

  if (parameter$family == "half_normal"){
    parameters_list <- list(heterogeneity_prior_family = parameter$family, std_dev = parameter$std_dev)
    label <- paste0("$\\tau \\sim HN(", round(as.numeric(parameters_list$std_dev),2), ")$")
  } else if (parameter$family == "inverse_gamma"){
    parameters_list <- list(heterogeneity_prior_family = parameter$family, alpha = parameter$alpha, beta = parameter$beta)
    label <- paste0("$\\tau \\sim IG(\\alpha = ",  round( as.numeric(parameters_list$alpha),2),", \\beta = ",  round( as.numeric(parameters_list$beta),2), ")$")
  } else if (parameter$family == "cauchy"){
    parameters_list <- list(heterogeneity_prior_family = parameter$family, location =  parameter$location, scale =  parameter$scale)
    label <- paste0("$\\tau \\sim Cauchy(x_0 = ",  round(as.numeric(parameter$location),2),"\\gamma = ", round( as.numeric(parameter$scale),2), ")$")
  } else {
    stop("Heterogeneity prior family not implemented.")
  }

  methods_dict[[method]][[key]][["range"]]

  return(filter)
}
