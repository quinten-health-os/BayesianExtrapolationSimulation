export_table <- function(data_table, ncollapses = 1, column_names = NULL, title, file_path, longtable = TRUE, remake_tables = TRUE, page_width = 5) {
  # Parameters:
  # page_width: Page width threshold in inches (default 6.5 inches for US Letter with 1-inch margins)

  if (file.exists(paste0(file_path, ".tex")) &&
      file.exists(paste0(file_path, ".pdf")) &&
      remake_tables == FALSE) {
    return()
  }

  if (is.null(column_names)) {
    column_names <- colnames(data_table)
  }

  # Escape '%' symbols in column names for LaTeX compatibility
  column_names <- gsub("%", "\\\\%", column_names)
  column_names <- lapply(column_names, function(x) paste0("\\textbf{", x, "}"))
  title <- gsub("%", "\\\\%", title)

  # Estimate the total width of the table
  avg_char_width <- 0.15  # Approximate width of a character in inches (adjust as needed)
  column_widths <- sapply(column_names, function(name) nchar(as.character(name)) * avg_char_width)
  total_table_width <- sum(column_widths)

  # Check if the table is too wide for the page
  use_landscape <-FALSE # total_table_width > page_width

  if (use_landscape){
    longtable = FALSE
  }

  # Generate alignment specification for centering columns
  align_spec <- rep("c", ncol(data_table))  # 'c' for centered alignment in LaTeX

  # Generate LaTeX table using kable and kableExtra
  data_kable <- data_table %>%
    knitr::kable(
      format = "latex",
      booktabs = TRUE,
      longtable = longtable,
      col.names = column_names,
      caption = title,
      escape = FALSE,
      digits = 3,
      align = align_spec
    ) %>%
    kableExtra::collapse_rows(columns = 1:ncollapses, valign = "middle") %>%
    kableExtra::kable_styling(
      latex_options = c("striped", "hold_position", "scale_down"),
      position = "center"
    )

  tex_file <- paste0(file_path, ".tex")

  # Write LaTeX code to file with the correct document structure
  cat(
    "\\documentclass{article}\n",
    "\\usepackage{geometry}\n",
    "\\geometry{margin=0.5in}\n",  # Adjusted margins to 0.5 inches
    "\\usepackage{longtable}\n",
    "\\usepackage{booktabs}\n",
    "\\usepackage{multirow}\n",
    "\\usepackage{graphicx}\n",
    "\\usepackage[table]{xcolor}\n",
    "\\usepackage{pdflscape}\n",  # Include pdflscape package
    "\\pagestyle{empty}\n",  # No page numbering
    "\\begin{document}\n",
    if (use_landscape) "\\begin{landscape}\n" else "",
    "\\centering\n",
    data_kable,
    if (use_landscape) "\n\\end{landscape}\n" else "",
    "\n\\end{document}\n",
    file = tex_file
  )

  # Compile the .tex file into a PDF
  tinytex::latexmk(tex_file)

  pdf_file <- paste0(file_path, ".pdf")

  # Safeguard against incorrect cropping
  if (file.exists(pdf_file)) {
    pdfcrop_command <- paste("pdfcrop", shQuote(pdf_file), shQuote(pdf_file))
    crop_result <- system(pdfcrop_command, intern = TRUE)
    if (length(crop_result) == 0) {
      warning("pdfcrop execution failed or returned no output. The PDF file might not be cropped correctly.")
    }
  } else {
    stop("The PDF file was not created successfully. Please check for LaTeX compilation errors.")
  }
}
