library(rmarkdown)
library(fs) # for file manipulation functions

# Get a list of all RMD files in the directory
rmd_files <- list.files(pattern = "\\.Rmd$", recursive = TRUE)

# Render each RMD file to HTML in its own directory
for (rmd_file in rmd_files) {
  tryCatch(
    {
      # Get the directory of the RMD file
      output_dir <- dirname(rmd_file)

      # Render the RMD file to HTML
      render(rmd_file,
        output_format = "html_document",
        output_dir = output_dir
      )

      cat("Knitting successful:", rmd_file, "\n")
    },
    error = function(e) {
      cat("Knitting failed for:", rmd_file, "\n")
      cat("Error message:", conditionMessage(e), "\n")
      # Continue to the next iteration of the loop
      next
    }
  )
}
