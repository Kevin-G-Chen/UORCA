#!/usr/bin/env Rscript

# download_metadata.R
# Description: Downloads GEO metadata for a given GEO accession number and saves it as CSV files.

# Load required libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(GEOquery))
suppressPackageStartupMessages(library(tidyverse))

# Define command-line options
option_list <- list(
  make_option(c("-g", "--geo_accession"), type = "character", default = NULL,
              help = "GEO accession number (e.g., GSE12345)", metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character", default = ".",
              help = "Directory where the metadata files will be saved", metavar = "character")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Function to log messages
log_message <- function(message, log_file) {
  cat(paste0(message, "\n"))
  if (!is.null(log_file)) {
    write(paste0(Sys.time(), " - ", message), file = log_file, append = TRUE)
  }
}

# Main execution
main <- function() {
  # Check if GEO accession is provided
  if (is.null(opt$geo_accession)) {
    print_help(opt_parser)
    stop("Error: GEO accession must be provided.", call. = FALSE)
  }
  
  # Ensure output directory exists
  if (!dir.exists(opt$output_dir)) {
    log_message(paste("Creating output directory:", opt$output_dir), NULL)
    dir.create(opt$output_dir, recursive = TRUE)
  } else {
    log_message(paste("Using existing output directory:", opt$output_dir), NULL)
  }
  
  # Download GEO metadata
  log_message(paste("Retrieving metadata for GEO accession:", opt$geo_accession), NULL)
  gds <- tryCatch({
    getGEO(opt$geo_accession, GSEMatrix = TRUE, AnnotGPL = TRUE)
  }, error = function(e) {
    stop(paste("Error: Failed to retrieve GEO data for accession", opt$geo_accession, "\n", e))
  })
  
  if (is.null(gds)) {
    stop(paste("Error: No GEO data found for accession", opt$geo_accession))
  }
  
  # Extract and save metadata
  invisible(purrr::map(seq_along(gds), function(i) {
    # Extract metadata for each ExpressionSet
    if (inherits(gds[[i]], "ExpressionSet")) {
      metadata <- pData(gds[[i]])
      
      # Extract name of the matrix and format the output file name
      matrix_name <- names(gds)[i]
      base_name <- sub("\\..*$", "", matrix_name)  # Remove everything after the first dot
      
      # Define the output file path
      output_file <- file.path(opt$output_dir, paste0(base_name, "_metadata.csv"))
      
      # Write metadata to file using readr's write_csv
      log_message(paste("Saving metadata to:", output_file), NULL)
      write_csv(metadata, file = output_file)
    } else {
      log_message(paste("Skipping unsupported object type at index", i), NULL)
    }
  }))
  
  log_message("Metadata downloaded and saved successfully!", NULL)
}

# Run the main function
main()