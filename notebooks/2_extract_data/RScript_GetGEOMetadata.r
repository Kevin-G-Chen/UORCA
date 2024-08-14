# Load required libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(GEOquery))
suppressPackageStartupMessages(library(tidyverse))

# Define command-line options
option_list <- list(
  make_option(c("-g", "--geo_accession"), type = "character", default = NULL,
              help = "GEO accession number (e.g., GSE12345)", metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character", default = ".",
              help = "Directory where the metadata file will be saved", metavar = "character")
)
# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if GEO accession is provided
if (is.null(opt$geo_accession)) {
  print_help(opt_parser)
  stop("Error: GEO accession must be provided", call. = FALSE)
}

if (!dir.exists(opt$output_dir)) {
  cat("Output directory does not exist. Creating directory:", opt$output_dir, "\n")
  dir.create(opt$output_dir, recursive = TRUE)
} else {
  cat("Output directory exists:", opt$output_dir, "\n")
}

# Download and save GEO metadata
cat("Retrieving metadata for GEO accession:", opt$geo_accession, "\n")
gds <- getGEO(opt$geo_accession, GSEMatrix = TRUE, AnnotGPL = TRUE)

# Extract metadata
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
    cat("Saving metadata to:", output_file, "\n")
    write_csv(metadata, file = output_file)
  } else {
    cat("Skipping unsupported object type at index", i, "\n")
  }
}))

cat("Metadata saved successfully!\n")