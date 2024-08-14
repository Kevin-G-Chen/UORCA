# Load required libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(GEOquery))

# Define command-line options
option_list <- list(
  make_option(c("-g", "--geo_accession"), type = "character", default = NULL,
              help = "GEO accession number", metavar = "character")
)

# Parse the options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if GEO accession is provided
if (is.null(opt$geo_accession)) {
  print_help(opt_parser)
  stop("Please specify a GEO accession number using the -g or --geo_accession option.", call. = FALSE)
}

# Get the GEO accession from the options
geo_accession <- opt$geo_accession

# Fetch the GEO data
gse <- getGEO(geo_accession, GSEMatrix = TRUE, getGPL = FALSE)

# Print the structure of the loaded object
print("Structure of the GEO object:\n")
str(gse)

# Check for expression data
if (length(gse) > 0 && !is.null(exprs(gse[[1]]))) {
  print("\nExpression data found. Dimensions:\n")
  print(dim(exprs(gse[[1]])))
  
  # Print the first few rows and columns of the expression data
  print("\nFirst few rows and columns of expression data:\n")
  print(exprs(gse[[1]])[1:5, 1:5])
} else {
  print("\nNo expression data found in this GEO accession.\n")
}

# Return the GEO object for further analysis if needed
return(gse)