# Load libraries silently
suppressMessages({
  suppressWarnings({
    library(edgeR)
    library(optparse)
  })
})

# Define command-line options
option_list <- list(
  make_option(c("-p", "--path"), type="character", default=NULL, 
              help="Path to the DGE RDS file", metavar="character")
)

# Parse command-line options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if the path is provided
if (is.null(opt$path)) {
  print_help(opt_parser)
  stop("Path to the DGE RDS file must be provided.", call.=FALSE)
}

# Step 1: Read in the DGEList object
DGE <- readRDS(opt$path)

# Step 2: Extract metadata
meta <- DGE$samples

# Step 3: Print column names and unique values
cat("Metadata columns and their unique values:\n")

for (col_name in colnames(meta)) {
  unique_values <- unique(meta[[col_name]])
  
  # Ignore the column if it has only one unique value
  if (length(unique_values) == 1) {
    next
  }
  
  # Ignore the column if all values are unique
  if (length(unique_values) == nrow(meta)) {
    next
  }
  
  # Check if any unique value exceeds 100 characters
  if (any(nchar(as.character(unique_values)) > 100)) {
    next  # Skip this column
  }
  
  # Print column name and unique values
  cat(sprintf("Column: %s\nValues: %s\n\n", col_name, paste(unique_values, collapse=", ")))
}