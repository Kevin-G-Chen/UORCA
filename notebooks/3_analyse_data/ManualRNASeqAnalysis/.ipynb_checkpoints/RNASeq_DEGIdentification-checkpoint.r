suppressMessages({
  suppressWarnings({
    library(tidyverse)
    library(edgeR)
    library(limma)
    library(tximport)
    library(DESeq2)
    library(janitor)
    library(viridis)
    library(optparse)
    library(jsonlite)  # For handling JSON strings
  })
})

# Define options for optparse
option_list <- list(
  make_option(c("--DGE"), type="character", default=NULL, help="Path to the DGE object (RDS file)", metavar="character"),
  make_option(c("--column"), type="character", default=NULL, help="Column to use in model.matrix", metavar="character"),
  make_option(c("--comparisons"), type="character", default=NULL, help="JSON string of comparisons", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if required options are provided
if (is.null(opt$DGE) | is.null(opt$column) | is.null(opt$comparisons)) {
  print_help(opt_parser)
  stop("Please provide all the required arguments (--DGE, --column, --comparisons)", call.=FALSE)
}

### Step 1 - Read in the DGEList object
DGE <- readRDS(opt$DGE)

### Step 2 - Clean and Parse Comparisons JSON

# Remove unnecessary characters and convert the string to valid JSON format
clean_comparisons <- opt$comparisons

# Remove the outer brackets and single quotes, and ensure it's a valid JSON array
clean_comparisons <- gsub("^\\[|\\]$", "", clean_comparisons)  # Remove the outer brackets
clean_comparisons <- gsub("^'|'$", "", clean_comparisons)  # Remove any leading and trailing single quotes
clean_comparisons <- gsub("', '", ",", clean_comparisons)  # Replace the separator with a comma
clean_comparisons <- gsub("\\\\", "", clean_comparisons)  # Remove escape backslashes

# Wrap in square brackets to form a valid JSON array
clean_comparisons <- paste0("[", clean_comparisons, "]")

# Parse the cleaned JSON string
comparisons <- fromJSON(clean_comparisons)

### Step 3 - Construct contrast matrix

# Use the column provided in the options
design <- model.matrix(data = DGE$samples, ~0 + get(opt$column))
colnames(design) <- str_remove_all(colnames(design), opt$column)

# Create a contrast matrix using the parsed comparisons
contrast.matrix <- makeContrasts(
  contrasts = setNames(lapply(comparisons, function(cmp) cmp$comparison), 
                       lapply(comparisons, function(cmp) cmp$name)),
  levels = colnames(design)
)

### Step 4 - Perform the DEG analysis

v <- voom(DGE, design)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrast.matrix)
efit <- eBayes(vfit)

### Step 5 - Record results of DEG analysis

contrasts <- colnames(contrast.matrix)

LFC.summary <- lapply(contrasts, function(x){
    topTable(efit, coef = x, number = Inf) %>% arrange(adj.P.Val)
})

# Example output: print the top results for the first contrast
print(LFC.summary[[1]])