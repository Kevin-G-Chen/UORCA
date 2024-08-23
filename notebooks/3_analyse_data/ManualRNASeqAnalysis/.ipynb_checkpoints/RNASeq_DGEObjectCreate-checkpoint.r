# Load libraries silently
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
  })
})

# Function to print messages with a timestamp
print_message <- function(message) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), message))
}

# Set up command-line options
option_list <- list(
  make_option(c("-d", "--directory"), type="character", default=NULL,
              help="Directory containing Kallisto abundance files", metavar="character"),
  make_option(c("-t", "--t2g"), type="character", default=NULL,
              help="Path to t2g.txt file", metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default=NULL,
              help="Path to metadata CSV file", metavar="character"),
  make_option(c("-c", "--clean_columns"), type="character", default=NULL,
              help="Columns to clean by removing spaces, separated by commas", metavar="character"),
  make_option(c("-g", "--group"), type="character", default="genotype_clean",
              help="Group column to use for filtering [default: %default]", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="./output/",
              help="Directory to save output files [default: %default]", metavar="character"),
  make_option(c("-r", "--geo_sra_mapping"), type="character", default=NULL,
              help="Path to GEO-SRA mapping file", metavar="character")
)

# Parse options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Ensure required arguments are provided
if (is.null(opt$directory) | is.null(opt$t2g) | is.null(opt$metadata) | is.null(opt$geo_sra_mapping)) {
  print_help(opt_parser)
  stop("Please supply the required arguments.", call.=FALSE)
}

# Ensure the output directory ends with a slash
if (!grepl("/$", opt$output)) {
  opt$output <- paste0(opt$output, "/")
}

# Create output directory if it doesn't exist
if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
}

# Step 1: Load GEO-SRA mapping data
print_message("Step 1: Loading GEO-SRA mapping data...")
geo_sra_map <- read_tsv(opt$geo_sra_mapping) %>%
  dplyr::rename(geo_accession = sample_ID)
print_message("GEO-SRA mapping data loaded.")

# Step 2: Load metadata
print_message("Step 2: Loading metadata...")
meta <- read_csv(opt$metadata)
print_message("Original metadata column names:")
print(colnames(meta))

meta <- meta %>%
  clean_names()
print_message("Cleaned metadata column names:")
print(colnames(meta))

# Step 3: Clean specified columns
if (!is.null(opt$clean_columns)) {
  print_message("Step 3: Cleaning specified columns...")
  cols_to_clean <- str_split(opt$clean_columns, ",")[[1]]
  for (col in cols_to_clean) {
    meta <- meta %>%
      mutate(!!sym(col) := str_remove_all(!!sym(col), " "))
  }
  print_message("Columns cleaned:")
  print(cols_to_clean)
}

# Step 4: Join metadata with GEO-SRA mapping data
print_message("Step 4: Joining metadata with GEO-SRA mapping data...")
meta <- left_join(meta, geo_sra_map, by = "geo_accession")
print_message("Metadata joined with GEO-SRA mapping data.")

# Step 5: List the abundance files
print_message("Step 5: Listing abundance files...")
files <- list.files(path = opt$directory,
                    pattern = "abundance.tsv",
                    recursive = TRUE,
                    full.names = TRUE)
print_message(sprintf("Abundance files found: %d files", length(files)))

file_data <- data.frame(path = files,
                        SRA_ID = basename(dirname(files)))

meta <- left_join(meta, file_data, by = "SRA_ID")
print_message("Metadata updated with abundance file paths.")

# Step 6: Load tx2gene data
print_message("Step 6: Loading tx2gene data...")
tx2gene <- read_tsv(opt$t2g, col_names = FALSE) %>%
  dplyr::select(1, 3) %>%
  drop_na()
print_message("tx2gene data loaded.")

# Step 7: Import the quantification files
print_message("Step 7: Importing quantification files...")
kallisto <- tximport(files = files,
                     type = "kallisto",
                     tx2gene = tx2gene,
                     ignoreAfterBar = TRUE,
                     countsFromAbundance = "lengthScaledTPM")
print_message("Quantification files imported.")

# Step 8: Create a DGEList object
print_message("Step 8: Creating DGEList object...")
DGE <- DGEList(counts = kallisto$counts, samples = meta)
print_message("DGEList object created.")

# Step 9: Filter based on expression
print_message("Step 9: Filtering based on expression...")
keep.exprs <- filterByExpr(DGE, group = DGE$samples[[opt$group]])
DGE.filtered <- DGE[keep.exprs, keep.lib.sizes = FALSE]
print_message(sprintf("Filtering complete. Number of genes retained: %d", sum(keep.exprs)))

# Step 10: Plot filtering results
print_message("Step 10: Plotting filtering results...")
L <- mean(DGE$samples$lib.size) * 1e-6
M <- median(DGE$samples$lib.size) * 1e-6

cpm <- cpm(DGE)
lcpm <- cpm(DGE, log = TRUE)

lcpm.cutoff <- log2(10/M + 2/L)
nsamples <- ncol(DGE)

col <- viridis(n = nsamples, direction = -1)

png(file.path(opt$output, "Filtering.png"))

par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.5), las=2, main="", xlab="")
title(main="A. Unfiltered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples) {
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
lcpm <- cpm(DGE.filtered, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.5), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples) {
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}

dev.off()
print_message("Filtering plot saved as Filtering.png.")

# Step 11: Normalize the DGEList object
print_message("Step 11: Normalizing the DGEList object...")
DGE.final <- calcNormFactors(DGE.filtered)
print_message("Normalization complete.")

# Step 12: Plot normalization results
print_message("Step 12: Plotting normalization results...")
png(file.path(opt$output, "Normalization.png"))

par(mfrow=c(1,2))
lcpm <- cpm(DGE.filtered, log=TRUE)
boxplot(lcpm, las=2, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
lcpm <- cpm(DGE.final, log=TRUE)
boxplot(lcpm, las=2, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")

dev.off()
print_message("Normalization plot saved as Normalization.png.")

# Step 13: Save the final DGEList object
print_message("Step 13: Saving the final DGEList object...")
saveRDS(DGE.final, file.path(opt$output, "DGE.RDS"))
print_message("Final DGEList object saved as DGE.RDS.")

print_message("Script execution completed successfully.")