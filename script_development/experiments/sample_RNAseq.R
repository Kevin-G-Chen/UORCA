user_lib <- Sys.getenv("R_LIBS_USER", unset="~/R/library")
if (!dir.exists(user_lib)) {
    dir.create(user_lib, recursive = TRUE)
}
.libPaths(c(user_lib, .libPaths()))

library(pacman)
# Load only the required packages (removing tidyverse)
p_load(edgeR, limma, tximport, gplots)

cat("=== R Script: edgeR/limma Analysis Start ===\n")

args <- commandArgs(trailingOnly = TRUE)

metadata_file <- args[1]
merged_group <- args[2]
output_dir <- args[3]
tx2gene_file <- args[4]

cat("Received command arguments:\n")
cat("  metadata_file =", metadata_file, "\n")
cat("  merged_group =", merged_group, "\n")
cat("  output_dir =", output_dir, "\n")
cat("  tx2gene_file =", tx2gene_file, "\n\n")

# Create DEG subfolder within the output directory
deg_dir <- file.path(output_dir, "DEG")
dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)
cat("Created DEG results directory at:", deg_dir, "\n")

cat("Step 1: Loading metadata from file...\n")
cat("Loading metadata from:", metadata_file, "\n")
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
if(nrow(metadata) == 0) {
  stop("Metadata is empty!")
}

if(!"abundance_file" %in% colnames(metadata)) {
  stop("Metadata must contain an 'abundance_file' column!")
}

# Load tx2gene mapping if provided
if(tx2gene_file != "NA" && file.exists(tx2gene_file)){
  cat("Step 2: Checking for and loading tx2gene mapping if provided...\n")
  cat("Loading tx2gene mapping from:", tx2gene_file, "\n")
  tx2gene <- read.delim(tx2gene_file, header = FALSE, stringsAsFactors = FALSE)
  # Select columns 1 and 3 (like dplyr::select(1,3))
  tx2gene <- tx2gene[, c(1, 3)]
  # Rename columns (like setNames)
  names(tx2gene) <- c("TXNAME", "GENEID")
  # Filter rows where GENEID is not NA and not empty (like dplyr::filter)
  tx2gene <- tx2gene[!is.na(tx2gene$GENEID) & tx2gene$GENEID != "", ]
  str(tx2gene)
  use_tx2gene <- TRUE
} else {
  cat("No valid tx2gene file provided. Proceeding without tx2gene mapping.\n")
  use_tx2gene <- FALSE
}

cat("Step 3: Importing Kallisto quantification data using tximport...\n")
cat("Importing Kallisto quantification data...\n")
if(use_tx2gene){
  kallisto <- tximport(metadata$abundance_file, type = "kallisto",
                       tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM",
                       ignoreAfterBar = TRUE)
} else {
  kallisto <- tximport(metadata$abundance_file, type = "kallisto",
                       txOut = TRUE, ignoreAfterBar = TRUE)
}

cat("Step 4: Creating DGEList with the imported counts...\n")
cat("Creating DGEList...\n")
DGE <- DGEList(counts = kallisto$counts)
# Combine DGE$samples with metadata (bind_cols replaced with cbind)
DGE$samples <- cbind(DGE$samples, metadata)

cat("Step 5: Normalizing the DGEList using calcNormFactors...\n")
cat("Normalizing DGEList...\n")
DGE.norm <- calcNormFactors(DGE)

cat("Step 6: Generating design matrix using grouping column:", merged_group, "\n")
cat("Creating design matrix using grouping column:", merged_group, "\n")
design <- model.matrix(as.formula(paste("~0 +", merged_group)), data = DGE.norm$samples)
colnames(design) <- sub(merged_group, "", colnames(design))
cat("Design matrix:\n")
print(design)

cat("Step 7: Performing voom transformation...\n")
cat("Performing voom transformation...\n")
v <- voom(DGE.norm, design, plot = FALSE)

cat("Step 8: Fitting linear model with limma...\n")
cat("Fitting linear model...\n")
fit <- lmFit(v, design)
fit <- eBayes(fit)


if(ncol(design) == 2){
  cat("Exactly two groups detected. Calculating contrast (group2 - group1)...\n")
  contrast_name <- paste(colnames(design)[2], "-", colnames(design)[1])
  contrast <- makeContrasts(diff = contrast_name, levels = design)
  cat("Contrast matrix:")
  print(contrast)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  cat("Top differential expression results for contrast:\n")
  top_results <- topTable(fit2, number = Inf)
  print(head(top_results))
  top_results$Gene <- rownames(top_results)
  top_results <- top_results[, c("Gene", setdiff(names(top_results), "Gene"))]
  print(head(top_results))
  write.csv(top_results, file = file.path(deg_dir, "DEG_results.csv"), row.names = FALSE)

} else {
  cat("Multiple groups detected. Generating top results for each coefficient...\n")
  for(i in 1:ncol(design)){
    coef_name <- colnames(design)[i]
    cat("Top results for", coef_name, ":\n")
    top_results <- topTable(fit, coef = i, number = Inf)

    top_results$Gene <- rownames(top_results)
    top_results <- top_results[, c("Gene", setdiff(names(top_results), "Gene"))]
    print(head(top_results))
    write.csv(top_results, file = file.path(deg_dir, paste0("DEG_results_", coef_name, ".csv")), row.names = FALSE)
  }
}

cat("Step 9: Saving DEG results and normalized DGEList...\n")
cat("Saving normalized DGEList object...\n")
saveRDS(DGE.norm, file = file.path(output_dir, "DGE_norm.RDS"))

cat("=== R Script: edgeR/limma Analysis Completed ===\n")
