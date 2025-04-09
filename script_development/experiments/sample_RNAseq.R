user_lib <- Sys.getenv("R_LIBS_USER", unset="~/R/library")
if (!dir.exists(user_lib)) {
    dir.create(user_lib, recursive = TRUE)
}
.libPaths(c(user_lib, .libPaths()))

library(pacman)
# Load only the required packages (removing tidyverse)
p_load(edgeR, limma, tximport, gplots)

# Create plot output directories and subdirectories
create_plot_output_dirs <- function(output_dir) {
  plots_dir <- file.path(output_dir, "plots")
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
  # Create subdirectories for different plot types
  for(subdir in c("mds", "filtering", "normalization", "voom", "ma", "volcano", "heatmap")) {
    dir.create(file.path(plots_dir, subdir), showWarnings = FALSE)
  }
  return(plots_dir)
}

# 1. MDS Plot
plot_mds <- function(dge, output_dir, color_by = NULL, shape_by = NULL) {
  plots_dir <- file.path(output_dir, "plots", "mds")
  pdf(file.path(plots_dir, "mds_plot_group.pdf"), width = 8, height = 6)
  if (is.null(color_by)) {
    if ("group" %in% colnames(dge$samples)) {
      color_by <- dge$samples$group
      colors <- rainbow(length(unique(color_by)))
      plotMDS(dge, col = colors[as.factor(color_by)], main = "MDS Plot - Colored by Group")
      legend("topright", legend = levels(as.factor(color_by)), col = colors, pch = 16)
    } else {
      plotMDS(dge, main = "MDS Plot")
    }
  } else {
    colors <- rainbow(length(unique(color_by)))
    plotMDS(dge, col = colors[as.factor(color_by)], main = "MDS Plot - Colored by Group")
    legend("topright", legend = levels(as.factor(color_by)), col = colors, pch = 16)
  }
  dev.off()
  if (!is.null(shape_by) || "lane" %in% colnames(dge$samples)) {
    if (is.null(shape_by) && "lane" %in% colnames(dge$samples)) { shape_by <- dge$samples$lane }
    pdf(file.path(plots_dir, "mds_plot_batch.pdf"), width = 8, height = 6)
    colors <- rainbow(length(unique(shape_by)))
    plotMDS(dge, col = colors[as.factor(shape_by)], dim = c(3,4),
            main = "MDS Plot (Dims 3,4) - Colored by Batch")
    legend("topright", legend = levels(as.factor(shape_by)), col = colors, pch = 16)
    dev.off()
  }
  return(paste0("MDS plots saved to ", plots_dir))
}

# 2. Filtering Plots (Density plots before/after filtering)
plot_filtering <- function(lcpm_pre, lcpm_post, output_dir) {
  plots_dir <- file.path(output_dir, "plots", "filtering")
  pdf(file.path(plots_dir, "filtering_density.pdf"), width = 10, height = 6)
  par(mfrow = c(1, 2))
  nsamples <- ncol(lcpm_pre)
  col <- rainbow(nsamples)
  plot(density(lcpm_pre[,1]), col = col[1], lwd = 2, ylim = c(0, 0.26),
       las = 2, main = "Before Filtering", xlab = "Log-CPM")
  for(i in 2:nsamples) { lines(density(lcpm_pre[,i]), col = col[i], lwd = 2) }
  plot(density(lcpm_post[,1]), col = col[1], lwd = 2, ylim = c(0, 0.26),
       las = 2, main = "After Filtering", xlab = "Log-CPM")
  for(i in 2:nsamples) { lines(density(lcpm_post[,i]), col = col[i], lwd = 2) }
  if (nsamples <= 12) {
    legend("topright", colnames(lcpm_post), col = col, lwd = 2, cex = 0.6)
  }
  dev.off()
  return(paste0("Filtering plots saved to ", plots_dir))
}

# 3. Normalization Boxplots
plot_normalization <- function(dge_raw, dge_norm, output_dir) {
  plots_dir <- file.path(output_dir, "plots", "normalization")
  lcpm_raw <- cpm(dge_raw, log = TRUE)
  lcpm_norm <- cpm(dge_norm, log = TRUE)
  pdf(file.path(plots_dir, "normalization_boxplots.pdf"), width = 10, height = 6)
  par(mfrow = c(1, 2))
  ylim <- range(c(lcpm_raw, lcpm_norm))
  boxplot(lcpm_raw, las = 2, col = rainbow(ncol(lcpm_raw)), 
          main = "Before Normalization", ylab = "Log-CPM", ylim = ylim)
  boxplot(lcpm_norm, las = 2, col = rainbow(ncol(lcpm_norm)), 
          main = "After Normalization", ylab = "Log-CPM", ylim = ylim)
  dev.off()
  return(paste0("Normalization plots saved to ", plots_dir))
}

# 4. Voom Plot
save_voom_plot <- function(v_object, output_dir) {
  plots_dir <- file.path(output_dir, "plots", "voom")
  pdf(file.path(plots_dir, "voom_mean_variance.pdf"), width = 8, height = 6)
  smoothScatter(v_object$voom.xy, xlab = "log2( count size + 0.5 )", 
                ylab = "Sqrt( standard deviation )", main = "Voom: Mean-Variance Trend")
  if (!is.null(v_object$voom.line)) {
    ord <- order(v_object$voom.xy[,1])
    lines(v_object$voom.xy[ord,1], v_object$voom.line[ord], col = "red", lwd = 2)
  }
  dev.off()
  pdf(file.path(plots_dir, "limma_mean_variance_fit.pdf"), width = 8, height = 6)
  plotSA(v_object, main = "Final Mean-Variance Trend")
  dev.off()
  return(paste0("Voom plots saved to ", plots_dir))
}

# 5. MA Plot
plot_ma <- function(fit_object, output_dir, coef = 1, status = NULL, 
                    top_genes = 10, highlight_color = "red") {
  plots_dir <- file.path(output_dir, "plots", "ma")
  contrast_name <- colnames(fit_object)[coef]
  if (is.null(status)) { status <- decideTests(fit_object)[, coef] }
  pdf(file.path(plots_dir, paste0("ma_plot_", gsub("[^a-zA-Z0-9]", "_", contrast_name), ".pdf")), width = 8, height = 6)
  plotMA(fit_object, coef = coef, status = status, 
         main = paste0("MA Plot: ", contrast_name), 
         xlab = "Average log-expression", ylab = "Log-fold change")
  if (!is.null(top_genes) && top_genes > 0) {
    top_indices <- order(fit_object$p.value[, coef])[1:min(top_genes, nrow(fit_object))]
    if (length(top_indices) > 0) {
      with(fit_object, points(Amean[top_indices], coefficients[top_indices, coef], 
                             col = highlight_color, pch = 16))
      with(fit_object, text(Amean[top_indices], coefficients[top_indices, coef], 
                          labels = rownames(fit_object)[top_indices], pos = 2, cex = 0.7))
    }
  }
  dev.off()
  return(paste0("MA plot for ", contrast_name, " saved to ", plots_dir))
}

# 6. Volcano Plot
plot_volcano <- function(fit_object, output_dir, coef = 1, status = NULL, top_genes = 10,
                         lfc_threshold = 1, pval_threshold = 0.05) {
  plots_dir <- file.path(output_dir, "plots", "volcano")
  contrast_name <- colnames(fit_object)[coef]
  if (is.null(status)) { status <- decideTests(fit_object)[, coef] }
  if (!is.null(fit_object$p.value)) { pvals <- fit_object$p.value[, coef]
  } else if (!is.null(fit_object$table)) { pvals <- fit_object$table[, "PValue"]
  } else { stop("Cannot find p-values in the fit object") }
  if (!is.null(fit_object$coefficients)) { lfc <- fit_object$coefficients[, coef]
  } else if (!is.null(fit_object$table)) { lfc <- fit_object$table[, "logFC"]
  } else { stop("Cannot find log-fold changes in the fit object") }
  df <- data.frame(logFC = lfc, negLogP = -log10(pvals), status = status)
  pdf(file.path(plots_dir, paste0("volcano_plot_", gsub("[^a-zA-Z0-9]", "_", contrast_name), ".pdf")), width = 8, height = 6)
  status_colors <- c("blue", "black", "red")
  point_colors <- status_colors[status + 2]
  plot(df$logFC, df$negLogP, xlab = "Log2 Fold Change", ylab = "-Log10 P-value",
       main = paste0("Volcano Plot: ", contrast_name), pch = 20, cex = 0.5, col = point_colors)
  abline(v = c(-lfc_threshold, lfc_threshold), lty = 2, col = "gray")
  abline(h = -log10(pval_threshold), lty = 2, col = "gray")
  if (!is.null(top_genes) && top_genes > 0) {
    top_indices <- order(pvals)[1:min(top_genes, length(pvals))]
    if (length(top_indices) > 0) {
      points(df$logFC[top_indices], df$negLogP[top_indices], col = "purple", pch = 16)
      text(df$logFC[top_indices], df$negLogP[top_indices], labels = rownames(df)[top_indices],
          pos = 2, cex = 0.7)
    }
  }
  legend("topleft", legend = c("Down", "Not significant", "Up"), col = status_colors, pch = 20)
  dev.off()
  return(paste0("Volcano plot for ", contrast_name, " saved to ", plots_dir))
}

# 7. Heatmap of Top DEGs
plot_deg_heatmap <- function(dge, fit_object, output_dir, coef = 1, n_genes = 100, clustering_method = "complete") {
  plots_dir <- file.path(output_dir, "plots", "heatmap")
  contrast_name <- colnames(fit_object)[coef]
  if (!is.null(fit_object$p.value) && !is.null(fit_object$lods)) {
    o <- order(fit_object$p.value[, coef])
  } else if (!is.null(fit_object$table)) {
    o <- order(fit_object$table$PValue)
  } else {
    stop("Cannot determine top genes from the fit object")
  }
  top_genes_idx <- head(o, n_genes)
  lcpm <- cpm(dge, log = TRUE)
  lcpm_top <- lcpm[top_genes_idx, ]
  if (!is.null(fit_object$genes)) {
    gene_symbols <- fit_object$genes$SYMBOL[top_genes_idx]
    if (is.null(gene_symbols)) { gene_symbols <- rownames(lcpm_top) }
  } else if (!is.null(dge$genes)) {
    gene_symbols <- dge$genes$SYMBOL[top_genes_idx]
    if (is.null(gene_symbols)) { gene_symbols <- rownames(lcpm_top) }
  } else {
    gene_symbols <- rownames(lcpm_top)
  }
  if (!is.null(dge$samples$group)) { groups <- dge$samples$group } else { groups <- colnames(lcpm_top) }
  lcpm_scaled <- t(scale(t(lcpm_top)))
  col_palette <- colorRampPalette(c("blue", "white", "red"))(100)
  pdf(file.path(plots_dir, paste0("heatmap_top", n_genes, "_", gsub("[^a-zA-Z0-9]", "_", contrast_name), ".pdf")),
      width = 10, height = max(8, n_genes/10))
  library(gplots)
  heatmap.2(lcpm_scaled, col = col_palette, labRow = gene_symbols, labCol = groups,
            trace = "none", density.info = "none", margin = c(8, 10), cexRow = 0.6,
            cexCol = 1, key.title = NA, key.xlab = "Z-score",
            main = paste0("Top ", n_genes, " DE Genes: ", contrast_name),
            hclustfun = function(x) hclust(x, method = clustering_method),
            scale = "none")
  dev.off()
  return(paste0("Heatmap of top ", n_genes, " DEGs for ", contrast_name, " saved to ", plots_dir))
}

cat("=== R Script: edgeR/limma Analysis Start ===\n")

args <- commandArgs(trailingOnly = TRUE)

metadata_file <- args[1]
merged_group <- args[2]
output_dir <- args[3]
tx2gene_file <- args[4]

plots_dir <- create_plot_output_dirs(output_dir)

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
lcpm_pre <- cpm(DGE, log = TRUE)

# Combine DGE$samples with metadata (bind_cols replaced with cbind)
DGE$samples <- cbind(DGE$samples, metadata)

# Assuming filtering is done here
# Example: keep.exprs <- filterByExpr(DGE)
# lcpm_post <- cpm(DGE[keep.exprs, ], log = TRUE)
# plot_filtering(lcpm_pre, lcpm_post, output_dir)

cat("Step 5: Normalizing the DGEList using calcNormFactors...\n")
cat("Normalizing DGEList...\n")
DGE.norm <- calcNormFactors(DGE)

dge_raw <- DGE  # or the appropriate variable holding raw counts
plot_normalization(dge_raw, DGE.norm, output_dir)

cat("Step 6: Generating design matrix using grouping column:", merged_group, "\n")
cat("Creating design matrix using grouping column:", merged_group, "\n")
design <- model.matrix(as.formula(paste("~0 +", merged_group)), data = DGE.norm$samples)
colnames(design) <- sub(merged_group, "", colnames(design))
cat("Design matrix:\n")
print(design)

plot_mds(DGE.norm, output_dir, color_by = DGE.norm$samples$group, shape_by = DGE.norm$samples$lane)

cat("Step 7: Performing voom transformation...\n")
cat("Performing voom transformation...\n")
v <- voom(DGE.norm, design, plot = FALSE)

save_voom_plot(v, output_dir)

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
  dt <- decideTests(fit2, p.value = 0.05)
  plot_ma(fit2, output_dir, coef = 1, status = dt)
  plot_volcano(fit2, output_dir, coef = 1, status = dt)
  plot_deg_heatmap(DGE.norm, fit2, output_dir, coef = 1, n_genes = 100)

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
    dt <- decideTests(fit, p.value = 0.05)[, i]
    plot_ma(fit, output_dir, coef = i, status = dt)
    plot_volcano(fit, output_dir, coef = i, status = dt)
    plot_deg_heatmap(DGE.norm, fit, output_dir, coef = i, n_genes = 100)
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
