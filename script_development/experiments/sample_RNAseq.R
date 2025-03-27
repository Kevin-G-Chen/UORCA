# This notebook will contain code pertaining to the original processing of the RNAseq analysis for BRAF samples.
# The overall design involves three separate analyses - CFC1, SNV2, and NS4. I actually don't know which is the "main" one of interest...

# %% Loading libraries
library(pacman)
p_load(
    tidyverse,
    edgeR,
    limma,
    DESeq2,
    org.Hs.eg.db,
    enrichplot,
    clusterProfiler,
    matrixStats,
    ComplexUpset,
    patchwork,
    DOSE,
    tximport,
    viridis
)
getwd()
# %% Load files
metadata <- read_csv("../input/2024_12_18_InitialInputs/metadata/RNASeq_Metadata_BRAF.csv")
kallisto_paths <- list.files(
    path = "../input/2024_12_18_InitialInputs/Kallisto/CFC1",
    recursive = TRUE,
    full.names = TRUE,
    pattern = "abundance.h5"
)
kallisto_filenames <- basename(dirname(kallisto_paths))
kallisto_df <- data.frame(
    path = kallisto_paths,
    file_name = kallisto_filenames
)
# Light processing on the metadata
metadata <- metadata %>%
    left_join(kallisto_df,
        by = "file_name"
    ) %>%
    mutate(
        geno_short = if_else(
            genotype == "WT/WT",
            "WT",
            "CFC1"
        ),
        merged_group = paste0(geno_short, "_", cell_type)
    ) %>%
    filter(!is.na(path))
str(metadata)

# %% Prepare the DGEList object
tx2gene <- read_csv(file = "../input/2024_12_18_InitialInputs/tx2gene_entrez_v38.csv", col_names = FALSE)
kallisto <- tximport(metadata$path,
    type = "kallisto",
    tx2gene = tx2gene,
    countsFromAbundance = "lengthScaledTPM",
    ignoreAfterBar = T
)
# %% Prepare DGEList object part 2
DGE <- DGEList(kallisto$counts)
DGE$samples <- bind_cols(DGE$samples, metadata)
DGE$genes <- bitr(
    geneID = rownames(DGE),
    fromType = "ENTREZID",
    toType = "SYMBOL",
    OrgDb = "org.Hs.eg.db",
    drop = FALSE
)
colnames(DGE$genes) <- c("GeneID", "Symbol")
identical(rownames(DGE), DGE$genes$GeneID) # Sanity check that the genes match

# Because they match we can save the output and be happy
saveRDS(DGE, "../scratch/R_outputs/InitialProcessing/CFC1/RawDGE.RDS")
# %% Now perform some processing

# Filtering
keep.exprs <- filterByExpr(DGE,
    group = DGE$samples$merged_group
)
DGE.filtered <- DGE[keep.exprs, keep.lib.sizes = FALSE]
dim(DGE)
dim(DGE.filtered)
# 16248 genes remain. Plot the results

dir.create("../results/2024_12_20_InitialProcessing/CFC1/QC",
    recursive = TRUE
)
L <- mean(DGE$samples$lib.size) * 1e-6
M <- median(DGE$samples$lib.size) * 1e-6

cpm <- cpm(DGE)
lcpm <- cpm(DGE, log = T)

lcpm.cutoff <- log2(10 / M + 2 / L)
nsamples <- ncol(DGE) # define number of samples

col <- inferno(n = nsamples, direction = -1)
pdf("../results/2024_12_20_InitialProcessing/CFC1/QC/Filtering.pdf",
    height = 7,
    width = 9
)

par(mfrow = c(1, 2)) # for seeing plots side by side

lcpm <- cpm(DGE, log = T)
plot(density(lcpm[, 1]), col = col[1], lwd = 2, ylim = c(0, 0.2), las = 2, main = "", xlab = "")
title(main = "A. Unfiltered data", xlab = "Log-cpm") # creating titles
abline(v = lcpm.cutoff, lty = 3) # add the vertical line to represent the lcpm cutoff
for (i in 2:nsamples) { # for each sample...
    den <- density(lcpm[, i]) # compute the Kernel density estimate?
    lines(den$x, den$y, col = col[i], lwd = 2)
}
lcpm <- cpm(DGE.filtered, log = TRUE)
plot(density(lcpm[, 1]), col = col[1], lwd = 2, ylim = c(0, 0.2), las = 2, main = "", xlab = "")
title(main = "B. Filtered data", xlab = "Log-cpm")
abline(v = lcpm.cutoff, lty = 3)
for (i in 2:nsamples) {
    den <- density(lcpm[, i])
    lines(den$x, den$y, col = col[i], lwd = 2)
}

dev.off()

# Normalisation

DGE.final <- calcNormFactors(DGE.filtered)
pdf("../results/2024_12_20_InitialProcessing/CFC1/QC/Normalisation.pdf",
    height = 7,
    width = 9
)
par(mfrow = c(1, 2))
lcpm <- cpm(DGE.filtered, log = TRUE)
boxplot(lcpm, las = 2, main = "")
title(main = "A. Unnormalised data", ylab = "Log-cpm")
lcpm <- cpm(DGE.final, log = TRUE)
boxplot(lcpm, las = 2, main = "")
title(main = "B. Normalised data", ylab = "Log-cpm")

dev.off()

saveRDS(
    DGE.final,
    "../scratch/R_outputs/InitialProcessing/CFC1/DGE_processed.RDS"
)



# %% Generate PCA plot
DGE.final <- readRDS("../scratch/R_outputs/InitialProcessing/CFC1/DGE_processed.RDS")

lcpm <- cpm(DGE.final, log = TRUE)

# Define colors using viridis palette for consistency
geno_colour <- DGE.final$samples$geno_short
colours <- c("red", "blue")
geno_colour <- colours[as.factor(geno_colour)]

genotypes <- c("CFC1", "WT")
pdf("../results/2024_12_20_InitialProcessing/CFC1/QC/PCA_plot.pdf",
    height = 7,
    width = 7
)

# Plot MDS
plotMDS(lcpm,
    col = geno_colour,
    labels = DGE.final$samples$sample_name,
    main = "PCA plot (CFC samples only)",
    cex = 0.7 # Reduce text size to 70% of default
)

# Add legend
legend("topright",
    legend = genotypes,
    col = colours,
    pch = 16,
    title = "Genotype"
)
dev.off()
# %% Produce CPM
cpm_towrite <- cpm(DGE.final) %>%
    as.data.frame()

colnames(cpm_towrite) <- DGE.final$samples$sample_name

Symbols <- bitr(rownames(cpm_towrite),
    fromType = "ENTREZID",
    toType = "SYMBOL",
    OrgDb = "org.Hs.eg.db",
    drop = FALSE
)
cpm_towrite$Symbol <- Symbols$SYMBOL

cpm_towrite <- cpm_towrite %>%
    rownames_to_column(var = "GeneID") %>%
    dplyr::select(Symbol, GeneID, everything())

write_csv(
    cpm_towrite,
    "../results/2024_12_20_InitialProcessing/CFC1/CPM.csv"
)

# %% Produce gene expression plots

# As a small sanity check, I will plot the expression of various iPSC and NPC cell markers.

SampleInfo <- DGE.final$samples %>%
    rownames_to_column(var = "SampleName") %>%
    dplyr::select(SampleName, sample_name, geno_short, cell_type, merged_group)

EntrezIDs <- c(
    "10763",
    "5080",
    "6656"
)
GeneData <- bitr(
    geneID = EntrezIDs,
    fromType = "ENTREZID",
    toType = "SYMBOL",
    OrgDb = "org.Hs.eg.db",
    drop = FALSE
)


cpm <- cpm(DGE.final) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ENTREZID") %>%
    filter(ENTREZID %in% EntrezIDs) %>%
    pivot_longer(
        cols = starts_with("Sample"),
        names_to = "SampleName",
        values_to = "CPM"
    )

data <- left_join(cpm, SampleInfo, by = "SampleName") %>%
    left_join(GeneData, by = "ENTREZID")
# %% Generation of expression plot
CreatePlot <- function(genename) {
    data_temp <- data %>%
        filter(ENTREZID == genename)

    Symbol <- GeneData %>%
        filter(ENTREZID == genename) %>%
        pull(SYMBOL)

    p <- ggplot(
        data = data_temp,
        aes(
            x = geno_short,
            y = CPM,
            fill = merged_group
        )
    ) +
        facet_grid(cols = vars(cell_type)) +
        geom_boxplot(outlier.shape = NULL) +
        scale_fill_manual(
            name = "Type of cell",
            values = c("pink", "#e41a1c", "lightblue", "#377eb8")
        ) +
        ggtitle(Symbol) +
        ylab("Counts per million (TMM normalised)") +
        xlab(NULL) +
        geom_point(aes(group = merged_group),
            position = position_dodge(width = 0.75)
        ) +
        theme(legend.position = "none")

    return(p)
}

plot_list <- lapply(EntrezIDs, FUN = CreatePlot)

design <- "
A
B
C
"

# probably could've just specified ncol = 1...

combined_plot <- wrap_plots(plot_list) +
    plot_layout(
        design = design,
        axis_titles = "collect"
    )

ggsave("../results/2024_12_20_InitialProcessing/CFC1/ExpressionPlots/NPC_Markers.pdf",
    plot = combined_plot,
    height = 14,
    width = 7
)
# %% DEG analysis
# I will now perform the DEG analysis

# First, construct the model matrix
design <- model.matrix(
    data = DGE.final$samples,
    ~ 0 + merged_group
)
colnames(design) <- str_remove_all(colnames(design), "merged_group")
str(design)
contrast.matrix <- makeContrasts(
    WT_differentiation = "WT_NPC - WT_IPSC",
    CFC1_differnetiation = "CFC1_NPC - CFC1_IPSC",
    DiffBetweenNPCs = "CFC1_NPC - WT_NPC",
    DiffofDiffns = "(CFC1_NPC - CFC1_IPSC) - (WT_NPC - WT_IPSC)",
    levels = colnames(design)
)
contrasts <- colnames(contrast.matrix)
v <- voom(DGE.final,
    design = design,
    plot = TRUE
)
vfit <- lmFit(
    v,
    design
)
vfit <- contrasts.fit(vfit,
    contrasts = contrast.matrix
)
efit <- eBayes(vfit)
plotSA(efit,
    main = "Mean-variance trend (using Empirical Bayes)"
)
summary(decideTests(efit))

LFC.summary <- sapply(contrasts, function(x) {
    lfc.list <- list()
    top <- topTable(efit,
        coef = x,
        number = Inf
    ) %>%
        list()

    lfc.list <- append(lfc.list, top)
})

saveRDS(
    LFC.summary,
    "../scratch/R_outputs/InitialProcessing/CFC1/LFC_summary.RDS"
)

# %% Save DEG CSVs

LFC.summary <- readRDS("../scratch/R_outputs/InitialProcessing/CFC1/LFC_summary.RDS")
contrasts <- names(LFC.summary)
sapply(contrasts, FUN = function(x) {
    a <- LFC.summary[[x]] %>%
        arrange(desc(logFC))

    write_csv(a, file = paste0("../results/2024_12_20_InitialProcessing/CFC1/DifferentialGeneExpression/", x, ".csv"))
})





# %% Enrichment analysis
contrasts <- colnames(contrast.matrix)

GSEA.lists <- sapply(contrasts, function(x) {
    gsea.list <- list()
    top <- topTable(efit,
        coef = x,
        number = Inf
    ) %>%
        arrange(desc(logFC))
    lfc <- top$logFC
    names(lfc) <- top$GeneID

    gsea.list <- append(gsea.list, list(lfc))
})
saveRDS(
    GSEA.lists,
    "./../scratch/R_outputs/InitialProcessing/CFC1/GSEALists.RDS"
)

# %% Performing enrichments
GSEA.Lists <- readRDS("../scratch/R_outputs/InitialProcessing/CFC1/GSEALists.RDS")

# %% Run DGN functions
str(GSEA.Lists)
DGN.unfiltered <- compareCluster(
    geneClusters = GSEA.Lists,
    fun = "gseDGN",
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    nPermSimple = 10000,
    seed = 42,
    eps = 0
) %>%
    setReadable("org.Hs.eg.db",
        keyType = "ENTREZID"
    )

DGN.p005 <- compareCluster(
    geneClusters = GSEA.Lists,
    fun = "gseDGN",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    nPermSimple = 10000,
    seed = 42,
    eps = 0
) %>%
    setReadable("org.Hs.eg.db",
        keyType = "ENTREZID"
    ) %>%
    pairwise_termsim(showCategory = Inf)

# Save these outputs
saveRDS(
    DGN.unfiltered,
    "../scratch/R_outputs/InitialProcessing/CFC1/DGN_unfiltered.RDS"
)
saveRDS(
    DGN.p005,
    "../scratch/R_outputs/InitialProcessing/CFC1/DGN_p005.RDS"
)

# %% Produce DGN CSVs

dir.create("../results/2024_12_20_InitialProcessing/CFC1/GSEA_Enrichments/DGN",
    recursive = TRUE
)

# Begin by creating CSVs

sapply(contrasts, FUN = function(x) {
    a <- DGN.unfiltered@compareClusterResult %>%
        filter(Cluster == x) %>%
        dplyr::select(!Cluster)
    b <- a %>%
        filter(p.adjust < 0.05)

    write_csv(a, file = paste0("../results/2024_12_20_InitialProcessing/CFC1/GSEA_Enrichments/DGN/DGN-", x, "_AllEnrichments.csv"))
    write_csv(b, file = paste0("../results/2024_12_20_InitialProcessing/CFC1/GSEA_Enrichments/DGN/DGN-", x, "_SignificantEnrichments.csv"))
})

# %% Run DO functions
str(GSEA.Lists)
DO.unfiltered <- compareCluster(
    geneClusters = GSEA.Lists,
    fun = "gseDO",
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    nPermSimple = 10000,
    seed = 42,
    eps = 0
) %>%
    setReadable("org.Hs.eg.db",
        keyType = "ENTREZID"
    )

DO.p005 <- compareCluster(
    geneClusters = GSEA.Lists,
    fun = "gseDO",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    nPermSimple = 10000,
    seed = 42,
    eps = 0
) %>%
    setReadable("org.Hs.eg.db",
        keyType = "ENTREZID"
    ) %>%
    pairwise_termsim(showCategory = Inf)

# Save these outputs
saveRDS(
    DO.unfiltered,
    "../scratch/R_outputs/InitialProcessing/CFC1/DO_unfiltered.RDS"
)
saveRDS(
    DO.p005,
    "../scratch/R_outputs/InitialProcessing/CFC1/DO_p005.RDS"
)

# Create directory and CSV outputs
dir.create("../results/2024_12_20_InitialProcessing/CFC1/GSEA_Enrichments/DO",
    recursive = TRUE
)

sapply(contrasts, FUN = function(x) {
    a <- DO.unfiltered@compareClusterResult %>%
        filter(Cluster == x) %>%
        dplyr::select(!Cluster)
    b <- a %>%
        filter(p.adjust < 0.05)

    write_csv(a, file = paste0("../results/2024_12_20_InitialProcessing/CFC1/GSEA_Enrichments/DO/DO-", x, "_AllEnrichments.csv"))
    write_csv(b, file = paste0("../results/2024_12_20_InitialProcessing/CFC1/GSEA_Enrichments/DO/DO-", x, "_SignificantEnrichments.csv"))
})

# %% Run GO functions
GO.unfiltered <- compareCluster(
    geneClusters = GSEA.Lists,
    fun = "gseGO",
    OrgDb = "org.Hs.eg.db",
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    nPermSimple = 10000,
    seed = 42,
    eps = 0,
    ont = "all"
) %>%
    setReadable("org.Hs.eg.db",
        keyType = "ENTREZID"
    )

GO.p005 <- compareCluster(
    geneClusters = GSEA.Lists,
    fun = "gseGO",
    pvalueCutoff = 0.05,
    OrgDb = "org.Hs.eg.db",
    pAdjustMethod = "BH",
    nPermSimple = 10000,
    seed = 42,
    eps = 0,
    ont = "all"
) %>%
    setReadable("org.Hs.eg.db",
        keyType = "ENTREZID"
    ) %>%
    pairwise_termsim(showCategory = Inf)

# Save these outputs
saveRDS(
    GO.unfiltered,
    "../scratch/R_outputs/InitialProcessing/CFC1/GO_unfiltered.RDS"
)
saveRDS(
    GO.p005,
    "../scratch/R_outputs/InitialProcessing/CFC1/GO_p005.RDS"
)

# Create directory structure for GO results
sapply(c("BP", "CC", "MF"), function(ont) {
    dir.create(paste0("../results/2024_12_20_InitialProcessing/CFC1/GSEA_Enrichments/GO/", ont),
        recursive = TRUE
    )
})

# Generate CSVs for each contrast and ontology
sapply(contrasts, FUN = function(x) {
    a <- GO.unfiltered@compareClusterResult %>%
        filter(Cluster == x)

    sapply(c("BP", "CC", "MF"), FUN = function(y) {
        data <- a %>%
            filter(ONTOLOGY == y) %>%
            dplyr::select(!c(Cluster, ONTOLOGY))
        data_sig <- data %>%
            filter(p.adjust < 0.05)

        write_csv(data,
            file = paste0(
                "../results/2024_12_20_InitialProcessing/CFC1/GSEA_Enrichments/GO/",
                y, "/GO_", y, "-", x, "_AllEnrichments.csv"
            )
        )
        write_csv(data_sig,
            file = paste0(
                "../results/2024_12_20_InitialProcessing/CFC1/GSEA_Enrichments/GO/",
                y, "/GO_", y, "-", x, "_SignificantEnrichments.csv"
            )
        )
    })
})

# %% KEGG
# Run KEGG enrichment analysis
KEGG.unfiltered <- compareCluster(
    geneClusters = GSEA.Lists,
    fun = "gseKEGG",
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    nPermSimple = 10000,
    seed = 42,
    eps = 0,
    organism = "hsa"
) %>%
    setReadable("org.Hs.eg.db",
        keyType = "ENTREZID"
    )

KEGG.p005 <- compareCluster(
    geneClusters = GSEA.Lists,
    fun = "gseKEGG",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    nPermSimple = 10000,
    seed = 42,
    eps = 0,
    organism = "hsa"
) %>%
    setReadable("org.Hs.eg.db",
        keyType = "ENTREZID"
    ) %>%
    pairwise_termsim(showCategory = Inf)

# Save these outputs
saveRDS(
    KEGG.unfiltered,
    "../scratch/R_outputs/InitialProcessing/CFC1/KEGG_unfiltered.RDS"
)
saveRDS(
    KEGG.p005,
    "../scratch/R_outputs/InitialProcessing/CFC1/KEGG_p005.RDS"
)

# Create directory and CSV outputs
dir.create("../results/2024_12_20_InitialProcessing/CFC1/GSEA_Enrichments/KEGG",
    recursive = TRUE
)

sapply(contrasts, FUN = function(x) {
    a <- KEGG.unfiltered@compareClusterResult %>%
        filter(Cluster == x) %>%
        dplyr::select(!Cluster)
    b <- a %>%
        filter(p.adjust < 0.05)

    write_csv(a, file = paste0("../results/2024_12_20_InitialProcessing/CFC1/GSEA_Enrichments/KEGG/KEGG-", x, "_AllEnrichments.csv"))
    write_csv(b, file = paste0("../results/2024_12_20_InitialProcessing/CFC1/GSEA_Enrichments/KEGG/KEGG-", x, "_SignificantEnrichments.csv"))
})

# %% Generate plots

ConstructPlots <- function(object, dir.name) {
    p <- dotplot(object,
        font.size = 16,
        color = "NES",
        label_format = 70,
        showCategory = 10
    ) +
        scale_fill_gradient2(
            low = "blue",
            mid = "white",
            high = "red"
        ) +
        theme(
            axis.text.x = element_text(
                angle = 45,
                hjust = 1,
                size = 14
            ),
            axis.text.y = element_text(size = 16)
        )

    ggsave(paste0("../results/2024_12_20_InitialProcessing/CFC1/GSEA_Enrichments/", dir.name, "/", dir.name, "_SummaryDotplot.pdf"),
        p,
        height = 8,
        width = 8
    )

    p2 <- emapplot(object,
        showCategory = 30,
        cex.params = list(
            category_label = 1,
            category_node = 3,
            line = 0.5
        ),
        pie.params = list(legend_n = 3)
    )

    ggsave(paste0("../results/2024_12_20_InitialProcessing/CFC1/GSEA_Enrichments/", dir.name, "/", dir.name, "_SummaryEnrichmentMap.pdf"),
        p2,
        height = 10,
        width = 10
    )
}

# %% Generate standard enrichment plots
ConstructPlots(object = DGN.p005, dir.name = "DGN")
ConstructPlots(object = KEGG.p005, dir.name = "KEGG")
ConstructPlots(object = DO.p005, dir.name = "DO")

# %% Generate GO-specific plots
ontologies <- c("BP", "CC", "MF")

sapply(ontologies, FUN = function(x) {
    temp <- GO.p005
    temp@compareClusterResult <- temp@compareClusterResult %>%
        filter(ONTOLOGY == x)

    p <- dotplot(temp,
        font.size = 16,
        color = "NES",
        label_format = 70,
        showCategory = 10
    ) +
        scale_fill_gradient2(
            low = "blue",
            mid = "white",
            high = "red"
        ) +
        theme(
            axis.text.x = element_text(
                angle = 45,
                hjust = 1,
                size = 14
            ),
            axis.text.y = element_text(size = 16)
        )

    ggsave(paste0("../results/2024_12_20_InitialProcessing/CFC1/GSEA_Enrichments/GO/", x, "/GO_", x, "-SummaryDotplot.pdf"),
        p,
        height = 10,
        width = 10
    )

    p2 <- emapplot(temp,
        showCategory = 30,
        cex.params = list(
            category_label = 1,
            category_node = 3,
            line = 0.5
        ),
        pie.params = list(legend_n = 3)
    )

    ggsave(paste0("../results/2024_12_20_InitialProcessing/CFC1/GSEA_Enrichments/GO/", x, "/GO_", x, "-SummaryEnrichmentMap.pdf"),
        p2,
        height = 10,
        width = 10
    )
})
