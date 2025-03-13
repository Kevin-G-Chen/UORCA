# %%
library(tidyverse)
library(limma)
library(edgeR)
library(tximport)
getwd()


# %%
tx2gene <- read_tsv("../../data/kallisto_indices/human/t2g.txt",
    col_names = FALSE
) %>%
    dplyr::select(1, 3) %>%
    drop_na()
files <- list.files(
    pattern = "abundance.tsv",
    path = "../3_analyse_data/Benchmark_Kallisto/GSE133702_Kallisto",
    recursive = TRUE,
    full.names = TRUE
)
kallist <- tximport(
    files = files,
    type = "kallisto",
    tx2gene = tx2gene,
    ignoreAfterBar = TRUE,
    countsFromAbundance = "lengthScaledTPM"
)
# %% Progress with metadata
meta <- read_csv("../2_extract_data/GSE133702_data/GSE133702_series_matrix_metadata.csv")
DGE <- DGEList(
    counts = kallist$counts,
    samples = meta
)
keep.exprs <- filterByExpr(DGE, group = DGE$samples$title)
DGE.filtered <- DGE
