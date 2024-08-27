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
    library(jsonlite)  # Add this library to handle JSON strings
  })
})

### Step 1 - Read in the DGEList object

DGE <- readRDS(DGE)

### Step 2 - Construct contrast matrix

design <- model.matrix(data = DGE$samples,
                       ~0 + genotype_clean)
colnames(design) <- str_remove_all(colnames(design),
                                   "genotype_clean")
contrast.matrix <- makeContrasts(
    KO = "GNASknockout - WT",
    levels = colnames(design))

### Step 3 - Perform the DEG analysis

v <- voom(DGE,
          design)
vfit <- lmFit(v,
              design)
vfit <- contrasts.fit(vfit,
                      contrast.matrix)
efit <- eBayes(vfit)

### Step 4 - Record results of DEG analysis

contrasts <- colnames(contrast.matrix)

LFC.summary <- sapply(contrasts, function(x){
    lfc.list <- list()
    top <- topTable(efit,
                    coef = x,
                    number = Inf) %>%
    list()
    lfc.list <- append(lfc.list, top)
    })

LFC.summary[[1]] %>%
arrange(adj.P.Val) # This is a dataframe