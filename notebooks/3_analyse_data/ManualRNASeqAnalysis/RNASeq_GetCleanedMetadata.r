DGE <- readRDS("path/to/DGE.RDS")

meta <- DGE$samples
print("Metadata columns: ", colnames(meta))