# 03_vst_and_pca.R
# Perform variance stabilizing transformation and PCA on DESeq2 dataset


suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
})


dds <- readRDS("data/processed/deseq/dds_vehRef.rds")

# Variance Stabilizing Transformation
vsd <- vst(dds)

if (!dir.exists("data/processed/deseq")){
    dir.create("data/processed/deseq", recursive = TRUE)
}
saveRDS(vsd, "data/processed/deseq/vsd.rds")

vsd <- readRDS("data/processed/deseq/vsd.rds")
coldata <- as.data.frame(colData(vsd))

pca <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))

out_fig_dir <- "results/pca/figures"
if (!dir.exists(out_fig_dir)){
    dir.create(out_fig_dir, recursive = TRUE)
}

p <- ggplot(pca, aes(PC1, PC2, color = group)) + 
    geom_point(size = 3) + 
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme_bw()
ggsave(file.path(out_fig_dir, "pca_group.png"), p, width = 5, height = 4)
cat("[deseq] VST and PCA completed. Saved VSD and PCA plot.\n")