# 04_vst_and_pca.R

library(DESeq2)
library(ggplot2)

dds <- readRDS("data/processed/deseq/dds_vehRef.rds")

# Variance Stabilizing Transformation
vsd <- vst(dds)

if (!dir.exists("results/deseq")){
    dir.create("results/deseq", recursive = TRUE)
}
saveRDS(vsd, "results/deseq/vsd.rds")

vsd <- readRDS("results/deseq/vsd.rds")
coldata <- as.data.frame(colData(vsd))

pca <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))

if (!dir.exists("results/figures")){
    dir.create("results/figures", recursive = TRUE)
}

p <- ggplot(pca, aes(PC1, PC2, color = group)) + 
    geom_point(size = 3) + 
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme_bw()

ggsave("results/figures/pca_group.png", p, width = 5, height = 4)