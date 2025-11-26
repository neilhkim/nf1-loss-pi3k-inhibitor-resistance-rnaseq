## 02_build_dds.R

library(tidyverse)
library(DESeq2)

# Load processed data
counts_mat <- readRDS("data/processed/counts_matrix.rds")
sample_info <- readRDS("data/metadata/sample_info.rds")

# Ensure ordering matches
counts_mat <- counts_mat[, sample_info$sample_id]

sample_info$group <- factor(sample_info$group)

# Relevel to set reference levels
# 1. Build DESeq2 dataset (Veh Reference) ----------

sample_info_vehRef <- sample_info

# Set CTRL_Veh as reference
sample_info_vehRef$group <- relevel(sample_info_vehRef$group, ref = "CTRL_Veh")

# Build DESeqDataSet
dds_vehRef <- DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData = sample_info_vehRef,
    design = ~ group
)

# Filter low counts (standard DESeq2 prefilter)
dds_vehRef <- dds_vehRef[rowSums(counts(dds_vehRef) >= 10) >= 3, ] # Keep genes with at least 10 counts in at least 3 samples

# Run DESeq and save results
dds_vehRef <- DESeq(dds_vehRef)

# Create output directory if it doesn't exist
if (!dir.exists("data/processed/deseq")) {
    dir.create("data/processed/deseq", recursive = TRUE)
}

saveRDS(dds_vehRef, "data/processed/deseq/dds_vehRef.rds")
cat("DESeq() finished and dds_vehRef.rds saved.\n")

 # 2. Build DESeq2 dataset (BYL Reference) ----------

sample_info_bylRef <- sample_info

# Set CTRL_BYL as reference
sample_info_bylRef$group <- relevel(sample_info_bylRef$group, ref = "CTRL_BYL")

# Build DESeqDataSet
dds_bylRef <- DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData = sample_info_bylRef,
    design = ~ group
)

# Filter low counts (standard DESeq2 prefilter)
dds_bylRef <- dds_bylRef[rowSums(counts(dds_bylRef) >= 10) >= 3, ] # Keep genes with at least 10 counts

# Run DESeq and save results
dds_bylRef <- DESeq(dds_bylRef)
saveRDS(dds_bylRef, "data/processed/deseq/dds_bylRef.rds")
cat("DESeq() finished and dds_bylRef.rds saved.\n")