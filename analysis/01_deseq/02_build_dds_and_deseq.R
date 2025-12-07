# 02_build_dds.R
# Build DESeq2 datasets and run DESeq

suppressPackageStartupMessages({
    library(tidyverse)
    library(DESeq2)
})
cat("[deseq] Starting DESeq2 dataset build and fit...\n")

# Load processed data
cat("[deseq] Loading counts matrix and sample metadata...\n")
counts_mat <- readRDS("data/processed/counts_matrix.rds")
sample_info <- readRDS("data/metadata/sample_info.rds")
cat("[deseq] Counts dims: ", paste(dim(counts_mat), collapse = " x "), "\n", sep = "")
cat("[deseq] Metadata rows: ", nrow(sample_info), "\n", sep = "")

# Ensure ordering matches
cat("[deseq] Ordering counts columns to match metadata sample_id...\n")
counts_mat <- counts_mat[, sample_info$sample_id]

sample_info$group <- factor(sample_info$group)
cat("[deseq] Groups: ", paste(levels(sample_info$group), collapse = ", "), "\n", sep = "")

# Relevel to set reference levels
# 1. Build DESeq2 dataset (Veh Reference) ----------

cat("[deseq] Building dataset with reference 'CTRL_Veh'...\n")
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
cat("[deseq] Prefiltered vehRef genes: ", nrow(dds_vehRef), "\n", sep = "")

# Run DESeq and save results
cat("[deseq] Running DESeq() for vehRef...\n")
dds_vehRef <- DESeq(dds_vehRef)

# Create output directory if it doesn't exist
if (!dir.exists("data/processed/deseq")) {
    dir.create("data/processed/deseq", recursive = TRUE)
}

saveRDS(dds_vehRef, "data/processed/deseq/dds_vehRef.rds")
cat("[deseq] Saved dds_vehRef to data/processed/deseq/dds_vehRef.rds\n")

 # 2. Build DESeq2 dataset (BYL Reference) ----------

cat("[deseq] Building dataset with reference 'CTRL_BYL'...\n")
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
cat("[deseq] Prefiltered bylRef genes: ", nrow(dds_bylRef), "\n", sep = "")

# Run DESeq and save results
cat("[deseq] Running DESeq() for bylRef...\n")
dds_bylRef <- DESeq(dds_bylRef)
saveRDS(dds_bylRef, "data/processed/deseq/dds_bylRef.rds")
cat("[deseq] Saved dds_bylRef to data/processed/deseq/dds_bylRef.rds\n")
cat("[deseq] All done.\n")