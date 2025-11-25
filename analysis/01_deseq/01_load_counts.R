## 01_load_counts.R
## Read raw counts and build sample metadata for GSE207514

library(tidyverse)

## 1. Read the count table ------

counts_path <- "data/raw/GSE207514_RNAseq_T47D_CTRL_NF1KO_annotated_count_table.txt"

# Read as a data frame
counts_raw <- read.delim(counts_path, check.names = TRUE)

cat("Column names in counts_raw:\n")
print(colnames(counts_raw))

# The first three columns are gene annotations:
# ENTREZID, SYMBOL, GENENAME
annotation_cols <- 1:3

gene_annot <- counts_raw[, annotation_cols, drop = FALSE]

# Use ENTREZID as identifier
gene_ids <- counts_raw$ENTREZID

# Counts matrix
counts_mat <- counts_raw[, -annotation_cols, drop = FALSE]
rownames(counts_mat) <- gene_ids

sample_ids <- colnames(counts_mat)

cat("\nSample IDs:\n")
print(sample_ids)

## 2. Build sample metadata ------

# Sample IDs are: 
# CTR_Veh_1, CTR_Veh_2, ..., 
# CTR_BYL_1, ..., 
# NF1_Veh_1, ..., 
# NF1_BYL_1, ...

sample_info <- tibble(
    sample_id = sample_ids
) %>%
    mutate(
        # Genotype: CTR vs NF1
        genotype = case_when(
            str_detect(sample_id, "^CTR") ~ "CTRL",
            str_detect(sample_id, "^NF1") ~ "NF1KO",
            TRUE ~ NA_character_
        ),
        # Treatment: Veh vs BYL
        treatment = case_when(
            str_detect(sample_id, "_Veh_") ~ "Veh",
            str_detect(sample_id, "_BYL_") ~ "BYL",
            TRUE ~ NA_character_
        ),
        # Group: e.g., CTRL_Veh, NF1KO_BYL, etc.
        group = paste(genotype, treatment, sep = "_")
    )

cat("\nSample metadata:\n")
print(sample_info)

## 2a. Simple sanity checks ----

if (any(is.na(sample_info$genotype))) {
    stop("Some samples have undefined genotype.")
}

if (any(is.na(sample_info$treatment))) {
    stop("Some samples have undefined treatment.")
}

## 3. Save outputs ------

# Create output folders if needed
if (!dir.exists("data/metadata")) {
    dir.create("data/metadata", recursive = TRUE)
}
if (!dir.exists("data/processed")) {
    dir.create("data/processed", recursive = TRUE)
}

# Save metadata
write_csv(sample_info, "data/metadata/sample_info.csv")
saveRDS(sample_info, "data/metadata/sample_info.rds")

# Save counts matrix and gene annotations
saveRDS(counts_mat, "data/processed/counts_matrix.rds")
saveRDS(gene_annot, "data/processed/gene_annotations.rds")

cat("\nData loading and processing complete. Outputs saved.\n")