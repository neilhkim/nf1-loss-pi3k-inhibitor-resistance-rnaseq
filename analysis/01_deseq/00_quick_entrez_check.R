#!/usr/bin/env Rscript
# 00_quick_entrez_check.R
# Minimal check: do our ENTREZ counts resemble the GEO ENTREZ counts?

suppressPackageStartupMessages({
  library(tidyverse)
})

cat("[quick] Loading pipeline ENTREZ counts...\n")
counts_entrez <- readRDS("data/processed/counts_matrix_entrez.rds")

cat("[quick] Pipeline ENTREZ matrix dimensions: ",
    nrow(counts_entrez), " genes x ", ncol(counts_entrez), " samples.\n")

cat("[quick] Loading GEO counts table...\n")
counts_geo <- read.delim(
  "data/raw/GSE207514_RNAseq_T47D_CTRL_NF1KO_annotated_count_table.txt",
  check.names = FALSE
)

cat("[quick] GEO matrix dimensions: ",
    nrow(counts_geo), " genes x ", ncol(counts_geo), " columns.\n")

# GEO has ENTREZID, SYMBOL, GENENAME plus count columns.
# We convert ENTREZID to character to match rownames(counts_entrez).
counts_geo <- counts_geo %>%
  mutate(ENTREZID = as.character(ENTREZID))

entrez_pipeline <- rownames(counts_entrez)
entrez_geo      <- counts_geo$ENTREZID

cat("[quick] Unique ENTREZ in pipeline: ", length(unique(entrez_pipeline)), "\n")
cat("[quick] Unique ENTREZ in GEO:      ", length(unique(entrez_geo)), "\n")

# 1. Find shared ENTREZ IDs
shared <- intersect(entrez_pipeline, entrez_geo)

cat("[quick] Shared ENTREZ IDs: ", length(shared), "\n")

if (length(shared) == 0) {
  stop("[quick] No shared ENTREZ IDs between pipeline and GEO. Something is wrong with mapping.")
}

# 2. Subset and align rows by ENTREZ ID
cat("[quick] Subsetting and aligning matrices by shared ENTREZ...\n")

# Pipeline side: subset rows by shared IDs and put them in the same order as 'shared'
ours_sub <- counts_entrez[shared, , drop = FALSE]

# GEO side: subset rows and reorder to match 'shared'
geo_sub <- counts_geo %>%
  filter(ENTREZID %in% shared) %>%
  arrange(match(ENTREZID, shared))

# Sanity check that row orders match
if (!all(geo_sub$ENTREZID == shared)) {
  stop("[quick] Row alignment by ENTREZID failed. Check matching logic.")
}

# 3. Drop annotation columns from GEO and keep only count columns
geo_count_cols <- setdiff(colnames(geo_sub), c("ENTREZID", "SYMBOL", "GENENAME"))
geo_counts_mat <- as.matrix(geo_sub[, geo_count_cols])

cat("[quick] GEO count columns used: ",
    paste(geo_count_cols, collapse = ", "), "\n")

# 4. Compute per-gene total counts in each source
cat("[quick] Computing per-gene total counts for shared genes...\n")

gene_totals_ours <- rowSums(ours_sub)
gene_totals_geo  <- rowSums(geo_counts_mat)

# 5. Correlation of gene totals
cor_val <- cor(gene_totals_ours, gene_totals_geo)

cat("[quick] Pearson correlation of gene totals (pipeline vs GEO): ",
    sprintf("%.4f", cor_val), "\n")

# Optional quick summary of ranges
cat("[quick] Pipeline gene total range: ",
    range(gene_totals_ours), "\n")
cat("[quick] GEO gene total range:      ",
    range(gene_totals_geo), "\n")

cat("[quick] Done.\n")
