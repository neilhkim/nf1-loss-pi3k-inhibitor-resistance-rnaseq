#!/usr/bin/env Rscript
# 07_convert_ensembl_to_entrez.R
# Convert pipeline Ensembl based counts_matrix_generated.rds
# to ENTREZID and SYMBOL based versions for DESeq2

suppressPackageStartupMessages({
  library(tidyverse)
  library(org.Hs.eg.db)
  # library(AnnotationDbi)  # removed to avoid select() conflicts
})

cat("[convert] Loading pipeline Ensembl-based counts (data/processed/counts_matrix_ensembl.rds)...\n")
counts_ens <- readRDS("data/processed/counts_matrix_ensembl.rds")

cat("[convert] Matrix dimensions: ",
    nrow(counts_ens), " genes x ", ncol(counts_ens), " samples.\n")

# rownames(counts_ens) should look like ENSG00000... with version:
#   ENSG00000290825.1, ENSG00000223972.6, ...
ensembl_ids <- rownames(counts_ens)

cat("[convert] Example rownames before stripping:\n")
print(head(ensembl_ids))

cat("[convert] Stripping Ensembl version numbers (e.g., .1, .6)...\n")
ensembl_stripped <- sub("\\.\\d+$", "", ensembl_ids)

cat("[convert] Example rownames after stripping:\n")
print(head(ensembl_stripped))

# Put counts into a data frame with Ensembl ID columns
counts_df <- counts_ens %>%
  as.data.frame() %>%
  rownames_to_column("ENSEMBL_ORIG") %>%
  mutate(ENSEMBL = ensembl_stripped)

cat("[convert] Mapping Ensembl → ENTREZID + SYMBOL using org.Hs.eg.db...\n")
map_df <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(counts_df$ENSEMBL),
  keytype = "ENSEMBL",
  columns = c("ENTREZID", "SYMBOL")
)

# Some Ensembl IDs have no mapping. Keep only those with ENTREZID.
counts_mapped <- counts_df %>%
  left_join(map_df, by = c("ENSEMBL" = "ENSEMBL")) %>%
  filter(!is.na(ENTREZID))

cat("[convert] After mapping, rows remaining (with ENTREZID): ",
    nrow(counts_mapped), "\n")

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)

# 1) ENTREZ based counts matrix ---------------------------------------------
counts_entrez <- counts_mapped %>%
  dplyr::select(ENTREZID, all_of(colnames(counts_ens))) %>%
  dplyr::group_by(ENTREZID) %>%
  dplyr::summarise(across(where(is.numeric), sum), .groups = "drop") %>%
  column_to_rownames("ENTREZID") %>%
  as.matrix()

saveRDS(counts_entrez, "data/processed/counts_matrix_entrez.rds")
cat("[convert] Saved ENTREZ-based matrix to data/processed/counts_matrix_entrez.rds\n")

# 2) ENTREZ based annotation table ------------------------------------------
gene_annot_entrez <- counts_mapped %>%
  dplyr::select(ENTREZID, SYMBOL) %>%
  distinct() %>%
  mutate(GENENAME = NA_character_) %>%
  arrange(as.numeric(ENTREZID))

saveRDS(gene_annot_entrez, "data/processed/gene_annotations_entrez.rds")
cat("[convert] Saved ENTREZ-based gene annotations to data/processed/gene_annotations_entrez.rds\n")

# 3) SYMBOL based counts matrix (optional) ----------------------------------
counts_symbol <- counts_mapped %>%
  dplyr::select(SYMBOL, all_of(colnames(counts_ens))) %>%
  filter(!is.na(SYMBOL)) %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarise(across(where(is.numeric), sum), .groups = "drop") %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()

saveRDS(counts_symbol, "data/processed/counts_matrix_symbol.rds")
cat("[convert] Saved SYMBOL-based matrix to data/processed/counts_matrix_symbol.rds\n")

cat("[convert] Done.\n")
