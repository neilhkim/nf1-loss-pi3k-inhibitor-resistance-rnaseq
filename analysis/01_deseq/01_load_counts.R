#!/usr/bin/env Rscript
## 01_load_counts.R
## Load counts and build sample metadata for GSE207514
## - Default: try pipeline-generated counts first; if not found, use GEO.
## - Options:
##     --source=auto      (default)
##     --source=geo       (force GEO; will prompt about overwriting existing file)

suppressPackageStartupMessages({
    library(tidyverse)
})

cat("[load] Starting count loading script.\n")

## ------------------------------------------------------------------
## 0. Parse command line arguments
## ------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
source_mode <- "auto"

for (a in args) {
    if (grepl("^--source=", a)) {
        source_mode <- sub("^--source=", "", a)
    }
}

if (!source_mode %in% c("auto", "geo")) {
    stop("[load] Invalid --source option. Use auto or geo.")
}

cat("[load] Source mode = ", source_mode, "\n")

## ------------------------------------------------------------------
## 1. Paths and GEO download helper
## ------------------------------------------------------------------

raw_dir <- "data/raw"
if (!dir.exists(raw_dir)) {
    dir.create(raw_dir, recursive = TRUE)
}

txt_path <- file.path(
    raw_dir,
    "GSE207514_RNAseq_T47D_CTRL_NF1KO_annotated_count_table.txt"
)
gz_path  <- paste0(txt_path, ".gz")

geo_urls <- c(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207514/suppl/GSE207514%5FRNAseq%5FT47D%5FCTRL%5FNF1KO%5Fannotated%5Fcount%5Ftable.txt.gz",
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE207514&format=file&file=GSE207514%5FRNAseq%5FT47D%5FCTRL%5FNF1KO%5Fannotated%5Fcount%5Ftable%2Etxt%2Egz"
)

download_geo_table <- function(force_download = FALSE, ask_overwrite = TRUE) {
    if (file.exists(txt_path)) {
        if (force_download) {
            do_overwrite <- TRUE
            if (ask_overwrite) {
                cat("[load] GEO txt exists at:\n        ", txt_path, "\n")
                cat("[load] Overwrite existing GEO txt? [y/N]: ")
                flush.console()
                ans <- tryCatch({
                    readLines(con = "stdin", n = 1)
                }, error = function(e) "")
                do_overwrite <- tolower(trimws(ans)) %in% c("y", "yes")
            }
            if (!do_overwrite) {
                message("[load] Keeping existing GEO txt. Skipping download.")
                return(invisible(NULL))
            }
        } else {
            message("[load] Found existing GEO txt table:\n        ", txt_path)
            return(invisible(NULL))
        }
    }

    message("[load] Downloading gzipped GEO table...")
    success <- FALSE
    for (u in geo_urls) {
        message("[load] Trying URL: ", u)
        res <- try(
            download.file(u, destfile = gz_path, mode = "wb", quiet = TRUE),
            silent = TRUE
        )
        if (!inherits(res, "try-error")) {
            success <- TRUE
            message("[load] Downloaded gzipped count table to:\n        ", gz_path)
            break
        }
    }
    if (!success) {
        stop("[load] Failed to download GEO count table. Please check your internet connection or download manually.")
    }

    message("[load] Reading gzipped table and writing plain text copy...")
    counts_tmp <- read.delim(gzfile(gz_path), check.names = TRUE)
    write.table(
        counts_tmp,
        file = txt_path,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    message("[load] Wrote plain text GEO table to:\n        ", txt_path)
}

## ------------------------------------------------------------------
## 2. Loaders for GEO and pipeline counts
## ------------------------------------------------------------------

load_geo_counts <- function(force_download = FALSE) {
    download_geo_table(force_download = force_download, ask_overwrite = TRUE)

    counts_raw <- read.delim(txt_path, check.names = TRUE)

    cat("[load] GEO table column names:\n")
    print(colnames(counts_raw))

    # First 3 columns are annotations: ENTREZID, SYMBOL, GENENAME
    annotation_cols <- 1:3
    gene_annot <- counts_raw[, annotation_cols, drop = FALSE]

    # Use ENTREZID as id
    gene_ids <- counts_raw$ENTREZID

    counts_mat <- counts_raw[, -annotation_cols, drop = FALSE]
    rownames(counts_mat) <- gene_ids

    list(counts_mat = counts_mat, gene_annot = gene_annot)
}

load_pipeline_counts <- function() {
    cp <- "data/processed/counts_matrix_entrez.rds"
    ga <- "data/processed/gene_annotations_entrez.rds"

    missing <- c()
    if (!file.exists(cp)) missing <- c(missing, cp)
    if (!file.exists(ga)) missing <- c(missing, ga)

    if (length(missing) > 0) {
        stop(
            "[load] Pipeline files not found:\n       ",
            paste(missing, collapse = "\n       ")
        )
    }

    counts_mat <- readRDS(cp)
    gene_annot <- readRDS(ga)

    cat("[load] Pipeline counts dimensions: ",
        nrow(counts_mat), " genes x ", ncol(counts_mat), " samples.\n")

    list(counts_mat = counts_mat, gene_annot = gene_annot)
}

## (Comparison logic removed; see scripts/counts_pipeline/08_compare_pipeline_vs_geo_visual.R)

## ------------------------------------------------------------------
## 4. Decide which source to use and load counts
## ------------------------------------------------------------------

cat("[load] Determining source of count matrix.\n")

counts_mat <- NULL
gene_annot <- NULL

if (source_mode == "geo") {
    cat("[load] Will use:  GEO (forced)\n")
    geo <- load_geo_counts(force_download = TRUE)
    counts_mat <- geo$counts_mat
    gene_annot <- geo$gene_annot

} else {  # auto
    cat("[load] Will use:  auto\n")
    cat("[load] Trying pipeline-generated counts first.\n")

    pipeline_try <- try(load_pipeline_counts(), silent = TRUE)
    if (!inherits(pipeline_try, "try-error")) {
        cat("[load] Successfully loaded pipeline counts. Using pipeline source.\n")
        counts_mat <- pipeline_try$counts_mat
        gene_annot <- pipeline_try$gene_annot
    } else {
        cat("[load] Pipeline counts not available. Falling back to GEO (with overwrite prompt if exists).\n")
        geo <- load_geo_counts(force_download = FALSE)
        counts_mat <- geo$counts_mat
        gene_annot <- geo$gene_annot
    }
}

if (is.null(counts_mat) || is.null(gene_annot)) {
    stop("[load] Failed to load counts and gene annotations from any source.")
}

cat("[load] Final counts dimensions: ",
    nrow(counts_mat), " genes x ", ncol(counts_mat), " samples.\n")

## ------------------------------------------------------------------
## 5. Build sample metadata (works for GEO and pipeline naming)
## ------------------------------------------------------------------

sample_ids <- colnames(counts_mat)

cat("[load] Sample IDs:\n")
print(sample_ids)

sample_info <- tibble(
    sample_id = sample_ids
) %>%
    mutate(
        # Genotype: CTRL vs NF1KO
        genotype = case_when(
            str_detect(sample_id, "^(CTR|CTRL)")   ~ "CTRL",
            str_detect(sample_id, "^(NF1|NF1KO)") ~ "NF1KO",
            TRUE ~ NA_character_
        ),
        # Treatment: Veh vs BYL
        treatment = case_when(
            str_detect(sample_id, "Veh") ~ "Veh",
            str_detect(sample_id, "BYL") ~ "BYL",
            TRUE ~ NA_character_
        ),
        # Group: e.g. CTRL_Veh, NF1KO_BYL
        group = paste(genotype, treatment, sep = "_")
    )

cat("\n[load] Sample metadata preview:\n")
print(sample_info)

## Simple sanity checks
if (any(is.na(sample_info$genotype))) {
    stop("[load] Some samples have undefined genotype.")
}

if (any(is.na(sample_info$treatment))) {
    stop("[load] Some samples have undefined treatment.")
}

## ------------------------------------------------------------------
## 6. Save outputs for downstream scripts
## ------------------------------------------------------------------

if (!dir.exists("data/metadata")) {
    dir.create("data/metadata", recursive = TRUE)
}
if (!dir.exists("data/processed")) {
    dir.create("data/processed", recursive = TRUE)
}

# Save metadata
readr::write_csv(sample_info, "data/metadata/sample_info.csv")
saveRDS(sample_info, "data/metadata/sample_info.rds")

# Save counts matrix and gene annotations in standard locations
saveRDS(counts_mat, "data/processed/counts_matrix.rds")
saveRDS(gene_annot, "data/processed/gene_annotations.rds")

cat("\n[load] Data loading and processing complete. Outputs saved.\n")
