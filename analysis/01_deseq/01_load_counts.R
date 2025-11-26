## 01_load_counts.R
## Download raw counts (if needed), read them, and build sample metadata for GSE207514

suppressPackageStartupMessages({
    library(tidyverse)
})

## 0. Paths and download helper ------------------------------------------

raw_dir <- "data/raw"
if (!dir.exists(raw_dir)) {
    dir.create(raw_dir, recursive = TRUE)
}

txt_path <- file.path(raw_dir, "GSE207514_RNAseq_T47D_CTRL_NF1KO_annotated_count_table.txt")
gz_path  <- paste0(txt_path, ".gz")

geo_urls <- c(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207514/suppl/GSE207514%5FRNAseq%5FT47D%5FCTRL%5FNF1KO%5Fannotated%5Fcount%5Ftable.txt.gz",
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE207514&format=file&file=GSE207514%5FRNAseq%5FT47D%5FCTRL%5FNF1KO%5Fannotated%5Fcount%5Ftable%2Etxt%2Egz"
)

download_counts_if_needed <- function() {
    if (file.exists(txt_path)) {
        message("Found existing count table: ", txt_path)
        return(invisible(NULL))
    }

    message("Count table not found at ", txt_path)
    message("Attempting to download from GEO and save under data/raw/")

    # Try to download gz if not present
    if (!file.exists(gz_path)) {
        success <- FALSE
        for (u in geo_urls) {
            message("Trying URL: ", u)
            res <- try(
                download.file(u, destfile = gz_path, mode = "wb", quiet = TRUE),
                silent = TRUE
            )
            if (!inherits(res, "try-error")) {
                success <- TRUE
                message("Downloaded gzipped count table to: ", gz_path)
                break
            }
        }
        if (!success) {
            stop("Failed to download count table from GEO. Please check your internet connection or download manually.")
        }
    } else {
        message("Found existing gzipped count table: ", gz_path)
    }

    # Read gzipped file and write plain txt for future runs
    message("Reading gzipped count table and writing plain text copy...")
    counts_tmp <- read.delim(gzfile(gz_path), check.names = TRUE)
    write.table(
        counts_tmp,
        file = txt_path,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    message("Wrote plain text count table to: ", txt_path)
}

## 1. Ensure counts file exists and load it -------------------------------

download_counts_if_needed()

counts_path <- txt_path

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

## 2. Build sample metadata -----------------------------------------------

# Sample IDs are:
# CTR_Veh_1, CTR_Veh_2, ...,
# CTR_BYL_1, ...,
# NF1_Veh_1, ...,
# NF1_BYL_1, ...

sample_info <- tibble(
    sample_id = sample_ids
) %>%
    mutate(
        # Genotype: CTRL vs NF1KO
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

## 2a. Simple sanity checks -----------------------------------------------

if (any(is.na(sample_info$genotype))) {
    stop("Some samples have undefined genotype.")
}

if (any(is.na(sample_info$treatment))) {
    stop("Some samples have undefined treatment.")
}

## 3. Save outputs --------------------------------------------------------

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