#!/usr/bin/env Rscript

# 06_summarize_counts.R
# Combine featureCounts output into a clean counts matrix
# using runs.tsv to map SRR run IDs to sample IDs (CTR_Veh_1, etc.)

suppressPackageStartupMessages({
    library(tidyverse)
})

cat("[summarize] Starting summarization of featureCounts matrix.\n")

# Input and output paths
fc_file <- "results/featurecounts/featurecounts_gene_counts.txt"
runs_tsv <- "data/metadata/runs.tsv"
out_rds <- "results/counts_matrix_generated.rds"
out_csv <- "results/counts_matrix_generated.csv"

# Step 1: read counts
cat("[summarize] Reading featureCounts file:\n", fc_file, "\n")
fc <- read.delim(fc_file, comment.char = "#", check.names = FALSE)

# Inspect structure
cat("[summarize] FeatureCounts table loaded. Dimensions: ",
    nrow(fc), " genes x ", ncol(fc), "columns.\n")

# Step 2: remove annotation columns
annot_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")

if (!all(annot_cols %in% colnames(fc))) {
    stop("[summarize] ERROR: annotation columns missing. Check input file.")
}

counts <- fc %>%
    select(-all_of(annot_cols))

genes <- fc$Geneid
rownames(counts) <- genes

# Step 3: read sample metadata and reorder columns
cat("[summarize] Reading run table:\n", runs_tsv, "\n")
runs <- read.delim(runs_tsv)

# Expect these columns
expected_runs_cols <- c("run_id", "sample_id", "group", "fastq")
if (!all(expected_runs_cols %in% colnames(runs))) {
    stop("[summarize] ERROR: runs.tsv missing required column: ",
        paste(expected_runs_cols, collapse = ", "))
}

# Column names in counts are BAM paths such as:
#   results/star/SRR19987593_Aligned.sortedByCoord.out.bam

bam_paths <- colnames(counts)

# Extract run_id (SRR accession) from each column
run_from_bam <- basename(bam_paths) %>%
    str_replace("_Aligned.sortedByCoord.out.bam$", "")

# Rename columns of counts to run_id so they match run$run_id
colnames(counts) <- run_from_bam

# Order runs_sub to match the order of columns in counts
runs_sub <- runs %>%
    filter(run_id %in% run_from_bam) %>%
    mutate(col_order = match(run_id, run_from_bam)) %>%
    arrange(col_order)

# Sanify check
if (!all(run_from_bam %in% runs_sub$run_id)) {
    missing <- setdiff(run_from_bam, runs_sub$run_id)
    stop("[summarize] ERROR: Some run_ids from featureCounts not found in runs.tsv:\n ",
        paste(missing, collapse = ", "))
}

# Reorder columns of counts to that same order
counts <- counts[, runs_sub$run_id]

# Rename columns to sample_id, so you get CTR_Veh_1 etc.
colnames(counts) <- runs_sub$sample_id

cat("[summarize] Counts matrix after mapping SRR -> sample_id: ",
    nrow(counts), " genes x ", ncol(counts), " samples.\n")

cat("[summarize] Sample columns:\n",
    paste(colnames(counts), collapse = ", "), "\n")


# Step 4: save results
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)

cat("[summarize] Saving counts matrix RDS to:\n", out_rds, "\n")
saveRDS(counts, out_rds)

cat("[summarize] Saving counts matrix CSV to:\n", out_csv, "\n")
write.csv(counts, out_csv, quote = FALSE)

cat("[summarize] Done. Final matrix: ",
    nrow(counts), " genes x ", ncol(counts), " samples.\n")
