#!/usr/bin/env Rscript
# 00_make_runs_tsv.R
# Create runs.tsv mapping SRR run IDs to sample IDs and metadata

library(tidyverse)

# Input file
sra <- read_csv("data/metadata/SraRunTable.csv", show_col_types = FALSE)

# Expected columns:
# Run, genotype, treatment, Sample Name, BioSample

# Clean treatment labels to match DESeq groups
sra <- sra %>%
    mutate(
        treatment_simple = case_when(
            treatment == "BYL719" ~ "BYL",
            treatment == "DMSO" ~ "Veh",
            TRUE ~ treatment
        ),
        group = paste(genotype, treatment_simple, sep = "_")
    )

# Assign replicate numbers per genotype-treatment group
sra <- sra %>%
    group_by(genotype, treatment_simple) %>%
    arrange(Run, .by_group = TRUE) %>%
    mutate(
        rep = row_number(),
        sample_id = paste0(genotype, "_", treatment_simple, "_rep", rep)
    ) %>%
    ungroup()

# Build final runs table
runs <- sra %>%
    transmute(
        run_id = Run,
        sample_id = sample_id,
        genotype = genotype,
        treatment = treatment_simple,
        group = group,
        fastq = file.path("dta/fastq", paste0(Run, ".fastq.gz"))
    )

# Write to TSV
write.table(
    runs,
    file = "data/metadata/runs.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat("Wrote data/metadata/runs.tsv with", nrow(runs), "entries.\n")