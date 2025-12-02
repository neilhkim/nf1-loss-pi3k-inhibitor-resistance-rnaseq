# 06_contrast_nf1ko_byl_vs_ctrl_byl.R
# Perform differential expression analysis for NF1KO_BYL vs CTRL_BYL contrast

library(DESeq2)
library(tidyverse)

dds <- readRDS("data/processed/deseq/dds_bylRef.rds")

res_nf1ko_byl <- results(
    dds,
    contrast = c("group", "NF1KO_BYL", "CTRL_BYL")
)

summary(res_nf1ko_byl)

# Export as CSV
if (!dir.exists("results/tables")) {
    dir.create("results/tables", recursive = TRUE)
}
res_tbl <- as.data.frame(res_nf1ko_byl) %>%
    tibble::rownames_to_column("gene")

write.csv(res_tbl, "results/tables/DEG_NF1KO_BYL_vs_CTRL_BYL.csv", row.names = FALSE)

cat("Differential expression results for NF1KO_BYL vs CTRL_BYL saved.\n")

# Filter significant DEGs
deg_nf1ko_byl <- res_tbl %>%
    filter(!is.na(padj)) %>%
    filter(padj < 0.05, abs(log2FoldChange) >= 1)

cat("Number of DEGs (padj < 0.05 & |log2FC| >= 1):", nrow(deg_nf1ko_byl), "\n")
# Count up vs down
deg_summary <- deg_nf1ko_byl %>%
    mutate(direction = if_else(log2FoldChange > 0, "up", "down")) %>%
    count(direction)

print(deg_summary)

# Save DEG table
write.csv(
    deg_nf1ko_byl,
    "results/tables/DEG_NF1KO_BYL_vs_CTRL_BYL_sig.csv",
    row.names = FALSE
)
cat("Filtered DEG table for NF1KO_BYL vs CTRL_BYL saved.\n")