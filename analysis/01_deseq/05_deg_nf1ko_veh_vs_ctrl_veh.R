# 05_contrast_nf1_veh_vs_ctrl_veh.R
# Perform differential expression analysis for NF1KO_Veh vs CTRL_Veh contrast

suppressPackageStartupMessages({
    library(DESeq2)
    library(tidyverse)
})

dds <- readRDS("data/processed/deseq/dds_vehRef.rds")

res_nf1ko_veh <- results(
    dds,
    contrast = c("group", "NF1KO_Veh", "CTRL_Veh")
)

summary(res_nf1ko_veh)

# Export as CSV
deg_tbl_dir <- "results/deg/tables"
if (!dir.exists(deg_tbl_dir)) {
    dir.create(deg_tbl_dir, recursive = TRUE)
}
res_tbl <- as.data.frame(res_nf1ko_veh) %>%
    tibble::rownames_to_column("gene")

write.csv(res_tbl, file.path(deg_tbl_dir, "DEG_NF1KO_Veh_vs_CTRL_Veh.csv"), row.names = FALSE)
cat("Differential expression results for NF1KO_Veh vs CTRL_Veh saved.\n")

# Filter significant DEGs
deg_nf1ko_veh <- res_tbl %>%
    filter(!is.na(padj)) %>%
    filter(padj < 0.05, abs(log2FoldChange) >= 1)
cat("Number of DEGs (padj < 0.05 & |log2FC| >= 1):", nrow(deg_nf1ko_veh), "\n")

# Count up vs down
deg_summary <- deg_nf1ko_veh %>%
    mutate(direction = if_else(log2FoldChange > 0, "up", "down")) %>%
    count(direction)
    
print(deg_summary)

# Save DEG table
write.csv(
    deg_nf1ko_veh,
    file.path(deg_tbl_dir, "DEG_NF1KO_Veh_vs_CTRL_Veh_sig.csv"),
    row.names = FALSE
)
cat("Filtered DEG table for NF1KO_Veh vs CTRL_Veh saved.\n")
