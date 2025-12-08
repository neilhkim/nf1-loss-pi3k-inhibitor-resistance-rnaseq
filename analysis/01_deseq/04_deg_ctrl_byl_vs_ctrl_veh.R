# 04_contrast_ctrl_byl_vs_ctrl_veh
# Perform differential expression analysis for CTRL_BYL vs CTRL_Veh contrast

suppressPackageStartupMessages({
    library(DESeq2)
    library(tidyverse)
})

dds <- readRDS("data/processed/deseq/dds_vehRef.rds")

res_ctrl_byl <- results(
    dds,
    contrast = c("group", "CTRL_BYL", "CTRL_Veh")
)

summary(res_ctrl_byl)

# Export as CSV
res_tbl <- as.data.frame(res_ctrl_byl) %>%
    tibble::rownames_to_column("gene")

deg_tbl_dir <- "results/deg/tables"
if (!dir.exists(deg_tbl_dir)){
    dir.create(deg_tbl_dir, recursive = TRUE)
}

write.csv(res_tbl, file.path(deg_tbl_dir, "DEG_CTRL_BYL_vs_CTRL_Veh.csv"), row.names = FALSE)
cat("Differential expression results for CTRL_BYL vs CTRL_Veh saved.\n")

# Filter significant DEGs
deg_ctrl_byl <- res_tbl %>%
    filter(!is.na(padj)) %>%
    filter(padj < 0.05, abs(log2FoldChange) >= 1)

cat("Number of DEGs (padj < 0.05 & |log2FC| >= 1):", nrow(deg_ctrl_byl), "\n")

# Count up vs down
deg_summary <- deg_ctrl_byl %>%
    mutate(direction = if_else(log2FoldChange > 0, "up", "down")) %>%
    count(direction)

print(deg_summary)

# Save DEG tale
write.csv(
    deg_ctrl_byl,
    file.path(deg_tbl_dir, "DEG_CTRL_BYL_vs_CTRL_Veh_sig.csv"),
    row.names = FALSE
)
cat("Filtered DEG table for CTRL_BYL vs CTRL_Veh saved.\n")