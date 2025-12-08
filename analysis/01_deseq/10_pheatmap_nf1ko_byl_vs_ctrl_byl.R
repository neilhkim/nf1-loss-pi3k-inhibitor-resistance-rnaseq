#!/usr/bin/env Rscript

# 10_pheatmap_nf1ko_byl_vs_ctrl_byl.R
# Heatmap of top DE genes for NF1KO_BYL vs CTRL_BYL comparison
# Using the same convention used in 07)volcano_detailed.R (no assumptions about IDs)

suppressPackageStartupMessages({
    library(DESeq2)
    library(tidyverse)
    library(pheatmap)
})

set.seed(123) # for reproducibility (pheatmap clustering)

# 1. Paths -----------

vsd_file <- "data/processed/deseq/vsd.rds"
deg_file <- "results/deg/tables/DEG_NF1KO_BYL_vs_CTRL_BYL.csv"
gene_annot_file <- "data/processed/gene_annotations_entrez.rds"
out_fig_dir <- "results/heatmap/figures"
if (!dir.exists(out_fig_dir)) dir.create(out_fig_dir, recursive = TRUE)
out_png <- file.path(out_fig_dir, "heatmap_NF1KO_BYL_vs_CTRL_BYL_top50.png")
out_pdf <- file.path(out_fig_dir, "heatmap_NF1KO_BYL_vs_CTRL_BYL_top50.pdf")

# 2. Load data -----------

message("Loading VST object...")
vsd <- readRDS(vsd_file)

deg <- read_csv(
    deg_file,
    show_col_types = FALSE
) %>%
    filter(
        !is.na(padj),
        is.finite(padj),
        is.finite(log2FoldChange)
    ) %>%
    mutate(
        gene = as.character(gene),
        padj = ifelse(padj == 0, 1e-320, padj)  # avoid -Inf in -log10
    )

gene_annot <- readRDS(gene_annot_file) %>%
    mutate(ENTREZID = as.character(ENTREZID)) %>%
    select(ENTREZID, SYMBOL)


# 3. Select top DE genes by adjusted p-value -----------

top_n <- 50

top_deg <- deg %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    slice_head(n = top_n)

gene_ids <- top_deg$gene

message("Selected top ", nrow(top_deg), " DE genes by adjusted p-value for heatmap.")

# 4. Build expression matrix for these genes ----------------------------
# ---- Important: rownames(vsd) are Entrez IDs and match top_deg$gene ---

mat <- assay(vsd)[gene_ids, , drop = FALSE]

# 5. Make nicer gene labels (SYMBOL when available) ---------------------

mat_labels <- tibble(gene = rownames(mat)) %>%
    left_join(gene_annot, by = c("gene" = "ENTREZID")) %>%
    mutate(pretty = if_else(!is.na(SYMBOL), SYMBOL, gene))

rownames(mat) <- mat_labels$pretty

# 6. Row-scale matrix --------------------------------

mat_scaled <- t(scale(t(mat)))

na_rows <- which(apply(mat_scaled, 1, function(x) any(is.na(x))))
if (length(na_rows) > 0) {
    mat_scaled <- mat_scaled[-na_rows, , drop = FALSE]
    message("Removed ", length(na_rows), " genes with NA values after scaling.")
}

# 7. Prepare annotation for columns (samples) ------------------------

annot_col <- as.data.frame(colData(vsd)) %>%
    select(genotype, treatment, group)

# Ensure factor levels encode the desired logical order
annot_col <- annot_col %>%
    mutate(
        genotype = factor(
            genotype,
            levels = c("CTRL", "NF1KO")
        ),
        treatment = factor(
            treatment,
            levels = c("Veh", "BYL")
        ),
        group = factor(
            group,
            levels = c(
                "CTRL_Veh",
                "CTRL_BYL",
                "NF1KO_Veh",
                "NF1KO_BYL"
            )
        )
    )


# Reorder samples in both mat_scaled and annot_col
# Order by genotype then treatment then group
sample_order <- order(
    annot_col$genotype,
    annot_col$treatment,
    annot_col$group
)

mat_scaled <- mat_scaled[, sample_order, drop = FALSE]
annot_col <- annot_col[sample_order, , drop = FALSE]

ann_colors <- list(
    genotype = c(CTRL = "#1f78b4", NF1KO = "#e31a1c"),
    treatment = c(Veh = "#33a02c", BYL = "#ff7f00"),
    group = c(
        CTRL_Veh = "#a6cee3",
        CTRL_BYL = "#1f78b4",
        NF1KO_Veh = "#fb9a99",
        NF1KO_BYL = "#e31a1c"
    )
)

# Only keep relevant levels
for (nm in names(ann_colors)) {
    ann_colors[[nm]] <- ann_colors[[nm]][
        intersect(names(ann_colors[[nm]]), unique(annot_col[[nm]]))
    ]
}

# 8. Plot heatmap and save (PNG and PDF) --------------------------------

png(out_png, width = 1800, height = 2200, res = 200)
pheatmap(
    mat_scaled,
    annotation_col = annot_col,
    annotation_colors = ann_colors,
    cluster_cols = FALSE,
    clustering_distance_rows = "correlation",
    clustering_method = "complete",
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 6,
    fontsize_col = 8,
    main = "Top 50 DE genes (NF1KO_BYL vs CTRL_BYL)"
)
dev.off()
message("Saved PNG heatmap to: ", out_png)

pdf(out_pdf, width = 8, height = 10)
pheatmap(
    mat_scaled,
    annotation_col        = annot_col,
    annotation_colors     = ann_colors,
    cluster_cols = FALSE,
    clustering_distance_rows = "correlation",
    clustering_method     = "complete",
    show_rownames         = TRUE,
    show_colnames         = TRUE,
    fontsize_row          = 6,
    fontsize_col          = 8,
    main                  = "Top 50 DE genes (NF1KO_BYL vs CTRL_BYL)"
)
dev.off()
message("Saved PDF heatmap to: ", out_pdf)