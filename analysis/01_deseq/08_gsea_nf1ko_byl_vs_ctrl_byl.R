# 09_gsea_nf1ko_byl_vs_ctrl_byl.R
# GSEA for NF1KO_BYL vs CTRL_BYL (resistance contrast)

suppressPackageStartupMessages({
    library(DESeq2)
    library(tidyverse)
    library(clusterProfiler)
    library(enrichplot)
    library(org.Hs.eg.db)
    library(msigdbr)
    library(dplyr)
})

set.seed(42)

# 0. Output folders ------------------

if (!dir.exists("results/tables")) {
    dir.create("results/tables", recursive = TRUE)
}

if (!dir.exists("results/figures")) {
    dir.create("results/figures", recursive = TRUE)
}

# 1. Load fitted DESeq2 object ------------------

dds <- readRDS("data/processed/deseq/dds_bylRef.rds")

# 2. Get shrunken log2FC for NF1KO_BYL vs CTRL_BYL ------------------

coef_name <- "group_NF1KO_BYL_vs_CTRL_BYL"

# Standard DESeq2 results (for reference)
res_nf1ko_byl <- results(dds, name = coef_name)

res_nf1ko_byl_shrunk <- lfcShrink(
    dds,
    coef = coef_name,
    type = "apeglm"
)

# Convert to tibble; rownames are ENTREZ IDs from 01_load_counts.R
res_tbl <- as.data.frame(res_nf1ko_byl_shrunk) %>%
    tibble::rownames_to_column("ENTREZID") %>%
    mutate(ENTREZID = as.character(ENTREZID)) %>%
    filter(!is.na(padj)) %>%
    filter(!is.na(log2FoldChange))

# Save the shrunken results table
write_csv(
    res_tbl,
    "results/tables/NF1KO_BYL_vs_CTRL_BYL_shrunk_results.csv"
)

cat("Shrunken DESeq2 results for NF1KO_BYL vs CTRL_BYL saved.\n")

# 3. Build ranked vector for GSEA ------------------

ranked_tbl <- res_tbl %>%
    arrange(desc(log2FoldChange))

ranked_vector <- ranked_tbl$log2FoldChange
names(ranked_vector) <- ranked_tbl$ENTREZID

# 4. Prepare Hallmark and KEGG gene sets (MSigDB) ------------------

hallmark <- msigdbr(
    species = "Homo sapiens",
    collection = "H" # Hallmark collection
) %>%
    mutate(ncbi_gene = as.character(ncbi_gene)) %>%
    dplyr::select(gs_name, ENTREZID = ncbi_gene)

cat("Loaded Hallmark gene sets.\n")

# Important note: use Reactome instead of KEGG (KEGG removed in msigdbr v10)
reactome_sets <- msigdbr(
    species = "Homo sapiens",
    collection = "C2",
    subcollection = "CP:REACTOME"   # CHANGED: Reactome subcollection instead of CP:KEGG
) %>%
    mutate(ncbi_gene = as.character(ncbi_gene)) %>%  # CHANGED: make IDs character
    dplyr::select(gs_name, ENTREZID = ncbi_gene)     # CHANGED: harmonize to ENTREZID

cat("Loaded Reactome gene sets.\n")

# 5. Run GSEA (Hallmark + Reactome)-------------------------------------

gsea_hallmark <- GSEA(
    ranked_vector,
    TERM2GENE = hallmark,
    pvalueCutoff = 0.05,
    verbose = FALSE
)

# NEW: Reactome GSEA to play the role KEGG would have played
gsea_reactome <- GSEA(
    ranked_vector,
    TERM2GENE = reactome_sets,
    pvalueCutoff = 0.05,
    verbose = FALSE
)

# Convert to data frame and save
gsea_h_tbl <- as.data.frame(gsea_hallmark)
gsea_r_tbl <- as.data.frame(gsea_reactome)

write_csv(
    gsea_h_tbl,
    "results/tables/GSEA_Hallmark_NF1KO_BYL_vs_CTRL_BYL.csv"
)
write_csv(
    gsea_r_tbl,
    "results/tables/GSEA_Reactome_NF1KO_BYL_vs_CTRL_BYL.csv"
)

cat("GSEA Hallmark and Reactome results for NF1KO_BYL vs CTRL_BYL saved.\n")

# 6. Plot top pathways ---------------

# Hallmark dotplot of top 20 enriched pathways
p_dot_hallmark <- dotplot(gsea_hallmark, showCategory = 20) +
    ggtitle("Hallmark GSEA: NF1KO_BYL vs CTRL_BYL")

ggsave(
    "results/figures/GSEA_Hallmark_NF1KO_BYL_vs_CTRL_BYL_dotplot.png",
    p_dot_hallmark,
    width = 8,
    height = 6,
    dpi = 300
)

cat("Hallmark dotplot saved to results/figures/GSEA_Hallmark_NF1KO_BYL_vs_CTRL_BYL_dotplot.png\n")

# KEGG dotplot of top 20 enriched pathways
p_dot_reactome <- dotplot(gsea_reactome, showCategory = 20) +
    ggtitle("Reactome GSEA: NF1KO_BYL vs CTRL_BYL")

ggsave(
    "results/figures/GSEA_Reactome_NF1KO_BYL_vs_CTRL_BYL_dotplot.png",
    p_dot_reactome,
    width = 8,
    height = 6,
    dpi = 300
)

cat("Reactome dotplot saved to results/figures/GSEA_Reactome_NF1KO_BYL_vs_CTRL_BYL_dotplot.png\n")

# 7. Show a specific pathway in detail ------------------

# Check gsea_h_tbl and gsea_k_tbl for available pathways.

# Hallmark: PI3K/AKT/mTOR and KRAS signaling
if ("HALLMARK_PI3K_AKT_MTOR_SIGNALING" %in% gsea_hallmark@result$ID) {
    p_pi3k <- gseaplot2(
        gsea_hallmark,
        geneSetID = "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
        title = "GSEA: PI3K/AKT/mTOR SIGNALING"
    )
    ggsave(
        "results/figures/GSEA_Hallmark_PI3K_AKT_MTOR_SIGNALING.png",
        p_pi3k,
        width = 8,
        height = 6,
        dpi = 300
    )
    cat("Saved Hallmark PI3K_AKT_MTOR_SIGNALING GSEA plot.\n")
}

if ("HALLMARK_KRAS_SIGNALING_UP" %in% gsea_hallmark@result$ID) {
    p_kras <- gseaplot2(
        gsea_hallmark,
        geneSetID = "HALLMARK_KRAS_SIGNALING_UP",
        title = "Hallmark GSEA: KRAS SIGNALING UP"
    )
    ggsave(
        "results/figures/GSEA_Hallmark_KRAS_SIGNALING_UP.png",
        p_kras,
        width = 8,
        height = 6,
        dpi = 300
    )
    cat("Saved Hallmark KRAS_SIGNALING_UP GSEA plot.\n")
}

# Reactome: PI3K-Akt and MAPK signaling pathway

reactome_res <- gsea_reactome@result

# Check if any Reactome pathways related to PI3K exist
if (!is.null(reactome_res) && any(grepl("PI3K", reactome_res$Description, ignore.case = TRUE))) {
    reactome_pi3k_row <- reactome_res %>%
        filter(grepl("PI3K", Description, ignore.case = TRUE)) %>%
        slice(1)  # Take the first match

    reactome_pi3k_id <- reactome_pi3k_row$ID

    p_reactome_pi3k <- gseaplot2(
        gsea_reactome,
        geneSetID = reactome_pi3k_id,
        title = paste0("Reactome GSEA: ", reactome_pi3k_row$Description)
    )

    ggsave(
        "results/figures/GSEA_Reactome_PI3K_Pathway.png",
        p_reactome_pi3k,
        width = 8,
        height = 6,
        dpi = 300
    )
    cat("Saved Reactome PI3K pathway GSEA plot.\n")
}

# Check if any Reactome pathways related to MAPK exist
if (!is.null(reactome_res) && any(grepl("MAPK", reactome_res$Description, ignore.case = TRUE))) {
    reactome_mapk_row <- reactome_res %>%
        filter(grepl("MAPK", Description, ignore.case = TRUE)) %>%
        slice(1)  # Take the first match

    reactome_mapk_id <- reactome_mapk_row$ID

    p_reactome_mapk <- gseaplot2(
        gsea_reactome,
        geneSetID = reactome_mapk_id,
        title = paste0("Reactome GSEA: ", reactome_mapk_row$Description)
    )

    ggsave(
        "results/figures/GSEA_Reactome_MAPK_Pathway.png",
        p_reactome_mapk,
        width = 8,
        height = 6,
        dpi = 300
    )

    cat("Saved Reactome MAPK pathway GSEA plot.\n")
}


# if ("KEGG_PI3K_AKT_SIGNALING_PATHWAY" %in% gsea_kegg@result$ID) {
#     p_kegg_pi3k <- gseaplot2(
#         gsea_kegg,
#         geneSetID = "KEGG_PI3K_AKT_SIGNALING_PATHWAY",
#         title = "GSEA: KEGG PI3K_AKT_SIGNALING_PATHWAY"
#     )
#     ggsave(
#         "results/figures/GSEA_KEGG_PI3K_AKT_SIGNALING_PATHWAY.png",
#         p_kegg_pi3k,
#         width = 8,
#         height = 6,
#         dpi = 300
#     )
#     cat("Saved KEGG PI3K_AKT_SIGNALING_PATHWAY GSEA plot.\n")
# }

# if ("KEGG_MAPK_SIGNALING_PATHWAY" %in% gsea_kegg@result$ID) {
#     p_kegg_mapk <- gseaplot2(
#         gsea_kegg,
#         geneSetID = "KEGG_MAPK_SIGNALING_PATHWAY",
#         title = "GSEA: KEGG MAPK_SIGNALING_PATHWAY"
#     )
#     ggsave(
#         "results/figures/GSEA_KEGG_MAPK_SIGNALING_PATHWAY.png",
#         p_kegg_mapk,
#         width = 8,
#         height = 6,
#         dpi = 300
#     )
#     cat("Saved KEGG MAPK_SIGNALING_PATHWAY GSEA plot.\n")
# }

cat("GSEA analysis for NF1KO_BYL vs CTRL_BYL completed.\n")