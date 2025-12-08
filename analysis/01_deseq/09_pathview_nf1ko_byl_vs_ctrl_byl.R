# 09_pathview_nf1ko_byl_vs_ctrl_byl.R
# Generate pathway maps using Pathview for significant pathways from GSEA
# Comparing NF1KO_BYL vs CTRL_BYL

suppressPackageStartupMessages({
    library(tidyverse)
    library(pathview)
})

# 1. Load shrunken log2FC results ------------------
res_tbl <- readr::read_csv(
    "results/gsea/tables/NF1KO_BYL_vs_CTRL_BYL_shrunk_results.csv",
    show_col_types = FALSE
)

# Row names from DESeq2 were ENTREZ IDs, which I stored in ENTREZID column
res_tbl <- res_tbl %>%
    filter(!is.na(log2FoldChange)) %>%
    mutate(ENTREZID = as.character(ENTREZID))

# Build named vector: log2FC with Entrez IDs as names
fc_vec <- res_tbl$log2FoldChange
names(fc_vec) <- res_tbl$ENTREZID

# 2. Output directory ------------------
outdir <- "results/pathview/figures"
if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
}

# Helper function to run pathview and move outputs --------
run_and_save_pathview <- function(fc_vec, pathway_id, label) {
    pv <- pathview(
        gene.data = fc_vec,
        pathway.id = pathway_id,
        species = "hsa",
        out.suffix = label,
        kegg.native = TRUE
    )

    # Two files are created: hsa<id>.<suffix>.png and .xml
    # Example: hsa04151.NF1KO_BYL_vs_CTRL_BYL.png

    pattern <- paste0("hsa", pathway_id, ".*", label, ".*\\.png$")
    files <- list.files(".", pattern = pattern, full.names = TRUE)

    if (length(files) == 0) {
        warning(paste("No PNG files found for pathway", pathway_id))
        return(invisible(NULL))
    }

    # Move files into results/pathview/figures/
    for (f in files) {
        new_name <- file.path(
            outdir,
            paste0("KEGG_", pathway_id, "_", label, ".png")
        )
        file.rename(f, new_name)
        message("Saved: ", new_name)
    }

    invisible(NULL)
}

# 3. PI3K-Akt signaling pathway (hsa04151) ------------------

run_and_save_pathview(
    fc_vec = fc_vec,
    pathway_id = "04151",
    label = "NF1KO_BYL_vs_CTRL_BYL_PI3K_Akt"
)

# 4. MAPK signaling pathway (hsa04010) ------------------

run_and_save_pathview(
    fc_vec = fc_vec,
    pathway_id = "04010",
    label = "NF1KO_BYL_vs_CTRL_BYL_MAPK"
)


cat("Pathview KEGG diagrams written to results/pathview.\n")