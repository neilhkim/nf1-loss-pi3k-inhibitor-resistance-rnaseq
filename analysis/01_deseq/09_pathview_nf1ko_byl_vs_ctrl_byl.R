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
        kegg.native = TRUE,
        kegg.dir = outdir
    )

    # Files are created in `outdir`: hsa<id>.<suffix>.png and .xml
    # Example: hsa04151.NF1KO_BYL_vs_CTRL_BYL.png

    # Pathview naming can vary: with out.suffix it may create both
    # base files (hsa<id>.png/xml) and suffixed ones.
    pattern_png_suffixed <- paste0("hsa", pathway_id, ".*", label, ".*\\.png$")
    pattern_xml_suffixed <- paste0("hsa", pathway_id, ".*", label, ".*\\.xml$")
    files_png <- list.files(outdir, pattern = pattern_png_suffixed, full.names = TRUE)
    files_xml <- list.files(outdir, pattern = pattern_xml_suffixed, full.names = TRUE)
    # Fallback to base filenames if suffixed not found
    if (length(files_png) == 0) {
        files_png <- list.files(outdir, pattern = paste0("^hsa", pathway_id, "\\.png$"), full.names = TRUE)
    }
    if (length(files_xml) == 0) {
        files_xml <- list.files(outdir, pattern = paste0("^hsa", pathway_id, "\\.xml$"), full.names = TRUE)
    }

    if (length(files_png) == 0) {
        warning(paste("No PNG files found for pathway", pathway_id))
        return(invisible(NULL))
    }

    # Rename files inside results/pathview/figures/ to a clearer scheme
    # PNG: KEGG_<pathway_id>_<label>.png
    # XML: KEGG_<pathway_id>_<label>.xml
    for (f in files_png) {
        new_png <- file.path(outdir, paste0("KEGG_", pathway_id, "_", label, ".png"))
        if (basename(f) != basename(new_png)) {
            file.rename(f, new_png)
        }
        message("Saved: ", new_png)
    }
    for (f in files_xml) {
        new_xml <- file.path(outdir, paste0("KEGG_", pathway_id, "_", label, ".xml"))
        if (basename(f) != basename(new_xml)) {
            file.rename(f, new_xml)
        }
        message("Saved: ", new_xml)
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