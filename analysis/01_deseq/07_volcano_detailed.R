# 07_volcano_detailed.R
# Improved volcano plots with gene labels and better aesthetics

suppressPackageStartupMessages({
    library(tidyverse)
    library(ggrepel)
})

# Prepare output directories
volcano_fig_dir <- "results/volcano/figures"
dir.create(volcano_fig_dir, showWarnings = FALSE, recursive = TRUE)

# Load gene annotation (ENTREZID -> SYMBOL)
gene_annot <- readRDS("data/processed/gene_annotations_entrez.rds") %>%
    mutate(ENTREZID = as.character(ENTREZID)) %>%
    dplyr::select(ENTREZID, SYMBOL)

# Helper: make a nice volcano plot -------------
make_volcano_labeled <- function(infile, 
                            contrast_label, 
                            outfile_png,
                            top_n = 15,
                            xlim_range = NULL, # e.g., c(-6, 6)
                            ylim_max = NULL # e.g., 250
                            ) {

    df <- read_csv(infile, show_col_types = FALSE) %>%
        filter(!is.na(padj), 
        is.finite(padj),
        is.finite(log2FoldChange))

    # gene column holds Entrez IDs from DESeq2 results
    df <- df %>%
        mutate(gene = as.character(gene)) %>%
        left_join(gene_annot, by = c("gene" = "ENTREZID")) %>%
        mutate(
            # use SYMBOL when available, else fallback to Entrez ID
            label = dplyr::coalesce(SYMBOL, gene),
            # sig = padj < 0.05 & abs(log2FoldChange) >= 1,
            direction = case_when(
                padj < 0.05 & log2FoldChange > 1 ~ "Up",
                padj < 0.05 & log2FoldChange < -1 ~ "Down",
                TRUE ~ "NS"
            ),
            neglog10padj = -log10(padj)
        )


    # Choose label candidates among significant genes with largest |LFC|
    label_df <- df %>%
        filter(direction != "NS") %>%
        arrange(desc(abs(log2FoldChange))) %>%
        slice_head(n = top_n)

    # Default axis limits if not provided
    if (is.null(xlim_range)) {
        max_abs_fc <- max(abs(df$log2FoldChange), na.rm = TRUE)
        xlim_range <- c(-max_abs_fc, max_abs_fc)
    }

    # # Top genes to label
    # top_up <- df %>% 
    #     filter(direction == "Up") %>%
    #     arrange(padj) %>%
    #     slice_head(n = 10)
    # top_down <- df %>%
    #     filter(direction == "Down") %>%
    #     arrange(padj) %>%
    #     slice_head(n = 10)

    # genes_to_label <- bind_rows(top_up, top_down)

    p <- ggplot(df, aes(x = log2FoldChange, y = neglog10padj)) +
        geom_point(aes(color = direction), alpha = 0.6, size = 1.2) +
        scale_color_manual(
            values = c( 
                "Up" = "firebrick",
                "Down" = "steelblue",
                "NS" = "grey70"
                )
        ) + 
        geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "grey50") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
        geom_text_repel(
            data = label_df,
            aes(label = label),
            size = 3,
            max.overlaps = Inf,
            box.padding = 0.5,
            point.padding = 0.2,
            force = 3,
            # direction = "y",
            force_pull = 0.2,
            min.segment.length = 0, # Always draw a segment
            segment.color = "grey50",
            segment.size = 0.3
        ) + 
        coord_cartesian(xlim = xlim_range, ylim = ylim_max) +
        labs(
            title = paste0("Volcano Plot (labeled): ", contrast_label),
            x = "log2 Fold Change",
            y = "-log10 adjusted P-value",
            color = "direction"
        ) +
        theme_minimal(base_size = 12)

    ggsave(outfile_png, p, width = 7, height = 6, dpi = 300)
    cat("Saved volcano plot to: ", outfile_png, "\n")
}

# -------------------------------------------
# Make labeled volcano for each contrast

# First, compute global axis range so all three plots share same scale

get_global_ranges <- function(files) {
    all_df <- purrr::map_df(files, ~ read_csv(.x, show_col_types = FALSE)) %>%
        # ensure numeric
        mutate(
            log2FoldChange = as.numeric(log2FoldChange),
            padj = as.numeric(padj)
        ) %>%
        filter(!is.na(padj), 
                is.finite(padj),
                is.finite(log2FoldChange)) %>%
        mutate(
            padj = ifelse(padj == 0, 1e-320, padj), # avoid -Inf
            neglog10padj = -log10(padj)
            )

    x_max <- max(abs(all_df$log2FoldChange), na.rm = TRUE)
    y_max <- max(all_df$neglog10padj, na.rm = TRUE) * 1.2 # add 5% headroom

    cat("Global x-axis limits: ", -x_max, " to ", x_max, "\n")
    cat("Global y-axis max: ", y_max, "\n")

    list(xlim = c(-x_max, x_max), ylim = c(0, y_max))
}

# Files for all contrasts
files <- c(
    "results/deg/tables/DEG_CTRL_BYL_vs_CTRL_Veh.csv",
    "results/deg/tables/DEG_NF1KO_Veh_vs_CTRL_Veh.csv",
    "results/deg/tables/DEG_NF1KO_BYL_vs_CTRL_BYL.csv"
)

ranges <- get_global_ranges(files)

make_volcano_labeled(
    infile = "results/deg/tables/DEG_CTRL_BYL_vs_CTRL_Veh.csv",
    contrast_label = "CTRL_BYL vs CTRL_Veh",
    outfile_png = file.path(volcano_fig_dir, "volcano_CTRL_BYL_vs_CTRL_Veh_labeled.png"),
    top_n = 12,
    xlim_range = ranges$xlim,
    ylim_max = ranges$ylim
)

make_volcano_labeled(
    infile = "results/deg/tables/DEG_NF1KO_Veh_vs_CTRL_Veh.csv",
    contrast_label = "NF1KO_Veh vs CTRL_Veh",
    outfile_png = file.path(volcano_fig_dir, "volcano_NF1KO_Veh_vs_CTRL_Veh_labeled.png"),
    top_n = 12,
    xlim_range = ranges$xlim,
    ylim_max = ranges$ylim
)

make_volcano_labeled(
    infile = "results/deg/tables/DEG_NF1KO_BYL_vs_CTRL_BYL.csv",
    contrast_label = "NF1KO_BYL vs CTRL_BYL",
    outfile_png = file.path(volcano_fig_dir, "volcano_NF1KO_BYL_vs_CTRL_BYL_labeled.png"),
    top_n = 12,
    xlim_range = ranges$xlim,
    ylim_max = ranges$ylim
)

message("All detailed volcano plots generated.")