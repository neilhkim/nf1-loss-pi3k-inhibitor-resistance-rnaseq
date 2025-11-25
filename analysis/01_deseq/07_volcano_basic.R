# 08_volcano_plots.R
# Make volcano plots for the three main contrasts

library(tidyverse)

# Create output folder if needed
if (!dir.exists("results/tables")) {
    dir.create("results/tables", recursive = TRUE)
}

# Helper function to build one volcano plot -----

make_volcano <- function(infile, contrast_label, outfile) {
    # infile: path to CSV with full DESeq2 results
    # contrast_label: nice label for the title
    # outfile: path to save PNG

    # 1. Read the table
    df <- read_csv(infile, show_col_types = FALSE)

    # 2. Keep rows with non-missing padj and finite values
    df <- df %>%
        filter(!is.na(padj)) %>%
        filter(is.finite(log2FoldChange), is.finite(pvalue))

    # 3. Define significance
    df <- df %>%
        mutate(
            sig = padj < 0.05 & abs(log2FoldChange) >= 1
        )


    # 4. Basic volcano plot
    p <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj))) + 
        geom_point(aes(alpha = sig), size = 1) +
        scale_alpha_manual(values = c(`FALSE` = 0.2, `TRUE` = 0.9)) +
        labs(
            title = paste0("Volcano Plot: ", contrast_label),
            x = "Log2 Fold Change",
            y = "-Log10 Adjusted P-value"
        ) + 
        theme_minimal()

    # 5. Save plot
    ggsave(
        filename = outfile,
        plot = p,
        width = 6,
        height = 5,
        dpi = 300
    )

    cat("Volcano plot saved to:", outfile, "\n")
}

# -----------------------------------------------------------------------
# 1) Drug effect in NF1 intact cells: CTRL_BYL vs CTRL_Veh
# -----------------------------------------------------------------------

make_volcano(
    infile = "results/tables/DEG_CTRL_BYL_vs_CTRL_Veh.csv",
    contrast_label = "CTRL_BYL vs CTRL_Veh",
    outfile = "results/figures/volcano_CTRL_BYL_vs_CTRL_Veh.png"
)

# -----------------------------------------------------------------------
# 2) NF1KO effect at baseline: NF1KO_Veh vs CTRL_Veh
# -----------------------------------------------------------------------

make_volcano(
    infile = "results/tables/DEG_NF1KO_Veh_vs_CTRL_Veh.csv",
    contrast_label = "NF1KO_Veh vs CTRL_Veh",
    outfile = "results/figures/volcano_NF1KO_Veh_vs_CTRL_Veh.png"
)

# -----------------------------------------------------------------------
# 3) NF1KO mediated resistance to BYL: NF1KO_BYL vs CTRL_BYL
# -----------------------------------------------------------------------

make_volcano(
    infile = "results/tables/DEG_NF1KO_BYL_vs_CTRL_BYL.csv",
    contrast_label = "NF1KO_BYL vs CTRL_BYL",
    outfile = "results/figures/volcano_NF1KO_BYL_vs_CTRL_BYL.png"
)

cat("All volcano plots generated.\n")