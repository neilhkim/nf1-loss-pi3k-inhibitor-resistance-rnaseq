#!/usr/bin/env Rscript
# 08_compare_pipeline_vs_geo_visual.R
# Visual comparison between this-project-generated counts (ENTREZ-based)
# and GEO-downloaded counts (ENTREZ-based), including:
# - Global QC plots (scatter, Bland–Altman, ratio, MA)
# - Condition-level trends for pathway/control genes
#   * One PNG per category
#   * One tall multi-row PNG combining all categories,
#     where each gene row has the same height across categories

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

cat("[compare] Starting targeted visual comparison: this-project-generated vs GEO-downloaded\n")

# -------------------------------------------------------------------
# 0. Helper: ensure GEO count table is present
# -------------------------------------------------------------------

raw_dir <- "data/raw"
if (!dir.exists(raw_dir)) {
  dir.create(raw_dir, recursive = TRUE)
}

geo_txt <- file.path(raw_dir, "GSE207514_RNAseq_T47D_CTRL_NF1KO_annotated_count_table.txt")
geo_gz  <- paste0(geo_txt, ".gz")

geo_urls <- c(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207514/suppl/GSE207514%5FRNAseq%5FT47D%5FCTRL%5FNF1KO%5Fannotated%5Fcount%5Ftable.txt.gz",
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE207514&format=file&file=GSE207514%5FRNAseq%5FT47D%5FCTRL%5FNF1KO%5Fannotated%5Fcount%5Ftable%2Etxt%2Egz"
)

download_geo_if_needed <- function() {
  if (file.exists(geo_txt)) {
    cat("[compare] Found existing GEO txt:\n   ", geo_txt, "\n")
    return(invisible(NULL))
  }

  cat("[compare] GEO txt not found. Attempting download...\n")

  if (!file.exists(geo_gz)) {
    success <- FALSE
    for (u in geo_urls) {
      cat("[compare]  Trying URL: ", u, "\n")
      res <- try(
        download.file(u, destfile = geo_gz, mode = "wb", quiet = TRUE),
        silent = TRUE
      )
      if (!inherits(res, "try-error")) {
        success <- TRUE
        cat("[compare]  Downloaded gzipped GEO table to: ", geo_gz, "\n")
        break
      }
    }
    if (!success) {
      stop("[compare] ERROR: Failed to download GEO count table. Please download manually.")
    }
  } else {
    cat("[compare] Found existing GEO gz:\n   ", geo_gz, "\n")
  }

  cat("[compare] Decompressing GEO table to plain txt...\n")
  geo_df <- read.delim(gzfile(geo_gz), check.names = FALSE)
  write.table(
    geo_df,
    file = geo_txt,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  cat("[compare] Wrote GEO txt table to:\n   ", geo_txt, "\n")
}

download_geo_if_needed()

# -------------------------------------------------------------------
# 1. Load this-project-generated ENTREZ matrix and GEO matrix
# -------------------------------------------------------------------

pipeline_path <- "data/processed/counts_matrix_entrez.rds"
if (!file.exists(pipeline_path)) {
  stop("[compare] ERROR: counts_matrix_entrez.rds not found at ", pipeline_path,
       "\n       Run 07_convert_ensembl_to_entrez.R first.")
}

counts_pipeline <- readRDS(pipeline_path)
cat("[compare] This-project-generated matrix dimensions: ",
    nrow(counts_pipeline), " x ", ncol(counts_pipeline), "\n")

# GEO table (ENTREZID, SYMBOL, GENENAME, count columns)
counts_geo_raw <- read.delim(geo_txt, check.names = FALSE)
cat("[compare] GEO matrix dimensions: ",
    nrow(counts_geo_raw), " x ", ncol(counts_geo_raw), "\n")

# Identify GEO count columns (samples)
geo_count_cols <- grep("^(CTR_|NF1_)", colnames(counts_geo_raw), value = TRUE)

counts_geo_mat <- as.matrix(counts_geo_raw[, geo_count_cols, drop = FALSE])
rownames(counts_geo_mat) <- as.character(counts_geo_raw$ENTREZID)

# -------------------------------------------------------------------
# 2. Align by shared ENTREZ IDs
# -------------------------------------------------------------------

shared_entrez <- intersect(rownames(counts_pipeline), rownames(counts_geo_mat))
cat("[compare] Shared ENTREZ IDs: ", length(shared_entrez), "\n")

pipeline_shared <- counts_pipeline[shared_entrez, , drop = FALSE]
geo_shared      <- counts_geo_mat[shared_entrez, , drop = FALSE]

# Output directory
out_dir <- "results/comparison"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------
# 3. Global comparisons: gene totals, Bland–Altman, ratios, MA plot
# -------------------------------------------------------------------

cat("[compare] Computing global summaries...\n")

pseudo <- 1

pipe_tot <- rowSums(pipeline_shared)
geo_tot  <- rowSums(geo_shared)

global_df <- tibble(
  ENTREZID    = shared_entrez,
  this_proj   = pipe_tot,
  geo         = geo_tot,
  log10_proj  = log10(pipe_tot + pseudo),
  log10_geo   = log10(geo_tot  + pseudo),
  mean_log10  = (log10(pipe_tot + pseudo) + log10(geo_tot + pseudo)) / 2,
  diff_log10  = log10(pipe_tot + pseudo) - log10(geo_tot + pseudo),
  log2_ratio  = log2((pipe_tot + pseudo) / (geo_tot + pseudo)),
  A_log2      = (log2(pipe_tot + pseudo) + log2(geo_tot + pseudo)) / 2,
  M_log2      = log2(pipe_tot + pseudo) - log2(geo_tot + pseudo)
)

# 3a. Scatter
p_scatter <- ggplot(global_df, aes(x = log10_geo, y = log10_proj)) +
  geom_point(alpha = 0.3, size = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    x = "GEO-downloaded gene total (log10)",
    y = "This-project-generated gene total (log10)",
    title = "Gene total comparison\n(this-project-generated vs GEO-downloaded)"
  ) +
  theme_bw()

ggsave(
  filename = file.path(out_dir, "pipeline_vs_geo_gene_totals.png"),
  plot = p_scatter, width = 6, height = 5, dpi = 300
)

# 3b. Bland–Altman (log10)
p_ba <- ggplot(global_df, aes(x = mean_log10, y = diff_log10)) +
  geom_point(alpha = 0.3, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Mean log10(gene total)",
    y = "Difference log10(this-project-generated) - log10(GEO-downloaded)",
    title = "Bland–Altman plot\n(this-project-generated vs GEO-downloaded)"
  ) +
  theme_bw()

ggsave(
  filename = file.path(out_dir, "pipeline_vs_geo_bland_altman.png"),
  plot = p_ba, width = 6, height = 5, dpi = 300
)

# 3c. Ratio distribution
p_ratio <- ggplot(global_df, aes(x = log2_ratio)) +
  geom_histogram(bins = 80) +
  labs(
    x = "log2(this-project-generated / GEO-downloaded)",
    y = "Number of genes",
    title = "Distribution of gene-total ratios\n(this-project-generated / GEO-downloaded)"
  ) +
  theme_bw()

ggsave(
  filename = file.path(out_dir, "pipeline_vs_geo_ratio_distribution.png"),
  plot = p_ratio, width = 6, height = 5, dpi = 300
)

# 3d. MA plot
p_ma <- ggplot(global_df, aes(x = A_log2, y = M_log2)) +
  geom_point(alpha = 0.3, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "A = mean log2(gene total)",
    y = "M = log2(this-project-generated) - log2(GEO-downloaded)",
    title = "MA plot\n(this-project-generated vs GEO-downloaded)"
  ) +
  theme_bw()

ggsave(
  filename = file.path(out_dir, "pipeline_vs_geo_MA_plot.png"),
  plot = p_ma, width = 6, height = 5, dpi = 300
)

# -------------------------------------------------------------------
# 4. Condition-level means for this-project-generated and GEO
# -------------------------------------------------------------------

cat("[compare] Computing condition-level means...\n")

# 4a. This-project-generated condition means (using runs.tsv)
runs <- read.delim("data/metadata/runs.tsv")

runs <- runs %>%
  filter(sample_id %in% colnames(counts_pipeline)) %>%
  mutate(
    condition = case_when(
      group == "CTRL_Veh"  ~ "CTR_Veh",
      group == "CTRL_BYL"  ~ "CTR_BYL",
      group == "NF1KO_Veh" ~ "NF1_Veh",
      group == "NF1KO_BYL" ~ "NF1_BYL",
      TRUE ~ NA_character_
    )
  )

if (any(is.na(runs$condition))) {
  stop("[compare] ERROR: some samples in runs.tsv could not be mapped to conditions.")
}

cond_groups <- split(runs$sample_id, runs$condition)

pipeline_df <- as.data.frame(pipeline_shared)

pipeline_means <- tibble(
  ENTREZID = rownames(pipeline_df),
  CTR_Veh  = rowMeans(pipeline_df[, cond_groups$CTR_Veh,  drop = FALSE]),
  CTR_BYL  = rowMeans(pipeline_df[, cond_groups$CTR_BYL,  drop = FALSE]),
  NF1_Veh  = rowMeans(pipeline_df[, cond_groups$NF1_Veh,  drop = FALSE]),
  NF1_BYL  = rowMeans(pipeline_df[, cond_groups$NF1_BYL,  drop = FALSE])
)

# 4b. GEO condition means
geo_cols <- colnames(geo_shared)

geo_means <- tibble(
  ENTREZID = rownames(geo_shared),
  CTR_Veh  = rowMeans(geo_shared[, grep("^CTR_Veh_", geo_cols),  drop = FALSE]),
  CTR_BYL  = rowMeans(geo_shared[, grep("^CTR_BYL_", geo_cols),  drop = FALSE]),
  NF1_Veh  = rowMeans(geo_shared[, grep("^NF1_Veh_", geo_cols),  drop = FALSE]),
  NF1_BYL  = rowMeans(geo_shared[, grep("^NF1_BYL_", geo_cols),  drop = FALSE])
)

# rename for join
pipeline_cond <- pipeline_means %>%
  rename(
    CTR_Veh_pipeline = CTR_Veh,
    CTR_BYL_pipeline = CTR_BYL,
    NF1_Veh_pipeline = NF1_Veh,
    NF1_BYL_pipeline = NF1_BYL
  )

geo_cond <- geo_means %>%
  rename(
    CTR_Veh_geo = CTR_Veh,
    CTR_BYL_geo = CTR_BYL,
    NF1_Veh_geo = NF1_Veh,
    NF1_BYL_geo = NF1_BYL
  )

merged_cond <- pipeline_cond %>%
  inner_join(geo_cond, by = "ENTREZID")

# -------------------------------------------------------------------
# 5. Add SYMBOL and define genes-of-interest with categories
# -------------------------------------------------------------------

annot_path <- "data/processed/gene_annotations_entrez.rds"
if (!file.exists(annot_path)) {
  stop("[compare] ERROR: gene_annotations_entrez.rds not found at ", annot_path,
       "\n       It should have been created by 07_convert_ensembl_to_entrez.R.")
}

gene_annot <- readRDS(annot_path) %>%
  mutate(
    ENTREZID = as.character(ENTREZID),
    SYMBOL   = as.character(SYMBOL)
  ) %>%
  select(ENTREZID, SYMBOL)

merged_cond <- merged_cond %>%
  mutate(ENTREZID = as.character(ENTREZID)) %>%
  inner_join(gene_annot, by = "ENTREZID") %>%
  relocate(SYMBOL, .after = ENTREZID)

# Genes-of-interest with categories
genes_of_interest <- tribble(
  ~SYMBOL,   ~category,
  # RAS / MAPK core and feedback
  "KRAS",   "RAS_MAPK",
  "HRAS",   "RAS_MAPK",
  "NRAS",   "RAS_MAPK",
  "BRAF",   "RAS_MAPK",
  "RAF1",   "RAS_MAPK",
  "MAP2K1", "RAS_MAPK",
  "MAP2K2", "RAS_MAPK",
  "MAPK1",  "RAS_MAPK",
  "MAPK3",  "RAS_MAPK",
  "DUSP4",  "RAS_MAPK",
  "DUSP6",  "RAS_MAPK",
  "SPRY2",  "RAS_MAPK",
  "ERRFI1", "RAS_MAPK",
  "EGR1",   "RAS_MAPK",
  "EGR2",   "RAS_MAPK",
  "JUN",    "RAS_MAPK",
  "FOS",    "RAS_MAPK",
  "FOSL1",  "RAS_MAPK",
  "FOSB",   "RAS_MAPK",

  # PI3K / AKT / mTOR inputs and FOXO axis
  "PIK3CA", "PI3K_AKT_mTOR",
  "PIK3CD", "PI3K_AKT_mTOR",
  "ERBB2",  "PI3K_AKT_mTOR",
  "ERBB3",  "PI3K_AKT_mTOR",
  "IGF1R",  "PI3K_AKT_mTOR",
  "AKT1",   "PI3K_AKT_mTOR",
  "AKT2",   "PI3K_AKT_mTOR",
  "AKT3",   "PI3K_AKT_mTOR",
  "FOXO1",  "PI3K_AKT_mTOR",
  "FOXO3",  "PI3K_AKT_mTOR",

  # FOXO / stress response genes
  "TXNIP",  "FOXO_response",
  "BNIP3",  "FOXO_response",
  "GADD45A","FOXO_response",
  "GADD45B","FOXO_response",
  "DDIT4",  "FOXO_response",
  "ATF3",   "FOXO_response",
  "SESN1",  "FOXO_response",
  "RHOB",   "FOXO_response",
  "PMAIP1","FOXO_response",
  "TRIB3",  "FOXO_response",

  # MYC / E2F / G2M proliferation programs
  "MYC",    "Proliferation_MYC_E2F_G2M",
  "E2F1",   "Proliferation_MYC_E2F_G2M",
  "E2F2",   "Proliferation_MYC_E2F_G2M",
  "E2F7",   "Proliferation_MYC_E2F_G2M",
  "CCND1",  "Proliferation_MYC_E2F_G2M",
  "CDK2",   "Proliferation_MYC_E2F_G2M",
  "CDK4",   "Proliferation_MYC_E2F_G2M",
  "MKI67",  "Proliferation_MYC_E2F_G2M",
  "TOP2A",  "Proliferation_MYC_E2F_G2M",
  "CDC20",  "Proliferation_MYC_E2F_G2M",
  "CDC25C","Proliferation_MYC_E2F_G2M",
  "BUB1",   "Proliferation_MYC_E2F_G2M",
  "AURKA",  "Proliferation_MYC_E2F_G2M",
  "AURKB",  "Proliferation_MYC_E2F_G2M",
  "PLK1",   "Proliferation_MYC_E2F_G2M",

  # Apoptosis and stress
  "BCL2L11","Apoptosis_Stress",
  "BAX",    "Apoptosis_Stress",
  "CASP3",  "Apoptosis_Stress",
  "CASP7",  "Apoptosis_Stress",
  "ATF4",   "Apoptosis_Stress",
  "CHAC1",  "Apoptosis_Stress",
  "HERPUD1","Apoptosis_Stress",

  # Housekeeping / reference genes
  "ACTB",   "Housekeeping",
  "GAPDH",  "Housekeeping",
  "RPL13A", "Housekeeping",
  "RPLP0",  "Housekeeping",
  "HPRT1",  "Housekeeping",
  "PPIA",   "Housekeeping",
  "TBP",    "Housekeeping",
  "TUBB",   "Housekeeping",
  "UBC",    "Housekeeping"
)

cat("[compare] Genes of interest:\n")
print(genes_of_interest$SYMBOL)

merged_interest <- merged_cond %>%
  inner_join(genes_of_interest, by = "SYMBOL")

cat("[compare] Found ", nrow(merged_interest), " rows for genes-of-interest.\n")

if (nrow(merged_interest) == 0) {
  cat("[compare] WARNING: None of the genes-of-interest were found in the shared ENTREZ set.\n")
  cat("[compare] Done (no trend plots).\n")
  quit(save = "no")
}

# -------------------------------------------------------------------
# 6. Long-format trends for plotting
# -------------------------------------------------------------------

trend_long <- merged_interest %>%
  select(
    SYMBOL, category,
    CTR_Veh_pipeline, CTR_BYL_pipeline, NF1_Veh_pipeline, NF1_BYL_pipeline,
    CTR_Veh_geo,      CTR_BYL_geo,      NF1_Veh_geo,      NF1_BYL_geo
  ) %>%
  pivot_longer(
    cols      = -c(SYMBOL, category),
    names_to  = c("condition", "source"),
    names_pattern = "^(.*)_(pipeline|geo)$",
    values_to = "expression"
  ) %>%
  mutate(
    source = recode(
      source,
      pipeline = "this-project-generated",
      geo      = "GEO-downloaded"
    ),
    condition = factor(
      condition,
      levels = c("CTR_Veh", "CTR_BYL", "NF1_Veh", "NF1_BYL")
    )
  ) %>%
  group_by(SYMBOL, source) %>%
  mutate(expr_rel = expression / mean(expression + 1)) %>%
  ungroup()

trend_long$category <- factor(
  trend_long$category,
  levels = c(
    "RAS_MAPK",
    "PI3K_AKT_mTOR",
    "FOXO_response",
    "Proliferation_MYC_E2F_G2M",
    "Apoptosis_Stress",
    "Housekeeping"
  )
)

# -------------------------------------------------------------------
# 7. One PNG per category (faceted by SYMBOL)
# -------------------------------------------------------------------

categories <- levels(trend_long$category)

for (cat_name in categories) {
  df_cat <- trend_long %>% filter(category == cat_name)

  if (nrow(df_cat) == 0) next

  n_genes <- df_cat %>% distinct(SYMBOL) %>% nrow()
  fig_width  <- max(7, min(0.9 * n_genes, 18))
  fig_height <- 5

  p_cat <- ggplot(
    df_cat,
    aes(
      x     = condition,
      y     = expr_rel,
      group = source,
      color = source
    )
  ) +
    geom_line(position = position_dodge(width = 0.25)) +
    geom_point(
      position = position_dodge(width = 0.25),
      size     = 1.2
    ) +
    facet_wrap(~ SYMBOL, scales = "free_y") +
    labs(
      title = paste0(
        "Relative expression trends in category: ",
        cat_name,
        "\n(this-project-generated vs GEO-downloaded)"
      ),
      x = "Condition",
      y = "Relative expression (per gene, normalized within each dataset)",
      color = "Data source"
    ) +
    scale_color_manual(
      values = c(
        "this-project-generated" = "#1b73ba",
        "GEO-downloaded"         = "#c93f3f"
      )
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text  = element_text(size = 9),
      plot.title  = element_text(size = 12, face = "bold"),
      legend.position = "top",
      legend.title    = element_text(size = 11, face = "bold"),
      legend.text     = element_text(size = 10)
    )

  out_file <- file.path(
    out_dir,
    paste0("pipeline_vs_geo_trend_", cat_name, ".png")
  )

  ggsave(
    filename = out_file,
    plot     = p_cat,
    width    = fig_width,
    height   = fig_height,
    dpi      = 300
  )

  cat("[compare] Saved category plot for", cat_name, "to:\n  ", out_file, "\n")
}

# -------------------------------------------------------------------
# 8. Multi-category tall figure:
#    split genes into rows of 5 per row, each row is one subplot.
#    This gives equal row height across categories.
# -------------------------------------------------------------------

cat("[compare] Building multi-category trend figure with per-row equal heights...\n")

genes_per_row <- 5L
plots_rows <- list()
row_heights <- c()
row_id <- 1L

for (cat_name in categories) {
  df_cat <- trend_long %>% filter(category == cat_name)
  if (nrow(df_cat) == 0) next

  genes_cat <- sort(unique(df_cat$SYMBOL))
  n_genes_cat <- length(genes_cat)

  # Category label row
  label_plot <- ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = cat_name,
      size  = 5,
      fontface = "bold",
      hjust = 0
    ) +
    theme_void() +
    theme(
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    )

  plots_rows[[row_id]] <- label_plot
  row_heights <- c(row_heights, 0.4)   # label row shorter but fixed
  row_id <- row_id + 1L

  # Split genes of this category into chunks of genes_per_row
  # Optional spacer between categories (same height as label row)

  for (start_idx in seq(1, n_genes_cat, by = genes_per_row)) {

    idx_end <- min(start_idx + genes_per_row - 1L, n_genes_cat)
    genes_chunk <- genes_cat[start_idx:idx_end]

    # Pad chunk with dummy genes so row always has exactly 5 facets
    missing_n <- genes_per_row - length(genes_chunk)
    if (missing_n > 0) {
        dummy_genes <- paste0("EMPTY_", seq_len(missing_n))
        genes_chunk <- c(genes_chunk, dummy_genes)
    }

    # Build data frame including dummy rows with NA expression
    df_row <- df_cat %>%
        filter(SYMBOL %in% genes_chunk) %>%
        mutate(SYMBOL = factor(SYMBOL, levels = genes_chunk))

    # Add dummy NA rows if needed
    if (missing_n > 0) {
        df_dummy <- expand.grid(
        SYMBOL    = dummy_genes,
        condition = levels(df_row$condition),
        source    = c("this-project-generated", "GEO-downloaded"),
        expr_rel  = NA_real_
        )
        df_row <- bind_rows(df_row, df_dummy)
    }

    p_row <- ggplot(
        df_row,
        aes(
            x     = condition,
            y     = expr_rel,
            group = source,
            color = source
        )
    ) +
    geom_line(position = position_dodge(width = 0.25), na.rm = TRUE) +
    geom_point(position = position_dodge(width = 0.25), size = 1.2, na.rm = TRUE) +
    facet_wrap(
        ~ SYMBOL,
        nrow   = 1,
        scales = "free_y"
    ) +
    labs(
        x = "Condition",
        y = "Relative expression"
    ) +
    scale_color_manual(
        values = c(
        "this-project-generated" = "#1b73ba",
        "GEO-downloaded"         = "#c93f3f"
        )
    ) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text  = element_text(size = 9),
        legend.position = "none",
        plot.margin = margin(t = 2, r = 5, b = 5, l = 5)
    )

    plots_rows[[row_id]] <- p_row
    row_heights <- c(row_heights, 1)
    row_id <- row_id + 1L
    }

    spacer <- plot_spacer() +
        theme(plot.margin = margin(t = 0, r = 5, b = 8, l = 5))
    plots_rows[[row_id]] <- spacer
    row_heights <- c(row_heights, 0.3)
    row_id <- row_id + 1L
}

# Remove trailing spacer if present
if (length(plots_rows) > 0 && inherits(plots_rows[[length(plots_rows)]], "patchwork")) {
  # not strictly reliable, keep as is
}

# Combine all rows
combined_trends <- wrap_plots(
  plots_rows,
  ncol    = 1,
  heights = row_heights,
  guides  = "collect"
) +
  plot_annotation(
    title = "Pathway-wide expression trends\nthis-project-generated vs GEO-downloaded",
    subtitle = "Each category has its own block; genes are split into rows of 5 per row"
  )

combined_trends <- combined_trends &
  theme(
    legend.position = "top",
    legend.title    = element_text(size = 11, face = "bold"),
    legend.text     = element_text(size = 10)
  )

out_big <- file.path(
  out_dir,
  "pipeline_vs_geo_trends_all_categories.png"
)

# Height scales with number of row plots so each row stays tall
total_rows <- length(plots_rows)
ggsave(
  filename = out_big,
  plot     = combined_trends,
  width    = 16,
  height   = max(4, 1.3 * total_rows),
  dpi      = 300
)

cat("[compare] Saved multi-category trend figure to:\n  ", out_big, "\n")
cat("[compare] All comparison plots written to:\n  ", out_dir, "\n")
cat("[compare] Done.\n")
