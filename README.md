# NF1 Loss Drives PI3Kα Inhibitor Resistance - RNA-seq Reanalysis (GSE207514)

This repo is an **independent reanalysis** of public RNA-seq data from **GSE207514**.

This project:

- **Demonstrates my skills** in RNA-seq analysis (DESeq2, GSEA, Pathview, visualization) on a real oncology dataset.
- **Independently reproduce and validate** a published resistance mechanism.  
- **Provide an enhanced analysis** than the original paper by:
  - Using a clean factor design with two reference levels (CTRL_Veh and CTRL_BYL).
  - Emphasizing KRAS/MAPK vs PI3K signaling using Hallmark + Reactome GSEA and KEGG overlays. {do: write about how the original analysis was so that I can call this, by comparison, an improvement}
  - Inspecting heatmaps of top DE genes under NF1KO + BYL to see how FOXO-like drug response programs change {do: if this is something the original authors did not do, make clear of that.}.

## Background

- PI3K/AKT/mTOR signaling is a key driver in hormone receptor-positive (HR+) / HER2-negative breast cancer.
- Activating mutations in PI3KCA is present in ~40% of cases {do: add ref}.
- Alpelisib (BYL719) is a PI3Kα inhibitor, which improves progression-free survival in patients.
- However, reisistance builds.
- NF1, a tumor supressor and RAS-GAP, when absent, has been identified to enable bypassing of PI3Kα inhibition by Alpelisib.
- {do: Add motivation for study. e.g., It would benefit to look at how XXXX under YYYY}

## Dataset

- RNA-seq data from GEO series GSE207514.
- Profiles T47D ER+ breast cancer cells harboring PI3KCA H1047R mutation.
- Includes 4 groups:
    - CTRL_Veh (control, vehicle)
    - CTRL_BYL (control, alpelisib)
    - NF1_Veh (NF1 knockout, vehicle)
    - **NF1_BYL (NF1 knockout, alpelisib)** 

## Primary Comparisons

I made three key comparisons to understand the transcriptomic effects:

1. Drug effect in control cells
- CTRL_BYL vs CTRL_Veh
- Question: What does alpelisib do in NF1-intact PI3KCA H1047R cells?

2. NF1 loss at baseline
- NF1_Veh vs CTRL_Veh
- Question: How does NF1 loss alone change signaling and transcription?

3. NF1-mediated resistance under drug (*main project goal*) 
- NF1_BYL vs CTRL_BYL
- Question: When treated with alpelisib, how do the drug-resistant NF1-KO cells' gene expression profiles differ from the drug-sensitive CTRL cells? 

## Methods Summary

1. **Load counts and metadata**  
   `01_load_counts.R` builds `counts_matrix.rds`, gene annotations, and `sample_info.rds`.

2. **Build DESeq2 objects with two reference groups**  
   `02_build_dds.R` creates `dds_vehRef` and `dds_bylRef` so that:
   - CTRL_BYL vs CTRL_Veh  
   - NF1KO_Veh vs CTRL_Veh  
   - NF1KO_BYL vs CTRL_BYL  
   are all cleanly available.

3. **PCA**  
   `04_vst_and_pca.R` runs VST and saves a PCA plot showing separation by NF1 status and drug.

4. **Differential expression**  
   `05_`, `06_`, `07_` scripts compute DEGs for the three main contrasts and save full and filtered tables.

5. **Volcano plots**  
   `08_volcano_plots.R` and `07_volcano_detailed.R` generate basic and labeled volcano plots with gene symbols.

6. **GSEA (Hallmark + Reactome)**  
   `09_gsea_nf1ko_byl_vs_ctrl_byl.R` uses **apeglm-shrunken** log2FC, Hallmark and Reactome gene sets, and saves dotplots and GSEA curves (for example KRAS_SIGNALING_UP).

7. **KEGG overlays (Pathview)**  
   `10_pathview_nf1ko_byl_vs_ctrl_byl.R` maps log2FC onto **PI3K-Akt** and **MAPK** pathways.

8. **Heatmap of top DE genes**  
   `10_pheatmap_nf1ko_byl_vs_ctrl_byl.R` plots a VST heatmap of the top 50 DE genes, with samples ordered as:  
   `CTRL_Veh → CTRL_BYL → NF1KO_Veh → NF1KO_BYL`.

Outputs are written to `results/tables/`, `results/figures/`, and `results/pathview/`.


## Key Results {do: clarify which are validation of original authors analysis/claim and which are mine - or add some new value}

From this reanalysis I confirm that:

- **NF1KO cells maintain KRAS/MAPK activity** under BYL719, with persistent cell cycle and DNA replication programs.  
- **CTRL cells** show a strong BYL719-induced FOXO/stress-like response that is **blunted** in NF1KO cells.  
- Pathview overlays highlight:
  - Downregulation of several PI3K inputs (for example EGFR, ERBB2, PIK3CD).
  - Upregulation of HRAS–RAF–ERK components in NF1KO + BYL.

This independently supports the published model that:

> **NF1 loss reduces dependence on PI3Kα and allows escape through Ras-driven parallel signaling.**

## How to run the main pieces

From the project root:

```bash
Rscript analysis/05_deseq2/01_load_counts.R
Rscript analysis/05_deseq2/02_build_dds.R
Rscript analysis/05_deseq2/04_vst_and_pca.R
Rscript analysis/05_deseq2/05_contrast_ctrl_byl_vs_ctrl_veh.R
Rscript analysis/05_deseq2/06_contrast_nf1_veh_vs_ctrl_veh.R
Rscript analysis/05_deseq2/07_contrast_nf1ko_byl_vs_ctrl_byl.R
Rscript analysis/05_deseq2/07_volcano_detailed.R
Rscript analysis/05_deseq2/09_gsea_nf1ko_byl_vs_ctrl_byl.R
Rscript analysis/05_deseq2/10_pathview_nf1ko_byl_vs_ctrl_byl.R
Rscript analysis/05_deseq2/10_pheatmap_nf1ko_byl_vs_ctrl_byl.R
```

## Packages

I use a Conda/Micromamba environment with R 4.4 and the following key packages:

DESeq2, clusterProfiler, msigdbr, enrichplot, pathview, org.Hs.eg.db, tidyverse, pheatmap.
