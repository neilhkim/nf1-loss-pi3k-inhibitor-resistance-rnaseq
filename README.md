# NF1 Loss Drives PI3Kα Inhibitor Resistance - RNA-seq Reanalysis (GSE207514)

This repo is an **independent reanalysis** of public RNA-seq data from **GSE207514**.
All wet-lab experiments were performed by the original study.

This project aims to:

- **Demonstrates my skills** in RNA-seq analysis (DESeq2, GSEA, Pathview, visualization) on a real oncology signaling dataset.
- **Reproduce and validate the authors' main conclusions** about NF1-loss-mediated resistance to the PI3Kα inhibitor alpelisib.
- **Provide an enhanced analysis** compared to the original study by:
  - Emphasizing KRAS/MAPK vs PI3K signaling using Hallmark + Reactome GSEA and KEGG overlays. 
    - *The original paper analyzed Hallmark and KEGG pathways, with only limited use of Reactome. In this reanalysis, I use the current msigdbr gene sets, which provide Hallmark and Reactome but not KEGG. So I perform systematic enrichment with Hallmark + Reactome, and include KEGG through Pathview overlays.*
  - Inspecting heatmaps of top DE genes to visualize FOXO-like drug response programs under Alpelisib (BYL719; a PI3Kα drug).
    - *This level of ranked-gene expression patterning is not shown in the original figures and adds a clear view of condition-specific expression shifts.*

## Background

- The PI3K/AKT/mTOR signaling axis is a key driver in hormone receptor-positive (HR+) / HER2-negative breast cancer.
- PI3KCA mutations is present in ~40% of cases.
- Alpelisib (BYL719; a PI3Kα inhibitor) improves progression-free survival in PIK3CA-mutant patients.
- However, reisistance builds.
- NF1, a tumor supressor and RAS-GAP, when absent, has been identified to enable bypassing of PI3Kα inhibition by Alpelisib.
- By reanalyzing this dataset, I hoped to indenpendently examine how NF1 loss reshapes drug response at the transcriptome level. In particular I hoped to examine:
    - which findings are reproducible under a slightly different analytical approach, and
    - how broad PI3K and MAPK pathway activities reorganize under NF1 loss.

## Dataset

- RNA-seq data from GEO series GSE207514.
- Profiles T47D ER+ breast cancer cells harboring PI3KCA H1047R mutation.
- Includes 4 groups:
    - CTRL_Veh (control, vehicle)
    - CTRL_BYL (control, alpelisib)
    - NF1_Veh (NF1 knockout, vehicle)
    - **NF1_BYL (NF1 knockout, alpelisib)** 

- This design allows separation of 
    - drug effects,
    - NF1 knockout effect, and
    - NF1-driven drug resistance

## Primary Comparisons

1. Drug effect in control cells
- CTRL_BYL vs CTRL_Veh
- Question: What does alpelisib do in NF1-intact PI3KCA H1047R cells?

2. NF1 loss at baseline
- NF1_Veh vs CTRL_Veh
- Question: How does NF1 loss alone change signaling and transcription?

3. NF1-mediated resistance under drug (*main project goal*) 
- NF1_BYL vs CTRL_BYL
- Question: When treated with alpelisib, how do the drug-resistant NF1-KO cells' gene expression profiles differ from the drug-sensitive CTRL cells? 

## Methods

The analysis is organized into two stages:

1. **`analysis/00_counts_pipeline/`** – scripts (to be added) that will regenerate the raw count matrix from FASTQ files.
2. **`analysis/01_deseq/`** – the full downstream DESeq2-based transcriptomic analysis used in this reanalysis.

---

### `analysis/00_counts_pipeline/` *(to be added)*

This folder will contain a complete count-generation pipeline that produces a raw gene-by-sample count matrix from FASTQ files.  
The intended workflow is:

1. **Quality control**
   - Run `fastqc` on FASTQs.

2. **Optional trimming**
   - Use tools such as `trimmomatic` if adapter or quality trimming is needed.

3. **Alignment with STAR**
   - Build a STAR genome index.
   - Align sequences with `STAR --runThreadN ...`.
   - Generate sorted BAM files per sample.

4. **Quantification with featureCounts**
   - Run `featureCounts -T ...` to assign aligned reads to genes.
   - Produce a unified raw count matrix.

5. **Save outputs**
   - `counts_matrix.raw.rds`
   - `gene_annotation.rds`

Although this reanalysis uses the count table provided by GSE207514, adding this pipeline will make the repository fully reproducible from raw FASTQ files.

---

### `analysis/01_deseq/` – Downstream RNA-seq analysis

This folder contains all steps from count loading to differential expression, visualization, and pathway interpretation.

#### 1. Load counts and metadata  
**`01_load_counts.R`**
- Reads the GEO count table.
- Separates annotation from counts.
- Builds `sample_info` (genotype, treatment, group).
- Saves:  
  - `counts_matrix.rds`  
  - `gene_annotations.rds`  
  - `sample_info.rds`

#### 2. Build DESeq2 objects and run DESeq  
**`02_build_dds_and_deseq.R`**
- Constructs two DESeq2 datasets:  
  - `dds_vehRef` (reference level: `CTRL_Veh`)  
  - `dds_bylRef` (reference level: `CTRL_BYL`)
- Filters low-count genes.
- Runs `DESeq()`.
- Saves both fitted objects.

#### 3. Variance stabilizing transform and PCA  
**`03_vst_pca.R`**
- Applies VST to `dds_vehRef`.
- Generates PCA plot colored by group.
- Saves `vsd.rds` and PCA figure.

#### 4. Differential expression  
- **`04_deg_ctrl_byl_vs_ctrl_veh.R`**
- **`05_deg_nf1ko_veh_vs_ctrl_veh.R`**
- **`06_deg_nf1ko_byl_vs_ctrl_byl.R`**

Each script:
- Extracts the appropriate contrast from DESeq2.
- Writes full and filtered DEG tables.
- Prints counts of up- and down-regulated genes.

#### 5. Volcano plots  
- **`07_volcano_basic.R`** – simple volcano plots for all contrasts.  
- **`07_volcano_detailed.R`** – labeled volcano plots with shared axes.

#### 6. GSEA (Hallmark + Reactome)  
**`08_gsea_nf1ko_byl_vs_ctrl_byl.R`**
- Uses **apeglm-shrunken** log2FC values.
- Runs GSEA with Hallmark and Reactome gene sets (msigdbr).
- Saves enrichment tables, dotplots, and selected GSEA curves.

#### 7. KEGG pathway overlays  
**`09_pathview_nf1ko_byl_vs_ctrl_byl.R`**
- Maps log2FC values onto KEGG PI3K-Akt and MAPK pathways using Pathview.
- Produces PNG pathway diagrams in `results/pathview/`.

#### 8. Heatmap of top DE genes  
**`10_pheatmap_nf1ko_byl_vs_ctrl_byl.R`**
- Selects top DE genes from the NF1KO_BYL vs CTRL_BYL comparison.
- Extracts VST expression values and performs row-scaling.
- Orders samples as: `CTRL_Veh → CTRL_BYL → NF1_Veh → NF1_BYL`.
- Saves heatmap (PNG + PDF).

## Key Results 

### Validation of the original analysis

- **NF1KO cells maintain KRAS/MAPK activity** even under PI3Kalpha inhibition (BYL).
- **CTRL cells** show a strong BYL719-induced FOXO/stress-like response that is **blunted** in NF1KO cells.  
- NF1KO cells fail to activate many of PI3K-inhibitor-response genes, matching the resistance phenotype described by the authors.

### Additional findings from this reanalysis

- Systematic Reactome enrichment (not emphasized in the original study) reveals more detailed pathway-level distinctions, such as persistent E2F/MYC and G2M checkpoint signals in NF1KO + BYL.
- The heatmap of the top 50 DE genes highlights a sharp separation between:
    - BYL-induced FOXO-like genes rising only in CTRL, and
    - MAPK-feedback genes staying high in NF1KO under drug.
    - This expression-pattern view is complementary to the original study’s pathway summaries.
- Pathview overlays provide node-level pathway regulations:
  - Downregulation of several PI3K inputs (e.g., ERBB2, PIK3CD).
  - Upregulation of HRAS–RAF–ERK components in NF1KO + BYL.

Together, these confirm and extend the interpretation that:
> NF1 loss reduces dependence on PI3Kα and enables escape through Ras-driven parallel signaling.

## How to run the main pieces

### Create env from environment.yml

```bash
micromamba env create -f environment.yml
micromamba activate nf1-rnaseq
```

### Run scripts

```bash
Rscript analysis/01_deseq/01_load_counts.R
Rscript analysis/01_deseq/02_build_dds_and_deseq.R
Rscript analysis/01_deseq/03_vst_pca.R
Rscript analysis/01_deseq/04_deg_ctrl_byl_vs_ctrl_veh.R
Rscript analysis/01_deseq/05_deg_nf1ko_veh_vs_ctrl_veh.R
Rscript analysis/01_deseq/06_deg_nf1ko_byl_vs_ctrl_byl.R
Rscript analysis/01_deseq/07_volcano_basic.R
Rscript analysis/01_deseq/07_volcano_detailed.R
Rscript analysis/01_deseq/08_gsea_nf1ko_byl_vs_ctrl_byl.R
Rscript analysis/01_deseq/09_pathview_nf1ko_byl_vs_ctrl_byl.R
Rscript analysis/01_deseq/10_pheatmap_nf1ko_byl_vs_ctrl_byl.R
```

## Packages

Conda/Micromamba environment with R 4.4.
Key packages: DESeq2, clusterProfiler, msigdbr, enrichplot, pathview, org.Hs.eg.db, tidyverse, pheatmap.