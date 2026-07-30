# Psoriasis meta-analysis — reproducible pipeline

Differential expression of **lesional psoriatic skin (PP) vs normal skin (NN)**,
built toward a STAT3-centred meta-analysis. This directory contains everything a
reviewer needs to regenerate our results from the raw public data, and everything
*we* need to rebuild our checkpoints after a disconnect.

- **Study:** SRP035988 (Li et al.), accessed via **recount3**
- **Annotation:** Gencode v26 (recount3 "G026"), genome build GRCh38/hg38
- **Design:** NN (n=83) vs PP (n=95); limma-voom moderated t-test
- **Headline result:** 3,477 genes at BH-FDR < 0.05 **and** |log2FC| > 1
  (1,511 up, 1,966 down). STAT3 up +1.16 log2FC (~2.2x) with its IL-17 / SOCS3
  neighbourhood coherently up.

---

## 1. What's in this directory

| File | Description |
|------|-------------|
| `psoriasis_pipeline.R`          | **Master script** — raw download → DE → figures, top to bottom |
| `dge_filt_norm.rds`             | Checkpoint 1: filtered + TMM-normalised DGEList (24,528 × 178) |
| `de_results_full.rds` / `.csv`  | Checkpoint 2: all 24,528 genes (logFC, p, adj.P.Val, qvalue, B) |
| `de_results_significant.csv`    | The 3,477 significant genes |
| `fig6_volcano.png`              | Volcano plot |
| `fig7_qvalue.png`               | p-value histogram + Storey q vs BH diagnostic |
| `whitepaper/`                   | White paper + method deep-dives (normalisation, DE) |

## 2. Software environment

Analysis runs in **R (>= 4.4)** with Bioconductor. We use a conda environment
(`psoriasis-r`) so the toolchain is pinned and portable.

### Option A — recreate the conda environment (recommended)
```bash
conda create -n psoriasis-r -c conda-forge -c bioconda \
  r-base=4.4 \
  bioconductor-recount3 bioconductor-edger bioconductor-limma \
  bioconductor-deseq2 bioconductor-qvalue \
  r-tidyverse r-ggrepel r-pheatmap r-patchwork
conda activate psoriasis-r
```

### Option B — install from within R
```r
install.packages("BiocManager")
BiocManager::install(c("recount3","edgeR","limma","DESeq2","qvalue"))
install.packages(c("dplyr","ggplot2","ggrepel","patchwork"))
```

**Core package versions used** (from `sessionInfo()`): recount3, edgeR, limma,
qvalue, ggplot2, ggrepel. Exact versions are recorded in the artifact environment
snapshot saved alongside the results.

## 3. How to run

From this directory, with the environment active:
```bash
Rscript psoriasis_pipeline.R
```
The script prints progress for each STEP and writes checkpoints as it goes.
Expected wall-clock: a few minutes (the recount3 download in STEP 1 dominates).

To run **interactively** (e.g. to resume from a checkpoint), open the script and
execute step by step. Because each STEP reloads its input from the previous
checkpoint (`readRDS(...)`), you can start at STEP 3 if `dge_filt_norm.rds`
already exists — no need to re-download.

## 4. How to verify the headline numbers

After a run, these should hold:
```r
de <- readRDS("de_results_full.rds")
nrow(de)                                        # 24528 genes tested
sum(de$adj.P.Val < 0.05 & abs(de$logFC) > 1)    # 3477 significant
de[de$gene == "STAT3", ]                         # logFC ~ +1.16, adj.P ~ 5.9e-52
```
The volcano (`fig6_volcano.png`) and the q-value diagnostic (`fig7_qvalue.png`)
should match the figures embedded in the white paper.

## 5. Multiple-testing note (BH vs Storey)

`adj.P.Val` is the **Benjamini-Hochberg** FDR (limma default). `qvalue` is the
**Storey** q-value, which additionally estimates pi0 — the proportion of truly
null genes — from the data (here **pi0 = 0.133**). Because q = pi0 x BH, the
q-values are uniformly smaller (more powerful), but after the |log2FC| > 1
filter both methods converge (3,477 vs 3,479), showing the hit list is robust to
the multiple-testing method. See `whitepaper/differential_expression_deep_dive.md`.

## 6. Provenance / caveats

- recount3 counts are base-pair coverage counts transformed with
  `transform_counts()`; gene rows are collapsed from versioned Ensembl IDs to
  HGNC symbols by summation.
- Group labels are parsed from `sra.sample_attributes` ("normal_skin" → NN,
  "Psoriasis_skin" → PP).
- This is one study; the multi-study meta-analysis (Fisher/Stouffer p-combination
  or inverse-variance effect-size pooling) is a later stage.
