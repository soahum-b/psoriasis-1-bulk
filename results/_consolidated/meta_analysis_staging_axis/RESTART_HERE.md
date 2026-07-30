# RESTART POINT — pathway analysis

Everything up to (but not including) pathway analysis is complete and checkpointed.
To resume this evening, reload the checkpoints below — no need to re-download or re-run DE.

## State as of this checkpoint
- **Study:** SRP035988, recount3, Gencode v26 / GRCh38
- **DE done:** limma-voom, NN (83) vs PP (95), 24,528 genes tested
- **Significant:** 3,477 genes (adj.P<0.05 & |log2FC|>1) — 1,511 up, 1,966 down
- **STAT3:** logFC +1.16, adj.P 5.85e-52 — up, with IL-17/SOCS3 neighbourhood

## Checkpoint files (in this directory)
| File | What it is | Reload with |
|------|-----------|-------------|
| `dge_filt_norm.rds`   | filtered + TMM DGEList (24,528 x 178) | `dge <- readRDS("dge_filt_norm.rds")` |
| `de_results_full.rds` | full DE table (BH + Storey q) | `de <- readRDS("de_results_full.rds")` |
| `gsea_ranks.rds`      | **named t-stat ranking vector, 24,528 genes, pre-sorted** | `ranks <- readRDS("gsea_ranks.rds")` |

## Environment
- conda env **psoriasis-r**: R 4.4.3, edgeR 4.4.0, limma 3.62.1, qvalue 2.38.0
- Full setup + package list in `README.md`; full pipeline in `psoriasis_pipeline.R`

## Next step: PATHWAY ANALYSIS (planned)
1. **GSEA (primary)** on `gsea_ranks.rds` (all genes ranked by limma t) using `fgsea`.
   Gene-set collections: **GO Biological Process, KEGG, Reactome, MSigDB Hallmark**.
   Focus pathways: IL-17 signalling, JAK-STAT, NF-kB, keratinocyte differentiation.
2. **ORA (cross-check)** on the 3,477-gene hit list (`de_results_significant.csv`) via
   `clusterProfiler::enrichGO` / `enrichKEGG`.
3. Packages to install at restart (not yet installed):
   `bioconductor-fgsea`, `bioconductor-clusterprofiler`, `bioconductor-msigdbr` (or `r-msigdbr`),
   `bioconductor-org.hs.eg.db`, `bioconductor-reactomepa`.
4. First figure target: GSEA enrichment plot for the top IL-17 / STAT3 pathway + a dot-plot of
   the top enriched sets.

## To restart cleanly
```r
setwd("/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3")
dge   <- readRDS("dge_filt_norm.rds")
de    <- readRDS("de_results_full.rds")
ranks <- readRDS("gsea_ranks.rds")
```
