# Loss-of-Y (transcriptional) — exploratory check in skin single-cell data

**Date:** 2026-07-17 · **Data:** `results_full/reference_scissor_full.rds` (GSE173706, 89,058 cells)
**Method:** Y-chromosome gene *expression* proxy (RNA-seq), NOT genomic mLOY.

## What this measures (and what it does not)
- **Measured:** in male donors, the fraction of cells with **zero detectable Y-gene expression**
  (RPS4Y1, DDX3Y, UTY, EIF1AY, KDM5D, USP9Y). This is an RNA-level readout — a transcriptional
  proxy for loss-of-Y.
- **NOT measured:** genomic mLOY (DNA-level clonal Y-loss vs a diploid baseline). That requires
  SNP array / WGS / genotyping on blood; not available here.
- Sex called per donor from mean Y-gene sum vs XIST: **21 male, 12 female** donors.

## Result — dominated by 10x dropout; no disease signal
Male cells n = 65,524. Overall Y-silent = 10.3%, BUT this collapses with sequencing depth:

| nFeature quartile | Y-silent % |
|---|---|
| Q1 (low)  | 31.5 |
| Q2        | 5.8 |
| Q3        | 2.5 |
| Q4 (high) | **1.4** |

→ 22× drop shallow→deep ⇒ apparent "loss of Y" is overwhelmingly **technical dropout**, not biology.

High-depth (Q4) male cells, Y-silent % by lesion tier: NN 0.6 → PN 3.3 → PP 1.2 — **not monotonic**
with disease; genuine mLOY pathology would predict PP-high. By celltype the "high Y-silent" lineages
(Melanocyte 10.3, T-cell 10.0, Myeloid 8.8) are the lowest-capture populations → capture artifact.

## Conclusion
**No detectable loss-of-Y-expression signal above the dropout floor in skin, and no association with
disease tier.** Consistent with the ruling that mLOY is a hematopoietic phenomenon: skin has only
sparse infiltrating immune cells, and RNA dropout in 10x cannot separate "cell lost Y" from "we
failed to capture Y transcripts." A real test needs blood (deep scRNA-seq or, for genomic mLOY,
DNA/genotyping). Clean negative — theory 2 should not be pursued in skin.
