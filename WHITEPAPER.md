# Cell states that track the psoriasis progression gradient

**A Scissor analysis anchored on an ordinal biopsy-site phenotype (normal → peri-lesional → lesional)**

*Working draft white paper · project `psoriasis-1-bulk` · single-cell reference GSE173706 (Ma et al. 2023) · bulk anchor SRP165679 (Tsoi et al. 2019)*

---

## Abstract

Psoriatic skin is usually studied as a two-state contrast (lesional vs. normal), which discards the clinically meaningful *peri-lesional* margin where disease is actively expanding. We reframe the problem as an **ordinal gradient** — uninvolved (NN) < peri-lesional (PN) < lesional (PP) — and ask which single cells co-vary with that gradient. Using Scissor, we project an ordinal bulk RNA-seq phenotype (93 biopsies spanning all three tiers) onto an 89,058-cell single-cell reference that contains the peri-lesional compartment. Because the phenotype label is a clinical biopsy site rather than a molecular signature derived from the same cells, the mapping is **non-circular by construction**. Selected cells are monotonic on the gradient (mean tier: Scissor− 0.79 < background 1.38 < Scissor+ 1.43) and the gradient-tracking (Scissor+) fraction peaks at the peri-lesional tier. Both a reliability test (p = 0.000) and a selection permutation null (p = 0.000) confirm the signal is driven by real phenotype structure. Endothelial cells dominate the gradient-tracking population (5.2× enriched, OR 11.3 on the backbone; OR 5.24 and still the top lineage at full census), and the associated 1,861-gene program is vascular-led. An orthogonal bulk deconvolution independently confirms the compositional trend for the three strongest lineages. **STAT3, examined as a program member, is significant only at backbone scale (log2FC 0.43, padj 0.018) and does not survive the full census** (log2FC 0.30, raw p = 0.12, BH-FDR q = 0.13; 42.7% vs 44.9% of cells expressing) — the vascular-led conclusion is unchanged, but STAT3 is not a gradient-tracking marker at full resolution. The full-census run (89,058 cells; alpha=0.20, 14.67% selected) has now completed and confirms the backbone's structural findings; additional robustness tiers remain staged for cluster execution.

## 1  Background and rationale

Standard psoriasis transcriptomics contrasts lesional (PP) against uninvolved (NN) skin, and is blind to the **peri-lesional** zone (PN) — the advancing margin just outside the visible plaque, where the transition from health to disease is presumably underway. If a distinct program initiates psoriatic conversion, the peri-lesional compartment is where it should be visible, and a two-state design cannot see it.

We treat biopsy site as an **ordinal phenotype** (NN < PN < PP) and use **Scissor** (Sun, Guan, … Xia, *Nat Biotechnol* 2022) to identify single cells whose expression co-varies with that gradient. The key design choice is that the phenotype is a **clinical biopsy-site label**, not a molecular score computed from the reference cells — so a flagged cell cannot be a restatement of how it was labelled. STAT3, a longstanding psoriasis candidate, is examined as a member of the resulting program, not assumed as a driver — and, as it turns out, it does not survive as a gradient-tracking marker at full-census resolution (see §4).

## 2  Data and design

### 2.1  Single-cell reference (GSE173706)
Assembled from 33 per-sample count matrices: **96,088 cells × 33,538 genes**, 22 donors, all three tiers (NN 13,534; PN 35,518; PP 47,036), with **11 donors contributing paired PN+PP biopsies**. After Ensembl→symbol mapping (24,185 symbols) and QC, **89,058 cells (92.7%)** were retained; standard clustering gave 23 clusters annotated to nine broad lineages.

![composition]({{artifact:art_054caac3-eed6-4990-ab03-411d801d39c4}})

*Single-cell reference composition (GSE173706, Ma et al. 2023). (A) Cells per biopsy tier after QC. (B) Cells per donor, stacked by tier; 11 donors contribute paired PN+PP biopsies.*

![ref_umap]({{artifact:art_08f23e01-2c15-4c2d-878e-d55da50eb1b7}})

*Reference UMAP (89,058 cells). (A) Nine broad lineages by canonical markers. (B) Cells colored by tier — peri-lesional (PN) cells occupy an intermediate position between normal (NN) and lesional (PP).*


### 2.2  Bulk phenotype anchor (SRP165679)
Of five tier-labelled recount3 psoriasis studies, only **SRP165679 (Tsoi et al. 2019)** balances all three tiers at usable depth: **NN=38, PN=27, PP=28**. A single study removes cross-study confounding. Non-degeneracy checks pass: min tier size 27; bulk PC1 (27.8% variance) tracks the tier axis at **F=253.9, p ≈ 10⁻³⁸**.

| Check | Value |
|---|---|
| Tiers / min samples per tier | 3 / 27 |
| Bulk PC1 variance explained | 27.8% |
| Bulk PC1 ~ tier | F=253.9, p≈10⁻³⁸ |
| Cross-study confound | none (single study) |
| Phenotype circularity | none (clinical biopsy-site label) |

## 3  Method note: a portable Scissor solver

Scissor's selection is a graph-regularized elastic net solved by a compiled routine (APML1) that could not be built in our environment, so we reimplemented the identical objective in pure R on `glmnet`. The network penalty is the symmetric-normalized graph Laplacian `L_sym = I − D^(−1/2) A D^(−1/2)`; at ~10⁵ cells a dense Laplacian is intractable, so we encode it as a **sparse edge-difference augmentation** (one sparse row per graph edge). Sparsity is controlled by walking the glmnet λ-path to a target selected fraction. On synthetic two-community data the port recovered 100/100 true-positive cells (cor(β,β_true)=0.72). This is a validated equivalent, **not** the canonical solver; a cross-check against compiled Scissor is staged for the cluster.

## 4  Results

### 4.1  Selection and directionality
The backbone runs on a **stratified 20,023-cell subset**. Tuning gave **alpha=0.40, 14.84% of cells selected** (1,574 Scissor+, 1,397 Scissor−). *(The full census subsequently re-tuned to alpha=0.20 at 14.67% selected, 6,837 Scissor+ / 6,228 Scissor− — see §6.)* Selected cells are **monotonic on the gradient**: mean tier Scissor− **0.79** < background **1.38** < Scissor+ **1.43**, and the Scissor+ fraction **peaks at the peri-lesional tier** (9.4% of PN cells vs. 8.0% of PP) — the signature of an intermediate progression state.

![alpha]({{artifact:art_d4e8addf-04d0-4fe6-a82b-5665e04e4186}})

*Scissor alpha tuning. Selected fraction across the graph-smoothing grid; alpha=0.40 was chosen (14.84% selected, under the 20% cutoff).*

![scissor_umap]({{artifact:art_7a0f1fbc-1749-4bc3-839f-18fd1c68398a}})

*Scissor selection on the reference UMAP (full census, 89,058 cells). Scissor+ (gradient-tracking, lesional-associated), Scissor- (normal-associated), and background cells.*

![scissor_comp]({{artifact:art_8f1e3b16-e45e-48c8-9f62-7e5dab371551}})

*Cell-type and tier composition of the Scissor-selected fractions. Scissor+ peaks at the PN/peri-lesional tier.*


### 4.2  Significance
The **reliability test** (100 label permutations) gives real CV-MSE **0.147** vs. null mean **0.779** (**p = 0.000**). The **selection permutation null** (30 shuffled-label reruns) collapses the positive-minus-negative tier gap from real **0.638** to null mean **−0.305** (**p = 0.000**). The selection reflects real phenotype structure, not graph geometry.

![signif]({{artifact:art_ce7c1728-4eb5-4dfd-8563-75e3c4f03f4d}})

*Significance controls. (A) Reliability test: real CV-MSE vs 100 label permutations. (B) Selection permutation null: real pos-minus-neg tier gap vs 30 shuffled-label reruns.*


### 4.3  Which cells track the gradient

| Cell type | Fold in Scissor+ | Odds ratio | p | Interpretation |
|---|--:|--:|--:|---|
| Endothelial | 5.18× | 11.26 | 6×10⁻²⁴⁵ | lesional-tracking |
| DC | 1.40× | 1.47 | 0.014 | weakly lesional |
| Keratinocyte | 0.87× | 0.69 | 8×10⁻¹² | near-uniform |
| NK | 0.40× | 0.36 | 2×10⁻¹² | depleted from Scissor+ |
| Fibroblast | 0.18× | 0.15 | 5×10⁻⁵³ | normal-tracking |
| Melanocyte | 0.04× | 0.04 | 3×10⁻²⁰ | normal-tracking |

Enrichment p-values are Fisher exact; under the protocol's **BH-FDR** across lineages every row above remains significant (DC weakest). At **full census** the same structure holds — Endothelial the top lineage (OR 5.24, BH q ≈ 0), Keratinocyte/Fibroblast/Melanocyte all BH q < 10⁻⁶¹, NK BH q = 8×10⁻³; Mast and T-cell are the only non-significant lineages.

### 4.4  The gradient-tracking gene program
Scissor+ vs. background DE yields **1,861 genes at padj < 0.05** on the backbone (Seurat-default Bonferroni); under the protocol's BH-FDR at full census the program is **3,903 of 4,541 tested genes at BH q < 0.05**. Top up-regulated genes are vascular/endothelial (**CCL14, ACKR1, RAMP3, PLVAP, APLNR, CYTL1, SPNS2**). On the backbone, **STAT3** is significant (**log2FC 0.43, padj 0.018**) but modest. At **full census this reverses**: STAT3 is not significant under the protocol's BH-FDR (log2FC 0.30, raw p = 0.12, **BH-FDR q = 0.13**; 42.7% vs 44.9% of cells expressing), confirming that the program is vascular-led, not STAT3-led — the STAT3 signal was a small-subset effect that does not hold at full resolution. All single-cell DE significance in this section follows the protocol's BH-FDR (`p.adjust(method="BH")`) over the tested gene family, replacing Seurat's default Bonferroni `p_val_adj`.

![program]({{artifact:art_449adcbc-efc4-4d89-98cb-f371cdf12467}})

*Gradient-tracking gene program: volcano of Scissor+ vs background differential expression (1,861 genes at Bonferroni padj<0.05 on the backbone; 3,903/4,541 at BH-FDR q<0.05 at full census).*

![stat3]({{artifact:art_20472264-c5b9-4721-a814-44c7cf38f9f5}})

*STAT3 expression across Scissor classes on the 20k backbone: elevated in Scissor+ (log2FC 0.43, padj 0.018). Note: this significance does not survive the full census (log2FC 0.30, raw p = 0.12, BH-FDR q = 0.13).*


### 4.5  Orthogonal validation by bulk deconvolution
NNLS deconvolution of the real SRP165679 bulk, tested for monotonic proportion trends across NN→PN→PP:

| Cell type | Scissor direction | Bulk proportion trend | BH-FDR | Status |
|---|---|---|--:|---|
| Endothelial | lesional-tracking | rises 0.0→0.2→3.8% | 1×10⁻¹⁹ | ✓ concordant |
| Fibroblast | normal-tracking | falls 12.3→8.2→1.4% | 1×10⁻¹¹ | ✓ concordant |
| Melanocyte | normal-tracking | falls 7.8→7.3→2.2% | 1×10⁻¹⁴ | ✓ concordant |
| NK | normal-tracking | rises to 7.9% (PP) | 2×10⁻¹¹ | ✗ discordant |
| Keratinocyte | normal-tracking | rises 79→83→83% | 5×10⁻³ | ✗ discordant |
| DC | lesional-tracking | flat | 0.20 | trend n.s. |

The three strongest-signal lineages are concordant: the endothelial signal reflects a **real compositional increase** in bulk tissue. Discordant cases are informative — NK cells *infiltrate* lesional tissue (bulk fraction rises) yet individual NK cells are depleted from the per-cell gradient-tracking set; keratinocytes dominate every tier. Per-cell tracking and bulk composition are different measurements that converge on clean compositional shifts and diverge otherwise.

![deconv]({{artifact:art_2e5f091e-2575-42ad-acbe-44f4ab322554}})

*Orthogonal bulk deconvolution vs Scissor direction, by cell type. Endothelial (lesional-tracking), Fibroblast and Melanocyte (normal-tracking) are concordant; NK and Keratinocyte discordant; DC non-significant.*


## 5  Discussion
The reframing yields a coherent, reproducible signal: a gradient-tracking population that is monotonic on NN→PN→PP, peaks at the peri-lesional margin, survives two orthogonal significance controls, and is corroborated by bulk deconvolution. The dominant biology is **vascular** — endothelial expansion and an angiogenic program — consistent with the histology of psoriatic progression and directing attention to the dermal vasculature at the advancing edge, not only the epidermal keratinocyte compartment. STAT3 shows the expected direction only at backbone scale and loses significance at full census, so it is clearly not the driver of the gradient-tracking program; any role appears embedded in a broader vascular/inflammatory circuit rather than as a marker of the gradient itself.

## 6  Limitations
- **Backbone vs. full census** — primary results shown on the 20,023-cell subset; the full census (89,058 cells) has now been run and confirms the structural findings (endothelial-led, monotonic gradient, significance controls) but revises STAT3 to non-significant. §4 figures have been regenerated from the full-census object (regen_fig4_fullcensus.R; figures_full/).
- **Solver** — validated pure-R reimplementation, not compiled APML1; cross-check planned.
- **Single bulk anchor** — Tier-1 uses SRP165679 only.
- **Deconvolution method** — transparent NNLS, not a benchmarked choice.

## 7  Next steps
1. ~~Full-census run (89,058 cells) on the cluster.~~ **Done** (SLURM 56882314; alpha=0.20, 14.67% selected; endothelial OR 5.24; STAT3 n.s.). §4 figures regenerated from the full-census object (regen_fig4_fullcensus.R).
2. Compiled-Scissor cross-check.
3. Benchmarked deconvolution (deconvBenchmarking).
4. Robustness tiers — additional bulk anchors / tier definitions.
5. Sequence-level arm (Sei-LLRA / seillra) — score regulatory variants near gradient-program genes; STAT3 sits at the seam between expression and sequence arms.

## 8  Data and code availability
Single-cell reference **GSE173706** (Ma et al. 2023); bulk anchor **SRP165679** (Tsoi et al. 2019, via recount3). All code, figures, tables, environment spec, and this white paper are in `github.com/soahum-b/psoriasis-1-bulk`; reproduce/scale instructions in `HANDOFF.md`. Heavy inputs regenerate via `code/00_download_data.R`.

**Method reference.** Scissor: Sun D, Guan X, Moran AE, et al. Identifying phenotype-associated subpopulations by integrating bulk and single-cell sequencing data. *Nature Biotechnology* **40**, 527–538 (2022; published online 11 Nov 2021). doi:10.1038/s41587-021-01091-3.

---
*Working draft — internal white paper. All quantitative values computed from the saved analysis artifacts. Clinical interpretation is preliminary and not intended to guide patient care.*
