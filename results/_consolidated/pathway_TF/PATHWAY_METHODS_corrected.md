# Pathway / enrichment analysis — corrected methodology & before/after

*Companion to METHODS.md (which documented the DE/meta-analysis backbone but not the
pathway tests). Scope: PP-vs-NN lesional-psoriasis meta-analysis, k=4. This note records
what the pathway layer did before, what was wrong with it, what it does now, and the
resulting change in conclusions.*

---

## 1. Summary

The gene-level DE and random-effects meta-analysis were rigorous. The **enrichment layer
built on top of them was not**, and in two places was statistically invalid. This note
documents a corrected re-run. The corrected results **do not overturn the biology** —
the same programmes come up — but they replace inflated / saturated statistics with valid
ones, and they reveal an honest ordering: the most *reproducible* signals across studies
are proliferation, RNA-processing/translation, and interferon; the IL6-JAK-STAT3 / IL-17 /
TNF-NF-κB inflammatory core is real and significant but more heterogeneous across cohorts
than the original mega-analysis suggested.

---

## 2. Four issues found, and the fix for each

### Issue 1 (High) — enrichment ranked on a mega-analysis, contradicting the DE model
- **Before:** GSEA/ORA were ranked on `de_results_full.rds`, a single pooled `limma` fit
  over all samples from all studies (a *mega-analysis*). In that fit STAT3 has t = 22.3,
  p = 9×10⁻⁵⁴.
- **Why wrong:** the project's DE deliberately uses *model-then-merge* (per-study limma →
  random-effects meta) and explicitly rejects *merge-then-model*, so that between-study
  heterogeneity stays modelled. GSEA is a competitive test whose null is sensitive to
  sample size and inter-gene correlation; ranking on a pooled fit inflates every pathway
  p-value and discards the heterogeneity the meta-analysis was built to capture.
- **Fix:** re-rank genes on the **random-effects meta statistic** z = logFC/SE (signed,
  precision-weighted) from `meta_de_PPvsNN_4study.csv`. STAT3 meta-z = 4.95 (primary k=4).
  Re-run `fgsea` (eps=0, seed=42) on Hallmark, Reactome, GO:BP, KEGG.
  Output: `gsea_meta_k4_primary.csv`, `gsea_meta_k4_extended.csv`.

### Issue 2 (High) — the CAMERA meta-combination was arithmetically broken
- **Before:** `meta_pathway_camera_PPvsNN.csv` reported p-values pinned at **1×10⁻³⁰⁸**
  (double-precision underflow floor) and combined Z ≈ 60 from **k=3** studies. Reading the
  code, this came from: per-study p floored at 1e-300; Z **capped at 38**; Stouffer
  pooling weighted by **√(library size)** (sequencing depth, not statistical precision).
  A Z of 60 from 3 studies is not an interpretable significance — it is the arithmetic of
  the caps colliding.
- **Fix:** re-run per-study CAMERA on the **full k=4 set** (matching the headline DE),
  compute the signed per-study z **in the log domain** (`qnorm(logp − log 2,
  lower.tail=FALSE, log.p=TRUE)`) so it never overflows, weight by **√N** (information),
  and combine with Stouffer — reporting the pooled p also in the log domain so it never
  underflows. A **Fisher** combination on the two-sided p is reported as a
  direction-agnostic sensitivity check, and a **direction-consistency flag** records
  whether all studies agree on sign. Result: max |Zc| = 16.8 (was ~60); smallest pooled
  p ≈ 5×10⁻⁶³ (was pinned at 1e-308). Output: `meta_pathway_camera_k4_corrected.csv`.

### Issue 3 (Medium) — CAMERA used a different, smaller cohort than the headline meta
- **Before:** CAMERA ran on k=3 (SRP035988, SRP165679, SRP126422), silently dropping
  SRP065812 — the 4th study in the primary meta. It also would have included all 36
  SRP065812 lesional samples.
- **Fix:** k=4, with SRP065812 restricted to the **18 pre-adalimumab** lesional samples
  (18 post-treatment dropped), so the per-study composition (83/95, 38/28, 4/4, 16/18)
  matches the DE exactly.

### Issue 4 (Medium) — per-collection multiple-testing and uncollapsed GO redundancy
- **Before:** Hallmark, Reactome, GO:BP, KEGG were each BH-corrected separately, and
  GO:BP parent/child terms inflated the apparent hit count — the "top pathways" list
  double-counted the same biology.
- **Fix:** (a) **joint BH** across all collections tested together (9,613 unique sets);
  (b) redundancy collapse via `fgsea::collapsePathways` (removes a set whose enrichment is
  explained by a more significant one, conditional on the leading-edge genes). Result:
  1,242 up-regulated significant sets → **308 non-redundant** main pathways, annotated by
  coarse biological theme. Output: `gsea_meta_primary_jointBH.csv`,
  `pathway_nonredundant_top.csv`.

### Issue 5 (Minor, noted not changed) — ssGSEA timing normalization
- `pathway_timing_stats.csv` uses per-set min–max normalization, which makes the
  "fraction of lesional reached at PN" sensitive to the single most extreme sample. Not
  invalid; flagged for a future robustness pass (e.g. rank-based or winsorized scaling).

---

## 3. Before/after — what actually changed in the conclusions

**Effect sizes were broadly right; significances were inflated.** GSEA NES rank
concordance between the mega and meta rankings is ρ = 0.877, so the *direction and
magnitude* of pathway enrichment were reliable. What changed is the **significance of the
heterogeneous inflammatory sets**, which the mega-analysis overstated by pooling samples
and ignoring between-study disagreement:

| Hallmark set (GSEA) | mega padj (before) | meta padj (after) | verdict |
|---|---|---|---|
| E2F_TARGETS / MYC / G2M / IFN | ~10⁻³⁴ to 10⁻⁴³ | ~10⁻³¹ to 10⁻³⁷ | robust, essentially unchanged |
| IL6_JAK_STAT3_SIGNALING | 6.4×10⁻¹¹ | 8.3×10⁻⁷ | still significant, de-inflated |
| INFLAMMATORY_RESPONSE | 1.0×10⁻¹⁵ | 4.9×10⁻⁴ | significant → weak |
| TNFA_SIGNALING_VIA_NFKB | 5.2×10⁻¹³ | 0.50 | **significant → not significant** |

The corrected competitive meta (CAMERA) tells the same story: the coherent up-regulated
cassettes (interferon-α/γ, E2F, MYC, G2M, allograft rejection, IL6-JAK-STAT3, complement)
are flagged **direction-consistent** (all 4 studies agree), whereas TNFA-NF-κB (2/4 up),
INFLAMMATORY_RESPONSE (3/4), and MTORC1 (3/4) are **mixed** — exactly the sets GSEA also
de-emphasized.

**STAT3 axis specifically (your lead interest):** after joint BH + redundancy collapse,
`HALLMARK_IL6_JAK_STAT3_SIGNALING` survives as a non-redundant representative at
padj = 4.4×10⁻⁶ (NES 2.24), together with a coherent NF-κB signaling cluster
(NF-κB in B cells, TNFR2 non-canonical NF-κB, NIK) and `KERATINOCYTE_DIFFERENTIATION`.
So the STAT3/JAK-STAT programme is a **robust, reproducible** part of the lesional
signature — it is simply not the single strongest signal once the statistics are honest;
proliferation and interferon rank above it.

---

## 3b. Meta-analysis vs mega-analysis — and three combination designs

**None of the corrected results is a mega-analysis.** A *mega-analysis* pools all samples
from all studies into one matrix and fits a single model (this was the bug: the old GSEA
ranked on `de_results_full.rds`, a pooled limma fit, which is why it inflated
significance). A *meta-analysis* analyzes each study separately and then combines the
per-study results, keeping between-study heterogeneity visible. Every corrected table is a
meta-analysis.

The pathway enrichment is now reported under **three meta-analytic combination designs**,
which agree closely (Hallmark rank concordance Spearman ρ = 0.956 between designs 1 and 3):

1. **GSEA on the gene-level random-effects meta statistic** (`gsea_meta_k4_primary.csv`).
   Genes are ranked by z = logFC/SE from the DerSimonian–Laird meta (per-study →
   inverse-variance effect-size pooling); one GSEA on that combined ranking. *This is the
   primary result.*
2. **Per-study competitive test combined by p-value — CAMERA**
   (`meta_pathway_camera_k4_corrected.csv`). CAMERA within each study, then Stouffer
   weighted-Z on the per-study p-values (+ Fisher cross-check).
3. **Per-study GSEA combined by p-value** (`gsea_perstudy_combined_k4.csv`). GSEA within
   each study (ranked on that study's own t-statistic), then the per-study pathway
   p-values combined by directional Stouffer (one-sided for up-regulation) and Fisher.
   *This is the most literal "combine every study by p-value" design at the pathway level,
   and serves as a direction/consistency cross-check.*

**Note on reading design 3:** because it combines four independent tests, its combined
p-values are systematically smaller than a single GSEA's (design 1) — it should be read
for **direction consistency and rank agreement**, not as "more significant." Its
`direction_consistent` and `nUp` (of k) columns are the key output: a set that is up in
all 4 studies (e.g. IL6_JAK_STAT3, 4/4) is trustworthy; a set up in only 2/4
(TNFA-NF-κB) is heterogeneous and stays non-significant under both meta designs.

Concordance across designs for the STAT3-relevant / heterogeneous sets:

| Hallmark set | nUp (of 4) | design 1 (meta-GSEA) padj | design 3 (per-study) FDR |
|---|---|---|---|
| IL6_JAK_STAT3_SIGNALING | 4/4 | 3.7×10⁻⁶ | 4.4×10⁻¹³ |
| INFLAMMATORY_RESPONSE | 3/4 | 1.5×10⁻³ | 1.8×10⁻¹² |
| TNFA_SIGNALING_VIA_NFKB | 2/4 | 0.58 | 0.29 |

The gene-level combination itself (from METHODS.md) uses inverse-variance effect-size
pooling (DL random-effects), not p-value combination — a stronger meta-analytic method
that retains direction, magnitude, and heterogeneity (I², τ²). If a strict end-to-end
p-value-combination design is ever required, design 3 is its pathway-level realization.

---

## 4. Cohort-naming clarification

The project's "5-study" / "7-study" labels refer to project-wide study count across all
three contrasts (PP-NN, PN-NN, PN-PP), **not** 7 studies in the PP-vs-NN contrast. For
PP-vs-NN specifically there are three meta tables, all with per-gene k ≤ 4:
- **k=3:** SRP035988, SRP165679, SRP126422 (`meta_de_PPvsNN.csv`)
- **primary k=4:** k=3 + SRP065812 (`meta_de_PPvsNN_4study.csv`) — STAT3 +1.06 — **headline**
- **extended k=4:** k=3 + SRP154474 (`meta_de_PPvsNN_7study.csv`) — STAT3 +1.21

The corrected GSEA was run on both k=4 cohorts (labelled `primary` / `extended` by actual
study composition, not the misleading "7"); CAMERA and the redundancy collapse used the
primary k=4.

---

## 5. Artifacts produced

| File | What it is |
|---|---|
| `gsea_meta_k4_primary.csv` | GSEA ranked on random-effects meta z, primary k=4 (+SRP065812) |
| `gsea_meta_k4_extended.csv` | same, extended k=4 (+SRP154474) |
| `gsea_mega_vs_meta_hallmark.csv` | before/after GSEA significance, Hallmark |
| `meta_pathway_camera_k4_corrected.csv` | corrected competitive meta (Stouffer + Fisher + direction flag) |
| `camera_old_vs_corrected.csv` | before/after CAMERA (Zc, p) |
| `camera_perstudy_k4_full.rds` | per-study CAMERA inputs (checkpoint) |
| `gsea_meta_primary_jointBH.csv` | all sets, joint BH across collections |
| `pathway_nonredundant_top.csv` | 308 non-redundant pathways after collapsePathways, themed |
| `gsea_perstudy_combined_k4.csv` | design 3: per-study GSEA combined by p-value (Stouffer + Fisher + nUp/direction) |
| `gsea_perstudy_raw_k4.csv` | per-study GSEA results (NES, pval) for all 4 studies, un-combined |
| `per_study_de_k4_primary.rds` | per-study gene-level DE, primary k=4 (checkpoint) |
| `fig_pathway_corrected.png` | consolidated 4-panel figure |
| `PATHWAY_METHODS_corrected.md` | this note |

## 6. Reproducibility
- Environment: `psoriasis-r` (R 4.5.3); `fgsea`, `msigdbr` 26.1.0, `limma`, `edgeR`,
  `recount3`. Python figure in env `python` with the `figure-style` skill.
- Gene sets: MSigDB via `msigdbr` — Hallmark (H), Reactome (C2 CP:REACTOME),
  GO:BP (C5 GO:BP), KEGG (C2 CP:KEGG_LEGACY).
- Per-study expression re-derived from recount3 (Gencode v26) with the same
  filterByExpr → TMM → voom pipeline as the DE, so CAMERA's voom matrices match the DE.
- Seeds: fgsea and collapsePathways run with set.seed(42).

## 7. Caveats
- CAMERA used its default preset inter-gene correlation rather than a per-set estimated
  one; the combined statistic is competitive-meta and should be read as ranking evidence,
  not an exact FWER. The GSEA meta ranking (fgsea on meta z) is the primary enrichment
  result; CAMERA is the competitive cross-check.
- The redundancy theme labels in `pathway_nonredundant_top.csv` are rule-based keyword
  assignments for readability, not an ontology mapping; 158/308 sets fall in "Other".
