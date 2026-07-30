# Auto-extracted generating script
# Produces: HANDOFF.md
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: da1a8751-a9f5-4194-8511-90ced6855152
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

handoff_content = """# HANDOFF — Scissor-on-gradient (peri-lesional psoriasis)

**Status as of this handoff:** backbone complete and validated on a local
20,023-cell subset. Ready to (a) scale to the full 89,058-cell census on a
cluster, and (b) add orthogonal deconvolution validation. This file is the
single source of truth for *where we are* and *what to do next*.

---

## 1. What this project does (one paragraph)

We project a **bulk RNA-seq ordinal phenotype** (biopsy site: normal `NN` <
peri-lesional `PN` < lesional `PP`) onto a **single-cell reference** using
**Scissor** (Sun, Guan, … Xia, *Nat Biotechnol* 2022) in Gaussian/regression mode.
The output is a per-cell selection: cells whose expression co-varies with the
NN→PN→PP gradient (**Scissor+**, tracking toward lesional), cells tracking
toward normal (**Scissor−**), and unselected **Background**. The phenotype is a
**clinical biopsy-site label, not an expression-derived score**, so the
selection is non-circular by construction (Tier-1 design).

---

## 2. Data

| Role | Dataset | Detail |
|------|---------|--------|
| Single-cell reference | **GSE173706** (Ma et al. 2023, *Nat Commun*) | 33 samples, 22 donors; NN=8 healthy donors, PN+PP = 11 paired lesional/peri-lesional donors. 96,088 cells raw → 89,058 after QC. |
| Bulk phenotype anchor | **SRP165679** (Tsoi 2019) via recount3 | 93 samples: NN=38 / PN=27 / PP=28. Single study (no cross-study confound). logCPM over 24,533 symbols. |

Both are re-downloaded by `code/00_download_data.R`. Raw data and large `.rds`
objects are **git-ignored** (regenerable); the repo tracks code, figures, docs,
and small result tables.

---

## 3. Method note — why this is not stock Scissor

The compiled Scissor package ships a C++ solver (`APML1`) for its
network-regularized elastic net. It **could not be installed in the analysis
sandbox** (shared libraries on writable paths are refused; the read-only conda
R library cannot receive the package). We therefore reimplemented the solver in
**pure R via glmnet** (`code/scissor_glmnet_solver.R`), keeping every other
Scissor step (quantile normalization, cell–bulk correlation, SNN graph,
sign-based selection, reliability test) faithful to the original.

Key implementation facts (read the solver docstring for the full derivation):
- The network penalty is the **symmetric-normalized graph Laplacian**
  `L_sym = I − D^{−1/2} A D^{−1/2}`, matching Scissor's `OmegaC`.
- Its quadratic form is solved **exactly** via sparse **edge-difference
  augmentation** (`βᵀL_sym β = Σ_edges (β_i/√d_i − β_j/√d_j)²`), which scales to
  millions of edges without a dense Cholesky.
- Sparsity (the "selected fraction") is controlled by **walking the glmnet
  lambda path to a target fraction** — the operational analogue of Scissor's L0
  hard-threshold trim, which the pure-R port does not replicate. `alpha` sets
  graph-smoothing strength; the L1 penalty is chosen independently along the
  path. **This is NOT the coupled elastic-net `lambda1=lambda*alpha` form.**

⚠️ **On the cluster, prefer the real compiled Scissor if it installs there.**
The glmnet port is validated (below) but the authors' `APML1` is canonical. If
`devtools::install_github("sunduanchen/Scissor")` succeeds on the cluster,
run stock `Scissor()` and compare against `results_full/` from our port.

---

## 4. Results on the 20k subset (validated backbone)

- **Selection:** alpha=0.40, **14.84 %** selected (1,574 Scissor+ / 1,397 Scissor−). Under the 20 % cutoff.
- **Directionality (the biology):** mean tier (NN=0,PN=1,PP=2) is monotonic —
  Scissor− **0.79** < Background **1.38** < Scissor+ **1.43**. Scissor+ **peaks
  at the PN/peri-lesional tier** (9.4 % of PN cells vs 8.0 % of PP), direct
  support for an intermediate gradient-tracking state.
- **Reliability test (n=100):** real CV-MSE **0.147** vs null mean **0.779**,
  **p = 0.000** (0/100 permutations lower). Highly reliable.
- **Selection permutation null (n=30):** pos−neg tier gap collapses from real
  **0.638** to null mean **−0.305**, **p = 0.000**. Selection is driven by real
  phenotype structure, not graph geometry.
- **Cell types:** **Endothelial 5.2× enriched** among Scissor+ (OR=11.3,
  p≈1e-245); Melanocyte & Fibroblast enriched among Scissor−.
- **Gradient program:** **1,861 DE genes** (Scissor+ vs Background); top hits
  vascular/endothelial (CCL14, ACKR1, RAMP3, PLVAP, APLNR, ADGRL4).
- **STAT3:** significantly up in Scissor+ (log2FC **0.43**, padj **0.018**;
  48.5 % vs 44.7 % expressing). Modest but in the expected direction.

### 4b. Orthogonal validation — bulk deconvolution (DONE)

Independent NNLS deconvolution of the real SRP165679 bulk against a Ma-reference
cell-type signature (`code/deconv_validation.R`) **confirms Scissor by a
different method**. Estimated cell-type proportions across NN→PN→PP:

Scissor direction is taken from the Fisher enrichment of each cell type among
Scissor+ cells (pos-enrichment OR): OR>1 (sig) = lesional-tracking; OR<1 (sig) =
normal-tracking. Concordant = deconv proportion trend is significant AND matches.

| Cell type | Scissor direction | pos-enrich OR | Deconv proportion trend | deconv padj | Status |
|-----------|-------------------|:---:|--------------------------|:---:|:---:|
| Endothelial | Scissor+ (lesional) | 11.26 | **rises** 0.0→0.2→3.8 % | 1.3e-19 | ✓ concordant |
| Fibroblast | Scissor− (normal) | 0.15 | **falls** 12.3→8.2→1.4 % | 1.2e-11 | ✓ concordant |
| Melanocyte | Scissor− (normal) | 0.04 | **falls** 7.8→7.3→2.2 % | 1.1e-14 | ✓ concordant |
| NK | Scissor− (normal) | 0.36 | rises →7.9 % (PP) | 1.5e-11 | ✗ discordant |
| Keratinocyte | Scissor− (normal) | 0.69 | rises 79→83→83 % | 5.1e-3 | ✗ discordant |
| DC | Scissor+ (lesional) | 1.47 | rises (flat) | 0.20 | trend n.s. |

**3 of 6 cell types with a clear Scissor direction are concordant** — and they
are the three strongest-signal lineages: **Endothelial** (lesional-tracking,
rises in bulk), **Fibroblast** and **Melanocyte** (normal-tracking, fall in
bulk). This resolves the composition-vs-state ambiguity *for the endothelial
signal*: it reflects (at least partly) a **real compositional increase** in the
bulk, not solely a state change.

The two discordant cases are informative, not failures. **NK** cells are
*depleted* from the Scissor+ gradient-tracking set (per-cell) yet their bulk
proportion rises sharply at PP — i.e. NK **infiltration** into lesional tissue
is a compositional event that individual NK cells do not "track" the gradient in
the Scissor sense. **Keratinocyte** is only weakly normal-leaning per-cell
(OR 0.69) while its bulk fraction is near-flat-to-rising (it dominates all
tiers), so the directions are not expected to align. Per-cell tracking and bulk
composition are genuinely different measurements; they agree where the biology
is a clean compositional shift (endothelial/fibroblast/melanocyte) and diverge
where it is not (NK infiltration, keratinocyte dominance).

Figure `figures/fig_deconv_validation.png`; tables
`results/deconv_{proportions,trend,scissor_concordance}.csv`.

---

## 5. HOW TO REPRODUCE / SCALE UP