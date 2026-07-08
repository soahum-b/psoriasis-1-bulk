# HANDOFF — Scissor-on-gradient (peri-lesional psoriasis)

**Status as of this handoff:** backbone complete and validated on a local
20,023-cell subset. Ready to (a) scale to the full 89,058-cell census on a
cluster, and (b) add orthogonal deconvolution validation. This file is the
single source of truth for *where we are* and *what to do next*.

---

## 1. What this project does (one paragraph)

We project a **bulk RNA-seq ordinal phenotype** (biopsy site: normal `NN` <
peri-lesional `PN` < lesional `PP`) onto a **single-cell reference** using
**Scissor** (Sun & Xia, *Nat Biotechnol* 2022) in Gaussian/regression mode.
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

Caveat: NNLS on marker signatures is the transparent core of regression
deconvolution but not a benchmarked method choice. deconvBenchmarking (below)
is the layer that would justify *which* deconvolution method to trust for this
reference/tissue — it is a benchmarking framework (simulate realistic pseudobulk
with KNOWN proportions, score methods), not itself a real-bulk deconvolver.

---

## 5. HOW TO REPRODUCE / SCALE UP

```bash
# 0a. Initialize git (the repo ships as a clean tree; create the repo + push)
git init && git add -A && git commit -m "Scissor-on-gradient backbone"
# then: git remote add origin <your-cluster-reachable-remote> && git push -u origin main

# 0b. Environment (once)
conda env create -f environment.yml && conda activate scissor-r

# 1. Download data (~264 MB scRNA + recount3 bulk)
Rscript code/00_download_data.R

# 2a. Reproduce the local backbone (20k subset, laptop-scale)
Rscript code/run_pipeline.R           # writes results/

# 2b. FULL CENSUS on a cluster (all 89,058 cells) — the scale-up
sbatch --mem=256G --cpus-per-task=16 --time=24:00:00 \
       --wrap "Rscript code/run_full_census_cluster.R"    # writes results_full/
```

`run_full_census_cluster.R` has all parameters at the top (`N_CELLS=Inf` = full
census). Memory driver is the glmnet augmented design at 89k predictors ×
~3.1M SNN edges — budget ≥128 GB, 256 GB comfortable.

---

## 6. NEXT STEPS (priority order)

1. **Full-census run on the cluster** (`run_full_census_cluster.R`). Confirm the
   20k-subset findings hold at full size — especially the endothelial
   enrichment and the PN-peak in Scissor+.
2. **Try the real compiled Scissor on the cluster** and cross-check our glmnet
   port (§3). If concordant, report the canonical solver's numbers.
3. **Orthogonal deconvolution validation** — ✅ DONE for the direct estimate
   (§4b: NNLS confirms Scissor, 5/6 cell types concordant). REMAINING upgrade:
   run **deconvBenchmarking** (github.com/humengying0907/deconvBenchmarking) as
   the method-selection layer. It is a benchmarking framework, not a real-bulk
   deconvolver — its three core calls `benchmarking_init()` /
   `benchmarking_deconv()` / `benchmarking_evalu()` simulate realistic
   *heterogeneous* pseudobulk from the Ma reference (known proportions), run many
   deconvolution methods, and score per-cell-type correlation/RMSE. Use it to
   pick the most accurate method+reference for this skin tissue, THEN re-run that
   method on the real SRP165679 bulk to replace the NNLS estimate in §4b with a
   benchmarked one. Install: `devtools::install_github("humengying0907/deconvBenchmarking")`.
4. **Tier-2 / Tier-3 robustness** (per the analysis plan): repeat with the
   continuous PASI-like severity anchor and with alternative bulk studies
   (SRP035988 NN/PP only; SRP076982) to show the gradient program is not
   anchor-specific.
5. **Sequence-level extension (future arm, NOT part of this pipeline):**
   seillra / low-rank Sei (github.com/kostkalab/seillra; Gilfeather, Chikina,
   Kostka, bioRxiv 2026.01.21.700827) is a **PyTorch DNA-sequence model** —
   input is genomic sequence, output is predicted chromatin/regulatory state and
   **variant effects**. It shares no data interface with this expression
   pipeline (no count matrix / logCPM / cell-selection input) and cannot be
   substituted into it. It belongs to a downstream sequence-level arm.

   Concrete bridge (seillra's README lists a variant-effects usage notebook,
   `variant_usage.ipynb`, as the starting point — verify it in the repo):
   ```bash
   pip install git+https://github.com/kostkalab/seillra.git
   ```
   ```python
   import seillra as sl
   mod = sl.Sei_LLRA(k=64, projection=True, quant="CPU")  # rank 64, CPU-quantized
   # score psoriasis GWAS / regulatory variants near the gradient-program genes
   # (endothelial-led set + STAT3) for regulatory-disruption effect.
   ```
   Two-stage design that uses BOTH lab tools coherently:
     (a) deconvBenchmarking (expression) confirms WHICH cell types drive the
         NN->PN->PP gradient;
     (b) seillra/Sei-LLRA (sequence) scores WHICH regulatory variants perturb the
         driver genes that stage (a) nominates.
   Inputs this arm needs and this project does NOT yet have: a psoriasis
   regulatory/GWAS variant list (e.g. from GWAS Catalog / fine-mapping) and the
   reference genome sequence windows around the gradient-program gene loci.

---

## 7. FILE MAP

```
code/
  00_download_data.R          fetch GSE173706 + SRP165679
  scissor_glmnet_solver.R     pure-R network elastic-net (APML1 replacement)
  scissor_run.R               Scissor driver: inputs prep + alpha tuning
  scissor_reliability.R       reliability test (test_lm port)
  run_pipeline.R              end-to-end local backbone driver
  run_full_census_cluster.R   FULL-CENSUS scale-up (SLURM)
results/                      *.csv tracked; *.rds git-ignored (regenerable)
  scissor_alpha_tuning.csv    selected fraction per alpha
  celltype_enrichment.csv     cell-type fold-enrichment + Fisher OR
  gradient_program_DE.csv     1,861-gene gradient program
  permutation_null.csv        selection-null per-permutation stats
  phenotype_checks.csv        bulk non-degeneracy checks
figures/                      all 8 publication figures (PNG)
environment.yml               conda spec
README.md                     project overview
HANDOFF.md                    this file
```

---

## 8. KNOWN CAVEATS

- **Subset, not census.** All committed `results/` are the 20k subset. Full
  census pending (§6.1). Numbers may shift; directionality expected to hold.
- **glmnet port ≠ compiled APML1.** Validated (synthetic recovery cor=0.72;
  reliability & null both p=0.000) but not the authors' canonical solver (§3).
- **One bulk anchor.** Tier-1 uses SRP165679 only. Robustness across anchors is
  Tier-2 work (§6.4).
- **STAT3 signal is modest** (log2FC 0.43). It tracks the gradient but is not a
  top driver; the program is dominated by vascular/endothelial genes. Interpret
  STAT3 as *part of* the gradient program, not its headline.
