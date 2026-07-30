# Auto-extracted generating script
# Produces: README.md
# Conda env: python   (run with this environment activated)
# Inputs (expected alongside / in data/): celltype_enrichment.csv, reliability_test.rds, scissor_result.rds, scissor_run.R, scissor_reliability.R, scissor_glmnet_solver.R, permutation_null.rds
# Source artifact version: 088da65c-e90a-4e4b-a4ec-b594f48d749f
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

## Headline result (20k-cell backbone)

- Gradient-tracking cells (**Scissor+**) are **monotonic on NN<PN<PP** and
  **peak at the peri-lesional (PN) tier** — the intermediate-state hypothesis.
- Both significance controls pass at **p = 0.000**: reliability test (real
  CV-MSE 0.147 vs null 0.779) and selection permutation null (directionality
  gap 0.638 → −0.305 under shuffled labels).
- **Endothelial cells** are the dominant gradient-trackers (5.2× enriched,
  OR 11.3); the 1,861-gene program is vascular-led. **STAT3** is significantly
  up in Scissor+ (log2FC 0.43, padj 0.018) as part of that program.

## Method

Faithful Scissor pipeline with the compiled `APML1` solver replaced by a
**pure-R glmnet network-regularized elastic net** (sparse edge augmentation of
the symmetric-normalized graph Laplacian). See `code/scissor_glmnet_solver.R`
and HANDOFF.md §3 for the derivation and the caveat about the canonical solver.

## Layout

`code/` scripts · `results/` tables (CSV tracked, RDS regenerable) ·
`figures/` 8 PNGs · `environment.yml` · `HANDOFF.md` (status & next steps).
"""

with open('README.md', 'w') as f:
    f.write(content)