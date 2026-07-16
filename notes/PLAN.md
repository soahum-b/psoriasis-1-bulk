# Plan (mirror of the live plan — for offline review)

Static snapshot. Update checkboxes as steps complete.

## Phase 1 — Bulk meta-analysis (STAT3 / IL-1β)  [largely done]
- [x] 7-study meta-analysis + DE (recount3)
- [x] IL-1β/STAT3 literature overview
- [x] PP-vs-NN focused contrast, peri-lesional excluded
- [x] Study-set expansion k=3 → k=4 (add SRP065812; exclude ERP110816)
- [x] k=4 forest plot + k=3-vs-k=4 comparison
- [x] IL1B_STAT3_expanded_analysis.md writeup section
- [ ] Healthy-vs-uninvolved control-type moderator (exploratory, under-powered)
- [ ] Fold Arm-A section into a meta-analysis whitepaper

## Phase 2 — Scissor single-cell arm
- [x] 20k-cell backbone (alpha tuning, significance controls, deconv validation)
- [x] Full-census run on CSB (job 56882314, 89,058 cells)
- [ ] **Retrieve + validate full-census numbers vs backbone**  ← in progress
- [ ] Update WHITEPAPER.md §4 with full-census figures; re-render HTML/PDF
- [ ] Compiled-Scissor (APML1) cross-check
- [ ] Rerun celltype-enrichment / phenotype-checks / deconv-concordance vs corrected anchor

## Phase 3 — Integration & downstream
- [ ] Benchmarked deconvolution (deconvBenchmarking) as Stage 1 of two-stage design
- [ ] Druggability / docking follow-up on target hits (STAT3 6NJS, JAK3 5LWM, RORC 5APH)
- [ ] Sequence-level arm (Sei-LLRA) — regulatory-variant scoring (GPU)
- [ ] STAT3 isoform ecological correlation (indirect, optional)
