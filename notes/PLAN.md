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
- [x] Retrieve + validate full-census numbers vs backbone (alpha=0.20, 14.67% selected; endothelial OR 5.24 top lineage; STAT3 n.s.)
- [x] Update WHITEPAPER.md §4 with full-census figures + BH-FDR correction (STAT3 q=0.13)
- [ ] Re-render WHITEPAPER HTML/PDF from updated .md
- [ ] Compiled-Scissor (APML1) cross-check
- [ ] Regenerate §4 figures from the 89k object (text/figure mismatch: figures still depict backbone)

## Phase 3 — Theory 1: IL-1-responsive vascular endothelium (→ atherosclerosis)
- [x] Reframe to DC→endothelial paracrine model; IL-1 receptor localization by cell type
- [x] IL1R1+ endothelium enrichment in Scissor+ set (OR 2.32) + downstream co-expression
- [x] Endothelium PP-vs-NN BH-FDR panel (IL1R1 early-high, down in lesional)
- [x] Figure + writeup (`fig_theory1_endothelial_IL1.png`, `notes/theory1_endothelial_IL1_and_blood_arm.md`)
- [ ] Literature integration into references.md (psoriatic march; IL-1 early; S1PR3-STAT3)

## Phase 4 — Blood / circulation arm (two-part story: skin + blood)
- [x] recount3 blood screen → known psoriasis-blood = GPP-only (EXCLUDED per user)
- [x] ERP110814 baseline feasibility (plaque psoriasis blood; wk-0 recoverable via SDRF)
- [ ] **Include ERP110814 baseline** (10 tx-naive PP blood; join SDRF time → recount3 ERR ids)
- [ ] Screen GEO/CELLxGENE for plaque-psoriasis PBMC single-cell + a healthy-blood control
- [ ] Systemic IL-1/myeloid-priming test in blood (closes the paracrine loop's source end)

## Phase 5 — Theory 2: sex / mosaic loss of Y (male severity)
- [x] Sex-inference feasibility in skin data (Y-genes vs XIST) — VIABLE
- [ ] Sex-stratified pass on skin data (feasible now; down payment on theory 2)
- [ ] mLOY proper — requires the blood arm (deferred per user)

## Phase 6 — Integration & downstream
- [ ] Benchmarked deconvolution (deconvBenchmarking) as Stage 1 of two-stage design
- [ ] Druggability / docking follow-up on target hits (STAT3 6NJS, JAK3 5LWM, RORC 5APH)
- [ ] Sequence-level arm (Sei-LLRA) — regulatory-variant scoring (GPU)
- [ ] STAT3 isoform ecological correlation (indirect, optional)
