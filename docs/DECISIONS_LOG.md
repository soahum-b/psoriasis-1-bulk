# Decisions Log (living)

Record every non-trivial choice. Tag SCOPE so a later choice does not silently override a
differently-scoped earlier one.

**Format:** `decision = value  [scope, kind, by:user|agent]  (why)`

## Log

### Scissor single-cell arm
- `phenotype = ordinal biopsy-site (NN<PN<PP)  [scissor, design, by:user]`  (reframe from 2-state; captures the advancing peri-lesional margin)
- `non_circularity = clinical biopsy-site label, not molecular score  [scissor, design, by:agent]`  (a flagged cell cannot be a restatement of how it was labelled)
- `bulk_anchor = SRP165679 only  [scissor, design, by:agent]`  (only recount3 study balancing all three tiers at usable depth)
- `solver = pure-R glmnet reimplementation, not compiled APML1  [scissor, method, by:agent]`  (APML1 wouldn't build; validated equivalent, canonical cross-check staged)
- `alpha = 0.40  [scissor, param, by:agent]`  (tuning grid; 14.84% cells selected on backbone, under 20% cutoff)
- `full_census = CPU partition (dept_cpu), not GPU  [scissor, compute, by:agent]`  (coordinate-descent elastic net is CPU-bound; glmnet has no GPU backend)

### Deconvolution / two-stage design
- `deconv_scissor = two-stage (deconv maps landscape → Scissor resolves sub-state)  [design, by:user]`  (composition vs state; Scissor+ & proportion-flat = state-driven = druggable-target home)
- `deconv_stage1 = benchmarked deconvolution (deconvBenchmarking) promoted to Stage 1  [design, by:agent]`  (staged upgrade: proportions as Scissor covariates)

### Bulk IL-1β/STAT3 meta-analysis arm
- `focus_contrast = PP-vs-NN, peri-lesional excluded  [bulk, design, by:user]`  (PN contrast underpowered: only k=2 studies, 1/49 panel genes significant)
- `study_constraint = must contain both NN and PP (NOT three-tier)  [bulk, ruling, by:agent]`  (corrects an earlier misstatement; SRP035988 has PN=0 yet qualifies)
- `expand_studies = add SRP065812 → k=4  [bulk, design, by:agent]`  (16 healthy NN + 18 pre-adalimumab PP; drop 18 post-tx)
- `exclude_ERP110816_from_PPvsNN  [bulk, ruling, by:agent]`  (skin, all on etanercept/anti-TNF; **CORRECTION**: timepoint IS recoverable — `Factor Value[time]` in the ArrayExpress SDRF E-MTAB-6556, dropped by recount3 — so week-0 baseline can be isolated; real reason it stays out of the PP-vs-NN meta is that baseline has only lesional+uninvolved, **no healthy NN**. Available for a PP-vs-PN sensitivity check only.)
- `pooled_normal_design = lesional vs all-normal pooled + control-type moderator + donor random effect  [bulk, design, by:user]`  (pooling healthy+uninvolved is not double-dipping; pre-specify one primary contrast)
- `whitepaper_split = k=4 IL-1β/STAT3 documented as a separate meta-analysis section, NOT folded into the Scissor WHITEPAPER  [doc, ruling, by:agent]`  (different arm; avoids muddling two stories)

### Theory 1 — IL-1-responsive vascular endothelium (psoriasis→atherosclerosis)
- `theory1_reframe = DC→endothelial paracrine model  [analysis, ruling, by:agent]`  (IL-1β is DC-produced, not endothelial (48% of DCs vs 0.2% endo); vessels carry the RECEPTOR (IL1R1 ~29% of endo). Source and responder are spatially separate.)
- `theory1_finding = IL1R1+ endothelium IS the gradient-tracking, STAT3/inflammasome-active population  [analysis, result, by:agent]`  (within endothelium, IL1R1+ cells 2.3× enriched in Scissor+ (Fisher OR 2.32, p 1.4e-38); IL1R1+ vs IL1R1- co-express more STAT3 66/40, IL6 31/13, NFKB1 28/13, CASP1 20/10, GSDMD 51/30, PYCARD 38/19. Figure: `figures/fig_theory1_endothelial_IL1.png`.)
- `theory1_direction = IL-1 responsiveness is early/constitutive, not late  [analysis, ruling, by:agent]`  (endothelial IL1R1 DOWN in lesional vs normal (log2FC −1.12, BH q≈0; tier 51%→29%→27%) — consistent with IL-1 as an initiating signal that desensitizes as lesions establish.)

### Circulation / blood arm ("two-part story": skin + blood) — DROPPED 2026-07-17
- `blood_bulk_arm = DROPPED  [blood, ruling, by:user]`  (attempted ERP110814 baseline (10 tx-naive plaque-psoriasis blood, wk-0 via SDRF) vs a healthy-blood reference. The only in-recount3 healthy whole-blood reference is GTEx, which is **perfectly study-batch confounded** with ERP110814: disease status is 100% aliased with lab/protocol (QMUL Globin-Zero vs GTEx PAXgene undepleted — globin fraction 0.05% vs 37%). Case/control differences are inseparable from batch; matching/covariates cannot fix a confounder collinear with the exposure. User: "double dipping … doesn't make sense." Correct call — dropped.)
- `exclude_GPP = do NOT include generalized pustular psoriasis  [blood, constraint, by:user]`  (still stands independently; GPP is a distinct subtype. recount3 psoriasis-blood = SRP173379/78/SRP132160, all GPP.)
- `recount3_is_bulk_only  [data, ruling, by:agent+user]`  (retained as a general data fact.)
- `blood_reopen_condition  [blood, note, by:agent]`  (if the blood question is ever revisited, it needs a SINGLE study containing both plaque-psoriasis and healthy blood processed together — internal case/control, no cross-study batch. Not available in recount3; would require GEO/ArrayExpress. Not pursued.)

### Theory 2 — sex / mosaic loss of Y (mLOY), male severity
- `theory2_sex_feasible_in_skin = YES  [analysis, ruling, by:agent]`  (sex cleanly inferable per donor from Y-genes RPS4Y1/DDX3Y/UTY/EIF1AY/KDM5D/USP9Y vs XIST; ~20 male / ~10 female donors → sex-stratified skin analysis feasible now.)
- `theory2_mLOY_needs_blood  [analysis, ruling, by:agent]`  (mLOY is hematopoietic Y-dropout vs diploid baseline — cannot be established from skin biopsy; requires the PBMC/blood arm. Deferred per user "we can test it later".)
- `LOY_skin_check = clean NEGATIVE  [analysis, result, by:agent]`  (2026-07-17 exploratory: RNA-level Y-gene-expression proxy in skin. 21 male donors / 65,524 male cells. Y-silent fraction is depth-driven (31.5% at low depth → 1.4% at high depth = dropout floor) and does NOT track disease tier (NN 0.6 → PN 3.3 → PP 1.2, non-monotonic). No loss-of-Y-expression signal above technical dropout. Note: RNA Y-gene EXPRESSION is measurable (this is what we did); GENOMIC mLOY still needs DNA/blood. See `notes/lossofY_skin_check.md`.)

### Isoforms
- `isoform_by_celltype = infeasible per-cell with GSE173706  [ruling, by:agent]`  (10x short-read cannot resolve STAT3α/β C-terminal splice; needs Smart-seq/long-read; indirect ecological correlation available as fallback)
