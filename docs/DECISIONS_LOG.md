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
- `exclude_ERP110816  [bulk, ruling, by:agent]`  (all patients on etanercept/anti-TNF confounds TNF→IL-1→STAT3; no timepoint field to isolate treatment-naive baseline)
- `pooled_normal_design = lesional vs all-normal pooled + control-type moderator + donor random effect  [bulk, design, by:user]`  (pooling healthy+uninvolved is not double-dipping; pre-specify one primary contrast)
- `whitepaper_split = k=4 IL-1β/STAT3 documented as a separate meta-analysis section, NOT folded into the Scissor WHITEPAPER  [doc, ruling, by:agent]`  (different arm; avoids muddling two stories)

### Isoforms
- `isoform_by_celltype = infeasible per-cell with GSE173706  [ruling, by:agent]`  (10x short-read cannot resolve STAT3α/β C-terminal splice; needs Smart-seq/long-read; indirect ecological correlation available as fallback)
