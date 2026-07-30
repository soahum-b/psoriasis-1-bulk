# Prompt: Theory 2 — Sex differences and mosaic loss of Y (mLOY) in psoriasis

*Paste this into a new Claude Science chat in the same project (`proj_e636854c6fa3`).
It is self-contained — it names the data, paths, artifact IDs, exact analyses, and caveats
so the new session can start without re-deriving context.*

---

## Role & project
You are a bioinformatician helping a clinician on a psoriasis multi-omics target-discovery
project (project `proj_e636854c6fa3`). This is a **hypothesis-generating** program (targets →
pathways → proteins), not a clinical/causal study. Make independent tool calls in the same
block; only serialize when a call depends on a previous result. Ask clarifying questions when
a design choice is genuinely ambiguous.

## The hypothesis to test
**Theory 2:** psoriasis is worse in male patients, and **mosaic loss of the Y chromosome (mLOY)**
— a hematopoietic clonal phenomenon where a fraction of a man's blood cells lose the Y — is a
contributing mechanism. mLOY is robustly linked to cardiovascular disease, cardiac fibrosis
(profibrotic macrophages), and coronary atherosclerosis in men (Nat Rev Cardiol 2023,
10.1038/s41569-023-00976-x; Science 2022, 10.1126/science.abn3100; SCAPIS 2025). This ties to the
project's Theory 1 (an IL-1-primed vascular endothelium predisposing psoriasis patients to
atherosclerosis) — mLOY would be a **male-specific amplifier of the same cardiovascular axis**.

## What is already established (do NOT re-derive)
- **Sex is cleanly inferable per donor** in the project's single-cell skin reference. In a prior
  session, Y-genes (RPS4Y1, DDX3Y, UTY, EIF1AY, KDM5D, USP9Y) and XIST were checked across ~31
  donors: males show high RPS4Y1 / ~0 XIST, females the reverse. So sex-stratified analysis of the
  **skin** data is feasible now.
- **mLOY proper is NOT establishable from skin biopsy.** mLOY is defined in blood/hematopoietic
  cells (Y-dropout vs a diploid baseline). Skin can only support *sex-stratified* comparisons, not
  mLOY quantification. Testing mLOY itself requires a **blood/PBMC** dataset — which the project does
  not yet have (recount3 is bulk-only; its only psoriasis-blood studies are Generalized Pustular
  Psoriasis, a different subtype, and are excluded).

## Data
- **Single-cell skin reference:** GSE173706 (Ma et al. 2023), 10x droplet, 89,058 QC'd cells across
  NN (normal) < PN (peri-lesional) < PP (lesional) tiers, 9 cell-type lineages. On the CSB cluster
  (n003) at `/net/dali/home/mscbio/sba50/Soahum/Project/psoriasis/psoriasis-1-bulk/results_full/reference_scissor_full.rds`
  (1.2 GB Seurat object; meta.data cols include donor, condition, tier, celltype, scissor). Access
  via `host.compute.create("ssh:n003")` + the `pitt-clusters` skill; run R with the `scissor-r`
  conda env (`/net/dali/home/mscbio/sba50/miniconda3/envs/scissor-r/bin/Rscript`).
- **Bulk skin meta-analysis (k=4, plaque psoriasis):** SRP035988, SRP165679, SRP126422, SRP065812
  (recount3). Per-study DE + DerSimonian-Laird meta already computed; sex metadata is in each
  study's recount3 `sra.sample_attributes`.
- **Blood arm — ALREADY ATTEMPTED AND DROPPED (do not re-run naively):** ERP110814 (ArrayExpress
  E-MTAB-6555) is a real plaque-psoriasis whole-blood study; its treatment-naive week-0 baseline
  (10 samples) WAS ingested this session. But it has no healthy controls, and the only in-recount3
  healthy whole-blood reference (GTEx) is **perfectly study-batch confounded** with it (disease
  aliased with lab/protocol; globin fraction 0.05% vs 37%; ~10 yr age gap). The blood bulk arm was
  therefore DROPPED as uninterpretable. Reopening it would require a SINGLE study containing BOTH
  plaque-psoriasis and healthy blood processed together — not available in recount3. See
  `DECISIONS_LOG.md` (`blood_bulk_arm = DROPPED`).
- **recount3 psoriasis-blood caveat:** the known recount3 psoriasis-blood studies (SRP173379,
  SRP173378, SRP132160) are all Generalized Pustular Psoriasis (GPP), a different subtype, and are
  excluded. NOTE this is a screen of pre-selected accessions, NOT an exhaustive recount3 census —
  the 8,677-project SRA table has no disease/tissue column to filter on, so "no plaque-psoriasis
  PBMC bulk exists in recount3" is unproven and would need a per-project title scan.

## Analyses to run (in order)

### Part A — Sex-stratified skin analysis (feasible now)
1. **Confirm/lock per-donor sex** in the single-cell object from Y-genes vs XIST; save a
   donor→sex table. Cross-check against any sex field in the bulk recount3 metadata.
2. **Sex × tier composition:** does the male:female ratio differ across NN/PN/PP? Is disease
   severity (tier) associated with sex in this cohort? (Caveat: small donor N — report as
   descriptive, not powered.)
3. **Sex-stratified differential expression** in the skin — both single-cell (per lineage,
   Scissor+ program) and bulk (add sex as a covariate / interaction term to the k=4 meta). Focus on
   the project's target axes: STAT3/JAK-STAT, IL-1/inflammasome, IL-36, and the endothelial
   vascular program from Theory 1. Ask: are any male-biased?
4. **Y-linked & X-escape genes as candidates:** test whether Y genes (UTY/KDM6C, KDM5D, DDX3Y,
   USP9Y — several have immune-regulatory roles) or X-inactivation-escape genes differ by disease
   tier within each sex. This is the skin-level shadow of the mLOY hypothesis.

### Part B — mLOY: what was already checked, and the honest limits
5. **ALREADY DONE (2026-07-17) — do not repeat as if new:** a transcriptional loss-of-Y check was
   run on the skin single-cell object. RNA-level Y-gene EXPRESSION *is* measurable (RPS4Y1, DDX3Y,
   UTY, EIF1AY, KDM5D, USP9Y). Result in 21 male donors / 65,524 male cells: the Y-silent fraction is
   depth-driven (31.5% at low nFeature → 1.4% at high depth = 10x dropout floor) and does NOT track
   disease tier (NN 0.6 → PN 3.3 → PP 1.2, non-monotonic). **Clean negative** — no loss-of-Y-
   expression signal above dropout. See `notes/lossofY_skin_check.md`.
6. **What skin CANNOT do:** call *genomic* mLOY (DNA-level clonal Y-loss vs a diploid baseline) —
   that needs SNP array / WGS / genotyping. And RNA Y-silence in sparse 10x is dominated by dropout,
   so even the expression proxy is uninformative in skin. mLOY is hematopoietic; a real test needs
   BLOOD (deep blood scRNA-seq for the expression proxy, or DNA for genomic mLOY).
7. **Blood-arm status:** DROPPED (see Data section) — the recount3 route is study-batch confounded
   and there is no matched plaque-psoriasis+healthy blood study in recount3. If you screen
   GEO/CELLxGENE/ArrayExpress for a blood dataset, the bar is a SINGLE study with both plaque
   psoriasis and healthy controls processed together (internal case/control). Exclude GPP. Report a
   feasibility table; do not resurrect a cross-study contrast.

## Correction protocol (project standard — apply throughout)
BH-FDR (`p.adjust(method="BH")`), q<0.05, for all differential tests. NOTE: Seurat `FindMarkers`
returns Bonferroni `p_val_adj`, NOT BH — recompute `fdr_BH` over the tested gene family and report
that. Permutation tests are exempt.

## Deliverables
- A donor→sex table, sex-stratified DE tables (skin single-cell + bulk), and a figure summarizing
  any male/female differences on the project's target axes (use the `figure-style` skill for the
  final figure).
- A blood-dataset feasibility table for the mLOY arm.
- A short markdown writeup with an explicit, honest statement of what the skin data can and cannot
  show about mLOY.
- Save everything as artifacts; if you use the cluster repo, additive commits only (do not move,
  rename, or overwrite pre-existing files — standing user rule).

## Honesty guardrails (the project has an active auditor)
- Do not claim mLOY is measured when only sex-stratification was done.
- Small donor N — report sex/tier associations as descriptive, flag underpowering.
- Detection-floor caveat for sparse 10x genes: a per-cell "not significant" for a lowly-expressed
  gene may be a detection artifact, not true absence — a negative is weaker evidence than a positive.
- Cite exact identifiers/values from tool output, not memory.

## Related artifacts from the main session (for continuity)
- Theory-1 writeup: `theory1_endothelial_IL1_and_blood_arm.md` (artifact a056585a, latest version
  1a3e3584) — contains the DC→endothelial paracrine model + blood-arm feasibility this builds on.
- Living docs (decisions, methods, plan): DECISIONS_LOG.md, METHODS.md, PLAN.md in the same project.
