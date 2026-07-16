# Theory 1: IL-1-responsive vascular endothelium and the psoriasis→atherosclerosis axis

*Single-cell full-census analysis (89,058 cells, GSE173706) + recount3 blood-arm feasibility.
All single-cell DE uses BH-FDR per protocol. Generated in-session; figures/tables saved as artifacts.*

## The hypothesis (clinician's framing, reframed by the data)

Original theory: *IL-1β is high on the vessels/circulation → drives the IL-1 pathway → increases
STAT3 and pro-inflammatory inflammasome activity*, and this vascular inflammation is the mechanism
**predisposing psoriasis patients to atherosclerosis** (the "psoriatic march").

The data required one correction to the source of IL-1β and one to the direction of the receptor,
but the **core mechanism is supported at the single-cell level**.

## What the data shows

**1. IL-1β is not made by endothelium — it is DC-restricted.**
IL1B is expressed in only 1.85% of all cells; 48.4% of dendritic cells express it vs <1% of every
other lineage (endothelial 0.2%). So the vessels do not *produce* IL-1β. This makes the mechanism a
**paracrine model**: DC- (and systemically-) produced IL-1β acting on an IL-1-responsive vascular
compartment — a stronger, more specific claim than co-localized production.

**2. The endothelium is IL-1-responsive.**
IL1R1 (the signaling receptor) is expressed in ~29% of endothelial cells — comparable to DCs (29%)
and fibroblasts (31%), and higher than keratinocytes (17%). The vascular compartment is wired to
receive IL-1.

**3. IL1R1⁺ endothelial cells ARE the progression-tracking, inflammatory population.**
Within endothelium, IL1R1⁺ cells are **2.3× enriched in the Scissor+ gradient-tracking set**
(Fisher OR 2.32, p = 1.4×10⁻³⁸; 40.5% of IL1R1⁺ vs 22.7% of IL1R1⁻ endothelial cells are
gradient-tracking). And IL1R1⁺ endothelial cells carry markedly more of the downstream program than
IL1R1⁻ endothelial cells:

| Program gene | IL1R1⁺ endo (% expr) | IL1R1⁻ endo (% expr) |
|---|--:|--:|
| STAT3  | 65.6 | 39.8 |
| GSDMD  | 50.9 | 30.4 |
| PYCARD | 38.2 | 18.8 |
| IL6    | 31.0 | 12.6 |
| NFKB1  | 28.1 | 12.8 |
| CASP1  | 19.5 |  9.7 |

This is the causal chain made visible: the endothelial cells expressing the IL-1 receptor are the
same ones running STAT3, IL-6, NF-κB and the pyroptotic/inflammasome machinery — the "IL-1-primed,
pro-atherogenic endothelium" the cardiovascular literature describes.

**4. Direction nuance — receptor is early-high, not late-high.**
Endothelium, PP vs NN (BH over a 16-gene IL-1/adhesion set): IL1R1 is *down* in lesional endothelium
(log2FC −1.12, BH q≈0), tracking a tier drop (51%→29%→27% expressing across NN→PN→PP). IL6
(−1.39), CCL2 (−0.64), ICAM1 (−0.44), STAT3 (−0.29), NFKB1 (−0.44) are also down in lesional vs
normal endothelium; PYCARD (+1.08, q=0.009) and IL1RN (+4.10, q=0.017) are up. Interpretation: the
vascular IL-1 response is an **early/constitutive** feature that is consumed or desensitized (classic
receptor downregulation under sustained ligand) as lesions establish — consistent with IL-1 as an
*initiating* rather than *maintaining* signal.

## Literature support

- **Psoriasis→atherosclerosis via endothelium** is well established: the psoriatic vascular
  endothelium adopts a pro-inflammatory, proatherogenic phenotype (IL-1, IL-6, TNF-α, VCAM1, ICAM1,
  E-selectin), unified clinically as the "psoriatic march" linking skin inflammation to
  cardiovascular disease; psoriasis carries ~50% increased CV risk. (Frontiers Med 2022;864185;
  PMC9744099; PMC5867651)
- **IL-1 as an early signal / IL-1↔STAT3**: IL-1β is elevated in lesional skin and IL-1R correlates
  with progression and treatment response (PMC6392027); IL-1β stimulates keratinocytes to make CCL20
  via STAT3+NF-κB; an S1PR3–Src–STAT3 axis drives early and prolonged STAT3 activation
  (Nature s41419-025-07358-w); IL-36 controls early IL-23 (PMC7190273).

The **cell-level receptor→program wiring** shown here (IL1R1⁺ endothelium = the STAT3/inflammasome-
active, gradient-tracking population) is a specific, less-explored angle within this established theme.

## Blood arm (the "two-part story": skin + circulation) — feasibility

The theory's IL-1β *source* is systemic/circulating, and the atherosclerosis link is measured in
blood — so a circulation arm is the right addition. mLOY (theory 2) is also a blood-cell phenomenon.
**recount3 psoriasis blood/PBMC feasibility:**

| Study | Tissue | Subtype | Groups | n | Verdict |
|---|---|---|---|--:|---|
| SRP173379 | Neutrophils (sorted) | GPP | 8 GPP / 11 healthy | 19 | Marginal — sorted, not PBMC |
| SRP173378 | Whole blood | GPP | 9 GPP / 7 healthy | 16 | Marginal — mixed blood |
| SRP132160 | PBMC | GPP | GPP only, acitretin-treated | 15 | Not usable alone — no healthy ctrl, drug-treated |

**Honest limit + scope of what was screened:** the three psoriasis blood accessions checked here
(SRP173379/SRP173378/SRP132160 — pre-selected from prior domain knowledge, *not* an exhaustive search)
are **all Generalized Pustular Psoriasis (GPP)** — a distinct, rarer subtype from the plaque psoriasis
(PsV) of the skin arm — and all small. This is *not* a comprehensive recount3 census: the full 8,677-
project SRA table has no disease/tissue columns to filter on, so a definitive "no plaque-psoriasis PBMC
bulk exists in recount3" claim would require a per-project title/metadata scan that has **not** been run.
What is established: the known/obvious psoriasis blood studies are GPP, so recount3 does not *readily*
supply a matched plaque-psoriasis PBMC bulk arm. For a proper two-part story the blood data most likely
needs to come from GEO/ArrayExpress plaque-psoriasis PBMC datasets (to be screened). The GPP studies
remain useful as a *pustular-subtype* sensitivity check. **Open task:** run a title-keyword scan across
the recount3 SRA table before concluding recount3 has no usable plaque-psoriasis blood study.

## Theory 2 (loss-of-Y, male severity) — feasibility in current data

Sex is cleanly inferable per donor in the skin single-cell reference (Y-genes RPS4Y1/DDX3Y/UTY/
EIF1AY/KDM5D/USP9Y strongly expressed in ~20 donors, XIST in the other ~10). So **sex-stratified**
analysis of the skin data is feasible now. **mLOY proper** (hematopoietic Y-dropout vs diploid
baseline) requires *blood* cells and cannot be established from skin biopsy — it needs the PBMC arm.
Literature grounding is strong: mLOY is robustly linked to cardiovascular disease, cardiac fibrosis
(profibrotic macrophages) and coronary atherosclerosis in men (Nature Rev Cardiol 2023; Science
abn3100; SCAPIS 2025).

## Status / next steps
- Theory-1 endothelial finding: **computed, rigorous, figure saved** (fig_theory1_endothelial_IL1.png).
- Blood arm: recount3 = GPP-only; **screen GEO/ArrayExpress for plaque-psoriasis PBMC** next.
- Theory 2: sex-stratified skin pass feasible now; mLOY needs blood data.
