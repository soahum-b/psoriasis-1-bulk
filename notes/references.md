# References (durable library)

Per paper: full citation + DOI, why it is in the project, and the SPECIFIC claim/
parameter/method taken from it.

## Papers

### Sun, Guan, Moran, … Xia 2022 — Scissor
- DOI: 10.1038/s41587-021-01091-3 (*Nature Biotechnology* 40, 527–538; online 11 Nov 2021)
- Why: the core method — identifying phenotype-associated single-cell subpopulations by integrating bulk + single-cell.
- Used: the network-regularized elastic-net objective (reimplemented on glmnet); Gaussian/continuous phenotype mode for the ordinal gradient. Citation verified vs CrossRef.

### Ma et al. 2023 — GSE173706 single-cell psoriasis
- DOI: (GEO accession GSE173706)
- Why: the single-cell reference for the Scissor arm.
- Used: 10x droplet scRNA-seq of NN/PN/PP skin (96,088 cells); note — 10x short-read, so cannot resolve STAT3 isoforms per cell.

### Tsoi et al. 2019 — SRP165679 bulk RNA-seq
- DOI: (recount3 project SRP165679)
- Why: the tier-balanced bulk phenotype anchor (Scissor arm) + a k=4 meta-analysis study (bulk arm).
- Used: NN=38 / PN=27 / PP=28 biopsies; the only recount3 study carrying BOTH healthy and uninvolved controls (anchors the control-type moderator).

### recount3 (Wilks et al. 2021)
- DOI: 10.1186/s13059-021-02533-6
- Why: uniform reprocessing of all bulk SRA studies used.
- Used: `create_rse_manual` gene-level counts (gencode v26) for all meta-analysis studies.

### Psoriatic march / endothelial dysfunction → atherosclerosis (Theory 1 context)
- Endothelial Dysfunction in Psoriasis: An Updated Review — *Front. Med.* 2022;9:864185.
- Psoriasis and Cardiovascular Disease: Novel Mechanisms — PMC9744099.
- Why: the psoriasis→atherosclerosis comorbidity; the vascular endothelium adopts a pro-inflammatory, proatherogenic phenotype (IL-1, IL-6, TNF-α, VCAM1/ICAM1/E-selectin); "psoriatic march."
- Used: framing for Arm-C — the cell-level IL1R1+ endothelium→STAT3/inflammasome wiring is a specific angle within this established theme.

### IL-1 as an early/initiating signal; IL-1↔STAT3
- A Critical Role of the IL-1β-IL-1R Signaling Pathway in Psoriasis — PMC6392027 (IL-1R correlates with progression/treatment response).
- S1PR3-driven feedback sustains STAT3 activation in psoriasis keratinocytes — *Cell Death & Dis.* s41419-025-07358-w (early + prolonged STAT3).
- IL-36 signaling in keratinocytes controls early IL-23 — PMC7190273.
- Used: supports "IL-1/IL-36 as early signal" and the receptor-early-high direction finding.

### Mosaic loss of Y (mLOY) & cardiovascular disease (Theory 2 context)
- Mosaic loss of chromosome Y and cardiovascular disease — *Nat. Rev. Cardiol.* 2023; DOI 10.1038/s41569-023-00976-x.
- Hematopoietic loss of Y → cardiac fibrosis and heart failure — *Science* 2022; DOI 10.1126/science.abn3100.
- Why/Used: grounds theory 2 (mLOY worse in males); mLOY is a hematopoietic phenomenon → motivates the blood arm.

> IL-36-dominant-axis / IL-1β-blockade-disappointing / STAT3-Th17-node summary remains in
> IL1B_STAT3_expanded_analysis.md (Arm-A) and theory1_endothelial_IL1_and_blood_arm.md (Arm-C).
