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

> Add IL-1β / IL-36 / STAT3 psoriasis literature entries here as the Arm-A writeup is finalised
> (IL-36 as dominant IL-1-family axis; IL-1β blockade disappointing in plaque psoriasis; STAT3 as
> Th17 convergence node — currently summarised in IL1B_STAT3_expanded_analysis.md).
