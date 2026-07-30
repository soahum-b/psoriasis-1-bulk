# Lead target selection — rationale

## Selection logic
From 2,379 significant meta-analysis DEGs (PP vs NN, FDR<0.05 & |log2FC|>1), 154 protein-coding
candidates were annotated for druggability. A composite score combined DE evidence (effect size,
significance), pathway centrality, and druggability (Open Targets small-molecule tractability bucket,
ChEMBL potent-inhibitor count, drug-candidate count, structural availability). Leads were then required
to sit on a psoriasis-relevant signalling axis AND be structurally workable for docking/MD.

The tractable, on-axis pool is dominated by the **JAK-STAT / Th17 axis** — the core psoriasis
inflammatory circuit: JAK3 -> STAT3 -> RORgt(RORC) -> IL-17. Three leads were carried forward:

## Selected leads

### STAT3 (designated lead)  — Transcription factor, JAK-STAT axis
- Meta-analysis: **up, log2FC +1.25, FDR 2.2e-7, concordant across all 3 studies**.
- Druggability: druggable family + high-quality ligand + ligand-bound structures; 112 potent
  (pChEMBL>=7) binders in ChEMBL, but **no approved/advanced small-molecule** — the classic
  "hard but real" TF target. Composite druggability rank 43/154 (lower by design — reflects TF difficulty).
- Structural handle: the **SH2 domain** (phosphotyrosine pocket) is the established small-molecule
  site (e.g. TTI-101/C188-9, napabucasin, stattic chemotypes).
- Why lead: central psoriasis node with prior isoform/splicing work in this project; the SH2 pocket
  is the docking/MD target.

### RORC / RORgt  — Nuclear receptor / transcription factor, Th17 master regulator
- Meta-analysis: **down, log2FC -2.47, FDR 9.3e-227** (very strong, 2 studies).
- Druggability: **best-in-class** among the leads — a genuine ligand-binding domain (LBD) pocket,
  9,766 potent inhibitors, Phase-1 small-molecule precedence. Composite druggability rank 7/154.
- Why lead: master Th17 transcription factor, directly upstream of IL-17; the most tractable TF on
  the axis and an ideal positive-control pocket for the docking pipeline.

### JAK3  — Tyrosine kinase, JAK-STAT axis
- Meta-analysis: **up, log2FC +1.55, FDR (highly significant), 3 studies**.
- Druggability: **approved-drug tractability** (tofacitinib/JAK-inhibitor class), 6,069 potent
  inhibitors, well-defined ATP pocket. Composite druggability rank 3/154.
- Why lead: the upstream kinase of the STAT3 module; a validated, structurally excellent kinase pocket
  that anchors the pipeline with a clinically de-risked target.

## Note on the highest raw-score targets
NOS2, MMP9/12, IDO1, PLK1/AURKA also score highly and are on-axis (innate/protease/immunometabolic/
proliferation). They are retained in the ranking table as secondary candidates but were not carried into
docking/MD to keep the structural work focused on the mechanistically-linked JAK-STAT/Th17 core.
