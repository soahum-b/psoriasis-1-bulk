# Resume: STAT3-led genome browser (NN→PN→PP)

**Session goal:** Build (1) static locus-track figures + (3) exon/junction-level STAT3 comparison,
leading with STAT3 (α vs β), comparing NN → PN → PP groups. Running as GENOMICS_ANALYST profile.

## Status: setup complete, first figure not yet rendered. Stopped for VPN/cluster issue.

## Decisions locked
- Track data source: **gene/exon-level from existing project counts** (user chose this), runs LOCALLY.
  No cluster needed for the browser. Cluster (ssh:n013) reserved for optional per-base BigWig later.
- Env: **psoriasis-r** — Gviz + EnsDb.Hsapiens.v86 + ensembldb + GenomicFeatures + org.Hs.eg.db INSTALLED & working.
- UCSC ideogram domain (hgdownload.soe.ucsc.edu) granted for IdeogramTrack.

## Data feeding the tracks (all artifacts, current project)
- **trend_SRP165679.rds** (vid d61263b9-9355-4e92-97d8-4987bcb38e86): per-group means mu_NN/mu_PN/mu_PP
  for all 24,640 genes + monotonicity flags. THE per-group signal source. 33-sample gradient cohort
  (8 NN, 11 PN, 14 PP; SRP165679).
- **key_genes_7study_PPvsNN.csv** (5fd6338b-66d8-45e2-bcc4-6f931f2b358d): pooled meta log2FC/FDR/I2/k.
- **stat3_gradient_isoform_vs_gene.csv** (087461f6-2a17-486a-bb9b-4f76963e8035): PSI-β per arm.
- sample_manifest.csv (a8220646-c817-4350-bfe9-f46609d614d3): NN/PN/PP labels.

## STAT3 structure (GRCh38, EnsDb v86) — WORKED OUT
- Gene: chr17:42,313,324-42,388,568, minus strand, 15 transcripts.
- STAT3α = ENST00000264657 (canonical, full TAD).
- STAT3β = ENST00000585517.
- α/β diverge at **exon 23**: α uses ENSE00001302308 (ends 42,348,539), β uses ENSE00002810203
  (ends 42,348,520) — 19-bp alternative splice → β truncated C-terminus (loses transactivation domain).
- Isoform-zoom window: ~chr17:42,348,300-42,348,600.

## Focus genes (per-group means confirmed, all monotonic NN→PP unless noted)
- JAK-STAT: STAT3(7.47→7.63→9.15), STAT1, JAK3, SOCS1/3 up; JAK1, JAK2 DOWN.
- IL-1/inflammasome: AIM2, CASP5, IL1B, GSDMD, PYCARD, NLRP3, IL1RN, CASP1 up; IL18, IL1A down.
- IL-36: IL36A, IL36G, IL36RN, IL36B up; IL1F10 down.
- Th17/output: DEFB4A, S100A8/9/7, PI3, LCN2, CCL20, IL17A, IL23A up; RORC DOWN. (IL17F absent from cohort.)
- KEY ISOFORM POINT: STAT3 gene rises monotonically but PSI-β fraction does NOT (NN 8.85% → PN 7.69% → PP 8.78%).

## NEXT STEPS (tomorrow)
1. Fix GeneRegionTrack: restrict to STAT3 window (grt built OK, 223 features/17 tx) — but I called
   feature(grt) which errored (no such slot). Use Gviz plotTracks directly; don't call feature().
2. Render STAT3 lead panel: Ideogram + GenomeAxis + GeneRegionTrack(α/β labeled) + per-group
   expression track (mu_NN/mu_PN/mu_PP as DataTrack barchart) + exon-23 isoform zoom inset.
3. Render focus-gene locus panels per module (transcript model + 3-group expression track).
4. Exon/junction STAT3 comparison figure (PSI-β per group).
5. save_artifacts all figures (PNG/PDF), embed inline.
