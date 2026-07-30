# Psoriasis focus-gene genome browser (NN → PN → PP)

Static locus-track figures comparing the STAT3-led focus modules across the
ordinal biopsy gradient: **NN** (normal) < **PN** (peri-lesional) < **PP** (lesional).

## Data provenance
- **Gene models / exon structure:** Ensembl v86 (GRCh38) via `EnsDb.Hsapiens.v86`, rendered with Gviz.
- **Per-group expression:** mean log2 expression per arm from `trend_SRP165679.rds`
  — the SRP165679 gradient cohort (8 NN, 11 PN, 14 PP; the one study carrying all three tiers).
- **STAT3 β-isoform fraction (PSI):** `stat3_gradient_isoform_vs_gene.csv`.
- Bars show expression **change vs NN (log2)** so direction is comparable across genes with different baselines.

## Figures
- **stat3_headline.png** — STAT3 four-panel lead: (a) full locus α/β models, (b) exon-23 divergence zoom,
  (c) gene expression NN→PN→PP, (d) β-isoform PSI per group.
- **stat3_lead.png / stat3_exon23_zoom.png / stat3_quant_panel.png** — the individual STAT3 panels.
- **module_jakstat.png** — STAT3, STAT1, JAK1, JAK3, TYK2, SOCS3.
- **module_il1_inflammasome.png** — IL1B, CASP1, CASP5, PYCARD, AIM2, NLRP3, GSDMD, IL18.
- **module_il36.png** — IL36A, IL36G, IL36RN.
- **module_th17_output.png** — IL17A, IL23A, RORC, S100A7/8/9, DEFB4A, CCL20.

## STAT3 α vs β (the isoform story)
- STAT3α = ENST00000264657 (canonical, full transactivation domain).
- STAT3β = ENST00000585517.
- On this minus-strand gene the two isoforms are identical across the body and diverge only at
  the **3′ coding (C-terminal) end**: β terminates at an alternative coding exon near
  chr17:42,316,830, while α splices past it and extends ~1.1 kb further to its terminal coding
  exon near 42,315,770. That extra α coding builds the C-terminal transactivation domain; β is
  the dominant-negative/regulatory isoform lacking it. (An earlier draft mislocated this to a
  19-bp exon-23 difference at the 42,348k end — that is the 5′/N-terminal end and was corrected.)
- **Key result:** STAT3 *gene* expression rises monotonically with lesion stage
  (log2 7.47 → 7.63 → 9.15), but the STAT3β *isoform fraction* does not (PSI 8.85% → 7.69% → 8.78%).
  The disease signal is increased transcription of the full-length activator, not an isoform switch.

## Biological directionality (per-group trend)
- **Up NN→PP:** STAT3, STAT1, JAK3, TYK2, SOCS3 (JAK-STAT); IL1B, CASP1/5, PYCARD, AIM2, NLRP3, GSDMD
  (inflammasome); IL36A/G/RN; IL17A, IL23A, S100A7/8/9, DEFB4A, CCL20 (Th17/output).
- **Down / non-monotonic:** JAK1 (constitutive), IL18, RORC (falls as its IL17A output peaks).

## Notes
- Ideogram track omitted (UCSC cytoband parsing failed under Gviz); a labeled genome axis is used instead.
- IL17F is absent from this cohort (fails expression filter; k=1 in the meta-analysis).
- Signal is gene/exon-level from existing counts, not per-base BigWig coverage.
