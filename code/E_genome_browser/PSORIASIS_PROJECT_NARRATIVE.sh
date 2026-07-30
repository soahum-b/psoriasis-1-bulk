# Auto-extracted generating script
# Produces: PSORIASIS_PROJECT_NARRATIVE.md
# Conda env: bash   (run with this environment activated)
# Inputs (expected alongside / in data/): stat3_gradient_isoform_vs_gene.csv, key_genes_7study_PPvsNN.csv, sample_manifest.csv, trend_SRP165679.rds, lead_pipeline_summary.csv, docking_ranked.csv
# Source artifact version: be496398-cd79-4c36-a10e-5f568571b521
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

content = """# Psoriasis, from expression to target: the story of the project

*A prose narrative of the whole program — every analytical arm, what it asked, what it found, and how the pieces fit together. Written for a reader who wants the argument, not the file listing. This is target-discovery and hypothesis-generating research on public data; it makes no clinical, diagnostic, or causal claims about individual patients.*

---

## The question we set out to answer

The project began with a clinician's question and a favored suspect. The question: in psoriasis, which genes are deranged in lesional skin, which biological pathways do they assemble into, and — the part that matters for therapy — which of the proteins in those pathways are actually druggable? The favored suspect was **STAT3**, a transcription factor sitting at the crossroads of the IL-6/JAK-STAT and IL-17 signaling that dermatology already associates with plaque psoriasis. STAT3 was the lead hypothesis, but the project was deliberately open to letting the data nominate other targets.

What follows is the arc the data actually drew. It is more interesting than a single-gene confirmation, because at several points the evidence complicated the STAT3 story rather than simply endorsing it — and the discipline of following those complications is where the scientific value lies.

---

## The material and how it was made comparable

Two kinds of public data anchor the work. For the bulk transcriptomics, we used **recount3** — a resource that reprocesses tens of thousands of SRA RNA-seq runs through one uniform pipeline against one annotation (Gencode v26, GRCh38). The reason this matters is not convenience. Meta-analysis of expression data has historically been undermined by having to reconcile different microarray platforms, different probe-to-gene mappings, different quantification choices; recount3 removes that entire class of problem by making the gene space identical across studies by construction. For the single-cell arm, we used one deeply annotated reference — GSE173706 (Ma et al., 2023) — comprising 33 count matrices from 22 donors, which after quality control yielded roughly 89,000 cells annotated to nine skin lineages.

Across the whole project, biopsy site is treated as an ordered series rather than a binary: **normal skin (NN)**, **peri-lesional skin (PN)** — the clinically uninvolved margin next to a plaque — and **lesional skin (PP)**, the plaque itself. Reading these three as a progression (NN → PN → PP) rather than just contrasting plaque against normal turned out to be one of the more productive decisions in the project, because it exposes *when* along the march to disease each program switches on.

---

## Arm A — The bulk meta-analysis, and the verdict on STAT3

The backbone of the project is a random-effects meta-analysis of lesional-versus-normal skin. We screened 23 candidate human psoriasis studies in recount3 against sample-level metadata and kept the four that were whole-skin biopsies containing both lesional and normal tissue, were treatment-free, and had adequate sequencing depth: SRP035988, SRP165679, SRP126422, and SRP065812, together contributing 141 normal and 145 lesional samples. One large study (ERP110816) was deliberately excluded despite its size because every patient was on etanercept — an anti-TNF drug that pharmacologically suppresses the very TNF → IL-1 → STAT3 axis under study — so including it would have been like measuring a fire after the sprinklers came on.

Each study was analyzed on its own first (filter low-expressed genes, TMM-normalize, then limma-voom to get a fold-change and its standard error), and only then were the four combined. This "model-then-merge" order is what keeps a meta-analysis honest: each study retains its own normalization and its own noise structure, and the disagreement between studies stays visible rather than being averaged away silently. The pooling itself used the closed-form DerSimonian-Laird random-effects estimator, vectorized across all ~24,000 genes and checked against the reference `metafor` implementation. Because four studies is a small number, the primary confidence intervals used the Hartung-Knapp-Sidik-Jonkman correction, which is honest about uncertainty when the study count is low; the simpler estimator was kept as a sensitivity check.

The verdict on the lead hypothesis was affirmative and robust. **STAT3 is reproducibly up in lesional skin** — a pooled log2 fold-change of roughly +1.06 to +1.25 depending on the study set, and critically, this holds up under leave-one-out testing including removal of the single largest ("anchor") study. The direction is not an artifact of one dominant cohort. STAT3 sits inside the tightest, most reproducible signal in the whole panel: the **IL-17 / interferon / IL-6–JAK-STAT3** inflammatory core, with STAT1 (+1.70), JAK3 (+1.12 in the four-study meta), and TYK2 (+0.23) all up, and JAK1 slightly *down* (−0.23) in a pattern consistent with a constitutive rather than induced role.

Two negative results from this arm are worth as much as the positive one, because they redirect the target thinking:

First, **the IL-1β ligand gene (IL1B) is the least reliable member of the panel** — nominally up (+1.09) but not significant, with heterogeneity so high (I² = 94%) that the fourth cohort actually pulls it downward. What *is* consistently up is the machinery around it — the inflammasome components CASP1, PYCARD, IL1RN, AIM2, CASP5. The lesson: any target hypothesis resting on IL-1β should rest on the assembled inflammasome and its downstream output, not on IL1B transcript abundance, which the data will not support.

Second, **IL-36 dwarfs IL-1β** within the IL-1 family (IL36A +9.0, IL36G +4.7). And the single strongest signal in the entire analysis is not a signaling gene at all but the shared downstream *output* — the antimicrobial and chemokine cassette S100A7/8/9, DEFB4A, LCN2, CCL20 — the common readout that every upstream arm converges on.

A methodological point that recurs in the target discussion: adding the fourth study shrank several inflated fold-changes toward more honest pooled values and exposed (rather than created) genuine between-study disagreement in a few genes like SOCS3 and IL18. That is the meta-analysis working as intended.

---

## The STAT3 isoform question — tested, and retired

STAT3 exists as two functionally opposite isoforms: the full-length **STAT3α**, which carries the C-terminal transactivation domain and drives transcription, and a shorter **STAT3β**, which lacks that domain and acts as a brake. An attractive hypothesis was that psoriasis involves a *shift in the ratio* — that the disease tips the balance toward the activating isoform. We tested this carefully, computing the β-isoform fraction (percent-spliced-in) from junction-level data across cohorts and comparing seven different statistical estimators.

The hypothesis did not survive. The apparent β-switch was real only in the single anchor cohort (p = 0.017 there) and **vanished the moment it was subjected to the honest tests** — it disappeared under the few-study-corrected interval, under both mega-analysis variants, and completely under anchor-drop (the combined p-value went from 0.045 to 0.91). The conclusion, stated plainly in the project record, is that the psoriasis STAT3 signal is **gene-level up-regulation of the full-length activator, not a redistribution between isoforms**. This is near the ceiling of what short-read bulk RNA-seq can resolve about splicing, and the result was correctly retired rather than published as a finding.

*(This is precisely the tension the genome browser built in the most recent session makes visual: the STAT3 gene climbs monotonically NN → PN → PP, log2 7.47 → 7.63 → 9.15, while the STAT3β isoform fraction stays essentially flat at ~8–9%. The picture and the statistics agree.)*

---

## Arm B — The staging axis: when each program switches on

Because one cohort (SRP165679) contains all three biopsy tiers with good depth, we could treat NN → PN → PP as an ordered trajectory and ask, for each gene, *how far along the march to lesional the peri-lesional margin has already traveled.* We fit a linear trend and a monotonicity test per gene and classified each by the fraction of its total lesional change already reached at the peri-lesional stage.

The findings turn a static contrast into a story about timing. **About 85% of the lesional program is strictly monotonic** across the three stages, and the peri-lesional margin sits, on average, about 16% of the way to full lesional — meaning clinically "uninvolved" skin is already an early molecular stage of disease, not a clean baseline. The programs do not all switch on together: **inflammation and interferon signaling are early** (already elevated in the peri-lesional margin), while **keratinocyte proliferation is late** (essentially flat at the margin, switching on only at the plaque). And STAT3 specifically shows **activity before transcription** — its regulon activity reaches ~30% of the lesional level at the margin while its transcript reaches only ~10% — placing STAT3 *activation* early in the sequence, upstream of its own transcriptional induction.

---

## Arm C — Single cells, and the moment STAT3 stopped being the whole answer

The bulk arms establish that STAT3 is up and that its activity comes on early. The single-cell arm asked a different and sharper question: *which cells actually drive the progression?* Using Scissor — a method that projects a bulk phenotype onto individual cells — we mapped the ordinal biopsy-site label (NN < PN < PP) onto all ~89,000 cells. This design is non-circular by construction, because the phenotype being projected is the *clinical* label, not any molecular score derived from the cells themselves. The underlying solver is a pure-R re-implementation of Scissor's network-regularized elastic net, and the result was guarded with two independent permutation nulls (a reliability test over 100 label shuffles and a selection null over 30 reruns, both p ≈ 0).

The result reframed the project. The cells that track the progression are overwhelmingly **endothelial** — the vasculature — with an odds ratio of 5.24 at full census and a vascular signature (CCL14, ACKR1, RAMP3, PLVAP). And **STAT3 did not survive as a marker of the gradient** at full census (BH q = 0.13). Notably, STAT3 *had* looked significant on a smaller 20,000-cell backbone, and the honest move was made: the effect was downgraded when the full census showed it to be a small-subset artifact.

This is not a contradiction of Arm A — it is the reconciliation. STAT3 is expressed nearly everywhere in the tissue and is genuinely up in bulk, but at single-cell resolution it is a **passenger in a vascular-led program** rather than the cell-type-specific engine of progression. An orthogonal deconvolution (NNLS of the real bulk against the single-cell reference) confirmed the endothelial signal is a **real compositional increase** — endothelial proportion rises 0 → 0.2 → 3.8% across the tiers — and not merely a change in cell state.

---

## Arm D — Theory-1: an IL-1-primed vascular circuit

Having found a vascular-led program, we asked what wires it. The answer assembled into a specific, literature-grounded model. **IL-1β is not made by the endothelium** — it is largely restricted to dendritic cells (expressed in 48% of DCs, only 0.2% of endothelial cells), which makes the signaling **paracrine**: the DC speaks, the vessel listens. And the vessel is built to listen — the IL-1 receptor IL1R1 is carried by ~29% of endothelial cells. The IL1R1⁺ endothelial cells are **2.3× enriched** among the progression-tracking cells (odds ratio 2.32, p ≈ 10⁻³⁸) and co-express markedly more STAT3, IL6, NFKB1, CASP1, GSDMD, and PYCARD — the receptor-to-program wiring made visible at single-cell resolution. The receptor is early and constitutive (it actually falls in fully lesional endothelium), consistent with IL-1 acting as an *initiating* signal. This arm connects the skin phenotype directly to the well-known **psoriatic-march** cardiovascular comorbidity, in which skin inflammation tracks with vascular disease.

This arm is framed, correctly, as exploratory but grounded — a specific hypothesis within an established theme, not an overclaim.

---

## The threads that were deliberately dropped

Part of the project's credibility is what it declined to claim. A **blood bulk arm** — comparing psoriasis blood against healthy blood using ERP110814 versus GTEx — was dropped because disease status was perfectly confounded with study and globin-prep batch: no statistical adjustment can rescue a confounder that is 100% collinear with the exposure. A **loss-of-Y-chromosome (mLOY)** thread, motivated by a male-severity angle, gave a clean negative in skin: the fraction of Y-silent cells was driven by sequencing depth, not disease tier, and genuine genomic mLOY would need blood or DNA rather than skin RNA. Both threads remain open only in the sense that they would require a dataset that recount3 does not contain — a single study with both plaque-psoriasis and healthy blood processed together.

---

## From pathways to proteins — the target and druggability layer

The pathway biology converged beautifully, but converging on "these pathways are up" is not the same as naming ranked, tractable drug targets. That final half-step was built out. Starting from **2,379 significant meta-analysis DEGs** (FDR < 0.05 and |log2FC| > 1), we triaged **154 protein-coding candidates** for druggability using Open Targets (genetic association, tractability bucket, known drugs) and ChEMBL (counts of potent inhibitors). The tractability distribution across those 154 is itself informative: 8 are targets of an approved drug, 14 are in advanced clinical development, and the large tail is discovery-stage or unknown.

Three leads were carried into structure-based analysis, and the choice among them makes an important point about **abundance versus druggability**:

- **STAT3** — up (+1.25), the designated lead, with 112 potent binders known but no clinical small molecule. Its structural problem is real: the SH2/core surface is shallow (fpocket druggability 0.001), which is exactly why the field has moved to degrader and bivalent chemistry for STAT3. It is best treated as a **degrader / molecular-dynamics target**, not a classical docking target.
- **RORC / RORγt** — the Th17 master transcription factor, strongly *down* in the data (−2.47) yet by far the most tractable pocket: a textbook druggable ligand-binding domain (fpocket 0.94) with nearly 10,000 known inhibitors. Being down-regulated does not disqualify it as a target — its pathway position and pocket quality do the arguing.
- **JAK3** — up (+1.54 in the three-study base set on which the druggability layer was built; +1.12 in the four-study meta), in the **approved-drug** class (tofacitinib), a validated kinase ATP site.

The direction-of-change caveat deserves emphasis because it is where naïve target-picking goes wrong: **TYK2 is essentially flat at the transcript level yet is the target of an approved psoriasis drug (deucravacitinib), and JAK1 is actually down.** The target argument has to be made on pathway position plus tractability, not on fold-change alone.

Because RORγt offered the best pocket, it anchored a demonstration virtual screen: a 151-compound RORγt-focused library (drawn from the ChEMBL bioactivity space, with built-in positive controls) was docked into the ligand-binding domain with AutoDock Vina, yielding affinities from −7.8 to −11.3 kcal/mol. The protocol was validated by redocking the crystal ligand (best-pose RMSD 2.15 Å — marginally above the strict 2.0 Å bar, but recovering the correct subpocket and orientation) and by the enrichment of known actives at the top of the ranking. Full-scale docking, a Boltz-2 GPU co-folding cross-check, and OpenMM molecular-dynamics systems were packaged for the cluster as ready-to-run deliverables.

---

## The genome browser — the modules made visual

Most recently, the project's gene-level findings were rendered as a **STAT3-led genome browser**: locus-track figures comparing the NN → PN → PP groups across the four focus modules — JAK-STAT (STAT3, STAT1, JAK1, JAK3, TYK2, SOCS3), the IL-1/inflammasome set, IL-36, and the Th17/S100 output — each gene shown as its Ensembl transcript model paired with its per-group expression change. The STAT3 headline figure pairs the α/β transcript models with a C-terminus zoom (on this minus-strand gene the two isoforms are identical across the body and diverge only at the 3′ coding end: β terminates at an alternative exon near chr17:42,316,830, while α splices past it and extends ~1.1 kb further to its terminal coding exon near 42,315,770 — the extra coding that builds the transactivation domain β lacks) and the quantitative panels showing gene expression rising while the β-isoform fraction stays flat. The browser makes the biological directionality legible at a glance, including the informative exceptions — JAK1 down, IL18 down, and RORC falling sharply even as its IL17A output peaks.

---

## What the project actually concluded

The honest one-sentence synthesis is not "STAT3 is the target." It is this: the lesional psoriasis program is a reproducible **IL-17 / interferon / IL-6–JAK-STAT3 inflammatory core with a dominant S100/defensin/chemokine output**, it comes on in a definite order (inflammation early, proliferation late), and at single-cell resolution it is **led by an IL-1-primed IL1R1⁺ endothelium** in which STAT3 and the inflammasome run downstream. STAT3 is a genuine, reproducible hub — but the integrated target is better described as the **IL-1-primed vascular circuit** than as any single gene. On the druggability axis, RORγt offers the cleanest pocket, JAK3 the most clinically validated one, and STAT3 a harder problem suited to degrader chemistry.

---

## Limitations, stated plainly

Every expression conclusion rests on mRNA, not protein; the staging arm's activity-before-transcription result is the closest we come to an activation-state readout, and a protein-level (pSTAT3) cross-check remains a genuine gap. The staging axis is well-powered in one cohort and supported by pooled trends in the others rather than independently replicated at full power. The single-cell arm rests on one reference dataset. The docking scores are relative rankings from a rigid-receptor method, not absolute affinities — which is why the orthogonal co-folding and MD packages exist. And the between-study heterogeneity is real and pervasive, appropriately handled by the random-effects model but a caveat on any single pooled magnitude.

---

## The threads still open

Four next steps are queued and feasible. A **drug-first cross-check** (CMap/LINCS) would test whether the psoriasis meta-signature is *reversed* by JAK/TYK2 inhibitors, nominating targets from perturbation data independently of expression. A **protein/pSTAT3 layer** would move the central claim off mRNA. A **sex-stratified re-analysis** of the skin arm is possible right now with existing data (sex is already called per donor, 21 male / 12 female) and connects to the male-severity thread. And the **blood / mLOY threads** stay open until a dataset pairing psoriasis and healthy blood in one batch becomes available — something recount3 does not currently provide.

---

*Every quantitative value above is drawn from the project's saved analysis artifacts and synthesis documents; interpretation is preliminary and for research use, not clinical guidance.*
"""

with open("PSORIASIS_PROJECT_NARRATIVE.md", "w") as f:
    f.write(content)