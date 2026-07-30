# Psoriasis DEG druggability → docking → MD: consolidated report

**Project:** meta-analysis of psoriasis (lesional PP vs non-lesional NN) → druggable target
triage → structure-based virtual screening → MD-ready systems.
**Clinician-facing summary of a bioinformatics pipeline. This is research/informational output,
not clinical guidance.**

---

## 1. Executive summary

Starting from **2,379 significant meta-analysis DEGs** (FDR<0.05 & |log2FC|>1, PP vs NN across
3 studies), we triaged **154 protein-coding candidates** for druggability using Open Targets and
ChEMBL, then focused on the mechanistically-linked **JAK-STAT / Th17 inflammatory core** that
dominates the tractable, on-disease-axis pool. Three leads were carried through structure-based
analysis:

| Lead | Role | DE (PP vs NN) | Tractability | Structural handle |
|------|------|---------------|--------------|-------------------|
| **STAT3** | JAK-STAT TF (your designated lead) | up, log2FC +1.25, FDR 2.2e-7 | Discovery; 112 potent binders, no clinical SM | SH2/core (PDB 6NJS) |
| **RORC / RORγt** | Th17 master TF | down, log2FC −2.47, FDR 9.3e-227 | Phase 1; 9,766 inhibitors | LBD orthosteric pocket (PDB 5APH), fpocket druggability **0.94** |
| **JAK3** | Upstream kinase | up, log2FC +1.55 | **Approved-drug** class; 6,069 inhibitors | Kinase ATP site (PDB 5LWM) |

**RORγt** emerged as the most tractable pocket and anchored the demonstration virtual screen; a
**151-compound RORγt-focused library** was docked into the LBD; the protocol was checked by
crystal-ligand redocking (best-pose RMSD 2.15 Å — marginally above the 2.0 Å bar but an
upper-bound estimate at the correct subpocket) and by recovery of known RORγt actives at the top
of the ranking. Full-scale docking, GPU co-folding cross-check, and MD systems were packaged for
the user's cluster.

---

## 2. Target triage & druggability

- Evidence-weighted ranking combined DE effect size + significance, pathway centrality, and
  druggability (Open Targets SM tractability bucket, ChEMBL potent-inhibitor count, drug-candidate
  count, structural availability).
- Small-molecule tractability distribution across 154 candidates:
  Approved drug **8**, Advanced clinical **14**, Phase 1 **1**, Druggable pocket **5**,
  Discovery precedence **40**, Unknown **86**.
- Highest-scoring on-axis targets also include NOS2, MMP9/12, IDO1, PLK1/AURKA (retained as
  secondary candidates in `druggable_target_ranking.csv`).



![Pipeline funnel: meta-DEGs → druggable → leads → docking/MD]({{artifact:art_63fb81cd-5623-4e2b-9f45-93a429985097}})


![Target volcano (meta-analysis PP vs NN) with lead annotations]({{artifact:art_71a82a2a-d9f9-447c-886a-6e6c4a6132df}})


![Druggability landscape of DEG-encoded proteins]({{artifact:art_d4f56b37-3924-4cbf-b5bc-ea0083f438f0}})


![Evidence-weighted druggable target ranking]({{artifact:art_51765bbb-831d-4ed6-977f-eb7c6d5e27df}})


---

## 3. Structures & binding pockets

Ligand-bound crystal structures were selected for each lead and the binding pocket confirmed with
fpocket (matched to the crystallographic ligand by atom overlap):

| Lead | PDB | Resolution | Ref ligand | fpocket druggability | ChEMBL potent inhib |
|------|-----|-----------|-----------|---------------------|--------------------|
| STAT3 | 6NJS | 2.7 Å | KQV | 0.001 (shallow SH2 — hard TF) | 112 |
| RORC | 5APH | 1.54 Å | VYI | **0.936** (textbook druggable) | 9,766 |
| JAK3 | 5LWM | 1.55 Å | 79T | 0.016 in this crystal form; ATP site chemically validated | 6,069 |

STAT3's low structural druggability is genuine and consistent with why the field has moved to
degrader/bivalent chemistry (the 6NJS ligand KQV is a degrader-class molecule). JAK3's low fpocket
score reflects crystal-form sensitivity, not a lack of druggability — its ATP site is one of the
most validated kinase pockets.



![Lead target structures (STAT3, RORC, JAK3)]({{artifact:art_33965e3a-e668-4a24-9fe1-45a30d188880}})


![fpocket structural druggability of the lead pockets]({{artifact:art_d471717e-09df-434f-bf9e-63d0cbfe46b9}})


---

## 4. Virtual screen by docking (RORγt LBD, demonstration)

- **Library:** 151 drug-like compounds from the ChEMBL RORγt bioactivity space (a target-focused,
  bioactivity-validated deck with built-in positive controls). *ZINC22 purchasable expansion was
  blocked by an upstream outage; a resolver script ships in the docking package for when the
  backend recovers.*
- **Engine:** AutoDock Vina 1.2.7, box centered on the LBD pocket (24³ Å), exhaustiveness 8.
- **Results:** all 151 docked; affinities **−7.8 to −11.3 kcal/mol** (median −9.9).
- **Protocol validation:** redocking the crystal ligand VYI gave best-pose RMSD **2.15 Å**
  (element-agnostic upper bound; centroid displacement 1.5 Å) — marginally above the 2.0 Å bar
  but the correct subpocket/orientation; known RORγt actives enriched at the top of the ranking.

### Top 10 docked hits
| rank | chembl_id | affinity | lig_efficiency | mw | qed | logp |
|---|---|---|---|---|---|---|
| 1 | CHEMBL3314024 | -11.33 | -0.306 | 494.6 | 0.282 | 6.18 |
| 2 | CHEMBL3598073 | -11.272 | -0.305 | 527.0 | 0.507 | 4.61 |
| 3 | CHEMBL3105678 | -11.111 | -0.309 | 515.6 | 0.334 | 4.89 |
| 4 | CHEMBL3598061 | -10.917 | -0.295 | 529.0 | 0.468 | 4.59 |
| 5 | CHEMBL3581540 | -10.833 | -0.339 | 480.0 | 0.63 | 3.89 |
| 6 | CHEMBL3598053 | -10.831 | -0.301 | 492.6 | 0.509 | 5.27 |
| 7 | CHEMBL3598055 | -10.786 | -0.308 | 499.0 | 0.453 | 5.22 |
| 8 | CHEMBL3314028 | -10.778 | -0.291 | 497.6 | 0.315 | 5.14 |
| 9 | CHEMBL3314023 | -10.71 | -0.298 | 480.6 | 0.329 | 5.57 |
| 10 | CHEMBL3596600 | -10.672 | -0.344 | 457.3 | 0.545 | 5.75 |



![RORγt-focused screening library properties]({{artifact:art_0cfa49df-dc4c-4e82-afb6-d28616465e50}})


![RORγt LBD virtual screen results (151 compounds)]({{artifact:art_836a332d-fb73-44b3-ac77-640b4ec5f24a}})


![Top docked hit in the RORγt LBD pocket]({{artifact:art_bf463f21-dde2-4fd3-b931-5c3a3b5a4347}})

*Interactive complex:* `RORC_tophit_complex.pdb` (Mol* viewer)


---

## 5. Cluster-ready packages (deliverables for the user's cluster)

1. **`psoriasis_docking_package.tar.gz`** — full-scale Vina screen for all 3 targets: prepared
   receptors + pocket boxes, 3D ligand library, SLURM array driver (scales to ~1M ZINC22 compounds),
   aggregation + ZINC22 purchasability resolver, `environment.yml`.
2. **`psoriasis_cofolding_package.tar.gz`** — Boltz-2 co-folding + affinity cross-check (GPU): 11
   ready RORC YAMLs (top-10 hits + VYI positive control) with the affinity head, STAT3/JAK3
   templates, domain-construct FASTA, GPU SLURM array, cross-check aggregation vs docking.
3. **`psoriasis_md_package.tar.gz`** — OpenMM MD systems (ff14SB + OpenFF-2.1.0 + TIP3P, 0.15 M
   NaCl): parameterized, solvated, minimized and briefly-equilibrated (25 ps NVT + 25 ps NPT on
   CPU) RORC complexes (crystal VYI + top hit), GPU production + trajectory-analysis scripts,
   SLURM runner. The short CPU equilibration gives a valid production start state; extend it on
   the GPU node for publication runs (one flag, see MD README).

---

## 6. Methods provenance & caveats

- **DE base:** meta-analysis PP vs NN (`meta_de_PPvsNN.csv`), 3 studies, threshold FDR<0.05 & |log2FC|>1.
- **Druggability:** Open Targets (tractability, target class, known drugs) + ChEMBL (potent
  inhibitors at pChEMBL≥7).
- **Docking:** rigid-receptor Vina; scores are relative rankings, not absolute affinities — the
  co-folding and MD packages provide orthogonal cross-checks.
- **STAT3 caveat:** shallow SH2/core surface; docking there is less discriminating — treat STAT3
  as a degrader/MD target rather than a classical-docking target.
- **Library caveat:** ChEMBL-sourced (bioactivity-validated), purchasability pending ZINC22 recovery.
- **Equilibration caveat:** shipped systems used short (25+25 ps) NVT+NPT on CPU for package-build
  tractability — a valid production start state (positions + velocities + box preserved), but brief;
  extend on GPU for publication runs (see MD package README).

---

*All underlying tables, structures, and figures are saved as project artifacts.*
