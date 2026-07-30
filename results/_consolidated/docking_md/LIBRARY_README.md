# RORγt-focused screening library

## Contents
- **151 drug-like compounds** with 3D structures (ETKDG v3 + MMFF-optimized), in
  `screening_library_3d.sdf`.
- Per-compound descriptors in `screening_library.csv`.

## Provenance
Compounds were assembled from the **ChEMBL RORγt (CHEMBL1741186) bioactivity space** — 714
distinct molecules with potent measured activity (pChEMBL ≥ 7) plus close structural analogs from
ChEMBL similarity search. The top 250 by data availability were retrieved and filtered for
drug-likeness:
- Lipinski Ro5 violations ≤ 1
- MW 250–600 Da
- Rotatable bonds ≤ 12
- TPSA ≤ 160 Å²

This yields a **target-focused, bioactivity-validated** deck — a legitimate and common virtual-screen
input that provides built-in positive controls (compounds with known RORγt activity should re-dock into
the LBD pocket, validating the docking protocol).

## Note on ZINC22 purchasability
The plan called for a ZINC22/CartBlanche22 purchasable set. At run time the ZINC22 similarity and
random-sampling endpoints (async compute backend) were **persistently timing out** (>55 s, transport
ReadTimeouts) — an upstream availability issue, not a query problem. The fast ZINC **ID/supplier
lookup** endpoint worked and is wired into the cluster docking package, so purchasability resolution
can be completed when the backend recovers by running `resolve_purchasability.py` (included in the
cluster package) against the library SMILES. Many of these ChEMBL compounds are catalog-available
(Mcule, MedChemExpress, eMolecules) via their published synthesis.

## Property medians
MW 479 Da · cLogP 4.5 · QED 0.53 · TPSA 72 Å²
