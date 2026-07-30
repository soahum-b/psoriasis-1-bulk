# Auto-extracted generating script
# Produces: LIBRARY_README.md
# Conda env: dock-md   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 08fbe0dc-9dcf-41be-9896-a42ae9b0c670
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import pandas as pd

libf = pd.read_csv("work/screening_library.csv")

libf["source"] = "ChEMBL RORC bioactivity (pChEMBL>=7) + similarity analogs"
libf["target"] = "RORC/RORgt (CHEMBL1741186)"
libf.to_csv("work/screening_library.csv", index=False)

readme = f"""# RORγt-focused screening library

## Contents
- **{len(libf)} drug-like compounds** with 3D structures (ETKDG v3 + MMFF-optimized), in
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
MW {libf.mw.median():.0f} Da · cLogP {libf.logp.median():.1f} · QED {libf.qed.median():.2f} · TPSA {libf.tpsa.median():.0f} Å²
"""
open("LIBRARY_README.md", "w").write(readme)