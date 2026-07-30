# Auto-extracted generating script
# Produces: screening_library.csv
# Conda env: dock-md   (run with this environment activated)
# Inputs (expected alongside / in data/): rorc_actives_smiles.json, chembl_simlib.json
# Source artifact version: 6a14ccbf-ab06-4b88-9ef6-8bfbeece98a9
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import json
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, QED, Crippen, AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize

recs = json.load(open("rorc_actives_smiles.json"))
simlib = json.load(open("chembl_simlib.json"))

# merge, dedupe by chembl_id
allc = {r["chembl_id"]: r for r in recs}
for s in simlib:
    allc.setdefault(s["chembl_id"], s)
print("combined unique compounds:", len(allc))

# Standardize + compute descriptors + drug-likeness filter
un = rdMolStandardize.Uncharger()
lfrag = rdMolStandardize.LargestFragmentChooser()
rows = []
for cid, r in allc.items():
    m = Chem.MolFromSmiles(r["smiles"])
    if m is None:
        continue
    m = lfrag.choose(m)
    m = un.uncharge(m)
    mw = Descriptors.MolWt(m)
    logp = Crippen.MolLogP(m)
    hbd = Descriptors.NumHDonors(m)
    hba = Descriptors.NumHAcceptors(m)
    tpsa = Descriptors.TPSA(m)
    rb = Descriptors.NumRotatableBonds(m)
    qed = QED.qed(m)
    ar = Descriptors.NumAromaticRings(m)
    ro5 = int(mw > 500) + int(logp > 5) + int(hbd > 5) + int(hba > 10)
    rows.append(dict(chembl_id=cid, smiles=Chem.MolToSmiles(m), mw=round(mw, 1), logp=round(logp, 2),
        hbd=hbd, hba=hba, tpsa=round(tpsa, 1), rot_bonds=rb, qed=round(qed, 3),
        arom_rings=ar, ro5_violations=ro5, max_phase=r.get("max_phase")))

lib = pd.DataFrame(rows)
# Drug-like filter: Ro5<=1 violation, MW 250-600, reasonable flexibility
libf = lib[(lib.ro5_violations <= 1) & (lib.mw.between(250, 600)) & (lib.rot_bonds <= 12) & (lib.tpsa <= 160)].copy()
libf = libf.sort_values("qed", ascending=False).reset_index(drop=True)
print(f"after drug-like filter: {len(libf)} / {len(lib)}")
libf.to_csv("screening_library.csv", index=False)