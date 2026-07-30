# Auto-extracted generating script
# Produces: screening_library_3d.sdf
# Conda env: dock-md   (run with this environment activated)
# Inputs (expected alongside / in data/): rorc_actives_smiles.json, chembl_simlib.json
# Source artifact version: 6e0d86a1-457f-44fb-8aaf-9e87976ef2d3
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import json, pandas as pd, numpy as np, os
from rdkit import Chem
from rdkit.Chem import Descriptors, QED, Crippen, AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize

recs = json.load(open("rorc_actives_smiles.json"))
simlib = json.load(open("chembl_simlib.json"))

allc = {r["chembl_id"]: r for r in recs}
for s in simlib:
    allc.setdefault(s["chembl_id"], s)

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
libf = lib[(lib.ro5_violations <= 1) & (lib.mw.between(250, 600)) & (lib.rot_bonds <= 12) & (lib.tpsa <= 160)].copy()
libf = libf.sort_values("qed", ascending=False).reset_index(drop=True)

os.makedirs("work/ligands_3d", exist_ok=True)

writer = Chem.SDWriter("work/ligands_3d/screening_library_3d.sdf")
ok = 0
fail = []
for _, r in libf.iterrows():
    m = Chem.MolFromSmiles(r["smiles"])
    if m is None:
        fail.append(r["chembl_id"])
        continue
    m = Chem.AddHs(m)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if AllChem.EmbedMolecule(m, params) != 0:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(m, params) != 0:
            fail.append(r["chembl_id"])
            continue
    try:
        AllChem.MMFFOptimizeMolecule(m, maxIters=500)
    except:
        pass
    m.SetProp("_Name", r["chembl_id"])
    m.SetProp("QED", str(r["qed"]))
    m.SetProp("MW", str(r["mw"]))
    writer.write(m)
    ok += 1
writer.close()
print(f"3D conformers generated: {ok}/{len(libf)}, failed: {len(fail)}")
print("SDF saved: work/ligands_3d/screening_library_3d.sdf")