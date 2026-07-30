# Auto-extracted generating script
# Produces: all_docked_poses.tar.gz
# Conda env: dock-md   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 9c8a0fce-0b95-48ac-ab99-8b869d64ae09
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import os
import glob
import json
import time
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, QED, Crippen, AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
from vina import Vina

BASE = "/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3"
WP = "/tmp/psor_pathway"
os.makedirs(WP, exist_ok=True)
os.makedirs("work", exist_ok=True)
os.makedirs("work/structures", exist_ok=True)
os.makedirs("work/ligands_3d", exist_ok=True)
os.makedirs("work/docking/poses", exist_ok=True)
os.makedirs("handoff", exist_ok=True)

# --- Load and prepare screening library ---
import json as _json

rorc_actives_smiles = _json.load(open("handoff/rorc_actives_smiles.json"))
chembl_simlib = _json.load(open("handoff/chembl_simlib.json"))

allc = {r["chembl_id"]: r for r in rorc_actives_smiles}
for s in chembl_simlib:
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
libf["source"] = "ChEMBL RORC bioactivity (pChEMBL>=7) + similarity analogs"
libf["target"] = "RORC/RORgt (CHEMBL1741186)"
libf.to_csv("work/screening_library.csv", index=False)

# --- Generate 3D conformers ---
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

# --- Prepare ligand PDBQTs ---
os.makedirs("work/ligands_3d/pdbqt", exist_ok=True)
os.system("mk_prepare_ligand.py -i work/ligands_3d/screening_library_3d.sdf --multimol_outdir work/ligands_3d/pdbqt 2>/dev/null")

# --- Run docking ---
rec = "work/structures/RORC_receptor.pdbqt"
lig_dir = "work/ligands_3d/pdbqt"
center = [-25.8, 5.4, -15.3]
box = [24, 24, 24]
ligs = sorted(glob.glob(f"{lig_dir}/*.pdbqt"))
results = []
t0 = time.time()
for i, lp in enumerate(ligs):
    name = os.path.basename(lp).replace(".pdbqt", "")
    try:
        v = Vina(sf_name='vina', cpu=4, verbosity=0)
        v.set_receptor(rec)
        v.set_ligand_from_file(lp)
        v.compute_vina_maps(center=center, box_size=box)
        v.dock(exhaustiveness=8, n_poses=5)
        e = v.energies(n_poses=1)
        top = float(e[0][0])
        out = f"work/docking/poses/{name}_out.pdbqt"
        v.write_poses(out, n_poses=1, overwrite=True)
        results.append({"ligand": name, "affinity": top})
    except Exception as ex:
        results.append({"ligand": name, "affinity": None, "error": str(ex)[:80]})
    if (i + 1) % 25 == 0:
        print(f"{i+1}/{len(ligs)} docked, {time.time()-t0:.0f}s elapsed", flush=True)

json.dump(results, open("work/docking/dock_results.json", "w"))

# --- Create archive ---
os.system("cd work/docking && tar -czf all_docked_poses.tar.gz poses/ 2>/dev/null")