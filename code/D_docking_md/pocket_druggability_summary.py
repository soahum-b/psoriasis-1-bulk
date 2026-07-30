# Auto-extracted generating script
# Produces: pocket_druggability_summary.csv
# Conda env: dock-md   (run with this environment activated)
# Inputs (expected alongside / in data/): druggable_target_ranking.csv, selected_structures.json
# Source artifact version: 2bff9ff2-6890-4f8d-952b-608cd39ee3e2
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import pandas as pd
import numpy as np
import json

sel = json.load(open("selected_structures.json"))

# pocket_report.csv is needed — rebuild from fpocket parsing logic
import re, glob

def parse_fpocket_info(info_path):
    pockets=[]
    cur=None
    for ln in open(info_path):
        m=re.match(r"Pocket (\d+) :", ln)
        if m:
            if cur: pockets.append(cur)
            cur={"pocket":int(m.group(1))}
        elif cur is not None and ":" in ln:
            k,v=ln.split(":",1)
            k=k.strip(); v=v.strip().rstrip()
            try: cur[k]=float(v.split()[0])
            except: pass
    if cur: pockets.append(cur)
    return pockets

def pocket_alpha_coords(vpath):
    c=[]
    for ln in open(vpath):
        if ln.startswith(("ATOM","HETATM")):
            c.append((float(ln[30:38]),float(ln[38:46]),float(ln[46:54])))
    return np.array(c) if c else None

def read_lig_atoms(pdb_path, ligname):
    a=[]
    for ln in open(pdb_path):
        if ln.startswith("HETATM") and ln[17:20].strip()==ligname:
            a.append((float(ln[30:38]),float(ln[38:46]),float(ln[46:54])))
    return np.array(a)

ligmap={"STAT3":("6NJS","KQV"),"RORC":("5APH","VYI"),"JAK3":("5LWM","79T")}

rows=[]; matched_pockets={}
for sym,(pdb,lig) in ligmap.items():
    ligA = read_lig_atoms(f"work/structures/{sym}_{pdb}.pdb", lig)
    base=f"work/structures/{sym}_apo_out"
    info=parse_fpocket_info(f"{base}/{sym}_apo_info.txt")
    best=None; best_overlap=-1
    for pk in info:
        pid=pk["pocket"]
        v=glob.glob(f"{base}/pockets/pocket{pid}_vert.pqr")
        ac = pocket_alpha_coords(v[0]) if v else None
        if ac is None: pk["overlap"]=0; pk["center"]=None; continue
        d = np.sqrt(((ligA[:,None,:]-ac[None,:,:])**2).sum(-1))
        overlap = int((d.min(1)<4.5).sum())
        pk["overlap"]=overlap; pk["center"]=ac.mean(0)
        if overlap>best_overlap: best_overlap=overlap; best=pk
    matched_pockets[sym]=best
    for pk in info:
        rows.append(dict(gene=sym, pocket=pk["pocket"],
            druggability=pk.get("Druggability Score"), fpocket_score=pk.get("Score"),
            volume=pk.get("Volume"), n_alpha=pk.get("Number of Alpha Spheres"),
            hydrophob=pk.get("Hydrophobicity score"), polarity=pk.get("Polarity score"),
            lig_atoms_in_pocket=pk.get("overlap",0),
            is_ligand_pocket=(best is not None and pk["pocket"]==best["pocket"])))

pdf=pd.DataFrame(rows)

summ=[]
for sym in ["STAT3","RORC","JAK3"]:
    sub=pdf[pdf.gene==sym]
    ligpk=sub[sub.is_ligand_pocket].iloc[0]
    bestpk=sub.loc[sub.druggability.idxmax()]
    ref_ctr = sel[sym]["pocket_center"]
    summ.append(dict(
        gene=sym, pdb=sel[sym]["pdb"], domain=sel[sym]["domain"], resolution_A=sel[sym]["res"],
        ref_ligand=sel[sym]["ligand"],
        docking_center_x=ref_ctr[0], docking_center_y=ref_ctr[1], docking_center_z=ref_ctr[2],
        ligand_pocket_id=int(ligpk.pocket), ligand_pocket_druggability=round(ligpk.druggability,3),
        ligand_pocket_volume=round(ligpk.volume,0), ligand_pocket_hydrophob=round(ligpk.hydrophob,1),
        best_pocket_druggability=round(bestpk.druggability,3),
        chembl_potent_inhibitors=None,
    ))
summ=pd.DataFrame(summ)
ann=pd.read_csv("druggable_target_ranking.csv")
summ["chembl_potent_inhibitors"]=summ.gene.map(ann.set_index("gene").n_potent_inhib)
summ.to_csv("pocket_druggability_summary.csv", index=False)
print(summ.to_string(index=False))