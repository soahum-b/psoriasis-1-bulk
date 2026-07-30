# Auto-extracted generating script
# Produces: structure_inventory.csv
# Conda env: dock-md   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: d2ef34c9-7cea-4a90-953d-78153a0fd7ba
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import os
import shutil
import json
import numpy as np
import pandas as pd

os.makedirs("work/structures", exist_ok=True)

selected = {
 "STAT3": {"pdb":"6NJS","ligand":"KQV","res":2.7,"domain":"SH2 / core","uniprot":"P40763"},
 "RORC":  {"pdb":"5APH","ligand":"VYI","res":1.54,"domain":"Ligand-binding domain (LBD)","uniprot":"P51449"},
 "JAK3":  {"pdb":"5LWM","ligand":"79T","res":1.55,"domain":"Kinase domain (ATP site)","uniprot":"P52333"},
}

# Download PDB files
import urllib.request
for sym, info in selected.items():
    pdb_id = info["pdb"]
    dest = f"/tmp/{pdb_id}.pdb"
    if not os.path.exists(dest):
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        urllib.request.urlretrieve(url, dest)
    shutil.copy(dest, f"work/structures/{sym}_{pdb_id}.pdb")

def parse_ligand_and_center(pdb_path, ligname):
    """Extract reference ligand coords + centroid + protein chains from a PDB."""
    coords=[]; lig_lines=[]; chains=set(); prot_atoms=0
    with open(pdb_path) as fh:
        for ln in fh:
            if ln.startswith(("ATOM","HETATM")):
                chains.add(ln[21])
            if ln.startswith("ATOM"): prot_atoms+=1
            if ln.startswith("HETATM") and ln[17:20].strip()==ligname:
                x,y,z=float(ln[30:38]),float(ln[38:46]),float(ln[46:54])
                coords.append((x,y,z)); lig_lines.append(ln)
    coords=np.array(coords)
    ctr = coords.mean(0)
    size = coords.max(0)-coords.min(0)
    return ctr, size, lig_lines, sorted(chains), prot_atoms

inv=[]
for sym,info in selected.items():
    p=f"work/structures/{sym}_{info['pdb']}.pdb"
    ctr,size,lig_lines,chains,natom = parse_ligand_and_center(p, info["ligand"])
    # write reference ligand out
    with open(f"work/structures/{sym}_ref_ligand_{info['ligand']}.pdb","w") as o:
        o.writelines(lig_lines); o.write("END\n")
    # docking box = ligand bbox + 8 A padding on each side
    box = size + 16.0
    info.update(dict(pocket_center=[round(float(c),2) for c in ctr],
                     ligand_bbox=[round(float(s),2) for s in size],
                     docking_box=[round(float(b),2) for b in box],
                     n_chains=len(chains), chains=chains, n_protein_atoms=natom))
    inv.append(dict(gene=sym, pdb=info["pdb"], ligand=info["ligand"], resolution_A=info["res"],
                    domain=info["domain"], uniprot=info["uniprot"],
                    pocket_center_x=info["pocket_center"][0], pocket_center_y=info["pocket_center"][1],
                    pocket_center_z=info["pocket_center"][2],
                    box_x=info["docking_box"][0], box_y=info["docking_box"][1], box_z=info["docking_box"][2],
                    n_chains=len(chains)))
    print(f"{sym} {info['pdb']} ({info['ligand']}): center={info['pocket_center']} box={info['docking_box']} chains={chains}")

inv_df=pd.DataFrame(inv)
inv_df.to_csv("work/structures/structure_inventory.csv", index=False)
sel_ser={k:{kk:(sorted(vv) if isinstance(vv,set) else vv) for kk,vv in v.items()} for k,v in selected.items()}
json.dump(sel_ser, open("work/structures/selected_structures.json","w"), indent=1)
print(inv_df.to_string(index=False))