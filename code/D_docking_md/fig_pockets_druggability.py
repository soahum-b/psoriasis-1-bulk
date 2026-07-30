# Auto-extracted generating script
# Produces: fig_pockets_druggability.png
# Conda env: dock-md   (run with this environment activated)
# Inputs (expected alongside / in data/): STAT3_6NJS.pdb, RORC_5APH.pdb, JAK3_5LWM.pdb, STAT3_apo.pdb, pocket39_vert.pqr, RORC_apo.pdb, pocket1_vert.pqr, JAK3_apo.pdb, pocket25_vert.pqr
# Source artifact version: 51f8e17c-b51c-4194-9cb6-a405b5c24f5a
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams.update({"font.size": 8, "axes.linewidth": 0.6})


def read_ca(path, chain=None):
    ca = {}
    for ln in open(path):
        if ln.startswith("ATOM") and ln[12:16].strip() == "CA":
            ch = ln[21]
            if chain and ch != chain:
                continue
            ca.setdefault(ch, []).append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    return {k: np.array(v) for k, v in ca.items()}


def read_lig_atoms(pdb_path, ligname):
    a = []
    for ln in open(pdb_path):
        if ln.startswith("HETATM") and ln[17:20].strip() == ligname:
            a.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    return np.array(a)


def pocket_alpha_coords(vpath):
    c = []
    for ln in open(vpath):
        if ln.startswith(("ATOM", "HETATM")):
            c.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    return np.array(c) if c else None


ligmap = {
    "STAT3": ("6NJS", "KQV", 39),
    "RORC": ("5APH", "VYI", 1),
    "JAK3": ("5LWM", "79T", 25),
}

sel = {
    "STAT3": {"pdb": "6NJS", "ligand": "KQV", "res": 2.7, "domain": "SH2 / core domain",
              "fpocket_druggability": None},
    "RORC": {"pdb": "5APH", "ligand": "VYI", "res": 1.54, "domain": "Ligand-binding domain (LBD)",
             "fpocket_druggability": None},
    "JAK3": {"pdb": "5LWM", "ligand": "79T", "res": 1.55, "domain": "Kinase domain (ATP site)",
             "fpocket_druggability": None},
}

apo_paths = {
    "STAT3": "STAT3_apo.pdb",
    "RORC": "RORC_apo.pdb",
    "JAK3": "JAK3_apo.pdb",
}

full_pdb_paths = {
    "STAT3": "STAT3_6NJS.pdb",
    "RORC": "RORC_5APH.pdb",
    "JAK3": "JAK3_5LWM.pdb",
}

pocket_pqr_paths = {
    "STAT3": "pocket39_vert.pqr",
    "RORC": "pocket1_vert.pqr",
    "JAK3": "pocket25_vert.pqr",
}

# pocket druggability summary (reconstructed from trace values)
summ_data = {
    "gene": ["STAT3", "RORC", "JAK3"],
    "ligand_pocket_druggability": [None, None, None],
    "chembl_potent_inhibitors": [112, 9766, 6069],
}

# Parse fpocket druggability from pocket report data by recomputing from pqr/apo files
import re


def parse_fpocket_info_from_lines(lines):
    pockets = []
    cur = None
    for ln in lines:
        m = re.match(r"Pocket (\d+) :", ln)
        if m:
            if cur:
                pockets.append(cur)
            cur = {"pocket": int(m.group(1))}
        elif cur is not None and ":" in ln:
            k, v = ln.split(":", 1)
            k = k.strip()
            v = v.strip().rstrip()
            try:
                cur[k] = float(v.split()[0])
            except:
                pass
    if cur:
        pockets.append(cur)
    return pockets


# We'll use the druggability values known from the trace output
# STAT3 pocket39 druggability, RORC pocket1 druggability, JAK3 pocket25 druggability
# From the trace: pocket report printed for each gene
# We set them from the pocket_druggability_summary known values
ligand_pocket_druggability = {
    "STAT3": 0.608,  # from trace pocket_report output for STAT3 is_ligand_pocket row
    "RORC": 0.844,
    "JAK3": 0.783,
}

summ = pd.DataFrame({
    "gene": ["STAT3", "RORC", "JAK3"],
    "ligand_pocket_druggability": [ligand_pocket_druggability["STAT3"],
                                    ligand_pocket_druggability["RORC"],
                                    ligand_pocket_druggability["JAK3"]],
    "chembl_potent_inhibitors": [112, 9766, 6069],
})

fig = plt.figure(figsize=(11, 6.6))
gs = fig.add_gridspec(2, 3, height_ratios=[1.5, 1.0], hspace=0.38, wspace=0.12)

for i, (sym, (pdb, lig, pid)) in enumerate(ligmap.items()):
    ca = read_ca(apo_paths[sym])
    ligA = read_lig_atoms(full_pdb_paths[sym], lig)
    pocket = pocket_alpha_coords(pocket_pqr_paths[sym])
    ax = fig.add_subplot(gs[0, i], projection='3d')
    for ch, arr in ca.items():
        for j in range(len(arr) - 1):
            ax.plot(arr[j:j+2, 0], arr[j:j+2, 1], arr[j:j+2, 2], color='0.72', lw=0.7)
    if pocket is not None:
        ax.scatter(pocket[:, 0], pocket[:, 1], pocket[:, 2], s=42, c='#f39c12', alpha=0.16, edgecolors='none')
    ax.scatter(ligA[:, 0], ligA[:, 1], ligA[:, 2], s=22, c='#c0392b', edgecolors='k', linewidths=0.3, zorder=10)
    ax.set_title(f"{sym} ({pdb})\n{sel[sym]['domain']}", fontsize=7.6)
    ax.set_axis_off()
    allp = np.vstack(list(ca.values()))
    mid = allp.mean(0)
    rng = (allp.max(0) - allp.min(0)).max() / 2 * 0.92
    ax.set_xlim(mid[0] - rng, mid[0] + rng)
    ax.set_ylim(mid[1] - rng, mid[1] + rng)
    ax.set_zlim(mid[2] - rng, mid[2] + rng)
    ax.view_init(elev=20, azim=50 + i * 35)

ax2 = fig.add_subplot(gs[1, :])
x = np.arange(3)
w = 0.36
drg = summ.set_index("gene").loc[["STAT3", "RORC", "JAK3"]]
b1 = ax2.bar(x - w / 2, drg.ligand_pocket_druggability, w, color='#4575b4')
ax2b = ax2.twinx()
b2 = ax2b.bar(x + w / 2, np.log10(drg.chembl_potent_inhibitors + 1), w, color='#d73027', alpha=0.85)
ax2.set_xticks(x)
ax2.set_xticklabels(["STAT3 (SH2/core)", "RORC (LBD)", "JAK3 (ATP site)"], fontsize=7.5)
ax2.set_ylabel("fpocket druggability", color='#4575b4', fontsize=7.5)
ax2b.set_ylabel("log$_{10}$(inhibitors+1)", color='#d73027', fontsize=7.5)
ax2.set_ylim(0, 1.05)
ax2.tick_params(axis='y', colors='#4575b4')
ax2b.tick_params(axis='y', colors='#d73027')
for xi, (g, r) in zip(x, drg.iterrows()):
    ax2.text(xi - w / 2, r.ligand_pocket_druggability + 0.02, f"{r.ligand_pocket_druggability:.2f}",
             ha='center', fontsize=6.5, color='#4575b4')
    ax2b.text(xi + w / 2, np.log10(r.chembl_potent_inhibitors + 1) + 0.05,
              f"{int(r.chembl_potent_inhibitors):,}", ha='center', fontsize=6.5, color='#d73027')
ax2.set_title("Structural vs chemical druggability of the lead pockets", fontsize=8.2)
ax2.legend(handles=[b1, b2],
           labels=['fpocket druggability (structural)', 'ChEMBL potent inhibitors (chemical)'],
           loc='upper center', bbox_to_anchor=(0.5, -0.16), frameon=False, fontsize=7, ncol=2)
fig.suptitle("Binding-pocket detection & structural druggability of psoriasis lead targets", fontsize=9.2, y=0.99)
fig.savefig("fig_pockets_druggability.png", dpi=190, bbox_inches='tight')
print("saved final pocket figure")