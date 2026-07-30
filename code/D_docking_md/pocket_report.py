# Auto-extracted generating script
# Produces: pocket_report.csv
# Conda env: dock-md   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: ae05ff5b-eced-4441-b405-5b6c483f6b98
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import pandas as pd
import numpy as np
import re
import glob
import json
import os

# skill:figure-style kernel.py (auto-injected on skill load)
META_GREY = "#888888"


def apply_figure_style(*, frame="open", font=None, sizes=(8, 7, 6), grid=False):
    import matplotlib as mpl
    if frame not in ("open", "boxed", "none"):
        raise ValueError(f"frame must be 'open'|'boxed'|'none', got {frame!r}")
    try:
        import os, sys, glob, matplotlib.font_manager as fm
        fdir = os.path.join(os.environ.get("CONDA_PREFIX") or sys.prefix, "fonts")
        if os.path.isdir(fdir):
            known = {f.fname for f in fm.fontManager.ttflist}
            for f in glob.glob(os.path.join(fdir, "*.ttf")):
                if f not in known:
                    fm.fontManager.addfont(f)
    except Exception:
        pass
    base, secondary, tick = sizes
    boxed = (frame == "boxed")
    rc = {
        "font.family": "sans-serif",
        "font.size": base,
        "axes.labelsize": base,
        "axes.titlesize": base,
        "legend.fontsize": secondary,
        "xtick.labelsize": tick,
        "ytick.labelsize": tick,
        "axes.linewidth": 0.6,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "axes.spines.top": boxed, "axes.spines.right": boxed,
        "axes.spines.left": frame != "none", "axes.spines.bottom": frame != "none",
        "axes.grid": bool(grid),
        "legend.frameon": False,
        "figure.dpi": 200,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.titleweight": "normal",
        "axes.titlelocation": "left",
        "axes.labelweight": "normal",
        "lines.linewidth": 1.2,
        "patch.linewidth": 0.6,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    }
    if font:
        rc["font.sans-serif"] = [font, "DejaVu Sans"]
    mpl.rcParams.update(rc)


BASE = "/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3"
WP = "/tmp/psor_pathway"
os.makedirs("work", exist_ok=True)
os.makedirs("work/structures", exist_ok=True)

# ---- Reproduce structure inventory (required for selected_structures.json) ----
import shutil

selected_raw = {
    "STAT3": {"pdb": "6NJS", "ligand": "KQV", "res": 2.7, "domain": "SH2 / core", "uniprot": "P40763"},
    "RORC":  {"pdb": "5APH", "ligand": "VYI", "res": 1.54, "domain": "Ligand-binding domain (LBD)", "uniprot": "P51449"},
    "JAK3":  {"pdb": "5LWM", "ligand": "79T", "res": 1.55, "domain": "Kinase domain (ATP site)", "uniprot": "P52333"},
}

for sym, info in selected_raw.items():
    src = f"/tmp/{info['pdb']}.pdb"
    dst = f"work/structures/{sym}_{info['pdb']}.pdb"
    if not os.path.exists(dst):
        shutil.copy(src, dst)

def parse_ligand_and_center(pdb_path, ligname):
    coords = []; lig_lines = []; chains = set(); prot_atoms = 0
    with open(pdb_path) as fh:
        for ln in fh:
            if ln.startswith(("ATOM", "HETATM")):
                chains.add(ln[21])
            if ln.startswith("ATOM"):
                prot_atoms += 1
            if ln.startswith("HETATM") and ln[17:20].strip() == ligname:
                x, y, z = float(ln[30:38]), float(ln[38:46]), float(ln[46:54])
                coords.append((x, y, z)); lig_lines.append(ln)
    coords = np.array(coords)
    ctr = coords.mean(0)
    size = coords.max(0) - coords.min(0)
    return ctr, size, lig_lines, sorted(chains), prot_atoms

selected = {}
inv = []
for sym, info in selected_raw.items():
    p = f"work/structures/{sym}_{info['pdb']}.pdb"
    ctr, size, lig_lines, chains, natom = parse_ligand_and_center(p, info["ligand"])
    with open(f"work/structures/{sym}_ref_ligand_{info['ligand']}.pdb", "w") as o:
        o.writelines(lig_lines); o.write("END\n")
    box = size + 16.0
    info_upd = dict(info)
    info_upd.update(dict(
        pocket_center=[round(float(c), 2) for c in ctr],
        ligand_bbox=[round(float(s), 2) for s in size],
        docking_box=[round(float(b), 2) for b in box],
        n_chains=len(chains), chains=sorted(chains), n_protein_atoms=natom
    ))
    selected[sym] = info_upd
    inv.append(dict(
        gene=sym, pdb=info["pdb"], ligand=info["ligand"], resolution_A=info["res"],
        domain=info["domain"], uniprot=info["uniprot"],
        pocket_center_x=info_upd["pocket_center"][0],
        pocket_center_y=info_upd["pocket_center"][1],
        pocket_center_z=info_upd["pocket_center"][2],
        box_x=info_upd["docking_box"][0],
        box_y=info_upd["docking_box"][1],
        box_z=info_upd["docking_box"][2],
        n_chains=len(chains)
    ))

inv_df = pd.DataFrame(inv)
inv_df.to_csv("work/structures/structure_inventory.csv", index=False)
sel_ser = {k: {kk: (sorted(vv) if isinstance(vv, set) else vv) for kk, vv in v.items()} for k, v in selected.items()}
json.dump(sel_ser, open("work/structures/selected_structures.json", "w"), indent=1)

# ---- Generate apo PDB files ----
for sym in ["STAT3", "RORC", "JAK3"]:
    info = selected_raw[sym]
    pdb = info["pdb"]
    apo_path = f"work/structures/{sym}_apo.pdb"
    with open(f"work/structures/{sym}_{pdb}.pdb") as fh, open(apo_path, "w") as out:
        for ln in fh:
            if ln.startswith(("ATOM", "TER", "HEADER")):
                out.write(ln)
        out.write("END\n")

# ---- Run fpocket ----
import subprocess
for sym in ["STAT3", "RORC", "JAK3"]:
    apo_path = f"work/structures/{sym}_apo.pdb"
    subprocess.run(["fpocket", "-f", apo_path], capture_output=True, cwd="work/structures")

# ---- Parse fpocket results ----
def parse_fpocket_info(info_path):
    pockets = []
    cur = None
    for ln in open(info_path):
        m = re.match(r"Pocket (\d+) :", ln)
        if m:
            if cur: pockets.append(cur)
            cur = {"pocket": int(m.group(1))}
        elif cur is not None and ":" in ln:
            k, v = ln.split(":", 1)
            k = k.strip(); v = v.strip().rstrip()
            try: cur[k] = float(v.split()[0])
            except: pass
    if cur: pockets.append(cur)
    return pockets

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

sel = json.load(open("work/structures/selected_structures.json"))
ligmap = {"STAT3": ("6NJS", "KQV"), "RORC": ("5APH", "VYI"), "JAK3": ("5LWM", "79T")}

rows = []; matched_pockets = {}
for sym, (pdb, lig) in ligmap.items():
    ligA = read_lig_atoms(f"work/structures/{sym}_{pdb}.pdb", lig)
    base = f"work/structures/{sym}_apo_out"
    info = parse_fpocket_info(f"{base}/{sym}_apo_info.txt")
    best = None; best_overlap = -1
    for pk in info:
        pid = pk["pocket"]
        v = glob.glob(f"{base}/pockets/pocket{pid}_vert.pqr")
        ac = pocket_alpha_coords(v[0]) if v else None
        if ac is None: pk["overlap"] = 0; pk["center"] = None; continue
        d = np.sqrt(((ligA[:, None, :] - ac[None, :, :]) ** 2).sum(-1))
        overlap = int((d.min(1) < 4.5).sum())
        pk["overlap"] = overlap; pk["center"] = ac.mean(0)
        if overlap > best_overlap: best_overlap = overlap; best = pk
    matched_pockets[sym] = best
    for pk in info:
        rows.append(dict(
            gene=sym, pocket=pk["pocket"],
            druggability=pk.get("Druggability Score"), fpocket_score=pk.get("Score"),
            volume=pk.get("Volume"), n_alpha=pk.get("Number of Alpha Spheres"),
            hydrophob=pk.get("Hydrophobicity score"), polarity=pk.get("Polarity score"),
            lig_atoms_in_pocket=pk.get("overlap", 0),
            is_ligand_pocket=(best is not None and pk["pocket"] == best["pocket"])))
    c = best["center"]
    print(f"{sym}: matched pocket #{best['pocket']} — {best_overlap}/{len(ligA)} ligand atoms overlap | "
          f"druggability={best.get('Druggability Score')} vol={best.get('Volume'):.0f} "
          f"hydrophob={best.get('Hydrophobicity score')} | center=({c[0]:.1f},{c[1]:.1f},{c[2]:.1f})")

pdf = pd.DataFrame(rows)
pdf.to_csv("work/pocket_report.csv", index=False)

for sym in sel:
    c = matched_pockets[sym]["center"]
    sel[sym]["fpocket_center"] = [round(float(x), 2) for x in c]
    sel[sym]["fpocket_druggability"] = matched_pockets[sym].get("Druggability Score")
    sel[sym]["fpocket_volume"] = matched_pockets[sym].get("Volume")
    sel[sym]["matched_pocket_id"] = matched_pockets[sym]["pocket"]
json.dump(sel, open("work/structures/selected_structures.json", "w"), indent=1)

print("\n=== Ligand-matched (docking) pockets ===")
print(pdf[pdf.is_ligand_pocket][["gene", "pocket", "druggability", "fpocket_score", "volume", "hydrophob", "polarity", "lig_atoms_in_pocket"]].to_string(index=False))