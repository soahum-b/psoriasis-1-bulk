# Auto-extracted generating script
# Produces: fig_lead_structures.png
# Conda env: dock-md   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 9f03a5c8-5558-4997-b727-ee07bfd5bd52
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

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


def read_ca_and_lig(pdb_path, ligname, chain=None):
    ca = {}  # chain -> list of (resi, x,y,z)
    lig = []
    with open(pdb_path) as fh:
        for ln in fh:
            if ln.startswith("ATOM") and ln[12:16].strip() == "CA":
                ch = ln[21]
                if chain and ch != chain: continue
                ca.setdefault(ch, []).append((int(ln[22:26]), float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
            if ln.startswith("HETATM") and ln[17:20].strip() == ligname:
                lig.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    return ca, np.array(lig)


mpl.rcParams.update({"font.size": 8, "axes.linewidth": 0.6, "figure.dpi": 120,
                     "font.family": "DejaVu Sans", "axes.titlesize": 9})

info_render = {
    "STAT3": {"pdb": "6NJS", "lig": "KQV", "chain": "A", "domain": "SH2 / core domain", "col": "#c0392b"},
    "RORC":  {"pdb": "5APH", "lig": "VYI", "chain": "A", "domain": "Ligand-binding domain", "col": "#2c7fb8"},
    "JAK3":  {"pdb": "5LWM", "lig": "79T", "chain": "A", "domain": "Kinase domain (ATP site)", "col": "#27ae60"},
}
fig = plt.figure(figsize=(11, 4.0))
for i, (sym, ri) in enumerate(info_render.items()):
    p = f"work/structures/{sym}_{ri['pdb']}.pdb"
    ca, lig = read_ca_and_lig(p, ri["lig"], chain=ri["chain"])
    ax = fig.add_subplot(1, 3, i + 1, projection='3d')
    for ch, atoms in ca.items():
        arr = np.array([(a[1], a[2], a[3]) for a in atoms])
        n = len(arr)
        for j in range(n - 1):
            ax.plot(arr[j:j+2, 0], arr[j:j+2, 1], arr[j:j+2, 2], color=plt.cm.Spectral(j / n), lw=1.4, solid_capstyle='round')
    if len(lig):
        ax.scatter(lig[:, 0], lig[:, 1], lig[:, 2], s=32, c=ri["col"], edgecolors='k', linewidths=0.4, depthshade=True, zorder=10)
        c = lig.mean(0)
    ax.set_title(f"{sym}  ({ri['pdb']})\n{ri['domain']}", fontsize=8.5)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    ax.grid(False); ax.set_axis_off()
    all_pts = np.vstack([np.array([(a[1], a[2], a[3]) for a in atoms]) for atoms in ca.values()])
    mid = all_pts.mean(0); rng = (all_pts.max(0) - all_pts.min(0)).max() / 2
    ax.set_xlim(mid[0] - rng, mid[0] + rng); ax.set_ylim(mid[1] - rng, mid[1] + rng); ax.set_zlim(mid[2] - rng, mid[2] + rng)
    ax.view_init(elev=18, azim=45 + i * 30)
fig.suptitle("Lead target structures — Cα ribbon (N→C rainbow) with pocket ligand (dark spheres)", fontsize=9, y=1.02)
fig.tight_layout()
fig.savefig("fig_lead_structures.png", dpi=190, bbox_inches='tight')
print("saved fig_lead_structures.png")