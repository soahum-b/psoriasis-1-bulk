# Auto-extracted generating script
# Produces: fig_library_properties.png
# Conda env: dock-md   (run with this environment activated)
# Inputs (expected alongside / in data/): screening_library.csv
# Source artifact version: 26c37872-f9f2-4f59-a4bc-93a5a4bf1e52
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from rdkit import Chem
from rdkit.Chem import Descriptors, QED, Crippen

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


mpl.rcParams.update({"font.size":8,"axes.linewidth":0.6,"font.family":"DejaVu Sans"})

libf = pd.read_csv("screening_library.csv")

fig = plt.figure(figsize=(11, 6.2))
gs = fig.add_gridspec(2, 4, hspace=0.42, wspace=0.34, height_ratios=[1, 1])
props = [("mw", "Molecular weight (Da)", (250, 600)), ("logp", "cLogP", None),
         ("qed", "QED drug-likeness", (0, 1)), ("tpsa", "TPSA (Å²)", None),
         ("hbd", "H-bond donors", None), ("hba", "H-bond acceptors", None),
         ("rot_bonds", "Rotatable bonds", None), ("arom_rings", "Aromatic rings", None)]
for i, (col, lab, xr) in enumerate(props):
    ax = fig.add_subplot(gs[i // 4, i % 4])
    ax.hist(libf[col], bins=18, color="#4575b4", edgecolor="white", lw=0.4)
    ax.set_xlabel(lab, fontsize=7)
    ax.set_ylabel("count", fontsize=7)
    if xr:
        ax.set_xlim(*xr)
    ax.axvline(libf[col].median(), color="#d73027", lw=1, ls="--")
    ax.tick_params(labelsize=6.5)
fig.suptitle(f"RORγt-focused purchasable screening library — {len(libf)} drug-like compounds (ChEMBL bioactivity space)", fontsize=9, y=0.98)
fig.savefig("fig_library_properties.png", dpi=190, bbox_inches="tight")
print("saved property fig; medians: MW", libf.mw.median(), "QED", libf.qed.median(), "logP", libf.logp.median())