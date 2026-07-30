# Auto-extracted generating script
# Produces: fig_pipeline_funnel.png
# Conda env: dock-md   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 28dbface-d129-452a-9e56-02f70a3a735f
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import FancyBboxPatch

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

fig,ax=plt.subplots(figsize=(8.5,5.6))
ax.axis("off")
stages=[
 ("Meta-analysis DEGs\n(PP vs NN, FDR<0.05 & |log2FC|>1)", "n = 2,379", "#34495e", 0.0),
 ("Protein-coding candidates\nannotated for druggability", "n = 154", "#2c7fb8", 0.12),
 ("Tractable & on psoriasis axis\n(≥ Discovery precedence)", "n = 28", "#16a085", 0.24),
 ("Structure-ready leads\n(pocket + ligand-bound PDB)", "STAT3 · RORC · JAK3  (n = 3)", "#e67e22", 0.36),
 ("RORC: docked (151-cmpd screen) + MD systems built\nSTAT3 · JAK3: docking + co-folding cluster packages", "", "#c0392b", 0.30),
]
y=len(stages); maxw=0.86
for i,(lab,n,col,taper) in enumerate(stages):
    w=maxw*(1-taper); x=(1-w)/2; yy=y-i-1
    box=FancyBboxPatch((x,yy),w,0.72, boxstyle="round,pad=0.01,rounding_size=0.02",
                       fc=col, ec="white", lw=1.5, alpha=0.92)
    ax.add_patch(box)
    ty = yy+0.36 if not n else yy+0.46
    ax.text(0.5, ty, lab, ha="center", va="center", color="white", fontsize=6.9, fontweight="bold")
    if n: ax.text(0.5, yy+0.15, n, ha="center", va="center", color="white", fontsize=6.8)
    if i<len(stages)-1:
        ax.annotate("", xy=(0.5,yy-0.02), xytext=(0.5,yy+0.02), arrowprops=dict(arrowstyle="-|>", color="0.4", lw=1.3))
ax.set_xlim(0,1); ax.set_ylim(-0.3,y)
ax.set_title("Psoriasis meta-analysis → druggable target triage → docking/MD pipeline", fontsize=9.5, pad=12)
ax.text(0.02, y-1.5, "Open Targets\n+ ChEMBL", fontsize=6, color="0.4", ha="left", style="italic")
ax.text(0.02, y-3.5, "fpocket +\nPDB structures", fontsize=6, color="0.4", ha="left", style="italic")
ax.text(0.985, y-4.55, "Vina · Boltz-2 · OpenMM", fontsize=6, color="0.4", ha="right", style="italic")
fig.savefig("fig_pipeline_funnel.png", dpi=190, bbox_inches="tight")