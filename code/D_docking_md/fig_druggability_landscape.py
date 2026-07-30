# Auto-extracted generating script
# Produces: fig_druggability_landscape.png
# Conda env: python   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: e38ec650-8f2a-4c4f-9670-8d36442bf997
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import pandas as pd
import numpy as np
import json
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from adjustText import adjust_text

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


apply_figure_style()

ann = pd.read_csv("work/druggability_annotation.csv")

order = ["Unknown","Druggable family","Druggable pocket","Discovery precedence","Phase 1","Advanced clinical","Approved drug"]
ann["ybase"] = ann.ot_sm_bucket.map({b:i for i,b in enumerate(order)})
rng = np.random.default_rng(0)
ann["yj"] = ann.ybase + rng.uniform(-0.26,0.26,len(ann))
ann["absLFC"] = ann.logFC.abs()
ann["msize"] = 18 + 3.2*np.sqrt(ann.n_potent_inhib.clip(lower=0))
col = np.where(ann.direction=="up","#c0392b","#2c7fb8")

fig, ax = plt.subplots(figsize=(8.4,5.8))
ax.scatter(ann.absLFC, ann.yj, s=ann.msize, c=col, alpha=0.55, lw=0.4, edgecolors="white", zorder=3)
lab=["RORC","JAK3","STAT3","IL17A","CXCR2","NOS2","IDO1","MMP9","MMP12","PLK1","AURKA","HRH3","CCL2","BIRC3","IL36G","LTB4R","PLA2G2A","CXCR1"]
labdf=ann[ann.gene.isin(lab)].copy()
texts=[ax.text(r.absLFC, r.yj, r.gene, fontsize=6.6, fontstyle='italic', zorder=6) for _,r in labdf.iterrows()]
adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='0.5', lw=0.45),
            expand_points=(1.5,1.7), force_text=(0.5,0.8), force_points=(0.4,0.6))
ax.set_yticks(range(len(order))); ax.set_yticklabels(order)
ax.set_xlabel("Meta-analysis |log$_2$FC| (PP vs NN)")
ax.set_ylabel("Small-molecule tractability (Open Targets)")
ax.set_title("Druggability landscape of psoriasis DEG-encoded proteins", fontsize=8.5)
ax.margins(x=0.05, y=0.06)
leg1=[Line2D([0],[0],marker='o',color='w',markerfacecolor='#c0392b',markersize=6,label='Up in lesion'),
      Line2D([0],[0],marker='o',color='w',markerfacecolor='#2c7fb8',markersize=6,label='Down in lesion')]
l1=ax.legend(handles=leg1, loc='center right', frameon=False, fontsize=6.5, title="DE direction", title_fontsize=6.5)
ax.add_artist(l1)
for n in [10,1000,9000]:
    ax.scatter([],[],s=18+3.2*np.sqrt(n),c='0.6',alpha=0.6,label=f'{n:,}')
ax.legend(loc='lower right', frameon=False, fontsize=6, scatterpoints=1, labelspacing=1.3,
          borderpad=0.9, title="ChEMBL potent\ninhibitors (pChEMBL≥7)", title_fontsize=6)
fig.tight_layout()
fig.savefig("fig_druggability_landscape.png", dpi=200)