# Auto-extracted generating script
# Produces: fig_target_volcano_annotated.png
# Conda env: python   (run with this environment activated)
# Inputs (expected alongside / in data/): target_universe.csv
# Source artifact version: 46e8873d-976c-4b02-8cf6-70bf0d3f56ab
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
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


sig = pd.read_csv("target_universe.csv")
meta = sig.copy()

apply_figure_style()

d = meta.copy()
d["neglog10FDR"] = -np.log10(d.FDR.clip(lower=1e-300))
d["sig"] = (d.FDR < 0.05) & (d.logFC.abs() > 1)

fig, ax = plt.subplots(figsize=(7.4, 5.6))
ns = d[~d.sig]; up = d[d.sig & (d.logFC > 0)]; dn = d[d.sig & (d.logFC < 0)]
ax.scatter(ns.logFC, ns.neglog10FDR, s=4, c=META_GREY, alpha=0.25, lw=0, rasterized=True)
ax.scatter(up.logFC, up.neglog10FDR, s=6, c="#c0392b", alpha=0.45, lw=0, rasterized=True, label=f"Up (n={len(up)})")
ax.scatter(dn.logFC, dn.neglog10FDR, s=6, c="#2c7fb8", alpha=0.45, lw=0, rasterized=True, label=f"Down (n={len(dn)})")
hi = ["STAT3", "IL17A", "RORC", "JAK2", "JAK1", "JAK3", "NFKBIA", "S100A9", "IL36A", "VNN3", "CXCL13", "MMP9", "PTGS2", "NOS2", "IDO1"]
present = d[d.gene.isin(hi)]
ax.scatter(present.logFC, present.neglog10FDR, s=30, facecolors="none", edgecolors="black", lw=0.9, zorder=5)
texts = [ax.text(r.logFC, r.neglog10FDR, r.gene, fontsize=6.8, fontstyle='italic', zorder=6) for _, r in present.iterrows()]
adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='0.45', lw=0.5),
            expand_points=(1.6, 1.8), force_text=(0.5, 0.8), force_points=(0.3, 0.5))
ax.axvline(1, color='k', lw=0.5, ls=':'); ax.axvline(-1, color='k', lw=0.5, ls=':'); ax.axhline(-np.log10(0.05), color='k', lw=0.5, ls=':')
ax.set_xlabel("Meta-analysis log$_2$FC (PP vs NN)"); ax.set_ylabel("$-$log$_{10}$ FDR")
ax.set_title("Psoriasis meta-analysis DEGs: STAT3 and signalling nodes among 2,379 hits", fontsize=8)
ax.margins(0.05); ax.legend(frameon=False, loc='upper left', fontsize=6.5)
fig.tight_layout()
fig.savefig("fig_target_volcano_annotated.png", dpi=200)