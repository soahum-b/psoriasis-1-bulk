# Auto-extracted generating script
# Produces: fig_il1b_stat3_NNvsPP.png
# Conda env: python   (run with this environment activated)
# Inputs (expected alongside / in data/): meta_de_PPvsNN.csv, meta_de_PNvsNN.csv, meta_de_PNvsPP.csv
# Source artifact version: 7b214293-212f-46c9-9325-3607cd09585d
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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


apply_figure_style()

contrasts = {
  "PPvsNN": "meta_de_PPvsNN.csv",
  "PNvsNN": "meta_de_PNvsNN.csv",
  "PNvsPP": "meta_de_PNvsPP.csv",
}
dfs = {c: pd.read_csv(v).set_index("gene") for c, v in contrasts.items()}

d = dfs["PPvsNN"]
mods = {
 "STAT3 / JAK-STAT": ["STAT3","STAT1","JAK3","TYK2","JAK1","SOCS3"],
 "IL-1$\\beta$ / inflammasome": ["IL1B","IL1RN","NLRP3","CASP1","PYCARD","AIM2","CASP5","GSDMD","IL1A","IL18"],
 "IL-36 family": ["IL36G","IL36A","IL36RN"],
 "Th17 axis": ["IL17A","IL23A","RORC"],
 "IL-1 output": ["S100A9","S100A8","S100A7","DEFB4A","CCL20","CXCL2","LCN2"],
}
rows = []
for m, gs in mods.items():
    for g in gs:
        if g in d.index:
            r = d.loc[g]
            rows.append((m, g, float(r["logFC"]), float(r["SE"]), float(r["FDR"]), int(r["k"]), float(r["I2"])))

F = pd.DataFrame(rows, columns=["module","gene","logFC","SE","FDR","k","I2"])
F["lo"] = F.logFC - 1.96 * F.SE
F["hi"] = F.logFC + 1.96 * F.SE

mod_list = list(mods.keys())
cmap = dict(zip(mod_list, [plt.cm.Dark2(i) for i in range(len(mod_list))]))

fig, ax = plt.subplots(figsize=(6.8, 9.5))
y = 0; sep = []
for m in mod_list:
    sub = F[F["module"]==m].sort_values("logFC")
    for _, r in sub.iterrows():
        c = cmap[m]; sig = r.FDR < 0.05
        ax.plot([r.lo, r.hi], [y, y], color=c, lw=1.6, alpha=0.9, zorder=2)
        ax.scatter([r.logFC], [y], s=46 if sig else 30, color=c,
                   edgecolor="black" if sig else "0.5", linewidth=0.9 if sig else 0.5,
                   zorder=3, alpha=1 if sig else 0.55)
        star = "***" if r.FDR < 1e-3 else "**" if r.FDR < 1e-2 else "*" if r.FDR < 0.05 else "n.s."
        het = " ‡" if r.I2 >= 75 else ""
        kflag = " ᵏ⁼²" if int(r.k) == 2 else ""
        ax.text(r.hi+0.2, y, f"{r.gene}  {star}{het}{kflag}", va="center", ha="left", fontsize=6.6,
                fontstyle="italic", color="black" if sig else "0.45")
        y += 1
    sep.append(y - 0.5); y += 0.8
ax.axvline(0, color="0.3", lw=1, ls="--", zorder=1)
for s in sep[:-1]:
    ax.axhline(s, color="0.9", lw=0.6)
y = 0
for m in mod_list:
    n = (F["module"]==m).sum()
    ax.text(-0.02, y+(n-1)/2, m, transform=ax.get_yaxis_transform(),
            ha="right", va="center", fontsize=7.6, fontweight="bold", color=cmap[m])
    y += n + 0.8
ax.set_yticks([]); ax.invert_yaxis()
ax.set_xlabel("meta-analysis log$_2$ fold-change, lesional (PP) vs normal (NN)   [95% CI]", fontsize=8)
ax.set_xlim(-4, 17); ax.margins(y=0.01)
ax.set_xticks([-4, -2, 0, 2, 4, 6, 8, 10, 12])
ax.set_title("IL-1$\\beta$/inflammasome and STAT3 modules in psoriatic lesions\nvs normal skin (peri-lesional excluded; up to 3 studies per gene)",
             fontsize=9, loc="left", pad=10)
ax.text(0.99, -0.045, "filled = FDR<0.05 · faded = n.s. · ‡ high heterogeneity (I²≥75%) · ᵏ⁼² only 2 studies (else 3)",
        transform=ax.transAxes, ha="right", va="top", fontsize=6)
for sp in ["top", "right", "left"]:
    ax.spines[sp].set_visible(False)
fig.subplots_adjust(left=0.19, right=0.985, top=0.92, bottom=0.06)
fig.savefig("fig_il1b_stat3_NNvsPP.png", dpi=300, bbox_inches="tight")