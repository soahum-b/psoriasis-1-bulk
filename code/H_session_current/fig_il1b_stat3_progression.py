# Auto-extracted generating script
# Produces: fig_il1b_stat3_progression.png
# Conda env: python   (run with this environment activated)
# Inputs (expected alongside / in data/): meta_de_PPvsNN.csv, meta_de_PNvsNN.csv, meta_de_PNvsPP.csv, WHITEPAPER.md
# Source artifact version: 193182a7-5452-4117-915a-20e2709e650b
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


contrasts = {
  "PPvsNN": "meta_de_PPvsNN.csv",
  "PNvsNN": "meta_de_PNvsNN.csv",
  "PNvsPP": "meta_de_PNvsPP.csv",
}
dfs = {c: pd.read_csv(v).set_index("gene") for c, v in contrasts.items()}

panel = {
 "IL-1 / inflammasome": ["IL1B","IL1A","IL1RN","NLRP3","NLRP1","AIM2","CASP1","PYCARD","IL18","IL1R1","IL1R2","IL1RAP","NLRC4","CASP5","GSDMD"],
 "IL-36 (family control)": ["IL36G","IL36A","IL36B","IL36RN","IL1F10"],
 "STAT program / JAK": ["STAT3","STAT1","JAK1","JAK2","JAK3","TYK2","SOCS3","IL6","IL6R","IL22","IL22RA1","OSMR"],
 "Th17 axis (context)": ["IL17A","IL17F","IL23A","IL23R","IL12B","RORC","IL17RA","TRAF3IP2"],
 "IL-1 downstream output": ["S100A7","S100A8","S100A9","CXCL1","CXCL2","CXCL8","CCL20","DEFB4A","LCN2"],
}

apply_figure_style()

order_groups = ["STAT program / JAK","IL-1 / inflammasome","IL-36 (family control)",
                "Th17 axis (context)","IL-1 downstream output"]
genes_ord, grp_ord = [], []
for grp in order_groups:
    for g in panel[grp]:
        present = any((g in dfs[c].index) for c in dfs)
        if present:
            genes_ord.append(g); grp_ord.append(grp)

contr_cols = ["PNvsNN","PPvsNN"]
col_labels = ["PN vs NN\n(peri-lesional)","PP vs NN\n(lesional)"]

M = np.full((len(genes_ord), len(contr_cols)), np.nan)
FDR = np.full_like(M, np.nan)
I2 = np.full_like(M, np.nan)
for i,g in enumerate(genes_ord):
    for j,c in enumerate(contr_cols):
        d=dfs[c]
        if g in d.index:
            M[i,j]=d.loc[g,"logFC"]; FDR[i,j]=d.loc[g,"FDR"]; I2[i,j]=d.loc[g,"I2"]

fig, ax = plt.subplots(figsize=(5.2, 11))
vmax=6
im = ax.imshow(M, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
ax.set_xticks(range(len(col_labels))); ax.set_xticklabels(col_labels, fontsize=8)
ax.set_yticks(range(len(genes_ord)))
ax.set_yticklabels([f"$\\it{{{g}}}$" for g in genes_ord], fontsize=7)
ax.xaxis.set_ticks_position("top"); ax.xaxis.set_label_position("top")

def star(p):
    if np.isnan(p): return ""
    return "***" if p<1e-3 else "**" if p<1e-2 else "*" if p<0.05 else "·"
for i in range(len(genes_ord)):
    for j in range(len(contr_cols)):
        if np.isnan(M[i,j]):
            ax.text(j,i,"—",ha="center",va="center",color="0.6",fontsize=7); continue
        val=M[i,j]; txtcol="white" if abs(val)>3 else "black"
        s=star(FDR[i,j])
        ax.text(j,i,f"{val:+.1f}{s}",ha="center",va="center",color=txtcol,fontsize=6.2)

boundaries=[]; last=grp_ord[0]; start=0
for i,gr in enumerate(grp_ord+[None]):
    if gr!=last:
        boundaries.append((start,i-1,last)); start=i; last=gr
for (s,e,name) in boundaries:
    if s>0: ax.axhline(s-0.5,color="0.15",lw=1.2)
    ax.text(1.62, (s+e)/2, name.replace(" / ","/\n").replace(" (","\n("),
            va="center",ha="left",fontsize=6.8,fontweight="bold",rotation=0)

cbar=fig.colorbar(im,ax=ax,fraction=0.038,pad=0.02,location="bottom",shrink=0.7)
cbar.set_label("meta-analysis log$_2$FC (up in →)",fontsize=7)
cbar.ax.tick_params(labelsize=6)
ax.set_title("IL-1$\\beta$/inflammasome and STAT3 modules along the\npsoriasis progression axis (bulk meta-analysis)",
             fontsize=9, pad=28, loc="left")
ax.text(0.0,-0.055,"*** FDR<0.001  ** <0.01  * <0.05  · n.s.   — not tested",
        transform=ax.transAxes,fontsize=6,va="top")
fig.subplots_adjust(left=0.22,right=0.72,top=0.90,bottom=0.09)
fig.savefig("fig_il1b_stat3_progression.png",dpi=300,bbox_inches="tight")
print("saved", M.shape)