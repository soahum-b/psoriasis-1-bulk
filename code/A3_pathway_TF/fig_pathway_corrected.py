# Auto-extracted generating script
# Produces: fig_pathway_corrected.png
# Conda env: python   (run with this environment activated)
# Inputs (expected alongside / in data/): figdat_A.csv, figdat_B.csv, figdat_C.csv, figdat_D.csv
# Source artifact version: 8db6a21c-4017-46ee-a18b-9ef91433789e
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

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

A = pd.read_csv("figdat_A.csv")
B = pd.read_csv("figdat_B.csv")
C = pd.read_csv("figdat_C.csv")
D = pd.read_csv("figdat_D.csv")

# theme -> color (outer hue families). Focal = inflammatory/JAK-STAT3 (project lead)
theme_order=["Proliferation / cell cycle","Interferon / antiviral","RNA processing / translation",
             "Immune effector / antimicrobial","Metabolism / MYC-mTOR",
             "Inflammatory signaling (IL/JAK-STAT/TNF)","Keratinocyte / epidermis","Other"]
palette={"Proliferation / cell cycle":"#4C72B0","Interferon / antiviral":"#55A868",
         "RNA processing / translation":"#8172B3","Immune effector / antimicrobial":"#937860",
         "Metabolism / MYC-mTOR":"#64B5CD","Inflammatory signaling (IL/JAK-STAT/TNF)":"#C44E52",
         "Keratinocyte / epidermis":"#CCB974","Other":"#B9B9B9"}

fig=plt.figure(figsize=(13.2,10.0))
gs=fig.add_gridspec(2,2,hspace=0.42,wspace=0.30,height_ratios=[1,1.05])
axA=fig.add_subplot(gs[0,0]); axB=fig.add_subplot(gs[0,1])
axC=fig.add_subplot(gs[1,0]); axD=fig.add_subplot(gs[1,1])

# ---- Panel A: mega vs meta de-inflation (Hallmark), slope/dumbbell ----
a=A.copy()
a["nlp_mega"]=-np.log10(a["padj_mega"].clip(lower=1e-320))
a["nlp_meta"]=-np.log10(a["padj_meta"].clip(lower=1e-320))
a=a.sort_values("nlp_mega",ascending=True).tail(14).reset_index(drop=True)
y=np.arange(len(a))
for i,r in a.iterrows():
    drop = r["nlp_mega"]-r["nlp_meta"]
    col = "#C44E52" if drop>3 else META_GREY
    axA.plot([r["nlp_meta"],r["nlp_mega"]],[i,i],color=col,lw=1.4,zorder=1,alpha=0.8)
axA.scatter(a["nlp_mega"],y,s=26,color=META_GREY,label="mega-analysis (pooled fit)",zorder=3)
axA.scatter(a["nlp_meta"],y,s=26,color="#4C72B0",label="random-effects meta (corrected)",zorder=3)
axA.axvline(-np.log10(0.05),ls=":",color="#999999",lw=0.9)
axA.set_yticks(y); axA.set_yticklabels(a["label"],fontsize=6)
axA.set_xlabel("$-\\log_{10}$ adjusted p (GSEA)")
axA.set_title("A   Mega-analysis inflated significance;\nmeta ranking corrects it",loc="left",fontsize=8.5)
axA.legend(frameon=False,fontsize=5.8,loc="lower right")
axA.margins(y=0.04); axA.text(0.98,0.02,"red = drops >3 log-units",transform=axA.transAxes,
                              ha="right",va="bottom",fontsize=5.5,color="#C44E52")

# ---- Panel B: corrected CAMERA, direction-consistency ----
b=B.sort_values("Zc").reset_index(drop=True)
yb=np.arange(len(b))
colb=["#4C72B0" if c else "#C44E52" for c in b["direction_consistent"]]
axB.barh(yb,b["Zc"],color=colb,height=0.72,zorder=2)
axB.set_yticks(yb); axB.set_yticklabels(b["label"],fontsize=6)
axB.set_xlabel("combined meta-Z (Stouffer, k=4)")
axB.axvline(0,color="#333",lw=0.8)
axB.set_title("B   Corrected competitive meta (CAMERA)\nno saturated statistics; all 4 studies",loc="left",fontsize=8.5)
axB.legend(handles=[Patch(color="#4C72B0",label="all studies agree (consistent)"),
                    Patch(color="#C44E52",label="mixed direction")],
           frameon=False,fontsize=5.8,loc="lower right")
axB.margins(y=0.03)

# ---- Panel C: non-redundant top pathways, colored by theme ----
c=C.sort_values("nlp",ascending=True).tail(18).reset_index(drop=True)
yc=np.arange(len(c))
cc=[palette.get(t,"#B9B9BF") for t in c["theme"]]
axC.scatter(c["nlp"],yc,s=(c["NES"]*22),color=cc,zorder=3,edgecolor="white",lw=0.4)
for i,r in c.iterrows():
    axC.plot([0,r["nlp"]],[i,i],color=cc[i],lw=0.7,alpha=0.45,zorder=1)
axC.set_yticks(yc); axC.set_yticklabels([l[:34] for l in c["label"]],fontsize=6)
axC.set_xlabel("$-\\log_{10}$ joint-adjusted p")
axC.set_title("C   Non-redundant top pathways\n(collapsePathways + joint BH; size $\\propto$ NES)",loc="left",fontsize=8.5)
axC.margins(y=0.03,x=0.08)

# ---- Panel D: theme composition of full non-redundant set ----
d=D.set_index("theme").reindex([t for t in theme_order if t in D["theme"].values]).reset_index()
yd=np.arange(len(d))[::-1]
axD.barh(yd,d["n"],color=[palette[t] for t in d["theme"]],height=0.72,zorder=2)
for i,(_,r) in enumerate(d.iterrows()):
    axD.text(r["n"]+2,yd[i],str(int(r["n"])),va="center",fontsize=6.5)
axD.set_yticks(yd); axD.set_yticklabels(d["theme"],fontsize=6.2)
axD.set_xlabel("non-redundant significant sets (n=308)")
axD.set_title("D   What the corrected enrichment says\nis reproducibly up in lesional skin",loc="left",fontsize=8.5)
axD.margins(x=0.12)

fig.savefig("fig_pathway_corrected.png",dpi=300,bbox_inches="tight")
print("saved", A.shape, B.shape, C.shape, D.shape)