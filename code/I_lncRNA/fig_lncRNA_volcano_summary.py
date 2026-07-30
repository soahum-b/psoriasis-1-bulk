# Auto-extracted generating script
# Produces: fig_lncRNA_volcano_summary.png
# Conda env: python   (run with this environment activated)
# Inputs (expected alongside / in data/): lncRNA_significant_ranked.csv, meta_de_PPvsNN_4study_biotyped.csv
# Source artifact version: 655fe386-10f6-4697-b4f9-0ac3437dfaca
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

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


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

apply_figure_style(sizes=(8,7,6))

bt = pd.read_csv("meta_de_PPvsNN_4study_biotyped.csv")
rk = pd.read_csv("lncRNA_significant_ranked.csv")
lnc = bt[(bt.is_lncRNA==True) & bt.FDR.notna()].copy()
lnc["neglog10FDR"] = -np.log10(lnc["FDR"].clip(lower=1e-300))
# merge tier onto lnc
lnc = lnc.merge(rk[["gene","confidence"]], on="gene", how="left")
lnc["tier"] = lnc["confidence"].fillna("NS")

# cap y for display; mark true infinities/large separately
YCAP = 60
lnc["y"] = lnc["neglog10FDR"].clip(upper=YCAP)
lnc["capped"] = lnc["neglog10FDR"] > YCAP

tier_color = {"NS":"#c9ccd1","C_significant_but_flagged":"#f0b0a0",
              "B_moderate":"#6ba3d6","A_high":"#c0392b"}
tier_label = {"A_high":"Tier A (k≥3, I²<50%)","B_moderate":"Tier B (k≥3, I²<75%)",
              "C_significant_but_flagged":"Tier C (flagged: low k / high I²)","NS":"Not significant"}
zord = {"NS":1,"C_significant_but_flagged":2,"B_moderate":3,"A_high":4}

fig = plt.figure(figsize=(11,5.2))
gs = fig.add_gridspec(1, 2, width_ratios=[2.15,1.0], wspace=0.30)
ax = fig.add_subplot(gs[0,0])

for t in ["NS","C_significant_but_flagged","B_moderate","A_high"]:
    d = lnc[lnc.tier==t]
    ax.scatter(d.logFC, d.y, s=14 if t!="NS" else 8,
               c=tier_color[t], alpha=0.55 if t=="NS" else 0.8,
               edgecolors="none", zorder=zord[t], label=tier_label[t], rasterized=(t=="NS"))
cap = lnc[lnc.capped]
ax.scatter(cap.logFC, [YCAP]*len(cap), marker="^", s=22,
           c=[tier_color[t] for t in cap.tier], edgecolors="k", linewidths=0.3, zorder=5)
ax.axhline(-np.log10(0.05), ls="--", lw=0.8, c="#555", zorder=0.5)
ax.text(4.8, -np.log10(0.05), "FDR 0.05", va="bottom", ha="right", fontsize=6, c="#555")
ax.axvline(0, ls="-", lw=0.6, c="#999", zorder=0.4)

labA = lnc[lnc.tier=="A_high"].reindex(lnc[lnc.tier=="A_high"].logFC.abs().sort_values(ascending=False).index).head(11)
labA = labA[~labA.gene.isin(["XIST","TSIX"])]
for _,rr in labA.iterrows():
    ax.annotate(rr.gene, (rr.logFC, rr.y), fontsize=5.5, fontstyle="italic",
                xytext=(3,3), textcoords="offset points", zorder=6)
for gene, tx, ty in [("XIST",-3.6,11),("TSIX",-1.6,4.2)]:
    rr = lnc[lnc.gene==gene].iloc[0]
    ax.annotate(rr.gene, (rr.logFC, rr.y), fontsize=6, fontstyle="italic", fontweight="bold",
                xytext=(tx,ty), textcoords="data", zorder=7,
                arrowprops=dict(arrowstyle="-", lw=0.5, color="#c0392b"))
ax.set_xlabel("Meta log$_2$FC (PP vs NN)")
ax.set_ylabel("−log$_{10}$ FDR  (capped at 60)")
ax.set_title("Differentially expressed lncRNAs in psoriatic (PP) vs normal (NN) skin", loc="left", fontsize=8)
ax.margins(0.04)
ax.legend(loc="upper right", frameon=False, fontsize=5.6, handletextpad=0.2, borderpad=0.2)

ax2 = fig.add_subplot(gs[0,1])
tiers = ["A_high","B_moderate","C_significant_but_flagged"]
up = [((rk.confidence==t)&(rk.direction=="up_in_PP")).sum() for t in tiers]
dn = [((rk.confidence==t)&(rk.direction=="down_in_PP")).sum() for t in tiers]
ypos = np.arange(len(tiers))[::-1]
ax2.barh(ypos, up, color="#c0392b", label="up in PP", height=0.6)
ax2.barh(ypos, [-x for x in dn], color="#2c7fb8", label="down in PP", height=0.6)
for i,(u,d) in enumerate(zip(up,dn)):
    y=ypos[i]
    ax2.text(u-6, y, str(u), va="center", ha="right", fontsize=6, color="white", fontweight="bold")
    ax2.text(-d+6, y, str(d), va="center", ha="left", fontsize=6, color="white", fontweight="bold")
ax2.axvline(0, c="k", lw=0.6)
ax2.set_yticks(ypos)
ax2.set_yticklabels(["Tier A\n(k≥3, I²<50%)","Tier B\n(k≥3, I²<75%)","Tier C\n(flagged)"], fontsize=6)
ax2.set_xlabel("Significant lncRNAs (FDR<0.05)")
ax2.set_title("1,866 significant lncRNAs\nby confidence & direction", loc="left", fontsize=8)
ax2.legend(loc="upper right", frameon=False, fontsize=5.6)
ax2.set_xlim(-720, 400)
for s in ["top","right"]: ax2.spines[s].set_visible(False)
for a,l in [(ax,"a"),(ax2,"b")]:
    a.text(-0.08,1.02,l,transform=a.transAxes,fontweight="bold",fontsize=11,va="bottom")

fig.savefig("fig_lncRNA_volcano_summary.png", dpi=300, bbox_inches="tight")