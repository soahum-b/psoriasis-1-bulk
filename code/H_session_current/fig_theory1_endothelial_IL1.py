# Auto-extracted generating script
# Produces: fig_theory1_endothelial_IL1.png
# Conda env: python   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 93ad7fe6-84c0-44d0-bdd3-a875cd9cad74
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


def panel_letter(ax, letter, dx=-0.18, dy=1.02, case="lower", fontsize=None):
    import matplotlib.pyplot as plt
    if fontsize is None:
        fontsize = plt.rcParams.get("font.size", 8) + 1
    s = letter.lower() if case == "lower" else letter.upper()
    ax.text(dx, dy, s, transform=ax.transAxes,
            fontweight="bold", fontsize=fontsize, va="bottom", ha="left")


import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
apply_figure_style()

progA = ["STAT3","IL6","NFKB1","CASP1","GSDMD","PYCARD"]
il1r1_pos = [0.656,0.310,0.281,0.195,0.509,0.382]
il1r1_neg = [0.398,0.126,0.128,0.097,0.304,0.188]
frac_pos = 637/(637+934); frac_neg = 868/(868+2951)
OR_B, p_B = 2.32, 1.36e-38

genesC = ["IL1RN","VCAM1","PYCARD","IL1B","GSDMD","IL1R2","CASP1","CXCL8",
          "STAT3","ICAM1","NFKB1","CCL2","SELE","IL1RAP","IL1R1","IL6"]
lfcC =   [4.10,0.79,1.08,0.31,0.17,0.02,0.53,0.65,-0.29,-0.44,-0.44,-0.64,-0.66,-0.62,-1.12,-1.39]
qC =     [0.017,0.216,0.009,0.375,0.375,0.088,0.764,0.503,1e-4,1e-4,1e-4,1e-4,0.070,0.375,1e-4,1e-4]
sigC=[q<0.05 for q in qC]

fig = plt.figure(figsize=(13,5.2))
gs = fig.add_gridspec(1,3, width_ratios=[1.15,0.8,1.25], wspace=0.5)

axA = fig.add_subplot(gs[0,0])
x=np.arange(len(progA)); w=0.38
axA.bar(x-w/2, np.array(il1r1_pos)*100, w, label="IL1R1⁺ endothelial", color="#c0392b")
axA.bar(x+w/2, np.array(il1r1_neg)*100, w, label="IL1R1⁻ endothelial", color="#95a5a6")
axA.set_xticks(x); axA.set_xticklabels(progA, rotation=45, ha="right")
axA.set_ylabel("% of cells expressing")
axA.set_title("Downstream program tracks with the\nIL-1 receptor on endothelium")
axA.legend(fontsize=7, frameon=False, loc="upper right")
panel_letter(axA,"A")

axB = fig.add_subplot(gs[0,1])
axB.bar([0,1],[frac_neg*100,frac_pos*100], color=["#95a5a6","#c0392b"], width=0.6)
axB.set_xticks([0,1]); axB.set_xticklabels(["IL1R1⁻","IL1R1⁺"])
axB.set_ylabel("% gradient-tracking (Scissor+)")
axB.set_title(f"IL1R1⁺ endothelium enriched\nin gradient set\n(OR {OR_B}, p={p_B:.0e})")
panel_letter(axB,"B")

axC = fig.add_subplot(gs[0,2])
order=np.argsort(lfcC)
yg=[genesC[i] for i in order]; yl=[lfcC[i] for i in order]; ys=[sigC[i] for i in order]
colors=[("#2471a3" if v<0 else "#c0392b") if s else "#d5dbdb" for v,s in zip(yl,ys)]
axC.barh(np.arange(len(yg)), yl, color=colors)
axC.set_yticks(np.arange(len(yg))); axC.set_yticklabels(yg, fontsize=8)
axC.axvline(0,color="k",lw=0.8)
axC.set_xlabel("log₂FC  (lesional PP vs normal NN)")
axC.set_title("Endothelial IL-1 axis + atherosclerosis\nmarkers, PP vs NN (BH<0.05 = colored)")
axC.margins(x=0.15)
panel_letter(axC,"C")

fig.suptitle("IL-1-responsive vascular endothelium carries the STAT3 / inflammasome program in psoriasis",
             fontsize=11, y=1.03)
fig.savefig("fig_theory1_endothelial_IL1.png", dpi=200, bbox_inches="tight")