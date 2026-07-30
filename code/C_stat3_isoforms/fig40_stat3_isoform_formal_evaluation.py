# Auto-extracted generating script
# Produces: fig40_stat3_isoform_formal_evaluation.png
# Conda env: python   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 40730fd4-670d-4558-a050-330705afcd54
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

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


def set_frame(ax, style="open"):
    show = {"open": (False, False, True, True),
            "boxed": (True, True, True, True),
            "none": (False, False, False, False)}[style]
    for side, vis in zip(("top", "right", "bottom", "left"), show):
        ax.spines[side].set_visible(vis)
        if vis:
            ax.spines[side].set_linewidth(0.6)
    ax.tick_params(direction="out", length=0 if style == "none" else 3, width=0.6)


def panel_letter(ax, letter, dx=-0.18, dy=1.02, case="lower", fontsize=None):
    import matplotlib.pyplot as plt
    if fontsize is None:
        fontsize = plt.rcParams.get("font.size", 8) + 1
    s = letter.lower() if case == "lower" else letter.upper()
    ax.text(dx, dy, s, transform=ax.transAxes,
            fontweight="bold", fontsize=fontsize, va="bottom", ha="left")


apply_figure_style()

# ---- Data ----
studies = ["SRP126422\n(4v4)", "SRP154474\n(n=17)", "SRP165679\n(deep, n=66)", "SRP035988\n(anchor, n=95v83)"]
shift   = [ 1.008, 2.209, -0.071, 1.167]
wp      = [ 0.86,  0.37,   0.70,  0.017]
pooled = [("DL random-effects", 0.868, 0.147, 1.588, 0.018),
          ("HKSJ (few-study honest)", 0.804, -0.440, 2.048, 0.132)]

ladder = [("Meta: DL random-effects", 0.018),
          ("Stouffer weighted-Z (\u221an)", 0.0446),
          ("Meta: HKSJ", 0.132),
          ("Stouffer unweighted", 0.124),
          ("Mega: depth-wtd linear", 0.211),
          ("Mega: quasibinomial (268 samp.)", 0.228),
          ("Fisher (unsigned)", 0.193)]
ladder.sort(key=lambda x: x[1])

loo = [("drop SRP165679", 0.0113),
       ("drop SRP126422", 0.0449),
       ("drop SRP154474", 0.0655),
       ("drop SRP035988 (anchor)", 0.9133)]
loo.sort(key=lambda x: x[1])

grp = ["NN\n(normal)", "PN\n(uninvolved)", "PP\n(lesional)"]
psi_beta = [8.85, 7.69, 8.78]
gene_lfc = [0.0, 0.156, 1.21]

FOCAL, MUTE, ALARM = "#1f4e79", "#9aa7b5", "#c1432c"

fig = plt.figure(figsize=(11.5, 8.6))
gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.34,
                      left=0.20, right=0.90, top=0.90, bottom=0.09)

axA = fig.add_subplot(gs[0, 0])
yA = np.arange(len(studies))
for i,(s,p) in enumerate(zip(shift, wp)):
    c = FOCAL if p < 0.05 else MUTE
    axA.scatter(s, yA[i], s=70, color=c, zorder=3, edgecolor="white", linewidth=0.8)
    axA.annotate(f"p={p:.3g}", (s, yA[i]), xytext=(0,9), textcoords="offset points",
                 fontsize=6, color="#333", ha="center")
yp = [-1.4, -2.2]
for (lab,est,lb,ub,p),yy in zip(pooled, yp):
    c = FOCAL if p < 0.05 else "#6d7f95"
    axA.plot([lb,ub],[yy,yy], color=c, lw=2, zorder=2)
    axA.scatter([est],[yy], marker="D", s=90, color=c, zorder=3, edgecolor="white", linewidth=0.8)
    axA.annotate(f"{lab}: +{est:.2f} pp, p={p:.3g}", (est,yy), xytext=(0,-13),
                 textcoords="offset points", fontsize=6, ha="center", color=c)
axA.axvline(0, color="#555", lw=0.9, ls="--", zorder=1)
axA.axhline(-0.6, color="#ccc", lw=0.8)
axA.set_yticks(list(yA)+yp)
axA.set_yticklabels(studies+["pooled (DL)","pooled (HKSJ)"], fontsize=6.5)
axA.set_xlabel("PP\u2212NN shift in STAT3\u03b2 PSI (percentage points)")
axA.set_title("Per-study \u03b2 shift: one cohort carries the signal", fontsize=8, loc="left")
axA.set_xlim(-0.7, 2.7); axA.set_ylim(-3.0, len(studies)-0.1)
set_frame(axA); panel_letter(axA, "a")

axB = fig.add_subplot(gs[0, 1])
labs=[l for l,_ in ladder]; ps=[p for _,p in ladder]; yB=np.arange(len(labs))[::-1]
cols=[FOCAL if p<0.05 else MUTE for p in ps]
axB.hlines(yB, 0.007, ps, color=cols, lw=1.2, alpha=0.6)
axB.scatter(ps, yB, s=60, color=cols, zorder=3, edgecolor="white", linewidth=0.7)
axB.axvline(0.05, color=ALARM, lw=1.1, ls="--")
axB.annotate("p = 0.05", (0.05, len(labs)-0.5), color=ALARM, fontsize=6.5, ha="center")
axB.set_xscale("log"); axB.set_xlim(0.007, 0.4)
axB.set_xticks([0.01,0.05,0.1,0.2]); axB.set_xticklabels(["0.01","0.05","0.1","0.2"])
axB.set_yticks(yB); axB.set_yticklabels(labs, fontsize=6.3)
axB.set_xlabel("p-value for PP\u2212NN \u03b2 shift")
axB.set_title("Significance depends on the estimator", fontsize=8, loc="left")
set_frame(axB); panel_letter(axB, "b")

axC = fig.add_subplot(gs[1, 0])
labs=[l for l,_ in loo]; ps=[p for _,p in loo]; yC=np.arange(len(labs))[::-1]
cols=[FOCAL if p<0.05 else MUTE for p in ps]
axC.hlines(yC, 0.007, ps, color=cols, lw=1.2, alpha=0.6)
axC.scatter(ps, yC, s=60, color=cols, zorder=3, edgecolor="white", linewidth=0.7)
for p,y in zip(ps,yC):
    axC.annotate(f"{p:.3g}", (p,y), xytext=(7,0), textcoords="offset points",
                 fontsize=6, va="center", color="#333")
axC.axvline(0.05, color=ALARM, lw=1.1, ls="--")
axC.set_xscale("log"); axC.set_xlim(0.007, 1.6)
axC.set_xticks([0.01,0.05,0.1,0.5,1.0]); axC.set_xticklabels(["0.01","0.05","0.1","0.5","1.0"])
axC.set_yticks(yC); axC.set_yticklabels(labs, fontsize=6.5)
axC.set_xlabel("Combined p (Stouffer weighted-Z), each dropping one study")
axC.set_title("Drop the anchor \u2192 signal collapses (p=0.91)", fontsize=8, loc="left")
set_frame(axC); panel_letter(axC, "c")

axD = fig.add_subplot(gs[1, 1])
xD=np.arange(3)
axD.plot(xD, psi_beta, "-o", color=MUTE, lw=2, ms=7, zorder=3, markeredgecolor="white")
axD.set_ylim(0, 10); axD.set_ylabel("STAT3\u03b2 PSI (%)", color="#57687a")
axD.tick_params(axis="y", labelcolor="#57687a")
axD.set_xticks(xD); axD.set_xticklabels(grp, fontsize=6.8); axD.set_xlim(-0.25, 2.25)
axD2=axD.twinx()
axD2.plot(xD, gene_lfc, "-s", color=FOCAL, lw=2.4, ms=7, zorder=4, markeredgecolor="white")
axD2.set_ylim(-0.15, 1.45); axD2.set_ylabel("STAT3 gene log$_2$FC vs NN", color=FOCAL)
axD2.tick_params(axis="y", labelcolor=FOCAL)
axD2.spines["top"].set_visible(False); axD.spines["top"].set_visible(False)
axD.annotate("isoform ratio: flat", (1.0, psi_beta[1]), xytext=(0,-20),
             textcoords="offset points", ha="center", color="#57687a", fontsize=6.8)
axD2.annotate("gene: +1.2 log$_2$FC", (2, gene_lfc[2]), xytext=(-6,-10),
              textcoords="offset points", ha="right", color=FOCAL, fontsize=6.8)
axD.set_title("Disease signal is gene up-regulation, not \u03b2 switch", fontsize=8, loc="left")
panel_letter(axD, "d")

fig.suptitle("Formal evaluation of the STAT3 \u03b1/\u03b2 isoform ratio in psoriasis (lesional PP vs normal NN)",
             fontsize=9.5, x=0.20, ha="left", y=0.965)
fig.savefig("fig40_stat3_isoform_formal_evaluation.png", dpi=300, bbox_inches="tight")