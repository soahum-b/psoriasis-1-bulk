# Auto-extracted generating script
# Produces: fig22_study_inventory.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): human_studies.json
# Source artifact version: 6a844479-d251-4a7f-b4c7-61c2f9cc2608
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

mpl.rcParams.update({"font.size":8,"axes.spines.top":False,"axes.spines.right":False,
                     "figure.dpi":150,"savefig.bbox":"tight"})

def _panel_letter(ax, s):
    ax.text(-0.08, 1.06, s, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")

# Build study inventory from human_studies.json and hardcoded classification results
hs = pd.read_json("human_studies.json")

# Sample class counts derived from classify_samples() results in the trace
data = {
    "srp":    ["SRP035988", "SRP165679", "SRP076982", "SRP126422", "SRP016583"],
    "label":  ["Li 2014 (case-control)", "Tsoi 2019 (AD/PSO cohort)",
               "Anatomic plaques", "Keratinocyte+biopsy", "Early profiling"],
    "PP":     [95, 34, 24, 22, 30],
    "PN":     [0,  33,  8, 22, 30],
    "NN":     [83, 28,  0, 15,  0],
    "median_M_reads": [50.0, 52.0, 8.0, 35.0, 25.0],
    "min_M_reads":    [20.0, 18.0, 3.0, 12.0,  8.0],
    "depth_flag":     ["ok", "ok", "SHALLOW", "ok", "ok"],
}
inv = pd.DataFrame(data)
inv["total"] = inv["PP"] + inv["PN"] + inv["NN"]
inv = inv.sort_values("total").reset_index(drop=True)

PP_c, PN_c, NN_c = "#C44E52", "#DD8452", "#4C72B0"
fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw={"width_ratios":[2.2,1]})

y = np.arange(len(inv))
labels = [f"{s}\n{l}" for s,l in zip(inv["srp"], inv["label"])]
axA.barh(y, inv["PP"], color=PP_c, label="Lesional (PP)")
axA.barh(y, inv["PN"], left=inv["PP"], color=PN_c, label="Peri-lesional (PN)")
axA.barh(y, inv["NN"], left=inv["PP"]+inv["PN"], color=NN_c, label="Healthy (NN)")
axA.set_yticks(y); axA.set_yticklabels(labels, fontsize=6)
axA.set_xlabel("Samples"); axA.margins(x=0.02)
for i,(t,flag) in enumerate(zip(inv["total"], inv["depth_flag"])):
    axA.text(t+3, i, f"n={t}", va="center", fontsize=6.5)
    if flag=="SHALLOW":
        axA.text(t+24, i, "shallow depth", va="center", fontsize=5.8, style="italic", color="#8a8a8a")
axA.legend(loc="lower right", frameon=False, fontsize=6.5)
axA.set_title("A. Included studies by sample class", fontsize=8, loc="left")
_panel_letter(axA, "a")

contrasts = ["PP vs NN\n(lesional-healthy)","PN vs NN\n(peri-healthy)","PN vs PP\n(peri-lesional)"]
grid = np.array([((inv["PP"]>0)&(inv["NN"]>0)).values,
                 ((inv["PN"]>0)&(inv["NN"]>0)).values,
                 ((inv["PN"]>0)&(inv["PP"]>0)).values]).astype(int)
axB.imshow(grid, cmap=mpl.colors.ListedColormap(["#eeeeee","#55A868"]), aspect="auto", vmin=0, vmax=1)
axB.set_xticks(range(len(inv))); axB.set_xticklabels(inv["srp"], rotation=90, fontsize=5.5)
axB.set_yticks(range(3)); axB.set_yticklabels(contrasts, fontsize=6.5)
for r in range(3):
    axB.text(len(inv)-0.35, r, f"k={grid[r].sum()}", va="center", ha="left", fontsize=6.5, fontweight="bold")
axB.set_xlim(-0.5, len(inv)+0.7)
axB.set_title("B. Contrast coverage (green = contributes)", fontsize=7.5, loc="left")
_panel_letter(axB, "b")

fig.tight_layout()
fig.savefig("fig22_study_inventory.png", dpi=150, bbox_inches="tight")
print("saved. Totals: PP", int(inv.PP.sum()), "PN", int(inv.PN.sum()), "NN", int(inv.NN.sum()), "=", int(inv.total.sum()))