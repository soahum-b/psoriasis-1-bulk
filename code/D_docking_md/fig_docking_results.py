# Auto-extracted generating script
# Produces: fig_docking_results.png
# Conda env: dock-md   (run with this environment activated)
# Inputs (expected alongside / in data/): screening_library.csv
# Source artifact version: 2a59b872-0d59-4be7-aee2-5cd6fb824649
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from rdkit import Chem

# Global state
mpl.rcParams.update({"font.size":8,"axes.linewidth":0.6,"font.family":"DejaVu Sans"})

# Load screening library
lib = pd.read_csv("screening_library.csv")

# Recreate docking results - need to compute from the data
# The docking ranked CSV was generated from dock_results.json merged with lib
# Since we don't have the docking results artifact, we need to reconstruct
# the docking_ranked.csv data. However, since dock_results.json is not available
# as an artifact, we use the screening_library as the base and note that
# docking_ranked.csv was saved in the trace.

# Load docking ranked results - reconstruct from the trace logic
# dock_results.json -> merged with lib -> valid (sorted by affinity)
# Since docking_ranked.csv was saved as work/docking/docking_ranked.csv
# and is not listed as an artifact, we must recreate from screening_library + docking scores
# The trace shows the docking was run via run_dock.py and results saved.
# We'll load the screening library and simulate the structure of valid df
# using what was saved.

import json
import os

# Try to load docking results if available
dock_results_path = "work/docking/dock_results.json"
docking_ranked_path = "work/docking/docking_ranked.csv"

if os.path.exists(docking_ranked_path):
    valid = pd.read_csv(docking_ranked_path)
else:
    res = json.load(open(dock_results_path))
    dock = pd.DataFrame(res)
    dock = dock.merge(lib, left_on="ligand", right_on="chembl_id", how="left")
    valid = dock[dock.affinity.notna()].copy().sort_values("affinity").reset_index(drop=True)
    valid["rank"] = np.arange(1, len(valid)+1)
    valid["heavy_atoms"] = valid.smiles.apply(
        lambda s: Chem.MolFromSmiles(s).GetNumHeavyAtoms() if Chem.MolFromSmiles(s) else np.nan)
    valid["lig_efficiency"] = (valid.affinity / valid.heavy_atoms).round(3)

fig = plt.figure(figsize=(11, 4.4))
gs = fig.add_gridspec(1, 3, wspace=0.32)
# (a) affinity distribution
ax = fig.add_subplot(gs[0])
ax.hist(valid.affinity, bins=22, color="#2c7fb8", edgecolor="white", lw=0.4)
ax.axvline(valid.affinity.median(), color="#d73027", lw=1.2, ls="--", label=f"median {valid.affinity.median():.1f}")
ax.axvline(-10, color="0.3", lw=0.8, ls=":", label="−10 kcal/mol")
ax.set_xlabel("Vina affinity (kcal/mol)"); ax.set_ylabel("compounds"); ax.legend(frameon=False, fontsize=6.5)
ax.set_title("(a) Docking score distribution", fontsize=8.2)
# (b) affinity vs QED (are top hits also drug-like?)
ax = fig.add_subplot(gs[1])
sc = ax.scatter(valid.affinity, valid.qed, c=valid.logp, cmap="viridis", s=22, edgecolors="k", lw=0.3)
ax.set_xlabel("Vina affinity (kcal/mol)"); ax.set_ylabel("QED")
top10 = valid.head(10)
ax.scatter(top10.affinity, top10.qed, facecolors="none", edgecolors="#d73027", s=55, lw=1.2, label="top 10")
ax.legend(frameon=False, fontsize=6.5, loc="lower left")
cb = fig.colorbar(sc, ax=ax, pad=0.02); cb.set_label("cLogP", fontsize=7); cb.ax.tick_params(labelsize=6)
ax.set_title("(b) Affinity vs drug-likeness", fontsize=8.2)
# (c) ligand efficiency vs MW
ax = fig.add_subplot(gs[2])
ax.scatter(valid.mw, valid.lig_efficiency, s=22, c="0.6", edgecolors="k", lw=0.3)
ax.scatter(top10.mw, top10.lig_efficiency, c="#d73027", s=30, edgecolors="k", lw=0.4, label="top 10")
ax.set_xlabel("Molecular weight (Da)"); ax.set_ylabel("Ligand efficiency (kcal/mol/heavy atom)")
ax.legend(frameon=False, fontsize=6.5); ax.set_title("(c) Ligand efficiency", fontsize=8.2)
fig.suptitle("RORγt LBD virtual screen (AutoDock Vina, 151 compounds) — CPU demo", fontsize=9, y=1.02)
fig.savefig("fig_docking_results.png", dpi=190, bbox_inches="tight")
print("saved docking results figure")
print("top-10 mean affinity:", round(valid.head(10).affinity.mean(), 2), "| all mean:", round(valid.affinity.mean(), 2))