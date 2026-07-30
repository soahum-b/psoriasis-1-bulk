# Auto-extracted generating script
# Produces: docking_ranked.csv
# Conda env: dock-md   (run with this environment activated)
# Inputs (expected alongside / in data/): screening_library.csv, dock_results.json
# Source artifact version: f97beaaa-b8c6-4bcd-91fb-3669eb9bb4e3
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import json
import pandas as pd
import numpy as np
from rdkit import Chem

res = json.load(open("dock_results.json"))
dock = pd.DataFrame(res)
lib = pd.read_csv("screening_library.csv")
dock = dock.merge(lib, left_on="ligand", right_on="chembl_id", how="left")
valid = dock[dock.affinity.notna()].copy().sort_values("affinity").reset_index(drop=True)
valid["rank"] = np.arange(1, len(valid)+1)
valid["heavy_atoms"] = valid.smiles.apply(lambda s: Chem.MolFromSmiles(s).GetNumHeavyAtoms() if Chem.MolFromSmiles(s) else np.nan)
valid["lig_efficiency"] = (valid.affinity / valid.heavy_atoms).round(3)
valid.to_csv("docking_ranked.csv", index=False)