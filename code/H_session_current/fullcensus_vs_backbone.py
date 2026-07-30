# Auto-extracted generating script
# Produces: fullcensus_vs_backbone.csv
# Conda env: python   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: ea8e21a8-1d8d-44c4-858b-fab6c62f92c8
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import pandas as pd

comp = pd.DataFrame([
 ("Cells analyzed","20,023 (subset)","89,058 (full census)"),
 ("Scissor+ / Scissor−","1,574 / 1,397","6,837 / 6,228"),
 ("% selected","14.84%","14.67%"),
 ("Reliability CV-MSE (real vs null)","0.147 vs 0.779","0.135 vs ~0.835"),
 ("Endothelial enrichment (OR)","11.26","5.24"),
 ("Endothelial Scissor+ fraction","(top lineage)","0.220 vs 0.051 bg"),
 ("STAT3 in Scissor+ program","log2FC 0.43, padj 0.018 (sig)","log2FC 0.30, p=0.12, padj=1 (n.s.)"),
 ("Dominant gradient-tracker","Endothelial","Endothelial (confirmed)"),
], columns=["metric","20k_backbone","full_census"])
comp.to_csv("fullcensus_vs_backbone.csv", index=False)

# insert a selected-alpha row near the top if not present
if "Selected alpha" not in comp["metric"].values:
    row = pd.DataFrame([("Selected alpha (tuning)","0.40","0.20")], columns=comp.columns)
    comp = pd.concat([comp.iloc[:1], row, comp.iloc[1:]], ignore_index=True)
comp.to_csv("fullcensus_vs_backbone.csv", index=False)
print(comp.to_string(index=False))