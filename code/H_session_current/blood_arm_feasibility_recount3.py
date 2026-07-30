# Auto-extracted generating script
# Produces: blood_arm_feasibility_recount3.csv
# Conda env: python   (run with this environment activated)
# Inputs (expected alongside / in data/): recount3_selection_2026-07-05 20_53_44.267187.csv
# Source artifact version: 45450d34-5da3-4142-8d9c-eb9e16adbe6d
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import pandas as pd

blood = pd.DataFrame([
 ("SRP173379","Neutrophils (sorted)","GPP","8 GPP / 11 healthy",19,"Sorted neutrophils — not PBMC; GPP subtype","Marginal"),
 ("SRP173378","Whole blood","GPP","9 GPP / 7 healthy",16,"Whole blood — mixed; GPP subtype","Marginal"),
 ("SRP132160","PBMC","GPP","GPP only, acitretin-treated",15,"PBMC but NO healthy controls + drug-treated","Not usable alone"),
], columns=["study","tissue","disease_subtype","groups","n","note","verdict"])
blood.to_csv("blood_arm_feasibility_recount3.csv", index=False)