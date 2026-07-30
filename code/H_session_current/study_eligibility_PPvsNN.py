# Auto-extracted generating script
# Produces: study_eligibility_PPvsNN.csv
# Conda env: python   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 8099b3fa-80df-413a-870b-649451ed2dba
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import pandas as pd
rows = [
 # currently used for PPvsNN
 ("SRP035988","Li 2014 case-control","USED (k=3)","healthy","83 NN + 95 PP","40.8M","clean healthy-vs-lesional; NN=healthy skin"),
 ("SRP165679","Tsoi 2019 AD/PSO","USED (k=3)","healthy + uninvolved","38 NN + 28 PP (+27 PN)","31.8M","has BOTH control types → anchors moderator"),
 ("SRP126422","Keratinocyte+biopsy","USED (k=3)","healthy","4 NN + 4 PP","26.8M","tiny but real skin"),
 # NEW clean
 ("SRP065812","lncRNA psor+healthy","ADD ✓","healthy","16 NN + 18 PP*","43.4M","*pre-adalimumab only; 18 post-tx dropped"),
 # judgment call
 ("ERP110816","QMUL etanercept","JUDGMENT ⚠","uninvolved only","86 UN + 92 PP","deep","ALL on etanercept (anti-TNF) — confounds TNF→IL-1→STAT3"),
 # excluded
 ("SRP026042","AhR agonist culture","EXCLUDE ✗","(in-vitro)","cultured explants","—","ex-vivo drug-perturbed; not comparable protocol"),
 ("SRP049599","JunB keratinocyte","EXCLUDE ✗","(in-vitro)","foreskin KC line","—","cell line + siRNA, no lesional tissue"),
 ("SRP093255","skin MSC","EXCLUDE ✗","(MSC)","6 MSC","—","mesenchymal stem cells, not skin biopsy"),
 ("SRP076982","anatomic plaques","EXCLUDE ✗","none","211 PP + 48 PN, 0 NN","4.0M","no normal + very shallow (~4M)"),
 ("SRP116922/130972/086612","T-cell / KC in-vitro","EXCLUDE ✗","(in-vitro)","—","—","explant/culture, not whole-skin"),
 ("SRP173379/173378/132160","GPP blood/PBMC/neut","EXCLUDE ✗","(blood)","—","—","blood compartment, not skin"),
]
E = pd.DataFrame(rows, columns=["study","label","decision","control_type","n (NN/PP)","depth","note"])
E.to_csv("study_eligibility_PPvsNN.csv", index=False)
print(E.to_string(index=False))