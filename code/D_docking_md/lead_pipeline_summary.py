# Auto-extracted generating script
# Produces: lead_pipeline_summary.csv
# Conda env: dock-md   (run with this environment activated)
# Inputs (expected alongside / in data/): druggable_target_ranking.csv, pocket_druggability_summary.csv, docking_ranked.csv, validation.json
# Source artifact version: 8e504728-b362-4469-9989-1a84520dba99
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import pandas as pd, json, numpy as np

rank=pd.read_csv("druggable_target_ranking.csv")
pocket=pd.read_csv("pocket_druggability_summary.csv")
dock=pd.read_csv("docking_ranked.csv")
val=json.load(open("validation.json"))

leads=["STAT3","RORC","JAK3"]
rows=[]
for g in leads:
    r=rank[rank.gene==g].iloc[0]
    p=pocket[pocket.gene==g].iloc[0]
    rows.append(dict(
        gene=g, direction=r.direction, meta_logFC=round(r.logFC,2), meta_FDR=f"{r.FDR:.1e}",
        protein_class=r.protein_class, disease_axis=r.disease_axis,
        ot_tractability=r.ot_sm_bucket, chembl_potent_inhib=int(r.n_potent_inhib),
        pdb=p.pdb, domain=p.domain, resolution_A=p.resolution_A,
        fpocket_druggability=p.ligand_pocket_druggability,
        docked=("yes (151-cmpd screen)" if g=="RORC" else "cluster package ready"),
    ))
summary=pd.DataFrame(rows)
summary.to_csv("lead_pipeline_summary.csv", index=False)