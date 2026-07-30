# Auto-extracted generating script
# Produces: druggability_annotation.csv
# Conda env: python   (run with this environment activated)
# Inputs (expected alongside / in data/): target_universe.csv, ot_targets.json, chembl_targets.json, chembl_bioactivity.json
# Source artifact version: 70d934da-ef2b-45db-982f-d354ae25af31
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import json
import numpy as np
import pandas as pd

sig = pd.read_csv("target_universe.csv")
ot = json.load(open("ot_targets.json"))
chembl = json.load(open("chembl_targets.json"))
bio = json.load(open("chembl_bioactivity.json"))

def tract_flags(t):
    """Summarize Open Targets small-molecule (SM) tractability into interpretable flags + a bucket score."""
    sm = {x["label"]: x["value"] for x in t["tractability"] if x["modality"]=="SM"}
    ab = {x["label"]: x["value"] for x in t["tractability"] if x["modality"]=="AB"}
    # SM clinical precedence bucket
    if sm.get("Approved Drug"): sm_bucket="Approved drug"
    elif sm.get("Advanced Clinical"): sm_bucket="Advanced clinical"
    elif sm.get("Phase 1 Clinical"): sm_bucket="Phase 1"
    elif sm.get("Structure with Ligand") or sm.get("High-Quality Ligand"): sm_bucket="Discovery precedence"
    elif sm.get("High-Quality Pocket") or sm.get("Med-Quality Pocket"): sm_bucket="Druggable pocket"
    elif sm.get("Druggable Family"): sm_bucket="Druggable family"
    else: sm_bucket="Unknown"
    ab_clin = ab.get("Approved Drug") or ab.get("Advanced Clinical") or ab.get("Phase 1 Clinical")
    return dict(
        sm_bucket=sm_bucket,
        sm_has_ligand=bool(sm.get("Structure with Ligand")),
        sm_hq_ligand=bool(sm.get("High-Quality Ligand")),
        sm_hq_pocket=bool(sm.get("High-Quality Pocket")),
        sm_mq_pocket=bool(sm.get("Med-Quality Pocket")),
        sm_druggable_family=bool(sm.get("Druggable Family")),
        ab_clinical=bool(ab_clin),
    )

# ordinal score for SM bucket (clinical precedence weighting)
bucket_rank = {"Approved drug":6,"Advanced clinical":5,"Phase 1":4,"Discovery precedence":3,
               "Druggable pocket":2,"Druggable family":1,"Unknown":0}

rows=[]
for sym, t in ot.items():
    tf = tract_flags(t)
    cls = "; ".join(sorted({c["label"] for c in t.get("targetClass",[]) if c.get("level")=="l1"})) or "Unclassified"
    cb = chembl.get(sym) or {}
    bb = bio.get(sym) or {}
    de = sig[sig.gene==sym]
    de = de.iloc[0] if len(de) else None
    rows.append(dict(
        gene=sym, ensembl=t["id"], approvedName=t.get("approvedName",""),
        protein_class=cls,
        logFC=(de.logFC if de is not None else np.nan),
        FDR=(de.FDR if de is not None else np.nan),
        direction=(de.direction if de is not None else ""),
        evidence_rank=(int(de.evidence_rank) if de is not None else np.nan),
        pathway_support=(int(de.pathway_support) if de is not None else 0),
        in_leading_edge=(bool(de.in_leading_edge) if de is not None else False),
        stat3_target=(bool(de.stat3_target) if de is not None else False),
        ot_sm_bucket=tf["sm_bucket"], ot_sm_bucket_rank=bucket_rank[tf["sm_bucket"]],
        sm_has_ligand=tf["sm_has_ligand"], sm_hq_ligand=tf["sm_hq_ligand"],
        sm_hq_pocket=tf["sm_hq_pocket"], sm_mq_pocket=tf["sm_mq_pocket"],
        sm_druggable_family=tf["sm_druggable_family"], ab_clinical=tf["ab_clinical"],
        ot_drug_candidates=t.get("drugAndClinicalCandidates",{}).get("count",0),
        chembl_target_id=cb.get("target_chembl_id"),
        n_potent_inhib=bb.get("n_potent_pchembl7",0),
        n_active_inhib=bb.get("n_active_pchembl6",0),
    ))
ann = pd.DataFrame(rows)
ann.to_csv("druggability_annotation.csv", index=False)