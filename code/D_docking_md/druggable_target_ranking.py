# Auto-extracted generating script
# Produces: druggable_target_ranking.csv
# Conda env: python   (run with this environment activated)
# Inputs (expected alongside / in data/): druggability_annotation.csv
# Source artifact version: e404033f-d152-49bb-a29f-a83ce87b108a
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import pandas as pd
import numpy as np

ann = pd.read_csv("druggability_annotation.csv")

# --- Component scores (each 0-1 normalized) ---
def mm(x):
    x=x.astype(float); lo,hi=np.nanmin(x),np.nanmax(x)
    return (x-lo)/(hi-lo) if hi>lo else x*0

# DE evidence: effect size + significance (already have)
ann["s_effect"] = mm(ann.logFC.abs())
ann["s_signif"] = mm(-np.log10(ann.FDR.clip(lower=1e-300)))
# pathway centrality
ann["s_pathway"] = mm(ann.pathway_support.astype(float))
# druggability: tractability bucket rank + chemical tractability (log inhibitors) + drug candidates
ann["s_tract"] = mm(ann.ot_sm_bucket_rank.astype(float))
ann["s_chem"]  = mm(np.log10(ann.n_potent_inhib.clip(lower=0)+1))
ann["s_drugs"] = mm(np.log10(ann.ot_drug_candidates.clip(lower=0)+1))
# structural availability (has ligand/pocket) — key for docking/MD feasibility
ann["s_struct"] = ((ann.sm_has_ligand.astype(int)+ann.sm_hq_ligand.astype(int)
                    +ann.sm_hq_pocket.astype(int)+ann.sm_mq_pocket.astype(int))/4.0)

# Composite (weights: DE 30%, pathway 15%, druggability 40%, structure 15%)
ann["druggability_score"] = (
    0.18*ann.s_effect + 0.12*ann.s_signif +
    0.15*ann.s_pathway +
    0.16*ann.s_tract + 0.12*ann.s_chem + 0.12*ann.s_drugs +
    0.15*ann.s_struct
)
ann = ann.sort_values("druggability_score", ascending=False).reset_index(drop=True)
ann["druggability_rank"] = np.arange(1,len(ann)+1)
ann.to_csv("druggable_target_ranking.csv", index=False)

# Disease-axis relevance
axis_map = {
 "Th17/IL-17 axis": ["RORC","RORA","IL17A","IL17F","IL17RA","CCL20","CCR6","IL23A","IL12B","IL23R"],
 "JAK-STAT": ["JAK1","JAK2","JAK3","TYK2","STAT1","STAT2","STAT3","IL6R","IL22RA1","IL22","IL19"],
 "NF-kB/TNF": ["NFKB1","NFKB2","RELA","RELB","REL","BIRC3","TNIP3","ZC3H12A","IRAK2","TNF"],
 "Innate/antimicrobial": ["S100A9","S100A7","S100A8","S100A12","DEFB4A","DEFB4B","IL36A","IL36G","IL36RN","LTF","PI3","NOS2"],
 "Chemokine/leukocyte": ["CXCR1","CXCR2","CXCL13","CCL2","CCL20","LTB4R"],
 "Tryptophan/immunometab": ["IDO1","KYNU","TDO2","VNN3"],
 "Eicosanoid/protease": ["PTGS2","ALOX5","ALOX15","MMP9","MMP1","MMP12","PLA2G2A","PLA2G4D","LTA4H"],
 "Proliferation": ["PLK1","AURKA","AURKB","CDK1","RRM2","TYMS","BIRC5"],
}
g2axis={}
for ax_,gs in axis_map.items():
    for g in gs: g2axis.setdefault(g,[]).append(ax_)
ann["disease_axis"] = ann.gene.map(lambda g: "; ".join(g2axis.get(g,[])) )
ann["on_psoriasis_axis"] = ann.disease_axis.str.len()>0

# Structural tractability for docking/MD (has a real pocket + ligand)
ann["struct_ready"] = (ann.sm_has_ligand | ann.sm_hq_pocket | ann.sm_mq_pocket)

LEADS = ["STAT3","RORC","JAK3"]
ann["is_lead"] = ann.gene.isin(LEADS)

ann.to_csv("druggable_target_ranking.csv", index=False)