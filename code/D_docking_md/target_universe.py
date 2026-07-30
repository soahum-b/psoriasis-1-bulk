# Auto-extracted generating script
# Produces: target_universe.csv
# Conda env: python   (run with this environment activated)
# Inputs (expected alongside / in data/): Recount-3/meta_de_PPvsNN.csv, Recount-3/clust_module_genes.csv, Recount-3/tf_activity_collectri.csv
# Source artifact version: 0a3b5d42-80ce-4919-b6a2-5b6b8e22e3d1
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import pandas as pd
import numpy as np
import os

os.makedirs("work", exist_ok=True)

# --- Meta-analysis DE (primary evidence base) ---
meta = pd.read_csv('/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3/meta_de_PPvsNN.csv')
# significant DEGs
sig = meta[(meta.FDR<0.05) & (meta.logFC.abs()>1)].copy()

# --- Pathway evidence ---
# leading edge sets need to be reconstructed from R outputs
# Since we don't have the RDS files directly accessible here, we use what's available
# We'll need the leading_edge_long.csv and stat3_targets.txt from the R step
# These were written to /tmp/psor_pathway in the trace
WP = "/tmp/psor_pathway"

le = pd.read_csv(f"{WP}/leading_edge_long.csv")
le_genes = set(le.gene)
le_by_gene = le.groupby("gene")["set"].apply(lambda s: ";".join(sorted(set(s))))
stat3_targets = set(open(f"{WP}/stat3_targets.txt").read().split())
tf = pd.read_csv('/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3/tf_activity_collectri.csv')
# module genes (co-expression module from clustering)
mod = pd.read_csv('/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3/clust_module_genes.csv')
mod_genes = set(mod.gene)

# flags
sig["in_leading_edge"] = sig.gene.isin(le_genes)
sig["leading_edge_sets"] = sig.gene.map(le_by_gene).fillna("")
sig["stat3_target"] = sig.gene.isin(stat3_targets)
sig["in_coexpr_module"] = sig.gene.isin(mod_genes)

# --- Composite evidence score (transparent, rank-based) ---
# effect size: |logFC|; significance: -log10(FDR); robustness: k studies, penalize high I2
sig["absLFC"] = sig.logFC.abs()
sig["neglog10FDR"] = -np.log10(sig.FDR.clip(lower=1e-300))
# robustness: reward more studies (k), penalize heterogeneity (I2 in %)
sig["robustness"] = sig.k * (1 - sig.I2/100.0)
# pathway support count (0-4): LE, STAT3 target, module, TF-axis pathway member
sig["pathway_support"] = (sig.in_leading_edge.astype(int) + sig.stat3_target.astype(int) + sig.in_coexpr_module.astype(int))

def z(x): 
    x=x.astype(float); 
    return (x-x.mean())/x.std(ddof=0)
sig["evidence_score"] = (
    0.35*z(sig.absLFC) + 0.30*z(sig.neglog10FDR) +
    0.15*z(sig.robustness) + 0.20*z(sig.pathway_support)
)
sig = sig.sort_values("evidence_score", ascending=False).reset_index(drop=True)
sig["evidence_rank"] = np.arange(1, len(sig)+1)

sig.to_csv("work/target_universe.csv", index=False)