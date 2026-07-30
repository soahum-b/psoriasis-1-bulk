# Auto-extracted generating script
# Produces: plan_druggability-triage-of-psoriasis-meta-de_bf412257.json
# Conda env: unknown   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 074fd09c-aab3-45e6-9097-e5daf7780210
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import json

plan = {
    "project": "Psoriasis meta-analysis — druggability triage",
    "evidence_base": "meta_de_PPvsNN.csv (23,409 genes; FDR<0.05 & |logFC|>1 → 2,379 candidates; k≥2 studies)",
    "headline_target": {
        "gene": "STAT3",
        "logFC": 1.246,
        "FDR": 2.20e-7,
        "k": 3,
        "I2": 91.1,
        "direction": "up"
    },
    "phases": [
        {
            "phase": 1,
            "name": "Breadth-first druggability triage",
            "inputs": "2,379 meta-DE genes",
            "tools": ["Open Targets (tractability buckets, known drugs)", "ChEMBL (bioactivities, mechanisms)", "PubChem / BindingDB", "FDA approvals", "ClinicalTrials.gov"],
            "output": "Ranked druggable-target shortlist (top ~10–20)"
        },
        {
            "phase": 2,
            "name": "Structure retrieval / prediction",
            "inputs": "Top-ranked targets from Phase 1",
            "tools": ["RCSB PDB", "AlphaFold DB"],
            "output": "PDB files for each target; pocket annotation"
        },
        {
            "phase": 3,
            "name": "Virtual screen — purchasable library",
            "inputs": "Top target pocket(s); ZINC22 purchasable compounds",
            "tools": ["ZINC22 (CartBlanche22)", "AutoDock Vina / smina (CPU docking)"],
            "output": "Ranked pose list; top-N hits by docking score"
        },
        {
            "phase": 4,
            "name": "MD simulation package preparation",
            "inputs": "Top docked pose(s) from Phase 3",
            "tools": ["OpenMM / GROMACS input preparation", "AMBER/CHARMM force fields"],
            "note": "No GPU on this machine — packages prepared for cluster execution",
            "output": "Ready-to-run cluster MD packages (topology, coordinates, run scripts)"
        }
    ],
    "compute_note": "No GPU / no remote compute configured locally. Docking uses CPU-friendly AutoDock Vina/smina. MD simulation packages will be prepared here and executed on user's cluster.",
    "de_filter": {
        "FDR_threshold": 0.05,
        "logFC_abs_threshold": 1.0,
        "min_studies_k": 2,
        "n_passing": 2379
    }
}

with open("plan_druggability-triage-of-psoriasis-meta-de_bf412257.json", "w") as f:
    json.dump(plan, f, indent=2)