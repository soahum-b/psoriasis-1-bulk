# Auto-extracted generating script
# Produces: plan_scissor-backbone-tier-1-ordinal-nn-pn-pp_689f69d8.json
# Conda env: unknown   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 08cb3291-6a22-456c-949b-1158636db681
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import json

plan = {
    "title": "Scissor-on-Gradient Backbone: Tier 1 Ordinal NN < PN < PP",
    "dataset": {
        "geo_accession": "GSE173706",
        "paper": "Ma et al. 2023, Nat Commun",
        "doi": "10.1038/s41467-023-39020-4",
        "raw_tar_url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE173nnn/GSE173706/suppl/GSE173706_RAW.tar",
        "raw_tar_size_mb": 264.1,
        "n_samples": 33,
        "sample_tiers": {
            "NS_normal_healthy": {"n": 8, "bulk_tier": "NN", "ordinal": 0},
            "PN_perilesional_uninvolved": {"n": 11, "bulk_tier": "PN", "ordinal": 1},
            "PP_lesional_psoriatic": {"n": 14, "bulk_tier": "PP", "ordinal": 2}
        },
        "format": "per_sample_csv_gz_count_matrices",
        "requires_cellranger": False
    },
    "bulk_anchor": {
        "source": "recount3_7study_meta_analysis",
        "phenotype": "ordinal_biopsy_site_label",
        "phenotype_values": {"NN": 0, "PN": 1, "PP": 2},
        "scissor_mode": "Gaussian",
        "circularity_tier": "Tier1_non_circular_by_construction"
    },
    "compute_assessment": {
        "needs_cluster": False,
        "reason": "264MB pre-quantified count matrices; Gaussian Scissor against ~100 bulk samples fits in local 64GB RAM",
        "cluster_note": "Pipeline structured as portable git repo; cluster dispatch possible via git pull + env build if permutation null scale-up needed"
    },
    "phases": [
        "Phase1: Download GSE173706_RAW.tar and parse sample metadata",
        "Phase2: Load 33 count matrices into Seurat, QC, normalize, integrate (Harmony)",
        "Phase3: Annotate cell types (keratinocyte, fibroblast, myeloid subtypes)",
        "Phase4: Run Gaussian Scissor with ordinal phenotype y={0,1,2} from bulk anchor",
        "Phase5: Identify intermediate-peaking cell states (PN-enriched)",
        "Phase6: Back-project to bulk signature, close loop against disease severity"
    ],
    "hypothesis": "Intermediate cell states (partial keratinocyte differentiation arrest, early-activation fibroblasts, pre-lesional myeloid activation) peak in PN and are missed by binary PP-vs-NN contrasts"
}

with open("plan_scissor-backbone-tier-1-ordinal-nn-pn-pp_689f69d8.json", "w") as f:
    json.dump(plan, f, indent=2)