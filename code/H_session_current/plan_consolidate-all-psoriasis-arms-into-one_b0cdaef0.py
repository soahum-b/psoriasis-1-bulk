# Auto-extracted generating script
# Produces: plan_consolidate-all-psoriasis-arms-into-one_b0cdaef0.json
# Conda env: unknown   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 6176d94e-2254-4c11-a403-ee13e52cee92
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import json

plan = {
    "project": "psoriasis-1-bulk",
    "goal": "Arrange the whole project as a reproducible project — consolidate all arms (bulk meta-analysis, Scissor, docking/MD, STAT3 isoforms, lncRNA, genome-browser, presentation) into the canonical cluster repo",
    "repo": "github.com/soahum-b/psoriasis-1-bulk",
    "cluster": "n013",
    "repo_path": "/net/dali/home/mscbio/sba50/Soahum/Project/psoriasis/psoriasis-1-bulk",
    "current_head": "115e663",
    "rule": "additive only — do not move, rename, or overwrite pre-existing files unless explicitly approved",
    "arms": [
        {
            "id": "A",
            "name": "Bulk meta-analysis (backbone)",
            "frame": "db606fdf",
            "n_artifacts": 104,
            "status": "SOLID — in repo",
            "key_content": "38 png, 31 rds, 16 csv, 11 md, 2 R"
        },
        {
            "id": "B",
            "name": "Scissor gradient single-cell",
            "frame": "689f69d8",
            "n_artifacts": 52,
            "status": "SOLID — in repo",
            "key_content": "16 rds, 9 R, 9 png, 9 csv"
        },
        {
            "id": "C",
            "name": "Theory-1: IL1R1+ endothelium / psoriatic march",
            "frame": "b0cdaef0",
            "n_artifacts": 34,
            "status": "SOLID exploratory — docs in repo, figure artifact only",
            "key_content": "18 md, 8 csv, 4 png"
        },
        {
            "id": "D",
            "name": "Docking + MD simulations",
            "frame": "bf412257",
            "n_artifacts": 95,
            "status": "ARTIFACTS ONLY — not in repo",
            "key_content": "12 pdf, 10 pdb, 2 sdf, 17 csv, 15 json, 10 png"
        },
        {
            "id": "E",
            "name": "lncRNA analysis",
            "frame": "9961df5f",
            "n_artifacts": 12,
            "status": "ARTIFACTS ONLY — figure in open workspace, not in repo",
            "key_content": "9 csv, 1 R, 1 png"
        },
        {
            "id": "F",
            "name": "Pathway / TF enrichment",
            "frame": "baba2263",
            "n_artifacts": 23,
            "status": "ARTIFACTS ONLY",
            "key_content": "13 csv, 8 png"
        },
        {
            "id": "G",
            "name": "STAT3 isoform analysis",
            "frame": "073dd8db",
            "n_artifacts": 5,
            "status": "RETIRED — correctly dropped; note only",
            "key_content": "2 md"
        },
        {
            "id": "H",
            "name": "Presentation / figure compile",
            "frame": "b16cdb24",
            "n_artifacts": 12,
            "status": "ARTIFACTS ONLY",
            "key_content": "8 png, 1 pdf, 1 pptx"
        },
        {
            "id": "I",
            "name": "Audit / peer review",
            "frame": "afa8da2d",
            "n_artifacts": 22,
            "status": "in repo as PROJECT_AUDIT.md",
            "key_content": "13 csv, 1 md"
        }
    ],
    "phases": [
        {
            "phase": 1,
            "label": "Repo structure + manifest",
            "description": "Add arm subdirs (arms/A_bulk, arms/B_scissor, arms/C_theory1, arms/D_docking, arms/E_lncRNA, arms/F_pathways, arms/H_presentation) under an arms/ top-level directory; write MANIFEST.md mapping every arm to its frame_id, artifact count, and status; commit additive.",
            "effort": "low",
            "reversible": True
        },
        {
            "phase": 2,
            "label": "lncRNA arm (E) — pull generating code + data",
            "description": "Fetch the lncRNA R script (artifact in frame 9961df5f), the two input CSVs (lncRNA_significant_ranked.csv, meta_de_PPvsNN_4study_biotyped.csv), and the figure. Stage under arms/E_lncRNA/. Also fix the open legend annotation on the figure. Commit.",
            "effort": "low",
            "reversible": True
        },
        {
            "phase": 3,
            "label": "Pathway / TF arm (F)",
            "description": "Retrieve artifact CSVs and figures from frame baba2263; add README explaining the enrichment method; stage under arms/F_pathways/.",
            "effort": "medium",
            "reversible": True
        },
        {
            "phase": 4,
            "label": "Docking + MD arm (D)",
            "description": "Retrieve PDB, SDF, CSV, PDF artifacts from frame bf412257; add a methods note; stage under arms/D_docking/. Note: raw MD trajectories are too large for git — reference by cluster path only.",
            "effort": "high",
            "reversible": True
        },
        {
            "phase": 5,
            "label": "Presentation / figure compile (H)",
            "description": "Retrieve final figures and .pptx from frame b16cdb24; stage under arms/H_presentation/.",
            "effort": "low",
            "reversible": True
        },
        {
            "phase": 6,
            "label": "Top-level README + navigation",
            "description": "Update README.md to give a one-paragraph-per-arm overview with paths to each arm dir. Update CHECKPOINT.md to reflect the consolidated layout.",
            "effort": "low",
            "reversible": True
        }
    ],
    "immediate_next": "Phase 1 + Phase 2 (lncRNA figure fix is already in progress — the open workspace annotation is on that figure)",
    "scope_approved_by_user": True,
    "note": "additive commits only throughout; no file moves without explicit approval; legend fix on fig_lncRNA_volcano_summary.png to happen in same session"
}

with open("plan_consolidate-all-psoriasis-arms-into-one_b0cdaef0.json", "w") as f:
    json.dump(plan, f, indent=2)

print(json.dumps(plan, indent=2)[:500])