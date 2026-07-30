# Auto-extracted generating script
# Produces: plan_multi-study-psoriasis-meta-analysis-per_db606fdf.json
# Conda env: unknown   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: a0408b80-29fe-4857-aa8d-429d8aec557a
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import json

plan = {
    "recommendation": "Meta-analysis: per-study DE, then statistical pooling (Recommended)",
    "frame_id": "db606fdf-5628-47d4-b5c3-aac2d0ebf657"
}

with open("plan_multi-study-psoriasis-meta-analysis-per_db606fdf.json", "w") as f:
    json.dump(plan, f, indent=2)