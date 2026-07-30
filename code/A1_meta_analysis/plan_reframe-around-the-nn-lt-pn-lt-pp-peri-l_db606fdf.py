# Auto-extracted generating script
# Produces: plan_reframe-around-the-nn-lt-pn-lt-pp-peri-l_db606fdf.json
# Conda env: python   (run with this environment activated)
# Inputs (expected alongside / in data/): plan_reframe-around-the-nn-lt-pn-lt-pp-peri-l_db606fdf.json
# Source artifact version: 67eb2bfd-3644-4ef6-8e24-1d641cd09d58
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import json
import os

p = host.artifact_path("b4d2c633-dd37-42cc-9d06-0d6b59e3d2f5")
plan = json.load(open(p))

steps = plan["phases"][0]["delegations"][0]["steps"]
# rewrite step 7 (index 6)
steps[6]["title"] = "Write a NEW staging-axis white paper (separate document)"
steps[6]["description"] = (
    "Create a NEW, standalone white paper (e.g. psoriasis_staging_axis_whitepaper.md) that leads "
    "with the peri-lesional NN<PN<PP molecular-staging axis as the central result, presents the "
    "STAT3/IL-17 confirmation as validation of the axis, and integrates the new gradient analyses "
    "(trend test, gradient-gene taxonomy, pathway-timing, early/PN-specific genes) plus the "
    "literature-positioning section. The original psoriasis_meta_analysis_whitepaper.md is left "
    "INTACT as the companion STAT3-led reference document - this is a second paper, not an overwrite. "
    "The new paper embeds the new figures with its own numbering and reference list, and cross-references "
    "the original where the shared meta-analysis methods/figures already live. Deliverable: "
    "psoriasis_staging_axis_whitepaper.md (new artifact)."
)
# also update desired_outputs wording
plan["desired_outputs"] = [o.replace("Reframed white-paper narrative leading with the gradient, STAT3/IL-17 as validation",
                                     "NEW standalone staging-axis white paper (original STAT3-led paper kept intact)") for o in plan["desired_outputs"]]

out = "plan_reframe-around-the-nn-lt-pn-lt-pp-peri-l_db606fdf.json"
json.dump(plan, open(out, "w"), indent=2)