# Auto-extracted generating script
# Produces: plan_theory-2-sex-stratified-psoriasis-skin-a_5a65de78.json
# Conda env: unknown   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 014ef68d-493c-4617-a8cf-daabd6e2bf42
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

# [lineage] reconstruction failed validation (prose) — raw producing cell below; see Execution Log for full trace
r = host.artifacts(limit=200)
print(r['count'])
for a in r['artifacts']:
    print(a['filename'], '|', a['content_type'], '|', a['size_bytes'], '|', a['latest_version_id'])
