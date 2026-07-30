# Auto-extracted generating script
# Produces: psoriasis_pipeline.R
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: b716656d-b7f5-4446-9ec1-bbf0e70f325a
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

# [lineage] reconstruction failed validation (max_tokens) — raw producing cell below; see Execution Log for full trace
SRC=/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3
WS=/Users/soahum/.claude-science/orgs/b2a12749-be9b-47f9-8aaf-51f78dadc583/workspaces/db606fdf-5628-47d4-b5c3-aac2d0ebf657
cp "$SRC"/psoriasis_pipeline.R "$WS"/
mkdir -p "$WS"/staging_axis/code
cp "$SRC"/staging_axis/code/staging_figures.R "$WS"/staging_axis/code/
echo copied; ls "$WS"/staging_axis/code/