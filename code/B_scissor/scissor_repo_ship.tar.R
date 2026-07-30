# Auto-extracted generating script
# Produces: scissor_repo_ship.tar.gz
# Conda env: scissor-r   (run with this environment activated)
# Inputs (expected alongside / in data/): deconv_scissor_concordance.csv, celltype_enrichment.csv, scissor_run.R, scissor_reliability.R, scissor_glmnet_solver.R, sample_manifest.csv, Recount-3/sample_classification.rds
# Source artifact version: cd755d82-4107-4890-8887-947d5b00a4f2
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

# [lineage] reconstruction failed validation (max_tokens) — raw producing cell below; see Execution Log for full trace
WS="/Users/soahum/.claude-science/orgs/b2a12749-be9b-47f9-8aaf-51f78dadc583/workspaces/689f69d8-f140-45ab-94f7-01ef3163cf0c/scissor_repo" && \
cp "$WS/code/run_pipeline.R" /tmp/scissor_ship/code/ && \
cp "$WS/.gitignore" /tmp/scissor_ship/ && \
cd /tmp && rm -f scissor_repo_ship.tar.gz && tar czf scissor_repo_ship.tar.gz -C /tmp/scissor_ship . && \
cp scissor_repo_ship.tar.gz "$WS/" && ls -lh "$WS/scissor_repo_ship.tar.gz" && \
echo "--- ship contents ---" && tar tzf scissor_repo_ship.tar.gz | sort | grep -vE '/$'