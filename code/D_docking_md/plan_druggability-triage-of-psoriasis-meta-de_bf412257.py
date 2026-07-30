# Auto-extracted generating script
# Produces: plan_druggability-triage-of-psoriasis-meta-de_bf412257.json
# Conda env: unknown   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 4f846efb-fd26-4613-980a-648ee262a00c
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

# [lineage] reconstruction failed validation (prose) — raw producing cell below; see Execution Log for full trace
cd "/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3" && echo "=== STAT3 (unquoted) ===" && awk -F',' '$1=="STAT3"' meta_de_PPvsNN.csv && echo && echo "=== how many pass FDR<0.05 & |logFC|>1 ===" && awk -F',' 'NR>1 && $7<0.05 && ($3>1 || $3< -1)' meta_de_PPvsNN.csv | wc -l && echo "=== how many pass FDR<0.05 & |logFC|>1 & k>=2 (in >=2 studies) ===" && awk -F',' 'NR>1 && $7<0.05 && ($3>1 || $3< -1) && $2>=2' meta_de_PPvsNN.csv | wc -l