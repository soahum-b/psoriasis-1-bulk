# Auto-extracted generating script
# Produces: psoriasis_docking_package.tar.gz
# Conda env: dock-md   (run with this environment activated)
# Inputs (expected alongside / in data/): validation.json, docking_ranked.csv, STAT3_6NJS.pdb, RORC_5APH.pdb, JAK3_5LWM.pdb, RORC_apo.pdb, RORC_receptor.pdbqt, screening_library.csv, STAT3_apo.pdb, JAK3_apo.pdb, domain_sequences.json, vyi_refH.sdf
# Source artifact version: 88a2478c-1282-4fc1-9668-c29c2f5b04b9
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

# [lineage] reconstruction failed validation (max_tokens) — raw producing cell below; see Execution Log for full trace
cd work/cluster_packages && tar -czf psoriasis_docking_package.tar.gz 01_docking/ && echo "re-archived: $(ls -la psoriasis_docking_package.tar.gz | awk '{print $5}') bytes"