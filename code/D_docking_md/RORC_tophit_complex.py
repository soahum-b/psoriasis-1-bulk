# Auto-extracted generating script
# Produces: RORC_tophit_complex.pdb
# Conda env: dock-md   (run with this environment activated)
# Inputs (expected alongside / in data/): RORC_apo.pdb
# Source artifact version: 2beb0615-a78b-47ec-9dd9-ca10e5778060
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import os

# Build RORC_tophit_complex.pdb from the receptor apo structure
# (ligand pose requires docking results not available as artifact)
with open('RORC_tophit_complex.pdb', 'w') as out:
    with open("RORC_apo.pdb") as f:
        for line in f:
            if line.startswith('ATOM'):
                out.write(line)
    out.write('TER\n')
    out.write('END\n')