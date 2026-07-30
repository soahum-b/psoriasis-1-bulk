# Auto-extracted generating script
# Produces: RORC_5APH.pdb
# Conda env: dock-md   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: c94075ca-4f37-4588-ad98-f96d4d5b88ab
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import shutil
import os

os.makedirs("work/structures", exist_ok=True)

shutil.copy("/tmp/5APH.pdb", "work/structures/RORC_5APH.pdb")