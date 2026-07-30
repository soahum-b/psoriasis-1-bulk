# Auto-extracted generating script
# Produces: STAT3_6NJS.pdb
# Conda env: dock-md   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: a223d16c-fa82-48d8-93f3-2b6a32092bf3
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import shutil
import os

os.makedirs("work/structures", exist_ok=True)

shutil.copy("/tmp/6NJS.pdb", "work/structures/STAT3_6NJS.pdb")