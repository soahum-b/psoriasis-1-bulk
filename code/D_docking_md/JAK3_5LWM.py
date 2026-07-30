# Auto-extracted generating script
# Produces: JAK3_5LWM.pdb
# Conda env: dock-md   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: c8912a6e-8970-48a8-97da-9d7c85b5d187
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import shutil
import os

os.makedirs("work/structures", exist_ok=True)

shutil.copy("/tmp/5LWM.pdb", "work/structures/JAK3_5LWM.pdb")