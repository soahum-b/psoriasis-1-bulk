# Auto-extracted generating script
# Produces: GenomeInfoDbData_1.2.13.tar.gz
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 20500bb4-c941-4eb2-8b04-eca4f1e0b81c
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import urllib.request

url = "https://bioconductor.org/packages/3.20/data/annotation/src/contrib/GenomeInfoDbData_1.2.13.tar.gz"
urllib.request.urlretrieve(url, "GenomeInfoDbData_1.2.13.tar.gz")