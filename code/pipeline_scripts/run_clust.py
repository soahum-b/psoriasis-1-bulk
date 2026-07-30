import numpy as np
if not hasattr(np, "in1d"):
    np.in1d = np.isin
# also common removed aliases clust may touch
for old,new in [("bool8","bool_"),("object0","object_"),("int0","intp"),("uint0","uintp"),("str0","str_")]:
    if not hasattr(np, old) and hasattr(np, new):
        setattr(np, old, getattr(np, new))
import sys
from clust.__main__ import main
sys.argv = ["clust", "clust_input", "-o", "clust_out"]
main()
