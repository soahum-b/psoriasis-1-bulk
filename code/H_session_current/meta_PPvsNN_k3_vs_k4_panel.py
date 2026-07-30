# Auto-extracted generating script
# Produces: meta_PPvsNN_k3_vs_k4_panel.csv
# Conda env: python   (run with this environment activated)
# Inputs (expected alongside / in data/): meta_de_PPvsNN.csv
# Source artifact version: 34f3b909-e18e-47cf-9746-25b87c801ebc
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import pandas as pd
import numpy as np

# Load k=3 meta-analysis results
m3 = pd.read_csv("meta_de_PPvsNN.csv").set_index("gene")

# Load k=4 meta-analysis results (created during the trace)
m4 = pd.read_csv("meta_de_PPvsNN_4study.csv").set_index("gene")

panel_all = {
 "STAT3 / JAK-STAT": ["STAT3","STAT1","JAK3","TYK2","JAK1","SOCS3"],
 "IL-1$\\beta$ / inflammasome": ["IL1B","IL1RN","NLRP3","CASP1","PYCARD","AIM2","CASP5","GSDMD","IL1A","IL18"],
 "IL-36 family": ["IL36G","IL36A","IL36RN"],
 "Th17 axis": ["IL17A","IL23A","RORC"],
 "IL-1 output": ["S100A9","S100A8","S100A7","DEFB4A","CCL20","CXCL2","LCN2"],
}

rows=[]
for mod,gs in panel_all.items():
    for g in gs:
        if g in m4.index and g in m3.index:
            rows.append((mod,g,
                m3.loc[g,"logFC"],m3.loc[g,"FDR"],int(m3.loc[g,"k"]),m3.loc[g,"I2"],
                m4.loc[g,"logFC"],m4.loc[g,"FDR"],int(m4.loc[g,"k"]),m4.loc[g,"I2"],m4.loc[g,"SE"]))
C=pd.DataFrame(rows,columns=["module","gene","lfc3","fdr3","k3","i23","lfc4","fdr4","k4","i24","se4"])
C["dlfc"]=C.lfc4-C.lfc3
def sig(p): return "***" if p<1e-3 else "**" if p<1e-2 else "*" if p<0.05 else "ns"
C["sig3"]=C.fdr3.apply(sig); C["sig4"]=C.fdr4.apply(sig)
C["flip"]=np.where((C.fdr3<0.05)!=(C.fdr4<0.05),"FLIP","")
C.to_csv("meta_PPvsNN_k3_vs_k4_panel.csv", index=False)