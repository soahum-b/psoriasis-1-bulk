# Auto-extracted generating script
# Produces: lead_network_centrality.csv
# Conda env: python   (run with this environment activated)
# Inputs (expected alongside / in data/): communities.csv, centrality.csv, string_network_full.json
# Source artifact version: 16e2bdf9-882e-4377-a52f-266ac854b552
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import json
import pandas as pd
import numpy as np
import networkx as nx

net = json.load(open("string_network_full.json"))
nodes = pd.DataFrame(net["nodes"])
edges = pd.DataFrame(net["edges"])

cent = pd.read_csv("centrality.csv")
comm = pd.read_csv("communities.csv")

g2c = dict(zip(comm.gene, comm.community))

comm_names = {
    0: "Immune / JAK-STAT / chemokine",
    1: "ECM / matrix remodeling",
    2: "Neuronal / Ca²⁺ signaling",
    3: "Cell cycle / proliferation",
    4: "Lipid / metabolism",
    5: "Interferon / antiviral"
}

G = nx.Graph()
for _, e in edges.iterrows():
    G.add_edge(e["a"], e["b"], weight=e["score"])

nd = nodes.copy()

leads = ["STAT3", "RORC", "JAK3"]

rows = []
for l in leads:
    r = cent[cent.gene == l]
    if len(r):
        b = r["betweenness"].values[0]
        d = int(nd[nd["query"] == l]["degree"].values[0])
        drank = int((nd["degree"] > d).sum()) + 1
        brank = int((cent["betweenness"] > b).sum()) + 1
        rows.append({
            "gene": l,
            "degree": d,
            "degree_rank": drank,
            "degree_pctile": round(100 * (nd["degree"] < d).mean()),
            "betweenness": round(b, 4),
            "betweenness_rank": brank,
            "community": comm_names.get(g2c.get(l), "—")
        })

lead_net = pd.DataFrame(rows)
lead_net.to_csv("lead_network_centrality.csv", index=False)
print(lead_net.to_string(index=False))