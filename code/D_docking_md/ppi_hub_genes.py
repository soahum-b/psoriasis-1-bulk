# Auto-extracted generating script
# Produces: ppi_hub_genes.csv
# Conda env: python   (run with this environment activated)
# Inputs (expected alongside / in data/): communities.csv, centrality.csv, string_network_full.json
# Source artifact version: 31c360bf-6cfd-4513-866b-9d9404c080e5
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import json
import pandas as pd
import numpy as np

net = json.load(open("string_network_full.json"))
nodes = pd.DataFrame(net["nodes"])
edges = pd.DataFrame(net["edges"])

meta = pd.read_csv("centrality.csv")

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

network_genes = meta.copy()

nd = nodes.merge(network_genes, left_on="query", right_on="gene", how="left")
nd_sorted = nd.sort_values("degree", ascending=False)

import networkx as nx
G = nx.Graph()
for _, e in edges.iterrows():
    G.add_edge(e["a"], e["b"], weight=e["score"])
cc = max(nx.connected_components(G), key=len)
Gc = G.subgraph(cc)

bc = nx.betweenness_centrality(Gc, seed=42)
deg = dict(G.degree())
cent = pd.DataFrame({"gene": list(bc.keys()), "betweenness": list(bc.values())})
cent["degree"] = cent["gene"].map(deg)
cent = cent.merge(network_genes[["gene", "direction", "logFC", "protein_class", "ot_sm_bucket", "evidence_rank"]], on="gene", how="left")

hub = cent.merge(nd[["query", "degree"]].rename(columns={"query": "gene", "degree": "degree_check"}), on="gene", how="left")
hub["community"] = hub["gene"].map(g2c)
hub["community_name"] = hub["community"].map(comm_names)
hub = hub.sort_values("degree", ascending=False)
hub_out = hub[["gene", "degree", "betweenness", "direction", "logFC", "protein_class", "ot_sm_bucket", "community_name", "evidence_rank"]]
hub_out.to_csv("ppi_hub_genes.csv", index=False)