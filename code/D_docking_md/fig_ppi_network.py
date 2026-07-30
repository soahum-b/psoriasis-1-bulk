# Auto-extracted generating script
# Produces: fig_ppi_network.png
# Conda env: python   (run with this environment activated)
# Inputs (expected alongside / in data/): communities.csv, centrality.csv, string_network_full.json, network_genes.csv
# Source artifact version: e3115920-bb27-438f-b407-bc00c53480f5
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import json
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.lines import Line2D
from adjustText import adjust_text

# skill:figure-style kernel.py (auto-injected on skill load)
META_GREY = "#888888"


def apply_figure_style(*, frame="open", font=None, sizes=(8, 7, 6), grid=False):
    import matplotlib as mpl
    if frame not in ("open", "boxed", "none"):
        raise ValueError(f"frame must be 'open'|'boxed'|'none', got {frame!r}")
    try:
        import os, sys, glob, matplotlib.font_manager as fm
        fdir = os.path.join(os.environ.get("CONDA_PREFIX") or sys.prefix, "fonts")
        if os.path.isdir(fdir):
            known = {f.fname for f in fm.fontManager.ttflist}
            for f in glob.glob(os.path.join(fdir, "*.ttf")):
                if f not in known:
                    fm.fontManager.addfont(f)
    except Exception:
        pass
    base, secondary, tick = sizes
    boxed = (frame == "boxed")
    rc = {
        "font.family": "sans-serif",
        "font.size": base,
        "axes.labelsize": base,
        "axes.titlesize": base,
        "legend.fontsize": secondary,
        "xtick.labelsize": tick,
        "ytick.labelsize": tick,
        "axes.linewidth": 0.6,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "axes.spines.top": boxed, "axes.spines.right": boxed,
        "axes.spines.left": frame != "none", "axes.spines.bottom": frame != "none",
        "axes.grid": bool(grid),
        "legend.frameon": False,
        "figure.dpi": 200,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.titleweight": "normal",
        "axes.titlelocation": "left",
        "axes.labelweight": "normal",
        "lines.linewidth": 1.2,
        "patch.linewidth": 0.6,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    }
    if font:
        rc["font.sans-serif"] = [font, "DejaVu Sans"]
    mpl.rcParams.update(rc)


apply_figure_style(sizes=(8,7,6))

# Load data
net = json.load(open("string_network_full.json"))
nodes = pd.DataFrame(net["nodes"])
edges = pd.DataFrame(net["edges"])
meta = pd.read_csv("network_genes.csv")
comm = pd.read_csv("communities.csv")
cent = pd.read_csv("centrality.csv")

nd = nodes.merge(meta, left_on="query", right_on="gene", how="left")
nd_sorted = nd.sort_values("degree", ascending=False)

g2c = dict(zip(comm.gene, comm.community))

comm_names = {0:"Immune / JAK-STAT / chemokine", 1:"ECM / matrix remodeling",
 2:"Neuronal / Ca²⁺ signaling", 3:"Cell cycle / proliferation",
 4:"Lipid / metabolism", 5:"Interferon / antiviral"}
big = [0,1,2,3,4,5]
palette = {0:"#c0392b", 1:"#2c7fb8", 2:"#8e6fb0", 3:"#e67e22", 4:"#16a085", 5:"#d4a017"}

def ncol(g):
    c = g2c.get(g)
    return palette.get(c, "#cfd4da")

G = nx.Graph()
for _, e in edges.iterrows():
    G.add_edge(e["a"], e["b"], weight=e["score"])
cc = max(nx.connected_components(G), key=len)
Gc = G.subgraph(cc)
print(f"graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges; largest CC: {len(cc)}")

bc = nx.betweenness_centrality(Gc, seed=42)
deg = dict(G.degree())

posG = nx.spring_layout(Gc, k=0.35, iterations=60, seed=42)
print("layout done for", len(posG), "nodes")

leads = ["STAT3","RORC","JAK3"]
top_hubs = nd_sorted["query"].head(6).tolist()
label_set = set(top_hubs) | set(leads)

fig = plt.figure(figsize=(13,6.2))
gs = fig.add_gridspec(1,2, width_ratios=[1.55,1.0], wspace=0.16)
axN = fig.add_subplot(gs[0]); axH = fig.add_subplot(gs[1])

# --- Panel A: network ---
for _, e in edges.iterrows():
    a, b = e["a"], e["b"]
    if a in posG and b in posG:
        x = [posG[a][0], posG[b][0]]; y = [posG[a][1], posG[b][1]]
        axN.plot(x, y, color="#d9dde3", lw=0.18, alpha=0.5, zorder=1)

degs = np.array([deg.get(g,1) for g in posG])
xs = [posG[g][0] for g in posG]; ys = [posG[g][1] for g in posG]
cols = [ncol(g) for g in posG]
sizes = [8+deg.get(g,1)*1.6 for g in posG]
axN.scatter(xs, ys, s=sizes, c=cols, edgecolors="white", linewidths=0.2, zorder=2, alpha=0.9)

for l in leads:
    if l in posG:
        x, y = posG[l]
        axN.scatter([x], [y], s=[8+deg.get(l,1)*1.6+40], facecolors="none",
                    edgecolors="black", linewidths=1.6, zorder=4)

txts = []
for g in label_set:
    if g in posG:
        x, y = posG[g]
        bold = g in leads
        txts.append(axN.text(x, y, g, fontsize=7 if bold else 6.2,
                    fontweight="bold" if bold else "normal",
                    color="black" if bold else "#333", zorder=5))
adjust_text(txts, ax=axN, arrowprops=dict(arrowstyle="-", color="0.5", lw=0.5),
            expand=(1.3,1.5))
axN.set_xticks([]); axN.set_yticks([])
for s in axN.spines.values(): s.set_visible(False)
axN.set_title("Psoriasis DEG protein–protein interaction network\n(STRING high-confidence ≥0.7; leads ringed)", fontsize=8.5)

leg = [Line2D([0],[0], marker='o', color='w', markerfacecolor=palette[c], markersize=7,
     label=f"{comm_names[c]}") for c in big]
axN.legend(handles=leg, loc="lower left", frameon=False, fontsize=5.8, handletextpad=0.2, labelspacing=0.25)

# --- Panel B: hub centrality ---
top15 = cent.sort_values("betweenness", ascending=False).head(15).copy()
top15 = top15[::-1]
ypos = np.arange(len(top15))
barcols = ["#c0392b" if g in leads else "#9aa6b5" for g in top15.gene]
axH.barh(ypos, top15.betweenness, color=barcols, edgecolor="white", height=0.72)
axH.set_yticks(ypos); axH.set_yticklabels(top15.gene, fontsize=6.5)
for i, g in enumerate(top15.gene):
    if g in leads:
        axH.get_yticklabels()[i].set_fontweight("bold")
axH.set_xlabel("Betweenness centrality (bottleneck score)", fontsize=7.5)
axH.set_title("Network-central hubs\n(bridges between functional modules)", fontsize=8.5)
for s in ["top","right"]: axH.spines[s].set_visible(False)

s3 = cent.sort_values("betweenness", ascending=False).reset_index(drop=True)
s3rank = s3[s3.gene=="STAT3"].index[0]+1
axH.annotate(f"STAT3 = #{s3rank} of {len(cent)}\n(designated lead)",
             xy=(top15[top15.gene=='STAT3'].betweenness.values[0], list(top15.gene).index("STAT3")),
             xytext=(0.45,0.5), textcoords="axes fraction", fontsize=6.5,
             arrowprops=dict(arrowstyle="->", color="#c0392b", lw=1))

fig.savefig("fig_ppi_network.png", dpi=200, bbox_inches="tight")
print("saved; STAT3 betweenness rank:", s3rank)