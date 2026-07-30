# Auto-extracted generating script
# Produces: 06_library.pdf
# Conda env: dock-md   (run with this environment activated)
# Inputs (expected alongside / in data/): screening_library.csv, style.css, master_refs.json, 06_library.md
# Source artifact version: ff2d12bb-8cf9-4acf-b3ce-9f0f44296caf
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import json, re, os
import markdown as md
import weasyprint
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

# Load master refs
master = json.load(open("master_refs.json"))

# Load screening library for the figure
libf = pd.read_csv("screening_library.csv")

# Generate fig_library_properties.png
mpl.rcParams.update({"font.size": 8, "axes.linewidth": 0.6, "font.family": "DejaVu Sans"})
import numpy as np

fig = plt.figure(figsize=(11, 6.2))
gs = fig.add_gridspec(2, 4, hspace=0.42, wspace=0.34, height_ratios=[1, 1])
props = [("mw", "Molecular weight (Da)", (250, 600)), ("logp", "cLogP", None),
         ("qed", "QED drug-likeness", (0, 1)), ("tpsa", "TPSA (Å²)", None),
         ("hbd", "H-bond donors", None), ("hba", "H-bond acceptors", None),
         ("rot_bonds", "Rotatable bonds", None), ("arom_rings", "Aromatic rings", None)]
for i, (col, lab, xr) in enumerate(props):
    ax = fig.add_subplot(gs[i // 4, i % 4])
    ax.hist(libf[col], bins=18, color="#4575b4", edgecolor="white", lw=0.4)
    ax.set_xlabel(lab, fontsize=7)
    ax.set_ylabel("count", fontsize=7)
    if xr:
        ax.set_xlim(*xr)
    ax.axvline(libf[col].median(), color="#d73027", lw=1, ls="--")
    ax.tick_params(labelsize=6.5)
fig.suptitle(
    f"RORγt-focused screening library — {len(libf)} drug-like compounds (ChEMBL bioactivity space; purchasability pending ZINC22)",
    fontsize=8.6, y=0.98)
fig_path = "fig_library_properties.png"
fig.savefig(fig_path, dpi=190, bbox_inches="tight")
plt.close()

fig_abs_path = os.path.abspath(fig_path)

# Load CSS
CSS = open("style.css").read()
CSS += "\n.figure { text-align:center; margin:14px 0; page-break-inside:avoid; }\n.figure img { max-width:100%; max-height:340px; border:0.6px solid #d0d8e4; }\n.figcap { font-size:8.5px; color:#666; font-style:italic; margin-top:3px; }\n"


def fmt_ref3(ref):
    auth = ref.get("first_author", "") or ""
    if auth and auth != "The UniProt Consortium":
        auth = f"{auth} et al"
    title = (ref.get("title", "") or "").rstrip(".")
    journal = ref.get("journal", "") or ""
    year = ref.get("year", "")
    doi = ref.get("DOI", "")
    s = f"{auth}. {title}." if auth else f"{title}."
    if journal:
        s += f" <i>{journal}</i>."
    if year:
        s += f" {year}."
    return s, doi


def render_paper3(mdpath, fig_paths):
    name = os.path.basename(mdpath).replace(".md", "")
    text = open(mdpath).read()
    cite_pat = re.compile(r'\[@([^\]]+)\]')
    order = []

    def collect(m):
        for k in [x.strip().lstrip('@').strip() for x in m.group(1).split(';')]:
            if k not in order:
                order.append(k)
        return m.group(0)

    cite_pat.sub(collect, text)
    numd = {k: i + 1 for i, k in enumerate(order)}

    def repl(m):
        nums = [str(numd[k]) for k in [x.strip().lstrip('@').strip() for x in m.group(1).split(';')] if k in numd]
        return "<sup>[" + ",".join(nums) + "]</sup>"

    body = cite_pat.sub(repl, text)
    reflines = []
    for k in order:
        if k in master:
            s, doi = fmt_ref3(master[k])
            reflines.append(f'<li>{s} doi:<a href="https://doi.org/{doi}">{doi}</a></li>')
    refs_html = "<ol class='refs'>\n" + "\n".join(reflines) + "\n</ol>"
    body = re.sub(r'##\s*References\s*$', '', body).rstrip()
    html_body = md.markdown(body, extensions=['tables', 'fenced_code', 'sane_lists'])
    if name in fig_paths:
        img = f'<div class="figure"><img src="file://{fig_paths[name]}"/><p class="figcap">Figure: key output of this pipeline step.</p></div>'
        if ">Deliverables<" in html_body:
            html_body = re.sub(r'(<h2>[^<]*Deliverables)', img + r'\1', html_body, count=1)
        else:
            html_body += img
    html_body += "\n<h2>References</h2>\n" + refs_html
    return name, html_body


fig_paths = {"06_library": fig_abs_path}

mdpath = "06_library.md"
name, html_body = render_paper3(mdpath, fig_paths)
full = f"<html><head><meta charset='utf-8'><style>{CSS}</style></head><body>{html_body}</body></html>"
weasyprint.HTML(string=full, base_url=".").write_pdf("06_library.pdf")
print("06_library.pdf written")