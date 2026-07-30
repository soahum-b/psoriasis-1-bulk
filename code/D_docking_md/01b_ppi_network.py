# Auto-extracted generating script
# Produces: 01b_ppi_network.pdf
# Conda env: dock-md   (run with this environment activated)
# Inputs (expected alongside / in data/): string_network_full.json, centrality.csv, lead_network_centrality.csv, ppi_hub_genes.csv, 01b_ppi_network.md, master_refs.json, fig_paths.json
# Source artifact version: 34422a8f-514d-4ba4-b946-117ec372768e
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import json
import os
import re
import markdown as md
import weasyprint

master = json.load(open("master_refs.json"))
fig_paths = json.load(open("fig_paths.json"))

CSS = """
@page { size: A4; margin: 2cm 2.2cm;
  @bottom-center { content: counter(page); font-size: 9px; color: #888; }
  @bottom-right { content: "Psoriasis druggability pipeline"; font-size: 8px; color: #aaa; } }
body { font-family: 'DejaVu Serif', Georgia, serif; font-size: 10.5px; line-height: 1.5; color: #1a1a1a; }
h1 { font-family: 'DejaVu Sans', Helvetica, sans-serif; font-size: 19px; color: #1a2e4a; border-bottom: 2.5px solid #2c5aa0; padding-bottom: 8px; margin-bottom: 4px; }
h2 { font-family: 'DejaVu Sans', Helvetica, sans-serif; font-size: 13.5px; color: #2c5aa0; margin-top: 18px; border-bottom: 0.7px solid #c8d4e6; padding-bottom: 3px; }
h3 { font-family: 'DejaVu Sans', Helvetica, sans-serif; font-size: 11.5px; color: #37527a; margin-top: 12px; }
p { margin: 6px 0; text-align: justify; }
blockquote { background:#eef3fa; border-left:4px solid #2c5aa0; margin:10px 0; padding:8px 12px; font-style:italic; color:#2a3a52; }
code { font-family: 'DejaVu Sans Mono', monospace; font-size: 9px; background: #f2f4f8; padding: 1px 3px; border-radius: 2px; color: #b0451a; }
pre { background: #f2f4f8; border-left: 3px solid #2c5aa0; padding: 8px 10px; font-size: 9px; border-radius: 3px; }
pre code { background: none; padding: 0; color: #222; }
table { border-collapse: collapse; width: 100%; margin: 10px 0; font-size: 9.3px; }
th { background: #2c5aa0; color: white; padding: 5px 7px; text-align: left; font-family: 'DejaVu Sans', sans-serif; }
td { padding: 4px 7px; border-bottom: 0.6px solid #d8e0ec; }
tr:nth-child(even) { background: #f6f8fb; }
sup { color: #2c5aa0; font-weight: bold; font-size: 8px; }
ol.refs { font-size: 9px; line-height: 1.45; padding-left: 18px; }
ol.refs li { margin-bottom: 5px; }
ol.refs a { color: #2c5aa0; text-decoration: none; word-break: break-all; }
strong { color: #12233d; }
.figure { text-align:center; margin:14px 0; page-break-inside:avoid; }
.figure img { max-width:100%; max-height:360px; border:0.6px solid #d0d8e4; }
.figcap { font-size:8.5px; color:#666; font-style:italic; margin-top:3px; }
"""

def fmt_ref(ref):
    if "_cite" in ref:
        return ref["_cite"].rstrip('. ') + ".", ref["DOI"]
    auth = ref.get("first_author", "") or ""
    if auth and auth != "The UniProt Consortium":
        auth = f"{auth} et al"
    title = (ref.get("title", "") or "").rstrip(".")
    journal = ref.get("journal", "") or ""
    year = ref.get("year", "")
    s = f"{auth}. {title}." if auth else f"{title}."
    if journal:
        s += f" <i>{journal}</i>."
    if year:
        s += f" {year}."
    return s, ref.get("DOI", "")

def render_paper(mdpath):
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
    numd = {k: i+1 for i, k in enumerate(order)}
    missing = [k for k in order if k not in master]
    def repl(m):
        nums = [str(numd[k]) for k in [x.strip().lstrip('@').strip() for x in m.group(1).split(';')] if k in numd]
        return "<sup>[" + ",".join(nums) + "]</sup>"
    body = cite_pat.sub(repl, text)
    reflines = []
    for k in order:
        if k in master:
            s, doi = fmt_ref(master[k])
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
    return name, html_body, missing

mp = "01b_ppi_network.md"
name, html, miss = render_paper(mp)
full = f"<html><head><meta charset='utf-8'><style>{CSS}</style></head><body>{html}</body></html>"
weasyprint.HTML(string=full, base_url=".").write_pdf("01b_ppi_network.pdf")
print(f"{name}: rendered, missing={miss}")