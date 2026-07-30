# Auto-extracted generating script
# Produces: 05_pockets.pdf
# Conda env: dock-md   (run with this environment activated)
# Inputs (expected alongside / in data/): pocket_druggability_summary.csv, pocket_report.csv, style.css, master_refs.json, 05_pockets.md
# Source artifact version: 4923ba94-4cbc-406e-a1ac-19e546c2676f
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import json, re, os
import markdown as md
import weasyprint

master = json.load(open("master_refs.json"))
master["uniprot_2023"]["first_author"] = "The UniProt Consortium"

fig_paths = {}  # no figure for 05_pockets paper

def fmt_ref3(ref):
    auth = ref.get("first_author", "") or ""
    if auth and auth != "The UniProt Consortium":
        auth = f"{auth} et al"
    title = (ref.get("title", "") or "").rstrip(".")
    journal = ref.get("journal", "") or ""
    year = ref.get("year", "")
    doi = ref.get("DOI", "")
    s = f"{auth}. {title}." if auth else f"{title}."
    if journal: s += f" <i>{journal}</i>."
    if year: s += f" {year}."
    return s, doi

def render_paper3(mdpath, name):
    text = open(mdpath).read()
    cite_pat = re.compile(r'\[@([^\]]+)\]')
    order = []
    def collect(m):
        for k in [x.strip().lstrip('@').strip() for x in m.group(1).split(';')]:
            if k not in order: order.append(k)
        return m.group(0)
    cite_pat.sub(collect, text)
    numd = {k: i+1 for i, k in enumerate(order)}
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

CSS = open("style.css").read()
CSS += "\n.figure { text-align:center; margin:14px 0; page-break-inside:avoid; }\n.figure img { max-width:100%; max-height:340px; border:0.6px solid #d0d8e4; }\n.figcap { font-size:8.5px; color:#666; font-style:italic; margin-top:3px; }\n"

name = "05_pockets"
mdpath = "05_pockets.md"
name_out, html_body = render_paper3(mdpath, name)
full = f"<html><head><meta charset='utf-8'><style>{CSS}</style></head><body>{html_body}</body></html>"
weasyprint.HTML(string=full, base_url=".").write_pdf("05_pockets.pdf")