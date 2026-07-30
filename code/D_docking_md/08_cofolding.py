# Auto-extracted generating script
# Produces: 08_cofolding.pdf
# Conda env: dock-md   (run with this environment activated)
# Inputs (expected alongside / in data/): style.css, master_refs.json, 08_cofolding.md, RORC_VYI_crystal_posctrl.yaml
# Source artifact version: ad2c68e2-2c48-4011-be39-a8946ee3159e
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import json, re, os
import markdown as md
import weasyprint

master = json.load(open("master_refs.json"))
master["uniprot_2023"]["first_author"] = "The UniProt Consortium"

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

CSS = open("style.css").read()
CSS += "\n.figure { text-align:center; margin:14px 0; page-break-inside:avoid; }\n.figure img { max-width:100%; max-height:340px; border:0.6px solid #d0d8e4; }\n.figcap { font-size:8.5px; color:#666; font-style:italic; margin-top:3px; }\n"

fig_paths = {}  # no figure for paper 08_cofolding

mdpath = "08_cofolding.md"
name = "08_cofolding"
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
html_body += "\n<h2>References</h2>\n" + refs_html

full = f"<html><head><meta charset='utf-8'><style>{CSS}</style></head><body>{html_body}</body></html>"
weasyprint.HTML(string=full, base_url=".").write_pdf("08_cofolding.pdf")