# Auto-extracted generating script
# Produces: whitepapers_psoriasis_pipeline.tar.gz
# Conda env: dock-md   (run with this environment activated)
# Inputs (expected alongside / in data/): 06_library.md, 09_md.md, 08_cofolding.md, 02_druggability_annotation.md, 05_pockets.md, 07_docking.md, 04_structures.md, 10_report.md, 03_scoring_leads.md, 10_report.pdf, 08_cofolding.pdf, 05_pockets.pdf, 09_md.pdf, 04_structures.pdf, 03_scoring_leads.pdf, 06_library.pdf, 02_druggability_annotation.pdf, 07_docking.pdf, 01b_ppi_network.md, 01_target_list.md
# Source artifact version: f6a6312b-72a4-4067-b405-9c12adf1f1ae
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import os
import re
import json
import shutil
import tarfile
import glob
import markdown as md
import weasyprint

# ── Load master refs from the BIBLIOGRAPHY.md artifact ──────────────────────
bib_text = open("01_target_list.md").read()
# Parse existing key->DOI from the bibliography text (reconstructed master)
master = {}
for line in bib_text.splitlines():
    m = re.match(r'\-\s*\*\*\[([^\]]+)\]\*\*\s*(.*?)\s*https://doi\.org/(\S+)', line)
    if m:
        key, cite, doi = m.group(1), m.group(2), m.group(3).rstrip('. ')
        master[key] = {"_cite": cite.strip(), "DOI": doi}

# Add network refs (structured, from the network refs verified in the trace)
network_refs = {
    "string_v12": {"DOI": "10.1093/nar/gkac1000", "first_author": "Szklarczyk", "year": 2023,
                   "title": "The STRING database in 2023: protein-protein association networks and functional enrichment analyses for any of 14 000 organisms", "journal": "Nucleic Acids Research"},
    "cytoscape": {"DOI": "10.1101/gr.1239303", "first_author": "Shannon", "year": 2003,
                  "title": "Cytoscape: A software environment for integrated models of biomolecular interaction networks", "journal": "Genome Research"},
    "freeman_betweenness": {"DOI": "10.2307/3033543", "first_author": "Freeman", "year": 1977,
                            "title": "A set of measures of centrality based on betweenness", "journal": "Sociometry"},
    "louvain": {"DOI": "10.1088/1742-5468/2008/10/P10008", "first_author": "Blondel", "year": 2008,
                "title": "Fast unfolding of communities in large networks", "journal": "Journal of Statistical Mechanics: Theory and Experiment"},
    "network_medicine": {"DOI": "10.1038/nrg2918", "first_author": "Barabasi", "year": 2011,
                         "title": "Network medicine: a network-based approach to human disease", "journal": "Nature Reviews Genetics"},
    "network_biology": {"DOI": "10.1038/nrg1272", "first_author": "Barabasi", "year": 2004,
                        "title": "Network biology: understanding the cell's functional organization", "journal": "Nature Reviews Genetics"},
    "lethality_centrality": {"DOI": "10.1038/35075138", "first_author": "Jeong", "year": 2001,
                             "title": "Lethality and centrality in protein networks", "journal": "Nature"},
}
for k, v in network_refs.items():
    master[k] = v

# ── Figure paths ─────────────────────────────────────────────────────────────
fig_for_paper = {
    "01_target_list": "fig_target_volcano_annotated.png",
    "01b_ppi_network": "fig_ppi_network.png",
    "02_druggability_annotation": "fig_druggability_landscape.png",
    "03_scoring_leads": "fig_target_ranking.png",
    "04_structures": "fig_lead_structures.png",
    "05_pockets": "fig_pockets_druggability.png",
    "06_library": "fig_library_properties.png",
    "07_docking": "fig_docking_results.png",
    "10_report": "fig_pipeline_funnel.png",
}
fig_paths = {}
for paper, fn in fig_for_paper.items():
    r = host.artifacts(filename=fn, exact=True)
    if r["artifacts"]:
        fig_paths[paper] = host.artifact_path(r["artifacts"][0]["latest_version_id"])

# ── CSS ──────────────────────────────────────────────────────────────────────
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

os.makedirs("work/whitepapers/md", exist_ok=True)
os.makedirs("work/whitepapers/pdf", exist_ok=True)
open("work/whitepapers/style.css", "w").write(CSS)

# ── Markdown source files ────────────────────────────────────────────────────
md_sources = {
    "01_target_list": "01_target_list.md",
    "01b_ppi_network": "01b_ppi_network.md",
    "02_druggability_annotation": "02_druggability_annotation.md",
    "03_scoring_leads": "03_scoring_leads.md",
    "04_structures": "04_structures.md",
    "05_pockets": "05_pockets.md",
    "06_library": "06_library.md",
    "07_docking": "07_docking.md",
    "08_cofolding": "08_cofolding.md",
    "09_md": "09_md.md",
    "10_report": "10_report.md",
}

# Copy md sources into working directory
for name, src_path in md_sources.items():
    shutil.copy(src_path, f"work/whitepapers/md/{name}.md")

# ── PDF artifacts (pre-rendered individual PDFs) ─────────────────────────────
pdf_artifacts = {
    "02_druggability_annotation": "02_druggability_annotation.pdf",
    "03_scoring_leads": "03_scoring_leads.pdf",
    "04_structures": "04_structures.pdf",
    "05_pockets": "05_pockets.pdf",
    "06_library": "06_library.pdf",
    "07_docking": "07_docking.pdf",
    "08_cofolding": "08_cofolding.pdf",
    "09_md": "09_md.pdf",
    "10_report": "10_report.pdf",
}

# Copy pre-rendered PDFs
for name, src_path in pdf_artifacts.items():
    shutil.copy(src_path, f"work/whitepapers/pdf/{name}.pdf")

# ── Render helper functions ───────────────────────────────────────────────────
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
    numd = {k: i + 1 for i, k in enumerate(order)}

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
    return name, html_body


# ── Render remaining PDFs not available as pre-rendered artifacts ─────────────
remaining_to_render = ["01_target_list", "01b_ppi_network"]
for name in remaining_to_render:
    mp = f"work/whitepapers/md/{name}.md"
    n, html = render_paper(mp)
    full = f"<html><head><meta charset='utf-8'><style>{CSS}</style></head><body>{html}</body></html>"
    weasyprint.HTML(string=full, base_url="work/whitepapers").write_pdf(f"work/whitepapers/pdf/{name}.pdf")

# ── Build combined PDF ────────────────────────────────────────────────────────
titles = [
    ("01_target_list", "Assembling the ranked DEG target list"),
    ("01b_ppi_network", "The protein–protein interaction network: a guided explainer"),
    ("02_druggability_annotation", "Druggability annotation via knowledge bases"),
    ("03_scoring_leads", "Scoring, ranking, and lead selection"),
    ("04_structures", "Retrieving structures for the lead targets"),
    ("05_pockets", "Binding-pocket detection & structural druggability"),
    ("06_library", "Preparing the purchasable screening library"),
    ("07_docking", "Virtual screening by molecular docking"),
    ("08_cofolding", "Structure-based cross-check by co-folding"),
    ("09_md", "Building MD-ready systems for the cluster"),
    ("10_report", "Compiling the report and deliverables"),
]

title_html = f"""<div style="page-break-after:always; text-align:center; padding-top:120px;">
<h1 style="font-size:26px; border:none; color:#1a2e4a;">Psoriasis Druggability \u2192 Docking \u2192 MD Pipeline</h1>
<h2 style="font-size:16px; border:none; color:#2c5aa0; margin-top:8px;">Methods White Papers</h2>
<p style="font-size:12px; color:#555; margin-top:30px;">Per-step methodological deep dives with references to the primary literature</p>
<p style="font-size:9px; color:#aaa; margin-top:80px;">All references verified against CrossRef \u00b7 Research/informational output, not clinical guidance</p></div>"""

labels = ["White Paper 1", "White Paper 1b", "White Paper 2", "White Paper 3", "White Paper 4",
          "White Paper 5", "White Paper 6", "White Paper 7", "White Paper 8", "White Paper 9", "White Paper 10"]
toc = "<div style='page-break-after:always;'><h1>Contents</h1><ol style='font-size:12px; line-height:2;'>"
for (name, t), lab in zip(titles, labels):
    toc += f"<li>{lab} \u2014 {t}</li>"
toc += "</ol></div>"

combined = title_html + toc
for name, _ in titles:
    n, html = render_paper(f"work/whitepapers/md/{name}.md")
    combined += f"<div style='page-break-before:always;'>{html}</div>"
fullc = f"<html><head><meta charset='utf-8'><style>{CSS}</style></head><body>{combined}</body></html>"
weasyprint.HTML(string=fullc, base_url="work/whitepapers").write_pdf(
    "work/whitepapers/pdf/00_ALL_whitepapers_combined.pdf")

# ── Build bibliography ────────────────────────────────────────────────────────
biblines = ["# Master Bibliography — Psoriasis Druggability Pipeline White Papers\n",
            "All references verified against CrossRef. Sorted alphabetically by citation key.\n"]
for k in sorted(master.keys()):
    s, doi = fmt_ref(master[k])
    s_plain = s.replace("<i>", "").replace("</i>", "")
    biblines.append(f"- **[{k}]** {s_plain} https://doi.org/{doi}")
open("work/whitepapers/BIBLIOGRAPHY.md", "w").write("\n".join(biblines))

# ── Bundle into tar.gz ────────────────────────────────────────────────────────
os.makedirs("work/whitepapers/bundle/markdown_sources", exist_ok=True)
for f in glob.glob("work/whitepapers/pdf/*.pdf"):
    shutil.copy(f, "work/whitepapers/bundle/")
shutil.copy("work/whitepapers/BIBLIOGRAPHY.md", "work/whitepapers/bundle/")
for f in glob.glob("work/whitepapers/md/*.md"):
    shutil.copy(f, "work/whitepapers/bundle/markdown_sources/")

with tarfile.open("whitepapers_psoriasis_pipeline.tar.gz", "w:gz") as t:
    t.add("work/whitepapers/bundle", arcname=".")

print("bundle:", os.path.getsize("whitepapers_psoriasis_pipeline.tar.gz") // 1024, "KB")
print("PDFs in bundle:", len(glob.glob("work/whitepapers/bundle/*.pdf")))