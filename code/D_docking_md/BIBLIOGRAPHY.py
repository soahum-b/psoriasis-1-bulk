# Auto-extracted generating script
# Produces: BIBLIOGRAPHY.md
# Conda env: dock-md   (run with this environment activated)
# Inputs (expected alongside / in data/): master_refs.json
# Source artifact version: 07c58133-db08-48c3-b848-dd32646d28fe
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import json
import re

master = json.load(open("master_refs.json"))

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

biblines = ["# Master Bibliography — Psoriasis Druggability Pipeline White Papers\n",
            "All references verified against CrossRef. Sorted alphabetically by citation key.\n"]
for k in sorted(master.keys()):
    s, doi = fmt_ref(master[k])
    s_plain = s.replace("<i>", "").replace("</i>", "")
    biblines.append(f"- **[{k}]** {s_plain} https://doi.org/{doi}")
open("BIBLIOGRAPHY.md", "w").write("\n".join(biblines))