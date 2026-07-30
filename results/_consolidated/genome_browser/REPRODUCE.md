# REPRODUCE — E_genome_browser

*12 artifacts staged from the artifact store into this directory. 10 have recoverable generating code in lineage.*

**Toolchain(s):** psoriasis-r (8), python (2)

**Files & provenance:**

| file | type | class | env | code_len |
|---|---|---|---|---|
| GENOME_BROWSER_README.md | markdown | authored_doc | python | 0 |
| PSORIASIS_PROJECT_NARRATIVE.md | markdown | code_recoverable | python | 19770 |
| PSORIASIS_PROJECT_NARRATIVE.pdf | pdf | code_recoverable | python | 21595 |
| RESUME_genome_browser.md | markdown | authored_doc | python | 0 |
| module_il1_inflammasome.png | png | code_recoverable | psoriasis-r | 2828 |
| module_il36.png | png | code_recoverable | psoriasis-r | 2740 |
| module_jakstat.png | png | code_recoverable | psoriasis-r | 2834 |
| module_th17_output.png | png | code_recoverable | psoriasis-r | 2924 |
| stat3_exon23_zoom.png | png | code_recoverable | psoriasis-r | 2141 |
| stat3_headline.png | png | code_recoverable | psoriasis-r | 4840 |
| stat3_lead.png | png | code_recoverable | psoriasis-r | 1612 |
| stat3_quant_panel.png | png | code_recoverable | psoriasis-r | 1858 |

**To regenerate a file:** fetch its generating code from artifact lineage (`host.lineage[version_id]['code']`), run it in the named conda env against the staged inputs. The version_id for each file is in `../../recoverability_audit_full.csv` (match on filename).