# REPRODUCE — A2_peer_review

*4 artifacts staged from the artifact store into this directory. 2 have recoverable generating code in lineage.*

**Toolchain(s):** psoriasis-r (2)

**Files & provenance:**

| file | type | class | env | code_len |
|---|---|---|---|---|
| fig31_stat3_hksj_loo.png | png | code_recoverable | psoriasis-r | 4809 |
| psoriasis_integrated_whitepaper.md | markdown | authored_doc | python | 0 |
| stat3_forest_robustness.csv | csv | code_recoverable | psoriasis-r | 1874 |
| stat3_robustness_section_11_5a.md | markdown | authored_doc | python | 0 |

**To regenerate a file:** fetch its generating code from artifact lineage (`host.lineage[version_id]['code']`), run it in the named conda env against the staged inputs. The version_id for each file is in `../../recoverability_audit_full.csv` (match on filename).