# REPRODUCE — C_stat3_isoforms

*3 artifacts staged from the artifact store into this directory. 2 have recoverable generating code in lineage.*

**Toolchain(s):** python (2)

**Files & provenance:**

| file | type | class | env | code_len |
|---|---|---|---|---|
| fig40_stat3_isoform_formal_evaluation.png | png | code_recoverable | python | 8323 |
| stat3_isoform_ab_formal_section.md | markdown | authored_doc | python | 0 |
| stat3_isoform_ab_formal_table.csv | csv | code_recoverable | python | 1671 |

**To regenerate a file:** fetch its generating code from artifact lineage (`host.lineage[version_id]['code']`), run it in the named conda env against the staged inputs. The version_id for each file is in `../../recoverability_audit_full.csv` (match on filename).