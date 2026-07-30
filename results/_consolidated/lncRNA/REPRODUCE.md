# REPRODUCE — I_lncRNA

*11 artifacts staged from the artifact store into this directory. 11 have recoverable generating code in lineage.*

**Toolchain(s):** psoriasis-r (8), python (2)

**Files & provenance:**

| file | type | class | env | code_len |
|---|---|---|---|---|
| coexpr_lnc_target_adjusted.csv | csv | code_recoverable | psoriasis-r | 5269 |
| coexpr_lnc_target_unadjusted.csv | csv | code_recoverable | psoriasis-r | 5393 |
| fig_lncRNA_volcano_summary.png | png | code_recoverable | python | 6317 |
| gencode_v26_annotation.csv | csv | code_recoverable | psoriasis-r | 415 |
| lncRNA_curated_STAT3_IL1_linkage.csv | csv | code_recoverable | python | 5433 |
| lncRNA_significant_ranked.csv | csv | code_recoverable | psoriasis-r | 2412 |
| lncRNA_target_linkage.csv | csv | code_recoverable | psoriasis-r | 6517 |
| logcpm_rebuild_summary.csv | csv | code_recoverable | psoriasis-r | 2613 |
| meta_de_PPvsNN_4study_biotyped.csv | csv | code_recoverable | psoriasis-r | 1706 |
| plan_full-lncrna-arm-for-the-psoriasis-pp-vs_9961df5f.json | json | code_recoverable | — | 251 |
| symbol_biotype_map.csv | csv | code_recoverable | psoriasis-r | 1318 |

**To regenerate a file:** fetch its generating code from artifact lineage (`host.lineage[version_id]['code']`), run it in the named conda env against the staged inputs. The version_id for each file is in `../../recoverability_audit_full.csv` (match on filename).