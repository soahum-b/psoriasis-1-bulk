# REPRODUCE — A3_pathway_TF

*14 artifacts staged from the artifact store into this directory. 13 have recoverable generating code in lineage.*

**Toolchain(s):** psoriasis-r (11), python (1)

**Files & provenance:**

| file | type | class | env | code_len |
|---|---|---|---|---|
| PATHWAY_METHODS_corrected.md | markdown | authored_doc | psoriasis-r | 0 |
| camera_old_vs_corrected.csv | csv | code_recoverable | psoriasis-r | 5667 |
| camera_perstudy_k4_full.rds | octet-stream | code_recoverable | psoriasis-r | 3690 |
| fig_pathway_corrected.png | png | code_recoverable | python | 6994 |
| gsea_mega_vs_meta_hallmark.csv | csv | code_recoverable | psoriasis-r | 1959 |
| gsea_meta_k4_extended.csv | csv | code_recoverable | psoriasis-r | 1698 |
| gsea_meta_k4_primary.csv | csv | code_recoverable | psoriasis-r | 1640 |
| gsea_meta_primary_jointBH.csv | csv | code_recoverable | psoriasis-r | 416 |
| gsea_perstudy_combined_k4.csv | csv | code_recoverable | psoriasis-r | 4261 |
| gsea_perstudy_raw_k4.csv | csv | code_recoverable | psoriasis-r | 3278 |
| meta_pathway_camera_k4_corrected.csv | csv | code_recoverable | psoriasis-r | 5210 |
| pathway_nonredundant_top.csv | csv | code_recoverable | psoriasis-r | 1868 |
| per_study_de_k4_primary.rds | octet-stream | code_recoverable | psoriasis-r | 2381 |
| plan_rebuild-the-psoriasis-pathway-enrichment_afa8da2d.json | json | code_recoverable | — | 461 |

**To regenerate a file:** fetch its generating code from artifact lineage (`host.lineage[version_id]['code']`), run it in the named conda env against the staged inputs. The version_id for each file is in `../../recoverability_audit_full.csv` (match on filename).