# REPRODUCE — A4_study_exclusion

*23 artifacts staged from the artifact store into this directory. 23 have recoverable generating code in lineage.*

**Toolchain(s):** psoriasis-r (23)

**Files & provenance:**

| file | type | class | env | code_len |
|---|---|---|---|---|
| fig32_stat3_forest_7study.png | png | code_recoverable | psoriasis-r | 5012 |
| fig33_robustness_7study.png | png | code_recoverable | psoriasis-r | 6386 |
| fig34_stat3_isoform_junction_7study.png | png | code_recoverable | psoriasis-r | 8625 |
| fig35_isoform_mega_vs_meta.png | png | code_recoverable | psoriasis-r | 7588 |
| fig36_isoform_pvalue_combination.png | png | code_recoverable | psoriasis-r | 3456 |
| fig37_stat3_gradient_isoform_vs_gene.png | png | code_recoverable | psoriasis-r | 3322 |
| fig38_stat3_vs_module_gradient.png | png | code_recoverable | psoriasis-r | 3420 |
| fig39_druggable_targets_gradient.png | png | code_recoverable | psoriasis-r | 3505 |
| key_genes_7study_PPvsNN.csv | csv | code_recoverable | psoriasis-r | 3464 |
| meta_7study_summary.csv | csv | code_recoverable | psoriasis-r | 3647 |
| meta_de_PPvsNN_7study.csv | csv | code_recoverable | psoriasis-r | 3319 |
| meta_de_results_7study.rds | octet-stream | code_recoverable | psoriasis-r | 3187 |
| per_study_de_7study.rds | octet-stream | code_recoverable | psoriasis-r | 2212 |
| stat3_5v7study_comparison.csv | csv | code_recoverable | psoriasis-r | 3503 |
| stat3_gradient_isoform_vs_gene.csv | csv | code_recoverable | psoriasis-r | 3848 |
| stat3_isoform_junction_7study.csv | csv | code_recoverable | psoriasis-r | 5394 |
| stat3_isoform_mega_vs_meta.csv | csv | code_recoverable | psoriasis-r | 6315 |
| stat3_isoform_pool_5v7.csv | csv | code_recoverable | psoriasis-r | 6219 |
| stat3_isoform_pvalue_combination.csv | csv | code_recoverable | psoriasis-r | 1490 |
| stat3_isoform_stouffer_inputs.csv | csv | code_recoverable | psoriasis-r | 851 |
| stat3_vs_module_gradient_position.csv | csv | code_recoverable | psoriasis-r | 1693 |
| targets_gradient_position.csv | csv | code_recoverable | psoriasis-r | 783 |
| targets_jak_tyk2_abundance.csv | csv | code_recoverable | psoriasis-r | 897 |

**To regenerate a file:** fetch its generating code from artifact lineage (`host.lineage[version_id]['code']`), run it in the named conda env against the staged inputs. The version_id for each file is in `../../recoverability_audit_full.csv` (match on filename).