# REPRODUCE — D_docking_md

*50 artifacts staged from the artifact store into this directory. 50 have recoverable generating code in lineage.*

**Toolchain(s):** dock-md (38), python (10)

**Files & provenance:**

| file | type | class | env | code_len |
|---|---|---|---|---|
| 00_ALL_whitepapers_combined.pdf | pdf | code_recoverable | dock-md | 9232 |
| 01_target_list.pdf | pdf | code_recoverable | dock-md | 5019 |
| 01b_ppi_network.pdf | pdf | code_recoverable | dock-md | 5041 |
| 02_druggability_annotation.pdf | pdf | code_recoverable | dock-md | 3406 |
| 03_scoring_leads.pdf | pdf | code_recoverable | dock-md | 3147 |
| 04_structures.pdf | pdf | code_recoverable | dock-md | 3156 |
| 05_pockets.pdf | pdf | code_recoverable | dock-md | 3092 |
| 06_library.pdf | pdf | code_recoverable | dock-md | 4745 |
| 07_docking.pdf | pdf | code_recoverable | dock-md | 3254 |
| 08_cofolding.pdf | pdf | code_recoverable | dock-md | 2560 |
| 09_md.pdf | pdf | code_recoverable | dock-md | 3112 |
| 10_report.pdf | pdf | code_recoverable | dock-md | 3055 |
| BIBLIOGRAPHY.md | markdown | code_recoverable | dock-md | 1171 |
| CONSOLIDATED_REPORT.md | markdown | code_recoverable | dock-md | 7338 |
| JAK3_5LWM.pdb | x-pdb | code_recoverable | dock-md | 133 |
| LIBRARY_README.md | markdown | code_recoverable | dock-md | 2027 |
| RORC_5APH.pdb | x-pdb | code_recoverable | dock-md | 133 |
| RORC_tophit_complex.pdb | x-pdb | code_recoverable | dock-md | 515 |
| STAT3_6NJS.pdb | x-pdb | code_recoverable | dock-md | 134 |
| all_docked_poses.tar.gz | gzip | code_recoverable | dock-md | 4327 |
| docking_ranked.csv | csv | code_recoverable | dock-md | 938 |
| druggability_annotation.csv | csv | code_recoverable | python | 3902 |
| druggable_target_ranking.csv | csv | code_recoverable | python | 2890 |
| fig_docking_results.png | png | code_recoverable | dock-md | 4032 |
| fig_druggability_landscape.png | png | code_recoverable | python | 4121 |
| fig_lead_structures.png | png | code_recoverable | dock-md | 4479 |
| fig_library_properties.png | png | code_recoverable | dock-md | 3405 |
| fig_pipeline_funnel.png | png | code_recoverable | dock-md | 3758 |
| fig_pockets_druggability.png | png | code_recoverable | dock-md | 7614 |
| fig_ppi_network.png | png | code_recoverable | python | 7090 |
| fig_target_ranking.png | png | code_recoverable | python | 10711 |
| fig_target_volcano_annotated.png | png | code_recoverable | python | 3799 |
| fig_tophit_pose.png | png | code_recoverable | dock-md | 2005 |
| lead_network_centrality.csv | csv | code_recoverable | python | 1855 |
| lead_pipeline_summary.csv | csv | code_recoverable | dock-md | 1498 |
| lead_rationale.md | markdown | code_recoverable | python | 2888 |
| plan_druggability-triage-of-psoriasis-meta-de_bf412257.json | json | code_recoverable | — | 567 |
| plan_druggability-triage-of-psoriasis-meta-de_bf412257.json | json | code_recoverable | — | 2287 |
| pocket_druggability_summary.csv | csv | code_recoverable | dock-md | 3758 |
| pocket_report.csv | csv | code_recoverable | dock-md | 9056 |
| ppi_hub_genes.csv | csv | code_recoverable | python | 2090 |
| psoriasis_cofolding_package.tar.gz | gzip | code_recoverable | dock-md | 8013 |
| psoriasis_docking_package.tar.gz | gzip | code_recoverable | dock-md | 291 |
| psoriasis_md_package.tar.gz | gzip | code_recoverable | dock-md | 13950 |
| screening_library.csv | csv | code_recoverable | dock-md | 2051 |
| screening_library_3d.sdf | x-mdl-sdfile | code_recoverable | dock-md | 2651 |
| structure_inventory.csv | csv | code_recoverable | dock-md | 3158 |
| target_universe.csv | csv | code_recoverable | python | 2354 |
| validation.json | json | code_recoverable | dock-md | 1878 |
| whitepapers_psoriasis_pipeline.tar.gz | gzip | code_recoverable | dock-md | 15978 |

**To regenerate a file:** fetch its generating code from artifact lineage (`host.lineage[version_id]['code']`), run it in the named conda env against the staged inputs. The version_id for each file is in `../../recoverability_audit_full.csv` (match on filename).