# Reproduce the current continuous-model workflow

## 1. Repository scope

This compendium documents the final six-program continuous model. It does not retain a discrete taxonomy. Final k is NOT_SELECTED and taxonomy is NOT_ASSIGNED.

## 2. Data access

Download raw or processed public data from GEO using the accessions in `DATA_SOURCES.md`. The repository supplies released derived inputs for scientific inspection and verification.

## 3. Environment setup

The validated R runtime was R 4.5.1. `environment/renv.lock` records 16 of the 28 direct R packages required by the public scripts. Exact versions for the remaining 12 packages are documented in `environment/README.md` and `environment/package_versions.csv` from retained execution records. Exact one-lock restoration is therefore partial and is not claimed.

```r
renv::restore(lockfile = "environment/renv.lock")
```

Python was used for final supplementary-figure engineering. Exact cross-platform reconstruction of every rendering environment is not claimed.

The retained figure scripts use configurable project and source-data roots. Some historical final-render scripts also expect raw GEO files, project-derived authorities, or journal-distributed supplementary material that is not duplicated here. These scripts are retained for provenance; the repository is not presented as a one-command publication-artwork build.

Use the relative paths in `source_data/SOURCE_DATA_INDEX.csv` to locate repository-distributed authorities. Entries explicitly marked as not included identify frozen archival or journal-distributed material and are not runtime claims.

## 4. Program definitions

The recorded membership table is `source_data/figure_rendering/Figure_1/canonical_programs_v2.csv`: six literature-informed curated composite programs, 22 memberships per program, 132 memberships, 121 unique genes, and 11 shared cross-program memberships. The filename is retained as an internal provenance identifier and does not denote a validated biological taxonomy.

## 5. Bulk scoring

Current GSE31312 primary analyses use the validated 498 x 6 historical score matrix. The primary score matrix is a validated transformation of the recorded historical ssGSEA score space; the original package environment cannot be reconstructed exactly. Consequently, exact regeneration of the historical matrix from raw expression is not claimed. Harmonized sensitivity scoring used GSVA 2.4.9 in ssGSEA mode across GSE31312 (n=498), GSE10846 unsorted bulk samples (n=414), and GSE181063 (n=1,310). Current summaries are under `source_data/final_sensitivity_summaries/`, and workflow scripts are under `code/continuous_model/`.

## 6. COO analyses

Use the validated GSE31312 score matrix, COO metadata, PCA coordinates, and Kruskal-Wallis result under `source_data/figure_rendering/Figure_2/`. The current Figure 2 script is under `code/table_figure_generation/`.

## 7. Candidate-k adjudication

Use `source_data/figure_rendering/Figure_3/`, `source_data/final_sensitivity_summaries/k2_k6_retention_and_continuity.csv`, and the current model-form adjudication scripts. The retention rules were prespecified before model-form adjudication. No candidate k = 2-6 was sufficiently stable under the full rule set; a possible k=2 structure remained stability-limited. Final k was not selected and no discrete taxonomy was assigned. Figure 3 is reproduced from the released derived inputs.

## 8. External continuous replication

The GSE10846 and GSE181063 score/correlation inputs and integration summaries are under `source_data/figure_rendering/Figure_4/`. The historical GSE10846 primary analysis used all 420 profiles; harmonized and de-overlap sensitivity analyses used the 414 unsorted bulk samples. Workflow scripts remain under `code/external_validation/` for provenance, and current public summaries are under `source_data/final_sensitivity_summaries/`.

## 9. Single-cell contextualization

GSE182434 analyses use the DLBCL-only subset of 14,368 retained cells from four patients. Validated UCell scores, UMAP coordinates, panel source tables, and gene-coverage records are under `source_data/single_cell/`. This is contextualization, not independent patient-level validation.

## 10. Spatial continuous analysis

GSE276542 workflow scripts are under `code/spatial/`. Final authority tables for spatial autocorrelation, bivariate association, method concordance, area scope, and Figure 6 are under `source_data/figure_rendering/Figure_6/`.

## 11. Spatial robustness analyses

Current matched-null, depth-adjusted, residual-permutation, matching-balance, and direct-source/core scripts are retained. Compact repository summaries are provided in `source_data/final_sensitivity_summaries/spatial_sensitivity_denominators.csv`; complete submission-facing tables are distributed with the journal materials. Matched-null support (42/54 evaluable combinations overall; 25/32 primary-plus-exploratory combinations) is distinct from depth-adjusted residual support (30/30 DLBCL-area combinations). Large detailed intermediates and replicate-level matched membership partitions are reserved for the associated archival release.

## 12. LR expression-only support

Supplementary Figure S6 reports expression support and local co-occurrence only. It does not establish signaling, receptor activation, direct interaction, causality, mediation, or confirmed communication.

## 13. Figure generation

Final scripts for Figures 1-6 are under `code/table_figure_generation/`. They are retained as provenance and may write generated outputs to caller-selected output directories. Publication-ready binaries are distributed by the journal and are not repository inputs.

## 14. Supplementary figure generation

Current generators for S1-S6 are under `code/supplementary_figure_generation/`. Several generators consume project-derived or journal-distributed authorities not duplicated in this repository; they are retained for provenance and methodological inspection.

## 15. Supplementary Tables

The compiled submission-facing Supplementary Tables workbook is distributed with the journal submission/publication and is not duplicated in this repository. Repository-distributed machine-readable authorities are indexed in `source_data/SOURCE_DATA_INDEX.csv`.

## 16. SHA verification

Recompute SHA-256 for every released file and compare it with `manifests/GITHUB_RELEASE_SHA256.csv`. `manifests/GITHUB_RELEASE_SHA256_SELF_CHECK.csv` records the independent release-time verification.
