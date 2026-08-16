# Reproduce the current continuous-model workflow

## 1. Repository scope

This compendium documents the final six-program continuous model. It does not retain a discrete taxonomy. Final k is NOT_SELECTED and taxonomy is NOT_ASSIGNED.

## 2. Data access

Download raw or processed public data from GEO using the accessions in `DATA_SOURCES.md`. The repository itself supplies released derived inputs needed for figure inspection and rendering.

## 3. Environment setup

The validated R runtime was R 4.5.1. `environment/renv.lock` records 16 of the 28 direct R packages required by the public scripts. Exact versions for the remaining 12 packages are documented in `environment/README.md` and `environment/package_versions.csv` from retained execution records. Exact one-lock restoration is therefore partial and is not claimed.

```r
renv::restore(lockfile = "environment/renv.lock")
```

Python was used for final supplementary-figure engineering. Exact cross-platform reconstruction of every rendering environment is not claimed; released figures and source workbooks remain directly inspectable.

The current figure scripts use three configurable roots: the repository root, `source_data/`, and `supplementary_tables/`. Some retained analysis and final-render scripts also expect raw GEO files or project-derived authorities that are not redistributed. For those components, use the released derived authority rather than treating the package as a one-command raw-to-figure build.

The source workbooks are immutable scientific artifacts. Embedded provenance-path strings record where analysis inputs were produced; they are not required to locate public inputs and do not replace the relative paths in `source_data/SOURCE_DATA_INDEX.csv`.

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

Current matched-null, depth-adjusted, residual-permutation, matching-balance, and direct-source/core scripts are retained. Their display-ready values are released in the S5 source workbook, the final Supplementary Tables, and `source_data/final_sensitivity_summaries/spatial_sensitivity_denominators.csv`. Matched-null support (42/54 evaluable combinations overall; 25/32 primary-plus-exploratory combinations) is distinct from depth-adjusted residual support (30/30 DLBCL-area combinations). Large detailed intermediates and replicate-level matched membership partitions are reserved for the associated Zenodo archival release.

## 12. LR expression-only support

Supplementary Figure S6 reports expression support and local co-occurrence only. It does not establish signaling, receptor activation, direct interaction, causality, mediation, or confirmed communication.

## 13. Figure generation

Final scripts for Figures 1-6 are under `code/table_figure_generation/`. Figures 1-5 use repository-configurable roots. Figure 6 consumes released project-derived inputs and is therefore documented as a derived-input workflow.

## 14. Supplementary figure generation

Current generators for S1-S6 are under `code/supplementary_figure_generation/`. The public package retains one final submission TIFF, one legend, and one source workbook for each figure. Several generators consume project-derived inputs not redistributed in the lightweight GitHub package; their released source workbooks provide the public reproduction inputs.

## 15. Supplementary Tables

The sole release workbook is [`supplementary_tables/DLBCL_continuous_model_Supplementary_Tables_FRONTIERS_FINAL_SUBMISSION_READY.xlsx`](supplementary_tables/DLBCL_continuous_model_Supplementary_Tables_FRONTIERS_FINAL_SUBMISSION_READY.xlsx). It contains 71 worksheets and no formulas.

## 16. SHA verification

Recompute SHA-256 for every released file and compare it with `manifests/GITHUB_RELEASE_SHA256.csv`. `manifests/GITHUB_RELEASE_SHA256_SELF_CHECK.csv` records the independent release-time verification.
