# Reproduce the current continuous-model workflow

## 1. Repository scope

This compendium documents the final six-program continuous model. It does not retain a discrete taxonomy. Final k is NOT_SELECTED and taxonomy is NOT_ASSIGNED.

## 2. Data access

Download raw or processed public data from GEO using the accessions in `DATA_SOURCES.md`. The repository itself supplies frozen derived inputs needed for figure inspection and rendering.

## 3. Environment setup

The validated R runtime was R 4.5.1. `environment/renv.lock` covers the core environment but not every direct package used across the current scripts. Additional version evidence is summarized in `environment/README.md` and `environment/package_versions.csv`. Exact one-lock restoration is therefore PARTIAL.

```r
renv::restore(lockfile = "environment/renv.lock")
```

Python was used for final supplementary-figure engineering. Exact cross-platform reconstruction of every rendering environment is not claimed; released figures and source workbooks remain directly inspectable.

The current figure scripts use three configurable roots: the repository root, `source_data/`, and `supplementary_tables/`. Some retained analysis and final-render scripts also expect raw GEO files or project-derived authorities that are not redistributed. For those components, use the released derived authority rather than treating the package as a one-command raw-to-figure build.

The source workbooks are immutable scientific artifacts. Embedded provenance-path strings record where frozen inputs were produced; they are not required to locate public inputs and do not replace the relative paths in `source_data/SOURCE_DATA_INDEX.csv`.

## 4. Program definitions

The canonical public definition is `source_data/figure_rendering/Figure_1/canonical_programs_v2.csv`: six programs, 22 memberships per program, 132 memberships, 121 unique genes, and 11 duplicated cross-program memberships.

## 5. Bulk scoring

Current GSE31312 analyses use the frozen 498 x 6 score authority. The primary score matrix is a validated transformation of a historical ssGSEA score space; the original historical package environment is not fully recoverable. Current score-freeze and validation scripts are under `code/continuous_model/`.

## 6. COO analyses

Use the frozen GSE31312 score matrix, COO metadata, PCA coordinates, and Kruskal-Wallis result under `source_data/figure_rendering/Figure_2/`. The current Figure 2 script is under `code/table_figure_generation/`.

## 7. Candidate-k adjudication

Use `source_data/figure_rendering/Figure_3/` and the current model-form adjudication scripts. Candidate k = 2-6 were checked against prespecified criteria; no candidate passed all criteria. Figure 3 is reproduced from the released derived authority.

## 8. External continuous replication

The GSE10846 and GSE181063 frozen score/correlation inputs and integration summaries are under `source_data/figure_rendering/Figure_4/`. Workflow scripts are under `code/external_validation/`.

## 9. Single-cell contextualization

GSE182434 analyses use the DLBCL-only subset of 14,368 retained cells from four patients. Frozen UCell scores, UMAP coordinates, panel source tables, and canonical-gene coverage are under `source_data/single_cell/`. This is contextualization, not independent patient-level validation.

## 10. Spatial continuous analysis

GSE276542 workflow scripts are under `code/spatial/`. Final authority tables for spatial autocorrelation, bivariate association, method concordance, area scope, and Figure 6 are under `source_data/figure_rendering/Figure_6/`.

## 11. Spatial robustness analyses

Current matched-null, depth-adjusted, residual-permutation, matching-balance, and direct-source/core scripts are retained. Their display-ready values are released in the S5 source workbook and final Supplementary Tables. Large detailed intermediates and replicate-level matched membership partitions are reserved for the associated Zenodo archival release.

## 12. LR expression-only support

Supplementary Figure S6 reports expression support and local co-occurrence only. It does not establish signaling, receptor activation, direct interaction, causality, mediation, or confirmed communication.

## 13. Figure generation

Final scripts for Figures 1-6 are under `code/table_figure_generation/`. Figures 1-5 use repository-configurable roots. Figure 6 consumes frozen project-derived authorities and is therefore documented as a released-derived-authority workflow.

## 14. Supplementary figure generation

Current generators for S1-S6 are under `code/supplementary_figure_generation/`. The public package retains one PDF, one legend, and one source workbook for each figure. Several generators consume project-derived authorities not redistributed in the lightweight GitHub package; their released source workbooks are the reproduction authority.

## 15. Supplementary Tables

The sole release workbook is [`supplementary_tables/DLBCL_continuous_model_Supplementary_Tables_FRONTIERS_FINAL_SUBMISSION_READY.xlsx`](supplementary_tables/DLBCL_continuous_model_Supplementary_Tables_FRONTIERS_FINAL_SUBMISSION_READY.xlsx). It contains 71 worksheets and no formulas.

## 16. SHA verification

Recompute SHA-256 for every released file and compare it with `manifests/GITHUB_RELEASE_SHA256.csv`. `manifests/GITHUB_RELEASE_SHA256_SELF_CHECK.csv` records the independent release-time verification.
