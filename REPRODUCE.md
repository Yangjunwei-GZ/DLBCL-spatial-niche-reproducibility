# Reproduce the current continuous-model workflow

## 1. Repository scope

This compendium documents the final six-program continuous model. It does not retain a discrete taxonomy. Final k is NOT_SELECTED and taxonomy is NOT_ASSIGNED.

## 2. Data access

Download raw or processed public data from GEO using the accessions in `DATA_SOURCES.md`. The repository itself supplies frozen derived inputs needed for figure inspection and rendering.

## 3. Environment setup

The validated R runtime was R 4.5.1. `environment/renv.lock` covers the core environment but does not contain every direct package used across the current scripts. Exact version evidence for the uncovered packages is retained in session and package records. Consequently, exact one-lock restoration is PARTIAL; do not treat `renv::restore()` alone as proof of the complete historical environment.

```r
renv::restore(lockfile = "environment/renv.lock")
```

Python was used for final supplementary-figure engineering. Exact cross-platform reconstruction of every rendering environment is not claimed; the frozen figures, source workbooks, scripts, and hashes are supplied.

## 4. Program definitions

The canonical authority is `source_data/figure_rendering/Figure_1/canonical_programs_v2.csv`: six programs, 22 memberships per program, 132 memberships, 121 unique genes, and 11 duplicated cross-program memberships.

## 5. Bulk scoring

Current GSE31312 analyses use the frozen 498 x 6 score authority. The primary score matrix is a validated transformation of a historical ssGSEA score space; the original historical package environment is not fully recoverable. Current score-freeze and validation scripts are under `code/continuous_model/`.

## 6. COO analyses

Use the frozen GSE31312 score matrix, COO metadata, PCA coordinates, and Kruskal-Wallis result under `source_data/figure_rendering/Figure_2/`. The current Figure 2 script is under `code/table_figure_generation/`.

## 7. Candidate-k adjudication

Use `source_data/figure_rendering/Figure_3/` and the current model-form adjudication scripts. Candidate k = 2-6 were checked against prespecified criteria; no candidate passed all criteria.

## 8. External continuous replication

The GSE10846 and GSE181063 frozen score/correlation inputs and integration summaries are under `source_data/figure_rendering/Figure_4/`. Workflow scripts are under `code/external_validation/`.

## 9. Single-cell contextualization

GSE182434 analyses use the DLBCL-only subset of 14,368 retained cells from four patients. Frozen UCell scores, UMAP coordinates, panel source tables, and canonical-gene coverage are under `source_data/single_cell/`. This is contextualization, not independent patient-level validation.

## 10. Spatial continuous analysis

GSE276542 workflow scripts are under `code/spatial/`. Final authority tables for spatial autocorrelation, bivariate association, method concordance, area scope, and Figure 6 are under `source_data/figure_rendering/Figure_6/`.

## 11. Spatial robustness analyses

Matched-random-gene-set summaries, depth-adjusted and residual-permutation outputs, matching-balance summaries, and direct-source/core sensitivity scripts are retained. Replicate-level matched membership partitions are intentionally excluded as redundant; seed and input registries remain under `environment/`.

## 12. LR expression-only support

Supplementary Figure S6 reports expression support and local co-occurrence only. It does not establish signaling, receptor activation, direct interaction, causality, mediation, or confirmed communication.

## 13. Figure generation

Final scripts for Figures 1-6 are under `code/table_figure_generation/`. Figure 1-5 repository-facing scripts use portable repository paths. Figure 6 consumes its frozen source authorities.

## 14. Supplementary figure generation

Current generators for S1-S6 are under `code/supplementary_figure_generation/`. The released canonical outputs, legends, source workbooks, QC records, and package manifests are under `supplementary_figures/S1` through `S6`.

## 15. Supplementary Tables

The sole release workbook is under `supplementary_tables/`. It contains 71 worksheets and no formulas.

## 16. SHA verification

Recompute SHA-256 for every file and compare it with `manifests/FINAL_GITHUB_REPOSITORY_SHA256_MANIFEST.csv`. `manifests/FINAL_GITHUB_MANIFEST_SELF_VERIFICATION.csv` records the release-time independent check.
