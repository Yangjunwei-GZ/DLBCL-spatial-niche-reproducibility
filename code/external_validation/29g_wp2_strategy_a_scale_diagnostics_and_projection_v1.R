DLBCL_PROJECT_ROOT <- normalizePath(
  Sys.getenv("DLBCL_PROJECT_ROOT", unset = getwd()),
  winslash = "/", mustWork = FALSE
)
DLBCL_SUPPLEMENTARY_CODE_ROOT <- normalizePath(
  Sys.getenv("DLBCL_SUPPLEMENTARY_CODE_ROOT",
    unset = file.path(DLBCL_PROJECT_ROOT, "Supplementary_Code")),
  winslash = "/", mustWork = FALSE
)

# Node 29g: enforce the complete HD4 scale gate for exploratory GSE10846 Strategy A.

if (!exists("WP2", inherits = TRUE)) source(file.path(DLBCL_SUPPLEMENTARY_CODE_ROOT, "external_validation", "29_common_wp2_functions_v1.R"))

wp2_node_29g <- function() {
  wp2_assert_token()
  primary <- utils::read.csv(file.path(WP2$stage, "03_scores/GSE10846_PRIMARY_SIX_PROGRAM_SCORES.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  sensitivity <- utils::read.csv(file.path(WP2$stage, "03_scores/GSE10846_SENSITIVITY_SIX_PROGRAM_SCORES.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  coverage <- utils::read.csv(file.path(WP2$stage, "WP2_REAL_CANONICAL_COVERAGE.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  coverage <- coverage[coverage$dataset_id == "GSE10846" & coverage$collapse_rule == "MEDIAN_ACROSS_PROBES_PER_GENE", , drop = FALSE]
  coverage <- coverage[order(coverage$program_order), , drop = FALSE]
  diagnostics <- wp2_strategy_a_diagnostics(primary, sensitivity, coverage)
  out <- diagnostics$diagnostics
  out$standardized_matrix_finite_pass <- diagnostics$matrix_pass
  out$covariance_matrix_finite_pass <- diagnostics$covariance_pass
  out$all_HD4_gates_pass <- diagnostics$status == "PASS"
  out$strategy_A_status <- diagnostics$status
  wp2_write_csv(out, file.path(WP2$stage, "WP2_GSE10846_STRATEGY_A_SCALE_DIAGNOSTICS.csv"))

  if (identical(diagnostics$status, "PASS")) {
    coords <- diagnostics$z %*% wp2_loadings()$matrix
    coord_df <- data.frame(
      sample = primary$sample,
      PC1_A = coords[, 1],
      PC2_A = coords[, 2],
      coordinate_scale = "GSE31312_TRAINING_MEAN_SD",
      analysis_role = "EXPLORATORY_ABSOLUTE_SCALE_SENSITIVITY",
      stringsAsFactors = FALSE
    )
    wp2_write_csv(coord_df, file.path(WP2$stage, "07_strategy_a/GSE10846_STRATEGY_A_PC_COORDINATES.csv"))
  }

  status_181063 <- data.frame(
    dataset_id = "GSE181063",
    strategy_A_status = "STRATEGY_A_NOT_AUTHORIZED_INCOMPLETE_CANONICAL_COVERAGE",
    attempted = FALSE,
    coordinate_file_created = FALSE,
    reason = "Amendment 023 HD5 prohibits Strategy A for GSE181063.",
    stringsAsFactors = FALSE
  )
  wp2_write_csv(status_181063, file.path(WP2$stage, "07_strategy_a/GSE181063_STRATEGY_A_STATUS.csv"))
  invisible(diagnostics$status)
}

if (!identical(Sys.getenv("WP2_NODE_SOURCE_ONLY", unset = ""), "TRUE")) wp2_node_29g()
