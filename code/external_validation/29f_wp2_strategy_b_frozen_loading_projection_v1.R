DLBCL_PROJECT_ROOT <- normalizePath(
  Sys.getenv("DLBCL_PROJECT_ROOT", unset = getwd()),
  winslash = "/", mustWork = FALSE
)
DLBCL_SUPPLEMENTARY_CODE_ROOT <- normalizePath(
  Sys.getenv("DLBCL_SUPPLEMENTARY_CODE_ROOT",
    unset = file.path(DLBCL_PROJECT_ROOT, "Supplementary_Code")),
  winslash = "/", mustWork = FALSE
)

# Node 29f: apply frozen WP1 PC loadings to cohort-relative standardized scores.

if (!exists("WP2", inherits = TRUE)) source(file.path(DLBCL_SUPPLEMENTARY_CODE_ROOT, "external_validation", "29_common_wp2_functions_v1.R"))

wp2_node_29f <- function() {
  wp2_assert_token()
  correlations <- list()
  variances <- list()
  for (dataset in c("GSE10846", "GSE181063")) {
    scores <- utils::read.csv(file.path(WP2$stage, "03_scores", paste0(dataset, "_PRIMARY_SIX_PROGRAM_SCORES.csv")), check.names = FALSE, stringsAsFactors = FALSE)
    result <- wp2_strategy_b(scores, dataset)
    wp2_write_csv(result$coordinates, file.path(WP2$stage, "06_strategy_b", paste0(dataset, "_STRATEGY_B_PC_COORDINATES.csv")))
    correlations[[dataset]] <- result$correlations
    variances[[dataset]] <- result$variance
  }
  wp2_write_csv(do.call(rbind, correlations), file.path(WP2$stage, "06_strategy_b/WP2_STRATEGY_B_PROGRAM_PC_CORRELATIONS.csv"))
  wp2_write_csv(do.call(rbind, variances), file.path(WP2$stage, "06_strategy_b/WP2_STRATEGY_B_PROJECTED_VARIANCE.csv"))
  invisible(TRUE)
}

if (!identical(Sys.getenv("WP2_NODE_SOURCE_ONLY", unset = ""), "TRUE")) wp2_node_29f()
