DLBCL_PROJECT_ROOT <- normalizePath(
  Sys.getenv("DLBCL_PROJECT_ROOT", unset = getwd()),
  winslash = "/", mustWork = FALSE
)
DLBCL_SUPPLEMENTARY_CODE_ROOT <- normalizePath(
  Sys.getenv("DLBCL_SUPPLEMENTARY_CODE_ROOT",
    unset = file.path(DLBCL_PROJECT_ROOT, "Supplementary_Code")),
  winslash = "/", mustWork = FALSE
)

# Node 29d: score QC and primary-versus-highest-MAD mapping sensitivity.

if (!exists("WP2", inherits = TRUE)) source(file.path(DLBCL_SUPPLEMENTARY_CODE_ROOT, "external_validation", "29_common_wp2_functions_v1.R"))

wp2_node_29d <- function() {
  wp2_assert_token()
  rows <- lapply(c("GSE10846", "GSE181063"), function(dataset) {
    primary <- utils::read.csv(file.path(WP2$stage, "03_scores", paste0(dataset, "_PRIMARY_SIX_PROGRAM_SCORES.csv")), check.names = FALSE, stringsAsFactors = FALSE)
    sensitivity <- utils::read.csv(file.path(WP2$stage, "03_scores", paste0(dataset, "_SENSITIVITY_SIX_PROGRAM_SCORES.csv")), check.names = FALSE, stringsAsFactors = FALSE)
    if (!identical(primary$sample, sensitivity$sample)) stop("Primary/sensitivity sample mismatch", call. = FALSE)
    wp2_score_qc(dataset, primary, sensitivity)
  })
  qc <- do.call(rbind, rows)
  wp2_write_csv(qc, file.path(WP2$stage, "04_qc/WP2_SCORE_QC_AND_MAPPING_SENSITIVITY.csv"))
  invisible(qc)
}

if (!identical(Sys.getenv("WP2_NODE_SOURCE_ONLY", unset = ""), "TRUE")) wp2_node_29d()
