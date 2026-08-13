DLBCL_PROJECT_ROOT <- normalizePath(
  Sys.getenv("DLBCL_PROJECT_ROOT", unset = getwd()),
  winslash = "/", mustWork = FALSE
)
DLBCL_SUPPLEMENTARY_CODE_ROOT <- normalizePath(
  Sys.getenv("DLBCL_SUPPLEMENTARY_CODE_ROOT",
    unset = file.path(DLBCL_PROJECT_ROOT, "Supplementary_Code")),
  winslash = "/", mustWork = FALSE
)

source(file.path(DLBCL_SUPPLEMENTARY_CODE_ROOT, "spatial", "wp3_common_v3.R"))

wp3_node_36f <- function(state) {
  method_scores <- list(PRIMARY_LOGNORMALIZE_UCELL = state$primary_scores,
    SENSITIVITY_SCT_UCELL = state$sct_scores, SENSITIVITY_LOGNORMALIZE_SSGSEA = state$ssgsea_scores)
  state$score_qc <- do.call(rbind, lapply(names(method_scores), function(method) {
    table <- method_scores[[method]]
    do.call(rbind, lapply(setdiff(names(table), "barcode"), function(program) data.frame(
      capture_area_id = state$capture_area_id, method = method, program_id = program,
      n_spots = nrow(table), finite_n = sum(is.finite(table[[program]])),
      minimum = min(table[[program]], na.rm = TRUE), median = median(table[[program]], na.rm = TRUE),
      maximum = max(table[[program]], na.rm = TRUE), stringsAsFactors = FALSE)))
  }))
  state$concordance <- rbind(
    wp3_score_concordance(state$primary_scores, state$sct_scores, state$capture_area_id, "PRIMARY_vs_SCT_UCell"),
    wp3_score_concordance(state$primary_scores, state$ssgsea_scores, state$capture_area_id, "PRIMARY_vs_ssGSEA"))
  state
}

if (sys.nframe() == 0L) wp3_node_cli("36f", wp3_node_36f)
