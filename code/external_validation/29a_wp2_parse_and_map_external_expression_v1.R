DLBCL_PROJECT_ROOT <- normalizePath(
  Sys.getenv("DLBCL_PROJECT_ROOT", unset = getwd()),
  winslash = "/", mustWork = FALSE
)
DLBCL_SUPPLEMENTARY_CODE_ROOT <- normalizePath(
  Sys.getenv("DLBCL_SUPPLEMENTARY_CODE_ROOT",
    unset = file.path(DLBCL_PROJECT_ROOT, "Supplementary_Code")),
  winslash = "/", mustWork = FALSE
)

# Node 29a: parse authoritative series matrices and explicitly join metadata.

if (!exists("WP2", inherits = TRUE)) source(file.path(DLBCL_SUPPLEMENTARY_CODE_ROOT, "external_validation", "29_common_wp2_functions_v1.R"))

wp2_node_29a <- function() {
  wp2_assert_token()
  specs <- list(
    GSE10846 = list(expression = WP2$paths$GSE10846_expression, metadata = WP2$paths$GSE10846_metadata, expected_n = 420L),
    GSE181063 = list(expression = WP2$paths$GSE181063_expression, metadata = WP2$paths$GSE181063_metadata, expected_n = 1310L)
  )
  audits <- list()
  sample_counts <- integer()
  for (dataset in names(specs)) {
    spec <- specs[[dataset]]
    parsed <- wp2_read_series_matrix(dataset, spec$expression, spec$metadata, spec$expected_n)
    wp2_save_rds(parsed, file.path(WP2$stage, "local_only", paste0(dataset, "_parsed_probe_expression.rds")))
    audits[[dataset]] <- parsed$audit
    sample_counts[[dataset]] <- ncol(parsed$expression)
    rm(parsed)
    gc()
  }
  wp2_write_csv(do.call(rbind, audits), file.path(WP2$stage, "01_parse_audit/WP2_EXPRESSION_PARSE_AUDIT.csv"))
  invisible(sample_counts)
}

if (!identical(Sys.getenv("WP2_NODE_SOURCE_ONLY", unset = ""), "TRUE")) wp2_node_29a()
