DLBCL_PROJECT_ROOT <- normalizePath(
  Sys.getenv("DLBCL_PROJECT_ROOT", unset = getwd()),
  winslash = "/", mustWork = FALSE
)
DLBCL_SUPPLEMENTARY_CODE_ROOT <- normalizePath(
  Sys.getenv("DLBCL_SUPPLEMENTARY_CODE_ROOT",
    unset = file.path(DLBCL_PROJECT_ROOT, "Supplementary_Code")),
  winslash = "/", mustWork = FALSE
)

# Node 29c: compute six continuous GSVA ssGSEA scores for both collapse rules.

if (!exists("WP2", inherits = TRUE)) source(file.path(DLBCL_SUPPLEMENTARY_CODE_ROOT, "external_validation", "29_common_wp2_functions_v1.R"))

wp2_node_29c <- function() {
  wp2_assert_token()
  coverage <- utils::read.csv(file.path(WP2$stage, "WP2_REAL_CANONICAL_COVERAGE.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  rules <- c(primary = "MEDIAN_ACROSS_PROBES_PER_GENE", sensitivity = "HIGHEST_MAD_PROBE_PER_GENE")
  expected_n <- c(GSE10846 = 420L, GSE181063 = 1310L)
  for (dataset in names(expected_n)) {
    for (rule_id in names(rules)) {
      expr_path <- file.path(WP2$stage, "local_only", paste0(dataset, "_", ifelse(rule_id == "primary", "primary_median", "sensitivity_highest_MAD"), "_gene_expression.rds"))
      expr <- readRDS(expr_path)
      cov_rows <- coverage[coverage$dataset_id == dataset & coverage$collapse_rule == rules[[rule_id]], , drop = FALSE]
      score <- wp2_score_programs(expr, cov_rows)
      if (nrow(score) != expected_n[[dataset]] || ncol(score) != 7L) stop("Unexpected score dimensions", call. = FALSE)
      out_name <- paste0(dataset, "_", toupper(rule_id), "_SIX_PROGRAM_SCORES.csv")
      wp2_write_csv(score, file.path(WP2$stage, "03_scores", out_name))
      rm(expr, score)
      gc()
    }
  }
  invisible(TRUE)
}

if (!identical(Sys.getenv("WP2_NODE_SOURCE_ONLY", unset = ""), "TRUE")) wp2_node_29c()
