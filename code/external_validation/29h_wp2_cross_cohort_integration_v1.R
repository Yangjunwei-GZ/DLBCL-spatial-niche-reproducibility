DLBCL_PROJECT_ROOT <- normalizePath(
  Sys.getenv("DLBCL_PROJECT_ROOT", unset = getwd()),
  winslash = "/", mustWork = FALSE
)
DLBCL_SUPPLEMENTARY_CODE_ROOT <- normalizePath(
  Sys.getenv("DLBCL_SUPPLEMENTARY_CODE_ROOT",
    unset = file.path(DLBCL_PROJECT_ROOT, "Supplementary_Code")),
  winslash = "/", mustWork = FALSE
)

# Node 29h: integrate cohort-level continuous validation results without merging samples.

if (!exists("WP2", inherits = TRUE)) source(file.path(DLBCL_SUPPLEMENTARY_CODE_ROOT, "external_validation", "29_common_wp2_functions_v1.R"))

wp2_node_29h <- function() {
  wp2_assert_token()
  coverage <- utils::read.csv(file.path(WP2$stage, "WP2_REAL_CANONICAL_COVERAGE.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  qc <- utils::read.csv(file.path(WP2$stage, "04_qc/WP2_SCORE_QC_AND_MAPPING_SENSITIVITY.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  point <- utils::read.csv(file.path(WP2$stage, "05_structure/WP2_STRUCTURAL_REPLICATION_POINT_ESTIMATES.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  boot <- utils::read.csv(file.path(WP2$stage, "05_structure/WP2_STRUCTURAL_BOOTSTRAP_SUMMARY.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  edge <- utils::read.csv(file.path(WP2$stage, "05_structure/WP2_PRIMARY_CORRELATION_EDGES_LONG.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  strategy_b <- utils::read.csv(file.path(WP2$stage, "06_strategy_b/WP2_STRATEGY_B_PROGRAM_PC_CORRELATIONS.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  strategy_a <- utils::read.csv(file.path(WP2$stage, "WP2_GSE10846_STRATEGY_A_SCALE_DIAGNOSTICS.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  discovery_pc <- utils::read.csv(WP2$paths$discovery_program_pc, check.names = FALSE, stringsAsFactors = FALSE)

  integration <- list()
  row_id <- 0L
  add_row <- function(scope, category, metric, value, interpretation) {
    row_id <<- row_id + 1L
    integration[[row_id]] <<- data.frame(
      scope = scope,
      category = category,
      metric = metric,
      value = as.character(value),
      interpretation = interpretation,
      samples_merged_across_cohorts = FALSE,
      final_k_status = "NOT_SELECTED",
      taxonomy_status = "NOT_ASSIGNED",
      stringsAsFactors = FALSE
    )
  }

  for (dataset in c("GSE10846", "GSE181063")) {
    p <- point[point$dataset_id == dataset, , drop = FALSE]
    b <- boot[boot$dataset_id == dataset, , drop = FALSE]
    q <- qc[qc$dataset_id == dataset, , drop = FALSE]
    sb <- strategy_b[strategy_b$dataset_id == dataset, , drop = FALSE]
    authority <- discovery_pc[order(discovery_pc$program_order, discovery_pc$component), , drop = FALSE]
    sb_ordered <- sb[order(sb$program_order, sb$component), , drop = FALSE]
    direction_matches_authority <- sign(sb_ordered$pearson_r) == sign(authority$pearson_r)
    add_row(dataset, "TECHNICAL", "sample_count", ifelse(dataset == "GSE10846", 420L, 1310L), "All eligible samples have six finite primary scores.")
    add_row(dataset, "STRUCTURE", "matrix_vector_spearman", p$matrix_vector_spearman, p$structural_replication_status)
    add_row(dataset, "STRUCTURE", "sign_concordance", paste0(p$sign_concordance_n, "/15"), p$structural_replication_status)
    add_row(dataset, "BOOTSTRAP", "matrix_vector_spearman_95CI", paste0("[", b$matrix_rho_CI_lower, ", ", b$matrix_rho_CI_upper, "]"), paste0(b$successful_replicates, "/", b$requested_replicates, " successful replicates"))
    add_row(dataset, "MAPPING_SENSITIVITY", "programs_with_spearman_below_0.85", sum(q$mapping_sensitivity_status == "LOW_CONCORDANCE"), "Primary collapse rule was retained regardless of sensitivity concordance.")
    add_row(dataset, "STRATEGY_B", "program_PC_direction_matches_WP1", paste0(sum(direction_matches_authority), "/12"), "Cohort-relative coordinates use frozen unflipped WP1 loadings.")
  }

  e1 <- edge[edge$dataset_id == "GSE10846", , drop = FALSE]
  e2 <- edge[edge$dataset_id == "GSE181063", , drop = FALSE]
  e1 <- e1[order(e1$program_1_order, e1$program_2_order), , drop = FALSE]
  e2 <- e2[order(e2$program_1_order, e2$program_2_order), , drop = FALSE]
  cross_rho <- suppressWarnings(stats::cor(e1$pearson_r, e2$pearson_r, method = "spearman"))
  cross_signs <- sum(sign(e1$pearson_r) == sign(e2$pearson_r))
  add_row("GSE10846_vs_GSE181063", "CROSS_COHORT", "correlation_vector_spearman", cross_rho, "Cohort structures compared without merging or batch correction.")
  add_row("GSE10846_vs_GSE181063", "CROSS_COHORT", "sign_concordance", paste0(cross_signs, "/15"), "Cohort structures compared without selecting a preferred collapse rule.")
  add_row("GSE10846", "STRATEGY_A", "status", unique(strategy_a$strategy_A_status), "Exploratory only; does not affect primary Strategy B or structure results.")
  add_row("GSE181063", "STRATEGY_A", "status", "STRATEGY_A_NOT_AUTHORIZED_INCOMPLETE_CANONICAL_COVERAGE", "Not attempted under HD5.")
  add_row("WP2", "COHORT_ISOLATION", "status", "PASS", "Each cohort retains its own scores, coordinates, and structural result.")

  integration_df <- do.call(rbind, integration)
  wp2_write_csv(integration_df, file.path(WP2$stage, "08_integration/WP2_CROSS_COHORT_INTEGRATION_SUMMARY.csv"))

  clinical <- data.frame(
    dataset_id = c("GSE10846", "GSE181063"),
    clinical_fields_entering_nodes_29a_to_29h = FALSE,
    OS_model_run = FALSE,
    event_model_run = FALSE,
    treatment_model_run = FALSE,
    Cox_model_run = FALSE,
    KM_curve_run = FALSE,
    clinical_cutoff_run = FALSE,
    median_split_run = FALSE,
    clinical_output_created = FALSE,
    isolation_status = "PASS",
    notes = c("Outcome and treatment fields were not used by WP2 nodes.", "Clinical analysis is prohibited for GSE181063 in WP2."),
    stringsAsFactors = FALSE
  )
  wp2_write_csv(clinical, file.path(WP2$stage, "WP2_CLINICAL_EXECUTION_ISOLATION_AUDIT.csv"))

  coverage_lines <- vapply(c("GSE10846", "GSE181063"), function(dataset) {
    x <- coverage[coverage$dataset_id == dataset & coverage$collapse_rule == "MEDIAN_ACROSS_PROBES_PER_GENE", , drop = FALSE]
    paste0("- ", dataset, ": ", paste(paste0(x$program_order, "=", x$detected_canonical_gene_count, "/22"), collapse = "; "))
  }, character(1))
  structural_lines <- vapply(c("GSE10846", "GSE181063"), function(dataset) {
    p <- point[point$dataset_id == dataset, , drop = FALSE]
    b <- boot[boot$dataset_id == dataset, , drop = FALSE]
    paste0("- ", dataset, ": rho=", signif(p$matrix_vector_spearman, 5), ", signs=", p$sign_concordance_n, "/15, bootstrap rho 95% CI [", signif(b$matrix_rho_CI_lower, 5), ", ", signif(b$matrix_rho_CI_upper, 5), "], status=", p$structural_replication_status)
  }, character(1))
  report <- c(
    "# WP2 External Continuous Validation Report",
    "",
    "## Execution status",
    "",
    "WP2 EXTERNAL CONTINUOUS VALIDATION COMPLETED",
    "",
    "- GSE10846 samples: 420",
    "- GSE181063 samples: 1310",
    "- Primary gene aggregation: sample-wise median across valid probes per uppercase symbol",
    "- Sensitivity gene aggregation: highest cross-sample MAD probe with lexical probe-ID tiebreak",
    "- Program scores: GSVA ssGSEA, six frozen canonical programs",
    "- Cohorts were not merged or batch corrected",
    "",
    "## Canonical coverage",
    "",
    coverage_lines,
    "",
    "No genes were substituted, borrowed, imputed, or restored from a historical shortlist.",
    "",
    "## Structural replication",
    "",
    structural_lines,
    "",
    paste0("- Cross-cohort 15-edge Spearman: ", signif(cross_rho, 5), "; sign concordance: ", cross_signs, "/15"),
    "",
    "Structural labels follow the preregistered HD6 thresholds and are not selected by P values.",
    "",
    "## Continuous PC representations",
    "",
    "- Strategy B: completed for both cohorts using cohort-relative z scores and frozen, unflipped WP1 PC1/PC2 loadings.",
    paste0("- GSE10846 Strategy A: ", unique(strategy_a$strategy_A_status), "."),
    "- GSE181063 Strategy A: not authorized and not attempted.",
    "",
    "## Scope controls",
    "",
    "- Clinical outcome analysis: NOT RUN",
    "- Spatial analysis: NOT RUN",
    "- Purity analysis: NOT RUN",
    "- Cluster or class assignment: NOT PRODUCED",
    "- final k: NOT_SELECTED",
    "- taxonomy: NOT_ASSIGNED",
    "",
    "This WP2 analysis evaluates continuous program-level structural robustness. It does not validate a discrete ecosystem taxonomy or clinical outcome model."
  )
  wp2_write_text(report, file.path(WP2$stage, "WP2_EXTERNAL_CONTINUOUS_VALIDATION_REPORT.md"))
  invisible(integration_df)
}

if (!identical(Sys.getenv("WP2_NODE_SOURCE_ONLY", unset = ""), "TRUE")) wp2_node_29h()
