DLBCL_PROJECT_ROOT <- normalizePath(
  Sys.getenv("DLBCL_PROJECT_ROOT", unset = getwd()),
  winslash = "/", mustWork = FALSE
)
DLBCL_SUPPLEMENTARY_CODE_ROOT <- normalizePath(
  Sys.getenv("DLBCL_SUPPLEMENTARY_CODE_ROOT",
    unset = file.path(DLBCL_PROJECT_ROOT, "Supplementary_Code")),
  winslash = "/", mustWork = FALSE
)

# Node 29v: validate WP2 outputs, protected paths, scope controls, and token cleanup.

if (!exists("WP2", inherits = TRUE)) source(file.path(DLBCL_SUPPLEMENTARY_CODE_ROOT, "external_validation", "29_common_wp2_functions_v1.R"))

wp2_node_29v <- function() {
  checks <- list()
  check_id <- 0L
  add_check <- function(category, check_name, passed, observed, expected) {
    check_id <<- check_id + 1L
    checks[[check_id]] <<- data.frame(
      category = category,
      check_name = check_name,
      status = if (isTRUE(passed)) "PASS" else "FAIL",
      observed = as.character(observed),
      expected = as.character(expected),
      stringsAsFactors = FALSE
    )
  }

  add_check("GIT", "branch", identical(trimws(system("git branch --show-current", intern = TRUE)), WP2$branch), trimws(system("git branch --show-current", intern = TRUE)), WP2$branch)
  head <- trimws(system("git rev-parse HEAD", intern = TRUE))
  add_check("GIT", "baseline_HEAD", identical(head, WP2$baseline), head, WP2$baseline)
  input_gate <- tryCatch(wp2_verify_input_registry(), error = identity)
  add_check("INPUT", "all_registered_hashes_sizes_and_paths", !inherits(input_gate, "error"), if (inherits(input_gate, "error")) conditionMessage(input_gate) else nrow(input_gate), "13 PASS")

  protected <- tryCatch(wp2_verify_protected_paths(), error = identity)
  add_check("INTEGRITY", "protected_05o_to_05y", !inherits(protected, "error"), if (inherits(protected, "error")) conditionMessage(protected) else nrow(protected), "3341 unchanged files")
  files_05t <- list.files(file.path(WP2$root, "revision_2026_reviewer_response/05t_stage4c2a_v4_real_continuation_attempt"), recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  files_05t <- files_05t[file.info(files_05t)$isdir %in% FALSE]
  bytes_05t <- sum(file.info(files_05t)$size)
  add_check("INTEGRITY", "05t_file_count", length(files_05t) == 14L, length(files_05t), 14L)
  add_check("INTEGRITY", "05t_total_bytes", identical(as.numeric(bytes_05t), 71318), bytes_05t, 71318)

  parse <- utils::read.csv(file.path(WP2$stage, "01_parse_audit/WP2_EXPRESSION_PARSE_AUDIT.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  expected_samples <- c(GSE10846 = 420L, GSE181063 = 1310L)
  parse <- parse[match(names(expected_samples), parse$dataset_id), , drop = FALSE]
  add_check("PARSE", "sample_counts", identical(as.integer(parse$expression_sample_count), unname(expected_samples)), paste(parse$expression_sample_count, collapse = ";"), "420;1310")
  add_check("PARSE", "metadata_join_and_numeric_values", all(parse$explicit_join_complete & parse$numeric_parse_complete & parse$nonfinite_expression_count == 0L & parse$all_NA_sample_count == 0L), paste(parse$explicit_join_complete, collapse = ";"), "all TRUE")
  add_check("PARSE", "no_expression_transformation", all(parse$transformation_applied == "NONE"), paste(parse$transformation_applied, collapse = ";"), "NONE")

  coverage <- utils::read.csv(file.path(WP2$stage, "WP2_REAL_CANONICAL_COVERAGE.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  add_check("COVERAGE", "coverage_rows", nrow(coverage) == 24L, nrow(coverage), 24L)
  add_check("COVERAGE", "canonical_count", all(coverage$canonical_gene_count == 22L), paste(unique(coverage$canonical_gene_count), collapse = ";"), 22L)
  add_check("COVERAGE", "minimum_18_of_22", all(coverage$detected_canonical_gene_count >= 18L & coverage$coverage_gate == "PASS"), min(coverage$detected_canonical_gene_count), ">=18")
  substituted <- coverage$substituted_genes
  substituted[is.na(substituted)] <- ""
  add_check("COVERAGE", "no_substitution", all(!nzchar(substituted)), sum(nzchar(substituted)), 0L)
  add_check("COVERAGE", "program_order", all(vapply(split(coverage, interaction(coverage$dataset_id, coverage$collapse_rule)), function(x) identical(as.integer(x$program_order), 1:6), logical(1))), "all groups inspected", "1:6")

  for (dataset in names(expected_samples)) {
    mapping <- utils::read.csv(file.path(WP2$stage, "02_mapping", paste0(dataset, "_GENE_PROBE_COLLAPSE_AUDIT.csv")), check.names = FALSE, stringsAsFactors = FALSE)
    add_check("MAPPING", paste0(dataset, "_primary_rule"), all(mapping$primary_rule == "MEDIAN_ACROSS_PROBES_PER_GENE"), paste(unique(mapping$primary_rule), collapse = ";"), "MEDIAN_ACROSS_PROBES_PER_GENE")
    add_check("MAPPING", paste0(dataset, "_sensitivity_rule"), all(mapping$sensitivity_rule == "HIGHEST_MAD_PROBE_LEXICAL_TIEBREAK"), paste(unique(mapping$sensitivity_rule), collapse = ";"), "HIGHEST_MAD_PROBE_LEXICAL_TIEBREAK")
    add_check("MAPPING", paste0(dataset, "_no_outcome_or_program_selection"), all(!mapping$outcome_used & !mapping$program_membership_used_for_selection), "all rows inspected", "FALSE/FALSE")
    for (rule in c("PRIMARY", "SENSITIVITY")) {
      score_path <- file.path(WP2$stage, "03_scores", paste0(dataset, "_", rule, "_SIX_PROGRAM_SCORES.csv"))
      score <- utils::read.csv(score_path, check.names = FALSE, stringsAsFactors = FALSE)
      score_matrix <- as.matrix(score[, WP2$programs, drop = FALSE])
      add_check("SCORES", paste0(dataset, "_", rule, "_dimensions"), nrow(score) == expected_samples[[dataset]] && ncol(score) == 7L, paste(nrow(score), ncol(score), sep = "x"), paste(expected_samples[[dataset]], 7L, sep = "x"))
      add_check("SCORES", paste0(dataset, "_", rule, "_unique_samples"), !anyDuplicated(score$sample), anyDuplicated(score$sample), 0L)
      add_check("SCORES", paste0(dataset, "_", rule, "_finite_positive_SD"), all(is.finite(score_matrix)) && all(apply(score_matrix, 2, stats::sd) > 0), sum(!is.finite(score_matrix)), "0 nonfinite and all SD>0")
    }
  }

  for (dataset in c("GSE31312", "GSE10846", "GSE181063")) {
    path <- file.path(WP2$stage, "05_structure", paste0(dataset, ifelse(dataset == "GSE31312", "_PEARSON_CORRELATION_MATRIX.csv", "_PRIMARY_PEARSON_CORRELATION_MATRIX.csv")))
    m <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    values <- as.matrix(m[, WP2$programs, drop = FALSE])
    add_check("STRUCTURE", paste0(dataset, "_matrix"), nrow(m) == 6L && all(is.finite(values)) && max(abs(values - t(values))) < 1e-12, paste(nrow(m), max(abs(values - t(values))), sep = ";"), "6 rows; symmetric; finite")
  }
  point <- utils::read.csv(file.path(WP2$stage, "05_structure/WP2_STRUCTURAL_REPLICATION_POINT_ESTIMATES.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  boot_rep <- utils::read.csv(file.path(WP2$stage, "05_structure/WP2_STRUCTURAL_BOOTSTRAP_REPLICATES.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  boot_sum <- utils::read.csv(file.path(WP2$stage, "05_structure/WP2_STRUCTURAL_BOOTSTRAP_SUMMARY.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  add_check("STRUCTURE", "point_estimate_rows", nrow(point) == 2L && all(is.finite(point$matrix_vector_spearman)) && all(point$sign_concordance_n >= 0L & point$sign_concordance_n <= 15L), nrow(point), 2L)
  add_check("BOOTSTRAP", "replicate_count", nrow(boot_rep) == 4000L, nrow(boot_rep), 4000L)
  add_check("BOOTSTRAP", "all_replicates_successful", all(boot_rep$success) && all(boot_sum$successful_replicates == 2000L), sum(!boot_rep$success), 0L)
  add_check("BOOTSTRAP", "finite_intervals", all(is.finite(unlist(boot_sum[, grep("CI_", names(boot_sum)), drop = FALSE]))), "all interval columns inspected", "finite")

  strategy_b_cor <- utils::read.csv(file.path(WP2$stage, "06_strategy_b/WP2_STRATEGY_B_PROGRAM_PC_CORRELATIONS.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  strategy_b_var <- utils::read.csv(file.path(WP2$stage, "06_strategy_b/WP2_STRATEGY_B_PROJECTED_VARIANCE.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  add_check("STRATEGY_B", "program_PC_rows", nrow(strategy_b_cor) == 24L && all(is.finite(strategy_b_cor$pearson_r)), nrow(strategy_b_cor), 24L)
  add_check("STRATEGY_B", "projected_variance", nrow(strategy_b_var) == 4L && all(strategy_b_var$projected_variance > 0), paste(strategy_b_var$projected_variance, collapse = ";"), "4 positive values")
  for (dataset in names(expected_samples)) {
    coords <- utils::read.csv(file.path(WP2$stage, "06_strategy_b", paste0(dataset, "_STRATEGY_B_PC_COORDINATES.csv")), check.names = FALSE, stringsAsFactors = FALSE)
    add_check("STRATEGY_B", paste0(dataset, "_coordinates"), nrow(coords) == expected_samples[[dataset]] && all(is.finite(as.matrix(coords[, c("PC1_B", "PC2_B")]))), nrow(coords), expected_samples[[dataset]])
  }

  diag_a <- utils::read.csv(file.path(WP2$stage, "WP2_GSE10846_STRATEGY_A_SCALE_DIAGNOSTICS.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  status_a <- unique(diag_a$strategy_A_status)
  a_file <- file.path(WP2$stage, "07_strategy_a/GSE10846_STRATEGY_A_PC_COORDINATES.csv")
  a_correct <- length(status_a) == 1L && ((status_a == "PASS" && file.exists(a_file)) || (status_a == "BLOCKED" && !file.exists(a_file)))
  add_check("STRATEGY_A", "GSE10846_gate_and_conditional_file", a_correct, paste(status_a, file.exists(a_file), sep = ";"), "PASS+present or BLOCKED+absent")
  status_181063 <- utils::read.csv(file.path(WP2$stage, "07_strategy_a/GSE181063_STRATEGY_A_STATUS.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  add_check("STRATEGY_A", "GSE181063_not_authorized", identical(status_181063$strategy_A_status, "STRATEGY_A_NOT_AUTHORIZED_INCOMPLETE_CANONICAL_COVERAGE") && !status_181063$attempted, paste(status_181063$strategy_A_status, status_181063$attempted, sep = ";"), "not authorized; FALSE")

  clinical <- utils::read.csv(file.path(WP2$stage, "WP2_CLINICAL_EXECUTION_ISOLATION_AUDIT.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  logical_cols <- names(clinical)[vapply(clinical, is.logical, logical(1))]
  add_check("SCOPE", "clinical_isolation", nrow(clinical) == 2L && all(!unlist(clinical[, logical_cols, drop = FALSE])) && all(clinical$isolation_status == "PASS"), paste(nrow(clinical), length(logical_cols), sep = ";"), "2 rows; all execution flags FALSE")
  integration <- utils::read.csv(file.path(WP2$stage, "08_integration/WP2_CROSS_COHORT_INTEGRATION_SUMMARY.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  add_check("SCOPE", "final_k_not_selected", all(integration$final_k_status == "NOT_SELECTED"), paste(unique(integration$final_k_status), collapse = ";"), "NOT_SELECTED")
  add_check("SCOPE", "taxonomy_not_assigned", all(integration$taxonomy_status == "NOT_ASSIGNED"), paste(unique(integration$taxonomy_status), collapse = ";"), "NOT_ASSIGNED")

  runtime_attempt1 <- utils::read.csv(file.path(WP2$stage, "logs/WP2_RUNTIME_AUDIT.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  continuation_path <- file.path(WP2$stage, "logs/WP2_CONTINUATION_RUNTIME_AUDIT.csv")
  runtime_continuation <- if (file.exists(continuation_path)) utils::read.csv(continuation_path, check.names = FALSE, stringsAsFactors = FALSE) else runtime_attempt1[0, , drop = FALSE]
  runtime <- rbind(runtime_attempt1, runtime_continuation)
  add_check("RUNTIME", "warning_count", all(runtime$warning_count == 0L), sum(runtime$warning_count), 0L)
  add_check("RUNTIME", "stack_imbalance", all(!runtime$stack_imbalance), any(runtime$stack_imbalance), FALSE)
  completed_initial <- all(runtime_attempt1$status[match(c("29a", "29b", "29c", "29d", "29e", "29f"), runtime_attempt1$node)] == "COMPLETED")
  completed_continuation <- all(runtime_continuation$status[match(c("29g", "29h"), runtime_continuation$node)] == "COMPLETED")
  add_check("RUNTIME", "node_completion_across_attempt_and_continuation", completed_initial && completed_continuation, paste(c(runtime_attempt1$status[1:6], runtime_continuation$status), collapse = ";"), "29a-29h COMPLETED across frozen attempt plus continuation")
  add_check("RUNTIME", "GSVA_version", identical(as.character(utils::packageVersion("GSVA")), "2.4.9"), as.character(utils::packageVersion("GSVA")), "2.4.9")
  add_check("RUNTIME", "other_R_processes", identical(wp2_other_r_process_count(), 0L), wp2_other_r_process_count(), 0L)
  add_check("TOKEN", "execution_token_cleared", identical(Sys.getenv(WP2$token_name, unset = ""), ""), Sys.getenv(WP2$token_name, unset = "<UNSET>"), "UNSET")

  manifest <- utils::read.csv(file.path(WP2$stage, "WP2_REAL_OUTPUT_MANIFEST.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  pending_ids <- c("VALIDATOR_RESULTS", "VALIDATOR_REPORT", "FINAL_MANIFEST")
  presence <- vapply(seq_len(nrow(manifest)), function(i) {
    if (manifest$output_id[[i]] %in% pending_ids) return(TRUE)
    path <- file.path(WP2$stage, manifest$relative_path[[i]])
    if (isTRUE(manifest$conditional[[i]]) && status_a == "BLOCKED") return(!file.exists(path))
    file.exists(path)
  }, logical(1))
  add_check("OUTPUT", "manifested_outputs_present", all(presence), paste(manifest$output_id[!presence], collapse = ";"), "all required outputs present")
  existing_csv <- list.files(WP2$stage, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  forbidden <- tryCatch(wp2_assert_no_forbidden_columns(existing_csv), error = identity)
  add_check("SCOPE", "forbidden_scientific_columns_absent", !inherits(forbidden, "error"), if (inherits(forbidden, "error")) conditionMessage(forbidden) else length(existing_csv), "no forbidden columns")

  results <- do.call(rbind, checks)
  overall <- if (all(results$status == "PASS")) "PASS" else "FAIL"
  wp2_write_csv(results, file.path(WP2$stage, "WP2_FINAL_VALIDATOR_RESULTS.csv"))
  wp2_write_text(c(
    paste0("WP2_FINAL_VALIDATOR_STATUS=", overall),
    paste0("checks_total=", nrow(results)),
    paste0("checks_passed=", sum(results$status == "PASS")),
    paste0("checks_failed=", sum(results$status == "FAIL")),
    paste0("warning_count=", sum(runtime$warning_count)),
    paste0("stack_imbalance=", any(runtime$stack_imbalance)),
    paste0("execution_token_after_validation=", ifelse(nzchar(Sys.getenv(WP2$token_name, unset = "")), "SET", "UNSET")),
    paste0("final_k_status=NOT_SELECTED"),
    paste0("taxonomy_status=NOT_ASSIGNED")
  ), file.path(WP2$stage, "WP2_FINAL_VALIDATOR_REPORT.txt"))

  files <- list.files(WP2$stage, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  files <- files[file.info(files)$isdir %in% FALSE]
  files <- files[wp2_norm(files) != wp2_norm(file.path(WP2$stage, "WP2_FINAL_FILE_MANIFEST.csv"))]
  file_manifest <- do.call(rbind, lapply(sort(files), function(path) data.frame(
    relative_path = substring(wp2_norm(path), nchar(wp2_norm(WP2$stage)) + 2L),
    file_size_bytes = as.numeric(file.info(path)$size),
    sha256 = wp2_sha(path),
    tracking_policy = ifelse(grepl("^local_only/", substring(wp2_norm(path), nchar(wp2_norm(WP2$stage)) + 2L)), "LOCAL_ONLY_IGNORED", "TRACK"),
    stringsAsFactors = FALSE
  )))
  wp2_write_csv(file_manifest, file.path(WP2$stage, "WP2_FINAL_FILE_MANIFEST.csv"))
  if (overall != "PASS") stop("WP2 final validator failed", call. = FALSE)
  invisible(results)
}

if (!identical(Sys.getenv("WP2_NODE_SOURCE_ONLY", unset = ""), "TRUE")) wp2_node_29v()
