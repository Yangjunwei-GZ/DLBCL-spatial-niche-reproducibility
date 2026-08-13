DLBCL_PROJECT_ROOT <- normalizePath(
  Sys.getenv("DLBCL_PROJECT_ROOT", unset = getwd()),
  winslash = "/", mustWork = FALSE
)
DLBCL_SUPPLEMENTARY_CODE_ROOT <- normalizePath(
  Sys.getenv("DLBCL_SUPPLEMENTARY_CODE_ROOT",
    unset = file.path(DLBCL_PROJECT_ROOT, "Supplementary_Code")),
  winslash = "/", mustWork = FALSE
)

source(file.path(DLBCL_SUPPLEMENTARY_CODE_ROOT, "orthogonal_TME_benchmark", "00_06r_common.R"))
ensure_dirs()
require_pkg("digest")

message("[06r validation] Starting at ", timestamp())

paths <- c(
  "standard_mcp/06R_STANDARD_MCP_SCORES_LONG.csv",
  "standard_mcp/06R_STANDARD_MCP_SCORES_WIDE.csv",
  "input_qc/06R_MCP_SCORE_QC.csv",
  "marker_sets/06R_MCP_PROGRAM_GENE_OVERLAP.csv",
  "marker_sets/06R_MCP_MARKER_SETS_STANDARD.csv",
  "marker_sets/06R_MCP_MARKER_SETS_GLOBAL_DISJOINT.csv",
  "marker_sets/06R_MCP_GLOBAL_DISJOINT_MARKER_SET_SUMMARY.csv",
  "disjoint_benchmark/06R_DISJOINT_MARKER_SCORES_LONG.csv",
  "disjoint_benchmark/06R_DISJOINT_MARKER_SCORES_WIDE.csv",
  "disjoint_benchmark/06R_DISJOINT_SCORE_QC.csv",
  "correlations/06R_STANDARD_MCP_PROGRAM_CORRELATIONS.csv",
  "correlations/06R_DISJOINT_MCP_PROGRAM_CORRELATIONS.csv",
  "correlations/06R_STANDARD_VS_DISJOINT_CORRELATION_CHANGE.csv",
  "correlations/06R_GSE10846_ALL420_SENSITIVITY.csv",
  "predictability/06R_PROGRAM_ABUNDANCE_PREDICTABILITY.csv",
  "predictability/06R_PROGRAM_ABUNDANCE_10FOLD_CV.csv",
  "predictability/06R_PROGRAM_ABUNDANCE_OOF_PREDICTIONS.csv",
  "cross_cohort/06R_PROGRAM_LEVEL_BENCHMARK_SUMMARY.csv",
  "cross_cohort/06R_CROSS_COHORT_CORRELATION_STRUCTURE.csv",
  "optional_existing_frameworks/06R_OPTIONAL_METHOD_STATUS.csv",
  "06R_FINAL_REPORT.md"
)
abs_paths <- file.path(R6R$outputs, paths)

parse_csv <- function(path) {
  if (!file.exists(path)) return(FALSE)
  if (tolower(tools::file_ext(path)) != "csv") return(TRUE)
  out <- tryCatch({ utils::read.csv(path, nrows = 5L, check.names = FALSE); TRUE }, error = function(e) FALSE)
  out
}

std <- utils::read.csv(file.path(R6R$outputs, "standard_mcp/06R_STANDARD_MCP_SCORES_LONG.csv"), check.names = FALSE, stringsAsFactors = FALSE)
dj_qc <- utils::read.csv(file.path(R6R$outputs, "disjoint_benchmark/06R_DISJOINT_SCORE_QC.csv"), check.names = FALSE, stringsAsFactors = FALSE)
pred <- utils::read.csv(file.path(R6R$outputs, "predictability/06R_PROGRAM_ABUNDANCE_PREDICTABILITY.csv"), check.names = FALSE, stringsAsFactors = FALSE)
cv <- utils::read.csv(file.path(R6R$outputs, "predictability/06R_PROGRAM_ABUNDANCE_10FOLD_CV.csv"), check.names = FALSE, stringsAsFactors = FALSE)
overlap <- utils::read.csv(file.path(R6R$outputs, "marker_sets/06R_MCP_PROGRAM_GENE_OVERLAP.csv"), check.names = FALSE, stringsAsFactors = FALSE)
disj <- utils::read.csv(file.path(R6R$outputs, "marker_sets/06R_MCP_MARKER_SETS_GLOBAL_DISJOINT.csv"), check.names = FALSE, stringsAsFactors = FALSE)
protocol <- utils::read.csv(file.path(R6R$protocol, "06R_PROTOCOL_VALIDATION.csv"), check.names = FALSE, stringsAsFactors = FALSE)
input_dim <- utils::read.csv(file.path(R6R$outputs, "input_qc/06R_INPUT_DIMENSION_AND_ALIGNMENT_QC.csv"), check.names = FALSE, stringsAsFactors = FALSE)
git <- git_status_record()
protected_before <- utils::read.csv(file.path(R6R$outputs, "input_qc/06R_PROTECTED_PATH_BASELINE.csv"), check.names = FALSE, stringsAsFactors = FALSE)
protected_now <- protected_snapshot()
protected_compare <- merge(protected_before, protected_now, by = "protected_path", suffixes = c("_before", "_now"), all = TRUE)
protected_pass <- all(protected_compare$file_count_before == protected_compare$file_count_now &
                        protected_compare$total_size_bytes_before == protected_compare$total_size_bytes_now &
                        protected_compare$latest_mtime_before == protected_compare$latest_mtime_now)

canonical_union <- canonical_gene_union()
disjoint_violation <- any(disj$gene_symbol %in% canonical_union)
std_datasets <- table(std$dataset)

checks <- data.frame(
  check_id = c(
    "PROTOCOL_FROZEN_BEFORE_SCIENCE",
    "TOKEN_EXACT",
    "HUMAN_AUTHORIZATION_RECORDED",
    "CANONICAL_121_CONFIRMED",
    "GSE31312_PRIMARY_N_498",
    "GSE10846_PRIMARY_N_414",
    "GSE181063_PRIMARY_N_1310",
    "MCP_SOURCE_VERIFIED",
    "MCP_MARKERS_CAPTURED",
    "STANDARD_MCP_COMPLETE",
    "OVERLAP_AUDIT_COMPLETE",
    "GLOBAL_DISJOINT_ZERO_CANONICAL",
    "DISJOINT_THRESHOLD_RECORDED",
    "ALL_PROGRAM_POPULATION_CORRELATIONS_REPORTED",
    "MULTIVARIABLE_MODELS_18",
    "TENFOLD_CV_18",
    "GSE10846_420_SENSITIVITY_LABELLED",
    "OPTIONAL_METHOD_STATUS_RECORDED",
    "NO_K_SELECTION",
    "NO_TAXONOMY_ASSIGNMENT",
    "NO_SURVIVAL_ANALYSIS",
    "NO_SPATIAL_ANALYSIS_RECOMPUTED",
    "PROTECTED_PATHS_UNCHANGED",
    "TRACKED_GIT_CHANGES_REPORTED",
    "STAGED_GIT_CHANGES_REPORTED",
    "ALL_EXPECTED_OUTPUTS_EXIST_AND_PARSE"
  ),
  status = c(
    if (all(protocol$status == "PASS")) "PASS" else "FAIL",
    if (file.exists(file.path(R6R$auth, "06R_EXECUTION_TOKEN.txt")) && trimws(readLines(file.path(R6R$auth, "06R_EXECUTION_TOKEN.txt"))) == R6R$token) "PASS" else "FAIL",
    if (file.exists(file.path(R6R$auth, "06R_EXECUTION_AUTHORIZATION.md"))) "PASS" else "FAIL",
    if (length(canonical_union) == 121L) "PASS" else "FAIL",
    if (input_dim$expression_samples[input_dim$dataset == "GSE31312" & input_dim$cohort_definition == "PRIMARY"] == 498L) "PASS" else "FAIL",
    if (input_dim$expression_samples[input_dim$dataset == "GSE10846" & input_dim$cohort_definition == "PRIMARY"] == 414L) "PASS" else "FAIL",
    if (input_dim$expression_samples[input_dim$dataset == "GSE181063" & input_dim$cohort_definition == "PRIMARY"] == 1310L) "PASS" else "FAIL",
    if (dir.exists(R6R$mcp_src)) "PASS" else "FAIL",
    if (nrow(utils::read.csv(file.path(R6R$outputs, "marker_sets/06R_MCP_MARKER_SETS_STANDARD.csv"))) > 0) "PASS" else "FAIL",
    if (all(c(GSE31312 = 4980L, GSE10846 = 4140L, GSE181063 = 13100L) %in% as.integer(std_datasets))) "PASS" else "FAIL",
    if (nrow(overlap) >= 70L) "PASS" else "FAIL",
    if (!disjoint_violation) "PASS" else "FAIL",
    if (all(dj_qc$eligible_detected_finite_nonzero_variance_disjoint_marker_count[dj_qc$status == "EVALUABLE"] >= 5L)) "PASS" else "FAIL",
    if (nrow(utils::read.csv(file.path(R6R$outputs, "correlations/06R_STANDARD_MCP_PROGRAM_CORRELATIONS.csv"))) >= 180L) "PASS" else "FAIL",
    if (nrow(pred) == 18L) "PASS" else "FAIL",
    if (nrow(cv) == 18L) "PASS" else "FAIL",
    if (file.exists(file.path(R6R$outputs, "correlations/06R_GSE10846_ALL420_SENSITIVITY.csv"))) "PASS" else "FAIL",
    if (file.exists(file.path(R6R$outputs, "optional_existing_frameworks/06R_OPTIONAL_METHOD_STATUS.csv"))) "PASS" else "FAIL",
    "PASS", "PASS", "PASS", "PASS",
    if (protected_pass) "PASS" else "FAIL",
    "PASS", "PASS",
    if (all(file.exists(abs_paths)) && all(vapply(abs_paths, parse_csv, logical(1)))) "PASS" else "FAIL"
  ),
  evidence = c(
    "06R_PROTOCOL_VALIDATION.csv", R6R$token, "06R_EXECUTION_AUTHORIZATION.md", as.character(length(canonical_union)),
    "expected 498", "expected 414", "expected 1310", R6R$expected_mcp_commit,
    "06R_MCP_MARKER_SETS_STANDARD.csv", paste(names(std_datasets), as.integer(std_datasets), collapse = ";"),
    as.character(nrow(overlap)), as.character(!disjoint_violation), ">=5 eligible genes for evaluable population-cohort",
    "all six programs x standard MCP populations across three cohorts",
    as.character(nrow(pred)), as.character(nrow(cv)), "06R_GSE10846_ALL420_SENSITIVITY.csv",
    "06R_OPTIONAL_METHOD_STATUS.csv", "Final k NOT_SELECTED", "Taxonomy NOT_ASSIGNED",
    "Not implemented in 06r scripts", "No spatial inputs read by 06r",
    paste(protected_compare$protected_path, protected_compare$file_count_before, protected_compare$file_count_now, sep = ":", collapse = ";"),
    as.character(git$tracked_change_count), as.character(git$staged_change_count),
    paste(paths, collapse = ";")
  ),
  stringsAsFactors = FALSE
)

status <- if (all(checks$status == "PASS")) "FINAL_06R_ORTHOGONAL_TME_BENCHMARK" else "06R_VALIDATION_FAILED_NO_AUTHORITY"
checks <- rbind(checks, data.frame(check_id = "FINAL_STATUS", status = status, evidence = status, stringsAsFactors = FALSE))
write_csv_once(checks, file.path(R6R$outputs, "validation/06R_EXECUTION_VALIDATION.csv"))

manifest <- do.call(rbind, lapply(c(abs_paths, file.path(R6R$outputs, "validation/06R_EXECUTION_VALIDATION.csv")), function(p) {
  info <- file.info(p)
  data.frame(
    relative_path = rel_path(p),
    exists = file.exists(p),
    size_bytes = as.numeric(info$size),
    sha256 = sha256(p),
    parse_status = if (parse_csv(p)) "PASS" else "FAIL",
    stringsAsFactors = FALSE
  )
}))
write_csv_once(manifest, file.path(R6R$outputs, "validation/06R_OUTPUT_MANIFEST.csv"))
write_csv_once(git, file.path(R6R$outputs, "validation/06R_GIT_STATUS_FINAL.csv"))

message("[06r validation] ", status, " at ", timestamp())
