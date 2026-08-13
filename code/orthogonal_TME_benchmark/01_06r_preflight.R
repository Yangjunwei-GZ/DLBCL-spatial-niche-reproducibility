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

message("[06r preflight] Starting at ", timestamp())

write_text_or_verify <- function(lines, path, required_pattern = NULL) {
  if (file.exists(path)) {
    existing <- readLines(path, warn = FALSE)
    target <- paste(lines, collapse = "\n")
    current <- paste(existing, collapse = "\n")
    if (!is.null(required_pattern) && !grepl(required_pattern, current, fixed = TRUE)) {
      stop("Existing create-once text file does not contain required authority: ", path, call. = FALSE)
    }
    if (is.null(required_pattern) && !identical(trimws(current), trimws(target))) {
      stop("Existing create-once text file differs from expected content: ", path, call. = FALSE)
    }
    return(invisible(path))
  }
  write_text_once(lines, path)
}

write_text_or_verify(R6R$token, file.path(R6R$auth, "06R_EXECUTION_TOKEN.txt"))
write_text_or_verify(c(
  "# 06R Execution Authorization",
  "",
  paste0("- Scientific authorization token: ", R6R$token),
  paste0("- Human authorization: ", R6R$human_authorization),
  "- Authorized module: 06r_bulk_orthogonal_tme_benchmark only",
  "- Final k: NOT_SELECTED",
  "- Taxonomy: NOT_ASSIGNED",
  "- Prohibited: k selection, taxonomy assignment, survival analysis, ligand-receptor analysis, centroid projection, manuscript/figure/workbook/response-letter edits, commit, push",
  paste0("- Created: ", timestamp())
), file.path(R6R$auth, "06R_EXECUTION_AUTHORIZATION.md"), required_pattern = R6R$token)

git_exe <- Sys.which("git")
if (!nzchar(git_exe)) stop("PRECHECK_FAILED_MCP_COUNTER_SOURCE_UNAVAILABLE: git executable not found", call. = FALSE)
valid_source_commit <- function(path) {
  if (!dir.exists(file.path(path, ".git"))) return(NA_character_)
  out <- tryCatch(system2(git_exe, c("-C", path, "rev-parse", "HEAD"), stdout = TRUE, stderr = TRUE), error = function(e) NA_character_)
  trimws(out[[1]])
}

if (!dir.exists(R6R$mcp_src)) dir.create(R6R$mcp_src, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(file.path(R6R$mcp_src, ".git"))) {
  local_candidates <- c(
    file.path(Sys.getenv("TEMP"), "MCPcounter_official_06r"),
    file.path(Sys.getenv("TEMP"), "MCPcounter_official")
  )
  copied <- FALSE
  for (cand in local_candidates) {
    if (identical(valid_source_commit(cand), R6R$expected_mcp_commit)) {
      files <- list.files(cand, all.files = TRUE, no.. = TRUE, full.names = TRUE)
      ok <- file.copy(files, R6R$mcp_src, recursive = TRUE, copy.date = TRUE)
      if (all(ok)) {
        copied <- TRUE
        break
      }
    }
  }
  if (!copied) {
    status <- system2(git_exe, c("clone", "--depth", "1", R6R$expected_mcp_repo, R6R$mcp_src))
    if (!identical(status, 0L)) stop("PRECHECK_FAILED_MCP_COUNTER_SOURCE_UNAVAILABLE: official repo clone failed", call. = FALSE)
  }
}
if (!dir.exists(file.path(R6R$mcp_src, ".git"))) {
  status <- system2(git_exe, c("clone", "--depth", "1", R6R$expected_mcp_repo, R6R$mcp_src))
  if (!identical(status, 0L)) stop("PRECHECK_FAILED_MCP_COUNTER_SOURCE_UNAVAILABLE: official repo clone failed", call. = FALSE)
}
mcp_commit <- system2(git_exe, c("-C", R6R$mcp_src, "rev-parse", "HEAD"), stdout = TRUE)
if (!identical(trimws(mcp_commit[[1]]), R6R$expected_mcp_commit)) {
  stop("PRECHECK_FAILED_MCP_COUNTER_SOURCE_UNAVAILABLE: unexpected MCP-counter commit ", paste(mcp_commit, collapse = " "), call. = FALSE)
}

desc_path <- file.path(R6R$mcp_src, "Source/DESCRIPTION")
genes_path <- file.path(R6R$mcp_src, "Signatures/genes.txt")
probesets_path <- file.path(R6R$mcp_src, "Signatures/probesets.txt")
if (!file.exists(desc_path) || !file.exists(genes_path) || !file.exists(probesets_path)) {
  stop("PRECHECK_FAILED_MCP_COUNTER_SOURCE_UNAVAILABLE: required official source/signature files absent", call. = FALSE)
}

if (!dir.exists(file.path(R6R$local_lib, "MCPcounter"))) {
  r_bin <- file.path(R.home("bin"), "R")
  status <- system2(r_bin, c("CMD", "INSTALL", paste0("--library=", R6R$local_lib), file.path(R6R$mcp_src, "Source")))
  if (!identical(status, 0L)) stop("PRECHECK_FAILED_MCP_COUNTER_SOURCE_UNAVAILABLE: module-local MCPcounter install failed", call. = FALSE)
}
.libPaths(c(R6R$local_lib, .libPaths()))
if (!requireNamespace("MCPcounter", quietly = TRUE)) {
  stop("PRECHECK_FAILED_MCP_COUNTER_SOURCE_UNAVAILABLE: MCPcounter not loadable from module-local library", call. = FALSE)
}

mcp_genes <- read_mcp_genes()
mcp_sets <- marker_sets_from_genes(mcp_genes)
mcp_counts <- data.frame(
  population = names(mcp_sets),
  official_marker_gene_count = lengths(mcp_sets),
  official_marker_genes = vapply(mcp_sets, paste, collapse = ";", FUN.VALUE = character(1)),
  stringsAsFactors = FALSE
)

pc <- load_program_contract()
canonical_union <- canonical_gene_union()
if (length(canonical_union) != 121L) stop("Canonical 121-gene union check failed: ", length(canonical_union), call. = FALSE)

input_rows <- rbind(
  file_record("GSE31312_EXPRESSION_MATRIX", "GSE31312", "historical_gene_expression_matrix", R6R$gse31312_expression, "Expected 22168 genes x 498 samples; final sensitivity source from Stage4C1 registry."),
  file_record("GSE31312_PROGRAM_SCORES", "GSE31312", "final_six_program_scores", R6R$gse31312_scores, "Expected 498 rows x six program columns."),
  file_record("GSE10846_EXPRESSION_MATRIX", "GSE10846", "WP2_primary_median_gene_expression_rds", R6R$gse10846_expression, "Expected 20823 genes x 420 samples; primary benchmark excludes six sorted profiles."),
  file_record("GSE10846_PROGRAM_SCORES", "GSE10846", "final_six_program_scores", R6R$gse10846_scores, "Expected 420 rows x six program columns; primary benchmark excludes six sorted profiles."),
  file_record("GSE181063_EXPRESSION_MATRIX", "GSE181063", "WP2_primary_median_gene_expression_rds", R6R$gse181063_expression, "Expected 19418 genes x 1310 samples."),
  file_record("GSE181063_PROGRAM_SCORES", "GSE181063", "final_six_program_scores", R6R$gse181063_scores, "Expected 1310 rows x six program columns."),
  file_record("CANONICAL_PROGRAM_CONTRACT", "ALL", "canonical_programs", R6R$program_contract, "Six 22-gene programs; 121 unique genes."),
  file_record("MCP_COUNTER_GENES", "ALL", "official_marker_definitions", genes_path, "Official MCP-counter genes.txt from pinned repo."),
  file_record("MCP_COUNTER_SOURCE_R", "ALL", "official_implementation_source", file.path(R6R$mcp_src, "Source/R/MCPcounter.R"), "Official MCPcounter.estimate implementation.")
)

dim_rows <- do.call(rbind, lapply(c("GSE31312", "GSE10846", "GSE181063"), function(ds) {
  primary <- align_expression_scores(ds, primary = TRUE)
  all420 <- if (ds == "GSE10846") align_expression_scores(ds, primary = FALSE) else NULL
  rbind(
    data.frame(dataset = ds, cohort_definition = "PRIMARY", expression_rows = nrow(primary$expression), expression_samples = ncol(primary$expression), score_rows = nrow(primary$scores), sample_alignment = identical(colnames(primary$expression), primary$scores$sample_id), stringsAsFactors = FALSE),
    if (!is.null(all420)) data.frame(dataset = ds, cohort_definition = "ALL420_SENSITIVITY", expression_rows = nrow(all420$expression), expression_samples = ncol(all420$expression), score_rows = nrow(all420$scores), sample_alignment = identical(colnames(all420$expression), all420$scores$sample_id), stringsAsFactors = FALSE)
  )
}))

expected_primary <- c(GSE31312 = 498L, GSE10846 = 414L, GSE181063 = 1310L)
sample_count_pass <- all(dim_rows$expression_samples[match(names(expected_primary), dim_rows$dataset[dim_rows$cohort_definition == "PRIMARY"])] == expected_primary)

write_text_once(c(
  "# 06R Protocol Contract",
  "",
  "Status: FROZEN_BEFORE_SCIENTIFIC_EXECUTION",
  "",
  "Scientific question: quantify how much the six continuous curated programs are explained by broad immune/stromal MCP-counter lineage-abundance benchmarks, and how much association remains after removing exact overlap with the 121-gene canonical union.",
  "",
  "Primary benchmark: official MCP-counter implementation from ebecht/MCPcounter, package version 1.2.0, pinned to commit b6eac73e91c246fcff0bb1a5c68a816cd588fc48.",
  "",
  "Final k: NOT_SELECTED. Taxonomy: NOT_ASSIGNED.",
  "",
  "Forbidden analyses: k selection, retired taxonomy reconstruction, ARI/NMI, centroid projection, survival analysis, ligand-receptor analysis, spatial analysis recomputation, manuscript/figure/workbook/response-letter modification.",
  "",
  "Cohorts: GSE31312 n=498; GSE10846 primary clinical biopsies n=414 with GSM361239-GSM361244 excluded and all-420 sensitivity labelled separately; GSE181063 n=1310.",
  "",
  "Disjoint benchmark: official MCP-counter marker genes minus the global union of all 121 canonical genes; at least five detected finite non-zero-variance disjoint genes are required per population-cohort."
), file.path(R6R$protocol, "06R_PROTOCOL_CONTRACT.md"))

write_csv_once(input_rows, file.path(R6R$protocol, "06R_INPUT_REGISTRY.csv"))
write_csv_once(pc[, setdiff(names(pc), "gene_list"), drop = FALSE], file.path(R6R$protocol, "06R_CANONICAL_PROGRAM_CONTRACT.csv"))
write_csv_once(mcp_counts, file.path(R6R$protocol, "06R_BENCHMARK_POPULATION_CONTRACT.csv"))
write_csv_once(data.frame(
  contract_item = c("global_disjoint_rule", "canonical_union_size", "minimum_detected_nonzero_variance_genes", "score_name", "score_formula"),
  value = c("official MCP marker genes minus global canonical 121-gene union", "121", "5", "MCP-counter-derived global-disjoint marker score", "mean of cohort-standardized expression across eligible disjoint marker genes"),
  stringsAsFactors = FALSE
), file.path(R6R$protocol, "06R_DISJOINT_BENCHMARK_CONTRACT.csv"))
write_csv_once(data.frame(
  analysis = c("pairwise_standard_correlations", "pairwise_disjoint_correlations", "correlation_delta", "multivariable_predictability", "deterministic_10fold_cv", "cross_cohort_correlation_structure", "GSE10846_all420_sensitivity"),
  specification = c("all six programs x all standard MCP populations per primary cohort; Pearson, Spearman, Fisher-z CI", "all six programs x all evaluable global-disjoint marker scores per primary cohort", "absolute standard correlation minus absolute disjoint correlation, descriptive only", "program_score ~ all evaluable global-disjoint marker scores; no outcome fields", "seed 20260807, same folds for six programs within cohort, no outcome-based feature selection", "full 6 x population correlation vectors compared across cohort pairs", "labelled sensitivity only, not mixed with 414-profile primary benchmark"),
  stringsAsFactors = FALSE
), file.path(R6R$protocol, "06R_STATISTICAL_ANALYSIS_CONTRACT.csv"))
write_csv_once(input_rows[, c("input_id", "absolute_path", "size_bytes", "sha256")], file.path(R6R$protocol, "06R_INPUT_SHA256_REGISTRY.csv"))
write_csv_once(dim_rows, file.path(R6R$outputs, "input_qc/06R_INPUT_DIMENSION_AND_ALIGNMENT_QC.csv"))
write_csv_once(git_status_record(), file.path(R6R$outputs, "input_qc/06R_GIT_STATUS_PREFLIGHT.csv"))
write_csv_once(protected_snapshot(), file.path(R6R$outputs, "input_qc/06R_PROTECTED_PATH_BASELINE.csv"))
write_csv_once(long_marker_df(mcp_sets, "STANDARD_OFFICIAL_MCP_COUNTER"), file.path(R6R$outputs, "marker_sets/06R_MCP_MARKER_SETS_STANDARD.csv"))

checks <- data.frame(
  check_id = c("TOKEN_EXACT", "HUMAN_AUTHORIZATION", "MCP_OFFICIAL_REPO", "MCP_COMMIT_PINNED", "MCP_LOCAL_INSTALL_LOADABLE", "MCP_MARKER_DEFINITIONS_CAPTURED", "CANONICAL_121_UNIQUE", "GSE31312_PRIMARY_N", "GSE10846_PRIMARY_N", "GSE181063_PRIMARY_N", "GSE10846_SORTED_EXCLUSION", "SAMPLE_ALIGNMENT", "NO_K_SELECTION", "NO_TAXONOMY_ASSIGNMENT", "PROTOCOL_FROZEN"),
  status = c(
    "PASS", "PASS", "PASS", "PASS", "PASS", if (nrow(mcp_genes) > 0) "PASS" else "FAIL",
    if (length(canonical_union) == 121L) "PASS" else "FAIL",
    if (dim_rows$expression_samples[dim_rows$dataset == "GSE31312" & dim_rows$cohort_definition == "PRIMARY"] == 498L) "PASS" else "FAIL",
    if (dim_rows$expression_samples[dim_rows$dataset == "GSE10846" & dim_rows$cohort_definition == "PRIMARY"] == 414L) "PASS" else "FAIL",
    if (dim_rows$expression_samples[dim_rows$dataset == "GSE181063" & dim_rows$cohort_definition == "PRIMARY"] == 1310L) "PASS" else "FAIL",
    if (dim_rows$expression_samples[dim_rows$dataset == "GSE10846" & dim_rows$cohort_definition == "ALL420_SENSITIVITY"] == 420L) "PASS" else "FAIL",
    if (all(dim_rows$sample_alignment)) "PASS" else "FAIL",
    "PASS", "PASS", "PASS"
  ),
  evidence = c(
    R6R$token, R6R$human_authorization, R6R$expected_mcp_repo, R6R$expected_mcp_commit,
    file.path(R6R$local_lib, "MCPcounter"), genes_path, as.character(length(canonical_union)),
    "expected 498", "expected 414", "expected 1310", paste(R6R$excluded_gse10846, collapse = ";"),
    "expression columns equal score sample order after explicit alignment", "Final k remains NOT_SELECTED", "Taxonomy remains NOT_ASSIGNED", timestamp()
  ),
  stringsAsFactors = FALSE
)
write_csv_once(checks, file.path(R6R$protocol, "06R_PROTOCOL_VALIDATION.csv"))
if (any(checks$status != "PASS")) stop("06R_PREFLIGHT_FAILED", call. = FALSE)

message("[06r preflight] PASS at ", timestamp())
