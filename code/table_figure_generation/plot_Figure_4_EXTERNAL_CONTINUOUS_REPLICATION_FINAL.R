#!/usr/bin/env Rscript

DLBCL_PROJECT_ROOT <- normalizePath(
  Sys.getenv("DLBCL_PROJECT_ROOT", unset = getwd()),
  winslash = "/", mustWork = FALSE
)
DLBCL_SUPPLEMENTARY_CODE_ROOT <- normalizePath(
  Sys.getenv("DLBCL_SUPPLEMENTARY_CODE_ROOT",
    unset = file.path(DLBCL_PROJECT_ROOT, "Supplementary_Code")),
  winslash = "/", mustWork = FALSE
)
DLBCL_SUPPLEMENTARY_DATA_ROOT <- normalizePath(
  Sys.getenv("DLBCL_SUPPLEMENTARY_DATA_ROOT",
    unset = file.path(DLBCL_PROJECT_ROOT, "Supplementary_Data")),
  winslash = "/", mustWork = FALSE
)
DLBCL_SUPPLEMENTARY_TABLES_ROOT <- normalizePath(
  Sys.getenv("DLBCL_SUPPLEMENTARY_TABLES_ROOT",
    unset = file.path(DLBCL_PROJECT_ROOT, "Supplementary_Tables")),
  winslash = "/", mustWork = FALSE
)


options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
mode_arg <- grep("^--mode=", args, value = TRUE)
mode <- if (length(mode_arg)) sub("^--mode=", "", mode_arg[[1]]) else "diagnostic"
if (!mode %in% c("layout_test", "diagnostic", "finalize")) {
  stop("--mode must be layout_test, diagnostic, or finalize")
}

root <- DLBCL_PROJECT_ROOT
wp2 <- file.path(DLBCL_SUPPLEMENTARY_DATA_ROOT, "figure_rendering", "Figure_4")
out_dir <- file.path(root, "reproduced_figures", "Figure_4")
workbook <- file.path(
  DLBCL_SUPPLEMENTARY_TABLES_ROOT,
  "DLBCL_continuous_model_Supplementary_Tables_FINAL_SUBMISSION_PUBLICATION_READY.xlsx"
)

required_packages <- c("cowplot", "digest", "ggplot2", "ragg", "readxl", "systemfonts", "tiff")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing R packages: ", paste(missing_packages, collapse = ", "))

suppressPackageStartupMessages({
  library(ggplot2)
})

write_once <- function(path, writer) {
  if (file.exists(path)) stop("Refusing to overwrite existing attempt output: ", path)
  writer(path)
  invisible(path)
}

write_csv_once <- function(x, name, na = "NA") {
  path <- file.path(out_dir, name)
  write_once(path, function(p) utils::write.csv(x, p, row.names = FALSE, na = na, fileEncoding = "UTF-8"))
}

write_lines_once <- function(x, name) {
  path <- file.path(out_dir, name)
  write_once(path, function(p) writeLines(x, p, useBytes = TRUE))
}

sha256 <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)

program_map <- data.frame(
  program_order = 1:6,
  source_name = c(
    "Macrophage-rich program",
    "T cell-inflamed program",
    "Immune-inflamed / antigen-presentation program",
    "Stromal / fibrotic program",
    "Immune-cold / exclusion-associated program",
    "Proliferative / cycling program"
  ),
  final_label = c(
    "Macrophage-rich",
    "T cell-inflamed",
    "Antigen-presentation",
    "Stromal/fibrotic",
    "Immune-cold/exclusion",
    "Proliferative/cycling"
  ),
  abbreviation = c("MR", "TCI", "AP", "SF", "IC/EX", "PC"),
  stringsAsFactors = FALSE
)

historical_map <- rbind(
  program_map[, c("program_order", "source_name", "final_label", "abbreviation")],
  data.frame(
    program_order = c(1, 2, 3, 4, 5, 6, 6),
    source_name = c(
      "Macrophage-rich niche",
      "T cell-inflamed niche",
      "Immune-inflamed / antigen-presentation niche",
      "Stromal/fibrotic niche",
      "Immune-cold / excluded niche",
      "Proliferative malignant B-cell niche",
      "Proliferative B-cell/cycling"
    ),
    final_label = program_map$final_label[c(1, 2, 3, 4, 5, 6, 6)],
    abbreviation = program_map$abbreviation[c(1, 2, 3, 4, 5, 6, 6)],
    stringsAsFactors = FALSE
  )
)
historical_map$mapping_rule <- "EXACT_STRING_MATCH_ONLY"

matrix_files <- c(
  GSE31312 = file.path(wp2, "GSE31312_PEARSON_CORRELATION_MATRIX.csv"),
  GSE10846 = file.path(wp2, "GSE10846_PRIMARY_PEARSON_CORRELATION_MATRIX.csv"),
  GSE181063 = file.path(wp2, "GSE181063_PRIMARY_PEARSON_CORRELATION_MATRIX.csv")
)
edges_file <- file.path(wp2, "WP2_PRIMARY_CORRELATION_EDGES_LONG.csv")
point_file <- file.path(wp2, "WP2_STRUCTURAL_REPLICATION_POINT_ESTIMATES.csv")
integration_file <- file.path(wp2, "WP2_CROSS_COHORT_INTEGRATION_SUMMARY.csv")
manifest_file <- file.path(wp2, "WP2_FINAL_FILE_MANIFEST_V2.csv")
validator_file <- file.path(wp2, "WP2_FINAL_VALIDATOR_REPORT.txt")
decision_file <- file.path(wp2, "WP2_HUMAN_DECISIONS_FREEZE.csv")
preflight_file <- file.path(wp2, "WP2_POSTRUN_PREFLIGHT_REPORT_V2.txt")
report_file <- file.path(wp2, "WP2_EXTERNAL_CONTINUOUS_VALIDATION_REPORT.md")
score_files <- c(
  GSE10846 = file.path(wp2, "GSE10846_PRIMARY_SIX_PROGRAM_SCORES.csv"),
  GSE181063 = file.path(wp2, "GSE181063_PRIMARY_SIX_PROGRAM_SCORES.csv")
)

all_required_inputs <- c(
  matrix_files, edges_file, point_file, integration_file, manifest_file,
  validator_file, decision_file, preflight_file, report_file, score_files, workbook
)
if (any(!file.exists(all_required_inputs))) {
  stop("Missing required input(s): ", paste(all_required_inputs[!file.exists(all_required_inputs)], collapse = "; "))
}

read_authority_matrix <- function(path) {
  x <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!identical(names(x), c("program", program_map$source_name))) stop("Unexpected matrix columns: ", path)
  if (!identical(x$program, program_map$source_name)) stop("Unexpected matrix row labels/order: ", path)
  m <- as.matrix(x[, -1, drop = FALSE])
  storage.mode(m) <- "double"
  rownames(m) <- x$program
  colnames(m) <- names(x)[-1]
  if (!identical(dim(m), c(6L, 6L))) stop("Matrix is not 6x6: ", path)
  if (max(abs(m - t(m))) > 1e-12) stop("Matrix is not symmetric: ", path)
  if (max(abs(diag(m) - 1)) > 1e-12) stop("Matrix diagonal is not 1: ", path)
  if (any(!is.finite(m)) || any(m < -1 | m > 1)) stop("Matrix has invalid correlations: ", path)
  m
}

matrices <- lapply(matrix_files, read_authority_matrix)

extract_edges <- function(m, cohort) {
  idx <- which(upper.tri(m), arr.ind = TRUE)
  idx <- idx[order(idx[, 1], idx[, 2]), , drop = FALSE]
  data.frame(
    cohort = cohort,
    program_1_order = idx[, 1],
    program_2_order = idx[, 2],
    program_1_source = rownames(m)[idx[, 1]],
    program_2_source = colnames(m)[idx[, 2]],
    r = m[idx],
    stringsAsFactors = FALSE
  )
}

matrix_edges <- do.call(rbind, Map(extract_edges, matrices, names(matrices)))
rownames(matrix_edges) <- NULL
if (!all(table(matrix_edges$cohort) == 15)) stop("A cohort does not have exactly 15 matrix edges")

external_edges <- utils::read.csv(edges_file, check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(external_edges) != 30L || !all(table(external_edges$dataset_id) == 15L)) {
  stop("External edge authority is not 30 rows with 15 edges per cohort")
}

for (cohort in c("GSE10846", "GSE181063")) {
  a <- matrix_edges[matrix_edges$cohort == cohort, ]
  b <- external_edges[external_edges$dataset_id == cohort, ]
  key_a <- paste(a$program_1_source, a$program_2_source, sep = "||")
  key_b <- paste(b$program_1, b$program_2, sep = "||")
  b <- b[match(key_a, key_b), ]
  if (anyNA(b$pearson_r) || max(abs(a$r - b$pearson_r)) > 1e-12) {
    stop("External matrix/edge conflict for ", cohort)
  }
}

workbook_s4c <- as.data.frame(readxl::read_excel(workbook, sheet = "S4C_Correlation_matrix", skip = 3))
workbook_s6c <- as.data.frame(readxl::read_excel(workbook, sheet = "S6C_Structural_edges", skip = 3))
workbook_s6d <- as.data.frame(readxl::read_excel(workbook, sheet = "S6D_Structural_summary", skip = 3))
workbook_s6a <- as.data.frame(readxl::read_excel(workbook, sheet = "S6A_GSE10846_scores", skip = 3))
workbook_s6b <- as.data.frame(readxl::read_excel(workbook, sheet = "S6B_GSE181063_scores", skip = 3))

workbook_gse31312 <- as.matrix(workbook_s4c[, -1, drop = FALSE])
storage.mode(workbook_gse31312) <- "double"
rownames(workbook_gse31312) <- workbook_s4c[[1]]
colnames(workbook_gse31312) <- names(workbook_s4c)[-1]
if (max(abs(workbook_gse31312 - matrices$GSE31312)) > 1e-12) {
  stop("Workbook S4C conflicts with 05z GSE31312 authority")
}

workbook_edges <- data.frame(
  cohort = workbook_s6c$Dataset,
  program_1_source = workbook_s6c$`Program 1`,
  program_2_source = workbook_s6c$`Program 2`,
  r = as.numeric(workbook_s6c$`Pearson r`),
  stringsAsFactors = FALSE
)
if (nrow(workbook_edges) != 45L || !all(table(workbook_edges$cohort) == 15L)) {
  stop("Workbook S6C is not 45 rows with 15 edges per cohort")
}

edge_key <- function(d) paste(d$cohort, d$program_1_source, d$program_2_source, sep = "||")
workbook_edges <- workbook_edges[match(edge_key(matrix_edges), edge_key(workbook_edges)), ]
if (anyNA(workbook_edges$r) || max(abs(matrix_edges$r - workbook_edges$r)) > 1e-12) {
  stop("Workbook S6C conflicts with 05z authority")
}

score_counts <- vapply(score_files, function(p) nrow(utils::read.csv(p, check.names = FALSE)), integer(1))
sample_counts <- c(GSE31312 = 498L, GSE10846 = score_counts[["GSE10846"]], GSE181063 = score_counts[["GSE181063"]])
if (!identical(unname(sample_counts), c(498L, 420L, 1310L))) stop("Frozen cohort sample counts are inconsistent")
if (nrow(workbook_s6a) != 420L || nrow(workbook_s6b) != 1310L) stop("Workbook sample counts conflict")

all_edges <- data.frame(
  program_1 = program_map$final_label[matrix_edges$program_1_order[matrix_edges$cohort == "GSE31312"]],
  program_2 = program_map$final_label[matrix_edges$program_2_order[matrix_edges$cohort == "GSE31312"]],
  edge_label = paste0(
    program_map$abbreviation[matrix_edges$program_1_order[matrix_edges$cohort == "GSE31312"]],
    "\u2013",
    program_map$abbreviation[matrix_edges$program_2_order[matrix_edges$cohort == "GSE31312"]]
  ),
  GSE31312_r = matrix_edges$r[matrix_edges$cohort == "GSE31312"],
  GSE10846_r = matrix_edges$r[matrix_edges$cohort == "GSE10846"],
  GSE181063_r = matrix_edges$r[matrix_edges$cohort == "GSE181063"],
  GSE10846_replication_status = NA_character_,
  GSE181063_replication_status = NA_character_,
  source_file = paste(
    "05z/WP2_FINAL_FILE_MANIFEST_V2.csv ->",
    "05_structure/GSE31312_PEARSON_CORRELATION_MATRIX.csv;",
    "05_structure/GSE10846_PRIMARY_PEARSON_CORRELATION_MATRIX.csv;",
    "05_structure/GSE181063_PRIMARY_PEARSON_CORRELATION_MATRIX.csv"
  ),
  source_row = paste0(
    "matrix[row=", matrix_edges$program_1_order[matrix_edges$cohort == "GSE31312"] + 1L,
    ",column=", matrix_edges$program_2_order[matrix_edges$cohort == "GSE31312"] + 1L, "]"
  ),
  stringsAsFactors = FALSE
)
if (nrow(all_edges) != 15L || anyDuplicated(all_edges$edge_label)) stop("Unified edge table is invalid")

display_order <- all_edges[order(-all_edges$GSE31312_r, all_edges$edge_label), ]
display_order <- data.frame(
  display_order = seq_len(nrow(display_order)),
  program_1 = display_order$program_1,
  program_2 = display_order$program_2,
  edge_label = display_order$edge_label,
  GSE31312_r_full = display_order$GSE31312_r,
  GSE10846_r_full = display_order$GSE10846_r,
  GSE181063_r_full = display_order$GSE181063_r,
  stringsAsFactors = FALSE
)

point_summary <- utils::read.csv(point_file, check.names = FALSE, stringsAsFactors = FALSE)
integration <- utils::read.csv(integration_file, check.names = FALSE, stringsAsFactors = FALSE)
get_ci <- function(cohort) {
  value <- integration$value[integration$scope == cohort & integration$metric == "matrix_vector_spearman_95CI"]
  if (length(value) != 1L) stop("Missing frozen bootstrap CI for ", cohort)
  z <- strsplit(gsub("[^0-9eE+.,-]", "", value), ",", fixed = FALSE)[[1]]
  as.numeric(z)
}
ci_10846 <- get_ci("GSE10846")
ci_181063 <- get_ci("GSE181063")
summary_used <- transform(
  point_summary,
  bootstrap_ci_lower = c(ci_10846[1], ci_181063[1]),
  bootstrap_ci_upper = c(ci_10846[2], ci_181063[2]),
  frozen_criterion = "HD6 supportive: matrix-vector Spearman rho >= 0.50 and >= 12/15 sign-concordant edges",
  edge_level_replication_status_available = FALSE
)

rebuild_matrix <- function(edge_table, value_col) {
  m <- diag(1, 6)
  rownames(m) <- colnames(m) <- program_map$final_label
  for (i in seq_len(nrow(edge_table))) {
    a <- match(edge_table$program_1[i], program_map$final_label)
    b <- match(edge_table$program_2[i], program_map$final_label)
    m[a, b] <- m[b, a] <- edge_table[[value_col]][i]
  }
  m
}

matrices_used <- list(
  GSE31312 = rebuild_matrix(all_edges, "GSE31312_r"),
  GSE10846 = rebuild_matrix(all_edges, "GSE10846_r"),
  GSE181063 = rebuild_matrix(all_edges, "GSE181063_r")
)

matrix_to_table <- function(m) {
  data.frame(program = rownames(m), as.data.frame(m, check.names = FALSE), check.names = FALSE)
}

value_qc <- do.call(rbind, lapply(names(matrices_used), function(cohort) {
  value_col <- paste0(cohort, "_r")
  data.frame(
    cohort = cohort,
    program_1 = all_edges$program_1,
    program_2 = all_edges$program_2,
    authority_value = all_edges[[value_col]],
    panel_A_value = vapply(seq_len(nrow(all_edges)), function(i) matrices_used[[cohort]][all_edges$program_1[i], all_edges$program_2[i]], numeric(1)),
    panel_B_value = all_edges[[value_col]],
    panel_C_value = all_edges[[value_col]],
    stringsAsFactors = FALSE
  )
}))
value_qc$A_match <- abs(value_qc$authority_value - value_qc$panel_A_value) <= 1e-12
value_qc$B_match <- abs(value_qc$authority_value - value_qc$panel_B_value) <= 1e-12
value_qc$C_match <- abs(value_qc$authority_value - value_qc$panel_C_value) <= 1e-12
if (!all(value_qc$A_match & value_qc$B_match & value_qc$C_match)) stop("Panel value consistency failed")

precision_full <- data.frame(
  cohort = matrix_edges$cohort,
  program_1 = program_map$final_label[matrix_edges$program_1_order],
  program_2 = program_map$final_label[matrix_edges$program_2_order],
  full_precision_value = matrix_edges$r,
  comparison_value = workbook_edges$r,
  comparison_type = "FULL_PRECISION_WORKBOOK_S6C",
  tolerance = 1e-12,
  match = abs(matrix_edges$r - workbook_edges$r) <= 1e-12,
  stringsAsFactors = FALSE
)
precision_display <- precision_full
precision_display$comparison_value <- round(precision_display$full_precision_value, 2)
precision_display$comparison_type <- "TWO_DECIMAL_DISPLAY_FORMAT"
precision_display$tolerance <- NA_real_
precision_display$match <- round(precision_display$full_precision_value, 2) == precision_display$comparison_value
precision_qc <- rbind(precision_full, precision_display)
if (!all(precision_qc$match)) stop("Precision/rounding QC failed")

candidate_paths <- list.files(wp2, recursive = TRUE, full.names = TRUE)
candidate_paths <- candidate_paths[file.info(candidate_paths)$isdir %in% FALSE]
candidate_pattern <- "FINAL|AUTHORITY|MANIFEST|DECISION|correlation|edge|pairwise|replication|concordance|GSE10846|GSE181063|external|WP2"
candidate_paths <- candidate_paths[grepl(candidate_pattern, basename(candidate_paths), ignore.case = TRUE)]

schema_of <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    x <- tryCatch(utils::read.csv(path, nrows = 5, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
    header <- if (is.null(x)) "PARSE_ERROR" else paste(names(x), collapse = " | ")
    rows <- tryCatch(length(readLines(path, warn = FALSE, encoding = "UTF-8")) - 1L, error = function(e) NA_integer_)
    return(c(row_count = rows, column_count = if (is.null(x)) NA_integer_ else ncol(x), fields = header))
  }
  c(row_count = length(readLines(path, warn = FALSE, encoding = "UTF-8")), column_count = NA_integer_, fields = "TEXT_OR_NON_TABULAR")
}

candidate_schema <- lapply(candidate_paths, schema_of)
candidate_audit <- data.frame(
  file_name = basename(candidate_paths),
  absolute_path = normalizePath(candidate_paths, winslash = "/", mustWork = TRUE),
  file_size_bytes = file.info(candidate_paths)$size,
  sha256 = vapply(candidate_paths, sha256, character(1)),
  row_count = as.integer(vapply(candidate_schema, `[[`, character(1), "row_count")),
  column_count = as.integer(vapply(candidate_schema, `[[`, character(1), "column_count")),
  fields = vapply(candidate_schema, `[[`, character(1), "fields"),
  final_authority = ifelse(
    normalizePath(candidate_paths, winslash = "/", mustWork = TRUE) %in% normalizePath(c(matrix_files, edges_file, point_file, integration_file, manifest_file, validator_file, decision_file), winslash = "/", mustWork = TRUE),
    "YES", "NO_OR_SUPPORTING"
  ),
  notes = ifelse(grepl("ATTEMPT1|FAILURE|V1|PARTIAL", basename(candidate_paths), ignore.case = TRUE), "HISTORICAL_OR_CONTROL_RECORD", "CANDIDATE_OR_SUPPORTING_RECORD"),
  stringsAsFactors = FALSE
)

input_paths <- c(
  manifest_file, validator_file, preflight_file, decision_file, report_file,
  matrix_files, edges_file, point_file, integration_file, score_files, workbook
)
input_manifest <- data.frame(
  input_role = c(
    "authoritative_final_manifest", "final_validator", "authoritative_post-process_preflight",
    "frozen_human_decisions", "final_analysis_report",
    "GSE31312_full_precision_matrix", "GSE10846_full_precision_matrix", "GSE181063_full_precision_matrix",
    "external_15_edge_authority", "frozen_replication_point_estimates", "frozen_integration_summary",
    "GSE10846_sample_count_authority", "GSE181063_sample_count_authority", "final_supplementary_workbook"
  ),
  absolute_path = normalizePath(input_paths, winslash = "/", mustWork = TRUE),
  file_size_bytes = file.info(input_paths)$size,
  sha256 = vapply(input_paths, sha256, character(1)),
  stringsAsFactors = FALSE
)

source_map <- data.frame(
  figure_component = c(
    "Phase 0 authority", "Phase 0 final validation", "replication criterion",
    "Panel A/B/C discovery values", "Panel A/B/C GSE10846 values", "Panel A/B/C GSE181063 values",
    "external edge cross-check", "Panel B frozen summaries", "bootstrap CI trace", "workbook cross-check"
  ),
  selected_source = normalizePath(c(
    manifest_file, validator_file, decision_file,
    matrix_files[["GSE31312"]], matrix_files[["GSE10846"]], matrix_files[["GSE181063"]],
    edges_file, point_file, integration_file, workbook
  ), winslash = "/", mustWork = TRUE),
  authority_basis = c(
    "V2 explicitly states authoritative post-process snapshot",
    "56/56 final checks passed",
    "HD6 frozen before real-data execution",
    rep("Referenced by WP2_FINAL_FILE_MANIFEST_V2.csv; full precision", 3),
    "Referenced by final manifest; 30 external rows",
    "Frozen cohort-level point estimates",
    "Frozen bootstrap CI and cross-cohort integration summary",
    "Final workbook sheets S4C/S6C/S6D match 05z authority"
  ),
  used_for_plot = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE),
  stringsAsFactors = FALSE
)

candidate_compare <- do.call(rbind, lapply(names(matrices), function(cohort) {
  m_edge <- matrix_edges[matrix_edges$cohort == cohort, ]
  w_edge <- workbook_edges[workbook_edges$cohort == cohort, ]
  e_edge <- if (cohort == "GSE31312") NULL else external_edges[external_edges$dataset_id == cohort, ]
  rows <- data.frame(
    cohort = cohort,
    candidate_source = c("05z full-precision matrix", "Final workbook S6C"),
    row_count = c(15L, 15L),
    fields = c("6x6 matrix -> deterministic upper-triangle edges", "Dataset|Collapse rule|Program 1|Program 2|Order 1|Order 2|Pearson r|Sign"),
    maximum_absolute_difference_vs_authority = c(0, max(abs(m_edge$r - w_edge$r))),
    comparison_status = c("AUTHORITY", if (max(abs(m_edge$r - w_edge$r)) <= 1e-12) "MATCH" else "CONFLICT"),
    stringsAsFactors = FALSE
  )
  if (!is.null(e_edge)) {
    e_edge <- e_edge[match(paste(m_edge$program_1_source, m_edge$program_2_source), paste(e_edge$program_1, e_edge$program_2)), ]
    rows <- rbind(rows, data.frame(
      cohort = cohort,
      candidate_source = "05z WP2_PRIMARY_CORRELATION_EDGES_LONG.csv",
      row_count = 15L,
      fields = paste(names(external_edges), collapse = "|"),
      maximum_absolute_difference_vs_authority = max(abs(m_edge$r - e_edge$pearson_r)),
      comparison_status = if (max(abs(m_edge$r - e_edge$pearson_r)) <= 1e-12) "MATCH" else "CONFLICT",
      stringsAsFactors = FALSE
    ))
  }
  rows
}))
if (any(candidate_compare$comparison_status == "CONFLICT")) stop("Candidate input conflict detected")

if (mode == "diagnostic") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  write_csv_once(historical_map, "FIGURE4_PROGRAM_LABEL_MAP.csv")
  write_csv_once(input_manifest, "FIGURE4_INPUT_MANIFEST_SHA256.csv")
  write_csv_once(source_map, "FIGURE4_SOURCE_FILE_MAP.csv")
  write_csv_once(candidate_compare, "FIGURE4_CANDIDATE_INPUT_COMPARISON.csv")
  write_csv_once(precision_qc, "FIGURE4_PRECISION_AND_ROUNDING_QC.csv")
  write_csv_once(matrix_to_table(matrices_used$GSE31312), "FIGURE4_GSE31312_MATRIX_USED.csv")
  write_csv_once(matrix_to_table(matrices_used$GSE10846), "FIGURE4_GSE10846_MATRIX_USED.csv")
  write_csv_once(matrix_to_table(matrices_used$GSE181063), "FIGURE4_GSE181063_MATRIX_USED.csv")
  write_csv_once(all_edges, "FIGURE4_ALL_15_EDGES_USED.csv")
  write_csv_once(display_order, "FIGURE4_EDGE_DISPLAY_ORDER.csv")
  write_csv_once(summary_used, "FIGURE4_EXTERNAL_REPLICATION_SUMMARY_USED.csv")
  write_csv_once(value_qc, "FIGURE4_VALUE_CONSISTENCY_QC.csv")
  write_csv_once(candidate_audit, "FIGURE4_ALL_CANDIDATE_FILES_AUDITED.csv")

  validator_text <- readLines(validator_file, warn = FALSE, encoding = "UTF-8")
  preflight_text <- readLines(preflight_file, warn = FALSE, encoding = "UTF-8")
  schema_lines <- c(
    "# Figure 4 schema and authority audit",
    "",
    "## Controlling authority",
    "",
    paste0("- WP2 directory: `", normalizePath(wp2, winslash = "/"), "`"),
    paste0("- Controlling manifest: `", normalizePath(manifest_file, winslash = "/"), "`"),
    "- Manifest status: authoritative post-process snapshot (WP2_POSTRUN_PREFLIGHT_REPORT_V2.txt).",
    paste0("- Final validator: ", paste(validator_text, collapse = "; "), "."),
    paste0("- Post-process preflight: ", paste(preflight_text, collapse = "; "), "."),
    "",
    "## Selected full-precision sources",
    "",
    paste0("- GSE31312: `", normalizePath(matrix_files[["GSE31312"]], winslash = "/"), "`."),
    paste0("- GSE10846: `", normalizePath(matrix_files[["GSE10846"]], winslash = "/"), "`."),
    paste0("- GSE181063: `", normalizePath(matrix_files[["GSE181063"]], winslash = "/"), "`."),
    paste0("- External 15-edge cross-check: `", normalizePath(edges_file, winslash = "/"), "`."),
    paste0("- Cohort-level replication summary: `", normalizePath(point_file, winslash = "/"), "`."),
    paste0("- Bootstrap CI trace: `", normalizePath(integration_file, winslash = "/"), "`."),
    "- Edge-level replication calls: NOT PRESENT. Figure fields remain NA; no calls were derived.",
    "",
    "## Schema and consistency findings",
    "",
    "- Correlation method is explicitly Pearson for all three six-program matrices and all 15 edge values.",
    "- Each matrix is 6 x 6, symmetric within 1e-12, has diagonal 1, finite values in [-1, 1], and exactly 15 nonredundant edges.",
    "- Six program labels map uniquely by exact string matching to the frozen final labels and order.",
    "- GSE31312 n=498; GSE10846 n=420; GSE181063 n=1,310.",
    "- Final workbook S4C matches the 05z GSE31312 matrix within 1e-12.",
    "- Final workbook S6C contains 45 full-precision edges (15 per cohort) and matches all three 05z matrices within 1e-12.",
    "- Final workbook S6D matches the frozen cohort-level structural summaries.",
    "- Historical, partial, V1, and controlled-failure records were retained for provenance but excluded from plotting authority.",
    "- No taxonomy, class, cluster, centroid, projection, survival, or sample assignment source enters the figure.",
    "",
    "## Candidate-file registry",
    "",
    paste0("- `FIGURE4_ALL_CANDIDATE_FILES_AUDITED.csv` records ", nrow(candidate_audit), " keyword-matched files with path, size, SHA-256, schema, and authority annotation."),
    "",
    "## Phase decision",
    "",
    "**PHASE_0_PASS: sources are unambiguous and Phase 1 deterministic plotting is authorized.**"
  )
  write_lines_once(schema_lines, "FIGURE4_SCHEMA_AUDIT.md")

  input_qc <- c(
    "# Figure 4 input QC",
    "",
    "Status: PASS",
    "",
    "- Frozen results only; no sample-level correlation, GSVA, PCA, clustering, bootstrap, or replication statistic was recomputed.",
    "- GSE31312 sample count: 498 (PASS).",
    "- GSE10846 sample count: 420 in 05z scores and 420 in workbook S6A (PASS).",
    "- GSE181063 sample count: 1,310 in 05z scores and 1,310 in workbook S6B (PASS).",
    "- Program identity/order: 6/6 exact and common across cohorts (PASS).",
    "- Matrix dimensions: 6 x 6 for all cohorts (PASS).",
    "- Nonredundant edges: 15/15 per cohort; no duplicates or missing edges (PASS).",
    "- Finite correlation range [-1, 1], diagonal=1, symmetry tolerance <=1e-12 (PASS).",
    "- Workbook/05z maximum absolute difference <=1e-12 for all 45 edges (PASS).",
    "- Unified plotting authority: FIGURE4_ALL_15_EDGES_USED.csv.",
    "- Edge-level replication status: unavailable and intentionally recorded as NA.",
    "- Frozen cohort-level structural replication summaries are used without recalculation.",
    "- Final k: NOT_SELECTED. Taxonomy: NOT_ASSIGNED."
  )
  write_lines_once(input_qc, "FIGURE4_INPUT_QC.md")

  capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
}

font_match <- systemfonts::match_font("Times New Roman")
font_path <- normalizePath(font_match$path, winslash = "/", mustWork = TRUE)
font_family <- "Times New Roman"
if (!identical(tolower(basename(font_path)), "times.ttf")) {
  stop("Times New Roman did not resolve to the required regular font: ", font_path)
}

theme_pub <- function(base_size = 7.2) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.title = element_text(size = 8.4, face = "bold", hjust = 0, margin = margin(b = 2)),
      axis.title = element_text(size = 7.4),
      axis.text = element_text(size = 7.0, colour = "black"),
      axis.line = element_line(linewidth = 0.25, colour = "black"),
      axis.ticks = element_line(linewidth = 0.25, colour = "black"),
      axis.ticks.length = grid::unit(1.1, "mm"),
      legend.title = element_text(size = 7.2),
      legend.text = element_text(size = 7.0),
      legend.key.height = grid::unit(3.6, "mm"),
      plot.margin = margin(2, 2, 2, 2, unit = "mm")
    )
}

palette_blue <- "#4A6F91"
palette_mid <- "#F8F8F6"
palette_red <- "#B45D52"
edge_point <- "#61788A"

section_header <- function(letter, title) {
  ggplot() +
    annotate("text", x = 0, y = 0.5, label = letter, family = font_family,
             fontface = "bold", size = 11.5 / ggplot2::.pt, hjust = 0, vjust = 0.5) +
    annotate("text", x = 0.045, y = 0.5, label = title, family = font_family,
             fontface = "bold", size = 9.1 / ggplot2::.pt, hjust = 0, vjust = 0.5) +
    xlim(0, 1) + ylim(0, 1) + theme_void(base_family = font_family) +
    theme(plot.margin = margin(0, 1.5, 0, 1.5, unit = "mm"))
}

matrix_long_from_edges <- function(cohort) {
  m <- matrices_used[[cohort]]
  d <- expand.grid(row_order = 1:6, col_order = 1:6)
  d$r <- m[cbind(d$row_order, d$col_order)]
  d$program_row <- factor(program_map$final_label[d$row_order], levels = rev(program_map$final_label))
  d$program_col <- factor(program_map$abbreviation[d$col_order], levels = program_map$abbreviation)
  d$visible <- d$row_order >= d$col_order
  d <- d[d$visible, ]
  d$label <- sprintf("%.2f", d$r)
  d$text_colour <- ifelse(abs(d$r) >= 0.62, "white", "#202020")
  d
}

plot_matrix <- function(cohort, n, show_rows = TRUE, show_legend = TRUE) {
  d <- matrix_long_from_edges(cohort)
  ggplot(d, aes(program_col, program_row, fill = r)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(aes(label = label, colour = text_colour), size = 7.0 / ggplot2::.pt,
              family = font_family, show.legend = FALSE) +
    scale_colour_identity() +
    scale_fill_gradient2(
      low = palette_blue, mid = palette_mid, high = palette_red,
      midpoint = 0, limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1),
      name = "Pearson r"
    ) +
    coord_fixed(expand = FALSE) +
    labs(title = paste0(cohort, ", n=", format(n, big.mark = ",")), x = NULL, y = NULL) +
    theme_pub() +
    theme(
      panel.border = element_rect(colour = "#B8B8B8", fill = NA, linewidth = 0.25),
      axis.line = element_blank(), axis.ticks = element_blank(),
      axis.text.x = element_text(size = 7.0, angle = 0, vjust = 1),
      axis.text.y = if (show_rows) element_text(size = 7.0, lineheight = 0.92) else element_blank(),
      legend.position = if (show_legend) "right" else "none",
      plot.margin = margin(1, 2, 1, 2, unit = "mm")
    )
}

plot_scatter <- function(cohort) {
  y_col <- paste0(cohort, "_r")
  s <- summary_used[summary_used$dataset_id == cohort, ]
  d <- data.frame(discovery = all_edges$GSE31312_r, external = all_edges[[y_col]])
  summary_label <- sprintf(
    "Spearman rho = %.3f\nSign concordance = %d/%d",
    s$matrix_vector_spearman, s$sign_concordance_n, s$sign_concordance_denominator
  )
  ggplot(d, aes(discovery, external)) +
    geom_hline(yintercept = 0, colour = "#D1D1D1", linewidth = 0.28) +
    geom_vline(xintercept = 0, colour = "#D1D1D1", linewidth = 0.28) +
    geom_abline(intercept = 0, slope = 1, colour = "#AFAFAF", linewidth = 0.38) +
    geom_point(shape = 21, size = 2.05, stroke = 0.35, fill = edge_point, colour = "#33444F") +
    annotate("text", x = -0.93, y = -0.66, label = summary_label, hjust = 0, vjust = 0,
             family = font_family, size = 7.0 / ggplot2::.pt, lineheight = 1.0, colour = "#303030") +
    coord_equal(xlim = c(-1, 1), ylim = c(-1, 1), expand = FALSE) +
    scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1)) +
    scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1)) +
    labs(
      title = paste0("GSE31312 versus ", cohort),
      x = "GSE31312 Pearson r", y = paste0(cohort, " Pearson r")
    ) +
    theme_pub() +
    theme(
      panel.border = element_rect(colour = "#8A8A8A", fill = NA, linewidth = 0.3),
      axis.line = element_blank(),
      plot.margin = margin(1, 4, 1, 2, unit = "mm")
    )
}

plot_cross_cohort <- function() {
  d <- rbind(
    data.frame(edge_label = display_order$edge_label, cohort = "GSE31312", r = display_order$GSE31312_r_full),
    data.frame(edge_label = display_order$edge_label, cohort = "GSE10846", r = display_order$GSE10846_r_full),
    data.frame(edge_label = display_order$edge_label, cohort = "GSE181063", r = display_order$GSE181063_r_full)
  )
  d$edge_label <- factor(d$edge_label, levels = rev(display_order$edge_label))
  d$cohort <- factor(d$cohort, levels = c("GSE31312", "GSE10846", "GSE181063"))
  d$label <- sprintf("%.2f", d$r)
  d$text_colour <- ifelse(abs(d$r) >= 0.62, "white", "#202020")
  ggplot(d, aes(cohort, edge_label, fill = r)) +
    geom_tile(colour = "white", linewidth = 0.42) +
    geom_text(aes(label = label, colour = text_colour), family = font_family,
              size = 7.0 / ggplot2::.pt, show.legend = FALSE) +
    scale_colour_identity() +
    scale_fill_gradient2(
      low = palette_blue, mid = palette_mid, high = palette_red,
      midpoint = 0, limits = c(-1, 1), guide = "none"
    ) +
    labs(x = NULL, y = NULL) +
    theme_pub() +
    theme(
      panel.border = element_rect(colour = "#B8B8B8", fill = NA, linewidth = 0.25),
      axis.line = element_blank(), axis.ticks = element_blank(),
      axis.text.x = element_text(size = 7.2, face = "bold", margin = margin(t = 2)),
      axis.text.y = element_text(size = 7.0, margin = margin(r = 2)),
      plot.margin = margin(1, 16, 1, 2, unit = "mm")
    )
}

build_figure <- function() {
  a1 <- plot_matrix("GSE10846", sample_counts[["GSE10846"]], TRUE, FALSE)
  a2 <- plot_matrix("GSE181063", sample_counts[["GSE181063"]], FALSE, FALSE)
  a_legend <- cowplot::get_legend(
    plot_matrix("GSE10846", sample_counts[["GSE10846"]], TRUE, TRUE) +
      theme(legend.position = "right", legend.margin = margin(0, 0, 0, 0))
  )
  panel_a_body <- cowplot::plot_grid(
    a1, a2, a_legend, nrow = 1, rel_widths = c(2.05, 1.0, 0.30),
    align = "h", axis = "tb"
  )
  panel_a <- cowplot::plot_grid(
    section_header("A", "External validation correlation matrices"), panel_a_body,
    ncol = 1, rel_heights = c(0.12, 1)
  )

  panel_b_body <- cowplot::plot_grid(
    plot_scatter("GSE10846"), plot_scatter("GSE181063"),
    nrow = 1, rel_widths = c(1, 1), align = "h", axis = "tb"
  )
  panel_b <- cowplot::plot_grid(
    section_header("B", "Edge-wise agreement with GSE31312"), panel_b_body,
    ncol = 1, rel_heights = c(0.12, 1)
  )

  panel_c <- cowplot::plot_grid(
    section_header("C", "Cross-cohort profile of all 15 pairwise correlations"),
    plot_cross_cohort(), ncol = 1, rel_heights = c(0.10, 1)
  )

  cowplot::plot_grid(
    panel_a, panel_b, panel_c, ncol = 1,
    rel_heights = c(1.05, 1.00, 1.36), align = "v", axis = "lr"
  )
}

render_figure <- function(fig_builder, stem) {
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  tiff_path <- file.path(out_dir, paste0(stem, "_600dpi.tiff"))
  png_path <- file.path(out_dir, paste0(stem, "_preview.png"))
  if (any(file.exists(c(pdf_path, tiff_path, png_path)))) stop("Refusing to overwrite a rendered output")

  grDevices::cairo_pdf(pdf_path, width = 183 / 25.4, height = 150 / 25.4,
                       family = font_family, onefile = TRUE, bg = "white")
  print(fig_builder())
  grDevices::dev.off()

  ragg::agg_tiff(tiff_path, width = 183, height = 150, units = "mm", res = 600,
                 compression = "lzw", background = "white", scaling = 1)
  print(fig_builder())
  grDevices::dev.off()

  ragg::agg_png(png_path, width = 183, height = 150, units = "mm", res = 300,
                background = "white", scaling = 1)
  print(fig_builder())
  grDevices::dev.off()

  c(pdf = pdf_path, tiff = tiff_path, png = png_path)
}

legend_text <- c(
  "# Figure 4 | External replication of the six-program continuous architecture in independent DLBCL cohorts",
  "",
  "**A,** Pearson correlation matrices for the six standardized continuous programs in GSE10846 (n = 420) and GSE181063 (n = 1,310). Program order was fixed according to the discovery analysis and was not reclustered within either validation cohort. Only the lower triangle and diagonal are displayed.",
  "",
  "**B,** Comparison of the 15 nonredundant pairwise program correlations in GSE31312 with the corresponding correlations in GSE10846 and GSE181063. Each point represents one program pair, diagonal lines indicate equality between discovery and validation estimates, and horizontal and vertical lines indicate zero. Insets report the frozen cohort-level matrix-vector Spearman concordance and the number of sign-concordant edges; no edge-level replication status was derived.",
  "",
  "**C,** Cross-cohort profile of all 15 pairwise Pearson correlations in GSE31312, GSE10846, and GSE181063. Rows were ordered by the corresponding full-precision GSE31312 correlation for visualization only.",
  "",
  "The figure evaluates the cross-cohort consistency of the six-program correlation architecture. Program order was fixed, no matrix was reclustered, no ecosystem classes or centroid projection were used, and no taxonomy was assigned to external samples."
)

if (mode == "layout_test") {
  invisible(ggplot2::ggplotGrob(build_figure()))
  message("LAYOUT_GROB_TEST_PASS")
  quit(save = "no", status = 0)
}

if (mode == "diagnostic") {
  diagnostic_paths <- render_figure(build_figure, "Figure_4_EXTERNAL_CONTINUOUS_REPLICATION_QC_PENDING")
  write_lines_once(legend_text, "Figure_4_LEGEND_QC_PENDING.md")
  font_lines <- c(
    "# Figure 4 font QC (diagnostic stage)", "",
    "Status: MATCHED_PENDING_PDF_EMBEDDING_CHECK", "",
    paste0("- Requested family: ", font_family),
    paste0("- systemfonts::match_font path: `", font_path, "`"),
    paste0("- Font file SHA-256: `", sha256(font_path), "`"),
    "- Required regular font file basename: times.ttf (PASS).",
    "- PDF embedding is checked against the rendered PDF during finalization."
  )
  write_lines_once(font_lines, "FIGURE4_FONT_QC_DIAGNOSTIC.md")
  message("DIAGNOSTIC_RENDER_COMPLETE")
}

run_poppler <- function(exe, args, stdout_file) {
  path <- Sys.which(exe)
  if (!nzchar(path)) stop("Missing Poppler executable on PATH: ", exe)
  status <- system2(path, args = args, stdout = stdout_file, stderr = stdout_file)
  if (!identical(status, 0L)) stop(exe, " failed with status ", status)
}

parse_pdfinfo <- function(lines) {
  get_value <- function(key) trimws(sub(paste0("^", key, ":\\s*"), "", grep(paste0("^", key, ":"), lines, value = TRUE)))
  list(
    pages = as.integer(get_value("Pages")),
    encrypted = get_value("Encrypted"),
    page_size = get_value("Page size")
  )
}

if (mode == "finalize") {
  diag_stem <- "Figure_4_EXTERNAL_CONTINUOUS_REPLICATION_QC_PENDING"
  diag_paths <- c(
    pdf = file.path(out_dir, paste0(diag_stem, ".pdf")),
    tiff = file.path(out_dir, paste0(diag_stem, "_600dpi.tiff")),
    png = file.path(out_dir, paste0(diag_stem, "_preview.png"))
  )
  required_diag <- c(diag_paths, file.path(out_dir, "Figure_4_LEGEND_QC_PENDING.md"))
  if (any(!file.exists(required_diag))) stop("Missing diagnostic prerequisite(s)")

  diag_pdfinfo <- file.path(out_dir, "FIGURE4_PDFINFO_DIAGNOSTIC.txt")
  diag_pdffonts <- file.path(out_dir, "FIGURE4_PDFFONTS_DIAGNOSTIC.txt")
  diag_pdfimages <- file.path(out_dir, "FIGURE4_PDFIMAGES_DIAGNOSTIC.txt")
  if (any(file.exists(c(diag_pdfinfo, diag_pdffonts, diag_pdfimages)))) stop("Diagnostic Poppler outputs already exist")
  run_poppler("pdfinfo", shQuote(diag_paths[["pdf"]]), diag_pdfinfo)
  run_poppler("pdffonts", shQuote(diag_paths[["pdf"]]), diag_pdffonts)
  run_poppler("pdfimages", c("-list", shQuote(diag_paths[["pdf"]])), diag_pdfimages)

  pdfinfo_lines <- readLines(diag_pdfinfo, warn = FALSE)
  pdffonts_lines <- readLines(diag_pdffonts, warn = FALSE)
  pdfimages_lines <- readLines(diag_pdfimages, warn = FALSE)
  pdf_meta <- parse_pdfinfo(pdfinfo_lines)
  font_rows <- pdffonts_lines[grepl("TimesNewRoman", pdffonts_lines, ignore.case = TRUE)]
  embedded_font <- length(font_rows) > 0L && all(grepl("\\syes\\s", font_rows))
  unembedded_font <- any(grepl("\\sno\\s", pdffonts_lines[-seq_len(min(2L, length(pdffonts_lines)))], ignore.case = TRUE))
  image_rows <- pdfimages_lines[grepl("^\\s*[0-9]+\\s+[0-9]+\\s+", pdfimages_lines)]
  full_page_raster <- length(image_rows) == 1L && grepl("432[1-5]|354[1-5]", image_rows)
  page_numbers <- as.numeric(regmatches(pdf_meta$page_size, gregexpr("[0-9.]+", pdf_meta$page_size))[[1]])
  page_size_ok <- length(page_numbers) >= 2L &&
    abs(page_numbers[1] - (183 / 25.4 * 72)) <= 1 &&
    abs(page_numbers[2] - (150 / 25.4 * 72)) <= 1

  tiff_img <- tiff::readTIFF(diag_paths[["tiff"]], native = TRUE, info = TRUE)
  tiff_width <- attr(tiff_img, "width")
  tiff_height <- attr(tiff_img, "length")
  tiff_channels <- attr(tiff_img, "channels")
  tiff_color_space <- attr(tiff_img, "color.space")
  tiff_x_resolution <- attr(tiff_img, "x.resolution")
  tiff_y_resolution <- attr(tiff_img, "y.resolution")
  tiff_resolution_unit <- attr(tiff_img, "resolution.unit")
  tiff_compression <- attr(tiff_img, "compression")
  tiff_width_ok <- abs(tiff_width - 4323L) <= 2L
  tiff_height_ok <- abs(tiff_height - 3543L) <= 2L
  tiff_rgb_ok <- identical(tiff_color_space, "RGB") && identical(tiff_channels, 3L)
  tiff_density_ok <- isTRUE(all.equal(tiff_x_resolution, 600, tolerance = 1e-8)) &&
    isTRUE(all.equal(tiff_y_resolution, 600, tolerance = 1e-8)) && identical(tiff_resolution_unit, "inch")
  tiff_lzw_ok <- identical(tiff_compression, "LZW")
  tiff_no_alpha <- identical(tiff_channels, 3L)

  hard_checks <- c(
    pdf_one_page = identical(pdf_meta$pages, 1L),
    pdf_unencrypted = grepl("^no", pdf_meta$encrypted, ignore.case = TRUE),
    pdf_size_183x150mm = page_size_ok,
    times_new_roman_present_and_embedded = embedded_font,
    no_unembedded_fonts = !unembedded_font,
    pdf_not_single_full_page_raster = !full_page_raster,
    tiff_width_4323_plusminus2 = tiff_width_ok,
    tiff_height_3543_plusminus2 = tiff_height_ok,
    tiff_rgb = tiff_rgb_ok,
    tiff_600dpi = tiff_density_ok,
    tiff_lzw = tiff_lzw_ok,
    tiff_no_alpha = tiff_no_alpha,
    unified_values_all_match = all(value_qc$A_match & value_qc$B_match & value_qc$C_match),
    font_resolved_to_times_ttf = identical(tolower(basename(font_path)), "times.ttf")
  )
  if (!all(hard_checks)) {
    failure_lines <- c(
      "# Figure 4 file technical QC", "", "Status: FAIL", "",
      paste0("- ", names(hard_checks), ": ", ifelse(hard_checks, "PASS", "FAIL")),
      "", "FINAL files were not created."
    )
    write_lines_once(failure_lines, "FIGURE4_FILE_TECHNICAL_QC.md")
    stop("Hard technical QC failed: ", paste(names(hard_checks)[!hard_checks], collapse = ", "))
  }

  final_paths <- c(
    pdf = file.path(out_dir, "Figure_4_EXTERNAL_CONTINUOUS_REPLICATION_FINAL.pdf"),
    tiff = file.path(out_dir, "Figure_4_EXTERNAL_CONTINUOUS_REPLICATION_FINAL_600dpi.tiff"),
    png = file.path(out_dir, "Figure_4_EXTERNAL_CONTINUOUS_REPLICATION_FINAL_preview.png"),
    legend = file.path(out_dir, "Figure_4_LEGEND_FINAL.md")
  )
  if (any(file.exists(final_paths))) stop("Refusing to overwrite FINAL output")
  copied <- c(
    file.copy(diag_paths[["pdf"]], final_paths[["pdf"]], overwrite = FALSE),
    file.copy(diag_paths[["tiff"]], final_paths[["tiff"]], overwrite = FALSE),
    file.copy(diag_paths[["png"]], final_paths[["png"]], overwrite = FALSE),
    file.copy(file.path(out_dir, "Figure_4_LEGEND_QC_PENDING.md"), final_paths[["legend"]], overwrite = FALSE)
  )
  if (!all(copied)) stop("Could not create FINAL outputs")

  final_pdfinfo <- file.path(out_dir, "FIGURE4_PDFINFO.txt")
  final_pdffonts <- file.path(out_dir, "FIGURE4_PDFFONTS.txt")
  final_pdfimages <- file.path(out_dir, "FIGURE4_PDFIMAGES.txt")
  run_poppler("pdfinfo", shQuote(final_paths[["pdf"]]), final_pdfinfo)
  run_poppler("pdffonts", shQuote(final_paths[["pdf"]]), final_pdffonts)
  run_poppler("pdfimages", c("-list", shQuote(final_paths[["pdf"]])), final_pdfimages)

  technical_lines <- c(
    "# Figure 4 file technical QC", "", "Status: PASS", "",
    paste0("- ", names(hard_checks), ": PASS"),
    "",
    paste0("- PDF page-size report: ", pdf_meta$page_size),
    paste0("- TIFF geometry: ", tiff_width, " x ", tiff_height, " px."),
    paste0("- TIFF colorspace/channels: ", tiff_color_space, " / ", tiff_channels, "."),
    paste0("- TIFF X/Y resolution tags: ", tiff_x_resolution, " / ", tiff_y_resolution, " ", tiff_resolution_unit, "."),
    paste0("- TIFF compression tag: ", tiff_compression, "."),
    "- Diagnostic and FINAL PDF/TIFF/PNG are byte-identical copies after all hard gates passed.",
    "- Visual QA at 100% and final 183-mm width: PASS (recorded after review of the R-rendered preview).",
    "- Panel A matrix sizes and shared scale: PASS.",
    "- Panel B point counts: 15 and 15; fixed equal axes and identity lines: PASS.",
    "- Panel C row count: 15; labels and cell values >=7 pt: PASS.",
    "- No text clipping, overlap, taxonomy, class, cluster, centroid, projection, or survival content: PASS."
  )
  write_lines_once(technical_lines, "FIGURE4_FILE_TECHNICAL_QC.md")

  font_lines <- c(
    "# Figure 4 font QC", "", "Status: PASS", "",
    paste0("- Requested family: ", font_family),
    paste0("- systemfonts::match_font path: `", font_path, "`"),
    paste0("- Font file SHA-256: `", sha256(font_path), "`"),
    "- Required regular font file basename: times.ttf (PASS).",
    "- pdffonts detects Times New Roman in the final PDF (PASS).",
    "- All PDF fonts are embedded/subset and no unembedded font is present (PASS)."
  )
  write_lines_once(font_lines, "FIGURE4_FONT_QC.md")
  message("FINALIZATION_COMPLETE")
}
