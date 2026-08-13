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


project_root <- DLBCL_PROJECT_ROOT
attempt_root <- file.path(
  project_root,
  "revision_2026_reviewer_response/05e_stage4_GSE31312_execution_attempt2"
)
source_root <- file.path(
  project_root,
  "revision_2026_reviewer_response/05d_stage4c1_score_path_split"
)
output_root <- file.path(attempt_root, "01_score_space_validation")

if (!identical(
  Sys.getenv("DLBCL_REVISION_ALLOW_SCIENCE"),
  "EXPLICITLY_APPROVED_FUTURE_STAGE"
)) stop("Scientific execution token is absent.", call. = FALSE)

write_csv_once <- function(value, name) {
  path <- file.path(output_root, name)
  if (file.exists(path)) stop("Refusing to overwrite: ", path, call. = FALSE)
  utils::write.csv(value, path, row.names = FALSE, na = "")
}

write_lines_once <- function(value, name) {
  path <- file.path(output_root, name)
  if (file.exists(path)) stop("Refusing to overwrite: ", path, call. = FALSE)
  writeLines(value, path, useBytes = TRUE)
}

spaces_path <- file.path(
  attempt_root,
  "score_spaces/STAGE4_HISTORICAL_SCORE_SPACES.rds"
)
if (!file.exists(spaces_path)) stop("Missing validated score spaces.", call. = FALSE)
spaces <- readRDS(spaces_path)
primary <- spaces$historical_untruncated_z
clipped <- spaces$historical_clipped_z
if (!identical(dim(primary), c(498L, 6L))) stop("Primary matrix dimension changed.", call. = FALSE)
if (!identical(dim(clipped), c(498L, 6L))) stop("Clipped matrix dimension changed.", call. = FALSE)
if (!identical(dimnames(primary), dimnames(clipped))) stop("Score-space dimnames differ.", call. = FALSE)
if (!all(is.finite(primary)) || !all(is.finite(clipped))) stop("Non-finite score detected.", call. = FALSE)

config <- utils::read.csv(
  file.path(source_root, "config/STAGE4_INPUT_CONFIG.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
raw_path <- config$absolute_path[
  config$input_id == "gse31312_historical_score_matrix"
]
raw <- utils::read.csv(raw_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!identical(dim(raw), c(6L, 499L))) stop("Raw historical score schema changed.", call. = FALSE)

primary_table <- data.frame(sample = rownames(primary), primary, check.names = FALSE)
untruncated_table <- primary_table
clipped_table <- data.frame(sample = rownames(clipped), clipped, check.names = FALSE)

summary_table <- do.call(rbind, lapply(seq_len(ncol(primary)), function(j) {
  raw_values <- as.numeric(raw[j, -1L])
  data.frame(
    program = colnames(primary)[j],
    raw_mean = mean(raw_values),
    raw_sd = stats::sd(raw_values),
    raw_median = stats::median(raw_values),
    raw_MAD = stats::mad(raw_values),
    raw_minimum = min(raw_values),
    raw_maximum = max(raw_values),
    untruncated_z_mean = mean(primary[, j]),
    untruncated_z_sd = stats::sd(primary[, j]),
    clipped_z_mean = mean(clipped[, j]),
    clipped_z_sd = stats::sd(clipped[, j]),
    stringsAsFactors = FALSE
  )
}))

matrix_to_table <- function(value) {
  data.frame(program = rownames(value), value, check.names = FALSE)
}

s6a_path <- file.path(attempt_root, "score_spaces/S6A_REPRODUCTION_CHECK.csv")
s6b_path <- file.path(attempt_root, "score_spaces/S6B_REPRODUCTION_CHECK.csv")
s6a <- utils::read.csv(s6a_path, stringsAsFactors = FALSE, check.names = FALSE)
s6b <- utils::read.csv(s6b_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(s6a$pass) || !all(s6b$pass)) stop("S6 reproduction gate failed.", call. = FALSE)

write_csv_once(raw, "GSE31312_historical_raw_score_matrix_6x498.csv")
write_csv_once(primary_table, "GSE31312_primary_score_matrix_498x6.csv")
write_csv_once(untruncated_table, "GSE31312_historical_untruncated_z.csv")
write_csv_once(clipped_table, "GSE31312_historical_clipped_z.csv")
write_csv_once(summary_table, "program_summary_statistics.csv")
write_csv_once(matrix_to_table(stats::cor(primary)), "program_correlation_matrix.csv")
write_csv_once(matrix_to_table(stats::cov(primary)), "program_covariance_matrix.csv")
write_csv_once(s6a, "S6A_reproduction_check.csv")
write_csv_once(s6b, "S6B_reproduction_check.csv")

access <- utils::read.csv(
  file.path(attempt_root, "score_spaces/STAGE4C1_INPUT_ACCESS_LOG.csv"),
  stringsAsFactors = FALSE
)
gene_expression_read <- any(access$input_id == "gse31312_gene_expression_matrix")
report <- c(
  "# Score Space Validation Report",
  "",
  "## Observed Validation",
  "",
  "- Raw historical input: 6 program rows and 498 sample columns; identifier column retained separately.",
  "- Primary historical untruncated z score space: 498 samples by 6 programs.",
  "- Historical clipped z score space: 498 samples by 6 programs, element-wise bounded to [-2, 2].",
  paste0("- Maximum absolute S6A reproduction difference: ", format(max(s6a$absolute_difference), scientific = TRUE, digits = 16), "."),
  paste0("- Maximum absolute S6B reproduction difference: ", format(max(s6b$absolute_difference), scientific = TRUE, digits = 16), "."),
  "- S6A and S6B tolerance: 1e-10; both checks passed.",
  paste0("- Maximum absolute untruncated column mean: ", format(max(abs(colMeans(primary))), scientific = TRUE, digits = 8), "."),
  paste0("- Maximum absolute untruncated column SD minus 1: ", format(max(abs(apply(primary, 2, stats::sd) - 1)), scientific = TRUE, digits = 8), "."),
  paste0("- Gene-level expression input read by script 15a: ", toupper(as.character(gene_expression_read)), "."),
  paste0("- GSVA loaded in this formatting process: ", toupper(as.character("GSVA" %in% loadedNamespaces())), "."),
  "",
  "These checks validate faithful reuse of the frozen historical score matrix. They do not constitute re-scoring or Stage 4C-2 sensitivity analysis."
)
write_lines_once(report, "SCORE_SPACE_VALIDATION_REPORT.md")

cat("SCORE_SPACE_VALIDATION_FORMATTED\n")
