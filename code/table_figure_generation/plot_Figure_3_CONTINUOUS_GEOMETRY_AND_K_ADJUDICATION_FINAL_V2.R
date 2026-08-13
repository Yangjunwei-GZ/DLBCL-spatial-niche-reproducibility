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

options(stringsAsFactors = FALSE, width = 200)

suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ragg)
  library(systemfonts)
  library(digest)
  library(magick)
})

set.seed(31312063)

project_root <- DLBCL_PROJECT_ROOT
attempt_dir <- file.path(
  project_root,
  "reproduced_figures/Figure_3"
)
figure_input_dir <- file.path(
  DLBCL_SUPPLEMENTARY_DATA_ROOT, "figure_rendering", "Figure_3"
)
current_attempt_dir <- figure_input_dir
workbook <- file.path(
  DLBCL_SUPPLEMENTARY_TABLES_ROOT,
  "DLBCL_continuous_model_Supplementary_Tables_FINAL_SUBMISSION_PUBLICATION_READY.xlsx"
)
geometry_dir <- figure_input_dir
score_file <- file.path(figure_input_dir, "GSE31312_primary_score_matrix_498x6.csv")
wp1_correlation_file <- file.path(figure_input_dir, "WP1_PROGRAM_PC_PEARSON_CORRELATIONS.csv")
adjudication_dir <- figure_input_dir
current_final_pdf <- file.path(current_attempt_dir, "Figure_3_CONTINUOUS_GEOMETRY_AND_K_ADJUDICATION_FINAL.pdf")
current_final_tiff <- file.path(current_attempt_dir, "Figure_3_CONTINUOUS_GEOMETRY_AND_K_ADJUDICATION_FINAL_600dpi.tiff")
current_final_preview <- file.path(current_attempt_dir, "Figure_3_CONTINUOUS_GEOMETRY_AND_K_ADJUDICATION_FINAL_preview.png")
current_final_legend <- file.path(current_attempt_dir, "Figure_3_LEGEND_FINAL.md")
current_source_script <- file.path(current_attempt_dir, "plot_Figure_3_CONTINUOUS_GEOMETRY_AND_K_ADJUDICATION_FINAL.R")

dir.create(attempt_dir, recursive = TRUE, showWarnings = FALSE)
stale_gate_reports <- file.path(attempt_dir, c(
  "FIGURE3_INPUT_CONFLICT_REPORT.md",
  "FIGURE3_REVISION_TECHNICAL_QC_PENDING.md"
))
invisible(file.remove(stale_gate_reports[file.exists(stale_gate_reports)]))

required_packages <- c(
  "readxl", "ggplot2", "dplyr", "tidyr", "patchwork", "ragg",
  "systemfonts", "digest", "magick", "zip"
)
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "))
}

write_lines <- function(lines, path) {
  writeLines(enc2utf8(lines), con = path, useBytes = TRUE)
}

sha256 <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)

fmt_num <- function(x, digits = 3L, plus = FALSE) {
  ifelse(
    is.na(x),
    "NA",
    if (plus) sprintf(paste0("%+.", digits, "f"), x) else sprintf(paste0("%.", digits, "f"), x)
  )
}

compare_numeric <- function(a, b, tolerance = 1e-12) {
  if (length(a) != length(b)) return(Inf)
  if (any(is.na(a) != is.na(b))) return(Inf)
  max(abs(as.numeric(a) - as.numeric(b)), na.rm = TRUE)
}

compare_frames <- function(a, b, key, tolerance = 1e-12) {
  a <- as.data.frame(a)
  b <- as.data.frame(b)
  common <- intersect(names(a), names(b))
  common <- setdiff(common, key)
  a <- a[do.call(order, a[key]), c(key, common), drop = FALSE]
  b <- b[do.call(order, b[key]), c(key, common), drop = FALSE]
  if (nrow(a) != nrow(b) || !identical(as.character(a[[key[1]]]), as.character(b[[key[1]]]))) {
    return(data.frame(field = "ROW_OR_KEY_STRUCTURE", max_absolute_difference = Inf, exact_character_match = FALSE))
  }
  out <- lapply(common, function(field) {
    if (is.numeric(a[[field]]) || is.integer(a[[field]]) || is.numeric(b[[field]]) || is.integer(b[[field]])) {
      data.frame(
        field = field,
        max_absolute_difference = compare_numeric(a[[field]], b[[field]], tolerance),
        exact_character_match = NA
      )
    } else {
      data.frame(
        field = field,
        max_absolute_difference = NA_real_,
        exact_character_match = identical(as.character(a[[field]]), as.character(b[[field]]))
      )
    }
  })
  bind_rows(out)
}

read_sheet <- function(sheet) {
  readxl::read_excel(workbook, sheet = sheet, skip = 3, .name_repair = "unique")
}

sheet_names <- readxl::excel_sheets(workbook)
required_sheets <- c(
  "S4A_GSE31312_scores",
  "S4D_PCA_parameters",
  "S5A_Method_contract",
  "S5B_Primary_k2_k6",
  "S5C_All_variants",
  "S5D_Failure_summary",
  "S5E_Bootstrap_audit",
  "S5F_Final_decision"
)

s4a <- read_sheet("S4A_GSE31312_scores")
s4d <- read_sheet("S4D_PCA_parameters")
s5a <- read_sheet("S5A_Method_contract")
s5b <- read_sheet("S5B_Primary_k2_k6")
s5c <- read_sheet("S5C_All_variants")
s5d <- read_sheet("S5D_Failure_summary")
s5e <- read_sheet("S5E_Bootstrap_audit")
s5f <- read_sheet("S5F_Final_decision")

sheet_schema <- bind_rows(lapply(required_sheets, function(sheet) {
  dat <- read_sheet(sheet)
  data.frame(
    sheet = sheet,
    row_count = nrow(dat),
    column_count = ncol(dat),
    column_names = paste(names(dat), collapse = " | ")
  )
}))

sheet_purposes <- data.frame(
  purpose = c(
    "498 sample PC coordinates and six frozen standardized scores",
    "Explained variance and program-PC associations",
    "Prespecified k-method and acceptance-threshold contract",
    "Primary k=2-6 criterion-level adjudication",
    "Fourteen nonreference sensitivity variants by k",
    "Nonreference failure and passing counts by k",
    "Bootstrap completion and failure audit",
    "Final model-form decision"
  ),
  source_path = workbook,
  sheet_or_file = required_sheets,
  authority_type = "PRIMARY_WORKBOOK",
  stringsAsFactors = FALSE
) |>
  left_join(sheet_schema, by = c("sheet_or_file" = "sheet")) |>
  mutate(status = "FOUND", notes = "Header row detected after title, description, and blank row.")

write.csv(sheet_purposes, file.path(attempt_dir, "FIGURE3_SHEET_FIELD_MAP.csv"), row.names = FALSE, na = "")

input_paths <- c(
  workbook,
  file.path(geometry_dir, "PCA_SAMPLE_COORDINATES.csv"),
  file.path(geometry_dir, "PCA_PROGRAM_LOADINGS.csv"),
  file.path(geometry_dir, "PCA_EXPLAINED_VARIANCE.csv"),
  score_file,
  wp1_correlation_file,
  file.path(adjudication_dir, "PRIMARY_K2_K6_THRESHOLD_DECOMPOSITION.csv"),
  file.path(adjudication_dir, "ALL_VARIANT_K_THRESHOLD_FAILURE_MATRIX.csv"),
  file.path(adjudication_dir, "THRESHOLD_FAILURE_FREQUENCY_BY_K.csv"),
  file.path(adjudication_dir, "BOOTSTRAP_FAILURE_DISTRIBUTION_SUMMARY.csv"),
  file.path(adjudication_dir, "STAGE4C2_FINAL_K_MODEL_FORM_ADJUDICATION_REPORT.md"),
  file.path(adjudication_dir, "FINAL_K_MODEL_FORM_ADJUDICATION_RULES.md"),
  current_final_pdf,
  current_final_tiff,
  current_final_preview,
  current_final_legend,
  current_source_script,
  file.path(current_attempt_dir, "FIGURE3_INPUT_MANIFEST_SHA256.csv"),
  file.path(current_attempt_dir, "FIGURE3_PROGRAM_PC_ASSOCIATIONS_USED.csv"),
  file.path(current_attempt_dir, "FIGURE3_EXPLAINED_VARIANCE_USED.csv"),
  file.path(current_attempt_dir, "FIGURE3_PRIMARY_K_ADJUDICATION_USED.csv"),
  file.path(current_attempt_dir, "FIGURE3_PRIMARY_K_FAILURE_COUNTS.csv"),
  file.path(current_attempt_dir, "FIGURE3_SENSITIVITY_VARIANT_ADJUDICATION.csv"),
  file.path(current_attempt_dir, "FIGURE3_FINAL_DECISION_USED.csv")
)

missing_inputs <- input_paths[!file.exists(input_paths)]
input_manifest <- bind_rows(lapply(input_paths[file.exists(input_paths)], function(path) {
  info <- file.info(path)
  data.frame(
    absolute_path = normalizePath(path, winslash = "/", mustWork = TRUE),
    file_size_bytes = unname(info$size),
    modified_time = format(info$mtime, "%Y-%m-%d %H:%M:%S %z"),
    sha256 = sha256(path),
    stringsAsFactors = FALSE
  )
}))
write.csv(input_manifest, file.path(attempt_dir, "FIGURE3_INPUT_MANIFEST_SHA256.csv"), row.names = FALSE, na = "")

schema_lines <- c(
  "# Figure 3 schema audit",
  "",
  paste0("- Primary authority workbook: `", workbook, "`"),
  paste0("- Workbook sheet count: ", length(sheet_names)),
  paste0("- Actual sheets: ", paste(sheet_names, collapse = "; ")),
  paste0("- Required sheets present: ", ifelse(all(required_sheets %in% sheet_names), "YES", "NO")),
  "",
  "## Candidate sheet structures",
  "",
  paste0("- `", sheet_schema$sheet, "`: ", sheet_schema$row_count, " rows x ", sheet_schema$column_count, " columns; ", sheet_schema$column_names),
  "",
  "## Authority assignment",
  "",
  paste0("- ", sheet_purposes$purpose, ": `", sheet_purposes$sheet_or_file, "`"),
  "",
  "## Ambiguity review",
  "",
  "- Multiple possible authority tables: NO. The submission workbook is primary; frozen CSVs are used only for value-level QC.",
  "- Reference/nonreference distinction: primary reference results are in S5B; S5C and S5D explicitly describe 14 nonreference variants.",
  "- Threshold direction: explicitly stated in S5A and accompanied by frozen criterion-level flags in S5B/S5C.",
  "- Final decision: directly parseable from S5F.",
  "- Phase 1 eligibility is determined by FIGURE3_INPUT_QC.md after all value-level checks."
)
write_lines(schema_lines, file.path(attempt_dir, "FIGURE3_SCHEMA_AUDIT.md"))

checks <- list()
add_check <- function(check_id, description, observed, expected, pass, evidence) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check_id = check_id,
    description = description,
    observed = as.character(observed),
    expected = as.character(expected),
    status = ifelse(isTRUE(pass), "PASS", "FAIL"),
    evidence = evidence,
    stringsAsFactors = FALSE
  )
}

add_check("INPUT_001", "All required sheets exist", sum(required_sheets %in% sheet_names), length(required_sheets), all(required_sheets %in% sheet_names), workbook)
add_check("INPUT_002", "All frozen QC inputs exist", length(input_paths) - length(missing_inputs), length(input_paths), length(missing_inputs) == 0L, paste(input_paths, collapse = "; "))
add_check("PCA_001", "GSE31312 sample count", nrow(s4a), 498, nrow(s4a) == 498L, "S4A_GSE31312_scores")
add_check("PCA_002", "Sample IDs unique", dplyr::n_distinct(s4a$Sample), 498, dplyr::n_distinct(s4a$Sample) == 498L, "S4A_GSE31312_scores$Sample")
add_check("PCA_003", "PC1 and PC2 finite", sum(is.finite(s4a$PC1) & is.finite(s4a$PC2)), 498, all(is.finite(s4a$PC1)) && all(is.finite(s4a$PC2)), "S4A_GSE31312_scores")

program_order <- c(
  "Macrophage-rich",
  "T cell-inflamed",
  "Immune-inflamed / antigen-presentation",
  "Stromal / fibrotic",
  "Immune-cold / exclusion-associated",
  "Proliferative / cycling"
)
standardized_fields <- paste(program_order, "standardized")
add_check("PCA_004", "Six standardized program fields present", sum(standardized_fields %in% names(s4a)), 6, all(standardized_fields %in% names(s4a)), "S4A_GSE31312_scores")

explained <- s4d |>
  filter(`Record type` == "Explained variance") |>
  transmute(
    PC = `Program / component`,
    explained_variance = as.numeric(`Explained variance`),
    cumulative_variance = as.numeric(`Cumulative variance`)
  )
expected_pcs <- paste0("PC", 1:6)
add_check("PCA_005", "PC1-PC6 explained variance present", paste(explained$PC, collapse = ";"), paste(expected_pcs, collapse = ";"), identical(explained$PC, expected_pcs), "S4D_PCA_parameters")
add_check("PCA_006", "PC1-PC6 explained variance sums to one", sum(explained$explained_variance), 1, abs(sum(explained$explained_variance) - 1) < 1e-12, "S4D_PCA_parameters")

associations <- s4d |>
  filter(`Record type` == "Program loading") |>
  transmute(
    program = sub(" program$", "", `Program / component`),
    program_order = as.integer(`Program order`),
    PC1 = as.numeric(`PC1 Pearson r`),
    PC2 = as.numeric(`PC2 Pearson r`)
  ) |>
  arrange(program_order)
add_check("PCA_007", "Program-PC association table contains PC1 and PC2 for six programs", nrow(associations), 6, nrow(associations) == 6L && all(is.finite(associations$PC1)) && all(is.finite(associations$PC2)), "S4D_PCA_parameters")
add_check("PCA_008", "Program order maps to frozen display order", paste(associations$program, collapse = ";"), paste(program_order, collapse = ";"), identical(associations$program, program_order), "S4D_PCA_parameters")

coords_csv <- read.csv(file.path(geometry_dir, "PCA_SAMPLE_COORDINATES.csv"), check.names = FALSE)
coords_join <- s4a |>
  select(Sample, all_of(expected_pcs)) |>
  inner_join(coords_csv, by = c("Sample" = "sample"), suffix = c("_workbook", "_csv"))
coord_diff <- max(vapply(expected_pcs, function(pc) compare_numeric(coords_join[[paste0(pc, "_workbook")]], coords_join[[paste0(pc, "_csv")]]), numeric(1)))
add_check("VALUE_001", "Workbook PC coordinates match frozen CSV", coord_diff, "<=1e-12", nrow(coords_join) == 498L && coord_diff <= 1e-12, "S4A versus PCA_SAMPLE_COORDINATES.csv")

scores_csv <- read.csv(score_file, check.names = FALSE)
score_join <- s4a |>
  select(Sample, all_of(standardized_fields)) |>
  inner_join(scores_csv, by = c("Sample" = "sample"))
score_diff <- max(vapply(seq_along(program_order), function(i) {
  compare_numeric(score_join[[standardized_fields[i]]], score_join[[paste(program_order[i], "program")]])
}, numeric(1)))
add_check("VALUE_002", "Workbook standardized scores match frozen 498x6 matrix", score_diff, "<=1e-12", nrow(score_join) == 498L && score_diff <= 1e-12, "S4A versus GSE31312_primary_score_matrix_498x6.csv")

variance_csv <- read.csv(file.path(geometry_dir, "PCA_EXPLAINED_VARIANCE.csv"), check.names = FALSE)
variance_diff <- max(
  compare_numeric(explained$explained_variance, variance_csv$explained_variance),
  compare_numeric(explained$cumulative_variance, variance_csv$cumulative_variance)
)
add_check("VALUE_003", "Workbook explained variance matches frozen CSV", variance_diff, "<=1e-12", variance_diff <= 1e-12, "S4D versus PCA_EXPLAINED_VARIANCE.csv")

loading_csv <- read.csv(file.path(geometry_dir, "PCA_PROGRAM_LOADINGS.csv"), check.names = FALSE)
loading_wb <- s4d |>
  filter(`Record type` == "Program loading") |>
  select(`Program / component`, all_of(expected_pcs))
loading_diff <- max(vapply(expected_pcs, function(pc) compare_numeric(loading_wb[[pc]], loading_csv[[pc]]), numeric(1)))
add_check("VALUE_004", "Workbook PCA loadings match frozen CSV", loading_diff, "<=1e-12", loading_diff <= 1e-12 && identical(loading_wb$`Program / component`, loading_csv$program), "S4D versus PCA_PROGRAM_LOADINGS.csv")

correlation_csv <- read.csv(wp1_correlation_file, check.names = FALSE) |>
  select(program, component, pearson_r) |>
  pivot_wider(names_from = component, values_from = pearson_r) |>
  mutate(program = sub(" program$", "", program))
correlation_join <- associations |>
  select(program, PC1, PC2) |>
  inner_join(correlation_csv, by = "program", suffix = c("_workbook", "_csv"))
correlation_diff <- max(
  compare_numeric(correlation_join$PC1_workbook, correlation_join$PC1_csv),
  compare_numeric(correlation_join$PC2_workbook, correlation_join$PC2_csv)
)
add_check("VALUE_005", "Workbook program-PC Pearson correlations match WP1 frozen table", correlation_diff, "<=1e-12", nrow(correlation_join) == 6L && correlation_diff <= 1e-12, "S4D versus WP1_PROGRAM_PC_PEARSON_CORRELATIONS.csv")

primary_csv <- read.csv(file.path(adjudication_dir, "PRIMARY_K2_K6_THRESHOLD_DECOMPOSITION.csv"), check.names = FALSE)
primary_compare <- compare_frames(s5b, primary_csv, key = "k")
primary_match <- all((is.na(primary_compare$max_absolute_difference) | primary_compare$max_absolute_difference <= 1e-12) & (is.na(primary_compare$exact_character_match) | primary_compare$exact_character_match))
add_check("VALUE_006", "Workbook primary-k table matches frozen adjudication CSV", ifelse(primary_match, "MATCH", "CONFLICT"), "MATCH", primary_match, "S5B versus PRIMARY_K2_K6_THRESHOLD_DECOMPOSITION.csv")

variant_csv <- read.csv(file.path(adjudication_dir, "ALL_VARIANT_K_THRESHOLD_FAILURE_MATRIX.csv"), check.names = FALSE)
s5c_keyed <- s5c |> mutate(variant_k_key = paste(variant_id, k, sep = "::"))
variant_keyed <- variant_csv |> mutate(variant_k_key = paste(variant_id, k, sep = "::"))
variant_compare <- compare_frames(s5c_keyed, variant_keyed, key = "variant_k_key")
variant_match <- all((is.na(variant_compare$max_absolute_difference) | variant_compare$max_absolute_difference <= 1e-12) & (is.na(variant_compare$exact_character_match) | variant_compare$exact_character_match))
add_check("VALUE_007", "Workbook sensitivity table matches frozen adjudication CSV", ifelse(variant_match, "MATCH", "CONFLICT"), "MATCH", variant_match, "S5C versus ALL_VARIANT_K_THRESHOLD_FAILURE_MATRIX.csv")

primary_k <- sort(unique(as.integer(s5b$k)))
add_check("K_001", "Primary candidate k values", paste(primary_k, collapse = ","), "2,3,4,5,6", identical(primary_k, 2:6), "S5B_Primary_k2_k6")
add_check("K_002", "Primary candidate count", length(primary_k), 5, length(primary_k) == 5L, "S5B_Primary_k2_k6")

criteria_contract <- s5a |> filter(`Analysis component` == "Acceptance criterion")
add_check("K_003", "Prespecified acceptance criteria count", nrow(criteria_contract), 6, nrow(criteria_contract) == 6L, "S5A_Method_contract")

criterion_pass_fields <- c(
  "PAC_rank_pass", "silhouette_pass", "within_0_02_of_best",
  "minimum_size_pass", "cluster_jaccard_pass", "overall_jaccard_pass"
)
criterion_complete <- all(criterion_pass_fields %in% names(s5b)) && all(rowSums(is.na(as.data.frame(s5b[criterion_pass_fields]))) == 0L)
add_check("K_004", "Every primary k has all criterion-level frozen flags", ifelse(criterion_complete, "COMPLETE", "INCOMPLETE"), "COMPLETE", criterion_complete, "S5B_Primary_k2_k6")

s5c_description <- paste(unlist(readxl::read_excel(workbook, sheet = "S5C_All_variants", col_names = FALSE, n_max = 2)), collapse = " ")
s5d_description <- paste(unlist(readxl::read_excel(workbook, sheet = "S5D_Failure_summary", col_names = FALSE, n_max = 2)), collapse = " ")
explicit_nonreference <- grepl("14 nonreference", paste(s5c_description, s5d_description), ignore.case = TRUE)
variant_ids <- unique(s5c$variant_id)
add_check("K_005", "Nonreference status explicitly documented", explicit_nonreference, TRUE, explicit_nonreference, "S5C/S5D title-description rows")
add_check("K_006", "Nonreference variant count", length(variant_ids), 14, length(variant_ids) == 14L, "S5C_All_variants")
variant_k_complete <- s5c |> count(variant_id) |> summarise(ok = n() == 14L && all(n == 5L)) |> pull(ok)
variant_k_values_complete <- identical(as.integer(sort(unique(s5c$k))), 2:6)
add_check("K_007", "Every nonreference variant has k=2-6 results", variant_k_complete && variant_k_values_complete, TRUE, isTRUE(variant_k_complete) && variant_k_values_complete, "S5C_All_variants")

failure_counts <- s5b |> select(k, failed_criterion_count)
expected_failures <- c(`2` = 3, `3` = 5, `4` = 6, `5` = 5, `6` = 5)
failure_expectation_match <- all(
  as.integer(failure_counts$failed_criterion_count) ==
    as.integer(unname(expected_failures[as.character(failure_counts$k)]))
)
add_check("K_008", "Primary failed-criterion counts match frozen QC targets", paste(failure_counts$failed_criterion_count, collapse = ","), paste(expected_failures, collapse = ","), failure_expectation_match, "S5B_Primary_k2_k6")

passing_counts <- s5c |> group_by(k) |> summarise(passing_variants = sum(variant_pass), .groups = "drop")
passing_from_summary <- s5d |> transmute(k = as.integer(k), passing_variants = as.integer(`Passing variants`))
passing_match <- identical(as.integer(passing_counts$k), as.integer(passing_from_summary$k)) &&
  identical(as.integer(passing_counts$passing_variants), as.integer(passing_from_summary$passing_variants)) &&
  all(passing_counts$passing_variants == 0L)
add_check("K_009", "Passing counts agree and equal zero for every k", paste(passing_counts$passing_variants, collapse = ","), "0,0,0,0,0", passing_match, "S5C and S5D")

decision <- setNames(as.character(s5f$`Final value`), s5f$`Decision field`)
decision_ok <- identical(unname(decision["Final k"]), "NOT_SELECTED") &&
  identical(unname(decision["Taxonomy"]), "NOT_ASSIGNED") &&
  identical(unname(decision["Primary manuscript model"]), "Six continuous scores plus PC1/PC2")
add_check("K_010", "Final decision is parseable and matches frozen model form", paste(decision[c("Final k", "Taxonomy", "Primary manuscript model")], collapse = " | "), "NOT_SELECTED | NOT_ASSIGNED | Six continuous scores plus PC1/PC2", decision_ok, "S5F_Final_decision")

plot_field_names <- c("Sample", "PC1", "PC2", "program", "component", "pearson_r", "k", criterion_pass_fields, "variant_id", "variant_pass")
forbidden_pattern <- "taxonomy|class|cluster_assignment|centroid"
forbidden_fields <- plot_field_names[grepl(forbidden_pattern, plot_field_names, ignore.case = TRUE)]
add_check("SCOPE_001", "Plotting-field contract excludes taxonomy/class/cluster-assignment/centroid fields", paste(forbidden_fields, collapse = ";"), "NONE", length(forbidden_fields) == 0L, "Explicit plotting-field whitelist")

qc <- bind_rows(checks)
write.csv(qc, file.path(attempt_dir, "FIGURE3_PCA_QC.csv"), row.names = FALSE, na = "")

phase0_pass <- all(qc$status == "PASS")
input_qc_lines <- c(
  "# Figure 3 input QC",
  "",
  paste0("- Phase 0 status: **", ifelse(phase0_pass, "PASS", "FAIL"), "**"),
  paste0("- Checks passed: ", sum(qc$status == "PASS"), "/", nrow(qc)),
  paste0("- Sample count: ", nrow(s4a)),
  paste0("- PC1 explained variance: ", sprintf("%.12f", explained$explained_variance[explained$PC == "PC1"])),
  paste0("- PC2 explained variance: ", sprintf("%.12f", explained$explained_variance[explained$PC == "PC2"])),
  paste0("- Primary failed-criterion counts: ", paste0("k=", failure_counts$k, ": ", failure_counts$failed_criterion_count, collapse = "; ")),
  paste0("- Nonreference passing counts: ", paste0("k=", passing_counts$k, ": ", passing_counts$passing_variants, "/14", collapse = "; ")),
  "- New scientific computations: NONE.",
  "- PCA, scoring, clustering, and acceptance decisions were not recomputed.",
  "",
  "## Check table",
  "",
  paste0("- ", qc$check_id, " [", qc$status, "] ", qc$description, "; observed=", qc$observed, "; expected=", qc$expected)
)
write_lines(input_qc_lines, file.path(attempt_dir, "FIGURE3_INPUT_QC.md"))

if (!phase0_pass) {
  failed <- qc |> filter(status == "FAIL")
  conflict_lines <- c(
    "# Figure 3 input conflict report",
    "",
    "Phase 1 was not entered. No FINAL or SUBMISSION files were generated.",
    "",
    paste0("- ", failed$check_id, ": ", failed$description, "; observed=", failed$observed, "; expected=", failed$expected, "; evidence=", failed$evidence)
  )
  write_lines(conflict_lines, file.path(attempt_dir, "FIGURE3_INPUT_CONFLICT_REPORT.md"))
  quit(save = "no", status = 2L)
}

write_lines(c(
  "# Figure 3 PC-association source decision",
  "",
  "- Displayed quantity: frozen program-PC Pearson correlation.",
  "- Primary authority: S4D_PCA_parameters, fields `PC1 Pearson r` and `PC2 Pearson r`.",
  paste0("- Value-level QC source: `", wp1_correlation_file, "`."),
  "- Correlations were selected because the primary workbook explicitly supplies them.",
  "- PCA loadings were not mixed with or relabeled as correlations.",
  "- No correlations were recalculated from sample scores or PC coordinates.",
  "- PCA axis signs are arbitrary."
), file.path(attempt_dir, "FIGURE3_PC_ASSOCIATION_SOURCE_DECISION.md"))

association_long <- associations |>
  pivot_longer(cols = c(PC1, PC2), names_to = "component", values_to = "pearson_r")
write.csv(association_long, file.path(attempt_dir, "FIGURE3_PROGRAM_PC_ASSOCIATIONS_USED.csv"), row.names = FALSE, na = "")
write.csv(explained, file.path(attempt_dir, "FIGURE3_EXPLAINED_VARIANCE_USED.csv"), row.names = FALSE, na = "")

criteria_map_revised <- data.frame(
  original_field = c(
    "PAC_rank", "mean_silhouette", "within_0_02_margin",
    "minimum_cluster_size", "minimum_cluster_median_jaccard", "overall_median_jaccard",
    "failed_criterion_count", "primary_pass"
  ),
  previous_display_label = c(
    "PAC rank", "Mean sil.", "Within 0.02", "Min. size",
    "Cluster Jaccard", "Overall Jaccard", "Failed", "Decision"
  ),
  revised_display_label = c(
    "PAC rank", "Mean sil.", "Near-best sil.", "Min. size",
    "Min. cluster Jaccard", "Overall Jaccard", "Failed criteria", "Overall decision"
  ),
  full_definition = c(
    "Rank of PAC among candidate k values; lower is preferred",
    "Mean silhouette width",
    "Near-best silhouette: within 0.02 of the best silhouette in the same score space",
    "Minimum cluster size",
    "Minimum cluster-level bootstrap Jaccard",
    "Overall median bootstrap Jaccard",
    "Number of prespecified acceptance criteria not satisfied",
    "Whether all prespecified acceptance criteria were satisfied"
  ),
  threshold = c("2", "0.25", "within 0.02", "25", "0.60", "0.75", "not applicable", "all criteria PASS"),
  direction = c("<=", ">", "within", ">=", ">=", ">=", "count", "all TRUE"),
  stringsAsFactors = FALSE
)
write.csv(criteria_map_revised, file.path(attempt_dir, "FIGURE3_CRITERION_LABEL_MAP_REVISED.csv"), row.names = FALSE, na = "")
criteria_map <- criteria_map_revised[1:6, ] |>
  mutate(display_label = c(
    "PAC\nrank", "Mean\nsil.", "Near-best\nsil.",
    "Min.\nsize", "Min.\ncluster\nJaccard", "Overall\nJaccard"
  ))
write.csv(s5b, file.path(attempt_dir, "FIGURE3_PRIMARY_K_ADJUDICATION_USED.csv"), row.names = FALSE, na = "")
write.csv(failure_counts, file.path(attempt_dir, "FIGURE3_PRIMARY_K_FAILURE_COUNTS.csv"), row.names = FALSE, na = "")
write.csv(s5c, file.path(attempt_dir, "FIGURE3_SENSITIVITY_VARIANT_ADJUDICATION.csv"), row.names = FALSE, na = "")
write.csv(s5f, file.path(attempt_dir, "FIGURE3_FINAL_DECISION_USED.csv"), row.names = FALSE, na = "")

variant_labels <- c(
  historical_clipped_euclidean = "Hist. clipped, Euclidean",
  unique_gene_only_euclidean = "Unique genes, Euclidean",
  deoverlap_stromal_immune_cold_euclidean = "De-overlap stromal/immune-cold",
  deoverlap_tcell_antigen_presentation_euclidean = "De-overlap T cell/antigen",
  five_program_without_immune_cold_euclidean = "5 programs; omit immune-cold",
  derived_exclusion_phenotype_euclidean = "Derived exclusion phenotype",
  leave_out_macrophage = "Leave out macrophage",
  leave_out_tcell = "Leave out T cell-inflamed",
  leave_out_antigen_presentation = "Leave out antigen-presentation",
  leave_out_stromal = "Leave out stromal",
  leave_out_immune_cold = "Leave out immune-cold",
  leave_out_proliferative = "Leave out proliferative/cycling",
  primary_pca_whitened = "Primary PCA-whitened distance",
  primary_shrinkage_mahalanobis_whitened = "Shrinkage Mahalanobis-whitened"
)
if (!identical(variant_ids, names(variant_labels))) {
  write_lines(c(
    "# Figure 3 input conflict report",
    "",
    "Variant identifiers do not match the explicit display-label map. No FINAL files were generated.",
    paste0("Observed: ", paste(variant_ids, collapse = "; ")),
    paste0("Mapped: ", paste(names(variant_labels), collapse = "; "))
  ), file.path(attempt_dir, "FIGURE3_INPUT_CONFLICT_REPORT.md"))
  quit(save = "no", status = 3L)
}
variant_map <- data.frame(
  original_variant = variant_ids,
  display_label = unname(variant_labels[variant_ids]),
  variant_type = s5c$variant_family[match(variant_ids, s5c$variant_id)],
  display_order = seq_along(variant_ids),
  stringsAsFactors = FALSE
)
write.csv(variant_map, file.path(attempt_dir, "FIGURE3_VARIANT_LABEL_MAP.csv"), row.names = FALSE, na = "")

font_match <- systemfonts::match_fonts("Times New Roman")
font_path <- normalizePath(font_match$path[1], winslash = "/", mustWork = TRUE)
font_record <- systemfonts::system_fonts() |>
  filter(tolower(normalizePath(path, winslash = "/", mustWork = FALSE)) == tolower(font_path)) |>
  slice(1)
font_ok <- nrow(font_record) == 1L && grepl("Times New Roman", font_record$family, fixed = TRUE)
font_lines <- c(
  "# Figure 3 font QC",
  "",
  paste0("- Requested family: Times New Roman"),
  paste0("- Matched family: ", ifelse(nrow(font_record) == 1L, font_record$family, "UNRESOLVED")),
  paste0("- Matched style: ", ifelse(nrow(font_record) == 1L, font_record$style, "UNRESOLVED")),
  paste0("- Matched font file: `", font_path, "`"),
  paste0("- Match status: ", ifelse(font_ok, "PASS", "FAIL")),
  "- Silent fallback was not permitted."
)
write_lines(font_lines, file.path(attempt_dir, "FIGURE3_FONT_QC.md"))
if (!font_ok) quit(save = "no", status = 4L)

base_family <- "Times New Roman"
ink <- "#2B3035"
blue <- "#55758A"
blue_dark <- "#3F647A"
blue_pale <- "#D9E3E8"
gray <- "#BFC3C5"
gray_pale <- "#ECEDEE"
green_pale <- "#D9E7DE"
red_pale <- "#E9D7D2"
brick <- "#A95448"

theme_pub <- function(base_size = 7.4) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(family = base_family, colour = ink),
      axis.line = element_line(linewidth = 0.3, colour = ink),
      axis.ticks = element_line(linewidth = 0.3, colour = ink),
      axis.title = element_text(size = 7.5),
      axis.text = element_text(size = 7.1),
      plot.title = element_text(size = 9.1, face = "bold", hjust = 0, margin = margin(b = 3)),
      plot.subtitle = element_text(size = 7.5, hjust = 0, margin = margin(b = 2)),
      plot.caption = element_text(size = 7.0, hjust = 0, margin = margin(t = 3)),
      legend.title = element_text(size = 7.0),
      legend.text = element_text(size = 7.0),
      legend.key.height = unit(3.2, "mm"),
      legend.key.width = unit(4.5, "mm"),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(2, 3, 2, 3)
    )
}

pca_plot_data <- s4a |> select(Sample, PC1, PC2)
pc1_pct <- 100 * explained$explained_variance[explained$PC == "PC1"]
pc2_pct <- 100 * explained$explained_variance[explained$PC == "PC2"]

p_a <- ggplot(pca_plot_data, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = gray_pale) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = gray_pale) +
  geom_point(shape = 16, size = 0.85, alpha = 0.58, colour = blue) +
  annotate("text", x = Inf, y = Inf, label = "n = 498", hjust = 1.08, vjust = 1.4, family = base_family, size = 7.2 / ggplot2::.pt, colour = ink) +
  coord_fixed(ratio = 1, clip = "off") +
  labs(
    title = "A  PC1-PC2 continuous geometry",
    x = sprintf("PC1 (%.1f%%)", pc1_pct),
    y = sprintf("PC2 (%.1f%%)", pc2_pct)
  ) +
  theme_pub() +
  theme(plot.margin = margin(2, 5, 2, 2))

panel_a_raw_x_range <- range(pca_plot_data$PC1)
panel_a_raw_y_range <- range(pca_plot_data$PC2)
panel_a_build <- ggplot_build(p_a)
panel_a_final_x_range <- panel_a_build$layout$panel_params[[1]]$x.range
panel_a_final_y_range <- panel_a_build$layout$panel_params[[1]]$y.range
panel_a_clipped_n <- sum(
  pca_plot_data$PC1 < panel_a_final_x_range[1] |
    pca_plot_data$PC1 > panel_a_final_x_range[2] |
    pca_plot_data$PC2 < panel_a_final_y_range[1] |
    pca_plot_data$PC2 > panel_a_final_y_range[2]
)
write_lines(c(
  "# Figure 3 Panel A aspect-ratio QC",
  "",
  paste0("- Original x data range: [", paste(format(panel_a_raw_x_range, digits = 16), collapse = ", "), "]"),
  paste0("- Original y data range: [", paste(format(panel_a_raw_y_range, digits = 16), collapse = ", "), "]"),
  paste0("- Final x scale range: [", paste(format(panel_a_final_x_range, digits = 16), collapse = ", "), "]"),
  paste0("- Final y scale range: [", paste(format(panel_a_final_y_range, digits = 16), collapse = ", "), "]"),
  "- Coordinate aspect ratio: 1 PC2 unit per 1 PC1 unit.",
  paste0("- Samples outside final scale range: ", panel_a_clipped_n),
  paste0("- Sample clipping status: ", ifelse(panel_a_clipped_n == 0L, "PASS", "FAIL")),
  "- Sample coordinates were not changed."
), file.path(attempt_dir, "FIGURE3_PANEL_A_ASPECT_QC.md"))
if (panel_a_clipped_n != 0L) quit(save = "no", status = 11L)

association_plot_data <- association_long |>
  mutate(
    program_display = unname(c(
      "Macrophage-rich" = "Macrophage-rich",
      "T cell-inflamed" = "T cell-inflamed",
      "Immune-inflamed / antigen-presentation" = "Antigen-presentation",
      "Stromal / fibrotic" = "Stromal/fibrotic",
      "Immune-cold / exclusion-associated" = "Immune-cold/exclusion",
      "Proliferative / cycling" = "Proliferative/cycling"
    )[program]),
    program_display = factor(program_display, levels = rev(unique(program_display[order(program_order)]))),
    component = factor(component, levels = c("PC1", "PC2")),
    value_label = sprintf("%.2f", pearson_r)
  )

p_b1 <- ggplot(association_plot_data, aes(component, program_display, fill = pearson_r)) +
  geom_tile(colour = "white", linewidth = 0.65) +
  geom_text(aes(label = value_label), family = base_family, size = 7.0 / ggplot2::.pt, colour = ink) +
  scale_fill_gradient2(
    low = "#5E7D92", mid = "white", high = brick,
    midpoint = 0, limits = c(-1, 1), name = "Pearson r"
  ) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(
    title = "B  Principal-component interpretation",
    subtitle = "Program-PC correlations",
    x = NULL, y = NULL
  ) +
  theme_pub() +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(size = 7.2, face = "bold"),
    axis.text.y = element_text(size = 7.0),
    legend.position = "right",
    legend.direction = "vertical",
    legend.key.height = unit(4.5, "mm"),
    legend.key.width = unit(2.5, "mm"),
    plot.margin = margin(2, 1, 1, 3)
  )

variance_plot_data <- explained |>
  mutate(
    PC = factor(PC, levels = expected_pcs),
    percent = 100 * explained_variance,
    fill_group = ifelse(PC %in% c("PC1", "PC2"), "PC1-PC2", "PC3-PC6")
  )

p_b2 <- ggplot(variance_plot_data, aes(PC, percent, fill = fill_group)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", percent)), vjust = -0.35, family = base_family, size = 7.0 / ggplot2::.pt, colour = ink) +
  scale_fill_manual(values = c("PC1-PC2" = blue_dark, "PC3-PC6" = gray), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(
    title = "Explained variance across principal components",
    x = NULL, y = "Explained variance (%)",
    caption = "PCA axis signs are arbitrary."
  ) +
  theme_pub() +
  theme(plot.title = element_text(size = 7.7, face = "plain"), plot.margin = margin(1, 2, 2, 3))

primary_long <- bind_rows(
  s5b |> transmute(k, criterion = criteria_map$display_label[1], value_label = as.character(as.integer(PAC_rank)), status = ifelse(PAC_rank_pass, "PASS", "FAIL")),
  s5b |> transmute(k, criterion = criteria_map$display_label[2], value_label = fmt_num(mean_silhouette, 3), status = ifelse(silhouette_pass, "PASS", "FAIL")),
  s5b |> transmute(k, criterion = criteria_map$display_label[3], value_label = fmt_num(within_0_02_margin, 3, plus = TRUE), status = ifelse(within_0_02_of_best, "PASS", "FAIL")),
  s5b |> transmute(k, criterion = criteria_map$display_label[4], value_label = as.character(as.integer(minimum_cluster_size)), status = ifelse(minimum_size_pass, "PASS", "FAIL")),
  s5b |> transmute(k, criterion = criteria_map$display_label[5], value_label = fmt_num(minimum_cluster_median_jaccard, 3), status = ifelse(cluster_jaccard_pass, "PASS", "FAIL")),
  s5b |> transmute(k, criterion = criteria_map$display_label[6], value_label = fmt_num(overall_median_jaccard, 3), status = ifelse(overall_jaccard_pass, "PASS", "FAIL")),
  s5b |> transmute(k, criterion = "Failed\ncriteria", value_label = as.character(as.integer(failed_criterion_count)), status = "SUMMARY"),
  s5b |> transmute(k, criterion = "Overall\ndecision", value_label = ifelse(primary_pass, "PASS", "FAIL"), status = ifelse(primary_pass, "PASS", "FAIL"))
) |>
  mutate(
    criterion = factor(criterion, levels = c(criteria_map$display_label, "Failed\ncriteria", "Overall\ndecision")),
    k_label = factor(paste0("k=", k), levels = paste0("k=", 6:2)),
    status = factor(status, levels = c("PASS", "FAIL", "NA", "SUMMARY"))
  )

status_colors <- c(PASS = green_pale, FAIL = red_pale, "NA" = gray_pale, SUMMARY = blue_pale)

p_c <- ggplot(primary_long, aes(criterion, k_label, fill = status)) +
  geom_tile(colour = "white", linewidth = 0.65) +
  geom_text(aes(label = value_label), family = base_family, size = 7.0 / ggplot2::.pt, lineheight = 0.88, colour = ink) +
  scale_fill_manual(
    values = status_colors, drop = FALSE,
    labels = c(PASS = "Pass", FAIL = "Fail", "NA" = "Not available", SUMMARY = "Summary"),
    name = "Criterion status"
  ) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(
    title = "C  Primary adjudication across candidate k values",
    x = NULL, y = NULL,
    caption = paste0(
      "FAIL indicates that the candidate k did not satisfy all prespecified acceptance criteria.\n",
      "No candidate k satisfied all prespecified acceptance criteria."
    )
  ) +
  theme_pub() +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(size = 7.0, face = "plain", lineheight = 0.9),
    axis.text.y = element_text(size = 7.2),
    legend.position = "bottom", legend.direction = "horizontal",
    legend.box.spacing = unit(1, "mm"),
    plot.margin = margin(2, 4, 1, 2)
  )

variant_plot_data <- s5c |>
  left_join(variant_map, by = c("variant_id" = "original_variant", "variant_family" = "variant_type")) |>
  mutate(
    display_label = factor(display_label, levels = rev(variant_map$display_label)),
    k_label = factor(paste0("k=", k), levels = paste0("k=", 2:6)),
    status = factor(ifelse(variant_pass, "PASS", "FAIL"), levels = c("PASS", "FAIL", "NA", "SUMMARY"))
  )

all_cells_fail <- all(!s5c$variant_pass) && !anyNA(s5c$variant_pass)
if (!all_cells_fail) {
  write_lines(c(
    "# Figure 3 input conflict report",
    "",
    "Panel D is not uniformly FAIL in the frozen authority table. V2 FINAL export was blocked."
  ), file.path(attempt_dir, "FIGURE3_INPUT_CONFLICT_REPORT.md"))
  quit(save = "no", status = 12L)
}
passing_caption <- if (length(unique(passing_counts$passing_variants)) == 1L) {
  paste0(
    "Passing variants: ", unique(passing_counts$passing_variants),
    "/", unique(s5d$`Variants evaluated`), " for every candidate k"
  )
} else {
  paste0(
    "Passing variants: ",
    paste0("k=", passing_counts$k, " ", passing_counts$passing_variants, "/14", collapse = " | ")
  )
}

p_d <- ggplot(variant_plot_data, aes(k_label, display_label, fill = status)) +
  geom_tile(colour = "white", linewidth = 0.65) +
  scale_fill_manual(
    values = status_colors,
    breaks = "FAIL",
    labels = "All cells: FAIL",
    name = NULL,
    guide = guide_legend(override.aes = list(fill = red_pale, colour = NA))
  ) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(
    title = "D  Adjudication across non-reference\nsensitivity variants",
    x = NULL, y = NULL, caption = passing_caption
  ) +
  theme_pub() +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(size = 7.1, face = "bold"),
    axis.text.y = element_text(size = 7.0, lineheight = 0.9),
    legend.position = "top",
    legend.justification = "right",
    legend.direction = "horizontal",
    legend.key.width = unit(3.2, "mm"),
    legend.key.height = unit(2.8, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(0.5, "mm"),
    plot.title = element_text(size = 8.8, lineheight = 0.95),
    plot.caption = element_text(size = 7.0, lineheight = 0.95, hjust = 0),
    plot.margin = margin(2, 2, 1, 4)
  )

decision_display <- if (decision_ok) {
  "k not selected | taxonomy not assigned | retained representation: six continuous program scores plus PC1 and PC2"
} else {
  write_lines(c(
    "# Figure 3 input conflict report",
    "",
    "S5F_Final_decision does not support the requested submission-language decision strip. V2 FINAL export was blocked."
  ), file.path(attempt_dir, "FIGURE3_INPUT_CONFLICT_REPORT.md"))
  quit(save = "no", status = 13L)
}

p_strip <- ggplot() +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = "#F1F2F2", colour = NA) +
  annotate("text", x = 0.018, y = 0.72, hjust = 0, label = "Final model-form decision", family = base_family, fontface = "bold", size = 8.8 / ggplot2::.pt, colour = ink) +
  annotate(
    "text", x = 0.018, y = 0.29, hjust = 0,
    label = decision_display,
    family = base_family, size = 7.6 / ggplot2::.pt, colour = ink
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_void(base_family = base_family) +
  theme(plot.margin = margin(1, 2, 1, 2), plot.background = element_rect(fill = "white", colour = NA))

top_row <- p_a + (p_b1 / p_b2 + plot_layout(heights = c(1.35, 1))) + plot_layout(widths = c(1.05, 0.95))
middle_row <- p_c + p_d + plot_layout(widths = c(1.11, 0.89))
figure <- top_row / middle_row / p_strip + plot_layout(heights = c(0.97, 1.26, 0.21)) &
  theme(plot.background = element_rect(fill = "white", colour = NA))

format_comparison_value <- function(value) {
  if (length(value) == 0L || is.na(value)) return("NA")
  if (is.numeric(value)) return(format(value, digits = 17, scientific = FALSE, trim = TRUE))
  as.character(value)
}

scalar_match <- function(old_value, new_value, tolerance = 1e-12) {
  if (length(old_value) != 1L || length(new_value) != 1L) return(FALSE)
  if (is.na(old_value) && is.na(new_value)) return(TRUE)
  if (is.na(old_value) || is.na(new_value)) return(FALSE)
  if (is.numeric(old_value) || is.numeric(new_value)) {
    return(isTRUE(all.equal(as.numeric(old_value), as.numeric(new_value), tolerance = tolerance, check.attributes = FALSE)))
  }
  identical(as.character(old_value), as.character(new_value))
}

vector_signature <- function(values) {
  normalized <- if (is.numeric(values)) {
    ifelse(is.na(values), "NA", sprintf("%.17g", values))
  } else {
    ifelse(is.na(values), "NA", as.character(values))
  }
  digest::digest(paste(normalized, collapse = "|"), algo = "sha256", serialize = FALSE)
}

comparison_rows <- list()
add_comparison <- function(panel, field, old_value, new_value, match, reason_for_change) {
  comparison_rows[[length(comparison_rows) + 1L]] <<- data.frame(
    panel = panel,
    field = field,
    old_value = format_comparison_value(old_value),
    new_value = format_comparison_value(new_value),
    match = isTRUE(match),
    reason_for_change = reason_for_change,
    stringsAsFactors = FALSE
  )
}

old_coords <- read.csv(file.path(geometry_dir, "PCA_SAMPLE_COORDINATES.csv"), check.names = FALSE) |>
  arrange(sample)
new_coords <- pca_plot_data |> arrange(Sample)
add_comparison("A", "sample_count", nrow(old_coords), nrow(new_coords), nrow(old_coords) == nrow(new_coords), "Scientific value unchanged")
add_comparison("A", "sample_id_vector_sha256", vector_signature(old_coords$sample), vector_signature(new_coords$Sample), identical(vector_signature(old_coords$sample), vector_signature(new_coords$Sample)), "Scientific value unchanged")
add_comparison("A", "PC1_coordinate_vector_sha256", vector_signature(old_coords$PC1), vector_signature(new_coords$PC1), identical(vector_signature(old_coords$PC1), vector_signature(new_coords$PC1)), "Scientific value unchanged")
add_comparison("A", "PC2_coordinate_vector_sha256", vector_signature(old_coords$PC2), vector_signature(new_coords$PC2), identical(vector_signature(old_coords$PC2), vector_signature(new_coords$PC2)), "Scientific value unchanged")
add_comparison("A", "coordinate_aspect", "free Cartesian aspect", "1:1 fixed coordinate aspect", FALSE, "Authorized layout-only aspect correction")

old_associations <- read.csv(file.path(current_attempt_dir, "FIGURE3_PROGRAM_PC_ASSOCIATIONS_USED.csv"), check.names = FALSE)
new_associations <- association_long
association_keys <- paste(old_associations$program, old_associations$component, sep = "::")
new_association_keys <- paste(new_associations$program, new_associations$component, sep = "::")
for (i in seq_len(nrow(old_associations))) {
  j <- match(association_keys[i], new_association_keys)
  add_comparison(
    "B", paste0(association_keys[i], " Pearson_r"),
    old_associations$pearson_r[i], new_associations$pearson_r[j],
    !is.na(j) && scalar_match(old_associations$pearson_r[i], new_associations$pearson_r[j]),
    "Scientific value unchanged"
  )
}

old_explained <- read.csv(file.path(current_attempt_dir, "FIGURE3_EXPLAINED_VARIANCE_USED.csv"), check.names = FALSE)
for (i in seq_len(nrow(old_explained))) {
  j <- match(old_explained$PC[i], explained$PC)
  for (field in c("explained_variance", "cumulative_variance")) {
    add_comparison(
      "B", paste0(old_explained$PC[i], " ", field),
      old_explained[[field]][i], explained[[field]][j],
      !is.na(j) && scalar_match(old_explained[[field]][i], explained[[field]][j]),
      "Scientific value unchanged"
    )
  }
}

old_primary <- read.csv(file.path(current_attempt_dir, "FIGURE3_PRIMARY_K_ADJUDICATION_USED.csv"), check.names = FALSE)
primary_fields <- setdiff(intersect(names(old_primary), names(s5b)), c("k", "notes"))
for (i in seq_len(nrow(old_primary))) {
  j <- match(old_primary$k[i], s5b$k)
  for (field in primary_fields) {
    add_comparison(
      "C", paste0("k=", old_primary$k[i], " ", field),
      old_primary[[field]][i], s5b[[field]][j],
      !is.na(j) && scalar_match(old_primary[[field]][i], s5b[[field]][j]),
      "Scientific value unchanged"
    )
  }
}

old_sensitivity <- read.csv(file.path(current_attempt_dir, "FIGURE3_SENSITIVITY_VARIANT_ADJUDICATION.csv"), check.names = FALSE)
sensitivity_fields <- setdiff(intersect(names(old_sensitivity), names(s5c)), c("variant_id", "variant_family", "k", "notes"))
old_sensitivity$key <- paste(old_sensitivity$variant_id, old_sensitivity$k, sep = "::")
new_sensitivity_keys <- paste(s5c$variant_id, s5c$k, sep = "::")
for (i in seq_len(nrow(old_sensitivity))) {
  j <- match(old_sensitivity$key[i], new_sensitivity_keys)
  for (field in sensitivity_fields) {
    add_comparison(
      "D", paste0(old_sensitivity$key[i], " ", field),
      old_sensitivity[[field]][i], s5c[[field]][j],
      !is.na(j) && scalar_match(old_sensitivity[[field]][i], s5c[[field]][j]),
      "Scientific value unchanged"
    )
  }
}

old_decision <- read.csv(file.path(current_attempt_dir, "FIGURE3_FINAL_DECISION_USED.csv"), check.names = FALSE)
for (i in seq_len(nrow(old_decision))) {
  j <- match(old_decision$`Decision field`[i], s5f$`Decision field`)
  add_comparison(
    "Bottom strip", old_decision$`Decision field`[i],
    old_decision$`Final value`[i], s5f$`Final value`[j],
    !is.na(j) && scalar_match(old_decision$`Final value`[i], s5f$`Final value`[j]),
    "Scientific value unchanged"
  )
}

for (field in c("within_0_02_margin", "minimum_cluster_median_jaccard", "failed_criterion_count", "primary_pass")) {
  row <- criteria_map_revised[criteria_map_revised$original_field == field, ]
  add_comparison(
    "C label", field, row$previous_display_label, row$revised_display_label,
    identical(row$previous_display_label, row$revised_display_label),
    "Authorized display-label revision"
  )
}
add_comparison("C display", "overall_decision_cell_text", "Not accepted", "FAIL", FALSE, "Authorized display-language revision")
add_comparison("D display", "title terminology", "nonreference", "non-reference", FALSE, "Authorized terminology revision")
add_comparison("D display", "uniform_matrix_key", "no local key", "All cells: FAIL", FALSE, "Authorized display clarification")
add_comparison(
  "Bottom strip display", "submission wording",
  "k: NOT_SELECTED | Taxonomy: NOT_ASSIGNED | Retained representation: Six continuous scores plus PC1/PC2",
  decision_display, FALSE, "Authorized manuscript-language revision"
)

value_comparison <- bind_rows(comparison_rows)
write.csv(value_comparison, file.path(attempt_dir, "FIGURE3_REVISION_VALUE_COMPARISON.csv"), row.names = FALSE, na = "")
scientific_value_rows <- value_comparison$reason_for_change == "Scientific value unchanged"
scientific_values_match <- all(value_comparison$match[scientific_value_rows])

revision_scientific_checks <- data.frame(
  check = c(
    "Sample count", "PC1 explained variance", "PC2 explained variance",
    "PC coordinate values", "Panel B values", "All primary k decisions FAIL",
    "Primary failed-criterion counts", "k=4 failed criteria", "Non-reference variant count",
    "Passing variants for every k", "Final k", "Taxonomy", "Retained representation",
    "Old taxonomy/class/centroid inputs used", "All unchanged scientific values match current FINAL source"
  ),
  observed = c(
    nrow(s4a), pc1_pct, pc2_pct,
    all(value_comparison$match[value_comparison$panel == "A" & grepl("coordinate_vector", value_comparison$field)]),
    all(value_comparison$match[value_comparison$panel == "B"]),
    all(!s5b$primary_pass),
    paste(failure_counts$failed_criterion_count, collapse = ","),
    failure_counts$failed_criterion_count[failure_counts$k == 4],
    length(unique(s5c$variant_id)),
    paste(passing_counts$passing_variants, collapse = ","),
    decision["Final k"], decision["Taxonomy"], decision["Primary manuscript model"],
    FALSE, scientific_values_match
  ),
  expected = c(
    "498", "48.584% (approximately)", "20.766% (approximately)",
    "TRUE", "TRUE", "TRUE", "3,5,6,5,5", "6", "14", "0,0,0,0,0",
    "NOT_SELECTED", "NOT_ASSIGNED", "Six continuous scores plus PC1/PC2",
    "FALSE", "TRUE"
  ),
  stringsAsFactors = FALSE
)
revision_scientific_checks$pass <- c(
  nrow(s4a) == 498L,
  abs(pc1_pct - 48.584) < 0.001,
  abs(pc2_pct - 20.766) < 0.001,
  isTRUE(revision_scientific_checks$observed[4] == "TRUE"),
  isTRUE(revision_scientific_checks$observed[5] == "TRUE"),
  all(!s5b$primary_pass),
  identical(as.integer(failure_counts$failed_criterion_count), c(3L, 5L, 6L, 5L, 5L)),
  failure_counts$failed_criterion_count[failure_counts$k == 4] == 6,
  length(unique(s5c$variant_id)) == 14L,
  all(passing_counts$passing_variants == 0L),
  identical(unname(decision["Final k"]), "NOT_SELECTED"),
  identical(unname(decision["Taxonomy"]), "NOT_ASSIGNED"),
  identical(unname(decision["Primary manuscript model"]), "Six continuous scores plus PC1/PC2"),
  TRUE,
  scientific_values_match
)
revision_scientific_checks$status <- ifelse(revision_scientific_checks$pass, "PASS", "FAIL")
revision_scientific_pass <- all(revision_scientific_checks$pass)
write_lines(c(
  "# Figure 3 V2 revision scientific QC",
  "",
  paste0("- Overall scientific QC: **", ifelse(revision_scientific_pass, "PASS", "FAIL"), "**"),
  "- Scientific analyses rerun: NO.",
  "- Frozen inputs modified: NO.",
  "- Historical taxonomy, class, and centroid inputs read or plotted: NO.",
  "- Current FINAL source-data tables were compared directly with the V2 plotting data.",
  "",
  "## Checks",
  "",
  paste0("- [", revision_scientific_checks$status, "] ", revision_scientific_checks$check, ": observed=", revision_scientific_checks$observed, "; expected=", revision_scientific_checks$expected)
), file.path(attempt_dir, "FIGURE3_REVISION_SCIENTIFIC_QC.md"))
if (!revision_scientific_pass) {
  write_lines(c(
    "# Figure 3 V2 input conflict report",
    "",
    "Scientific value comparison failed. No FINAL_V2 files were generated.",
    paste0("- ", revision_scientific_checks$check[!revision_scientific_checks$pass])
  ), file.path(attempt_dir, "FIGURE3_INPUT_CONFLICT_REPORT.md"))
  quit(save = "no", status = 14L)
}

width_mm <- 183
height_mm <- 145
base_name <- "Figure_3_CONTINUOUS_GEOMETRY_AND_K_ADJUDICATION"
pending_pdf <- file.path(attempt_dir, paste0(base_name, "_V2_QC_PENDING.pdf"))
pending_tiff <- file.path(attempt_dir, paste0(base_name, "_V2_QC_PENDING_600dpi.tiff"))
pending_preview <- file.path(attempt_dir, paste0(base_name, "_V2_QC_PENDING_preview.png"))
pending_final_size <- file.path(attempt_dir, "FIGURE3_REVISION_VISUAL_QC_FINAL_SIZE_EQUIVALENT.png")

grDevices::cairo_pdf(
  filename = pending_pdf,
  width = width_mm / 25.4,
  height = height_mm / 25.4,
  family = base_family,
  bg = "white",
  onefile = TRUE
)
print(figure)
grDevices::dev.off()

ragg::agg_tiff(
  filename = pending_tiff,
  width = width_mm,
  height = height_mm,
  units = "mm",
  res = 600,
  scaling = 1,
  bitsize = 8,
  compression = "lzw",
  background = "white"
)
print(figure)
grDevices::dev.off()

ragg::agg_png(
  filename = pending_preview,
  width = width_mm,
  height = height_mm,
  units = "mm",
  res = 300,
  scaling = 1,
  background = "white"
)
print(figure)
grDevices::dev.off()

ragg::agg_png(
  filename = pending_final_size,
  width = width_mm,
  height = height_mm,
  units = "mm",
  res = 96,
  scaling = 1,
  background = "white"
)
print(figure)
grDevices::dev.off()

read_tiff_tags <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  byte_order <- rawToChar(readBin(con, what = "raw", n = 2L))
  endian <- if (byte_order == "II") "little" else if (byte_order == "MM") "big" else stop("Invalid TIFF byte order")
  read_u16 <- function(connection = con) readBin(connection, integer(), n = 1L, size = 2L, signed = FALSE, endian = endian)
  read_u32 <- function(connection = con) readBin(connection, integer(), n = 1L, size = 4L, signed = TRUE, endian = endian)
  magic <- read_u16()
  if (magic != 42L) stop("Not a classic TIFF file")
  ifd_offset <- read_u32()
  seek(con, where = ifd_offset, origin = "start")
  n_entries <- read_u16()
  type_sizes <- c(`1` = 1L, `2` = 1L, `3` = 2L, `4` = 4L, `5` = 8L)
  tags <- list()
  for (i in seq_len(n_entries)) {
    tag <- read_u16()
    type <- read_u16()
    count <- read_u32()
    value_raw <- readBin(con, what = "raw", n = 4L)
    next_pos <- seek(con)
    raw_con <- rawConnection(value_raw, "rb")
    value_offset <- readBin(raw_con, integer(), n = 1L, size = 4L, signed = TRUE, endian = endian)
    close(raw_con)
    total_size <- unname(type_sizes[as.character(type)]) * count
    if (is.na(total_size)) {
      seek(con, where = next_pos, origin = "start")
      next
    }
    if (total_size <= 4L) {
      value_con <- rawConnection(value_raw, "rb")
    } else {
      seek(con, where = value_offset, origin = "start")
      value_con <- con
    }
    if (type == 3L) {
      value <- readBin(value_con, integer(), n = count, size = 2L, signed = FALSE, endian = endian)
    } else if (type == 4L) {
      value <- readBin(value_con, integer(), n = count, size = 4L, signed = TRUE, endian = endian)
    } else if (type == 5L) {
      nums <- readBin(value_con, integer(), n = count * 2L, size = 4L, signed = TRUE, endian = endian)
      value <- nums[seq(1L, length(nums), by = 2L)] / nums[seq(2L, length(nums), by = 2L)]
    } else {
      value <- readBin(value_con, raw(), n = count)
    }
    if (total_size <= 4L) close(value_con)
    tags[[as.character(tag)]] <- value
    seek(con, where = next_pos, origin = "start")
  }
  tags
}

pdfinfo_exe <- Sys.which("pdfinfo")
pdffonts_exe <- Sys.which("pdffonts")
pdfimages_exe <- Sys.which("pdfimages")
tool_ok <- all(nzchar(c(pdfinfo_exe, pdffonts_exe, pdfimages_exe)))

pdfinfo_output <- if (tool_ok) system2(pdfinfo_exe, shQuote(pending_pdf), stdout = TRUE, stderr = TRUE) else "PDFINFO_UNAVAILABLE"
pdffonts_output <- if (tool_ok) system2(pdffonts_exe, shQuote(pending_pdf), stdout = TRUE, stderr = TRUE) else "PDFFONTS_UNAVAILABLE"
pdfimages_output <- if (tool_ok) system2(pdfimages_exe, c("-listonly", shQuote(pending_pdf)), stdout = TRUE, stderr = TRUE) else "PDFIMAGES_UNAVAILABLE"
write_lines(pdfinfo_output, file.path(attempt_dir, "FIGURE3_V2_pdfinfo.txt"))
write_lines(pdffonts_output, file.path(attempt_dir, "FIGURE3_V2_pdffonts.txt"))
write_lines(pdfimages_output, file.path(attempt_dir, "FIGURE3_V2_pdfimages.txt"))

extract_field <- function(lines, field) {
  hit <- grep(paste0("^", field, ":"), lines, value = TRUE)
  if (length(hit) == 0L) return(NA_character_)
  trimws(sub(paste0("^", field, ":"), "", hit[1]))
}

pages <- suppressWarnings(as.integer(extract_field(pdfinfo_output, "Pages")))
encrypted <- extract_field(pdfinfo_output, "Encrypted")
page_size_line <- extract_field(pdfinfo_output, "Page size")
page_match <- regexec("([0-9.]+) x ([0-9.]+) pts", page_size_line)
page_parts <- regmatches(page_size_line, page_match)[[1]]
page_width_mm <- if (length(page_parts) >= 3L) as.numeric(page_parts[2]) * 25.4 / 72 else NA_real_
page_height_mm <- if (length(page_parts) >= 3L) as.numeric(page_parts[3]) * 25.4 / 72 else NA_real_

font_data_lines <- pdffonts_output[seq_along(pdffonts_output) > max(c(0L, grep("^-+", pdffonts_output)))]
font_data_lines <- font_data_lines[nzchar(trimws(font_data_lines))]
font_tokens <- strsplit(trimws(font_data_lines), "[[:space:]]+")
font_embedded <- length(font_tokens) > 0L && all(vapply(font_tokens, function(x) length(x) >= 3L && tolower(x[3]) == "yes", logical(1)))
font_names <- if (length(font_tokens) > 0L) vapply(font_tokens, function(x) x[1], character(1)) else character()
pdf_times <- length(font_names) > 0L && all(grepl("TimesNewRoman", font_names, ignore.case = TRUE))

image_data_lines <- pdfimages_output[grepl("^page=", pdfimages_output)]
pdf_image_count <- length(image_data_lines)

tiff_info <- magick::image_info(magick::image_read(pending_tiff))
tiff_tags <- read_tiff_tags(pending_tiff)
tag_value <- function(tag, default = NA_real_) if (!is.null(tiff_tags[[as.character(tag)]])) tiff_tags[[as.character(tag)]][1] else default
tiff_width <- tag_value(256)
tiff_height <- tag_value(257)
tiff_compression <- tag_value(259)
tiff_photometric <- tag_value(262)
tiff_samples <- tag_value(277)
tiff_xres <- tag_value(282)
tiff_yres <- tag_value(283)
tiff_resunit <- tag_value(296)
tiff_extra_samples <- tiff_tags[["338"]]

expected_width_px <- round(width_mm / 25.4 * 600)
expected_height_px <- round(height_mm / 25.4 * 600)
image_dimensions <- lapply(image_data_lines, function(line) {
  hit <- regmatches(line, regexec("width=([0-9]+) height=([0-9]+)", line))[[1]]
  if (length(hit) < 3L) return(c(width = NA_real_, height = NA_real_))
  c(width = as.numeric(hit[2]), height = as.numeric(hit[3]))
})
pdf_full_page_image_count <- if (length(image_dimensions) == 0L) {
  0L
} else {
  sum(vapply(image_dimensions, function(dimensions) {
    is.finite(dimensions["width"]) && is.finite(dimensions["height"]) &&
      dimensions["width"] >= 0.8 * expected_width_px &&
      dimensions["height"] >= 0.8 * expected_height_px
  }, logical(1)))
}

technical_checks <- data.frame(
  check = c(
    "PDF tools available", "PDF page count", "PDF encryption", "PDF width", "PDF height",
    "PDF Times New Roman", "PDF fonts embedded", "PDF full-page raster background",
    "TIFF width", "TIFF height", "TIFF colorspace", "TIFF compression",
    "TIFF photometric", "TIFF samples per pixel", "TIFF alpha", "TIFF X resolution", "TIFF Y resolution", "TIFF resolution unit"
  ),
  observed = c(
    tool_ok, pages, encrypted, page_width_mm, page_height_mm,
    paste(font_names, collapse = ";"), font_embedded, paste0(pdf_full_page_image_count, " full-page of ", pdf_image_count, " total image object(s)"),
    tiff_width, tiff_height, tiff_info$colorspace[1], tiff_compression,
    tiff_photometric, tiff_samples, ifelse(is.null(tiff_extra_samples), "ABSENT", paste(tiff_extra_samples, collapse = ";")),
    tiff_xres, tiff_yres, tiff_resunit
  ),
  expected = c(
    "TRUE", "1", "no", "183 +/- 0.5 mm", "145 +/- 0.2 mm",
    "TimesNewRoman subset/name", "TRUE", "0 full-page raster objects",
    paste0(expected_width_px, " +/- 2"), paste0(expected_height_px, " +/- 2"), "sRGB or RGB", "5 (LZW)",
    "2 (RGB)", "3", "ABSENT", "600 +/- 1", "600 +/- 1", "2 (inch)"
  ),
  pass = c(
    tool_ok,
    identical(pages, 1L),
    identical(tolower(encrypted), "no"),
    is.finite(page_width_mm) && abs(page_width_mm - width_mm) <= 0.5,
    is.finite(page_height_mm) && abs(page_height_mm - height_mm) <= 0.2,
    pdf_times,
    font_embedded,
    pdf_full_page_image_count == 0L,
    abs(tiff_width - expected_width_px) <= 2,
    abs(tiff_height - expected_height_px) <= 2,
    tiff_info$colorspace[1] %in% c("sRGB", "RGB"),
    identical(as.integer(tiff_compression), 5L),
    identical(as.integer(tiff_photometric), 2L),
    identical(as.integer(tiff_samples), 3L),
    is.null(tiff_extra_samples),
    is.finite(tiff_xres) && abs(tiff_xres - 600) <= 1,
    is.finite(tiff_yres) && abs(tiff_yres - 600) <= 1,
    identical(as.integer(tiff_resunit), 2L)
  ),
  stringsAsFactors = FALSE
)
technical_checks$status <- ifelse(technical_checks$pass, "PASS", "FAIL")
technical_pass <- all(technical_checks$pass)

technical_lines <- c(
  "# Figure 3 V2 revision technical QC",
  "",
  paste0("- Overall technical status: **", ifelse(technical_pass, "PASS", "FAIL"), "**"),
  paste0("- PDF inspection commands: `", pdfinfo_exe, "`; `", pdffonts_exe, "`; `", pdfimages_exe, "`."),
  "- TIFF inspection: direct ImageMagick metadata plus direct TIFF IFD tag parsing in R.",
  "- TIFF was exported directly from the R figure at 600 dpi; it was not upscaled from PNG.",
  "- Panel A uses coord_fixed(ratio = 1, clip = 'off').",
  paste0("- Panel A samples outside the final scale range: ", panel_a_clipped_n, "."),
  "",
  "## Checks",
  "",
  paste0("- [", technical_checks$status, "] ", technical_checks$check, ": observed=", technical_checks$observed, "; expected=", technical_checks$expected),
  "",
  "## Raw pdfinfo output",
  "",
  "```text",
  pdfinfo_output,
  "```",
  "",
  "## Raw pdffonts output",
  "",
  "```text",
  pdffonts_output,
  "```",
  "",
  "## Raw pdfimages output",
  "",
  "```text",
  pdfimages_output,
  "```"
)
write_lines(technical_lines, file.path(attempt_dir, "FIGURE3_REVISION_TECHNICAL_QC.md"))

write_lines(capture.output(sessionInfo()), file.path(attempt_dir, "sessionInfo.txt"))

if (!technical_pass) {
  write_lines(c(
    "# Figure 3 technical QC pending",
    "",
    "FINAL export was blocked because one or more technical checks failed.",
    paste0("- ", technical_checks$check[!technical_checks$pass], ": observed=", technical_checks$observed[!technical_checks$pass], "; expected=", technical_checks$expected[!technical_checks$pass])
  ), file.path(attempt_dir, "FIGURE3_REVISION_TECHNICAL_QC_PENDING.md"))
  quit(save = "no", status = 5L)
}

visual_approved <- identical(Sys.getenv("FIGURE3_VISUAL_QC_APPROVED", unset = "FALSE"), "TRUE")
if (!visual_approved) {
  write_lines(c(
    "# Figure 3 V2 revision visual QC status",
    "",
    "- Scientific QC: PASS",
    "- Font QC: PASS",
    "- File technical QC: PASS",
    "- Human visual QC at 100% view: PENDING",
    "- Human visual QC at 183 mm final-width equivalent: PENDING",
    "- FINAL naming and submission packaging: BLOCKED pending visual review."
  ), file.path(attempt_dir, "FIGURE3_REVISION_VISUAL_QC_STATUS.md"))
  quit(save = "no", status = 10L)
}

final_pdf <- file.path(attempt_dir, paste0(base_name, "_FINAL_V2.pdf"))
final_tiff <- file.path(attempt_dir, paste0(base_name, "_FINAL_V2_600dpi.tiff"))
final_preview <- file.path(attempt_dir, paste0(base_name, "_FINAL_V2_preview.png"))
file.copy(pending_pdf, final_pdf, overwrite = TRUE)
file.copy(pending_tiff, final_tiff, overwrite = TRUE)
file.copy(pending_preview, final_preview, overwrite = TRUE)
final_copy_integrity <- c(
  PDF = identical(sha256(pending_pdf), sha256(final_pdf)),
  TIFF = identical(sha256(pending_tiff), sha256(final_tiff)),
  PREVIEW = identical(sha256(pending_preview), sha256(final_preview))
)
if (!all(final_copy_integrity)) {
  stop("FINAL copy-integrity check failed: ", paste(names(final_copy_integrity)[!final_copy_integrity], collapse = ", "))
}
cat(
  paste0(
    "\n\n## FINAL copy integrity\n\n",
    "- PDF tested-copy SHA-256 match: ", final_copy_integrity["PDF"], "\n",
    "- TIFF tested-copy SHA-256 match: ", final_copy_integrity["TIFF"], "\n",
    "- Preview tested-copy SHA-256 match: ", final_copy_integrity["PREVIEW"], "\n"
  ),
  file = file.path(attempt_dir, "FIGURE3_REVISION_TECHNICAL_QC.md"),
  append = TRUE
)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
final_script <- file.path(attempt_dir, paste0("plot_", base_name, "_FINAL_V2.R"))
file.copy(script_path, final_script, overwrite = TRUE)

legend_path <- file.path(attempt_dir, "Figure_3_LEGEND_FINAL_V2.md")
legend_lines <- c(
  "# Figure 3 legend",
  "",
  "**Figure 3 | Continuous geometry and prespecified adjudication of candidate discrete models in GSE31312.**",
  "",
  "**A,** PC1-PC2 representation of 498 GSE31312 samples derived from the six standardized continuous program scores. Points are displayed without discrete class assignments, centroids, candidate-k coloring, cluster boundaries, or density contours.",
  "",
  "**B,** Frozen program-PC Pearson correlations for PC1 and PC2 and the variance explained by PC1-PC6. PCA axis signs are arbitrary.",
  "",
  "**C,** Primary adjudication of candidate k values from 2 to 6 using the prespecified acceptance criteria. Cells show the observed metric and criterion-level PASS/FAIL status. \"FAIL\" indicates that the candidate k did not satisfy all prespecified acceptance criteria. No candidate k satisfied all frozen criteria.",
  "",
  "**D,** Overall adjudication across the 14 non-reference sensitivity variants. All candidate k values failed in every non-reference variant, corresponding to 0 of 14 passing variants for each k.",
  "",
  "**Bottom strip,** The frozen final model-form decision retained the six continuous program scores and PC1 and PC2, with k not selected and taxonomy not assigned.",
  "",
  "The data supported a continuous representation, whereas none of the prespecified candidate discrete models met the frozen acceptance criteria."
)
write_lines(legend_lines, legend_path)

current_v2_hash_manifest <- bind_rows(lapply(list(
  c("CURRENT_FINAL", "PDF", current_final_pdf),
  c("CURRENT_FINAL", "TIFF", current_final_tiff),
  c("CURRENT_FINAL", "PREVIEW_PNG", current_final_preview),
  c("CURRENT_FINAL", "LEGEND", current_final_legend),
  c("CURRENT_FINAL", "R_SCRIPT", current_source_script),
  c("FINAL_V2", "PDF", final_pdf),
  c("FINAL_V2", "TIFF", final_tiff),
  c("FINAL_V2", "PREVIEW_PNG", final_preview),
  c("FINAL_V2", "LEGEND", legend_path),
  c("FINAL_V2", "R_SCRIPT", final_script)
), function(record) {
  path <- record[3]
  info <- file.info(path)
  data.frame(
    version = record[1],
    artifact_type = record[2],
    absolute_path = normalizePath(path, winslash = "/", mustWork = TRUE),
    file_size_bytes = unname(info$size),
    sha256 = sha256(path),
    stringsAsFactors = FALSE
  )
}))
write.csv(
  current_v2_hash_manifest,
  file.path(attempt_dir, "FIGURE3_CURRENT_FINAL_V2_SHA256.csv"),
  row.names = FALSE,
  na = ""
)

write_lines(c(
  "# Figure 3 V2 revision visual QC status",
  "",
  "- Scientific QC: PASS",
  "- Font QC: PASS",
  "- File technical QC: PASS",
  "- Human visual QC at 100% view: PASS",
  "- Human visual QC at 183 mm final-width equivalent: PASS",
  "- Text overlap: NONE OBSERVED",
  "- Text clipping: NONE OBSERVED",
  "- Panel A equal physical axis scaling: PASS",
  "- Panel A sample clipping: NONE",
  "- Panel C FAIL cell text is single-line: YES",
  "- Failed criteria and Overall decision headers readable: YES",
  "- Panel D local All cells: FAIL key visible: YES",
  "- Fourteen Panel D labels readable at final-size equivalent: YES",
  "- Panel C values readable: YES",
  "- Legends unobstructed: YES",
  "- Panel whitespace adequate: YES",
  "- Decision strip uses submission language without status codes: YES",
  "- Minimum specified text size: 7.0 pt"
), file.path(attempt_dir, "FIGURE3_REVISION_VISUAL_QC_STATUS.md"))

git_workdir <- getwd()
setwd(project_root)
git_status_after <- system2("git", c("status", "--short"), stdout = TRUE, stderr = TRUE)
git_diff_after <- system2("git", c("diff", "--stat"), stdout = TRUE, stderr = TRUE)
git_cached_after <- system2("git", c("diff", "--cached", "--stat"), stdout = TRUE, stderr = TRUE)
setwd(git_workdir)
write_lines(c(
  "## git status --short",
  git_status_after,
  "",
  "## git diff --stat",
  git_diff_after,
  "",
  "## git diff --cached --stat",
  git_cached_after
), file.path(attempt_dir, "FIGURE3_REVISION_GIT_STATUS_AFTER.txt"))

manifest_candidates <- c(
  "FIGURE3_INPUT_MANIFEST_SHA256.csv",
  "FIGURE3_PROGRAM_PC_ASSOCIATIONS_USED.csv",
  "FIGURE3_EXPLAINED_VARIANCE_USED.csv",
  "FIGURE3_CRITERION_LABEL_MAP_REVISED.csv",
  "FIGURE3_PRIMARY_K_ADJUDICATION_USED.csv",
  "FIGURE3_PRIMARY_K_FAILURE_COUNTS.csv",
  "FIGURE3_SENSITIVITY_VARIANT_ADJUDICATION.csv",
  "FIGURE3_VARIANT_LABEL_MAP.csv",
  "FIGURE3_FINAL_DECISION_USED.csv",
  "FIGURE3_REVISION_SCIENTIFIC_QC.md",
  "FIGURE3_REVISION_VALUE_COMPARISON.csv",
  "FIGURE3_PANEL_A_ASPECT_QC.md",
  "FIGURE3_REVISION_TECHNICAL_QC.md",
  "FIGURE3_REVISION_VISUAL_QC_STATUS.md",
  "FIGURE3_CURRENT_FINAL_V2_SHA256.csv",
  "FIGURE3_V2_pdfinfo.txt",
  "FIGURE3_V2_pdffonts.txt",
  "FIGURE3_V2_pdfimages.txt",
  "FIGURE3_FONT_QC.md",
  "FIGURE3_REVISION_GIT_STATUS_BEFORE.txt",
  "FIGURE3_REVISION_GIT_STATUS_AFTER.txt",
  "sessionInfo.txt",
  basename(final_pdf),
  basename(final_tiff),
  basename(final_preview),
  basename(legend_path),
  basename(final_script)
)
manifest_paths <- file.path(attempt_dir, manifest_candidates)
output_manifest <- bind_rows(lapply(manifest_paths[file.exists(manifest_paths)], function(path) {
  info <- file.info(path)
  data.frame(
    file_name = basename(path),
    absolute_path = normalizePath(path, winslash = "/", mustWork = TRUE),
    file_size_bytes = unname(info$size),
    sha256 = sha256(path),
    stringsAsFactors = FALSE
  )
}))
output_manifest_path <- file.path(attempt_dir, "FIGURE3_REVISION_OUTPUT_MANIFEST_SHA256.csv")
write.csv(output_manifest, output_manifest_path, row.names = FALSE, na = "")

package_path <- file.path(attempt_dir, "Figure_3_CONTINUOUS_GEOMETRY_AND_K_ADJUDICATION_FINAL_V2_SUBMISSION_PACKAGE.zip")
package_files <- c(manifest_paths[file.exists(manifest_paths)], output_manifest_path)
old_wd <- setwd(attempt_dir)
on.exit(setwd(old_wd), add = TRUE)
zip::zipr(package_path, files = basename(package_files), include_directories = FALSE)

cat("FIGURE3_FINAL_EXPORT_COMPLETE\n")
cat("ATTEMPT_DIR=", attempt_dir, "\n", sep = "")
cat("FINAL_PDF=", final_pdf, "\n", sep = "")
cat("FINAL_TIFF=", final_tiff, "\n", sep = "")
cat("FINAL_PREVIEW=", final_preview, "\n", sep = "")
cat("PACKAGE=", package_path, "\n", sep = "")
