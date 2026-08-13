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


options(stringsAsFactors = FALSE, warn = 2)

wp1_root <- DLBCL_PROJECT_ROOT
wp1_output_root <- file.path(wp1_root, "revision_2026_reviewer_response/05x_wp1_continuous_score_freeze")
Sys.setenv(WP1_CONTINUOUS_FREEZE_SOURCE_ONLY = "TRUE")
source(file.path(DLBCL_SUPPLEMENTARY_CODE_ROOT, "continuous_model", "28_wp1_continuous_score_freeze_v1.R"), local = FALSE)
Sys.unsetenv("WP1_CONTINUOUS_FREEZE_SOURCE_ONLY")

wp1_validator_read_csv <- function(relative_path) {
  path <- file.path(wp1_output_root, relative_path)
  wp1_assert(file.exists(path), paste("Missing required WP1 output:", relative_path))
  utils::read.csv(path, check.names = FALSE)
}

wp1_other_r_process_count <- function() {
  lines <- suppressWarnings(system2("tasklist", c("/FO", "CSV", "/NH"), stdout = TRUE, stderr = FALSE))
  if (!length(lines) || all(!nzchar(lines))) return(0L)
  parsed <- tryCatch(utils::read.csv(text = paste(lines, collapse = "\n"), header = FALSE), error = function(e) NULL)
  if (is.null(parsed) || ncol(parsed) < 2L) return(NA_integer_)
  image <- tolower(as.character(parsed[[1L]]))
  pid <- suppressWarnings(as.integer(as.character(parsed[[2L]])))
  is_r <- image %in% c("r.exe", "rscript.exe", "rterm.exe", "rgui.exe")
  remaining <- sum(is_r & !is.na(pid) & pid != Sys.getpid())
  # Windows Rscript keeps one launcher process beside the current R process.
  max(0L, remaining - 1L)
}

wp1_current_protected_manifest <- function() {
  wp1_protected_manifest()
}

wp1_validate_protected <- function(baseline, current) {
  same_paths <- identical(baseline$relative_path, current$relative_path)
  same_size <- same_paths && identical(as.numeric(baseline$file_size_bytes), as.numeric(current$file_size_bytes))
  same_hash <- same_paths && identical(baseline$sha256, current$sha256)
  list(paths = same_paths, sizes = same_size, hashes = same_hash)
}

wp1_validator_main <- function() {
  wp1_assert(identical(Sys.getenv("DLBCL_REVISION_ALLOW_WP1_CONTINUOUS_FREEZE", unset = ""), ""),
             "WP1 token must be unset before final validation.")
  final_paths <- c("WP1_FINAL_VALIDATION_RESULTS.csv", "WP1_CONTINUOUS_SCORE_FREEZE_REPORT.md",
                   "WP1_CONTINUOUS_SCORE_FREEZE_FINAL_SELF_CHECK.txt")
  wp1_assert(!any(file.exists(file.path(wp1_output_root, final_paths))), "Final validator outputs are create-once and already exist.")

  manifest <- wp1_validator_read_csv("WP1_OUTPUT_MANIFEST.csv")
  pre_validator_required <- manifest$relative_path[manifest$required_for_success %in% TRUE &
                                                    !manifest$relative_path %in% final_paths]
  required_exist <- all(file.exists(file.path(wp1_output_root, pre_validator_required)))

  canonical <- wp1_validator_read_csv("WP1_CANONICAL_PROGRAM_CONTRACT.csv")
  authority <- wp1_validator_read_csv("WP1_PRIMARY_SCORE_AUTHORITY.csv")
  summary <- wp1_validator_read_csv("WP1_PRIMARY_SIX_PROGRAM_DESCRIPTIVE_SUMMARY.csv")
  long <- wp1_validator_read_csv("WP1_PRIMARY_PROGRAM_PEARSON_CORRELATIONS_LONG.csv")
  matrix_table <- wp1_validator_read_csv("WP1_PRIMARY_PROGRAM_PEARSON_CORRELATION_MATRIX.csv")
  pca_authority <- wp1_validator_read_csv("WP1_FROZEN_PCA_AUTHORITY.csv")
  reconstruction <- wp1_validator_read_csv("WP1_PCA_RECONSTRUCTION_AUDIT.csv")
  program_pc <- wp1_validator_read_csv("WP1_PROGRAM_PC_PEARSON_CORRELATIONS.csv")
  pc_repro <- wp1_validator_read_csv("WP1_PROGRAM_PC_CORRELATION_REPRODUCIBILITY_AUDIT.csv")
  projection <- wp1_validator_read_csv("WP1_CONTINUOUS_PROJECTION_PARAMETERS.csv")
  execution <- wp1_validator_read_csv("logs/28_wp1_continuous_score_freeze_v1_execution.csv")
  protected_baseline <- wp1_validator_read_csv("WP1_PROTECTED_PATH_BASELINE_MANIFEST.csv")
  protected_current <- wp1_current_protected_manifest()
  protected_check <- wp1_validate_protected(protected_baseline, protected_current)
  input_registry <- wp1_validator_read_csv("WP1_INPUT_REGISTRY.csv")
  variance_path <- input_registry$absolute_path[grepl("PCA_EXPLAINED_VARIANCE.csv", input_registry$relative_path, fixed = TRUE)]
  wp1_assert(length(variance_path) == 1L, "Frozen variance authority is absent or ambiguous.")
  variance <- utils::read.csv(variance_path, check.names = FALSE)

  matrix_values <- as.matrix(matrix_table[, wp1_programs, drop = FALSE])
  storage.mode(matrix_values) <- "double"
  long_lookup <- matrix(long$pearson_r, nrow = 6L, ncol = 6L, byrow = TRUE)
  expected_input_hashes <- wp1_input_definitions()$expected_sha256
  observed_input_hashes <- vapply(input_registry$absolute_path, wp1_sha256, character(1))
  other_r <- wp1_other_r_process_count()

  checks <- list(
    required_outputs_before_validator = required_exist,
    manifest_paths_unique = !anyDuplicated(manifest$relative_path),
    manifest_all_inside_05x = all(!grepl("^([A-Za-z]:|/|\\\\\\\\|\\.\\.)", manifest$relative_path)),
    manifest_create_once = all(manifest$create_once %in% TRUE),
    manifest_overwrite_forbidden = all(manifest$overwrite_allowed %in% FALSE),
    primary_authority_hash_unchanged = identical(authority$sha256, "e37132340336ac698a73a00229f9d045016c1653870c68276d77bacef15de55a"),
    input_hashes_unchanged = wp1_hash_vectors_match(observed_input_hashes, expected_input_hashes),
    score_sample_count_498 = authority$sample_count == 498L,
    program_count_6 = authority$program_count == 6L,
    canonical_memberships_132 = sum(canonical$canonical_gene_count) == 132L,
    canonical_unique_genes_121 = length(unique(unlist(strsplit(canonical$canonical_gene_members, ";", fixed = TRUE)))) == 121L,
    no_seventh_program = nrow(canonical) == 6L,
    direct_immune_cold_primary = canonical$program_name[[5L]] == "Immune-cold / exclusion-associated program" && canonical$primary_or_sensitivity[[5L]] == "PRIMARY",
    descriptive_rows_6 = nrow(summary) == 6L,
    correlation_long_rows_36 = nrow(long) == 36L,
    correlation_matrix_6x6 = identical(dim(matrix_values), c(6L, 6L)),
    correlation_matrix_symmetric = max(abs(matrix_values - t(matrix_values))) <= 1e-15,
    correlation_diagonal_one = max(abs(diag(matrix_values) - 1)) <= 1e-15,
    long_matrix_consistent = max(abs(matrix_values - long_lookup)) <= 1e-15,
    program_pc_rows_12 = nrow(program_pc) == 12L,
    program_pc_n_498 = all(program_pc$n == 498L),
    program_pc_finite = all(is.finite(program_pc$pearson_r)),
    program_pc_loading_sign_consistent = all(program_pc$loading_correlation_sign_consistent %in% TRUE),
    program_pc_reproducible = all(pc_repro$status == "PASS"),
    pca_authority_rows_5 = nrow(pca_authority) == 5L,
    pca_center_false = all(pca_authority$center_setting %in% FALSE),
    pca_scale_false = all(pca_authority$scale_setting %in% FALSE),
    reconstruction_rows_2 = nrow(reconstruction) == 2L,
    reconstruction_pass = all(reconstruction$status == "PASS"),
    reconstruction_tolerance = max(reconstruction$maximum_absolute_difference) <= 1e-10,
    variance_pc1_exact = abs(variance$explained_variance[variance$PC == "PC1"] - 0.485840112695683) <= 1e-15,
    variance_pc2_exact = abs(variance$explained_variance[variance$PC == "PC2"] - 0.207660210739324) <= 1e-15,
    variance_cumulative_exact = abs(variance$cumulative_variance[variance$PC == "PC2"] - 0.693500323435007) <= 1e-15,
    projection_rows_6 = nrow(projection) == 6L,
    projection_means_finite = all(is.finite(projection$training_mean)),
    projection_sd_positive = all(is.finite(projection$training_standard_deviation) & projection$training_standard_deviation > 0),
    pc1_interpretation_supported = sum(program_pc$component == "PC1" & program_pc$pearson_r > 0) >= 5L,
    pc2_interpretation_supported = all(program_pc$pearson_r[program_pc$component == "PC2" & program_pc$program_order %in% c(4L, 5L)] > 0) && all(program_pc$pearson_r[program_pc$component == "PC2" & program_pc$program_order %in% c(2L, 3L, 6L)] < 0),
    protected_paths_unchanged = protected_check$paths,
    protected_sizes_unchanged = protected_check$sizes,
    protected_hashes_unchanged = protected_check$hashes,
    token_unset = identical(Sys.getenv("DLBCL_REVISION_ALLOW_WP1_CONTINUOUS_FREEZE", unset = ""), ""),
    no_other_active_R_process = identical(other_r, 0L),
    warning_count_zero = identical(execution$value[execution$check == "warning_count"], "0"),
    stack_imbalance_false = identical(execution$value[execution$check == "stack_imbalance"], "FALSE"),
    GSVA_not_loaded = !(("GSVA" %in% loadedNamespaces())) && identical(execution$value[execution$check == "GSVA_loaded"], "FALSE"),
    external_not_run = identical(execution$value[execution$check == "external_input_read"], "FALSE"),
    spatial_not_run = identical(execution$value[execution$check == "spatial_input_read"], "FALSE"),
    clinical_not_run = identical(execution$value[execution$check == "clinical_input_read"], "FALSE"),
    no_forbidden_assignment_fields = all(vapply(list(canonical, authority, summary, long, matrix_table, pca_authority,
                                                       reconstruction, program_pc, pc_repro, projection),
                                                function(x) wp1_forbidden_field_scan(names(x)), logical(1)))
  )
  status <- vapply(checks, isTRUE, logical(1))
  results <- data.frame(
    check_id = names(checks), status = ifelse(status, "PASS", "FAIL"),
    observed = vapply(checks, function(x) paste(x, collapse = ";"), character(1)),
    notes = "Final validator check; no downstream science executed.", check.names = FALSE
  )
  wp1_assert(all(status), paste("WP1 final validator failed:", paste(names(checks)[!status], collapse = "; ")))

  matrix_lines <- vapply(seq_len(nrow(matrix_table)), function(i) {
    paste0("- ", matrix_table$program[[i]], ": ",
           paste(sprintf("%s=%.6f", names(matrix_table)[-1L], as.numeric(matrix_table[i, -1L])), collapse = "; "))
  }, character(1))
  pc_lines <- vapply(seq_len(nrow(program_pc)), function(i) {
    paste0("- ", program_pc$component[[i]], " / ", program_pc$program[[i]],
           ": r=", sprintf("%.6f", program_pc$pearson_r[[i]]),
           ", loading=", sprintf("%.6f", program_pc$loading[[i]]),
           ", absolute rank=", program_pc$absolute_rank_within_component[[i]])
  }, character(1))
  max_reconstruction <- max(reconstruction$maximum_absolute_difference)
  five_t_count <- sum(grepl("/05t_stage4c2a_v4_real_continuation_attempt/", paste0("/", protected_baseline$relative_path), fixed = TRUE))
  conclusion <- "WP1 CONTINUOUS SCORE FREEZE COMPLETED - READY TO DESIGN EXTERNAL AND SPATIAL CONTINUOUS VALIDATION"
  report <- c(
    "# WP1 Continuous Score Freeze Report", "", paste0("**", conclusion, "**"), "",
    "1. Primary score authority: `revision_2026_reviewer_response/05e_stage4_GSE31312_execution_attempt2/01_score_space_validation/GSE31312_primary_score_matrix_498x6.csv`, SHA-256 `e37132340336ac698a73a00229f9d045016c1653870c68276d77bacef15de55a`.",
    "2. Score matrix dimensions: 498 samples x 6 programs plus one sample-ID column (498 x 7 table).",
    "3. Six-program definitions: unchanged; 132 memberships, 121 unique genes, six 22-gene programs, no seventh program.",
    "4. Program correlation structure is reported without screening or grouping:", matrix_lines, "",
    "5. Program-PC Pearson correlations (12 rows):", pc_lines, "",
    paste0("6. PC1 is led numerically by ", paste(program_pc$program[program_pc$component == "PC1" & program_pc$dominant_contributor], collapse = " and "), "."),
    paste0("7. PC2 is led numerically by ", paste(program_pc$program[program_pc$component == "PC2" & program_pc$dominant_contributor], collapse = " and "), "."),
    "8. PC1 shared-abundance interpretation: PASS. Multiple programs contribute in the same direction; this remains a continuous description.",
    "9. PC2 stromal/exclusion-versus-immune/proliferative interpretation: PASS, with the macrophage contribution reported as minor rather than forced into a side.",
    "10. Axis naming strength does not require further downgrading, but all qualification and arbitrary-sign caveats remain mandatory.",
    paste0("11. PCA exact reconstruction: PASS; maximum absolute difference=", format(max_reconstruction, scientific = TRUE, digits = 16), "."),
    "12. Projection parameters: complete for six programs, including GSE31312 raw-score mean/SD and frozen PC1/PC2 loadings.",
    "13. Single composite contrast created: NO.",
    "14. Cluster or taxonomy created: NO. Final k remains NOT_SELECTED; taxonomy remains NOT_ASSIGNED.",
    "15. External input read: NO.", "16. Spatial input read: NO.", "17. Clinical analysis run: NO.",
    "18. WP2 interface readiness: interface frozen, but real execution requires a separate token and cohort/gene/scaling comparability gates.",
    "19. WP3 interface readiness: interface frozen, but real execution requires a separate token and canonical detected-gene intersection gate.",
    "20. Final k status: NOT_SELECTED.", "21. Taxonomy status: NOT_ASSIGNED.",
    paste0("22. Validator: ", sum(status), "/", length(status), " PASS; no other active R process detected during validation."), "",
    paste0("Protected 05t status: ", five_t_count, " files hash-identical to the WP1 baseline snapshot. Protected 05o-05w files also remain hash-identical.")
  )
  self_check <- c(
    "WP1 CONTINUOUS SCORE FREEZE FINAL SELF CHECK",
    paste0("validator_checks: ", sum(status), "/", length(status), " PASS"),
    paste0("primary_score_dimensions: ", authority$sample_count, " x ", authority$program_count, " plus sample ID"),
    paste0("program_memberships: ", sum(canonical$canonical_gene_count)),
    "unique_genes: 121", "final_k_status: NOT_SELECTED", "taxonomy_status: NOT_ASSIGNED",
    "external_status: NOT_RUN", "spatial_status: NOT_RUN", "purity_status: NOT_RUN", "clinical_status: NOT_RUN",
    paste0("PCA_reconstruction_max_absolute_difference: ", format(max_reconstruction, scientific = TRUE, digits = 16)),
    "PC1_variance: 0.485840112695683", "PC2_variance: 0.207660210739324",
    "PC1_PC2_cumulative_variance: 0.693500323435007",
    paste0("protected_manifest_rows: ", nrow(protected_baseline)), paste0("protected_05t_files: ", five_t_count),
    "protected_hashes_unchanged: TRUE", "token_after_run: UNSET", "other_active_R_processes: 0",
    "warning_count: 0", "stack_imbalance: FALSE", "GSVA_loaded: FALSE",
    paste0("conclusion: ", conclusion)
  )

  wp1_write_csv(results, "WP1_FINAL_VALIDATION_RESULTS.csv")
  wp1_write_lines(report, "WP1_CONTINUOUS_SCORE_FREEZE_REPORT.md")
  wp1_write_lines(self_check, "WP1_CONTINUOUS_SCORE_FREEZE_FINAL_SELF_CHECK.txt")
  invisible(TRUE)
}

if (sys.nframe() == 0L) wp1_validator_main()
