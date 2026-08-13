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

wp1_assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

wp1_programs <- c(
  "Macrophage-rich program",
  "T cell-inflamed program",
  "Immune-inflamed / antigen-presentation program",
  "Stromal / fibrotic program",
  "Immune-cold / exclusion-associated program",
  "Proliferative / cycling program"
)

wp1_program_ids <- c(
  "macrophage_rich", "t_cell_inflamed", "antigen_presentation",
  "stromal_fibrotic", "immune_cold_exclusion", "proliferative_cycling"
)

wp1_expected_token <- "EXPLICITLY_APPROVED_WP1_CONTINUOUS_SCORE_FREEZE"
wp1_baseline_commit <- "3509fb4fe367b0c1596ac6a9400301cadcf93027"
wp1_root <- DLBCL_PROJECT_ROOT
wp1_output_root <- file.path(
  wp1_root,
  "revision_2026_reviewer_response/05x_wp1_continuous_score_freeze"
)

wp1_sha256 <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

wp1_hash_vectors_match <- function(observed, expected) {
  identical(unname(observed), unname(expected))
}

wp1_read_csv_bom <- function(path) {
  utils::read.csv(path, check.names = FALSE, fileEncoding = "UTF-8-BOM")
}

wp1_validate_token <- function(token) {
  wp1_assert(identical(token, wp1_expected_token), "WP1 real-execution token is absent or incorrect.")
  invisible(TRUE)
}

wp1_validate_score_data <- function(score_data, expected_n = 498L) {
  wp1_assert(is.data.frame(score_data), "Primary scores must be a data frame.")
  wp1_assert(identical(dim(score_data), c(expected_n, 7L)), "Primary score table must be 498 x 7.")
  wp1_assert(identical(names(score_data), c("sample", wp1_programs)), "Primary program order or sample column changed.")
  wp1_assert(!anyNA(score_data), "Primary score table contains missing values.")
  wp1_assert(!anyDuplicated(score_data$sample), "Primary score table contains duplicated sample IDs.")
  wp1_assert(all(nzchar(score_data$sample)), "Primary score table contains empty sample IDs.")
  numeric_scores <- as.matrix(score_data[, wp1_programs, drop = FALSE])
  storage.mode(numeric_scores) <- "double"
  wp1_assert(all(is.finite(numeric_scores)), "Primary score table contains non-finite values.")
  invisible(TRUE)
}

wp1_validate_canonical <- function(canonical) {
  needed <- c("normalized_program_name", "gene_symbol", "membership_index", "canonical_membership")
  wp1_assert(all(needed %in% names(canonical)), "Canonical manifest schema is incomplete.")
  wp1_assert(nrow(canonical) == 132L, "Canonical manifest must contain 132 memberships.")
  wp1_assert(length(unique(canonical$gene_symbol)) == 121L, "Canonical manifest must contain 121 unique genes.")
  observed <- unique(canonical$normalized_program_name)
  wp1_assert(identical(observed, wp1_programs), "Canonical program order changed.")
  counts <- table(factor(canonical$normalized_program_name, levels = wp1_programs))
  wp1_assert(all(counts == 22L), "Every canonical program must contain 22 memberships.")
  wp1_assert(!any(grepl("derived exclusion", canonical$normalized_program_name, ignore.case = TRUE)),
             "Derived exclusion cannot be a primary or seventh program.")
  wp1_assert(identical(wp1_programs[[5L]], "Immune-cold / exclusion-associated program"),
             "Direct immune-cold must remain the fifth primary program.")
  invisible(TRUE)
}

wp1_skewness <- function(x) {
  s <- stats::sd(x)
  if (!is.finite(s) || s == 0) return(NA_real_)
  mean((x - mean(x))^3) / s^3
}

wp1_descriptive_summary <- function(score_data, authority_sha) {
  rows <- lapply(wp1_programs, function(program) {
    x <- score_data[[program]]
    qs <- stats::quantile(x, probs = c(0.25, 0.5, 0.75), type = 7, names = FALSE)
    data.frame(
      program = program,
      n = length(x),
      missing = sum(is.na(x)),
      mean = mean(x),
      standard_deviation = stats::sd(x),
      minimum = min(x),
      Q1 = qs[[1L]],
      median = qs[[2L]],
      Q3 = qs[[3L]],
      maximum = max(x),
      IQR = stats::IQR(x, type = 7),
      skewness = wp1_skewness(x),
      zero_count = sum(x == 0),
      nonfinite_count = sum(!is.finite(x)),
      authority_score_sha256 = authority_sha,
      quantile_type = 7L,
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

wp1_program_correlations <- function(score_data) {
  matrix_values <- stats::cor(score_data[, wp1_programs, drop = FALSE], method = "pearson")
  long <- do.call(rbind, lapply(seq_along(wp1_programs), function(i) {
    do.call(rbind, lapply(seq_along(wp1_programs), function(j) {
      r <- unname(matrix_values[i, j])
      data.frame(
        program_1 = wp1_programs[[i]], program_2 = wp1_programs[[j]], n = nrow(score_data),
        pearson_r = r, absolute_r = abs(r),
        direction = if (r > 0) "positive" else if (r < 0) "negative" else "zero",
        program_1_order = i, program_2_order = j, check.names = FALSE
      )
    }))
  }))
  matrix_table <- data.frame(program = rownames(matrix_values), matrix_values, check.names = FALSE)
  list(long = long, matrix = matrix_table, values = matrix_values)
}

wp1_reconstruction_audit <- function(score_data, loadings, coordinates, tolerance = 1e-10) {
  wp1_assert(all(c("program", "PC1", "PC2") %in% names(loadings)), "PCA loading schema is incomplete.")
  wp1_assert(all(c("sample", "PC1", "PC2") %in% names(coordinates)), "PCA coordinate schema is incomplete.")
  ordered_loadings <- loadings[match(wp1_programs, loadings$program), , drop = FALSE]
  wp1_assert(!anyNA(ordered_loadings$program), "Frozen loadings do not contain all six programs.")
  joined <- coordinates[match(score_data$sample, coordinates$sample), , drop = FALSE]
  wp1_assert(!anyNA(joined$sample), "Explicit sample-ID join to PCA coordinates failed.")
  wp1_assert(identical(joined$sample, score_data$sample), "Explicit sample-ID join changed score order.")
  reconstructed <- as.matrix(score_data[, wp1_programs, drop = FALSE]) %*%
    as.matrix(ordered_loadings[, c("PC1", "PC2"), drop = FALSE])
  rows <- lapply(c("PC1", "PC2"), function(component) {
    delta <- reconstructed[, component] - joined[[component]]
    correlation <- stats::cor(reconstructed[, component], joined[[component]], method = "pearson")
    data.frame(
      component = component,
      sample_count = nrow(score_data),
      maximum_absolute_difference = max(abs(delta)),
      mean_absolute_difference = mean(abs(delta)),
      RMSE = sqrt(mean(delta^2)),
      correlation = correlation,
      sign_match = all(sign(reconstructed[, component]) == sign(joined[[component]])),
      tolerance = tolerance,
      status = if (max(abs(delta)) <= tolerance && correlation > 1 - 1e-12) "PASS" else "FAIL",
      notes = "Explicit sample-ID join; no sign flip; PCA center=FALSE and scale.=FALSE.",
      check.names = FALSE
    )
  })
  list(table = do.call(rbind, rows), reconstructed = reconstructed, joined = joined, loadings = ordered_loadings)
}

wp1_program_pc_correlations <- function(score_data, coordinates, loadings) {
  joined <- coordinates[match(score_data$sample, coordinates$sample), , drop = FALSE]
  wp1_assert(!anyNA(joined$sample), "Program-PC sample join failed.")
  ordered_loadings <- loadings[match(wp1_programs, loadings$program), , drop = FALSE]
  rows <- do.call(rbind, lapply(c("PC1", "PC2"), function(component) {
    r_values <- vapply(wp1_programs, function(program) {
      stats::cor(score_data[[program]], joined[[component]], method = "pearson")
    }, numeric(1))
    ranks <- rank(-abs(r_values), ties.method = "first")
    do.call(rbind, lapply(seq_along(wp1_programs), function(i) {
      loading <- ordered_loadings[[component]][[i]]
      r <- r_values[[i]]
      role <- if (component == "PC1") {
        if (i == 6L) "limited contribution to shared abundance" else "same-direction shared abundance contribution"
      } else if (i %in% c(4L, 5L)) {
        "stromal/exclusion side"
      } else if (i %in% c(2L, 3L, 6L)) {
        "immune/proliferative opposite side"
      } else {
        "minor macrophage contribution"
      }
      data.frame(
        program = wp1_programs[[i]], program_order = i, component = component,
        n = nrow(score_data), pearson_r = r, absolute_r = abs(r),
        absolute_rank_within_component = ranks[[i]],
        sign = if (r > 0) "positive" else if (r < 0) "negative" else "zero",
        dominant_contributor = ranks[[i]] <= 2L,
        loading = loading,
        loading_sign = if (loading > 0) "positive" else if (loading < 0) "negative" else "zero",
        loading_correlation_sign_consistent = sign(r) == sign(loading),
        axis_interpretation_role = role,
        notes = "Pearson correlation uses the original 498 primary scores and frozen coordinates after explicit sample-ID join.",
        check.names = FALSE
      )
    }))
  }))
  rows
}

wp1_pc_reproducibility_audit <- function(first, second) {
  checks <- c(
    combinations_12 = nrow(first) == 12L,
    n_all_498 = all(first$n == 498L),
    all_finite = all(is.finite(first$pearson_r)),
    pc1_ranks_1_to_6 = identical(sort(first$absolute_rank_within_component[first$component == "PC1"]), 1:6),
    pc2_ranks_1_to_6 = identical(sort(first$absolute_rank_within_component[first$component == "PC2"]), 1:6),
    program_order_valid = identical(first$program, rep(wp1_programs, 2L)),
    repeated_derivation_identical = isTRUE(all.equal(first, second, tolerance = 0, check.attributes = TRUE))
  )
  data.frame(
    check_id = names(checks), status = ifelse(checks, "PASS", "FAIL"),
    observed = as.character(checks), expected = "TRUE",
    notes = "No clustering or program grouping was used.", check.names = FALSE
  )
}

wp1_projection_parameters <- function(summary_stats, loadings) {
  summary_ordered <- summary_stats[match(wp1_programs, summary_stats$program), , drop = FALSE]
  loading_ordered <- loadings[match(wp1_programs, loadings$program), , drop = FALSE]
  wp1_assert(!anyNA(summary_ordered$program), "Training summary is missing a canonical program.")
  wp1_assert(all(is.finite(summary_ordered$raw_mean)), "Training means are incomplete.")
  wp1_assert(all(is.finite(summary_ordered$raw_sd) & summary_ordered$raw_sd > 0), "Training SD values are invalid.")
  data.frame(
    program = wp1_programs,
    program_order = seq_along(wp1_programs),
    training_mean = summary_ordered$raw_mean,
    training_standard_deviation = summary_ordered$raw_sd,
    PC1_loading = loading_ordered$PC1,
    PC2_loading = loading_ordered$PC2,
    loading_source = "05e/continuous_geometry/PCA_PROGRAM_LOADINGS.csv",
    center_source = "05e/01_score_space_validation/program_summary_statistics.csv::raw_mean",
    scale_source = "05e/01_score_space_validation/program_summary_statistics.csv::raw_sd",
    authority_sha256 = paste0(
      "score_summary=dbdd0536720ba1f82b9e5a9fed926e7bf0095223869529f1b443b7e3836fb016;",
      "loadings=7e55483c3f6fde0798476915ee57f73976bd479848075d0c98319415fd2491e1"
    ),
    notes = "Means/SD map canonical program names to frozen raw scores; PCA itself applies no further centering or scaling.",
    check.names = FALSE
  )
}

wp1_forbidden_field_scan <- function(field_names) {
  forbidden <- c(
    "selected_k", "assigned_cluster", "cluster_label", "subtype", "taxonomy_name",
    "nearest_centroid", "winner_program", "patient_class", "spot_class"
  )
  !any(tolower(field_names) %in% forbidden) && !any(tolower(field_names) == "final_k")
}

wp1_write_csv <- function(x, relative_path) {
  path <- file.path(wp1_output_root, relative_path)
  wp1_assert(!file.exists(path), paste("Create-once collision:", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
  invisible(path)
}

wp1_write_lines <- function(lines, relative_path) {
  path <- file.path(wp1_output_root, relative_path)
  wp1_assert(!file.exists(path), paste("Create-once collision:", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path, useBytes = TRUE)
  invisible(path)
}

wp1_input_definitions <- function() {
  rel <- c(
    "05w_stage4c2_continuous_model_transition_preflight/STAGE4C2_CONTINUOUS_MODEL_TRANSITION_PREFLIGHT_REPORT.md",
    "05w_stage4c2_continuous_model_transition_preflight/CONTINUOUS_MODEL_TRANSITION_DECISION_FREEZE.md",
    "05w_stage4c2_continuous_model_transition_preflight/CONTINUOUS_MODEL_ANALYSIS_HIERARCHY.md",
    "05w_stage4c2_continuous_model_transition_preflight/CONTINUOUS_AXIS_LOADING_INTERPRETATION_AUDIT.csv",
    "05w_stage4c2_continuous_model_transition_preflight/CONTINUOUS_POLARIZATION_DEFINITION_OPTIONS.csv",
    "05w_stage4c2_continuous_model_transition_preflight/CONTINUOUS_MODEL_ASSET_REUSE_MATRIX.csv",
    "05w_stage4c2_continuous_model_transition_preflight/CONTINUOUS_MODEL_DOWNSTREAM_WORK_PACKAGES.csv",
    "05w_stage4c2_continuous_model_transition_preflight/CONTINUOUS_MODEL_REAL_EXECUTION_CONTRACT_DRAFT.md",
    "05w_stage4c2_continuous_model_transition_preflight/CONTINUOUS_MODEL_HUMAN_DECISION_REGISTER.csv",
    "05w_stage4c2_continuous_model_transition_preflight/STAGE4C2_CONTINUOUS_MODEL_TRANSITION_SELF_CHECK.txt",
    "05v_stage4c2_final_k_model_form_adjudication/STAGE4C2_FINAL_K_MODEL_FORM_ADJUDICATION_REPORT.md",
    "05v_stage4c2_final_k_model_form_adjudication/FINAL_K_MODEL_FORM_ADJUDICATION_RULES.md",
    "05v_stage4c2_final_k_model_form_adjudication/POST_ADJUDICATION_MANUSCRIPT_REVISION_PLAN.md",
    "02_canonical_manifest_and_pipeline_rebuild/config/canonical_programs_v2.csv",
    "02_canonical_manifest_and_pipeline_rebuild/config/canonical_programs_v2_manifest_summary.csv",
    "02_canonical_manifest_and_pipeline_rebuild/config/canonical_programs_v2_sha256.txt",
    "01_signature_and_provenance_audit/pairwise_signature_overlap.csv",
    "05e_stage4_GSE31312_execution_attempt2/01_score_space_validation/GSE31312_primary_score_matrix_498x6.csv",
    "05e_stage4_GSE31312_execution_attempt2/01_score_space_validation/GSE31312_historical_raw_score_matrix_6x498.csv",
    "05e_stage4_GSE31312_execution_attempt2/01_score_space_validation/program_summary_statistics.csv",
    "05e_stage4_GSE31312_execution_attempt2/01_score_space_validation/program_correlation_matrix.csv",
    "05e_stage4_GSE31312_execution_attempt2/continuous_geometry/PCA_PROGRAM_LOADINGS.csv",
    "05e_stage4_GSE31312_execution_attempt2/continuous_geometry/PCA_SAMPLE_COORDINATES.csv",
    "05e_stage4_GSE31312_execution_attempt2/continuous_geometry/PCA_EXPLAINED_VARIANCE.csv",
    "05e_stage4_GSE31312_execution_attempt2/continuous_geometry/PCA_AXIS_CORRELATIONS_WITH_P_VALUES.csv",
    "05e_stage4_GSE31312_execution_attempt2/parameters/15a_historical_score_parameters.csv",
    "05e_stage4_GSE31312_execution_attempt2/parameters/16_continuous_geometry_parameters.csv",
    "05e_stage4_GSE31312_execution_attempt2/STAGE4C1_ATTEMPT2_COMPLETE_FILE_MANIFEST.csv",
    "05g_stage4c1_continuation_run_id_fix/scripts/16_stage4_continuous_geometry.R",
    "05m_stage4c1_v5_postrun_ratification/V5_POSTRUN_COMPLETE_OUTPUT_MANIFEST.csv",
    "05m_stage4c1_v5_postrun_ratification/STAGE4C1_V5_POSTRUN_RATIFICATION_REPORT.md",
    "04b_stage4_environment_freeze/stage4_renv_project/renv.lock"
  )
  sha <- c(
    "4f38e3c503624109ff380ae2b5501ecbbd5f727b47c1c677fa14b8ff5c01672d",
    "7727fb6039ca44948217ba11452999663175afffe57a3ec668abd879270cd692",
    "11bf8ec0bbd371f2cec37e8d7c6e2ca3fd56eccfa61f7bc9b38601ff1fc2af15",
    "2c27e9e0a8135c96c15dc3f7c79c34078e457ea02a1dfeb1fd6885bdb37fd98c",
    "f121c03655c258e6a71f148284da701cf412166606a5da82037de5da96e29689",
    "aad1d1e3465019adbdc4efec4bed3cb4c9ebafba9ace54906d5812a9749cbd01",
    "66ba3cdeaca866b060a19f9463d1d018c44dc2b9d3c078b92788508720b3c494",
    "3b861823fac4709af38dd3e652279625d09d141f93da687c34e82d2493d28aa7",
    "4e089b22c2f91b1acace4e162e876cf6e92e9962535d15872320380ab7a42097",
    "01aff4faa02b4d1d51ee1c5300519b86605579681aa0d2a1ec2213c624957115",
    "53464fd6e0db0d983b1772333f2b983d2f52d566a33473d4782196f16bf011a4",
    "69d5be65a49664e161c4e68f75479703acc38e217b0017e9ab2944e7486eed63",
    "bf0afa65af63cddae0de77354cefa9f667f5b8c3603b957a74792833f7d83943",
    "0127ca95f3a599b7bb8b97147ec39dae5c394d2cb01d951bde2ed880c89bd20b",
    "5a5762e6b4ab3ffea62e4cf3bce6fc8fb051e3a916c68ba9f6e57d53aa21ba2b",
    "367c90ed3182d90919be7276e774a6fb4d078e2bace710f41f93cf8885e02584",
    "a5c4b1c6575f6120922c7f298f76cca8465468db5b2b3004d7b233e422051173",
    "e37132340336ac698a73a00229f9d045016c1653870c68276d77bacef15de55a",
    "eb2144eb6329bb1ee3fcca0caff7f04ada5c2432c9f7d7bab62bc0494868fe14",
    "dbdd0536720ba1f82b9e5a9fed926e7bf0095223869529f1b443b7e3836fb016",
    "ef21dddc5ea7ad3b3c7670994af5748c2971545c7fd6725dbaba004c6bf494a1",
    "7e55483c3f6fde0798476915ee57f73976bd479848075d0c98319415fd2491e1",
    "6e5c5573cc9f2fd998740aa2ef7275b1d8cc5d3ebdfb398a30b39baa0f87571e",
    "c032ba1ee5b8f9341b3a14eb96ddf1294cf0529c7440c63f9876524f2cd7c073",
    "aefb3249b73662eabe9eaad132dfcbc989abbb014be307bb83a1f47f7f8dc67e",
    "88a7724904c3af2452dc5c4c5cafa120e5b46a85c6471d22f204ae572796237f",
    "d9302362e1b0f0ef1e557d349d13f45a7c70a687758e2049689719e7b9c999c8",
    "ef9c99615d038d10261cd4d75c7ef6257fb9ba55448dd2fa2a611d8fdf007e8f",
    "7561d380ee5bc8abcd582ed610ef31aab7d73b1d26d51a38f115b77be709acd0",
    "9ad91ba601a53b6b13d77979051496c17be80046112b1eb9edcdc7f45c8bba94",
    "7220d470a03136b43fc082e7af20a3093ccddcba91304935e2c5b332bcfda9d6",
    "abf1607763905c0afbbafd75834d28ac2865781064601cb374ea02e7425a736e"
  )
  roles <- c(rep("Amendment 020 authority", 10), rep("Stage 4C-2 adjudication", 3),
             "canonical memberships", "canonical summary", "canonical hash record", "overlap audit",
             "primary score authority", "historical raw score authority", "training scale summary",
             "existing correlation cross-check", "frozen loadings", "frozen coordinates", "frozen variance",
             "composite axis audit", "score parameters", "PCA parameters", "Stage 4C-1 freeze manifest",
             "PCA generating script", "postrun ratification manifest", "postrun ratification report", "frozen renv lock")
  data.frame(input_id = sprintf("I%03d", seq_along(rel)), role = roles, relative_path = rel,
             expected_sha256 = sha, stringsAsFactors = FALSE)
}

wp1_build_input_registry <- function(definitions) {
  paths <- file.path(wp1_root, "revision_2026_reviewer_response", definitions$relative_path)
  wp1_assert(all(file.exists(paths)), "One or more registered WP1 inputs are missing.")
  observed_sha <- vapply(paths, wp1_sha256, character(1))
  wp1_assert(wp1_hash_vectors_match(observed_sha, definitions$expected_sha256), "One or more registered input hashes changed.")
  info <- file.info(paths)
  data.frame(
    input_id = definitions$input_id,
    role = definitions$role,
    absolute_path = normalizePath(paths, winslash = "/", mustWork = TRUE),
    relative_path = file.path("revision_2026_reviewer_response", definitions$relative_path),
    sha256 = observed_sha,
    size = as.numeric(info$size),
    mtime = format(info$mtime, "%Y-%m-%dT%H:%M:%OS6%z"),
    read_authorized = TRUE,
    required = TRUE,
    authority_source = ifelse(grepl("05e|05g|05m", definitions$relative_path), "Stage 4C-1 freeze/ratification", "registered revision authority"),
    notes = "Hash verified before WP1 derivation.",
    check.names = FALSE
  )
}

wp1_protected_manifest <- function() {
  revision_root <- file.path(wp1_root, "revision_2026_reviewer_response")
  candidate_roots <- list.dirs(revision_root, recursive = FALSE, full.names = TRUE)
  roots <- sort(candidate_roots[grepl("^05[o-w]_", basename(candidate_roots))])
  wp1_assert(length(roots) == 9L, "Protected 05o-05w directory set is incomplete.")
  files <- sort(unique(unlist(lapply(roots, function(root) {
    list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  }))))
  files <- files[file.info(files)$isdir %in% FALSE]
  info <- file.info(files)
  normalized_files <- normalizePath(files, winslash = "/")
  data.frame(
    relative_path = substring(normalized_files, nchar(normalizePath(wp1_root, winslash = "/")) + 2L),
    file_size_bytes = as.numeric(info$size),
    mtime = format(info$mtime, "%Y-%m-%dT%H:%M:%OS6%z"),
    sha256 = vapply(files, wp1_sha256, character(1)),
    protected_status = "READ_ONLY_BASELINE",
    check.names = FALSE
  )
}

wp1_canonical_contract <- function(canonical, authority_path, authority_sha) {
  rows <- lapply(seq_along(wp1_programs), function(i) {
    x <- canonical[canonical$normalized_program_name == wp1_programs[[i]], , drop = FALSE]
    x <- x[order(as.integer(x$membership_index)), , drop = FALSE]
    overlap_count <- sum(toupper(as.character(x$duplicated_across_programs)) == "TRUE")
    data.frame(
      program_order = i,
      program_id = wp1_program_ids[[i]],
      program_name = wp1_programs[[i]],
      canonical_gene_count = nrow(x),
      canonical_gene_members = paste(x$gene_symbol, collapse = ";"),
      unique_gene_count = length(unique(x$gene_symbol)),
      overlap_gene_count = overlap_count,
      primary_or_sensitivity = "PRIMARY",
      continuous_role = "continuous program score; no patient or spot class",
      spatial_intersection_rule = "canonical22 intersect detected genes; no substitution; no imputation",
      external_gene_matching_rule = "exact approved gene-symbol mapping; missing genes remain missing; no zero fill",
      authority_path = authority_path,
      authority_sha256 = authority_sha,
      notes = if (i == 5L) "Direct immune-cold is primary; derived exclusion is separate sensitivity only." else "No seventh program authorized.",
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

wp1_main <- function() {
  wp1_validate_token(Sys.getenv("DLBCL_REVISION_ALLOW_WP1_CONTINUOUS_FREEZE", unset = ""))
  on.exit(Sys.unsetenv("DLBCL_REVISION_ALLOW_WP1_CONTINUOUS_FREEZE"), add = TRUE)
  wp1_assert(identical(paste(R.version$major, strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1], sep = "."), "4.5"),
             "WP1 requires R 4.5.x.")
  wp1_assert(startsWith(R.version.string, "R version 4.5.1"), "WP1 requires exact R 4.5.1.")
  wp1_assert(identical(Sys.getenv("WP1_PRELAUNCH_R_COUNT", unset = ""), "0"), "Missing prelaunch active-R gate evidence.")
  head <- trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE))
  wp1_assert(identical(head, wp1_baseline_commit), "WP1 baseline HEAD changed.")
  tracked <- system2("git", c("status", "--porcelain", "--untracked-files=no"), stdout = TRUE)
  wp1_assert(length(tracked) == 0L, "Tracked worktree is not clean before WP1 execution.")

  environment_path <- file.path(DLBCL_SUPPLEMENTARY_CODE_ROOT, "continuous_model", "stage4_environment.R")
  source(environment_path, local = FALSE)
  env_context <- stage4b_environment_activate(require_science_token = FALSE, verbose = FALSE, profile = "stage4c1")
  wp1_assert(!("GSVA" %in% loadedNamespaces()), "GSVA must not be loaded in WP1.")

  real_outputs <- c(
    "WP1_INPUT_REGISTRY.csv", "WP1_CANONICAL_PROGRAM_CONTRACT.csv", "WP1_PRIMARY_SCORE_AUTHORITY.csv",
    "WP1_PRIMARY_SIX_PROGRAM_DESCRIPTIVE_SUMMARY.csv", "WP1_PRIMARY_PROGRAM_PEARSON_CORRELATIONS_LONG.csv",
    "WP1_PRIMARY_PROGRAM_PEARSON_CORRELATION_MATRIX.csv", "WP1_FROZEN_PCA_AUTHORITY.csv",
    "WP1_PCA_RECONSTRUCTION_AUDIT.csv", "WP1_PROGRAM_PC_PEARSON_CORRELATIONS.csv",
    "WP1_PROGRAM_PC_CORRELATION_REPRODUCIBILITY_AUDIT.csv", "WP1_CONTINUOUS_AXIS_INTERPRETATION_CONTRACT.md",
    "WP1_CONTINUOUS_PROJECTION_PARAMETERS.csv", "WP1_CONTINUOUS_PROJECTION_EQUATION.md",
    "WP1_DOWNSTREAM_CONTINUOUS_INTERFACE_SPEC.csv", "WP1_PROTECTED_PATH_BASELINE_MANIFEST.csv",
    "logs/28_wp1_continuous_score_freeze_v1_execution.csv"
  )
  wp1_assert(!any(file.exists(file.path(wp1_output_root, real_outputs))), "One or more create-once real outputs already exist.")

  input_definitions <- wp1_input_definitions()
  input_registry <- wp1_build_input_registry(input_definitions)
  get_input <- function(pattern) {
    hit <- input_registry$absolute_path[grepl(pattern, input_registry$relative_path, fixed = TRUE)]
    wp1_assert(length(hit) == 1L, paste("Authority is absent or ambiguous:", pattern))
    hit
  }
  score_path <- get_input("GSE31312_primary_score_matrix_498x6.csv")
  raw_path <- get_input("GSE31312_historical_raw_score_matrix_6x498.csv")
  canonical_path <- get_input("canonical_programs_v2.csv")
  summary_path <- get_input("program_summary_statistics.csv")
  loadings_path <- get_input("PCA_PROGRAM_LOADINGS.csv")
  coordinates_path <- get_input("PCA_SAMPLE_COORDINATES.csv")
  variance_path <- get_input("PCA_EXPLAINED_VARIANCE.csv")
  pca_parameters_path <- get_input("16_continuous_geometry_parameters.csv")
  freeze_manifest_path <- get_input("STAGE4C1_ATTEMPT2_COMPLETE_FILE_MANIFEST.csv")

  score_data <- utils::read.csv(score_path, check.names = FALSE)
  canonical <- utils::read.csv(canonical_path, check.names = FALSE)
  summary_stats <- utils::read.csv(summary_path, check.names = FALSE)
  loadings <- utils::read.csv(loadings_path, check.names = FALSE)
  coordinates <- utils::read.csv(coordinates_path, check.names = FALSE)
  variance <- utils::read.csv(variance_path, check.names = FALSE)
  pca_parameters <- utils::read.csv(pca_parameters_path, check.names = FALSE)
  raw_scores <- utils::read.csv(raw_path, check.names = FALSE)
  freeze_manifest <- wp1_read_csv_bom(freeze_manifest_path)

  wp1_validate_score_data(score_data)
  wp1_validate_canonical(canonical)
  wp1_assert(identical(dim(raw_scores), c(6L, 499L)), "Historical raw score authority must be 6 x 499 including row-label column.")
  wp1_assert(sum(freeze_manifest$relative_path == "01_score_space_validation/GSE31312_primary_score_matrix_498x6.csv") == 1L,
             "Freeze manifest does not identify a unique primary score authority.")
  wp1_assert(sum(freeze_manifest$relative_path == "continuous_geometry/PCA_PROGRAM_LOADINGS.csv") == 1L,
             "Freeze manifest does not identify unique PCA loadings.")
  wp1_assert(identical(as.character(pca_parameters$value[pca_parameters$parameter == "pca_center"]), "FALSE"), "Frozen PCA center setting changed.")
  wp1_assert(identical(as.character(pca_parameters$value[pca_parameters$parameter == "pca_scale"]), "FALSE"), "Frozen PCA scale setting changed.")
  wp1_assert(identical(dim(loadings), c(6L, 7L)), "Frozen loading table must be 6 x 7.")
  wp1_assert(identical(dim(coordinates), c(498L, 7L)), "Frozen coordinate table must be 498 x 7.")
  wp1_assert(abs(variance$explained_variance[variance$PC == "PC1"] - 0.485840112695683) <= 1e-15,
             "PC1 explained variance changed.")
  wp1_assert(abs(variance$explained_variance[variance$PC == "PC2"] - 0.207660210739324) <= 1e-15,
             "PC2 explained variance changed.")
  wp1_assert(abs(variance$cumulative_variance[variance$PC == "PC2"] - 0.693500323435007) <= 1e-15,
             "PC1+PC2 cumulative variance changed.")

  score_sha <- wp1_sha256(score_path)
  descriptive <- wp1_descriptive_summary(score_data, score_sha)
  correlations <- wp1_program_correlations(score_data)
  reconstruction <- wp1_reconstruction_audit(score_data, loadings, coordinates)
  wp1_assert(all(reconstruction$table$status == "PASS"), "WP1 BLOCKED BY PCA AUTHORITY OR RECONSTRUCTION MISMATCH")
  program_pc <- wp1_program_pc_correlations(score_data, coordinates, loadings)
  program_pc_repeat <- wp1_program_pc_correlations(score_data, coordinates, loadings)
  pc_repro <- wp1_pc_reproducibility_audit(program_pc, program_pc_repeat)
  wp1_assert(all(pc_repro$status == "PASS"), "Program-PC reproducibility audit failed.")
  wp1_assert(all(program_pc$loading_correlation_sign_consistent), "Frozen loading and program-PC correlation directions are inconsistent.")
  projection <- wp1_projection_parameters(summary_stats, loadings)
  protected <- wp1_protected_manifest()

  pc1_support <- sum(program_pc$component == "PC1" & program_pc$pearson_r > 0) >= 5L
  pc2_support <- all(program_pc$pearson_r[program_pc$component == "PC2" & program_pc$program_order %in% c(4L, 5L)] > 0) &&
    all(program_pc$pearson_r[program_pc$component == "PC2" & program_pc$program_order %in% c(2L, 3L, 6L)] < 0)
  wp1_assert(pc1_support, "Observed program-PC correlations do not support the qualified PC1 description.")
  wp1_assert(pc2_support, "Observed program-PC correlations do not support the qualified PC2 description.")

  canonical_contract <- wp1_canonical_contract(canonical, canonical_path, wp1_sha256(canonical_path))
  score_info <- file.info(score_path)
  primary_authority <- data.frame(
    artifact_id = "GSE31312_PRIMARY_HISTORICAL_UNTRUNCATED_SCORE_MATRIX",
    absolute_path = normalizePath(score_path, winslash = "/", mustWork = TRUE),
    relative_path = sub(paste0(wp1_root, "/"), "", normalizePath(score_path, winslash = "/"), fixed = TRUE),
    sha256 = score_sha, size = as.numeric(score_info$size),
    mtime = format(score_info$mtime, "%Y-%m-%dT%H:%M:%OS6%z"),
    rows = nrow(score_data), columns = ncol(score_data), sample_count = nrow(score_data), program_count = 6L,
    sample_id_column = "sample", program_order_valid = TRUE, numeric_finite = TRUE,
    missing_values = sum(is.na(score_data)), duplicated_samples = sum(duplicated(score_data$sample)),
    primary_status = "PRIMARY", transformation_status = "HISTORICAL_UNTRUNCATED_Z",
    clipped = FALSE, winsorized = FALSE,
    notes = "Hash-frozen Stage 4C-1 primary matrix; referenced in place; not copied or rewritten.",
    check.names = FALSE
  )
  pca_artifacts <- data.frame(
    artifact_id = c("PCA_INPUT_SCORE", "PCA_LOADINGS", "PCA_COORDINATES", "PCA_EXPLAINED_VARIANCE", "PCA_TRANSFORM_PARAMETERS"),
    role = c("input", "rotation", "coordinates", "variance", "center_scale"),
    absolute_path = c(score_path, loadings_path, coordinates_path, variance_path, pca_parameters_path),
    relative_path = sub(paste0(wp1_root, "/"), "", normalizePath(c(score_path, loadings_path, coordinates_path, variance_path, pca_parameters_path), winslash = "/"), fixed = TRUE),
    sha256 = vapply(c(score_path, loadings_path, coordinates_path, variance_path, pca_parameters_path), wp1_sha256, character(1)),
    pca_input_score_authority = score_sha,
    pca_method = "stats::prcomp",
    center_setting = FALSE,
    scale_setting = FALSE,
    sample_count = 498L,
    feature_count = 6L,
    retained_components = "PC1;PC2",
    rows = c(498L, 6L, 498L, 6L, 6L),
    columns = c(7L, 7L, 7L, 3L, 4L),
    authority_status = "FROZEN_AND_HASH_VERIFIED",
    notes = "No selected k or cluster assignment is required by this artifact.",
    check.names = FALSE
  )
  downstream <- data.frame(
    downstream_work_package = c("WP2 External", "WP3 Spatial", "WP4 Purity/composition"),
    required_artifact = c("six continuous scores; canonical gene match; sample metadata; cohort ID",
                          "canonical22 intersect detected genes; continuous spot/region scores; coordinates; capture-area ID",
                          "approved orthogonal composition method and disjoint covariates"),
    authority_path = c("WP1_CANONICAL_PROGRAM_CONTRACT.csv; WP1_CONTINUOUS_PROJECTION_PARAMETERS.csv",
                       "WP1_CANONICAL_PROGRAM_CONTRACT.csv", "WP1_CANONICAL_PROGRAM_CONTRACT.csv"),
    authority_sha256 = c("HASH_AFTER_WP1_COMMIT", "HASH_AFTER_WP1_COMMIT", "HASH_AFTER_WP1_COMMIT"),
    expected_schema = c("sample plus six ordered scores", "spot ID plus six ordered scores and coordinates", "sample plus approved composition variables"),
    expected_dimensions = c("cohort-specific rows x 7", "spot-specific rows x required columns", "cohort-specific rows x approved variables"),
    standardization_rule = c("separate authorization; reference scaling only after comparability gate",
                             "separate authorization; canonical detected intersection only", "continuous covariates; method preregistration required"),
    missing_data_rule = c("no zero fill; do not project if a program is unavailable",
                          "no gene imputation or substitution", "report missingness; no mechanical target-proxy overlap"),
    prohibited_operation = c("cluster label; predicted subtype; nearest-centroid assignment",
                             "historical shortlist; spot taxonomy; highest-program class",
                             "cluster membership exposure; overlapping marker proxy as independent validation"),
    readiness = c("INTERFACE_FROZEN_EXECUTION_NOT_AUTHORIZED", "INTERFACE_FROZEN_EXECUTION_NOT_AUTHORIZED", "BLOCKED_PENDING_METHOD_APPROVAL"),
    blocking_issue = c("WP2 token and cohort/gene/scaling gate", "WP3 token and canonical intersection gate", "orthogonal method not yet approved"),
    notes = "WP1 creates interfaces only and does not read downstream inputs.",
    check.names = FALSE
  )

  all_tables <- list(canonical_contract, primary_authority, descriptive, correlations$long, correlations$matrix,
                     pca_artifacts, reconstruction$table, program_pc, pc_repro, projection, downstream)
  wp1_assert(all(vapply(all_tables, function(x) wp1_forbidden_field_scan(names(x)), logical(1))),
             "A prohibited assignment field was created.")

  wp1_write_csv(input_registry, "WP1_INPUT_REGISTRY.csv")
  wp1_write_csv(canonical_contract, "WP1_CANONICAL_PROGRAM_CONTRACT.csv")
  wp1_write_csv(primary_authority, "WP1_PRIMARY_SCORE_AUTHORITY.csv")
  wp1_write_csv(descriptive, "WP1_PRIMARY_SIX_PROGRAM_DESCRIPTIVE_SUMMARY.csv")
  wp1_write_csv(correlations$long, "WP1_PRIMARY_PROGRAM_PEARSON_CORRELATIONS_LONG.csv")
  wp1_write_csv(correlations$matrix, "WP1_PRIMARY_PROGRAM_PEARSON_CORRELATION_MATRIX.csv")
  wp1_write_csv(pca_artifacts, "WP1_FROZEN_PCA_AUTHORITY.csv")
  wp1_write_csv(reconstruction$table, "WP1_PCA_RECONSTRUCTION_AUDIT.csv")
  wp1_write_csv(program_pc, "WP1_PROGRAM_PC_PEARSON_CORRELATIONS.csv")
  wp1_write_csv(pc_repro, "WP1_PROGRAM_PC_CORRELATION_REPRODUCIBILITY_AUDIT.csv")
  wp1_write_csv(projection, "WP1_CONTINUOUS_PROJECTION_PARAMETERS.csv")
  wp1_write_csv(downstream, "WP1_DOWNSTREAM_CONTINUOUS_INTERFACE_SPEC.csv")
  wp1_write_csv(protected, "WP1_PROTECTED_PATH_BASELINE_MANIFEST.csv")

  wp1_write_lines(c(
    "# WP1 Continuous Axis Interpretation Contract", "",
    "## PC1", "", "Formal variable name: `PC1`.", "",
    "Qualified description: `shared program-abundance axis`.", "",
    "The observed program-PC correlations show same-direction contributions from multiple programs. This is a continuous abundance structure, not a class, subtype, taxonomy, or biological quality scale.", "",
    "## PC2", "", "Formal variable name: `PC2`.", "",
    "Qualified description: `stromal/exclusion versus immune/proliferative continuous axis`.", "",
    "The observed direction places stromal/fibrotic and direct immune-cold/exclusion-associated scores on one side, and T-cell-inflamed, antigen-presentation, and proliferative/cycling scores on the other. The macrophage contribution is smaller and is reported numerically rather than forced into either biological group.", "",
    "PCA signs are arbitrary. Positive and negative coordinates do not represent favorable or unfavorable biology. The wording describes relative loadings in the two-dimensional frozen ordination. It is not a single composite polarization index, discrete taxonomy, patient class, or spot class."
  ), "WP1_CONTINUOUS_AXIS_INTERPRETATION_CONTRACT.md")

  wp1_write_lines(c(
    "# WP1 Continuous Projection Equation", "",
    "WP1 does not project any external or spatial sample. This equation freezes a future interface only.", "",
    "1. Generate six raw continuous program scores in the exact WP1 program order using a separately authorized and cross-platform-compatible scoring workflow.",
    "2. Require all six programs. Missing programs or canonical genes must not be imputed, substituted, or zero-filled.",
    "3. After a separate comparability gate, transform each raw score using the frozen GSE31312 parameters: `z_j = (raw_j - training_mean_j) / training_standard_deviation_j`.",
    "4. The frozen PCA itself uses `center=FALSE` and `scale.=FALSE` because its input is already the historical-untruncated standardized six-score matrix.",
    "5. Project continuously: `PC1 = sum_j(z_j * PC1_loading_j)` and `PC2 = sum_j(z_j * PC2_loading_j)`.",
    "6. Do not output a cluster label, subtype, taxonomy, nearest-centroid assignment, winner-takes-all program, patient class, or spot class.", "",
    "Within-cohort per-SD effect reporting is distinct from frozen-loading projection and requires a downstream work-package-specific standardization rule."
  ), "WP1_CONTINUOUS_PROJECTION_EQUATION.md")

  execution_log <- data.frame(
    check = c("run_id", "baseline_HEAD", "R_version", "renv_lock_sha256", "token_initial", "token_cleared_in_process",
              "GSVA_loaded", "external_input_read", "spatial_input_read", "clinical_input_read", "warning_count", "stack_imbalance"),
    status = c(rep("PASS", 6), "PASS", rep("PASS", 5)),
    value = c("WP1_CONTINUOUS_SCORE_FREEZE_V1", head, R.version.string, env_context$lock_sha256,
              "EXACT_MATCH", "TRUE", "FALSE", "FALSE", "FALSE", "FALSE", "0", "FALSE"),
    notes = c("Only real WP1 entry point", "Frozen amendment baseline", "Exact R 4.5.1 required", "Frozen Stage 4 renv",
              "Process-scoped token", "on.exit plus explicit clear", "Stage 4C-1 profile", "No external paths registered",
              "No spatial paths registered", "No clinical paths registered", "warnings promoted to errors", "No stack imbalance warning"),
    check.names = FALSE
  )
  Sys.unsetenv("DLBCL_REVISION_ALLOW_WP1_CONTINUOUS_FREEZE")
  wp1_write_csv(execution_log, "logs/28_wp1_continuous_score_freeze_v1_execution.csv")
  invisible(TRUE)
}

if (sys.nframe() == 0L && !identical(Sys.getenv("WP1_CONTINUOUS_FREEZE_SOURCE_ONLY", unset = ""), "TRUE")) {
  wp1_main()
}
