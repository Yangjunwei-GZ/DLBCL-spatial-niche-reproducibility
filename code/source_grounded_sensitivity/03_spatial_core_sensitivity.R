DLBCL_PROJECT_ROOT <- normalizePath(
  Sys.getenv("DLBCL_PROJECT_ROOT", unset = getwd()),
  winslash = "/", mustWork = FALSE
)
DLBCL_SUPPLEMENTARY_CODE_ROOT <- normalizePath(
  Sys.getenv("DLBCL_SUPPLEMENTARY_CODE_ROOT",
    unset = file.path(DLBCL_PROJECT_ROOT, "Supplementary_Code")),
  winslash = "/", mustWork = FALSE
)

source(file.path(DLBCL_SUPPLEMENTARY_CODE_ROOT, "source_grounded_sensitivity", "01_common.R"))
o6_require_token(); o6_set_library()
o6_require_packages(c("Matrix", "SeuratObject", "Seurat", "UCell", "BiocParallel", "digest"))

sets <- o6_gene_sets(); areas <- o6_area_roles(); paths <- o6_authority_paths()
seed_registry <- utils::read.csv(file.path(O6_AMENDMENT, "06O_MASTER_SEED_REGISTRY.csv"), check.names = FALSE)
fdr_contract <- utils::read.csv(file.path(O6_AMENDMENT, "06O_FDR_FAMILY_CONTRACT.csv"), check.names = FALSE)
full_stats <- utils::read.csv(paths$moran_geary, check.names = FALSE)

gate_rows <- list(); detection_rows <- list(); concordance_rows <- list()
stat_rows <- list(); permutation_rows <- list()

for (area_index in seq_len(nrow(areas))) {
  area <- areas$capture_area_id[area_index]
  role <- areas$analysis_role[area_index]
  message("Spatial area ", area_index, "/", nrow(areas), ": ", area)
  counts <- o6_read_spatial_counts(area)
  full_detected <- lapply(sets$full, intersect, y = rownames(counts))
  full_scores <- o6_ucell(counts, full_detected, area)

  score_obj <- paths$registry[paths$registry$capture_area_id == area &
    paths$registry$object_type == "primary_program_score_matrix", , drop = FALSE]
  edge_obj <- paths$registry[paths$registry$capture_area_id == area &
    paths$registry$object_type == "adjacency_edge_list", , drop = FALSE]
  if (nrow(score_obj) != 1L || nrow(edge_obj) != 1L) stop("Spatial authority registry ambiguity: ", area)
  authority <- utils::read.csv(score_obj$path, check.names = FALSE)
  expected_programs <- unname(O6_PROGRAM_IDS)
  area_id_match <- identical(area, score_obj$capture_area_id)
  barcode_set_match <- setequal(full_scores$barcode, authority$barcode)
  barcode_order_match <- identical(full_scores$barcode, authority$barcode)
  program_names_match <- identical(names(full_scores)[-1L], names(authority)[-1L]) &&
    identical(names(full_scores)[-1L], expected_programs)
  dimensions_match <- identical(dim(full_scores), dim(authority))
  na_positions_match <- dimensions_match && identical(is.na(as.matrix(full_scores[, -1L, drop = FALSE])),
    is.na(as.matrix(authority[, -1L, drop = FALSE])))
  max_abs <- if (dimensions_match && program_names_match && barcode_order_match) {
    max(abs(as.matrix(full_scores[, -1L, drop = FALSE]) - as.matrix(authority[, -1L, drop = FALSE])), na.rm = TRUE)
  } else NA_real_
  gate_pass <- area_id_match && barcode_set_match && barcode_order_match && program_names_match &&
    dimensions_match && na_positions_match && is.finite(max_abs) && max_abs <= 1e-12
  gate_rows[[length(gate_rows) + 1L]] <- data.frame(capture_area_id = area, area_role = role,
    authority_path = normalizePath(score_obj$path, winslash = "/"), authority_sha256 = o6_sha256(score_obj$path),
    regenerated_spot_count = nrow(full_scores), authority_spot_count = nrow(authority),
    regenerated_program_count = ncol(full_scores) - 1L, authority_program_count = ncol(authority) - 1L,
    capture_area_identifier_exact = area_id_match, barcode_set_exact = barcode_set_match,
    barcode_order_exact = barcode_order_match, program_names_and_order_exact = program_names_match,
    dimensions_exact = dimensions_match, NA_positions_exact = na_positions_match,
    maximum_absolute_difference = max_abs, tolerance = 1e-12,
    authority_gate_status = if (gate_pass) "PASS" else "HARD_STOP", stringsAsFactors = FALSE)

  core_detected <- lapply(sets$core, intersect, y = rownames(counts))
  for (program_id in unname(O6_PROGRAM_IDS)) {
    core_genes <- if (program_id %in% names(sets$core)) sets$core[[program_id]] else character()
    detected <- if (program_id %in% names(core_detected)) core_detected[[program_id]] else character()
    seed_exists <- any(seed_registry$capture_area_id == area & seed_registry$program == o6_program_name(program_id))
    eligible <- gate_pass && length(core_genes) > 0L && length(detected) >= 5L && seed_exists
    reason <- if (!gate_pass) "FULL_REGENERATION_AUTHORITY_HARD_STOP" else if (!length(core_genes)) {
      "IMMUNE_COLD_HAS_NO_DIRECT_SOURCE_CORE"
    } else if (length(detected) < 5L) "FEWER_THAN_FIVE_DETECTED_CORE_GENES" else if (!seed_exists) {
      "NOT_IN_FROZEN_39_COMBINATION_REGISTRY"
    } else ""
    detection_rows[[length(detection_rows) + 1L]] <- data.frame(capture_area_id = area, area_role = role,
      program = o6_program_name(program_id), program_id = program_id,
      full_canonical_gene_count = length(sets$full[[program_id]]), full_detected_gene_count = length(full_detected[[program_id]]),
      full_detected_genes = paste(full_detected[[program_id]], collapse = ";"),
      full_missing_genes = paste(setdiff(sets$full[[program_id]], full_detected[[program_id]]), collapse = ";"),
      core_gene_count = length(core_genes), core_detected_gene_count = length(detected),
      core_detection_fraction = if (length(core_genes)) length(detected) / length(core_genes) else NA_real_,
      core_detected_genes = paste(detected, collapse = ";"),
      core_missing_genes = paste(setdiff(core_genes, detected), collapse = ";"),
      substitution_performed = FALSE, core_eligibility = if (eligible) "EVALUABLE" else "NOT_EVALUABLE",
      exclusion_reason = reason, stringsAsFactors = FALSE)
  }

  if (!gate_pass) {
    rm(counts, full_scores, authority); gc(verbose = FALSE)
    next
  }

  eligible_programs <- unname(O6_PROGRAM_IDS)[vapply(unname(O6_PROGRAM_IDS), function(id) {
    id %in% names(core_detected) && length(core_detected[[id]]) >= 5L &&
      any(seed_registry$capture_area_id == area & seed_registry$program == o6_program_name(id))
  }, logical(1))]
  core_scores <- o6_ucell(counts, core_detected[eligible_programs], area)
  edges <- utils::read.csv(edge_obj$path, check.names = FALSE)
  edge <- o6_edge_index(full_scores$barcode, edges)

  for (program_id in unname(O6_PROGRAM_IDS)) {
    evaluable <- program_id %in% eligible_programs
    full <- full_scores[[program_id]]
    core <- if (evaluable) core_scores[[program_id]] else rep(NA_real_, length(full))
    concordance_rows[[length(concordance_rows) + 1L]] <- data.frame(capture_area_id = area,
      area_role = role, program = o6_program_name(program_id), program_id = program_id,
      evaluation_status = if (evaluable) "EVALUABLE" else "NOT_EVALUABLE",
      spot_count = length(full), full_detected_gene_count = length(full_detected[[program_id]]),
      core_detected_gene_count = if (program_id %in% names(core_detected)) length(core_detected[[program_id]]) else 0L,
      full_mean = mean(full), full_sd = stats::sd(full),
      core_mean = if (evaluable) mean(core) else NA_real_, core_sd = if (evaluable) stats::sd(core) else NA_real_,
      full_zero_variance = stats::sd(full) == 0,
      core_zero_variance = if (evaluable) stats::sd(core) == 0 else NA,
      Pearson = if (evaluable) o6_safe_cor(full, core, "pearson") else NA_real_,
      Spearman = if (evaluable) o6_safe_cor(full, core, "spearman") else NA_real_,
      missingness_status = if (evaluable && !anyNA(c(full, core))) "COMPLETE" else if (evaluable) "MISSING_VALUES" else "NOT_APPLICABLE",
      stringsAsFactors = FALSE)
    if (!evaluable) next

    full_row <- full_stats[full_stats$capture_area_id == area & full_stats$program_id == program_id, , drop = FALSE]
    if (nrow(full_row) != 1L) stop("Missing final full-program spatial comparator: ", area, "/", program_id)
    endpoint_results <- list()
    for (endpoint in c("CORE_MORAN", "CORE_GEARY")) {
      analysis_type <- paste0(endpoint, "_PERMUTATION_ENDPOINT")
      sr <- seed_registry[seed_registry$capture_area_id == area & seed_registry$program == o6_program_name(program_id) &
        seed_registry$analysis_type == analysis_type, , drop = FALSE]
      if (nrow(sr) != 1L) stop("Missing unique endpoint seed: ", area, "/", program_id, "/", endpoint)
      result <- o6_permutation(core, edge, endpoint, as.integer(sr$derived_seed), 9999L)
      if (result$n != 9999L || !result$all_finite) stop("Incomplete/nonfinite permutation endpoint")
      endpoint_results[[endpoint]] <- result
      permutation_rows[[length(permutation_rows) + 1L]] <- data.frame(
        capture_area_id = area, area_role = role, program = o6_program_name(program_id), program_id = program_id,
        endpoint = endpoint, master_seed = O6_MASTER_SEED, derived_seed = as.integer(sr$derived_seed),
        requested_permutations = 9999L, completed_permutations = result$n,
        observed_statistic = result$observed, null_mean = result$null_mean, null_sd = result$null_sd,
        null_minimum = result$null_min, null_maximum = result$null_max, all_null_statistics_finite = result$all_finite,
        empirical_p_value = result$p, permutation_scheme = "SCORE_LABEL_PERMUTATION_ON_FIXED_FINAL_EDGE_GRAPH",
        execution_status = "COMPLETED", stringsAsFactors = FALSE)
    }
    cm <- endpoint_results$CORE_MORAN$observed; cg <- endpoint_results$CORE_GEARY$observed
    fm <- as.numeric(full_row$Moran_I); fg <- as.numeric(full_row$Geary_C)
    moran_ratio <- if (is.finite(fm) && fm > 0) cm / fm else NA_real_
    geary_den <- 1 - fg
    geary_ratio <- if (is.finite(geary_den) && geary_den > 0) (1 - cg) / geary_den else NA_real_
    stat_rows[[length(stat_rows) + 1L]] <- data.frame(capture_area_id = area, area_role = role,
      program = o6_program_name(program_id), program_id = program_id, spot_count = length(core),
      edge_count = nrow(edges), core_gene_count = length(sets$core[[program_id]]),
      detected_core_gene_count = length(core_detected[[program_id]]),
      full_Moran_I = fm, core_Moran_I = cm, core_minus_full_Moran = cm - fm,
      Moran_effect_retention_ratio = moran_ratio,
      Moran_retention_status = if (is.finite(moran_ratio)) "APPLICABLE" else "NOT_APPLICABLE",
      Moran_direction_concordant = sign(cm) == sign(fm), Moran_empirical_p = endpoint_results$CORE_MORAN$p,
      full_Geary_C = fg, core_Geary_C = cg, core_minus_full_Geary = cg - fg,
      Geary_spatial_departure_retention = geary_ratio,
      Geary_retention_status = if (is.finite(geary_ratio)) "APPLICABLE" else "NOT_APPLICABLE",
      Geary_departure_direction_concordant = sign(1 - cg) == sign(1 - fg),
      Geary_empirical_p = endpoint_results$CORE_GEARY$p,
      permutations_per_endpoint = 9999L, evaluation_status = "EVALUABLE", stringsAsFactors = FALSE)
  }
  rm(counts, full_scores, core_scores, authority, edges, edge); gc(verbose = FALSE)
}

gate <- do.call(rbind, gate_rows); detection <- do.call(rbind, detection_rows)
concordance <- do.call(rbind, concordance_rows); stats_out <- do.call(rbind, stat_rows)
permutations <- do.call(rbind, permutation_rows)

fdr_test_rows <- list(); fdr_family_rows <- list()
for (i in seq_len(nrow(fdr_contract))) {
  fc <- fdr_contract[i, ]
  if (fc$role_family == "EXPLORATORY_ANTIGEN") {
    fdr_family_rows[[length(fdr_family_rows) + 1L]] <- data.frame(row_type = "FAMILY_SUMMARY",
      family_id = fc$family_id, endpoint = fc$endpoint, role_family = fc$role_family,
      family_label = fc$family_label, capture_area_id = "", program = "", program_id = "",
      raw_p_value = NA_real_, BH_FDR = NA_real_, expected_test_count = as.integer(fc$expected_test_count),
      observed_test_count = 0L, family_status = "EMPTY_NOT_APPLICABLE", stringsAsFactors = FALSE)
    next
  }
  p_col <- if (fc$endpoint == "CORE_MORAN") "Moran_empirical_p" else "Geary_empirical_p"
  idx <- which(stats_out$area_role == fc$role_family)
  adjusted <- stats::p.adjust(stats_out[[p_col]][idx], method = "BH")
  if (length(idx) != as.integer(fc$expected_test_count)) stop("FDR family count mismatch: ", fc$family_id)
  if (fc$endpoint == "CORE_MORAN") stats_out$Moran_BH_FDR[idx] <- adjusted else stats_out$Geary_BH_FDR[idx] <- adjusted
  fdr_test_rows[[length(fdr_test_rows) + 1L]] <- data.frame(row_type = "TEST",
    family_id = fc$family_id, endpoint = fc$endpoint, role_family = fc$role_family,
    family_label = fc$family_label, capture_area_id = stats_out$capture_area_id[idx],
    program = stats_out$program[idx], program_id = stats_out$program_id[idx],
    raw_p_value = stats_out[[p_col]][idx], BH_FDR = adjusted,
    expected_test_count = as.integer(fc$expected_test_count), observed_test_count = length(idx),
    family_status = "COMPLETED", stringsAsFactors = FALSE)
  fdr_family_rows[[length(fdr_family_rows) + 1L]] <- data.frame(row_type = "FAMILY_SUMMARY",
    family_id = fc$family_id, endpoint = fc$endpoint, role_family = fc$role_family,
    family_label = fc$family_label, capture_area_id = "", program = "", program_id = "",
    raw_p_value = NA_real_, BH_FDR = NA_real_, expected_test_count = as.integer(fc$expected_test_count),
    observed_test_count = length(idx), family_status = "COMPLETED", stringsAsFactors = FALSE)
}
fdr_out <- do.call(rbind, c(fdr_test_rows, fdr_family_rows))

o6_write_csv_once(gate, file.path(O6_OUTPUTS, "spatial/FULL_REGENERATION_AUTHORITY_CHECK.csv"))
o6_write_csv_once(detection, file.path(O6_OUTPUTS, "spatial/CORE_GENE_DETECTION_BY_AREA.csv"))
o6_write_csv_once(concordance, file.path(O6_OUTPUTS, "spatial/FULL_CORE_SPOT_SCORE_CONCORDANCE.csv"))
o6_write_csv_once(stats_out, file.path(O6_OUTPUTS, "spatial/CORE_MORAN_GEARY.csv"))
o6_write_csv_once(permutations, file.path(O6_OUTPUTS, "spatial/CORE_SCORE_LABEL_PERMUTATIONS.csv"))
o6_write_csv_once(fdr_out, file.path(O6_OUTPUTS, "spatial/06O_SENSITIVITY_FDR.csv"))
cat("SPATIAL_EXECUTION=PASS\n")
