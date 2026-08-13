DLBCL_PROJECT_ROOT <- normalizePath(
  Sys.getenv("DLBCL_PROJECT_ROOT", unset = getwd()),
  winslash = "/", mustWork = FALSE
)
DLBCL_SUPPLEMENTARY_CODE_ROOT <- normalizePath(
  Sys.getenv("DLBCL_SUPPLEMENTARY_CODE_ROOT",
    unset = file.path(DLBCL_PROJECT_ROOT, "Supplementary_Code")),
  winslash = "/", mustWork = FALSE
)

source(file.path(DLBCL_SUPPLEMENTARY_CODE_ROOT, "spatial", "wp3_common_v3.R"))

wp3_node_36v <- function(state, output_root = WP3_FUTURE_OUTPUT) {
  wp3_validate_interrupted_scene()
  wp3_validate_complete_v2_caparea1()
  results <- state$area_results
  if (is.null(results) || length(results) != 9L) stop("V3 validator requires nine authoritative areas")
  expected_spots <- c(Cap.area1 = 3384L, Cap.area2 = 3318L, Cap.area3 = 1197L,
    Cap.area4 = 1723L, Cap.area5 = 1827L, Cap.area6 = 4992L, Cap.area7 = 4992L,
    Cap.area8 = 4674L, Cap.area9 = 4923L)
  contract <- wp3_read_area_contract()
  area_short <- vapply(results, function(z) contract$capture_area_short[
    match(z$capture_area_id, contract$capture_area_id)], character(1))
  if (anyNA(area_short) || anyDuplicated(area_short)) stop("Area identity is incomplete or duplicated")
  objects <- unlist(lapply(results, function(z) z[c("eligibility", "primary_scores", "sct_scores",
    "ssgsea_scores", "score_qc", "concordance", "edges", "autocorrelation", "bivariate", "pc")]),
    recursive = FALSE)
  invisible(lapply(objects, function(x) if (is.data.frame(x)) wp3_forbidden_field_check(x)))
  score_complete <- function(z, field, short) {
    programs <- wp3_read_program_contract()$program_id
    expected_columns <- if (field == "ssgsea_scores") programs else paste0(programs, "_UCell")
    nrow(z[[field]]) == expected_spots[[short]] &&
      identical(names(z[[field]]), c("barcode", expected_columns)) &&
      all(is.finite(as.matrix(z[[field]][, -1L, drop = FALSE])))
  }
  pc_ok <- vapply(seq_along(results), function(i) {
    if (area_short[[i]] %in% c("Cap.area1", "Cap.area2", "Cap.area3")) !is.null(results[[i]]$pc)
    else is.null(results[[i]]$pc)
  }, logical(1))
  concordance <- do.call(rbind, lapply(results, `[[`, "concordance"))
  concordance_counts <- table(concordance$comparison)
  equivalence <- read.csv(file.path(WP3_CONTINUATION_CONTRACT_DIR,
    "WP3_CAPAREA1_DENSE_EDGE_FULL_EQUIVALENCE.csv"), check.names = FALSE)
  completion_audit <- read.csv(file.path(WP3_CONTINUATION_CONTRACT_DIR,
    "WP3_CAPTURE_AREA_COMPLETION_AUDIT_AFTER_INTERRUPTION.csv"), check.names = FALSE)
  source_registry <- state$source_registry
  checks <- data.frame(
    check = c("frozen_interrupted_scene", "nine_unique_areas", "one_authoritative_source_per_area",
      "caparea1_reused_not_rerun", "partial_v2_not_reused", "eligibility_54", "spot_counts_exact",
      "primary_scores_complete", "sct_scores_complete", "ssgsea_scores_complete",
      "concordance_54_per_comparison", "adjacency_nine", "moran_geary_54", "bivariate_135",
      "pc_only_authorized_areas", "matrix_class_dgCMatrix", "deprecated_warning_zero",
      "other_warning_zero", "stack_imbalance_false", "edge_equivalence_pass",
      "permutations_9999", "seed_unchanged", "fdr_family_present", "no_forbidden_fields",
      "output_root_authorized", "final_k_not_selected", "taxonomy_not_assigned"),
    status = c(TRUE, length(unique(area_short)) == 9L,
      nrow(source_registry) == 9L && !anyDuplicated(source_registry$capture_area),
      source_registry$source[source_registry$capture_area == "Cap.area1"] == "V2_REUSED_COMPLETE",
      all(source_registry$source[source_registry$capture_area != "Cap.area1"] == "V3_RECOMPUTED"),
      sum(vapply(results, function(z) nrow(z$eligibility), integer(1))) == 54L,
      all(vapply(seq_along(results), function(i) nrow(results[[i]]$primary_scores) == expected_spots[[area_short[[i]]]], logical(1))),
      all(vapply(seq_along(results), function(i) score_complete(results[[i]], "primary_scores", area_short[[i]]), logical(1))),
      all(vapply(seq_along(results), function(i) score_complete(results[[i]], "sct_scores", area_short[[i]]), logical(1))),
      all(vapply(seq_along(results), function(i) score_complete(results[[i]], "ssgsea_scores", area_short[[i]]), logical(1))),
      setequal(names(concordance_counts), c("PRIMARY_vs_SCT_UCell", "PRIMARY_vs_ssGSEA")) && all(concordance_counts == 54L),
      all(vapply(results, function(z) nrow(z$edges) > 0L, logical(1))),
      sum(vapply(results, function(z) nrow(z$autocorrelation), integer(1))) == 54L,
      sum(vapply(results, function(z) nrow(z$bivariate), integer(1))) == 135L,
      all(pc_ok), all(vapply(results, function(z) identical(z$matrix_class, "dgCMatrix"), logical(1))),
      sum(vapply(results, `[[`, integer(1), "sparse_coercion_warning_count")) == 0L,
      sum(vapply(results, `[[`, integer(1), "other_warning_count")) == 0L,
      !any(vapply(results, `[[`, logical(1), "stack_imbalance")),
      nrow(equivalence) == 27L && all(equivalence$status == "PASS"),
      all(do.call(rbind, lapply(results, `[[`, "autocorrelation"))$permutations == 9999L) &&
        all(do.call(rbind, lapply(results, `[[`, "bivariate"))$permutations == 9999L),
      WP3_SEED == 20260730L,
      all(nzchar(do.call(rbind, lapply(results, `[[`, "autocorrelation"))$bh_family)) &&
        all(nzchar(do.call(rbind, lapply(results, `[[`, "bivariate"))$bh_family)),
      TRUE,
      identical(tolower(normalizePath(output_root, winslash = "/", mustWork = FALSE)),
        tolower(normalizePath(WP3_FUTURE_OUTPUT, winslash = "/", mustWork = FALSE))),
      TRUE, TRUE), stringsAsFactors = FALSE)
  if (completion_audit$completion_status[completion_audit$capture_area == "Cap.area1"] != "COMPLETE") {
    checks$status[checks$check == "caparea1_reused_not_rerun"] <- FALSE
  }
  if (!all(checks$status)) stop("WP3 V3 resume validator failed: ",
    paste(checks$check[!checks$status], collapse = ","))
  checks
}

if (sys.nframe() == 0L) wp3_node_cli("36v", wp3_node_36v)
