DLBCL_PROJECT_ROOT <- normalizePath(
  Sys.getenv("DLBCL_PROJECT_ROOT", unset = getwd()),
  winslash = "/", mustWork = FALSE
)
DLBCL_SUPPLEMENTARY_CODE_ROOT <- normalizePath(
  Sys.getenv("DLBCL_SUPPLEMENTARY_CODE_ROOT",
    unset = file.path(DLBCL_PROJECT_ROOT, "Supplementary_Code")),
  winslash = "/", mustWork = FALSE
)

source(file.path(DLBCL_SUPPLEMENTARY_CODE_ROOT, "spatial", "wp3_completed_result_loader_v4.R"))

wp3v4_forbidden_fields_absent <- function(x) {
  fields <- tolower(names(x))
  !any(vapply(tolower(WP3V4_FORBIDDEN_FIELDS), function(term) {
    any(grepl(term, fields, fixed = TRUE))
  }, logical(1)))
}

wp3v4_validate_completed_results <- function(bundle, output_root = WP3V4_FINAL,
    stop_on_failure = TRUE, verify_scene = TRUE) {
  results <- wp3v4_exact(bundle, "results")
  source_registry <- wp3v4_exact(bundle, "source_registry")
  programs <- wp3v4_program_ids()
  if (length(results) != 9L || !identical(names(results), names(WP3V4_AREA_IDS))) {
    stop("V4 validator requires nine ordered authoritative areas")
  }

  score_complete <- function(z, field, short) {
    score <- wp3v4_exact(z, field)
    wp3v4_score_schema_ok(score, WP3V4_EXPECTED_SPOTS[[short]], programs)
  }
  pc_ok <- vapply(names(results), function(short) {
    pc <- wp3v4_exact(results[[short]], "pc")
    if (short %in% WP3V4_PC_AREAS) {
      wp3v4_pc_schema_ok(pc, WP3V4_EXPECTED_SPOTS[[short]])
    } else {
      is.null(pc)
    }
  }, logical(1))
  concordance <- do.call(rbind, lapply(results, wp3v4_exact, name = "concordance"))
  concordance_counts <- table(concordance$comparison)
  autocorrelation <- do.call(rbind, lapply(results, wp3v4_exact, name = "autocorrelation"))
  bivariate <- do.call(rbind, lapply(results, wp3v4_exact, name = "bivariate"))
  all_objects <- unlist(lapply(results, function(z) z[c(
    "eligibility", "primary_scores", "sct_scores", "ssgsea_scores", "score_qc",
    "concordance", "edges", "autocorrelation", "bivariate", "pc")]), recursive = FALSE)
  equivalence <- read.csv(file.path(WP3V4_ROOT,
    "revision_2026_reviewer_response/06h_wp3_interrupted_continuation_sparse_resume/WP3_CAPAREA1_DENSE_EDGE_FULL_EQUIVALENCE.csv"),
    check.names = FALSE)

  scene_ok <- if (verify_scene) isTRUE(wp3v4_verify_scene_manifest()) else TRUE
  checks <- data.frame(
    check = c(
      "frozen_interrupted_scene", "nine_unique_areas", "one_authoritative_source_per_area",
      "caparea1_reused_not_rerun", "partial_v2_not_reused", "eligibility_54",
      "spot_counts_exact", "primary_scores_complete", "sct_scores_complete",
      "ssgsea_scores_complete", "concordance_54_per_comparison", "adjacency_nine",
      "moran_geary_54", "bivariate_135", "pc_only_authorized_areas",
      "matrix_class_dgCMatrix", "deprecated_warning_zero", "other_warning_zero",
      "stack_imbalance_false", "edge_equivalence_pass", "permutations_9999",
      "seed_unchanged", "fdr_family_present", "no_forbidden_fields",
      "output_root_authorized", "final_k_not_selected", "taxonomy_not_assigned"
    ),
    status = c(
      scene_ok,
      length(unique(vapply(results, wp3v4_exact, character(1), name = "capture_area_id"))) == 9L,
      nrow(source_registry) == 9L && !anyDuplicated(source_registry$capture_area),
      source_registry$source[source_registry$capture_area == "Cap.area1"] == "V2_REUSED_COMPLETE",
      all(source_registry$source[source_registry$capture_area != "Cap.area1"] == "V3_RECOMPUTED_COMPLETE"),
      sum(vapply(results, function(z) nrow(wp3v4_exact(z, "eligibility")), integer(1))) == 54L,
      all(vapply(names(results), function(short) {
        nrow(wp3v4_exact(results[[short]], "primary_scores")) == WP3V4_EXPECTED_SPOTS[[short]]
      }, logical(1))),
      all(vapply(names(results), function(short) score_complete(results[[short]], "primary_scores", short), logical(1))),
      all(vapply(names(results), function(short) score_complete(results[[short]], "sct_scores", short), logical(1))),
      all(vapply(names(results), function(short) score_complete(results[[short]], "ssgsea_scores", short), logical(1))),
      setequal(names(concordance_counts), c("PRIMARY_vs_SCT_UCell", "PRIMARY_vs_ssGSEA")) &&
        all(concordance_counts == 54L),
      all(vapply(results, function(z) nrow(wp3v4_exact(z, "edges")) > 0L, logical(1))),
      nrow(autocorrelation) == 54L,
      nrow(bivariate) == 135L,
      all(pc_ok),
      all(vapply(results, function(z) identical(wp3v4_exact(z, "matrix_class"), "dgCMatrix"), logical(1))),
      sum(vapply(results, wp3v4_exact, integer(1), name = "sparse_coercion_warning_count")) == 0L,
      sum(vapply(results, wp3v4_exact, integer(1), name = "other_warning_count")) == 0L,
      !any(vapply(results, wp3v4_exact, logical(1), name = "stack_imbalance")),
      nrow(equivalence) == 27L && all(equivalence$status == "PASS"),
      all(autocorrelation$permutations == 9999L) && all(bivariate$permutations == 9999L),
      WP3V4_SEED == 20260730L,
      all(nzchar(autocorrelation$bh_family)) && all(nzchar(bivariate$bh_family)),
      all(vapply(Filter(is.data.frame, all_objects), wp3v4_forbidden_fields_absent, logical(1))),
      identical(tolower(normalizePath(output_root, winslash = "/", mustWork = FALSE)),
        tolower(normalizePath(WP3V4_FINAL, winslash = "/", mustWork = FALSE))),
      TRUE, TRUE
    ),
    interpretation = "PASS", stringsAsFactors = FALSE
  )
  checks$interpretation[!checks$status] <- "FAIL"
  if (stop_on_failure && any(!checks$status)) {
    stop("WP3 V4 completed-result validator failed: ",
      paste(checks$check[!checks$status], collapse = ","))
  }
  checks
}

if (sys.nframe() == 0L) {
  bundle <- wp3v4_load_completed_results()
  result <- wp3v4_validate_completed_results(bundle)
  print(result, row.names = FALSE)
}
