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
source(file.path(DLBCL_SUPPLEMENTARY_CODE_ROOT, "spatial", "wp3_common_v3.R"))
for (script in c("36a_wp3_load_and_validate_spatial_inputs_v3.R",
  "36b_wp3_program_area_eligibility_v3.R", "36c_wp3_primary_lognormalize_ucell_v3.R",
  "36d_wp3_sct_ucell_sensitivity_v3.R", "36e_wp3_ssgsea_sensitivity_v3.R",
  "36f_wp3_score_qc_concordance_v3.R", "36g_wp3_visium_neighbors_v3.R",
  "36h_wp3_spatial_autocorrelation_edge_v3.R",
  "36i_wp3_bivariate_spatial_association_edge_v3.R", "36j_wp3_exploratory_pc_v3.R",
  "36k_wp3_capture_area_synthesis_v3.R", "36v_wp3_resume_validator_v3.R")) {
  source(file.path(DLBCL_SUPPLEMENTARY_CODE_ROOT, "spatial", script))
}

wp3_single_file <- function(directory, pattern) {
  paths <- list.files(directory, pattern = pattern, full.names = TRUE)
  if (length(paths) != 1L) stop("Expected one frozen file for pattern: ", pattern)
  paths[[1L]]
}

wp3_load_v2_caparea1 <- function() {
  root <- WP3_INTERRUPTED_OUTPUT
  read_part <- function(dir, pattern = "[.]csv$") read.csv(
    wp3_single_file(file.path(root, dir), pattern), check.names = FALSE, stringsAsFactors = FALSE)
  input_qc <- read_part("01_input_qc", "Cap[.]area1.*input_qc[.]csv$")
  list(capture_area_id = input_qc$capture_area_id[[1L]], area_role = "CONTEXT_ONLY",
    input_qc = input_qc,
    eligibility = read_part("02_eligibility", "Cap[.]area1.*eligibility[.]csv$"),
    primary_scores = read_part("03_primary_ucell", "Cap[.]area1.*scores[.]csv$"),
    sct_scores = read_part("04_sct_ucell", "Cap[.]area1.*scores[.]csv$"),
    ssgsea_scores = read_part("05_ssgsea", "Cap[.]area1.*scores[.]csv$"),
    score_qc = read_part("06_score_qc", "Cap[.]area1.*qc[.]csv$"),
    concordance = read_part("07_method_concordance", "Cap[.]area1.*concordance[.]csv$"),
    edges = read_part("08_adjacency", "Cap[.]area1.*edges[.]csv$"),
    autocorrelation = read_part("09_moran", "Cap[.]area1.*autocorrelation[.]csv$"),
    bivariate = read_part("11_bivariate_moran", "Cap[.]area1.*bivariate[.]csv$"),
    pc = read_part("12_exploratory_pc", "Cap[.]area1.*pc[.]csv$"),
    matrix_class = "dgCMatrix", sparse_coercion_warning_count = 0L,
    other_warning_count = 0L, stack_imbalance = FALSE)
}

wp3_write_v3_area_stage <- function(state, stage) {
  wp3_write_csv_once(state$input_qc, file.path(stage, "01_input_qc.csv"))
  wp3_write_csv_once(state$eligibility, file.path(stage, "02_eligibility.csv"))
  wp3_write_csv_once(state$primary_scores, file.path(stage, "03_primary_ucell.csv"))
  wp3_write_csv_once(state$sct_scores, file.path(stage, "04_sct_ucell.csv"))
  wp3_write_csv_once(state$ssgsea_scores, file.path(stage, "05_ssgsea.csv"))
  wp3_write_csv_once(state$score_qc, file.path(stage, "06_score_qc.csv"))
  wp3_write_csv_once(state$concordance, file.path(stage, "07_method_concordance.csv"))
  wp3_write_csv_once(state$edges, file.path(stage, "08_adjacency.csv"))
  wp3_write_csv_once(state$autocorrelation, file.path(stage, "09_moran_geary.csv"))
  wp3_write_csv_once(state$bivariate, file.path(stage, "10_bivariate_moran.csv"))
  if (!is.null(state$pc)) wp3_write_csv_once(state$pc, file.path(stage, "11_exploratory_pc.csv"))
  expected <- if (is.null(state$pc)) 10L else 11L
  files <- list.files(stage, full.names = TRUE)
  if (length(files) != expected || any(file.info(files)$size <= 0)) stop("Area staging output is incomplete")
  marker <- data.frame(capture_area_id = state$capture_area_id, status = "COMPLETE",
    completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE), nodes = "36a-36j",
    permutations = 9999L, engine = "UNIQUE_UNDIRECTED_EDGE_LIST", stringsAsFactors = FALSE)
  wp3_write_csv_once(marker, file.path(stage, "AREA_COMPLETE_V3.csv"))
  invisible(TRUE)
}

wp3_load_v3_area <- function(area_dir) {
  read_part <- function(name) read.csv(file.path(area_dir, name), check.names = FALSE,
    stringsAsFactors = FALSE)
  input_qc <- read_part("01_input_qc.csv")
  pc_path <- file.path(area_dir, "11_exploratory_pc.csv")
  list(capture_area_id = input_qc$capture_area_id[[1L]],
    area_role = wp3_read_area_contract()$analysis_role[
      match(input_qc$capture_area_id[[1L]], wp3_read_area_contract()$capture_area_id)],
    input_qc = input_qc, eligibility = read_part("02_eligibility.csv"),
    primary_scores = read_part("03_primary_ucell.csv"), sct_scores = read_part("04_sct_ucell.csv"),
    ssgsea_scores = read_part("05_ssgsea.csv"), score_qc = read_part("06_score_qc.csv"),
    concordance = read_part("07_method_concordance.csv"), edges = read_part("08_adjacency.csv"),
    autocorrelation = read_part("09_moran_geary.csv"), bivariate = read_part("10_bivariate_moran.csv"),
    pc = if (file.exists(pc_path)) read_part("11_exploratory_pc.csv") else NULL,
    matrix_class = "dgCMatrix", sparse_coercion_warning_count = 0L,
    other_warning_count = 0L, stack_imbalance = FALSE)
}

wp3_resume_orchestrator_v3 <- function(output_root = WP3_FUTURE_OUTPUT) {
  wp3_require_r(); wp3_require_token("36")
  on.exit(Sys.unsetenv(WP3_TOKEN_NAME), add = TRUE)
  frozen <- wp3_validate_interrupted_scene()
  wp3_validate_complete_v2_caparea1()
  if (dir.exists(output_root) || file.exists(output_root)) stop("continuation_v3 create-once refusal")
  plan <- read.csv(file.path(WP3_CONTINUATION_CONTRACT_DIR, "WP3_RESUME_AREA_PLAN.csv"),
    check.names = FALSE, stringsAsFactors = FALSE)
  audit <- read.csv(file.path(WP3_CONTINUATION_CONTRACT_DIR,
    "WP3_CAPTURE_AREA_COMPLETION_AUDIT_AFTER_INTERRUPTION.csv"), check.names = FALSE)
  equivalence <- read.csv(file.path(WP3_CONTINUATION_CONTRACT_DIR,
    "WP3_CAPAREA1_DENSE_EDGE_FULL_EQUIVALENCE.csv"), check.names = FALSE)
  performance <- read.csv(file.path(WP3_CONTINUATION_CONTRACT_DIR,
    "WP3_EDGE_ENGINE_PERFORMANCE_BENCHMARK.csv"), check.names = FALSE)
  if (nrow(plan) != 9L || plan$capture_area[[1L]] != "Cap.area1" ||
      plan$old_outputs_reused[[1L]] != "TRUE" || any(plan$old_outputs_reused[-1L] != "FALSE")) {
    stop("Resume plan violates complete/partial reuse rules")
  }
  if (nrow(equivalence) != 27L || any(equivalence$status != "PASS") ||
      !all(performance$performance_gate == "PASS")) stop("Sparse equivalence/performance gate failed")
  runtime_packages <- wp3_require_runtime_packages(); wp3_mkdir_once(output_root)
  dir.create(file.path(output_root, "protocol"), recursive = TRUE)
  dir.create(file.path(output_root, "logs"), recursive = TRUE)
  dir.create(file.path(output_root, "checkpoints"), recursive = TRUE)
  dir.create(file.path(output_root, "staging"), recursive = TRUE)
  dir.create(file.path(output_root, "areas"), recursive = TRUE)
  protocol <- data.frame(protocol = "WP3_RESUME_V3", seed = WP3_SEED, permutations = 9999L,
    engine = "UNIQUE_UNDIRECTED_EDGE_LIST", completed_v2_areas_reused = 1L,
    partial_v2_outputs_reused = FALSE, final_k = "NOT_SELECTED", taxonomy = "NOT_ASSIGNED")
  wp3_write_csv_once(protocol, file.path(output_root, "protocol", "WP3_RESUME_PROTOCOL_V3.csv"))
  startup <- data.frame(R_version = paste(R.version$major, R.version$minor, sep = "."),
    seed = WP3_SEED, input_root = file.path(WP3_ROOT, "00_raw_data/GSE276542/standard_10x"),
    output_root = output_root, token_status = "VALID_NOT_RECORDED",
    started_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
  wp3_write_csv_once(startup, file.path(output_root, "logs", "WP3_RESUME_STARTUP_AUDIT_V3.csv"))
  wp3_write_csv_once(runtime_packages, file.path(output_root, "logs", "WP3_RESUME_RUNTIME_PACKAGE_RESOLUTION_V3.csv"))
  wp3_write_csv_once(frozen, file.path(output_root, "protocol", "WP3_FROZEN_V2_REFERENCE.csv"))
  wp3_write_csv_once(plan, file.path(output_root, "protocol", "WP3_RESUME_AREA_PLAN_RUNTIME.csv"))
  areas <- wp3_read_area_contract(); rerun <- plan$rerun_required == "TRUE"
  for (row in which(rerun)) {
    id <- plan$capture_area_id[[row]]; stage <- file.path(output_root, "staging", id)
    final <- file.path(output_root, "areas", id); started <- Sys.time()
    if (dir.exists(stage) || dir.exists(final)) stop("Area create-once refusal: ", id)
    dir.create(stage, recursive = TRUE)
    outcome <- tryCatch({
      state <- list(capture_area_id = id)
      state <- wp3_node_36a(state); state <- wp3_node_36b(state); state <- wp3_node_36c(state)
      state <- wp3_node_36d(state); state <- wp3_node_36e(state); state <- wp3_node_36f(state)
      state <- wp3_node_36g(state); state <- wp3_node_36h(state); state <- wp3_node_36i(state)
      state <- wp3_node_36j(state)
      wp3_write_v3_area_stage(state, stage)
      if (!file.rename(stage, final)) stop("Atomic area finalization failed: ", id)
      list(status = "COMPLETED", error = "")
    }, error = function(e) list(status = "FAILED", error = conditionMessage(e)))
    checkpoint <- data.frame(capture_area_id = id, status = outcome$status,
      started_at = format(started, tz = "UTC", usetz = TRUE),
      finished_at = format(Sys.time(), tz = "UTC", usetz = TRUE), error = outcome$error)
    wp3_write_csv_once(checkpoint, file.path(output_root, "checkpoints", paste0(id, "_STATUS_V3.csv")))
    if (outcome$status != "COMPLETED") stop("V3 resume stopped at ", id, ": ", outcome$error)
    rm(state); gc()
  }
  results <- list(wp3_load_v2_caparea1())
  for (id in plan$capture_area_id[rerun]) results[[length(results) + 1L]] <-
    wp3_load_v3_area(file.path(output_root, "areas", id))
  names(results) <- plan$capture_area
  source_registry <- data.frame(capture_area = plan$capture_area,
    capture_area_id = plan$capture_area_id,
    source = c("V2_REUSED_COMPLETE", rep("V3_RECOMPUTED", 8L)), stringsAsFactors = FALSE)
  aggregate <- wp3_node_36k(list(area_results = results, source_registry = source_registry))
  dir.create(file.path(output_root, "combined"), recursive = TRUE)
  wp3_write_csv_once(source_registry, file.path(output_root, "combined", "WP3_COMBINED_AREA_REGISTRY_V3.csv"))
  wp3_write_csv_once(do.call(rbind, lapply(results, `[[`, "autocorrelation")),
    file.path(output_root, "combined", "WP3_COMBINED_MORAN_GEARY_V3.csv"))
  wp3_write_csv_once(do.call(rbind, lapply(results, `[[`, "bivariate")),
    file.path(output_root, "combined", "WP3_COMBINED_BIVARIATE_V3.csv"))
  pc_registry <- do.call(rbind, lapply(seq_along(results), function(i) if (!is.null(results[[i]]$pc))
    data.frame(capture_area = names(results)[[i]], source = source_registry$source[[i]],
      rows = nrow(results[[i]]$pc), stringsAsFactors = FALSE) else NULL))
  wp3_write_csv_once(pc_registry, file.path(output_root, "combined", "WP3_COMBINED_PC_REGISTRY_V3.csv"))
  wp3_write_csv_once(aggregate$across_area_summary,
    file.path(output_root, "combined", "WP3_CAPTURE_AREA_SYNTHESIS_V3.csv"))
  aggregate$source_registry <- source_registry
  validator <- wp3_node_36v(aggregate, output_root)
  wp3_write_csv_once(validator, file.path(output_root, "WP3_RESUME_FINAL_VALIDATOR_V3.csv"))
  scope <- data.frame(scope = c("patient", "specimen", "region", "OS", "purity",
    "composition_adjustment", "malignant_assignment", "cluster", "class", "taxonomy"), executed = FALSE)
  wp3_write_csv_once(scope, file.path(output_root, "15_scope_isolation", "WP3_RESUME_SCOPE_ISOLATION_V3.csv"))
  completion <- data.frame(status = "COMPLETED", capture_areas = 9L, reused_v2_areas = 1L,
    recomputed_v3_areas = 8L, final_validator = all(validator$status),
    final_k = "NOT_SELECTED", taxonomy = "NOT_ASSIGNED")
  wp3_write_csv_once(completion, file.path(output_root, "WP3_RESUME_COMPLETION_REPORT_V3.csv"))
  invisible(aggregate)
}

if (sys.nframe() == 0L) wp3_resume_orchestrator_v3()
