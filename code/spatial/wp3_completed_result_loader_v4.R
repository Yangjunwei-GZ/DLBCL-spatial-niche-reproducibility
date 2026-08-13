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

WP3V4_ROOT <- DLBCL_PROJECT_ROOT
WP3V4_DIR <- file.path(WP3V4_ROOT,
  "revision_2026_reviewer_response/06i_wp3_final_validator_contract_fix")
WP3V4_SCENE <- file.path(WP3V4_ROOT,
  "revision_2026_reviewer_response/06c_wp3_real_spatial_continuous_analysis")
WP3V4_V2 <- file.path(WP3V4_SCENE, "continuation_v2")
WP3V4_V3 <- file.path(WP3V4_SCENE, "continuation_v3")
WP3V4_FINAL <- file.path(WP3V4_V3, "finalization_v2")
WP3V4_PROGRAM_CONTRACT <- file.path(WP3V4_ROOT,
  "revision_2026_reviewer_response/05x_wp1_continuous_score_freeze/WP1_CANONICAL_PROGRAM_CONTRACT.csv")
WP3V4_SCENE_MANIFEST <- file.path(WP3V4_DIR,
  "WP3_COMPLETED_PREVALIDATOR_SCENE_MANIFEST.csv")
WP3V4_SEED <- 20260730L

WP3V4_AREA_IDS <- c(
  Cap.area1 = "GSM8500534_Cap.area1_LN_V1",
  Cap.area2 = "GSM8500535_Cap.area2_tonsil_V1",
  Cap.area3 = "GSM8500536_Cap.area3_DLBCL_V1",
  Cap.area4 = "GSM8500537_Cap.area4_DLBCL_V2",
  Cap.area5 = "GSM8500538_Cap.area5_DLBCL_V2",
  Cap.area6 = "GSM8500539_Cap.area6_DLBCL_V2",
  Cap.area7 = "GSM8500540_Cap.area7_DLBCL_V2",
  Cap.area8 = "GSM8500541_Cap.area8_LN_V2",
  Cap.area9 = "GSM8500542_Cap.area9_tonsil_V2"
)
WP3V4_EXPECTED_SPOTS <- c(
  Cap.area1 = 3384L, Cap.area2 = 3318L, Cap.area3 = 1197L,
  Cap.area4 = 1723L, Cap.area5 = 1827L, Cap.area6 = 4992L,
  Cap.area7 = 4992L, Cap.area8 = 4674L, Cap.area9 = 4923L
)
WP3V4_PC_AREAS <- c("Cap.area1", "Cap.area2", "Cap.area3")
WP3V4_FORBIDDEN_FIELDS <- c(
  "patient", "specimen", "region", "treatment", "relapse", "OS",
  "purity", "composition", "malignant", "cluster", "class", "taxonomy"
)

wp3v4_exact <- function(x, name) x[[name, exact = TRUE]]

wp3v4_program_ids <- function() {
  contract <- read.csv(WP3V4_PROGRAM_CONTRACT, check.names = FALSE)
  contract <- contract[order(contract$program_order), , drop = FALSE]
  if (nrow(contract) != 6L || anyDuplicated(contract$program_id)) {
    stop("Canonical WP1 program contract is invalid")
  }
  contract$program_id
}

wp3v4_read_csv <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 0) stop("Missing or empty file: ", path)
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

wp3v4_single_file <- function(directory, pattern) {
  paths <- list.files(directory, pattern = pattern, full.names = TRUE)
  if (length(paths) != 1L) stop("Expected exactly one file for pattern: ", pattern)
  paths[[1L]]
}

wp3v4_score_schema_ok <- function(x, expected_rows, programs = wp3v4_program_ids()) {
  is.data.frame(x) && nrow(x) == expected_rows &&
    identical(names(x), c("barcode", programs)) &&
    !anyNA(x$barcode) && !anyDuplicated(x$barcode) &&
    all(is.finite(as.matrix(x[, programs, drop = FALSE])))
}

wp3v4_pc_schema_ok <- function(x, expected_rows) {
  required <- c("barcode", "PC1", "PC2", "pc_role")
  is.data.frame(x) && nrow(x) == expected_rows && all(required %in% names(x)) &&
    !anyNA(x$barcode) && !anyDuplicated(x$barcode) &&
    all(is.finite(as.matrix(x[, c("PC1", "PC2"), drop = FALSE])))
}

wp3v4_load_v2_area1 <- function() {
  read_part <- function(dir, pattern) wp3v4_read_csv(
    wp3v4_single_file(file.path(WP3V4_V2, dir), pattern))
  input_qc <- read_part("01_input_qc", "Cap[.]area1.*input_qc[.]csv$")
  list(
    capture_area_id = input_qc$capture_area_id[[1L]],
    area_role = "CONTEXT_ONLY", source = "V2_REUSED_COMPLETE",
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
    pc_artifact = NULL, matrix_class = "dgCMatrix",
    sparse_coercion_warning_count = 0L, other_warning_count = 0L,
    stack_imbalance = FALSE
  )
}

wp3v4_load_v3_area <- function(short) {
  id <- unname(WP3V4_AREA_IDS[[short]])
  directory <- file.path(WP3V4_V3, "areas", id)
  marker <- wp3v4_read_csv(file.path(directory, "AREA_COMPLETE_V3.csv"))
  if (nrow(marker) != 1L || marker$status[[1L]] != "COMPLETE") {
    stop("Area completion marker failed: ", short)
  }
  read_part <- function(name) wp3v4_read_csv(file.path(directory, name))
  input_qc <- read_part("01_input_qc.csv")
  raw_pc_path <- file.path(directory, "11_exploratory_pc.csv")
  raw_pc <- if (file.exists(raw_pc_path)) read_part("11_exploratory_pc.csv") else NULL
  pc_is_authoritative <- short %in% WP3V4_PC_AREAS &&
    wp3v4_pc_schema_ok(raw_pc, WP3V4_EXPECTED_SPOTS[[short]])
  artifact <- NULL
  if (!is.null(raw_pc) && !pc_is_authoritative) {
    status_value <- if (ncol(raw_pc) == 1L && nrow(raw_pc) == 1L) as.character(raw_pc[[1L]][[1L]]) else ""
    artifact <- data.frame(
      capture_area = short, file = raw_pc_path, rows = nrow(raw_pc),
      columns = paste(names(raw_pc), collapse = ";"), pc_status = status_value,
      artifact_type = "TECHNICAL_STATUS_SENTINEL_MISNAMED_AS_PC",
      authoritative_pc_result = FALSE, excluded_from_final_registry = TRUE,
      stringsAsFactors = FALSE
    )
  }
  eligibility <- read_part("02_eligibility.csv")
  list(
    capture_area_id = input_qc$capture_area_id[[1L]],
    area_role = eligibility$area_role[[1L]], source = "V3_RECOMPUTED_COMPLETE",
    input_qc = input_qc, eligibility = eligibility,
    primary_scores = read_part("03_primary_ucell.csv"),
    sct_scores = read_part("04_sct_ucell.csv"),
    ssgsea_scores = read_part("05_ssgsea.csv"),
    score_qc = read_part("06_score_qc.csv"),
    concordance = read_part("07_method_concordance.csv"),
    edges = read_part("08_adjacency.csv"),
    autocorrelation = read_part("09_moran_geary.csv"),
    bivariate = read_part("10_bivariate_moran.csv"),
    pc = if (pc_is_authoritative) raw_pc else NULL,
    pc_artifact = artifact, matrix_class = "dgCMatrix",
    sparse_coercion_warning_count = 0L, other_warning_count = 0L,
    stack_imbalance = FALSE
  )
}

wp3v4_load_completed_results <- function() {
  results <- list(Cap.area1 = wp3v4_load_v2_area1())
  for (short in names(WP3V4_AREA_IDS)[-1L]) results[[short]] <- wp3v4_load_v3_area(short)
  artifacts <- do.call(rbind, Filter(Negate(is.null), lapply(results, wp3v4_exact, name = "pc_artifact")))
  if (is.null(artifacts)) artifacts <- data.frame()
  source_registry <- data.frame(
    capture_area = names(results),
    capture_area_id = vapply(results, wp3v4_exact, character(1), name = "capture_area_id"),
    source = vapply(results, wp3v4_exact, character(1), name = "source"),
    stringsAsFactors = FALSE
  )
  list(results = results, pc_artifacts = artifacts, source_registry = source_registry,
    raw_expression_read = FALSE, area_node_calls = 0L, score_computations = 0L,
    spatial_statistic_computations = 0L)
}

wp3v4_verify_scene_manifest <- function() {
  if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required for SHA-256 verification")
  manifest <- wp3v4_read_csv(WP3V4_SCENE_MANIFEST)
  if (nrow(manifest) != 138L || any(manifest$overwrite_allowed)) {
    stop("Completed prevalidator manifest contract failed")
  }
  observed <- list.files(WP3V4_SCENE, recursive = TRUE, full.names = TRUE,
    include.dirs = FALSE, no.. = TRUE)
  observed <- observed[!grepl("(^|/)finalization_v2(/|$)", gsub("\\\\", "/", observed))]
  relative <- substring(normalizePath(observed, winslash = "/", mustWork = TRUE),
    nchar(normalizePath(WP3V4_SCENE, winslash = "/", mustWork = TRUE)) + 2L)
  if (!setequal(relative, manifest$relative_path)) stop("Frozen 138-file set changed")
  for (i in seq_len(nrow(manifest))) {
    path <- file.path(WP3V4_SCENE, manifest$relative_path[[i]])
    info <- file.info(path)
    if (!identical(unname(info$size), as.numeric(manifest$size[[i]]))) stop("Size changed: ", path)
    if (digest::digest(path, file = TRUE, algo = "sha256") != manifest$sha256[[i]]) stop("SHA changed: ", path)
    observed_mtime <- format(info$mtime, "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")
    if (observed_mtime != manifest$mtime[[i]]) stop("Mtime changed: ", path)
  }
  TRUE
}
