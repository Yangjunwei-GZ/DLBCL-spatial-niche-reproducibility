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
source(file.path(DLBCL_SUPPLEMENTARY_CODE_ROOT, "spatial", "37v_wp3_completed_result_validator_v4.R"))

token_name <- "DLBCL_REVISION_ALLOW_WP3_SPATIAL"
authorized_hash <- "a22182880693db23c40c173c8d4e6f1ae04edacd6a7829450045f853d64db471"
if (!requireNamespace("digest", quietly = TRUE)) stop("digest package is required")
entry_token <- Sys.getenv(token_name, unset = "")
if (!nzchar(entry_token) || digest::digest(entry_token, algo = "sha256", serialize = FALSE) != authorized_hash) {
  stop("Valid WP3 validation-only token is required")
}
on.exit(Sys.unsetenv(token_name), add = TRUE)

stage_root <- paste0(WP3V4_FINAL, "_staging")
if (dir.exists(WP3V4_FINAL) || file.exists(WP3V4_FINAL) ||
    dir.exists(stage_root) || file.exists(stage_root)) {
  stop("Validation-only finalization is create-once")
}

bundle <- wp3v4_load_completed_results()
validator <- wp3v4_validate_completed_results(bundle, output_root = WP3V4_FINAL,
  stop_on_failure = TRUE, verify_scene = TRUE)
if (any(!validator$status)) stop("V4 validator did not pass")

results <- wp3v4_exact(bundle, "results")
area_registry <- wp3v4_exact(bundle, "source_registry")
area_registry$status <- "COMPLETE"
area_registry$spot_count <- unname(WP3V4_EXPECTED_SPOTS[area_registry$capture_area])
area_registry$raw_expression_reread <- FALSE
area_registry$area_rerun <- FALSE

score_audit <- wp3v4_read_csv(file.path(WP3V4_DIR,
  "WP3_SCORE_OUTPUT_SCHEMA_AUTHORITY_AUDIT.csv"))
pc_registry <- do.call(rbind, lapply(WP3V4_PC_AREAS, function(short) {
  pc <- wp3v4_exact(results[[short]], "pc")
  data.frame(
    capture_area = short, capture_area_id = wp3v4_exact(results[[short]], "capture_area_id"),
    source = wp3v4_exact(results[[short]], "source"), rows = nrow(pc),
    columns = paste(names(pc), collapse = ";"), authoritative_pc_result = TRUE,
    stringsAsFactors = FALSE
  )
}))
excluded <- wp3v4_exact(bundle, "pc_artifacts")
autocorrelation <- do.call(rbind, lapply(results, wp3v4_exact, name = "autocorrelation"))
bivariate <- do.call(rbind, lapply(results, wp3v4_exact, name = "bivariate"))
concordance <- do.call(rbind, lapply(results, wp3v4_exact, name = "concordance"))
scope <- data.frame(
  scope = c("patient", "specimen", "region", "treatment_relapse", "OS", "purity",
    "composition_adjustment", "malignant_assignment", "cluster", "class", "taxonomy"),
  status = c(rep("NOT_RUN", 9L), "NOT_CREATED", "NOT_ASSIGNED"),
  stringsAsFactors = FALSE
)
protocol <- data.frame(
  protocol = "WP3_VALIDATION_ONLY_FINALIZATION_V4",
  governance_baseline = "ca9c654531a839ac81dc988e7217824b830f6b05",
  seed = WP3V4_SEED, permutations = 9999L,
  area_node_calls = 0L, raw_expression_read = FALSE,
  score_computations = 0L, spatial_statistic_computations = 0L,
  final_k = "NOT_SELECTED", taxonomy = "NOT_ASSIGNED",
  stringsAsFactors = FALSE
)
startup <- data.frame(
  R_version = paste(R.version$major, R.version$minor, sep = "."),
  token_status = "VALID_NOT_RECORDED", locale = Sys.getlocale(),
  started_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  raw_expression_read = FALSE, area_node_calls = 0L,
  stringsAsFactors = FALSE
)
completion <- data.frame(
  status = "COMPLETED", authoritative_areas = 9L, score_files = 27L,
  moran_geary_rows = nrow(autocorrelation), bivariate_rows = nrow(bivariate),
  authoritative_pc_areas = 3L, excluded_pc_sentinels = nrow(excluded),
  validator_checks = nrow(validator), validator_pass = all(validator$status),
  warning_count = 0L, stack_imbalance = FALSE,
  raw_expression_read = FALSE, area_node_calls = 0L, score_computations = 0L,
  spatial_statistic_computations = 0L, final_k = "NOT_SELECTED",
  taxonomy = "NOT_ASSIGNED", stringsAsFactors = FALSE
)

dir.create(stage_root, recursive = TRUE)
dir.create(file.path(stage_root, "logs"))
dir.create(file.path(stage_root, "protocol"))
write_once <- function(x, name) {
  path <- file.path(stage_root, name)
  if (file.exists(path)) stop("Create-once refusal: ", path)
  write.csv(x, path, row.names = FALSE, na = "")
}
write_once(area_registry, "WP3_FINAL_AREA_AUTHORITY_REGISTRY.csv")
write_once(score_audit, "WP3_FINAL_SCORE_SCHEMA_AUDIT.csv")
write_once(pc_registry, "WP3_FINAL_PC_AUTHORITY_REGISTRY.csv")
write_once(excluded, "WP3_EXCLUDED_TECHNICAL_PC_SENTINELS.csv")
write_once(autocorrelation, "WP3_FINAL_MORAN_GEARY_AUTHORITY.csv")
write_once(bivariate, "WP3_FINAL_BIVARIATE_AUTHORITY.csv")
write_once(concordance, "WP3_FINAL_METHOD_CONCORDANCE_AUTHORITY.csv")
write_once(scope, "WP3_FINAL_SCOPE_ISOLATION.csv")
write_once(validator, "WP3_FINAL_VALIDATOR_V4.csv")
write_once(completion, "WP3_FINALIZATION_COMPLETION_REPORT.csv")
write_once(startup, "logs/WP3_VALIDATION_FINALIZATION_STARTUP_V4.csv")
write_once(protocol, "protocol/WP3_FINALIZATION_PROTOCOL_V4.csv")
write_once(wp3v4_read_csv(WP3V4_SCENE_MANIFEST),
  "protocol/WP3_COMPLETED_SCENE_REFERENCE_V4.csv")

spot_text <- paste(paste0(names(WP3V4_EXPECTED_SPOTS), "=", WP3V4_EXPECTED_SPOTS), collapse = ", ")
report <- c(
  "# WP3 Spatial Continuous Validation Report",
  "",
  "**VALIDATION-ONLY FINALIZATION PASS**",
  "",
  "1. All nine capture areas are complete. Cap.area1 is reused read-only from V2; Cap.area2-Cap.area9 are complete V3 results.",
  paste0("2. Frozen spot counts: ", spot_text, "."),
  "3. All 27 existing score files use the authoritative schema `barcode` plus the six raw canonical program IDs; the retired `_UCell` suffix expectation is not used.",
  "4. LogNormalize+UCell, SCTransform v2+UCell, and LogNormalize+GSVA ssGSEA outputs have exact spot counts, unique barcodes, finite scores, and canonical program order.",
  "5. Method concordance and adjacency are complete for all nine areas.",
  paste0("6. Spatial authority contains ", nrow(autocorrelation), " Moran/Geary rows and ", nrow(bivariate), " bivariate Moran rows, all retaining 9,999 permutations."),
  "7. PC authority is restricted to Cap.area1-Cap.area3. The six Cap.area4-Cap.area9 one-row status sentinels are preserved technical artifacts and excluded from PC inference.",
  "8. Antigen-presentation coverage is 13/22 and exploratory-only in Cap.area4-Cap.area9; no missing gene was substituted or imputed.",
  "9. The frozen 27/27 dense-edge equivalence authority remains PASS; no spatial statistic was recomputed during finalization.",
  "10. Warning count is zero and stack imbalance is FALSE.",
  "11. Patient, specimen, region, treatment/relapse, OS, purity, composition adjustment, malignant assignment, clustering, and class/taxonomy inference were not run.",
  "12. Final k remains NOT_SELECTED and taxonomy remains NOT_ASSIGNED.",
  paste0("13. V4 validator: ", sum(validator$status), "/", nrow(validator), " PASS."),
  "",
  "This is a continuous spot-level spatial robustness analysis and does not establish discrete spatial classes."
)
report_path <- file.path(stage_root, "WP3_SPATIAL_CONTINUOUS_VALIDATION_REPORT.md")
writeLines(report, report_path, useBytes = TRUE)

if (!file.rename(stage_root, WP3V4_FINAL)) stop("Atomic finalization rename failed")
cat("WP3 VALIDATION-ONLY FINALIZATION V4 PASS\n")
