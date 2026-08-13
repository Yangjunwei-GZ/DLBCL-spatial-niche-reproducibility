DLBCL_PROJECT_ROOT <- normalizePath(
  Sys.getenv("DLBCL_PROJECT_ROOT", unset = getwd()),
  winslash = "/", mustWork = FALSE
)
DLBCL_SUPPLEMENTARY_CODE_ROOT <- normalizePath(
  Sys.getenv("DLBCL_SUPPLEMENTARY_CODE_ROOT",
    unset = file.path(DLBCL_PROJECT_ROOT, "Supplementary_Code")),
  winslash = "/", mustWork = FALSE
)

options(stringsAsFactors = FALSE, warn = 1)

ROOT <- DLBCL_PROJECT_ROOT
REV <- file.path(ROOT, "revision_2026_reviewer_response")
P06 <- file.path(REV, "06p_wp3_matched_null_and_depth_sensitivity")
PROTO <- file.path(P06, "00_protocol_freeze")
AMEND <- file.path(PROTO, "amendments")
OUT <- file.path(P06, "01_execution_outputs")
RUN <- file.path(OUT, "run_control")
ATT4 <- file.path(RUN, "attempt_004")
VAL <- file.path(OUT, "validation")
LOG <- file.path(P06, "02_execution_logs", "05_06p_partitioned_membership_resume.log")

TOKEN <- "AUTHORIZE_06P_MATCHED_NULL_AND_DEPTH_EXECUTION_SEED_20260804"
O001_HASH <- "f3f557ceb7e6baf3c0f3f2e076e22523d338e72d2522e08cab6a38266ea32147"

.libPaths(c(
  file.path(REV, "06b_wp3b_spatial_scope_method_resolution/.wp3_r_library"),
  file.path(REV, "04b_stage4_environment_freeze/stage4_renv_project/renv/library/windows/R-4.5/x86_64-w64-mingw32"),
  .libPaths()
))

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Missing required package: digest")
}

dir.create(VAL, recursive = TRUE, showWarnings = FALSE)

final_paths <- c(
  file.path(VAL, "06P_EXECUTION_VALIDATION.csv"),
  file.path(VAL, "06P_OUTPUT_MANIFEST.csv"),
  file.path(OUT, "06P_FINAL_REPORT.md")
)
if (any(file.exists(final_paths))) {
  stop("create-once refusal for final validation artifact(s): ",
       paste(final_paths[file.exists(final_paths)], collapse = "; "))
}

sha_file <- function(path) digest::digest(path, file = TRUE, algo = "sha256")
read_csv <- function(path) read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
read_text <- function(path) readLines(path, warn = FALSE, encoding = "UTF-8")
read_utf16_text <- function(path) readLines(file(path, encoding = "UTF-16LE"), warn = FALSE)
csv_ok <- function(path) {
  tryCatch({
    read.csv(path, nrows = 5L, check.names = FALSE, stringsAsFactors = FALSE)
    TRUE
  }, error = function(e) FALSE)
}
truthy <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x %in% c("TRUE", "T", "1", "PASS", "VALID", "COMPLETE")
}
add_check <- function(rows, check_id, check, ok, detail) {
  rows[[length(rows) + 1L]] <- data.frame(
    check_id = check_id,
    check = check,
    status = if (isTRUE(ok)) "PASS" else "FAIL",
    detail = paste(detail, collapse = ""),
    stringsAsFactors = FALSE
  )
  rows
}

protocol_validation <- read_csv(file.path(PROTO, "06P_PROTOCOL_VALIDATION.csv"))
planned_outputs <- read_csv(file.path(PROTO, "06P_PLANNED_OUTPUT_REGISTRY.csv"))
seed_registry <- read_csv(file.path(PROTO, "06P_SEED_REGISTRY.csv"))
amd1 <- read_csv(file.path(AMEND, "AMENDMENT_001_06P_PREFLIGHT_IMPLEMENTATION_CLARIFICATION", "06P_AMENDMENT_VALIDATION.csv"))
amd2 <- read_csv(file.path(AMEND, "AMENDMENT_002_06P_ROOT_OUTPUT_TARGET_CLARIFICATION", "06P_AMENDMENT_002_VALIDATION.csv"))
amd3 <- read_csv(file.path(AMEND, "AMENDMENT_003_06P_PREEXISTING_CONTROL_OUTPUT_STATE_CLARIFICATION", "06P_AMENDMENT_003_VALIDATION.csv"))
amd4 <- read_csv(file.path(AMEND, "AMENDMENT_004_06P_PARTITIONED_MEMBERSHIP_SALVAGE_AND_RESUME", "06P_AMENDMENT_004_VALIDATION.csv"))
preflight4 <- read_csv(file.path(ATT4, "06P_PREFLIGHT_VALIDATION_ATTEMPT004.csv"))
planned_state4 <- read_csv(file.path(ATT4, "06P_PLANNED_OUTPUT_STATE_CHECKS_ATTEMPT004.csv"))
input_crosswalk4 <- read_csv(file.path(ATT4, "06P_INPUT_AUTHORITY_CROSSWALK_ATTEMPT004.csv"))

membership_validation <- read_csv(file.path(OUT, "matched_null", "MATCHED_RANDOM_SET_MEMBERSHIP_VALIDATION.csv"))
partition_validation <- read_csv(file.path(RUN, "membership_salvage_001", "06P_PARTITION_FILE_VALIDATION.csv"))
reused_checks <- read_csv(file.path(RUN, "membership_salvage_001", "06P_REUSED_SCIENTIFIC_OUTPUT_CHECKS.csv"))
matching_qc <- read_csv(file.path(OUT, "matched_null", "MATCHING_QC.csv"))
matched_null <- read_csv(file.path(OUT, "matched_null", "MATCHED_NULL_MORAN_GEARY.csv"))
empirical <- read_csv(file.path(OUT, "matched_null", "MATCHED_NULL_EMPIRICAL_TESTS.csv"))
spot_depth <- read_csv(file.path(OUT, "gene_summaries", "SPOT_DEPTH_COVARIATES.csv"))
gene_summary <- read_csv(file.path(OUT, "gene_summaries", "AREA_GENE_MATCHING_SUMMARIES.csv"))
depth_models <- read_csv(file.path(OUT, "depth", "DEPTH_MODEL_DIAGNOSTICS.csv"))
depth_spatial <- read_csv(file.path(OUT, "depth", "DEPTH_RESIDUAL_SPATIAL_STATISTICS.csv"))
depth_perm <- read_csv(file.path(OUT, "depth", "DEPTH_RESIDUAL_PERMUTATION_SUMMARY.csv"))
fdr <- read_csv(file.path(OUT, "fdr", "06P_SEPARATE_FDR_RESULTS.csv"))

token_value <- trimws(paste(read_text(file.path(P06, "00_execution_authorization", "06P_EXECUTION_TOKEN.txt")), collapse = ""))
continuation_auth <- paste(read_text(file.path(ATT4, "06P_CONTINUATION_AUTHORIZATION.md")), collapse = "\n")
execution_auth <- paste(read_text(file.path(P06, "00_execution_authorization", "06P_EXECUTION_AUTHORIZATION.md")), collapse = "\n")
log_text <- if (file.exists(LOG)) paste(read_utf16_text(LOG), collapse = "\n") else ""
git_tracked <- system2("git", c("-C", shQuote(REV), "diff", "--name-only"), stdout = TRUE)
git_staged <- system2("git", c("-C", shQuote(REV), "diff", "--cached", "--name-only"), stdout = TRUE)

scientific_csv_paths <- file.path(OUT, c(
  "gene_summaries/AREA_GENE_MATCHING_SUMMARIES.csv",
  "gene_summaries/SPOT_DEPTH_COVARIATES.csv",
  "matched_null/MATCHED_RANDOM_SET_MEMBERSHIP.csv",
  "matched_null/MATCHING_QC.csv",
  "matched_null/MATCHED_NULL_MORAN_GEARY.csv",
  "matched_null/MATCHED_NULL_EMPIRICAL_TESTS.csv",
  "depth/DEPTH_MODEL_DIAGNOSTICS.csv",
  "depth/DEPTH_RESIDUAL_SPATIAL_STATISTICS.csv",
  "depth/DEPTH_RESIDUAL_PERMUTATION_SUMMARY.csv",
  "fdr/06P_SEPARATE_FDR_RESULTS.csv",
  "failures/06P_FAILURE_LOG.csv"
))

seed_unique <- length(unique(seed_registry$derived_seed))
matched_seed_count <- sum(seed_registry$analysis_type == "MATCHED_GENE_SET_REPLICATE")
moran_seed_count <- sum(seed_registry$analysis_type == "DEPTH_RESIDUAL_MORAN_PERMUTATION_ENDPOINT")
geary_seed_count <- sum(seed_registry$analysis_type == "DEPTH_RESIDUAL_GEARY_PERMUTATION_ENDPOINT")
membership_metric <- function(id) {
  z <- membership_validation[membership_validation$check_id == id, , drop = FALSE]
  if (nrow(z) != 1L) return(NA_integer_)
  as.integer(z$observed)
}
membership_pass <- function(id) {
  z <- membership_validation[membership_validation$check_id == id, , drop = FALSE]
  nrow(z) == 1L && identical(z$status, "PASS")
}
fdr_family_table <- aggregate(raw_p ~ endpoint_type + role_family, fdr, length)
names(fdr_family_table)[3] <- "observed_count"
expected_fdr_strings <- as.vector(outer(
  c("MATCHED_NULL_MORAN", "MATCHED_NULL_GEARY", "DEPTH_ADJUSTED_MORAN", "DEPTH_ADJUSTED_GEARY"),
  c("PRIMARY_DLBCL 26", "EXPLORATORY_ANTIGEN 6", "CONTEXT_ONLY 22"),
  paste
))
observed_fdr_strings <- paste(fdr_family_table$endpoint_type, fdr_family_table$role_family, fdr_family_table$observed_count)

rows <- list()
rows <- add_check(rows, "06P-FV001", "01_execution_outputs exists", dir.exists(OUT), OUT)
rows <- add_check(rows, "06P-FV002", "Exact execution token", identical(token_value, TOKEN), token_value)
rows <- add_check(rows, "06P-FV003", "Human authorization recorded", grepl("PROVIDED", continuation_auth, fixed = TRUE) && grepl("PROVIDED", execution_auth, fixed = TRUE), "authorization text contains PROVIDED")
rows <- add_check(rows, "06P-FV004", "Original protocol validation PASS", all(protocol_validation$status == "PASS"), paste(sum(protocol_validation$status == "PASS"), "of", nrow(protocol_validation)))
rows <- add_check(rows, "06P-FV005", "Amendments 001-004 validation PASS", all(amd1$status == "PASS") && all(amd2$status == "PASS") && all(amd3$status == "PASS") && all(amd4$status == "PASS"), paste(c(nrow(amd1), nrow(amd2), nrow(amd3), nrow(amd4)), collapse = "/"))
rows <- add_check(rows, "06P-FV006", "Attempt 004 preflight structured validation PASS", all(preflight4$status == "PASS"), paste(sum(preflight4$status == "PASS"), "of", nrow(preflight4)))
rows <- add_check(rows, "06P-FV007", "Raw component authority 81/81 PASS", all(input_crosswalk4$attempt004_component_status == "81_OF_81_PASS") && all(input_crosswalk4$effective_raw_input_status == "PASS"), paste(input_crosswalk4$attempt004_component_status, collapse = ";"))
rows <- add_check(rows, "06P-FV008", "Seed registry 54108 unique split", nrow(seed_registry) == 54108L && seed_unique == 54108L && matched_seed_count == 54000L && moran_seed_count == 54L && geary_seed_count == 54L, paste("rows", nrow(seed_registry), "unique", seed_unique, "matched", matched_seed_count, "moran", moran_seed_count, "geary", geary_seed_count))
rows <- add_check(rows, "06P-FV009", "Planned output preflight states nonblocking", nrow(planned_state4) == 16L && all(planned_state4$effective_status == "PASS") && !any(truthy(planned_state4$scientific_execution_blocking)), paste(nrow(planned_state4), "planned outputs"))
rows <- add_check(rows, "06P-FV010", "O001 preserved control output unchanged", sha_file(file.path(RUN, "06P_INPUT_SHA256_REGISTRY.csv")) == O001_HASH, O001_HASH)
rows <- add_check(rows, "06P-FV011", "Membership authority intact", all(membership_validation$status == "PASS") && all(partition_validation$scientific_validation == "PASS") && all(truthy(partition_validation$sha256_matches)), paste(nrow(partition_validation), "partitions;", sum(grepl("^SALVAGED__", partition_validation$batch_id)), "salvaged;", sum(grepl("^GENERATED__", partition_validation$batch_id)), "generated"))
rows <- add_check(rows, "06P-FV012", "Validated membership key count", membership_pass("COMPOSITE_KEY_COUNT") && membership_metric("COMPOSITE_KEY_COUNT") == 54000L && membership_pass("COMPOSITE_KEY_UNIQUENESS") && membership_pass("AREA_PROGRAM_COMBINATIONS") && membership_metric("AREA_PROGRAM_COMBINATIONS") == 54L && membership_pass("REPLICATES_PER_COMBINATION"), paste("keys", membership_metric("COMPOSITE_KEY_COUNT"), "combinations", membership_metric("AREA_PROGRAM_COMBINATIONS")))
rows <- add_check(rows, "06P-FV013", "HLA-DRB1 exclusion PASS", any(reused_checks$check_id == "AREA1_HLA_DRB1_EXCLUDED" & reused_checks$status == "PASS"), "AREA1_HLA_DRB1_EXCLUDED")
rows <- add_check(rows, "06P-FV014", "Spot depth covariates generated", nrow(spot_depth) == 31030L && all(c("nCount_Spatial", "nFeature_Spatial") %in% names(spot_depth)) && all(is.finite(spot_depth$nCount_Spatial)) && all(spot_depth$nCount_Spatial >= 0) && all(is.finite(spot_depth$nFeature_Spatial)) && all(spot_depth$nFeature_Spatial > 0), paste("rows", nrow(spot_depth)))
rows <- add_check(rows, "06P-FV015", "Gene matching summaries generated", nrow(gene_summary) == 162315L && length(unique(gene_summary$capture_area_id)) == 9L, paste("rows", nrow(gene_summary)))
rows <- add_check(rows, "06P-FV016", "Matched-null QC 54/54 valid 1000", nrow(matching_qc) == 54L && all(matching_qc$attempted_sets == 1000L) && all(matching_qc$valid_sets == 1000L) && all(matching_qc$failed_sets == 0L) && all(matching_qc$null_validity_decision == "VALID"), paste("valid sets", sum(matching_qc$valid_sets)))
rows <- add_check(rows, "06P-FV017", "Matched-null result CSVs complete", nrow(matched_null) == 54000L && all(matched_null$validity_status == "VALID") && nrow(empirical) == 54L && all(empirical$valid_null_count == 1000L) && all(empirical$null_validity_status == "VALID"), paste(nrow(matched_null), "null rows;", nrow(empirical), "empirical rows"))
rows <- add_check(rows, "06P-FV018", "Depth models 54/54 valid", nrow(depth_models) == 54L && all(depth_models$model_status == "VALID") && all(depth_models$design_matrix_rank == 3L) && all(truthy(depth_models$residual_complete)), paste("valid", sum(depth_models$model_status == "VALID"), "rank deficient", sum(depth_models$model_status == "RANK_DEFICIENT")))
rows <- add_check(rows, "06P-FV019", "Depth residual Moran/Geary endpoints complete", nrow(depth_spatial) == 54L && all(depth_spatial$evaluation_status == "VALID") && nrow(depth_perm) == 108L && sum(depth_perm$endpoint == "DEPTH_ADJUSTED_MORAN") == 54L && sum(depth_perm$endpoint == "DEPTH_ADJUSTED_GEARY") == 54L && all(depth_perm$status == "COMPLETE") && all(depth_perm$permutations == 9999L), paste("spatial", nrow(depth_spatial), "permutation", nrow(depth_perm), "permutations", paste(unique(depth_perm$permutations), collapse = ";")))
rows <- add_check(rows, "06P-FV020", "FDR results complete", nrow(fdr) == 216L && length(unique(paste(fdr$endpoint_type, fdr$role_family))) == 12L && setequal(observed_fdr_strings, expected_fdr_strings), paste("rows", nrow(fdr), "families", length(unique(paste(fdr$endpoint_type, fdr$role_family)))))
rows <- add_check(rows, "06P-FV021", "Scientific output CSVs parse", all(file.exists(scientific_csv_paths)) && all(vapply(scientific_csv_paths, csv_ok, logical(1))), paste(sum(vapply(scientific_csv_paths, function(p) file.exists(p) && csv_ok(p), logical(1))), "of", length(scientific_csv_paths)))
rows <- add_check(rows, "06P-FV022", "Execution log reached FDR then halted in validator construction", file.exists(LOG) && grepl("FDR_START", log_text, fixed = TRUE) && grepl("arguments imply differing number of rows: 25, 22", log_text, fixed = TRUE), "scientific phases completed; prior halt was validation row construction")
rows <- add_check(rows, "06P-FV023", "Raw matrices not read by this final validator", TRUE, "validation-only script reads existing outputs and logs only")
rows <- add_check(rows, "06P-FV024", "Git tracked and staged changes zero before final artifact creation", length(git_tracked) == 0L && length(git_staged) == 0L, paste("tracked", length(git_tracked), "staged", length(git_staged)))

validation <- do.call(rbind, rows)
if (anyDuplicated(validation$check_id) || any(!nzchar(validation$check_id))) {
  stop("Validation check identifiers are not unique and nonempty")
}
if (!all(validation$status %in% c("PASS", "FAIL", "NOT_APPLICABLE"))) {
  stop("Unexpected validation status")
}
if (any(validation$status != "PASS")) {
  stop("Final validation failed before create-once output: ",
       paste(validation$check_id[validation$status != "PASS"], collapse = "; "))
}

authority_decision <- "FINAL_06P_AUTHORITY"
controlling_status <- "FINAL_06P_AUTHORITY_ASSIGNED"

write.csv(validation, final_paths[1], row.names = FALSE, na = "", fileEncoding = "UTF-8")

report <- c(
  "# 06p Final Report",
  "",
  paste("Controlling status:", controlling_status),
  paste("Authority decision:", authority_decision),
  "",
  "## Validation-Only Continuation",
  "The final validator was run after the authorized Rscript completed all scientific output phases and halted while constructing the original final validation data frame.",
  "This continuation reads existing protocol, run-control, log, and scientific-output records only. It does not read raw matrix values and does not rerun membership generation, matched-null scoring, depth models, permutations, FDR, or any scientific analysis.",
  "",
  "## Execution Summary",
  paste("Execution token:", TOKEN),
  "Human authorization: PROVIDED",
  paste("Scientific-result CSVs present:", sum(file.exists(scientific_csv_paths)), "of", length(scientific_csv_paths)),
  paste("Matched-null combinations:", nrow(matching_qc), "of 54"),
  paste("Matched-null valid sets:", sum(matching_qc$valid_sets)),
  paste("Matched-null result rows:", nrow(matched_null)),
  paste("Depth models:", sum(depth_models$model_status == "VALID"), "valid,", sum(depth_models$model_status == "RANK_DEFICIENT"), "rank-deficient"),
  paste("Depth residual endpoints:", nrow(depth_perm), "with permutations", paste(unique(depth_perm$permutations), collapse = ";")),
  paste("FDR result rows:", nrow(fdr)),
  "",
  "## Preflight Reconciliation Note",
  "The Attempt 004 markdown report line `Amendment 003 validation: FAIL` is reporting text only. The structured Amendment 003 validation CSV and Attempt 004 preflight validation CSV are PASS, and no scientific parameter requires amendment.",
  "",
  "## Final Validation",
  paste("Validation checks:", sum(validation$status == "PASS"), "of", nrow(validation), "PASS"),
  paste("Final authority status:", authority_decision),
  "",
  "Negative, failed, and heterogeneous findings are retained in the scientific CSV outputs.",
  "No causal spatial organization, cell-cell communication, patient-level replication, independent normal-control validation, spatial cell-state reconstruction, DLBCL specificity of context-area results, or complete technical-confounding-removal claim is made."
)
writeLines(report, final_paths[3], useBytes = TRUE)

manifest_paths <- c(
  file.path(P06, "00_execution_authorization", "06P_EXECUTION_TOKEN.txt"),
  file.path(P06, "00_execution_authorization", "06P_EXECUTION_AUTHORIZATION.md"),
  file.path(PROTO, "06P_PROTOCOL_VALIDATION.csv"),
  file.path(PROTO, "06P_PLANNED_OUTPUT_REGISTRY.csv"),
  file.path(PROTO, "06P_SEED_REGISTRY.csv"),
  file.path(ATT4, "06P_PREFLIGHT_VALIDATION_ATTEMPT004.csv"),
  file.path(ATT4, "06P_CONTINUATION_AUTHORIZATION.md"),
  file.path(RUN, "06P_INPUT_SHA256_REGISTRY.csv"),
  file.path(RUN, "membership_salvage_001", "06P_PARTITION_FILE_VALIDATION.csv"),
  file.path(RUN, "membership_salvage_001", "06P_REUSED_SCIENTIFIC_OUTPUT_CHECKS.csv"),
  scientific_csv_paths,
  final_paths[1],
  final_paths[2],
  final_paths[3]
)
manifest_exists <- file.exists(manifest_paths)
is_self_manifest <- normalizePath(manifest_paths, winslash = "/", mustWork = FALSE) ==
  normalizePath(final_paths[2], winslash = "/", mustWork = FALSE)

manifest <- data.frame(
  path = manifest_paths,
  exists = manifest_exists | is_self_manifest,
  size_bytes = ifelse(manifest_exists & !is_self_manifest, as.numeric(file.info(manifest_paths)$size), NA_real_),
  last_write_time = ifelse(manifest_exists & !is_self_manifest, as.character(file.info(manifest_paths)$mtime), NA_character_),
  sha256 = ifelse(manifest_exists & !is_self_manifest, vapply(manifest_paths, function(p) if (file.exists(p)) sha_file(p) else NA_character_, character(1)), NA_character_),
  parse_status = ifelse(grepl("\\.csv$", manifest_paths, ignore.case = TRUE),
                        ifelse(manifest_exists & vapply(manifest_paths, function(p) if (file.exists(p)) csv_ok(p) else FALSE, logical(1)), "PASS", ifelse(is_self_manifest, "SELF_REFERENCE", "MISSING")),
                        "NOT_CSV"),
  creation_origin = ifelse(manifest_paths %in% c(final_paths[1], final_paths[2], final_paths[3]), "VALIDATION_ONLY_CONTINUATION", "PREEXISTING_AUTHORIZED_OUTPUT"),
  overwrite_performed = "FALSE",
  stringsAsFactors = FALSE
)
write.csv(manifest, final_paths[2], row.names = FALSE, na = "", fileEncoding = "UTF-8")

message("06P_FINAL_VALIDATION_ONLY_STATUS=", controlling_status)
