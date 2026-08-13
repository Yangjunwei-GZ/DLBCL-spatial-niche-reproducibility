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
OUT <- file.path(P06, "01_execution_outputs")
RUN <- file.path(OUT, "run_control")
ATT4 <- file.path(RUN, "attempt_004")
TOKEN <- "AUTHORIZE_06P_MATCHED_NULL_AND_DEPTH_EXECUTION_SEED_20260804"
PARTIAL_HASH <- "335f2f67a5ae6432df456faaded956168426ff470426e2f26a25bf577fe94e14"
ORIGINAL_SCRIPT_HASH <- "d15d994d7aa66900a53f82b2fc7c263a0006d93fc8d2f2ee5cc220739829925c"
O001_HASH <- "f3f557ceb7e6baf3c0f3f2e076e22523d338e72d2522e08cab6a38266ea32147"

for (d in c("matched_null", "depth", "fdr", "failures", "validation")) {
  dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
}

write_csv_once <- function(x, path) {
  if (file.exists(path)) stop("create-once refusal: ", path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  invisible(path)
}
write_text_once <- function(lines, path) {
  if (file.exists(path)) stop("create-once refusal: ", path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path, useBytes = TRUE)
  invisible(path)
}
sha_file <- function(path) digest::digest(path, file = TRUE, algo = "sha256")
bool <- function(x) as.logical(toupper(as.character(x)))
csv_ok <- function(path) tryCatch({ read.csv(path, nrows = 5L, check.names = FALSE); TRUE }, error = function(e) FALSE)
log_msg <- function(...) message(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\t", paste(..., collapse = " "))

source(file.path(DLBCL_SUPPLEMENTARY_CODE_ROOT, "spatial", "wp3_common_v3.R"))
wp3_require_r()
.libPaths(c(file.path(REV, "06b_wp3b_spatial_scope_method_resolution/.wp3_r_library"),
            file.path(REV, "04b_stage4_environment_freeze/stage4_renv_project/renv/library/windows/R-4.5/x86_64-w64-mingw32"),
            .Library))
required <- c("Matrix", "SeuratObject", "Seurat", "BiocParallel", "UCell", "digest")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing required package(s): ", paste(missing, collapse = ","))
library(Matrix)
library(BiocParallel)

partial_path <- file.path(OUT, "matched_null/MATCHED_RANDOM_SET_MEMBERSHIP.csv")
original_script <- file.path(P06, "01_execution_scripts/04_06p_scientific_execution.R")
partition_manifest_path <- file.path(OUT, "matched_null/MATCHED_RANDOM_SET_MEMBERSHIP_PARTITION_MANIFEST.csv")
membership_validation_path <- file.path(OUT, "matched_null/MATCHED_RANDOM_SET_MEMBERSHIP_VALIDATION.csv")
batch_manifest_path <- file.path(RUN, "06P_BATCH_MANIFEST.csv")
spot_depth_path <- file.path(OUT, "gene_summaries/SPOT_DEPTH_COVARIATES.csv")
gene_summary_path <- file.path(OUT, "gene_summaries/AREA_GENE_MATCHING_SUMMARIES.csv")

if (!file.exists(file.path(ATT4, "06P_CONTINUATION_AUTHORIZATION.md"))) stop("Continuation authorization missing")
if (!identical(readLines(file.path(P06, "00_execution_authorization/06P_EXECUTION_TOKEN.txt"), warn = FALSE), TOKEN)) stop("Execution token mismatch")
if (sha_file(partial_path) != PARTIAL_HASH) stop("Interrupted membership file changed")
if (sha_file(original_script) != ORIGINAL_SCRIPT_HASH) stop("Original scientific script changed")
if (sha_file(file.path(RUN, "06P_INPUT_SHA256_REGISTRY.csv")) != O001_HASH) stop("O001 changed")
amd4 <- read.csv(file.path(PROTO, "amendments/AMENDMENT_004_06P_PARTITIONED_MEMBERSHIP_SALVAGE_AND_RESUME/06P_AMENDMENT_004_VALIDATION.csv"), check.names = FALSE)
if (!all(amd4$status == "PASS")) stop("Amendment 004 validation is not PASS")
membership_validation <- read.csv(membership_validation_path, check.names = FALSE)
if (!all(membership_validation$status == "PASS")) stop("Partition membership validation is not PASS")

blocked_outputs <- c("matched_null/MATCHING_QC.csv", "matched_null/MATCHED_NULL_MORAN_GEARY.csv",
  "matched_null/MATCHED_NULL_EMPIRICAL_TESTS.csv", "depth/DEPTH_MODEL_DIAGNOSTICS.csv",
  "depth/DEPTH_RESIDUAL_SPATIAL_STATISTICS.csv", "depth/DEPTH_RESIDUAL_PERMUTATION_SUMMARY.csv",
  "fdr/06P_SEPARATE_FDR_RESULTS.csv", "failures/06P_FAILURE_LOG.csv",
  "validation/06P_EXECUTION_VALIDATION.csv", "validation/06P_OUTPUT_MANIFEST.csv", "06P_FINAL_REPORT.md")
preexisting <- blocked_outputs[file.exists(file.path(OUT, blocked_outputs))]
if (length(preexisting)) stop("Downstream create-once output already exists: ", paste(preexisting, collapse = ";"))

inputs <- read.csv(file.path(PROTO, "06P_INPUT_REGISTRY.csv"), check.names = FALSE)
seeds <- read.csv(file.path(PROTO, "06P_SEED_REGISTRY.csv"), check.names = FALSE)
program_contract <- wp3_read_program_contract()
program_ids <- program_contract$program_id
program_names <- setNames(program_contract$program_name, program_contract$program_id)
program_genes <- setNames(program_contract$genes, program_contract$program_id)
program_order_by_id <- setNames(program_contract$program_order, program_contract$program_id)
canonical_121 <- unique(unlist(program_genes))
observed_authority <- read.csv(file.path(REV, "06c_wp3_real_spatial_continuous_analysis/continuation_v3/finalization_v2/WP3_FINAL_MORAN_GEARY_AUTHORITY.csv"), check.names = FALSE)
spot_depth <- read.csv(spot_depth_path, check.names = FALSE)
gene_summary <- read.csv(gene_summary_path, check.names = FALSE)
partition_manifest <- read.csv(partition_manifest_path, check.names = FALSE)
batch_manifest <- read.csv(batch_manifest_path, check.names = FALSE)
if (nrow(spot_depth) != 31030L || anyDuplicated(paste(spot_depth$capture_area_id, spot_depth$barcode))) stop("Spot-depth reusable output failed validation")
if (nrow(gene_summary) != 162315L || length(unique(gene_summary$capture_area_id)) != 9L) stop("Gene-summary reusable output failed validation")
if (nrow(partition_manifest) != 487L || !all(partition_manifest$validation_status == "PASS")) stop("Partition manifest failed validation")
if (nrow(batch_manifest) != 487L || !all(batch_manifest$validation_status == "PASS")) stop("Batch manifest failed validation")

role_family_for <- function(area_id, program_order, analysis_type) {
  x <- seeds[seeds$capture_area_id == area_id & seeds$program_order == program_order &
               seeds$analysis_type == analysis_type, , drop = FALSE]
  u <- unique(x$fdr_family)
  if (length(u) != 1L) stop("Nonunique role family for ", area_id, " program_order ", program_order)
  u
}

normalize_lognormalize <- function(counts) {
  totals <- Matrix::colSums(counts)
  norm <- counts
  norm@x <- log1p(norm@x / rep.int(totals, diff(norm@p)) * 10000)
  norm
}
score_sets <- function(norm_matrix, feature_sets, rank_cache = NULL) {
  UCell::ScoreSignatures_UCell(matrix = norm_matrix, features = feature_sets, precalc.ranks = rank_cache,
                               maxRank = 1500, name = "", missing_genes = "impute",
                               BPPARAM = BiocParallel::SerialParam())
}

log_msg("READ_MEMBERSHIP_PARTITIONS_START")
membership <- do.call(rbind, lapply(partition_manifest$output_file, read.csv, check.names = FALSE))
membership$replicate_id <- as.integer(membership$replicate_id)
membership$seed <- as.integer(membership$seed)
valid_membership <- membership[membership$replicate_validity == "VALID", , drop = FALSE]
final_keys <- unique(membership[, c("capture_area_id", "program_id", "replicate_id", "seed", "role_family")])
final_keys$key <- paste(final_keys$capture_area_id, final_keys$program_id, final_keys$replicate_id, sep = "||")
expected <- seeds[seeds$analysis_type == "MATCHED_GENE_SET_REPLICATE", , drop = FALSE]
expected$program_id <- program_ids[match(expected$program_order, program_contract$program_order)]
expected$replicate_id <- as.integer(expected$replicate_or_endpoint)
expected$key <- paste(expected$capture_area_id, expected$program_id, expected$replicate_id, sep = "||")
seed_expected <- setNames(as.integer(expected$derived_seed), expected$key)
if (nrow(final_keys) != 54000L || !setequal(final_keys$key, expected$key) || anyDuplicated(final_keys$key)) stop("Final membership key validation failed")
if (!all(final_keys$seed == seed_expected[final_keys$key])) stop("Final membership seed validation failed")
if (any(membership$capture_area_id == "GSM8500534_Cap.area1_LN_V1" &
        membership$program_id == "antigen_presentation" &
        membership$observed_canonical_gene == "HLA-DRB1")) stop("Area1 HLA-DRB1 exclusion failed")
if (any(valid_membership$selected_matched_gene %in% canonical_121)) stop("Canonical selected gene violation")
log_msg("READ_MEMBERSHIP_PARTITIONS_COMPLETE", "rows", nrow(membership), "keys", nrow(final_keys))

log_msg("LOAD_RAW_AUTHORITIES_START")
areas <- list(); score_checks <- list(); spatial_checks <- list()
for (i in seq_len(nrow(inputs))) {
  area_id <- inputs$capture_area_id[i]
  log_msg("LOAD_AREA", area_id)
  structure <- wp3_validate_area_structure(area_id)
  counts <- wp3_read_real_counts(structure)
  score_authority <- read.csv(inputs$primary_score_path[i], check.names = FALSE)
  counts <- counts[, score_authority$barcode, drop = FALSE]
  norm_counts <- normalize_lognormalize(counts)
  edge_index <- wp3_edge_index(score_authority$barcode, read.csv(inputs$adjacency_path[i], check.names = FALSE))
  eligibility <- read.csv(inputs$eligibility_path[i], check.names = FALSE)
  gene_sets <- setNames(strsplit(eligibility$detected_genes, ";", fixed = TRUE), eligibility$program_id)[program_ids]
  scores_repro <- wp3_primary_ucell(counts, gene_sets, area_id)$scores
  names(scores_repro) <- sub("_UCell$", "", names(scores_repro))
  scores_repro <- scores_repro[match(score_authority$barcode, scores_repro$barcode), , drop = FALSE]
  for (pid in program_ids) {
    maxdiff <- max(abs(score_authority[[pid]] - scores_repro[[pid]]), na.rm = TRUE)
    score_checks[[length(score_checks) + 1L]] <- data.frame(capture_area_id = area_id, program_id = pid,
      max_abs_difference = maxdiff, tolerance = 1e-12, status = if (is.finite(maxdiff) && maxdiff <= 1e-12) "PASS" else "FAIL")
    raw_m <- wp3_moran_edge(score_authority[[pid]], edge_index)
    raw_g <- wp3_geary_edge(score_authority[[pid]], edge_index)
    auth <- observed_authority[observed_authority$capture_area_id == area_id & observed_authority$program_id == pid, , drop = FALSE]
    spatial_checks[[length(spatial_checks) + 1L]] <- data.frame(capture_area_id = area_id, program_id = pid,
      reproduced_Moran_I = raw_m, authority_Moran_I = auth$Moran_I, reproduced_Geary_C = raw_g, authority_Geary_C = auth$Geary_C,
      moran_abs_difference = abs(raw_m - auth$Moran_I), geary_abs_difference = abs(raw_g - auth$Geary_C),
      status = if (abs(raw_m - auth$Moran_I) <= 1e-12 && abs(raw_g - auth$Geary_C) <= 1e-12) "PASS" else "FAIL")
  }
  areas[[area_id]] <- list(norm = norm_counts, scores = score_authority, edge_index = edge_index, area_role = inputs$area_role[i])
}
score_checks <- do.call(rbind, score_checks)
spatial_checks <- do.call(rbind, spatial_checks)
write_csv_once(score_checks, file.path(OUT, "gene_summaries/OBSERVED_SCORE_AUTHORITY_CHECK_RESUME.csv"))
write_csv_once(spatial_checks, file.path(OUT, "gene_summaries/OBSERVED_SPATIAL_STATISTIC_AUTHORITY_CHECK_RESUME.csv"))
if (!all(score_checks$status == "PASS") || !all(spatial_checks$status == "PASS")) stop("Observed authority gate failed")

log_msg("MATCHED_NULL_ALL_START")
matching_qc_rows <- list(); matched_null_rows <- list(); empirical_rows <- list(); failure_rows <- list()
for (area_id in names(areas)) {
  area <- areas[[area_id]]
  rank_cache <- UCell::StoreRankings_UCell(area$norm, maxRank = 1500, BPPARAM = BiocParallel::SerialParam())
  for (pid in program_ids) {
    sub <- membership[membership$capture_area_id == area_id & membership$program_id == pid, , drop = FALSE]
    sub_valid <- sub[sub$replicate_validity == "VALID", , drop = FALSE]
    rep_ids <- sort(unique(sub_valid$replicate_id))
    role <- unique(sub$role_family)
    if (length(role) != 1L) stop("Nonunique membership role family")
    log_msg("MATCHED_NULL_START", area_id, pid, "sets", length(rep_ids))
    replicate_stats <- do.call(rbind, lapply(split(sub, sub$replicate_id), function(z) {
      data.frame(replicate_id = z$replicate_id[1], seed = z$seed[1], set_size = nrow(z),
        mean_matching_distance = mean(z$matching_distance, na.rm = TRUE),
        max_matching_distance = max(z$matching_distance, na.rm = TRUE),
        fallback_50_rows = sum(z$pool_size_used == 50, na.rm = TRUE),
        fallback_100_rows = sum(z$pool_size_used == 100, na.rm = TRUE),
        fallback_200_rows = sum(z$pool_size_used == 200, na.rm = TRUE),
        replicate_validity = if (all(z$replicate_validity == "VALID")) "VALID" else "MATCH_FAILURE")
    }))
    valid_count <- length(rep_ids)
    matching_qc_rows[[length(matching_qc_rows) + 1L]] <- data.frame(
      capture_area_id = area_id, role_family = role, program_id = pid,
      attempted_sets = 1000L, valid_sets = valid_count, failed_sets = 1000L - valid_count,
      valid_fraction = valid_count / 1000, mean_matching_distance = mean(replicate_stats$mean_matching_distance, na.rm = TRUE),
      max_matching_distance = max(replicate_stats$max_matching_distance, na.rm = TRUE),
      fallback_50_rows = sum(replicate_stats$fallback_50_rows, na.rm = TRUE),
      fallback_100_rows = sum(replicate_stats$fallback_100_rows, na.rm = TRUE),
      fallback_200_rows = sum(replicate_stats$fallback_200_rows, na.rm = TRUE),
      null_validity_decision = if (valid_count >= 950L) "VALID" else "FAILED_VALID_SET_THRESHOLD")
    if (valid_count < 950L) {
      failure_rows[[length(failure_rows) + 1L]] <- data.frame(capture_area_id = area_id, program_id = pid,
        stage = "MATCHED_NULL", status = "FAILED_ENDPOINT", reason = paste("valid sets", valid_count, "< 950"))
      next
    }
    feature_sets <- split(sub_valid$selected_matched_gene, sub_valid$replicate_id)
    score_mat <- score_sets(area$norm, feature_sets, rank_cache)
    if (!is.data.frame(score_mat)) score_mat <- as.data.frame(score_mat, check.names = FALSE)
    if (nrow(score_mat) != ncol(area$norm)) score_mat <- as.data.frame(t(as.matrix(score_mat)), check.names = FALSE)
    names(score_mat) <- sub("_UCell$", "", names(score_mat))
    rownames(score_mat) <- colnames(area$norm)
    start_index <- length(matched_null_rows)
    for (rid in names(feature_sets)) {
      vals <- score_mat[[rid]]
      matched_null_rows[[length(matched_null_rows) + 1L]] <- data.frame(
        capture_area_id = area_id, role_family = role, program_id = pid,
        replicate_id = as.integer(rid), matched_set_seed = unique(sub_valid$seed[sub_valid$replicate_id == as.integer(rid)]),
        set_size = length(feature_sets[[rid]]), Moran_I = wp3_moran_edge(vals, area$edge_index),
        Geary_C = wp3_geary_edge(vals, area$edge_index), validity_status = "VALID")
    }
    null_df <- do.call(rbind, matched_null_rows[(start_index + 1L):length(matched_null_rows)])
    auth <- observed_authority[observed_authority$capture_area_id == area_id & observed_authority$program_id == pid, , drop = FALSE]
    empirical_rows[[length(empirical_rows) + 1L]] <- data.frame(
      capture_area_id = area_id, role_family = role, program_id = pid,
      observed_Moran_I = auth$Moran_I, observed_Geary_C = auth$Geary_C, valid_null_count = nrow(null_df),
      null_Moran_mean = mean(null_df$Moran_I), null_Moran_median = median(null_df$Moran_I), null_Moran_sd = sd(null_df$Moran_I),
      null_Moran_q025 = as.numeric(quantile(null_df$Moran_I, 0.025)), null_Moran_q975 = as.numeric(quantile(null_df$Moran_I, 0.975)),
      Moran_observed_percentile = mean(null_df$Moran_I <= auth$Moran_I),
      Moran_empirical_p = (1 + sum(null_df$Moran_I >= auth$Moran_I)) / (1 + nrow(null_df)),
      Moran_standardized_effect = if (sd(null_df$Moran_I) > 0) (auth$Moran_I - mean(null_df$Moran_I)) / sd(null_df$Moran_I) else NA_real_,
      null_Geary_mean = mean(null_df$Geary_C), null_Geary_median = median(null_df$Geary_C), null_Geary_sd = sd(null_df$Geary_C),
      null_Geary_q025 = as.numeric(quantile(null_df$Geary_C, 0.025)), null_Geary_q975 = as.numeric(quantile(null_df$Geary_C, 0.975)),
      Geary_observed_percentile = mean(null_df$Geary_C <= auth$Geary_C),
      Geary_empirical_p = (1 + sum(null_df$Geary_C <= auth$Geary_C)) / (1 + nrow(null_df)),
      Geary_standardized_effect = if (sd(null_df$Geary_C) > 0) (mean(null_df$Geary_C) - auth$Geary_C) / sd(null_df$Geary_C) else NA_real_,
      null_validity_status = "VALID")
    log_msg("MATCHED_NULL_COMPLETE", area_id, pid, "valid", valid_count)
  }
}
matching_qc <- do.call(rbind, matching_qc_rows)
matched_null <- do.call(rbind, matched_null_rows)
empirical <- do.call(rbind, empirical_rows)
write_csv_once(matching_qc, file.path(OUT, "matched_null/MATCHING_QC.csv"))
write_csv_once(matched_null, file.path(OUT, "matched_null/MATCHED_NULL_MORAN_GEARY.csv"))
write_csv_once(empirical, file.path(OUT, "matched_null/MATCHED_NULL_EMPIRICAL_TESTS.csv"))

log_msg("DEPTH_ALL_START")
depth_model_rows <- list(); depth_spatial_rows <- list(); depth_perm_rows <- list()
for (area_id in names(areas)) {
  area <- areas[[area_id]]
  depth <- spot_depth[spot_depth$capture_area_id == area_id, , drop = FALSE]
  for (pid in program_ids) {
    program_order <- program_order_by_id[[pid]]
    combo_role <- role_family_for(area_id, program_order, "DEPTH_RESIDUAL_MORAN_PERMUTATION_ENDPOINT")
    dat <- data.frame(score = area$scores[[pid]], nCount = depth$nCount_Spatial, nFeature = depth$nFeature_Spatial)
    fit <- lm(score ~ log1p(nCount) + log1p(nFeature), data = dat)
    X <- model.matrix(fit); rank <- qr(X)$rank; status <- if (rank < ncol(X)) "RANK_DEFICIENT" else "VALID"
    res <- residuals(fit); co <- summary(fit)$coefficients
    getc <- function(row, col) if (row %in% rownames(co)) co[row, col] else NA_real_
    depth_model_rows[[length(depth_model_rows) + 1L]] <- data.frame(
      capture_area_id = area_id, role_family = combo_role, program_id = pid,
      intercept = getc("(Intercept)", "Estimate"), beta_log1p_nCount = getc("log1p(nCount)", "Estimate"),
      beta_log1p_nFeature = getc("log1p(nFeature)", "Estimate"), se_log1p_nCount = getc("log1p(nCount)", "Std. Error"),
      se_log1p_nFeature = getc("log1p(nFeature)", "Std. Error"), t_log1p_nCount = getc("log1p(nCount)", "t value"),
      t_log1p_nFeature = getc("log1p(nFeature)", "t value"), p_log1p_nCount = getc("log1p(nCount)", "Pr(>|t|)"),
      p_log1p_nFeature = getc("log1p(nFeature)", "Pr(>|t|)"), r_squared = summary(fit)$r.squared,
      adjusted_r_squared = summary(fit)$adj.r.squared, residual_sd = sigma(fit),
      raw_cor_nCount = cor(dat$score, log1p(dat$nCount)), raw_cor_nFeature = cor(dat$score, log1p(dat$nFeature)),
      residual_cor_nCount = cor(res, log1p(dat$nCount)), residual_cor_nFeature = cor(res, log1p(dat$nFeature)),
      depth_covariate_correlation = cor(log1p(dat$nCount), log1p(dat$nFeature)),
      design_matrix_rank = rank, condition_number = kappa(X), residual_complete = all(is.finite(res)), model_status = status)
    if (status == "VALID") {
      auth <- observed_authority[observed_authority$capture_area_id == area_id & observed_authority$program_id == pid, , drop = FALSE]
      sm <- as.integer(seeds[seeds$capture_area_id == area_id & seeds$program_order == program_order & seeds$analysis_type == "DEPTH_RESIDUAL_MORAN_PERMUTATION_ENDPOINT", "derived_seed"])
      sg <- as.integer(seeds[seeds$capture_area_id == area_id & seeds$program_order == program_order & seeds$analysis_type == "DEPTH_RESIDUAL_GEARY_PERMUTATION_ENDPOINT", "derived_seed"])
      log_msg("DEPTH_ENDPOINT_START", area_id, pid)
      res_m <- wp3_moran_edge(res, area$edge_index); res_g <- wp3_geary_edge(res, area$edge_index)
      pm <- wp3_permutation_test_edge(res, area$edge_index, "moran", 9999L, sm)
      pg <- wp3_permutation_test_edge(res, area$edge_index, "geary", 9999L, sg)
      depth_spatial_rows[[length(depth_spatial_rows) + 1L]] <- data.frame(
        capture_area_id = area_id, role_family = combo_role, program_id = pid,
        raw_Moran_I = auth$Moran_I, residual_Moran_I = res_m, residual_Moran_empirical_p = pm$permutation_p,
        Moran_absolute_change = res_m - auth$Moran_I, Moran_relative_change = if (auth$Moran_I != 0) (res_m - auth$Moran_I) / abs(auth$Moran_I) else NA_real_,
        Moran_effect_retention_ratio = if (auth$Moran_I > 0) res_m / auth$Moran_I else NA_real_,
        raw_Geary_C = auth$Geary_C, residual_Geary_C = res_g, residual_Geary_empirical_p = pg$permutation_p,
        Geary_absolute_change = res_g - auth$Geary_C, raw_Geary_spatial_departure_strength = 1 - auth$Geary_C,
        residual_Geary_spatial_departure_strength = 1 - res_g,
        Geary_retention = if ((1 - auth$Geary_C) > 0) (1 - res_g) / (1 - auth$Geary_C) else NA_real_,
        raw_residual_score_correlation = cor(dat$score, res), evaluation_status = "VALID")
      depth_perm_rows[[length(depth_perm_rows) + 1L]] <- data.frame(
        capture_area_id = area_id, role_family = combo_role, program_id = pid,
        endpoint = c("DEPTH_ADJUSTED_MORAN", "DEPTH_ADJUSTED_GEARY"),
        seed = c(sm, sg), permutations = 9999L, empirical_p = c(pm$permutation_p, pg$permutation_p), status = "COMPLETE")
      log_msg("DEPTH_ENDPOINT_COMPLETE", area_id, pid)
    }
  }
}
depth_models <- do.call(rbind, depth_model_rows)
depth_spatial <- do.call(rbind, depth_spatial_rows)
depth_perm <- do.call(rbind, depth_perm_rows)
write_csv_once(depth_models, file.path(OUT, "depth/DEPTH_MODEL_DIAGNOSTICS.csv"))
write_csv_once(depth_spatial, file.path(OUT, "depth/DEPTH_RESIDUAL_SPATIAL_STATISTICS.csv"))
write_csv_once(depth_perm, file.path(OUT, "depth/DEPTH_RESIDUAL_PERMUTATION_SUMMARY.csv"))

log_msg("FDR_START")
expected_counts <- data.frame(role_family = c("PRIMARY_DLBCL", "EXPLORATORY_ANTIGEN", "CONTEXT_ONLY"),
                              expected_test_count = c(26L, 6L, 22L))
bh <- function(p) { out <- rep(NA_real_, length(p)); idx <- which(is.finite(p)); out[idx] <- p.adjust(p[idx], "BH"); out }
fdr_rows <- list()
add_fdr <- function(endpoint, table, p_col) {
  for (role in expected_counts$role_family) {
    sub <- table[table$role_family == role, , drop = FALSE]
    fdr_rows[[length(fdr_rows) + 1L]] <<- data.frame(
      endpoint_type = endpoint, role_family = role, family_label = paste("06P", endpoint, role, "BH", sep = "_"),
      expected_test_count = expected_counts$expected_test_count[expected_counts$role_family == role],
      observed_valid_test_count = sum(is.finite(sub[[p_col]])),
      capture_area_id = sub$capture_area_id, program_id = sub$program_id,
      raw_p = sub[[p_col]], BH_FDR = bh(sub[[p_col]]),
      evaluability = ifelse(is.finite(sub[[p_col]]), "EVALUABLE", "NOT_EVALUABLE"),
      exclusion_failure_reason = "")
  }
}
add_fdr("MATCHED_NULL_MORAN", empirical, "Moran_empirical_p")
add_fdr("MATCHED_NULL_GEARY", empirical, "Geary_empirical_p")
add_fdr("DEPTH_ADJUSTED_MORAN", depth_spatial, "residual_Moran_empirical_p")
add_fdr("DEPTH_ADJUSTED_GEARY", depth_spatial, "residual_Geary_empirical_p")
fdr_out <- do.call(rbind, fdr_rows)
write_csv_once(fdr_out, file.path(OUT, "fdr/06P_SEPARATE_FDR_RESULTS.csv"))
failure_log <- if (length(failure_rows)) do.call(rbind, failure_rows) else data.frame(capture_area_id = character(), program_id = character(), stage = character(), status = character(), reason = character())
write_csv_once(failure_log, file.path(OUT, "failures/06P_FAILURE_LOG.csv"))

scientific_csvs <- c(spot_depth_path, gene_summary_path, partition_manifest_path, membership_validation_path,
  file.path(OUT, "matched_null/MATCHING_QC.csv"), file.path(OUT, "matched_null/MATCHED_NULL_MORAN_GEARY.csv"),
  file.path(OUT, "matched_null/MATCHED_NULL_EMPIRICAL_TESTS.csv"), file.path(OUT, "depth/DEPTH_MODEL_DIAGNOSTICS.csv"),
  file.path(OUT, "depth/DEPTH_RESIDUAL_SPATIAL_STATISTICS.csv"), file.path(OUT, "depth/DEPTH_RESIDUAL_PERMUTATION_SUMMARY.csv"),
  file.path(OUT, "fdr/06P_SEPARATE_FDR_RESULTS.csv"), file.path(OUT, "failures/06P_FAILURE_LOG.csv"))
fdr_counts <- aggregate(raw_p ~ endpoint_type + role_family, fdr_out, length)
names(fdr_counts)[3] <- "observed_count"
expected_fdr_strings <- as.vector(outer(c("MATCHED_NULL_MORAN", "MATCHED_NULL_GEARY", "DEPTH_ADJUSTED_MORAN", "DEPTH_ADJUSTED_GEARY"),
                                        c("PRIMARY_DLBCL 26", "EXPLORATORY_ANTIGEN 6", "CONTEXT_ONLY 22"), paste))
git_tracked <- system2("git", c("-C", shQuote(REV), "diff", "--name-only"), stdout = TRUE)
git_staged <- system2("git", c("-C", shQuote(REV), "diff", "--cached", "--name-only"), stdout = TRUE)
validation <- data.frame(
  check = c("interrupted_membership_hash_unchanged", "original_script_hash_unchanged", "amendment004_pass",
            "membership_validation_pass", "salvage_5999_replicates", "salvage_129978_rows", "replacement_48001",
            "final_54000_keys", "one_thousand_keys_per_combination", "hla_drb1_excluded_area1_antigen",
            "exact_seed_identity", "no_canonical_selected_gene", "spot_depth_31030", "gene_summary_162315",
            "matching_qc_complete", "all_54_depth_models", "all_valid_endpoints_complete", "permutations_9999",
            "exactly_12_fdr_families", "expected_26_6_22_counts", "all_output_csvs_parse",
            "o001_unchanged", "tracked_git_changes_zero", "staged_git_changes_zero", "no_commit_or_push"),
  status = c(sha_file(partial_path) == PARTIAL_HASH, sha_file(original_script) == ORIGINAL_SCRIPT_HASH,
             all(amd4$status == "PASS"), all(membership_validation$status == "PASS"),
             as.integer(membership_validation$observed[membership_validation$check_id == "SALVAGED_REPLICATE_COUNT"]) == 5999L,
             as.integer(membership_validation$observed[membership_validation$check_id == "SALVAGED_ROW_COUNT"]) == 129978L,
             as.integer(membership_validation$observed[membership_validation$check_id == "REPLACEMENT_MISSING_COUNT"]) == 48001L,
             nrow(final_keys) == 54000L, all(table(paste(final_keys$capture_area_id, final_keys$program_id)) == 1000L),
             !any(membership$capture_area_id == "GSM8500534_Cap.area1_LN_V1" & membership$program_id == "antigen_presentation" & membership$observed_canonical_gene == "HLA-DRB1"),
             all(final_keys$seed == seed_expected[final_keys$key]), !any(valid_membership$selected_matched_gene %in% canonical_121),
             nrow(spot_depth) == 31030L, nrow(gene_summary) == 162315L, nrow(matching_qc) == 54L,
             nrow(depth_models) == 54L, nrow(depth_perm) == 108L && all(depth_perm$status == "COMPLETE"),
             all(depth_perm$permutations == 9999L), length(unique(paste(fdr_out$endpoint_type, fdr_out$role_family))) == 12L,
             setequal(paste(fdr_counts$endpoint_type, fdr_counts$role_family, fdr_counts$observed_count), expected_fdr_strings),
             all(vapply(scientific_csvs, csv_ok, logical(1))), sha_file(file.path(RUN, "06P_INPUT_SHA256_REGISTRY.csv")) == O001_HASH,
             length(git_tracked) == 0L, length(git_staged) == 0L, TRUE),
  detail = c(PARTIAL_HASH, ORIGINAL_SCRIPT_HASH, paste(sum(amd4$status == "PASS"), "of", nrow(amd4)),
             paste(sum(membership_validation$status == "PASS"), "of", nrow(membership_validation)),
             membership_validation$observed[membership_validation$check_id == "SALVAGED_REPLICATE_COUNT"],
             membership_validation$observed[membership_validation$check_id == "SALVAGED_ROW_COUNT"],
             membership_validation$observed[membership_validation$check_id == "REPLACEMENT_MISSING_COUNT"],
             nrow(final_keys), paste(range(table(paste(final_keys$capture_area_id, final_keys$program_id))), collapse = "-"),
             "0", "derived_seed equality", "valid rows checked", nrow(spot_depth), nrow(gene_summary),
             nrow(matching_qc), nrow(depth_models), nrow(depth_perm), paste(unique(depth_perm$permutations), collapse = ";"),
             length(unique(paste(fdr_out$endpoint_type, fdr_out$role_family))),
             paste(paste(fdr_counts$endpoint_type, fdr_counts$role_family, fdr_counts$observed_count, sep = ":"), collapse = ";"),
             paste(sum(vapply(scientific_csvs, csv_ok, logical(1))), "of", length(scientific_csvs)),
             O001_HASH, length(git_tracked), length(git_staged), "No git commit or push performed by script"),
  stringsAsFactors = FALSE)
validation$status_text <- ifelse(validation$status, "PASS", "FAIL")
authority_decision <- if (all(validation$status)) "FINAL_06P_AUTHORITY" else "VALIDATION_FAILED_NO_AUTHORITY"
write_csv_once(validation, file.path(OUT, "validation/06P_EXECUTION_VALIDATION.csv"))
write_text_once(c("# 06p Final Report", "",
                  paste("Authority decision:", authority_decision),
                  paste("Amendment 004 status: PASS"),
                  paste("Salvaged replicates:", membership_validation$observed[membership_validation$check_id == "SALVAGED_REPLICATE_COUNT"]),
                  paste("Salvaged rows:", membership_validation$observed[membership_validation$check_id == "SALVAGED_ROW_COUNT"]),
                  paste("Replacement/missing replicates:", membership_validation$observed[membership_validation$check_id == "REPLACEMENT_MISSING_COUNT"]),
                  paste("Final composite keys:", nrow(final_keys)),
                  paste("Total retained spots:", nrow(spot_depth)),
                  paste("Area-gene summary rows:", nrow(gene_summary)),
                  paste("Matched-null valid sets:", sum(matching_qc$valid_sets)),
                  paste("Matching failures:", sum(matching_qc$failed_sets)),
                  paste("Depth models:", nrow(depth_models)),
                  paste("Depth endpoints:", nrow(depth_perm)),
                  paste("FDR families:", length(unique(paste(fdr_out$endpoint_type, fdr_out$role_family)))),
                  "",
                  "Negative, failed, and heterogeneous findings are retained in the scientific CSV outputs.",
                  "No causal spatial organization, cell-cell communication, patient-level replication, independent normal-control validation, complete technical-confounding-removal, or DLBCL-specific context-area claim is made.",
                  paste("Final authority status:", authority_decision)),
                file.path(OUT, "06P_FINAL_REPORT.md"))
manifest_paths <- c(file.path(RUN, "06P_INPUT_SHA256_REGISTRY.csv"), batch_manifest_path,
                    file.path(RUN, "membership_salvage_001/06P_INTERRUPTED_MEMBERSHIP_FILE_AUDIT.csv"),
                    file.path(RUN, "membership_salvage_001/06P_REPLICATE_LEVEL_SALVAGE_AUDIT.csv"),
                    file.path(RUN, "membership_salvage_001/06P_SALVAGED_REPLICATE_KEY_REGISTRY.csv"),
                    file.path(RUN, "membership_salvage_001/06P_REPLACEMENT_AND_MISSING_REPLICATE_PLAN.csv"),
                    file.path(RUN, "membership_salvage_001/06P_MEMBERSHIP_SALVAGE_REPORT.md"),
                    scientific_csvs,
                    file.path(OUT, "validation/06P_EXECUTION_VALIDATION.csv"),
                    file.path(OUT, "validation/06P_OUTPUT_MANIFEST.csv"),
                    file.path(OUT, "06P_FINAL_REPORT.md"))
manifest_exists <- file.exists(manifest_paths)
manifest_exists[basename(manifest_paths) == "06P_OUTPUT_MANIFEST.csv"] <- TRUE
manifest <- data.frame(path = manifest_paths, exists = manifest_exists,
  size_bytes = ifelse(manifest_exists & basename(manifest_paths) != "06P_OUTPUT_MANIFEST.csv", as.numeric(file.info(manifest_paths)$size), NA_real_),
  sha256 = ifelse(manifest_exists & basename(manifest_paths) != "06P_OUTPUT_MANIFEST.csv",
                  vapply(manifest_paths, function(p) if (file.exists(p) && basename(p) != "06P_OUTPUT_MANIFEST.csv") sha_file(p) else NA_character_, character(1)), NA_character_),
  creation_origin = ifelse(basename(manifest_paths) == "06P_INPUT_SHA256_REGISTRY.csv", "ATTEMPT_001_PREFLIGHT", "AMENDMENT_004_CONTROLLED_RESUME"),
  final_execution_action = ifelse(basename(manifest_paths) == "06P_INPUT_SHA256_REGISTRY.csv", "PRESERVED_AND_REUSED", "CREATED_ONCE"),
  overwrite_performed = "FALSE",
  effective_input_authority = ifelse(basename(manifest_paths) == "06P_INPUT_SHA256_REGISTRY.csv", "81_OF_81_COMPONENT_SHA256_PASS", ""),
  stringsAsFactors = FALSE)
write_csv_once(manifest, file.path(OUT, "validation/06P_OUTPUT_MANIFEST.csv"))
log_msg("06P_PARTITIONED_RESUME_COMPLETE", authority_decision)
