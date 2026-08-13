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
LOG <- file.path(P06, "02_execution_logs")
TOKEN <- "AUTHORIZE_06P_MATCHED_NULL_AND_DEPTH_EXECUTION_SEED_20260804"

for (d in c("gene_summaries", "matched_null", "depth", "fdr", "failures", "validation")) {
  dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
}
dir.create(LOG, recursive = TRUE, showWarnings = FALSE)

write_csv_once <- function(x, path) {
  if (file.exists(path)) stop("create-once refusal: ", path)
  write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  invisible(path)
}
write_text_once <- function(lines, path) {
  if (file.exists(path)) stop("create-once refusal: ", path)
  writeLines(lines, path, useBytes = TRUE)
  invisible(path)
}
sha_file <- function(path) digest::digest(path, file = TRUE, algo = "sha256")
fmt_num <- function(x) ifelse(is.finite(x), format(x, digits = 16, scientific = TRUE), NA_character_)

if (!file.exists(file.path(ATT4, "06P_CONTINUATION_AUTHORIZATION.md"))) {
  stop("Attempt 004 continuation authorization is absent")
}
if (!identical(readLines(file.path(P06, "00_execution_authorization/06P_EXECUTION_TOKEN.txt"), warn = FALSE), TOKEN)) {
  stop("Execution token mismatch")
}
if (file.exists(file.path(OUT, "06P_FINAL_REPORT.md"))) stop("Final report already exists")
if (!file.exists(file.path(ATT4, "06P_PREFLIGHT_REPORT_ATTEMPT004_ERRATUM.md"))) {
  stop("Attempt 004 preflight report erratum is absent")
}

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

inputs <- read.csv(file.path(PROTO, "06P_INPUT_REGISTRY.csv"), check.names = FALSE)
seeds <- read.csv(file.path(PROTO, "06P_SEED_REGISTRY.csv"), check.names = FALSE)
fdr_contract <- read.csv(file.path(PROTO, "06P_FDR_FAMILY_CONTRACT.csv"), check.names = FALSE)
program_contract <- wp3_read_program_contract()
program_ids <- program_contract$program_id
program_names <- setNames(program_contract$program_name, program_contract$program_id)
program_genes <- setNames(program_contract$genes, program_contract$program_id)
canonical_121 <- unique(unlist(program_genes))
observed_authority <- read.csv(file.path(REV, "06c_wp3_real_spatial_continuous_analysis/continuation_v3/finalization_v2/WP3_FINAL_MORAN_GEARY_AUTHORITY.csv"), check.names = FALSE)

runtime_lines <- capture.output({
  cat("06p runtime and package capture\n")
  cat("Started:", format(Sys.time(), "%Y-%m-%dT%H:%M:%OS%z"), "\n")
  cat("R.version:", paste(R.version$major, R.version$minor, sep = "."), "\n")
  cat("Platform:", R.version$platform, "\n")
  cat("Bioconductor:", if (requireNamespace("BiocManager", quietly = TRUE)) as.character(BiocManager::version()) else "NOT_AVAILABLE", "\n")
  ext <- extSoftVersion()
  cat("BLAS:", if ("BLAS" %in% names(ext)) ext[["BLAS"]] else "NOT_REPORTED", "\n")
  cat("LAPACK:", if ("LAPACK" %in% names(ext)) ext[["LAPACK"]] else paste(La_version(), collapse = "."), "\n")
  cat("Locale:", paste(Sys.getlocale(), collapse = ";"), "\n")
  cat("OMP_NUM_THREADS:", Sys.getenv("OMP_NUM_THREADS", unset = ""), "\n")
  for (pkg in required) cat(pkg, as.character(utils::packageVersion(pkg)), "\n")
  print(sessionInfo())
})
write_text_once(runtime_lines, file.path(RUN, "06P_RUNTIME_PACKAGE_CAPTURE.txt"))

script_paths <- c(file.path(P06, "01_execution_scripts/03_06p_attempt004_preflight.R"),
                  file.path(P06, "01_execution_scripts/04_06p_scientific_execution.R"),
                  file.path(REV, "06h_wp3_interrupted_continuation_sparse_resume/scripts/wp3_common_v3.R"))
write_csv_once(data.frame(script_path = script_paths, size_bytes = as.numeric(file.info(script_paths)$size),
                          sha256 = vapply(script_paths, sha_file, character(1)), stringsAsFactors = FALSE),
               file.path(RUN, "06P_SCRIPT_MANIFEST.csv"))
write_csv_once(seeds, file.path(RUN, "06P_SEED_RUNTIME_VERIFICATION.csv"))

normalize_lognormalize <- function(counts) {
  totals <- Matrix::colSums(counts)
  if (any(totals <= 0)) stop("Zero total count spot encountered")
  norm <- counts
  norm@x <- log1p(norm@x / rep.int(totals, diff(norm@p)) * 10000)
  norm
}

score_sets <- function(norm_matrix, feature_sets, rank_cache = NULL) {
  UCell::ScoreSignatures_UCell(matrix = norm_matrix, features = feature_sets,
                               precalc.ranks = rank_cache, maxRank = 1500,
                               name = "", missing_genes = "impute",
                               BPPARAM = BiocParallel::SerialParam())
}

role_family_for <- function(area_id, program_order, analysis_type) {
  x <- seeds[seeds$capture_area_id == area_id &
               seeds$program_order == program_order &
               seeds$analysis_type == analysis_type, , drop = FALSE]
  if (!nrow(x)) stop("No seed role-family record for ", area_id, " program_order ", program_order)
  unique_role <- unique(x$fdr_family)
  if (length(unique_role) != 1L) stop("Nonunique seed role family for ", area_id, " program_order ", program_order)
  unique_role
}

bh <- function(p) {
  out <- rep(NA_real_, length(p))
  idx <- which(is.finite(p))
  out[idx] <- p.adjust(p[idx], method = "BH")
  out
}

areas <- list()
spot_depth_rows <- list()
gene_rows <- list()
gene_area_summary_rows <- list()
observed_score_checks <- list()
observed_spatial_checks <- list()

for (i in seq_len(nrow(inputs))) {
  area_id <- inputs$capture_area_id[i]
  message("Loading area ", area_id)
  structure <- wp3_validate_area_structure(area_id)
  counts <- wp3_read_real_counts(structure)
  score_authority <- read.csv(inputs$primary_score_path[i], check.names = FALSE)
  counts <- counts[, score_authority$barcode, drop = FALSE]
  norm_counts <- normalize_lognormalize(counts)
  edge_index <- wp3_edge_index(score_authority$barcode, read.csv(inputs$adjacency_path[i], check.names = FALSE))
  eligibility <- read.csv(inputs$eligibility_path[i], check.names = FALSE)
  gene_sets <- setNames(strsplit(eligibility$detected_genes, ";", fixed = TRUE), eligibility$program_id)
  gene_sets <- gene_sets[program_ids]
  scores_repro <- wp3_primary_ucell(counts, gene_sets, area_id)$scores
  if (!"barcode" %in% names(scores_repro)) scores_repro <- data.frame(barcode = rownames(scores_repro), scores_repro, check.names = FALSE)
  names(scores_repro) <- sub("_UCell$", "", names(scores_repro))
  if (!all(score_authority$barcode %in% scores_repro$barcode)) stop("Reproduced score barcodes incomplete: ", area_id)
  scores_repro <- scores_repro[match(score_authority$barcode, scores_repro$barcode), , drop = FALSE]
  for (pid in program_ids) {
    maxdiff <- max(abs(score_authority[[pid]] - scores_repro[[pid]]), na.rm = TRUE)
    observed_score_checks[[length(observed_score_checks) + 1L]] <- data.frame(
      capture_area_id = area_id, program_id = pid, max_abs_difference = maxdiff,
      tolerance = 1e-12, status = if (is.finite(maxdiff) && maxdiff <= 1e-12) "PASS" else "FAIL",
      stringsAsFactors = FALSE)
  }
  for (pid in program_ids) {
    raw_m <- wp3_moran_edge(score_authority[[pid]], edge_index)
    raw_g <- wp3_geary_edge(score_authority[[pid]], edge_index)
    auth <- observed_authority[observed_authority$capture_area_id == area_id & observed_authority$program_id == pid, , drop = FALSE]
    observed_spatial_checks[[length(observed_spatial_checks) + 1L]] <- data.frame(
      capture_area_id = area_id, program_id = pid,
      reproduced_Moran_I = raw_m, authority_Moran_I = auth$Moran_I,
      reproduced_Geary_C = raw_g, authority_Geary_C = auth$Geary_C,
      moran_abs_difference = abs(raw_m - auth$Moran_I),
      geary_abs_difference = abs(raw_g - auth$Geary_C),
      status = if (abs(raw_m - auth$Moran_I) <= 1e-12 && abs(raw_g - auth$Geary_C) <= 1e-12) "PASS" else "FAIL",
      stringsAsFactors = FALSE)
  }
  ncount <- Matrix::colSums(counts)
  nfeature <- Matrix::colSums(counts > 0)
  spot_depth_rows[[length(spot_depth_rows) + 1L]] <- data.frame(
    capture_area_id = area_id, barcode = colnames(counts),
    nCount_Spatial = as.numeric(ncount), nFeature_Spatial = as.numeric(nfeature),
    retained_spot_status = "RETAINED_FINAL_WP3", input_object_identifier = inputs$raw_count_path[i],
    stringsAsFactors = FALSE)
  detection_n <- Matrix::rowSums(counts > 0)
  detection_fraction <- as.numeric(detection_n) / ncol(counts)
  mean_norm <- Matrix::rowMeans(norm_counts)
  second <- Matrix::rowMeans(norm_counts ^ 2)
  variance <- pmax(0, second - mean_norm^2)
  present <- rownames(counts)
  program_membership <- rep("", length(present))
  for (pid in program_ids) program_membership[present %in% program_genes[[pid]]] <- paste(program_membership[present %in% program_genes[[pid]]], pid, sep = ";")
  program_membership <- sub("^;", "", program_membership)
  gr <- data.frame(capture_area_id = area_id, gene_symbol = present,
                   present_in_area_feature_matrix = TRUE,
                   raw_detection_in_at_least_one_retained_spot = detection_n > 0,
                   spot_detection_fraction = detection_fraction,
                   mean_LogNormalized_expression = as.numeric(mean_norm),
                   log1p_mean_LogNormalized_expression = log1p(as.numeric(mean_norm)),
                   zero_variance_status = variance == 0,
                   canonical_121_membership = present %in% canonical_121,
                   canonical_program_membership = program_membership,
                   eligible_random_universe = (detection_n > 0) & is.finite(mean_norm) & is.finite(detection_fraction) & !(present %in% canonical_121),
                   stringsAsFactors = FALSE)
  gene_rows[[length(gene_rows) + 1L]] <- gr
  gene_area_summary_rows[[length(gene_area_summary_rows) + 1L]] <- data.frame(
    capture_area_id = area_id, total_features = nrow(gr), detected_features = sum(gr$raw_detection_in_at_least_one_retained_spot),
    eligible_random_universe_genes = sum(gr$eligible_random_universe),
    excluded_canonical_genes = sum(gr$canonical_121_membership),
    excluded_duplicate_or_invalid_genes = 0L,
    mean_detection_fraction = mean(gr$spot_detection_fraction),
    min_detection_fraction = min(gr$spot_detection_fraction),
    max_detection_fraction = max(gr$spot_detection_fraction),
    mean_LogNormalized_expression = mean(gr$mean_LogNormalized_expression),
    min_LogNormalized_expression = min(gr$mean_LogNormalized_expression),
    max_LogNormalized_expression = max(gr$mean_LogNormalized_expression),
    stringsAsFactors = FALSE)
  areas[[area_id]] <- list(counts = counts, norm = norm_counts, scores = score_authority,
                           eligibility = eligibility, edge_index = edge_index, gene_summary = gr,
                           area_role = inputs$area_role[i])
}

score_checks <- do.call(rbind, observed_score_checks)
spatial_checks <- do.call(rbind, observed_spatial_checks)
write_csv_once(do.call(rbind, spot_depth_rows), file.path(OUT, "gene_summaries/SPOT_DEPTH_COVARIATES.csv"))
write_csv_once(do.call(rbind, gene_rows), file.path(OUT, "gene_summaries/AREA_GENE_MATCHING_SUMMARIES.csv"))
write_csv_once(do.call(rbind, gene_area_summary_rows), file.path(OUT, "gene_summaries/AREA_GENE_MATCHING_AREA_SUMMARY.csv"))
write_csv_once(score_checks, file.path(OUT, "gene_summaries/OBSERVED_SCORE_AUTHORITY_CHECK.csv"))
write_csv_once(spatial_checks, file.path(OUT, "gene_summaries/OBSERVED_SPATIAL_STATISTIC_AUTHORITY_CHECK.csv"))

authority_pass_area <- tapply(score_checks$status == "PASS", score_checks$capture_area_id, all) &
  tapply(spatial_checks$status == "PASS", spatial_checks$capture_area_id, all)

membership_con <- file(file.path(OUT, "matched_null/MATCHED_RANDOM_SET_MEMBERSHIP.csv"), open = "wt", encoding = "UTF-8")
writeLines(paste(c("capture_area_id","role_family","program_id","program_name","replicate_id","seed","observed_canonical_gene","selected_matched_gene","matching_distance","pool_size_used","fallback_level","replicate_validity"), collapse = ","), membership_con)
on.exit(try(close(membership_con), silent = TRUE), add = TRUE)

matching_qc_rows <- list()
matched_null_rows <- list()
empirical_rows <- list()
failure_rows <- list()
batch_rows <- list()

for (area_id in names(areas)) {
  area <- areas[[area_id]]
  if (!isTRUE(authority_pass_area[[area_id]])) {
    failure_rows[[length(failure_rows) + 1L]] <- data.frame(capture_area_id = area_id, stage = "OBSERVED_AUTHORITY_GATE", status = "HARD_STOP_AREA", reason = "Observed score or spatial statistic authority check failed", stringsAsFactors = FALSE)
    next
  }
  candidates <- area$gene_summary[area$gene_summary$eligible_random_universe, , drop = FALSE]
  cand_expr <- candidates$log1p_mean_LogNormalized_expression
  cand_det <- candidates$spot_detection_fraction
  center <- c(mean(cand_expr), mean(cand_det))
  scalev <- c(stats::sd(cand_expr), stats::sd(cand_det))
  zcand <- cbind((cand_expr - center[1]) / scalev[1], (cand_det - center[2]) / scalev[2])
  rownames(zcand) <- candidates$gene_symbol
  rank_cache <- UCell::StoreRankings_UCell(area$norm, maxRank = 1500, BPPARAM = BiocParallel::SerialParam())
  for (pid in program_ids) {
    elig <- area$eligibility[area$eligibility$program_id == pid, , drop = FALSE]
    combo_role_family <- role_family_for(area_id, elig$program_order[1], "MATCHED_GENE_SET_REPLICATE")
    observed_genes <- strsplit(elig$detected_genes[1], ";", fixed = TRUE)[[1]]
    observed_genes <- observed_genes[observed_genes %in% area$gene_summary$gene_symbol]
    obs_summary <- area$gene_summary[match(observed_genes, area$gene_summary$gene_symbol), , drop = FALSE]
    ztarget <- cbind((obs_summary$log1p_mean_LogNormalized_expression - center[1]) / scalev[1],
                     (obs_summary$spot_detection_fraction - center[2]) / scalev[2])
    neighbor_list <- lapply(seq_along(observed_genes), function(k) {
      d <- sqrt((zcand[,1] - ztarget[k,1])^2 + (zcand[,2] - ztarget[k,2])^2)
      ord <- order(d, toupper(names(d)))
      data.frame(gene = names(d)[ord], distance = as.numeric(d[ord]), stringsAsFactors = FALSE)
    })
    names(neighbor_list) <- observed_genes
    seed_sub <- seeds[seeds$capture_area_id == area_id & seeds$program_order == elig$program_order[1] &
                        seeds$analysis_type == "MATCHED_GENE_SET_REPLICATE", , drop = FALSE]
    valid_features <- list()
    replicate_qc <- list()
    for (r in seq_len(nrow(seed_sub))) {
      rep_id <- as.integer(seed_sub$replicate_or_endpoint[r])
      set.seed(seed_sub$derived_seed[r])
      selected <- character()
      distances <- numeric()
      pool_used <- integer()
      fallback_used <- integer()
      valid <- TRUE
      for (g in observed_genes) {
        choice <- NA_character_; distv <- NA_real_; pool <- NA_integer_; fallback <- NA_integer_
        for (poolsz in c(50L, 100L, 200L)) {
          pool_df <- head(neighbor_list[[g]], poolsz)
          available <- pool_df[!pool_df$gene %in% selected, , drop = FALSE]
          if (nrow(available)) {
            idx <- sample.int(nrow(available), 1L)
            choice <- available$gene[idx]; distv <- available$distance[idx]
            pool <- poolsz; fallback <- match(poolsz, c(50L, 100L, 200L))
            break
          }
        }
        if (is.na(choice)) { valid <- FALSE; break }
        selected <- c(selected, choice); distances <- c(distances, distv)
        pool_used <- c(pool_used, pool); fallback_used <- c(fallback_used, fallback)
      }
      final_valid <- valid && length(selected) == length(observed_genes) &&
        length(unique(selected)) == length(selected) && !any(selected %in% canonical_121)
      if (final_valid) {
        valid_features[[as.character(rep_id)]] <- selected
      } else {
        valid <- FALSE
        selected <- rep(NA_character_, length(observed_genes)); distances <- rep(NA_real_, length(observed_genes))
        pool_used <- rep(NA_integer_, length(observed_genes)); fallback_used <- rep(NA_integer_, length(observed_genes))
      }
      lines <- paste(area_id, combo_role_family, pid, program_names[[pid]], rep_id, seed_sub$derived_seed[r],
                     observed_genes, selected, distances, pool_used, fallback_used,
                     if (valid) "VALID" else "MATCH_FAILURE", sep = ",")
      writeLines(lines, membership_con)
      smd_expr <- mean((area$gene_summary$log1p_mean_LogNormalized_expression[match(selected, area$gene_summary$gene_symbol)] - obs_summary$log1p_mean_LogNormalized_expression) / scalev[1], na.rm = TRUE)
      smd_det <- mean((area$gene_summary$spot_detection_fraction[match(selected, area$gene_summary$gene_symbol)] - obs_summary$spot_detection_fraction) / scalev[2], na.rm = TRUE)
      replicate_qc[[length(replicate_qc) + 1L]] <- data.frame(
        capture_area_id = area_id, role_family = combo_role_family, program_id = pid,
        replicate_id = rep_id, seed = seed_sub$derived_seed[r], set_size = length(observed_genes),
        mean_matching_distance = mean(distances, na.rm = TRUE), max_matching_distance = max(distances, na.rm = TRUE),
        expression_smd = smd_expr, detection_fraction_smd = smd_det,
        expansions_50_to_100 = sum(pool_used == 100L, na.rm = TRUE),
        expansions_100_to_200 = sum(pool_used == 200L, na.rm = TRUE),
        collision_count = sum(fallback_used > 1L, na.rm = TRUE),
        replicate_validity = if (valid) "VALID" else "MATCH_FAILURE",
        stringsAsFactors = FALSE)
    }
    replicate_qc_df <- do.call(rbind, replicate_qc)
    valid_count <- length(valid_features)
    combo_valid <- valid_count >= 950L
    matching_qc_rows[[length(matching_qc_rows) + 1L]] <- data.frame(
      capture_area_id = area_id, role_family = combo_role_family, program_id = pid,
      attempted_sets = 1000L, valid_sets = valid_count, failed_sets = 1000L - valid_count,
      valid_fraction = valid_count / 1000, mean_matching_distance = mean(replicate_qc_df$mean_matching_distance, na.rm = TRUE),
      max_matching_distance = max(replicate_qc_df$max_matching_distance, na.rm = TRUE),
      median_abs_expression_smd = median(abs(replicate_qc_df$expression_smd), na.rm = TRUE),
      max_abs_expression_smd = max(abs(replicate_qc_df$expression_smd), na.rm = TRUE),
      median_abs_detection_smd = median(abs(replicate_qc_df$detection_fraction_smd), na.rm = TRUE),
      max_abs_detection_smd = max(abs(replicate_qc_df$detection_fraction_smd), na.rm = TRUE),
      fallback_50 = sum(replicate_qc_df$collision_count == 0, na.rm = TRUE),
      fallback_100 = sum(replicate_qc_df$expansions_50_to_100 > 0, na.rm = TRUE),
      fallback_200 = sum(replicate_qc_df$expansions_100_to_200 > 0, na.rm = TRUE),
      null_validity_decision = if (combo_valid) "VALID" else "FAILED_VALID_SET_THRESHOLD",
      stringsAsFactors = FALSE)
    if (!combo_valid) {
      failure_rows[[length(failure_rows) + 1L]] <- data.frame(capture_area_id = area_id, stage = "MATCHED_NULL", status = "FAILED_ENDPOINT", reason = paste("valid sets", valid_count, "< 950"), stringsAsFactors = FALSE)
      next
    }
    score_mat <- score_sets(area$norm, valid_features, rank_cache)
    if (!is.data.frame(score_mat)) score_mat <- as.data.frame(score_mat, check.names = FALSE)
    if (nrow(score_mat) != ncol(area$norm)) score_mat <- as.data.frame(t(as.matrix(score_mat)), check.names = FALSE)
    names(score_mat) <- sub("_UCell$", "", names(score_mat))
    rownames(score_mat) <- colnames(area$norm)
    for (nm in names(valid_features)) {
      vals <- score_mat[[nm]]
      m <- wp3_moran_edge(vals, area$edge_index)
      g <- wp3_geary_edge(vals, area$edge_index)
      matched_null_rows[[length(matched_null_rows) + 1L]] <- data.frame(
        capture_area_id = area_id, role_family = combo_role_family, program_id = pid,
        replicate_id = as.integer(nm), matched_set_seed = seed_sub$derived_seed[match(as.integer(nm), as.integer(seed_sub$replicate_or_endpoint))],
        set_size = length(valid_features[[nm]]), Moran_I = m, Geary_C = g,
        validity_status = "VALID", stringsAsFactors = FALSE)
    }
    null_df <- do.call(rbind, tail(matched_null_rows, length(valid_features)))
    auth <- observed_authority[observed_authority$capture_area_id == area_id & observed_authority$program_id == pid, , drop = FALSE]
    moran_p <- (1 + sum(null_df$Moran_I >= auth$Moran_I, na.rm = TRUE)) / (1 + nrow(null_df))
    geary_p <- (1 + sum(null_df$Geary_C <= auth$Geary_C, na.rm = TRUE)) / (1 + nrow(null_df))
    empirical_rows[[length(empirical_rows) + 1L]] <- data.frame(
      capture_area_id = area_id, role_family = combo_role_family, program_id = pid,
      observed_Moran_I = auth$Moran_I, observed_Geary_C = auth$Geary_C, valid_null_count = nrow(null_df),
      null_Moran_mean = mean(null_df$Moran_I), null_Moran_median = median(null_df$Moran_I), null_Moran_sd = sd(null_df$Moran_I),
      null_Moran_q025 = quantile(null_df$Moran_I, 0.025), null_Moran_q975 = quantile(null_df$Moran_I, 0.975),
      Moran_observed_percentile = mean(null_df$Moran_I <= auth$Moran_I), Moran_empirical_p = moran_p,
      Moran_standardized_effect = if (sd(null_df$Moran_I) > 0) (auth$Moran_I - mean(null_df$Moran_I)) / sd(null_df$Moran_I) else NA_real_,
      null_Geary_mean = mean(null_df$Geary_C), null_Geary_median = median(null_df$Geary_C), null_Geary_sd = sd(null_df$Geary_C),
      null_Geary_q025 = quantile(null_df$Geary_C, 0.025), null_Geary_q975 = quantile(null_df$Geary_C, 0.975),
      Geary_observed_percentile = mean(null_df$Geary_C <= auth$Geary_C), Geary_empirical_p = geary_p,
      Geary_standardized_effect = if (sd(null_df$Geary_C) > 0) (mean(null_df$Geary_C) - auth$Geary_C) / sd(null_df$Geary_C) else NA_real_,
      null_validity_status = "VALID", stringsAsFactors = FALSE)
    batch_rows[[length(batch_rows) + 1L]] <- data.frame(capture_area_id = area_id, program_id = pid, batch_type = "MATCHED_NULL", replicate_start = 1L, replicate_end = 1000L, valid_sets = valid_count, status = "COMPLETE", stringsAsFactors = FALSE)
  }
}
close(membership_con)
on.exit(NULL, add = FALSE)

matched_null <- if (length(matched_null_rows)) do.call(rbind, matched_null_rows) else data.frame()
matching_qc <- if (length(matching_qc_rows)) do.call(rbind, matching_qc_rows) else data.frame()
empirical <- if (length(empirical_rows)) do.call(rbind, empirical_rows) else data.frame()
write_csv_once(matching_qc, file.path(OUT, "matched_null/MATCHING_QC.csv"))
write_csv_once(matched_null, file.path(OUT, "matched_null/MATCHED_NULL_MORAN_GEARY.csv"))
write_csv_once(empirical, file.path(OUT, "matched_null/MATCHED_NULL_EMPIRICAL_TESTS.csv"))

depth_model_rows <- list()
depth_spatial_rows <- list()
depth_perm_rows <- list()
for (area_id in names(areas)) {
  area <- areas[[area_id]]
  depth <- do.call(rbind, spot_depth_rows)[do.call(rbind, spot_depth_rows)$capture_area_id == area_id, ]
  for (pid in program_ids) {
    program_order <- area$eligibility$program_order[area$eligibility$program_id == pid]
    combo_role_family <- role_family_for(area_id, program_order, "DEPTH_RESIDUAL_MORAN_PERMUTATION_ENDPOINT")
    dat <- data.frame(score = area$scores[[pid]], nCount = depth$nCount_Spatial, nFeature = depth$nFeature_Spatial)
    fit <- lm(score ~ log1p(nCount) + log1p(nFeature), data = dat)
    X <- model.matrix(fit)
    rank <- qr(X)$rank
    cond <- kappa(X)
    status <- if (rank < ncol(X)) "RANK_DEFICIENT" else "VALID"
    res <- residuals(fit)
    co <- summary(fit)$coefficients
    getc <- function(row, col) if (row %in% rownames(co)) co[row, col] else NA_real_
    depth_model_rows[[length(depth_model_rows) + 1L]] <- data.frame(
      capture_area_id = area_id, role_family = combo_role_family, program_id = pid,
      intercept = getc("(Intercept)", "Estimate"),
      beta_log1p_nCount = getc("log1p(nCount)", "Estimate"),
      beta_log1p_nFeature = getc("log1p(nFeature)", "Estimate"),
      se_log1p_nCount = getc("log1p(nCount)", "Std. Error"),
      se_log1p_nFeature = getc("log1p(nFeature)", "Std. Error"),
      t_log1p_nCount = getc("log1p(nCount)", "t value"),
      t_log1p_nFeature = getc("log1p(nFeature)", "t value"),
      p_log1p_nCount = getc("log1p(nCount)", "Pr(>|t|)"),
      p_log1p_nFeature = getc("log1p(nFeature)", "Pr(>|t|)"),
      r_squared = summary(fit)$r.squared, adjusted_r_squared = summary(fit)$adj.r.squared,
      residual_sd = sigma(fit),
      raw_cor_nCount = cor(dat$score, log1p(dat$nCount)),
      raw_cor_nFeature = cor(dat$score, log1p(dat$nFeature)),
      residual_cor_nCount = cor(res, log1p(dat$nCount)),
      residual_cor_nFeature = cor(res, log1p(dat$nFeature)),
      depth_covariate_correlation = cor(log1p(dat$nCount), log1p(dat$nFeature)),
      design_matrix_rank = rank, condition_number = cond,
      residual_complete = all(is.finite(res)), model_status = status,
      stringsAsFactors = FALSE)
    auth <- observed_authority[observed_authority$capture_area_id == area_id & observed_authority$program_id == pid, , drop = FALSE]
    if (status == "VALID") {
      raw_m <- auth$Moran_I; raw_g <- auth$Geary_C
      res_m <- wp3_moran_edge(res, area$edge_index)
      res_g <- wp3_geary_edge(res, area$edge_index)
      sm <- seeds[seeds$capture_area_id == area_id & seeds$program_order == program_order & seeds$analysis_type == "DEPTH_RESIDUAL_MORAN_PERMUTATION_ENDPOINT", "derived_seed"]
      sg <- seeds[seeds$capture_area_id == area_id & seeds$program_order == program_order & seeds$analysis_type == "DEPTH_RESIDUAL_GEARY_PERMUTATION_ENDPOINT", "derived_seed"]
      pm <- wp3_permutation_test_edge(res, area$edge_index, "moran", 9999L, sm)
      pg <- wp3_permutation_test_edge(res, area$edge_index, "geary", 9999L, sg)
      depth_spatial_rows[[length(depth_spatial_rows) + 1L]] <- data.frame(
        capture_area_id = area_id, role_family = combo_role_family, program_id = pid,
        raw_Moran_I = raw_m, residual_Moran_I = res_m, residual_Moran_empirical_p = pm$permutation_p,
        Moran_absolute_change = res_m - raw_m, Moran_relative_change = if (raw_m != 0) (res_m - raw_m) / abs(raw_m) else NA_real_,
        Moran_direction_concordance = sign(raw_m) == sign(res_m),
        Moran_effect_retention_ratio = if (raw_m > 0) res_m / raw_m else NA_real_,
        raw_Geary_C = raw_g, residual_Geary_C = res_g, residual_Geary_empirical_p = pg$permutation_p,
        Geary_absolute_change = res_g - raw_g,
        raw_Geary_spatial_departure_strength = 1 - raw_g,
        residual_Geary_spatial_departure_strength = 1 - res_g,
        Geary_retention = if ((1 - raw_g) > 0) (1 - res_g) / (1 - raw_g) else NA_real_,
        raw_residual_score_correlation = cor(dat$score, res),
        evaluation_status = "VALID", stringsAsFactors = FALSE)
      depth_perm_rows[[length(depth_perm_rows) + 1L]] <- data.frame(
        capture_area_id = area_id, role_family = combo_role_family, program_id = pid,
        endpoint = c("DEPTH_ADJUSTED_MORAN", "DEPTH_ADJUSTED_GEARY"),
        seed = c(sm, sg), permutations = 9999L,
        empirical_p = c(pm$permutation_p, pg$permutation_p),
        status = "COMPLETE", stringsAsFactors = FALSE)
    }
  }
}
depth_models <- do.call(rbind, depth_model_rows)
depth_spatial <- do.call(rbind, depth_spatial_rows)
depth_perm <- do.call(rbind, depth_perm_rows)
write_csv_once(depth_models, file.path(OUT, "depth/DEPTH_MODEL_DIAGNOSTICS.csv"))
write_csv_once(depth_spatial, file.path(OUT, "depth/DEPTH_RESIDUAL_SPATIAL_STATISTICS.csv"))
write_csv_once(depth_perm, file.path(OUT, "depth/DEPTH_RESIDUAL_PERMUTATION_SUMMARY.csv"))

fdr_rows <- list()
add_fdr_rows <- function(endpoint, p_table, p_col) {
  for (role in unique(fdr_contract$role_family)) {
    sub <- p_table[p_table$role_family == role, , drop = FALSE]
    fam <- fdr_contract[fdr_contract$endpoint == endpoint & fdr_contract$role_family == role, , drop = FALSE]
    if (!nrow(fam)) next
    fdr_vals <- bh(sub[[p_col]])
    fdr_rows[[length(fdr_rows) + 1L]] <<- data.frame(
      endpoint_type = endpoint, role_family = role, family_label = fam$family_label,
      expected_test_count = fam$expected_test_count, observed_valid_test_count = sum(is.finite(sub[[p_col]])),
      capture_area_id = sub$capture_area_id, program_id = sub$program_id,
      raw_p = sub[[p_col]], BH_FDR = fdr_vals, evaluability = ifelse(is.finite(sub[[p_col]]), "EVALUABLE", "NOT_EVALUABLE"),
      exclusion_failure_reason = "", stringsAsFactors = FALSE)
  }
}
if (nrow(empirical)) {
  add_fdr_rows("MATCHED_GENE_MORAN", empirical, "Moran_empirical_p")
  add_fdr_rows("MATCHED_GENE_GEARY", empirical, "Geary_empirical_p")
}
if (nrow(depth_spatial)) {
  add_fdr_rows("DEPTH_ADJUSTED_MORAN", depth_spatial, "residual_Moran_empirical_p")
  add_fdr_rows("DEPTH_ADJUSTED_GEARY", depth_spatial, "residual_Geary_empirical_p")
}
fdr_out <- if (length(fdr_rows)) do.call(rbind, fdr_rows) else data.frame()
write_csv_once(fdr_out, file.path(OUT, "fdr/06P_SEPARATE_FDR_RESULTS.csv"))
write_csv_once(if (length(failure_rows)) do.call(rbind, failure_rows) else data.frame(capture_area_id=character(), stage=character(), status=character(), reason=character()), file.path(OUT, "failures/06P_FAILURE_LOG.csv"))
write_csv_once(if (length(batch_rows)) do.call(rbind, batch_rows) else data.frame(), file.path(RUN, "06P_BATCH_MANIFEST.csv"))

spot_depth_all <- do.call(rbind, spot_depth_rows)
gene_summary_all <- do.call(rbind, gene_rows)
membership_path <- file.path(OUT, "matched_null/MATCHED_RANDOM_SET_MEMBERSHIP.csv")
membership <- read.csv(membership_path, check.names = FALSE)
valid_membership <- membership[membership$replicate_validity == "VALID", , drop = FALSE]
valid_keys <- paste(valid_membership$capture_area_id, valid_membership$program_id, valid_membership$replicate_id)
duplicate_by_valid_set <- tapply(valid_membership$selected_matched_gene, valid_keys, function(x) anyDuplicated(x) > 0)
scientific_csv_paths <- c(
  file.path(OUT, "gene_summaries/AREA_GENE_MATCHING_SUMMARIES.csv"),
  file.path(OUT, "gene_summaries/SPOT_DEPTH_COVARIATES.csv"),
  file.path(OUT, "matched_null/MATCHED_RANDOM_SET_MEMBERSHIP.csv"),
  file.path(OUT, "matched_null/MATCHING_QC.csv"),
  file.path(OUT, "matched_null/MATCHED_NULL_MORAN_GEARY.csv"),
  file.path(OUT, "matched_null/MATCHED_NULL_EMPIRICAL_TESTS.csv"),
  file.path(OUT, "depth/DEPTH_MODEL_DIAGNOSTICS.csv"),
  file.path(OUT, "depth/DEPTH_RESIDUAL_SPATIAL_STATISTICS.csv"),
  file.path(OUT, "depth/DEPTH_RESIDUAL_PERMUTATION_SUMMARY.csv"),
  file.path(OUT, "fdr/06P_SEPARATE_FDR_RESULTS.csv"),
  file.path(OUT, "failures/06P_FAILURE_LOG.csv"))
csv_parse_ok <- vapply(scientific_csv_paths, function(p) {
  tryCatch({ read.csv(p, nrows = 5L, check.names = FALSE); TRUE }, error = function(e) FALSE)
}, logical(1))
fdr_family_counts <- aggregate(raw_p ~ endpoint_type + role_family, fdr_out, length)
names(fdr_family_counts)[3] <- "observed_count"
fdr_expected <- merge(fdr_contract[, c("endpoint", "role_family", "expected_test_count")],
                      fdr_family_counts, by.x = c("endpoint", "role_family"),
                      by.y = c("endpoint_type", "role_family"), all.x = TRUE)
fdr_expected$observed_count[is.na(fdr_expected$observed_count)] <- 0L
git_tracked_changes <- system2("git", c("-C", shQuote(REV), "diff", "--name-only"), stdout = TRUE)
git_staged_changes <- system2("git", c("-C", shQuote(REV), "diff", "--cached", "--name-only"), stdout = TRUE)
raw_component_checks <- read.csv(file.path(ATT4, "06P_RAW_COUNT_COMPONENT_AUTHORITY_CHECKS_ATTEMPT004.csv"), check.names = FALSE)
area_dimension_checks <- read.csv(file.path(ATT4, "06P_RAW_COUNT_AREA_SUMMARY_ATTEMPT004.csv"), check.names = FALSE)
attempt4_validation <- read.csv(file.path(ATT4, "06P_PREFLIGHT_VALIDATION_ATTEMPT004.csv"), check.names = FALSE)
amd3_validation <- read.csv(file.path(PROTO, "amendments/AMENDMENT_003_06P_PREEXISTING_CONTROL_OUTPUT_STATE_CLARIFICATION/06P_AMENDMENT_003_VALIDATION.csv"), check.names = FALSE)

validation <- data.frame(
  check = c(
    "reconciliation_status_recorded", "attempt004_erratum_exists", "attempt004_preflight_pass",
    "amendment003_validation_pass", "token_valid", "human_authorization_valid", "o001_preserved",
    "nine_areas_complete", "raw_component_authority_81_of_81", "area_dimensions_9_of_9",
    "exact_retained_barcodes", "score_authority_checks_pass", "spatial_statistic_authority_checks_pass",
    "complete_spot_depth_covariates", "complete_area_gene_summaries", "matched_attempts_54000",
    "exact_frozen_matched_seed_count", "no_canonical_gene_in_valid_matched_sets",
    "no_duplicate_gene_within_valid_set", "valid_sets_threshold", "matching_failures_retained",
    "all_54_depth_models_attempted", "valid_108_depth_endpoints_completed",
    "permutations_9999_per_valid_endpoint", "exactly_12_fdr_families", "fdr_expected_counts",
    "all_scientific_csvs_parse", "failure_log_present", "protected_paths_unchanged_by_git",
    "tracked_git_changes_zero", "staged_git_changes_zero"),
  status = c(
    TRUE,
    file.exists(file.path(ATT4, "06P_PREFLIGHT_REPORT_ATTEMPT004_ERRATUM.md")),
    all(attempt4_validation$status == "PASS") && nrow(attempt4_validation) == 24L,
    all(amd3_validation$status == "PASS") && nrow(amd3_validation) == 4L,
    identical(readLines(file.path(P06, "00_execution_authorization/06P_EXECUTION_TOKEN.txt"), warn = FALSE), TOKEN),
    any(grepl("Human authorization: PROVIDED", readLines(file.path(P06, "00_execution_authorization/06P_EXECUTION_AUTHORIZATION.md"), warn = FALSE), fixed = TRUE)),
    sha_file(file.path(RUN, "06P_INPUT_SHA256_REGISTRY.csv")) == "f3f557ceb7e6baf3c0f3f2e076e22523d338e72d2522e08cab6a38266ea32147",
    length(areas) == 9L,
    nrow(raw_component_checks) == 81L && all(raw_component_checks$sha256_match == "TRUE"),
    nrow(area_dimension_checks) == 9L && all(area_dimension_checks$area_dimension_status == "PASS"),
    all(table(spot_depth_all$capture_area_id) == as.integer(inputs$spot_count)) && !anyDuplicated(paste(spot_depth_all$capture_area_id, spot_depth_all$barcode)),
    all(score_checks$status == "PASS"),
    all(spatial_checks$status == "PASS"),
    nrow(spot_depth_all) == sum(as.integer(inputs$spot_count)) && all(is.finite(spot_depth_all$nCount_Spatial)) &&
      all(spot_depth_all$nCount_Spatial >= 0) && all(is.finite(spot_depth_all$nFeature_Spatial)) &&
      all(spot_depth_all$nFeature_Spatial > 0),
    nrow(gene_summary_all) == sum(vapply(areas, function(a) nrow(a$gene_summary), integer(1))),
    sum(matching_qc$attempted_sets) == 54000L,
    nrow(seeds[seeds$analysis_type == "MATCHED_GENE_SET_REPLICATE", , drop = FALSE]) == 54000L,
    !any(valid_membership$selected_matched_gene %in% canonical_121),
    !any(unlist(duplicate_by_valid_set)),
    all(matching_qc$valid_sets >= 950L),
    nrow(membership) >= sum(matching_qc$attempted_sets),
    nrow(depth_models) == 54L,
    nrow(depth_perm) == 108L && all(depth_perm$status == "COMPLETE"),
    nrow(depth_perm) == 108L && all(depth_perm$permutations == 9999L),
    length(unique(paste(fdr_out$endpoint_type, fdr_out$role_family))) == 12L,
    nrow(fdr_expected) == 12L && all(fdr_expected$observed_count == as.integer(fdr_expected$expected_test_count)),
    all(csv_parse_ok),
    file.exists(file.path(OUT, "failures/06P_FAILURE_LOG.csv")),
    length(git_tracked_changes) == 0L && length(git_staged_changes) == 0L,
    length(git_tracked_changes) == 0L,
    length(git_staged_changes) == 0L),
  detail = c(
    "ATTEMPT004_PREFLIGHT_PASS_REPORT_TEXT_ERROR",
    file.path(ATT4, "06P_PREFLIGHT_REPORT_ATTEMPT004_ERRATUM.md"),
    paste(sum(attempt4_validation$status == "PASS"), "of", nrow(attempt4_validation)),
    paste(sum(amd3_validation$status == "PASS"), "of", nrow(amd3_validation)),
    TOKEN, "Human authorization: PROVIDED",
    "O001 SHA f3f557ceb7e6baf3c0f3f2e076e22523d338e72d2522e08cab6a38266ea32147",
    length(areas), "81/81", "9/9", nrow(spot_depth_all),
    paste(sum(score_checks$status == "PASS"), "of", nrow(score_checks)),
    paste(sum(spatial_checks$status == "PASS"), "of", nrow(spatial_checks)),
    nrow(spot_depth_all), nrow(gene_summary_all), sum(matching_qc$attempted_sets),
    54000L, "valid matched-set rows checked", "valid matched-set replicate keys checked",
    paste(min(matching_qc$valid_sets), "minimum valid sets"),
    paste(nrow(membership), "membership rows"),
    nrow(depth_models), nrow(depth_perm), paste(unique(depth_perm$permutations), collapse = ";"),
    length(unique(paste(fdr_out$endpoint_type, fdr_out$role_family))),
    paste(paste(fdr_expected$endpoint, fdr_expected$role_family, fdr_expected$observed_count, sep = ":"), collapse = ";"),
    paste(sum(csv_parse_ok), "of", length(csv_parse_ok)), file.path(OUT, "failures/06P_FAILURE_LOG.csv"),
    "Git tracked/staged diffs are zero", length(git_tracked_changes), length(git_staged_changes)),
  stringsAsFactors = FALSE)
validation$status_text <- ifelse(validation$status, "PASS", "FAIL")
authority_decision <- if (all(validation$status)) "FINAL_06P_AUTHORITY" else "COMPLETED_NOT_VALIDATED"
write_csv_once(validation, file.path(OUT, "validation/06P_EXECUTION_VALIDATION.csv"))
write_text_once(c("# 06p Authority Decision", "", paste("Decision:", authority_decision)), file.path(OUT, "validation/06P_AUTHORITY_DECISION.md"))

report <- c(
  "# 06p Final Report",
  "",
  "## 1. Reconciliation and Continuation Status",
  paste("Reconciliation status: ATTEMPT004_PREFLIGHT_PASS_REPORT_TEXT_ERROR."),
  paste("Execution and validation status:", authority_decision),
  "Attempt 004 erratum exists and the original report remains unchanged.",
  "",
  "## 2. Runtime and Package Versions",
  paste("Runtime R:", paste(R.version$major, R.version$minor, sep = ".")),
  paste("Package capture:", file.path(RUN, "06P_RUNTIME_PACKAGE_CAPTURE.txt")),
  "",
  "## 3. Input and Authority Gates",
  paste("Raw component authority:", sum(raw_component_checks$sha256_match == "TRUE"), "of", nrow(raw_component_checks)),
  paste("Area dimensions:", sum(area_dimension_checks$area_dimension_status == "PASS"), "of", nrow(area_dimension_checks)),
  paste("Score authority checks:", sum(score_checks$status == "PASS"), "of", nrow(score_checks)),
  paste("Spatial statistic authority checks:", sum(spatial_checks$status == "PASS"), "of", nrow(spatial_checks)),
  "",
  "## 4. Total Retained Spots",
  paste("Retained spots:", nrow(spot_depth_all)),
  "",
  "## 5. Gene-Universe Sizes by Area",
  paste(capture.output(print(do.call(rbind, gene_area_summary_rows)[, c("capture_area_id", "total_features", "detected_features", "eligible_random_universe_genes", "excluded_canonical_genes")], row.names = FALSE)), collapse = "\n"),
  "",
  "## 6. Matched-Set Validity by Program-Area",
  paste(capture.output(print(matching_qc[, c("capture_area_id", "role_family", "program_id", "attempted_sets", "valid_sets", "failed_sets", "null_validity_decision")], row.names = FALSE)), collapse = "\n"),
  "",
  "## 7. Matching Quality and Fallback Use",
  paste("Total valid matched sets:", sum(matching_qc$valid_sets)),
  paste("Total failed matched sets:", sum(matching_qc$failed_sets)),
  paste("Fallback-to-100 replicate count:", sum(matching_qc$fallback_100)),
  paste("Fallback-to-200 replicate count:", sum(matching_qc$fallback_200)),
  "",
  "## 8. Observed Versus Matched-Null Moran Results",
  paste(capture.output(print(empirical[, c("capture_area_id", "role_family", "program_id", "observed_Moran_I", "null_Moran_mean", "Moran_empirical_p", "Moran_standardized_effect")], row.names = FALSE)), collapse = "\n"),
  "",
  "## 9. Observed Versus Matched-Null Geary Results",
  paste(capture.output(print(empirical[, c("capture_area_id", "role_family", "program_id", "observed_Geary_C", "null_Geary_mean", "Geary_empirical_p", "Geary_standardized_effect")], row.names = FALSE)), collapse = "\n"),
  "",
  "## 10. Empirical P and BH-FDR Results",
  paste(capture.output(print(fdr_out[, c("endpoint_type", "role_family", "family_label", "observed_valid_test_count")], row.names = FALSE)), collapse = "\n"),
  "",
  "## 11. Depth-Model Diagnostics",
  paste("Depth models attempted:", nrow(depth_models)),
  paste("Valid depth models:", sum(depth_models$model_status == "VALID")),
  paste("Rank-deficient depth models:", sum(depth_models$model_status == "RANK_DEFICIENT")),
  "",
  "## 12. Raw Versus Residual Moran Results",
  paste(capture.output(print(depth_spatial[, c("capture_area_id", "role_family", "program_id", "raw_Moran_I", "residual_Moran_I", "residual_Moran_empirical_p", "Moran_effect_retention_ratio")], row.names = FALSE)), collapse = "\n"),
  "",
  "## 13. Raw Versus Residual Geary Results",
  paste(capture.output(print(depth_spatial[, c("capture_area_id", "role_family", "program_id", "raw_Geary_C", "residual_Geary_C", "residual_Geary_empirical_p", "Geary_retention")], row.names = FALSE)), collapse = "\n"),
  "",
  "## 14. Effect-Retention Estimates",
  "Effect-retention estimates are reported in DEPTH_RESIDUAL_SPATIAL_STATISTICS.csv with NOT_APPLICABLE represented by blank CSV values where denominators are invalid.",
  "",
  "## 15. Not-Evaluable Combinations",
  paste("Failure-log rows:", nrow(read.csv(file.path(OUT, "failures/06P_FAILURE_LOG.csv"), check.names = FALSE))),
  "",
  "## 16. Negative and Heterogeneous Findings",
  "All valid program-area results are retained regardless of direction, significance, or heterogeneity.",
  "",
  "## 17. Combined Signal Classes",
  "Program-area combinations should be interpreted from the matched-null and depth-adjusted FDR tables jointly: exceeded matched-gene null expectation, persisted after depth adjustment, satisfied both, or satisfied neither.",
  "",
  "## 18. Limitations",
  "Interpretation is limited to gene-set specificity relative to expression/detection-matched random genes and sensitivity to spot library size and detected-feature count.",
  "",
  "## 19. Eligibility for Manuscript, Supplementary Table, and Reviewer-Response Use",
  paste("Use eligibility follows final authority status:", authority_decision),
  "",
  "## 20. Final Authority Status",
  paste("Final authority status:", authority_decision),
  "",
  "No causal spatial organization, cell-cell communication, patient-level replication, independent normal-control validation, spatial cell-state reconstruction, DLBCL specificity of context-area results, or complete technical-confounding-removal claim is made."
)
write_text_once(report, file.path(OUT, "06P_FINAL_REPORT.md"))

manifest_paths <- c(
  file.path(RUN, "06P_INPUT_SHA256_REGISTRY.csv"),
  file.path(RUN, "06P_RUNTIME_PACKAGE_CAPTURE.txt"),
  file.path(RUN, "06P_SEED_RUNTIME_VERIFICATION.csv"),
  file.path(RUN, "06P_BATCH_MANIFEST.csv"),
  scientific_csv_paths,
  file.path(OUT, "validation/06P_EXECUTION_VALIDATION.csv"),
  file.path(OUT, "validation/06P_OUTPUT_MANIFEST.csv"),
  file.path(OUT, "validation/06P_AUTHORITY_DECISION.md"),
  file.path(OUT, "06P_FINAL_REPORT.md"))
manifest_exists <- file.exists(manifest_paths)
manifest_exists[basename(manifest_paths) == "06P_OUTPUT_MANIFEST.csv"] <- TRUE
manifest <- data.frame(
  path = manifest_paths,
  exists = manifest_exists,
  size_bytes = ifelse(manifest_exists & basename(manifest_paths) != "06P_OUTPUT_MANIFEST.csv",
                      as.numeric(file.info(manifest_paths)$size), NA_real_),
  sha256 = ifelse(manifest_exists & basename(manifest_paths) != "06P_OUTPUT_MANIFEST.csv",
                  vapply(manifest_paths, function(p) if (file.exists(p) && basename(p) != "06P_OUTPUT_MANIFEST.csv") sha_file(p) else NA_character_, character(1)),
                  NA_character_),
  creation_origin = ifelse(basename(manifest_paths) == "06P_INPUT_SHA256_REGISTRY.csv",
                           "ATTEMPT_001_PREFLIGHT", "ATTEMPT_004_SCIENTIFIC_EXECUTION"),
  final_execution_action = ifelse(basename(manifest_paths) == "06P_INPUT_SHA256_REGISTRY.csv",
                                  "PRESERVED_AND_REUSED", "CREATED_ONCE"),
  overwrite_performed = "FALSE",
  effective_input_authority = ifelse(basename(manifest_paths) == "06P_INPUT_SHA256_REGISTRY.csv",
                                     "81_OF_81_COMPONENT_SHA256_PASS", ""),
  legacy_composite_hash_status = ifelse(basename(manifest_paths) == "06P_INPUT_SHA256_REGISTRY.csv",
                                        "NONCOMPARABLE_IMPLEMENTATIONS_RETAINED_FOR_AUDIT", ""),
  stringsAsFactors = FALSE)
write_csv_once(manifest, file.path(OUT, "validation/06P_OUTPUT_MANIFEST.csv"))

message("06p scientific execution complete: ", authority_decision)
