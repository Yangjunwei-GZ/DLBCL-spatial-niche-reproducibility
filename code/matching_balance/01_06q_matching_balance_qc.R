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
Q06 <- file.path(REV, "06q_wp3_matched_null_balance_qc")
PROTO <- file.path(Q06, "00_protocol")
OUT <- file.path(Q06, "02_execution_outputs")
LOG_DIR <- file.path(Q06, "03_execution_logs")
LOG <- file.path(LOG_DIR, "01_06q_matching_balance_qc.log")

if (!requireNamespace("digest", quietly = TRUE)) stop("Missing required package: digest")

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
if (file.exists(LOG)) stop("create-once refusal: log exists")

log_msg <- function(...) {
  line <- paste(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), paste(..., collapse = " "), sep = "\t")
  cat(line, "\n", file = LOG, append = TRUE)
  message(line)
}
sha_file <- function(path) digest::digest(path, file = TRUE, algo = "sha256")
read_csv <- function(path) read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
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
truthy <- function(x) toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1", "PASS", "VALID", "COMPLETE")
finite_num <- function(x) is.finite(as.numeric(x))
qnum <- function(x, p) as.numeric(quantile(x, probs = p, na.rm = TRUE, names = FALSE, type = 7))
safe_cor <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  as.numeric(cor(x[ok], y[ok], method = method))
}
pooled_smd <- function(a, b) {
  psd <- sqrt((stats::sd(a)^2 + stats::sd(b)^2) / 2)
  if (!is.finite(psd) || psd == 0) return(NA_real_)
  (mean(b) - mean(a)) / psd
}
add_validation <- function(rows, check_id, check, ok, detail) {
  rows[[length(rows) + 1L]] <- data.frame(
    check_id = check_id,
    check = check,
    status = if (isTRUE(ok)) "PASS" else "FAIL",
    detail = paste(detail, collapse = ""),
    stringsAsFactors = FALSE
  )
  rows
}

outputs <- file.path(OUT, c(
  "06Q_REPLICATE_LEVEL_MATCHING_BALANCE.csv",
  "06Q_AREA_PROGRAM_MATCHING_BALANCE_SUMMARY.csv",
  "06Q_PROGRAM_LEVEL_MATCHING_BALANCE_SUMMARY.csv",
  "06Q_DLBCL_MATCHING_BALANCE_SUMMARY.csv",
  "06Q_BALANCE_DISTANCE_CORRELATIONS.csv",
  "06Q_ANTIGEN_PRESENTATION_BALANCE_AUDIT.csv",
  "06Q_EXECUTION_VALIDATION.csv",
  "06Q_MATCHING_BALANCE_QC_REPORT.md"
))
if (any(file.exists(outputs))) stop("create-once refusal: output exists: ", paste(outputs[file.exists(outputs)], collapse = "; "))

log_msg("06Q_START")

hash_registry <- read_csv(file.path(PROTO, "06Q_INPUT_HASH_REGISTRY.csv"))
hash_registry$current_sha256 <- vapply(hash_registry$path, sha_file, character(1))
hash_registry$current_size_bytes <- as.numeric(file.info(hash_registry$path)$size)
hash_registry$hash_status <- ifelse(hash_registry$current_sha256 == hash_registry$sha256 &
                                      hash_registry$current_size_bytes == hash_registry$size_bytes, "PASS", "FAIL")
if (!all(hash_registry$hash_status == "PASS")) stop("Registered input hash check failed")

p06_final_report <- paste(readLines(file.path(P06, "01_execution_outputs", "06P_FINAL_REPORT.md"), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
p06_validation <- read_csv(file.path(P06, "01_execution_outputs", "validation", "06P_EXECUTION_VALIDATION.csv"))
membership_validation <- read_csv(file.path(P06, "01_execution_outputs", "matched_null", "MATCHED_RANDOM_SET_MEMBERSHIP_VALIDATION.csv"))
partition_manifest <- read_csv(file.path(P06, "01_execution_outputs", "matched_null", "MATCHED_RANDOM_SET_MEMBERSHIP_PARTITION_MANIFEST.csv"))
batch_manifest <- read_csv(file.path(P06, "01_execution_outputs", "run_control", "06P_BATCH_MANIFEST.csv"))
matching_qc <- read_csv(file.path(P06, "01_execution_outputs", "matched_null", "MATCHING_QC.csv"))
gene_summary <- read_csv(file.path(P06, "01_execution_outputs", "gene_summaries", "AREA_GENE_MATCHING_SUMMARIES.csv"))
input_registry <- read_csv(file.path(P06, "00_protocol_freeze", "06P_INPUT_REGISTRY.csv"))
seed_registry <- read_csv(file.path(P06, "00_protocol_freeze", "06P_SEED_REGISTRY.csv"))

partition_manifest$current_sha256 <- vapply(partition_manifest$output_file, sha_file, character(1))
partition_manifest$current_size_bytes <- as.numeric(file.info(partition_manifest$output_file)$size)
partition_hash_pass <- partition_manifest$current_sha256 == partition_manifest$output_sha256 &
  partition_manifest$current_size_bytes == partition_manifest$output_size_bytes

program_map <- unique(partition_manifest[, c("program_id", "program_name")])
program_map <- program_map[order(program_map$program_id), , drop = FALSE]
program_map$program_order <- match(program_map$program_id, c("macrophage_rich", "t_cell_inflamed", "antigen_presentation", "stromal_fibrotic", "immune_cold_exclusion", "proliferative_cycling"))
program_map <- program_map[order(program_map$program_order), , drop = FALSE]

seed_matched <- seed_registry[seed_registry$analysis_type == "MATCHED_GENE_SET_REPLICATE", , drop = FALSE]
seed_matched$program_order <- as.integer(seed_matched$program_order)
seed_matched$replicate_id <- as.integer(seed_matched$replicate_or_endpoint)
seed_matched$program_id <- program_map$program_id[match(seed_matched$program_order, program_map$program_order)]
seed_matched$key <- paste(seed_matched$capture_area_id, seed_matched$program_id, seed_matched$replicate_id, sep = "||")
seed_expected <- setNames(as.integer(seed_matched$derived_seed), seed_matched$key)

pre_rows <- list()
pre_rows <- add_validation(pre_rows, "06Q-PRE001", "06p FINAL_06P_AUTHORITY present",
                           grepl("FINAL_06P_AUTHORITY", p06_final_report, fixed = TRUE) && all(p06_validation$status == "PASS"),
                           "06p final authority report and validation")
pre_rows <- add_validation(pre_rows, "06Q-PRE002", "Membership validation PASS",
                           all(membership_validation$status == "PASS"), paste(sum(membership_validation$status == "PASS"), "of", nrow(membership_validation)))
pre_rows <- add_validation(pre_rows, "06Q-PRE003", "487 partitions present unchanged",
                           nrow(partition_manifest) == 487L && all(partition_hash_pass) && all(partition_manifest$validation_status == "PASS"),
                           paste(nrow(partition_manifest), "partitions"))
pre_rows <- add_validation(pre_rows, "06Q-PRE004", "6 salvaged and 481 generated partitions",
                           sum(partition_manifest$partition_origin == "SALVAGED_FROM_INTERRUPTED_FILE") == 6L &&
                             sum(partition_manifest$partition_origin == "GENERATED_BY_CONTROLLED_RESUME") == 481L,
                           paste(sum(partition_manifest$partition_origin == "SALVAGED_FROM_INTERRUPTED_FILE"), "salvaged;", sum(partition_manifest$partition_origin == "GENERATED_BY_CONTROLLED_RESUME"), "generated"))
pre_rows <- add_validation(pre_rows, "06Q-PRE005", "AREA_GENE_MATCHING_SUMMARIES rows",
                           nrow(gene_summary) == 162315L, nrow(gene_summary))
pre_rows <- add_validation(pre_rows, "06Q-PRE006", "MATCHING_QC complete",
                           nrow(matching_qc) == 54L && all(matching_qc$attempted_sets == 1000L) && all(matching_qc$valid_sets == 1000L),
                           paste(nrow(matching_qc), "rows;", sum(matching_qc$valid_sets), "valid sets"))
pre_validation <- do.call(rbind, pre_rows)
if (any(pre_validation$status != "PASS")) stop("Pre-execution authority check failed: ", paste(pre_validation$check_id[pre_validation$status != "PASS"], collapse = "; "))

log_msg("READ_PARTITIONS_START", nrow(partition_manifest))
membership_parts <- lapply(partition_manifest$output_file, read_csv)
membership <- do.call(rbind, membership_parts)
rm(membership_parts)
membership$replicate_id <- as.integer(membership$replicate_id)
membership$seed <- as.integer(membership$seed)
membership$key <- paste(membership$capture_area_id, membership$program_id, membership$replicate_id, sep = "||")

key_table <- unique(membership[, c("capture_area_id", "role_family", "program_id", "program_name", "replicate_id", "seed", "key")])
key_counts <- table(paste(key_table$capture_area_id, key_table$program_id, sep = "||"))
missing_reps <- any(unlist(tapply(key_table$replicate_id, paste(key_table$capture_area_id, key_table$program_id, sep = "||"), function(x) !identical(sort(x), 1:1000))))
if (nrow(key_table) != 54000L || anyDuplicated(key_table$key) || length(key_counts) != 54L || any(key_counts != 1000L) || missing_reps) {
  stop("Composite key authority failed after partition read")
}
if (!all(key_table$seed == seed_expected[key_table$key])) stop("Frozen seed identity failed")
if (any(membership$capture_area_id == "GSM8500534_Cap.area1_LN_V1" &
        membership$program_id == "antigen_presentation" &
        membership$observed_canonical_gene == "HLA-DRB1")) {
  stop("Area1 antigen-presentation HLA-DRB1 exclusion failed")
}
log_msg("READ_PARTITIONS_COMPLETE", nrow(membership), "membership rows", nrow(key_table), "keys")

expr_var <- "log1p_mean_LogNormalized_expression"
det_var <- "spot_detection_fraction"
gene_summary[[expr_var]] <- as.numeric(gene_summary[[expr_var]])
gene_summary[[det_var]] <- as.numeric(gene_summary[[det_var]])
gene_summary$eligible_random_universe <- truthy(gene_summary$eligible_random_universe)
gene_summary$canonical_121_membership <- truthy(gene_summary$canonical_121_membership)

area_ref_rows <- list()
for (area_id in unique(gene_summary$capture_area_id)) {
  area_genes <- gene_summary[gene_summary$capture_area_id == area_id, , drop = FALSE]
  universe <- area_genes[area_genes$eligible_random_universe & !area_genes$canonical_121_membership &
                           finite_num(area_genes[[expr_var]]) & finite_num(area_genes[[det_var]]), , drop = FALSE]
  area_ref_rows[[length(area_ref_rows) + 1L]] <- data.frame(
    capture_area_id = area_id,
    area_expression_reference_mean = mean(universe[[expr_var]]),
    area_expression_reference_sd = stats::sd(universe[[expr_var]]),
    area_detection_reference_mean = mean(universe[[det_var]]),
    area_detection_reference_sd = stats::sd(universe[[det_var]]),
    eligible_random_universe_size = nrow(universe),
    stringsAsFactors = FALSE
  )
}
area_refs <- do.call(rbind, area_ref_rows)
if (any(!is.finite(area_refs$area_expression_reference_sd) | area_refs$area_expression_reference_sd <= 0 |
        !is.finite(area_refs$area_detection_reference_sd) | area_refs$area_detection_reference_sd <= 0)) {
  stop("Area reference SD check failed")
}

gene_lookup <- gene_summary[, c("capture_area_id", "gene_symbol", expr_var, det_var), drop = FALSE]
names(gene_lookup) <- c("capture_area_id", "gene_symbol", "expression_value", "detection_value")
gene_lookup$key <- paste(gene_lookup$capture_area_id, gene_lookup$gene_symbol, sep = "||")
expr_by_key <- setNames(gene_lookup$expression_value, gene_lookup$key)
det_by_key <- setNames(gene_lookup$detection_value, gene_lookup$key)
ref_by_area <- split(area_refs, area_refs$capture_area_id)
area_role <- setNames(input_registry$area_role, input_registry$capture_area_id)

log_msg("REPLICATE_BALANCE_START")
membership$observed_key <- paste(membership$capture_area_id, membership$observed_canonical_gene, sep = "||")
membership$matched_key <- paste(membership$capture_area_id, membership$selected_matched_gene, sep = "||")
membership$observed_expression <- expr_by_key[membership$observed_key]
membership$matched_expression <- expr_by_key[membership$matched_key]
membership$observed_detection <- det_by_key[membership$observed_key]
membership$matched_detection <- det_by_key[membership$matched_key]

split_membership <- split(membership, membership$key)
rep_rows <- vector("list", length(split_membership))
idx <- 0L
for (key in names(split_membership)) {
  z <- split_membership[[key]]
  idx <- idx + 1L
  ref <- ref_by_area[[z$capture_area_id[1]]]
  obs_expr <- unique(z[, c("observed_canonical_gene", "observed_expression")])
  obs_det <- unique(z[, c("observed_canonical_gene", "observed_detection")])
  mat_expr <- z$matched_expression
  mat_det <- z$matched_detection
  invalid <- any(!is.finite(obs_expr$observed_expression)) || any(!is.finite(obs_det$observed_detection)) ||
    any(!is.finite(mat_expr)) || any(!is.finite(mat_det))
  exp_smd <- (mean(mat_expr) - mean(obs_expr$observed_expression)) / ref$area_expression_reference_sd
  det_smd <- (mean(mat_det) - mean(obs_det$observed_detection)) / ref$area_detection_reference_sd
  exp_pool <- pooled_smd(obs_expr$observed_expression, mat_expr)
  det_pool <- pooled_smd(obs_det$observed_detection, mat_det)
  pool_zero <- (!is.finite(exp_pool)) || (!is.finite(det_pool))
  status <- if (invalid || !is.finite(exp_smd) || !is.finite(det_smd)) "INVALID_INPUT" else if (pool_zero) "POOLED_SD_ZERO" else "COMPLETE"
  rep_rows[[idx]] <- data.frame(
    capture_area_id = z$capture_area_id[1],
    role_family = z$role_family[1],
    program_id = z$program_id[1],
    program_name = z$program_name[1],
    replicate_id = z$replicate_id[1],
    seed = z$seed[1],
    observed_set_size = nrow(obs_expr),
    matched_set_size = length(mat_expr),
    observed_expression_mean = mean(obs_expr$observed_expression),
    matched_expression_mean = mean(mat_expr),
    observed_detection_mean = mean(obs_det$observed_detection),
    matched_detection_mean = mean(mat_det),
    area_expression_reference_mean = ref$area_expression_reference_mean,
    area_expression_reference_sd = ref$area_expression_reference_sd,
    area_detection_reference_mean = ref$area_detection_reference_mean,
    area_detection_reference_sd = ref$area_detection_reference_sd,
    expression_balance_smd = exp_smd,
    abs_expression_balance_smd = abs(exp_smd),
    detection_balance_smd = det_smd,
    abs_detection_balance_smd = abs(det_smd),
    expression_pooled_smd = exp_pool,
    abs_expression_pooled_smd = abs(exp_pool),
    detection_pooled_smd = det_pool,
    abs_detection_pooled_smd = abs(det_pool),
    mean_matching_distance = mean(z$matching_distance, na.rm = TRUE),
    max_matching_distance = max(z$matching_distance, na.rm = TRUE),
    pool_50_count = sum(z$pool_size_used == 50, na.rm = TRUE),
    pool_100_count = sum(z$pool_size_used == 100, na.rm = TRUE),
    pool_200_count = sum(z$pool_size_used == 200, na.rm = TRUE),
    replicate_balance_status = status,
    notes = if (status == "COMPLETE") "" else status,
    stringsAsFactors = FALSE
  )
}
rep_balance <- do.call(rbind, rep_rows)
rep_balance <- rep_balance[order(rep_balance$capture_area_id, rep_balance$program_id, rep_balance$replicate_id), , drop = FALSE]
log_msg("REPLICATE_BALANCE_COMPLETE", nrow(rep_balance))

combo_split <- split(rep_balance, paste(rep_balance$capture_area_id, rep_balance$program_id, sep = "||"))
combo_rows <- lapply(combo_split, function(z) {
  data.frame(
    capture_area_id = z$capture_area_id[1],
    role_family = z$role_family[1],
    program_id = z$program_id[1],
    program_name = z$program_name[1],
    observed_set_size = z$observed_set_size[1],
    replicate_count = nrow(z),
    mean_expression_smd = mean(z$expression_balance_smd),
    median_expression_smd = median(z$expression_balance_smd),
    sd_expression_smd = stats::sd(z$expression_balance_smd),
    mean_abs_expression_smd = mean(z$abs_expression_balance_smd),
    median_abs_expression_smd = median(z$abs_expression_balance_smd),
    p90_abs_expression_smd = qnum(z$abs_expression_balance_smd, 0.90),
    p95_abs_expression_smd = qnum(z$abs_expression_balance_smd, 0.95),
    max_abs_expression_smd = max(z$abs_expression_balance_smd),
    prop_abs_expression_smd_le_0_05 = mean(z$abs_expression_balance_smd <= 0.05),
    prop_abs_expression_smd_le_0_10 = mean(z$abs_expression_balance_smd <= 0.10),
    prop_abs_expression_smd_le_0_20 = mean(z$abs_expression_balance_smd <= 0.20),
    mean_detection_smd = mean(z$detection_balance_smd),
    median_detection_smd = median(z$detection_balance_smd),
    sd_detection_smd = stats::sd(z$detection_balance_smd),
    mean_abs_detection_smd = mean(z$abs_detection_balance_smd),
    median_abs_detection_smd = median(z$abs_detection_balance_smd),
    p90_abs_detection_smd = qnum(z$abs_detection_balance_smd, 0.90),
    p95_abs_detection_smd = qnum(z$abs_detection_balance_smd, 0.95),
    max_abs_detection_smd = max(z$abs_detection_balance_smd),
    prop_abs_detection_smd_le_0_05 = mean(z$abs_detection_balance_smd <= 0.05),
    prop_abs_detection_smd_le_0_10 = mean(z$abs_detection_balance_smd <= 0.10),
    prop_abs_detection_smd_le_0_20 = mean(z$abs_detection_balance_smd <= 0.20),
    median_abs_expression_pooled_smd = median(z$abs_expression_pooled_smd, na.rm = TRUE),
    p95_abs_expression_pooled_smd = qnum(z$abs_expression_pooled_smd, 0.95),
    max_abs_expression_pooled_smd = max(z$abs_expression_pooled_smd, na.rm = TRUE),
    median_abs_detection_pooled_smd = median(z$abs_detection_pooled_smd, na.rm = TRUE),
    p95_abs_detection_pooled_smd = qnum(z$abs_detection_pooled_smd, 0.95),
    max_abs_detection_pooled_smd = max(z$abs_detection_pooled_smd, na.rm = TRUE),
    mean_matching_distance = mean(z$mean_matching_distance),
    median_matching_distance = median(z$mean_matching_distance),
    max_matching_distance = max(z$max_matching_distance),
    pool_50_total = sum(z$pool_50_count),
    pool_100_total = sum(z$pool_100_count),
    pool_200_total = sum(z$pool_200_count),
    complete_replicates = sum(z$replicate_balance_status == "COMPLETE"),
    noncomplete_replicates = sum(z$replicate_balance_status != "COMPLETE"),
    stringsAsFactors = FALSE
  )
})
combo_summary <- do.call(rbind, combo_rows)
combo_summary <- combo_summary[order(combo_summary$capture_area_id, combo_summary$program_id), , drop = FALSE]

program_split <- split(rep_balance, rep_balance$program_id)
program_rows <- lapply(program_split, function(z) {
  data.frame(
    program_id = z$program_id[1],
    program_name = z$program_name[1],
    combinations = length(unique(paste(z$capture_area_id, z$program_id, sep = "||"))),
    replicates = nrow(z),
    median_abs_expression_smd = median(z$abs_expression_balance_smd),
    p95_abs_expression_smd = qnum(z$abs_expression_balance_smd, 0.95),
    max_abs_expression_smd = max(z$abs_expression_balance_smd),
    median_abs_detection_smd = median(z$abs_detection_balance_smd),
    p95_abs_detection_smd = qnum(z$abs_detection_balance_smd, 0.95),
    max_abs_detection_smd = max(z$abs_detection_balance_smd),
    prop_both_abs_smd_le_0_05 = mean(z$abs_expression_balance_smd <= 0.05 & z$abs_detection_balance_smd <= 0.05),
    prop_both_abs_smd_le_0_10 = mean(z$abs_expression_balance_smd <= 0.10 & z$abs_detection_balance_smd <= 0.10),
    prop_both_abs_smd_le_0_20 = mean(z$abs_expression_balance_smd <= 0.20 & z$abs_detection_balance_smd <= 0.20),
    median_matching_distance = median(z$mean_matching_distance),
    max_matching_distance = max(z$max_matching_distance),
    stringsAsFactors = FALSE
  )
})
program_summary <- do.call(rbind, program_rows)
program_summary <- program_summary[match(program_map$program_id, program_summary$program_id), , drop = FALSE]

dlbcl_areas <- names(area_role)[area_role == "PRIMARY_DLBCL"]
dlbcl <- rep_balance[rep_balance$capture_area_id %in% dlbcl_areas, , drop = FALSE]
dlbcl_summary <- data.frame(
  area_role = "PRIMARY_DLBCL",
  areas = length(unique(dlbcl$capture_area_id)),
  combinations = length(unique(paste(dlbcl$capture_area_id, dlbcl$program_id, sep = "||"))),
  replicates = nrow(dlbcl),
  median_abs_expression_smd = median(dlbcl$abs_expression_balance_smd),
  p95_abs_expression_smd = qnum(dlbcl$abs_expression_balance_smd, 0.95),
  max_abs_expression_smd = max(dlbcl$abs_expression_balance_smd),
  median_abs_detection_smd = median(dlbcl$abs_detection_balance_smd),
  p95_abs_detection_smd = qnum(dlbcl$abs_detection_balance_smd, 0.95),
  max_abs_detection_smd = max(dlbcl$abs_detection_balance_smd),
  prop_both_abs_smd_le_0_05 = mean(dlbcl$abs_expression_balance_smd <= 0.05 & dlbcl$abs_detection_balance_smd <= 0.05),
  prop_both_abs_smd_le_0_10 = mean(dlbcl$abs_expression_balance_smd <= 0.10 & dlbcl$abs_detection_balance_smd <= 0.10),
  prop_both_abs_smd_le_0_20 = mean(dlbcl$abs_expression_balance_smd <= 0.20 & dlbcl$abs_detection_balance_smd <= 0.20),
  median_matching_distance = median(dlbcl$mean_matching_distance),
  max_matching_distance = max(dlbcl$max_matching_distance),
  stringsAsFactors = FALSE
)

cor_pairs <- list(
  median_abs_expression_smd = combo_summary$median_abs_expression_smd,
  median_abs_detection_smd = combo_summary$median_abs_detection_smd,
  p95_abs_expression_smd = combo_summary$p95_abs_expression_smd,
  p95_abs_detection_smd = combo_summary$p95_abs_detection_smd
)
cor_rows <- list()
for (metric in names(cor_pairs)) {
  for (method in c("pearson", "spearman")) {
    cor_rows[[length(cor_rows) + 1L]] <- data.frame(
      balance_metric = metric,
      distance_metric = "mean_matching_distance",
      correlation_method = toupper(method),
      n_combinations = nrow(combo_summary),
      correlation = safe_cor(cor_pairs[[metric]], combo_summary$mean_matching_distance, method),
      interpretation_scope = "DESCRIPTIVE_QC_ONLY",
      stringsAsFactors = FALSE
    )
  }
}
distance_cor <- do.call(rbind, cor_rows)

antigen <- combo_summary[combo_summary$program_id == "antigen_presentation", , drop = FALSE]
antigen$hla_drb1_inclusion_status <- ifelse(antigen$capture_area_id == "GSM8500534_Cap.area1_LN_V1", "EXCLUDED", "NOT_APPLICABLE_AREA")
antigen_audit <- antigen[, c(
  "capture_area_id", "role_family", "program_id", "program_name", "observed_set_size",
  "median_abs_expression_smd", "p95_abs_expression_smd", "max_abs_expression_smd",
  "median_abs_detection_smd", "p95_abs_detection_smd", "max_abs_detection_smd",
  "mean_matching_distance", "max_matching_distance"
), drop = FALSE]
antigen_keys <- paste(antigen_audit$capture_area_id, antigen_audit$program_id, sep = "||")
antigen_audit$prop_both_abs_smd_le_0_10 <- vapply(antigen_keys, function(k) {
  z <- rep_balance[paste(rep_balance$capture_area_id, rep_balance$program_id, sep = "||") == k, , drop = FALSE]
  mean(z$abs_expression_balance_smd <= 0.10 & z$abs_detection_balance_smd <= 0.10)
}, numeric(1))
antigen_audit$prop_both_abs_smd_le_0_20 <- vapply(antigen_keys, function(k) {
  z <- rep_balance[paste(rep_balance$capture_area_id, rep_balance$program_id, sep = "||") == k, , drop = FALSE]
  mean(z$abs_expression_balance_smd <= 0.20 & z$abs_detection_balance_smd <= 0.20)
}, numeric(1))
antigen_audit$hla_drb1_inclusion_status <- ifelse(antigen_audit$capture_area_id == "GSM8500534_Cap.area1_LN_V1", "EXCLUDED", "NOT_INCLUDED_IN_AREA_PROGRAM")

write_csv_once(rep_balance, outputs[1])
write_csv_once(combo_summary, outputs[2])
write_csv_once(program_summary, outputs[3])
write_csv_once(dlbcl_summary, outputs[4])
write_csv_once(distance_cor, outputs[5])
write_csv_once(antigen_audit, outputs[6])

post_hash_registry <- hash_registry
post_hash_registry$post_sha256 <- vapply(post_hash_registry$path, sha_file, character(1))
post_hash_registry$post_size_bytes <- as.numeric(file.info(post_hash_registry$path)$size)
partition_manifest$post_sha256 <- vapply(partition_manifest$output_file, sha_file, character(1))
partition_unchanged_after <- partition_manifest$post_sha256 == partition_manifest$output_sha256

tracked_all <- system2("git", c("-C", shQuote(REV), "diff", "--name-only"), stdout = TRUE)
staged_all <- system2("git", c("-C", shQuote(REV), "diff", "--cached", "--name-only"), stdout = TRUE)
tracked_06p <- system2("git", c("-C", shQuote(REV), "diff", "--name-only", "--", "06p_wp3_matched_null_and_depth_sensitivity"), stdout = TRUE)
tracked_06o <- system2("git", c("-C", shQuote(REV), "diff", "--name-only", "--", "06o_source_grounded_program_sensitivity"), stdout = TRUE)
protected_paths <- c("05t_stage4c2a_v4_real_continuation_attempt", "06k_wp3_figure7_manuscript_integration", "source_snapshot/manuscript", "source_snapshot/workbooks")
tracked_protected <- unlist(lapply(protected_paths, function(p) system2("git", c("-C", shQuote(REV), "diff", "--name-only", "--", p), stdout = TRUE)))

validation_rows <- list()
validation_rows <- add_validation(validation_rows, "06Q-V001", "06p FINAL_06P_AUTHORITY present",
                                  grepl("FINAL_06P_AUTHORITY", p06_final_report, fixed = TRUE) && all(p06_validation$status == "PASS"),
                                  "FINAL_06P_AUTHORITY")
validation_rows <- add_validation(validation_rows, "06Q-V002", "06p input hashes unchanged",
                                  all(post_hash_registry$post_sha256 == post_hash_registry$sha256) &&
                                    all(post_hash_registry$post_size_bytes == post_hash_registry$size_bytes),
                                  paste(sum(post_hash_registry$post_sha256 == post_hash_registry$sha256), "of", nrow(post_hash_registry)))
validation_rows <- add_validation(validation_rows, "06Q-V003", "partition manifest unchanged",
                                  identical(hash_registry$current_sha256[hash_registry$input_id == "06Q-IN004"], hash_registry$sha256[hash_registry$input_id == "06Q-IN004"]) &&
                                    identical(hash_registry$current_sha256[hash_registry$input_id == "06Q-IN006"], hash_registry$sha256[hash_registry$input_id == "06Q-IN006"]),
                                  "partition manifest and batch manifest hashes match registry")
validation_rows <- add_validation(validation_rows, "06Q-V004", "487 partitions unchanged",
                                  nrow(partition_manifest) == 487L && all(partition_unchanged_after),
                                  paste(nrow(partition_manifest), "partitions"))
validation_rows <- add_validation(validation_rows, "06Q-V005", "54000 replicate keys",
                                  nrow(key_table) == 54000L && !anyDuplicated(key_table$key),
                                  nrow(key_table))
validation_rows <- add_validation(validation_rows, "06Q-V006", "54000 replicate-level QC rows",
                                  nrow(rep_balance) == 54000L, nrow(rep_balance))
validation_rows <- add_validation(validation_rows, "06Q-V007", "54 combination-summary rows",
                                  nrow(combo_summary) == 54L, nrow(combo_summary))
validation_rows <- add_validation(validation_rows, "06Q-V008", "6 program-summary rows",
                                  nrow(program_summary) == 6L, nrow(program_summary))
validation_rows <- add_validation(validation_rows, "06Q-V009", "all nine antigen combinations present",
                                  nrow(antigen_audit) == 9L, nrow(antigen_audit))
validation_rows <- add_validation(validation_rows, "06Q-V010", "exact frozen seeds preserved",
                                  all(key_table$seed == seed_expected[key_table$key]), "seed equality")
validation_rows <- add_validation(validation_rows, "06Q-V011", "expression SMD finite where area reference SD positive",
                                  all(is.finite(rep_balance$expression_balance_smd)), "all finite")
validation_rows <- add_validation(validation_rows, "06Q-V012", "detection SMD finite where area reference SD positive",
                                  all(is.finite(rep_balance$detection_balance_smd)), "all finite")
validation_rows <- add_validation(validation_rows, "06Q-V013", "no random-set membership changed",
                                  all(partition_unchanged_after), "partition hashes unchanged")
validation_rows <- add_validation(validation_rows, "06Q-V014", "no empirical P recomputed",
                                  TRUE, "06q does not read or write matched-null empirical P output")
validation_rows <- add_validation(validation_rows, "06Q-V015", "no FDR recomputed",
                                  TRUE, "06q does not read or write FDR output")
validation_rows <- add_validation(validation_rows, "06Q-V016", "no 06p file modified",
                                  all(post_hash_registry$post_sha256 == post_hash_registry$sha256) && length(tracked_06p) == 0L,
                                  paste("tracked_06p_diff", length(tracked_06p)))
validation_rows <- add_validation(validation_rows, "06Q-V017", "no 06o file modified",
                                  length(tracked_06o) == 0L, paste("tracked_06o_diff", length(tracked_06o)))
validation_rows <- add_validation(validation_rows, "06Q-V018", "no manuscript figure workbook response-letter file modified",
                                  length(tracked_protected) == 0L, paste("tracked_protected_diff", length(tracked_protected)))
validation_rows <- add_validation(validation_rows, "06Q-V019", "tracked Git changes reported",
                                  TRUE, paste("tracked_total", length(tracked_all)))
validation_rows <- add_validation(validation_rows, "06Q-V020", "staged Git changes reported",
                                  TRUE, paste("staged_total", length(staged_all)))
validation_rows <- add_validation(validation_rows, "06Q-V021", "no commit or push",
                                  TRUE, "script performs no commit or push")
validation <- do.call(rbind, validation_rows)
validation$authority_status <- ifelse(all(validation$status == "PASS"), "06Q_MATCHING_BALANCE_QC_PASS", "06Q_MATCHING_BALANCE_QC_FAIL")
write_csv_once(validation, outputs[7])

overall_median_expr <- median(rep_balance$abs_expression_balance_smd)
overall_p95_expr <- qnum(rep_balance$abs_expression_balance_smd, 0.95)
overall_median_det <- median(rep_balance$abs_detection_balance_smd)
overall_p95_det <- qnum(rep_balance$abs_detection_balance_smd, 0.95)
prop_both_005 <- mean(rep_balance$abs_expression_balance_smd <= 0.05 & rep_balance$abs_detection_balance_smd <= 0.05)
prop_both_010 <- mean(rep_balance$abs_expression_balance_smd <= 0.10 & rep_balance$abs_detection_balance_smd <= 0.10)
prop_both_020 <- mean(rep_balance$abs_expression_balance_smd <= 0.20 & rep_balance$abs_detection_balance_smd <= 0.20)
weak <- combo_summary[order(pmax(combo_summary$p95_abs_expression_smd, combo_summary$p95_abs_detection_smd), decreasing = TRUE), ][1:min(10L, nrow(combo_summary)), ]

report <- c(
  "# 06q Matched-Null Balance QC Report",
  "",
  "## 1. Purpose",
  "Balance was quantified post hoc using the same two variables and area-specific standardization used for matched-set construction.",
  "",
  "## 2. Input Authority",
  "Inputs are read-only from completed 06p authority. Membership is read only from validated partition files, not from the interrupted monolithic membership file.",
  "",
  "## 3. Exact SMD Definitions",
  "Primary expression and detection SMDs are matched-set mean z-score minus observed canonical-set mean z-score using the area-specific eligible random-gene universe as reference.",
  "Secondary pooled-SD SMDs use sqrt((sd_observed^2 + sd_matched^2) / 2); pooled-SD-zero cases are reported as NA and not replaced.",
  "",
  "## 4. Why Area-Reference SMD Is Primary",
  "The area-reference denominator maps the balance diagnostic directly to the frozen 06p matching space and does not re-standardize each observed/matched comparison separately.",
  "",
  "## 5. Replicate-Level Balance Results",
  paste("Replicate-level rows:", nrow(rep_balance)),
  paste("Complete replicate QC rows:", sum(rep_balance$replicate_balance_status == "COMPLETE")),
  paste("Noncomplete replicate QC rows:", sum(rep_balance$replicate_balance_status != "COMPLETE")),
  "",
  "## 6. Overall Expression-Balance Distribution",
  paste("Median absolute expression SMD:", signif(overall_median_expr, 6)),
  paste("P95 absolute expression SMD:", signif(overall_p95_expr, 6)),
  "",
  "## 7. Overall Detection-Balance Distribution",
  paste("Median absolute detection SMD:", signif(overall_median_det, 6)),
  paste("P95 absolute detection SMD:", signif(overall_p95_det, 6)),
  "",
  "## 8. Program-Level Differences",
  paste(capture.output(print(program_summary[, c("program_id", "median_abs_expression_smd", "p95_abs_expression_smd", "median_abs_detection_smd", "p95_abs_detection_smd")], row.names = FALSE)), collapse = "\n"),
  "",
  "## 9. DLBCL-Only Balance Results",
  paste(capture.output(print(dlbcl_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## 10. Antigen-Presentation Diagnostic",
  paste(capture.output(print(antigen_audit[, c("capture_area_id", "median_abs_expression_smd", "p95_abs_expression_smd", "median_abs_detection_smd", "p95_abs_detection_smd", "prop_both_abs_smd_le_0_10", "prop_both_abs_smd_le_0_20", "hla_drb1_inclusion_status")], row.names = FALSE)), collapse = "\n"),
  "",
  "## 11. Relationship With Matching Distance",
  paste(capture.output(print(distance_cor, row.names = FALSE)), collapse = "\n"),
  "",
  "## 12. Descriptive Thresholds",
  paste("Both primary |SMD| <= 0.05:", signif(prop_both_005, 6)),
  paste("Both primary |SMD| <= 0.10:", signif(prop_both_010, 6)),
  paste("Both primary |SMD| <= 0.20:", signif(prop_both_020, 6)),
  "",
  "## 13. Combinations With Weaker Balance",
  paste(capture.output(print(weak[, c("capture_area_id", "program_id", "p95_abs_expression_smd", "p95_abs_detection_smd", "max_abs_expression_smd", "max_abs_detection_smd")], row.names = FALSE)), collapse = "\n"),
  "",
  "## 14. Limitations",
  "These are descriptive balance diagnostics only. No post-hoc SMD cutoff is used to remove replicates or invalidate combinations.",
  "",
  "## 15. 06p Scientific Result Preservation",
  "No 06p scientific result, matched-set membership, replicate validity, empirical P value, FDR value, or FINAL_06P_AUTHORITY assignment was changed."
)
write_text_once(report, outputs[8])

if (any(validation$status != "PASS")) stop("06q validation failed: ", paste(validation$check_id[validation$status != "PASS"], collapse = "; "))
log_msg("06Q_COMPLETE", "06Q_MATCHING_BALANCE_QC_PASS")
