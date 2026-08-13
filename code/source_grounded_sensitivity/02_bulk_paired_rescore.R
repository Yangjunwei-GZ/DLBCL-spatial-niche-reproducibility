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
o6_require_token(); o6_set_library(); o6_require_packages(c("GSVA", "BiocParallel", "Matrix", "digest"))
set.seed(O6_MASTER_SEED)

sets <- o6_gene_sets()
datasets <- list(
  GSE31312 = list(path = file.path(O6_PROJECT, "00_raw_data/GSE31312/GSE31312_series_matrix.txt.gz"),
    annotation = file.path(O6_PROJECT, "04_tables/GSE10846/GSE10846_GPL570_probe_to_symbol_mapping.csv"),
    expected_n = 498L, collapse = "highest_mean", platform = "GPL570"),
  GSE10846 = list(path = file.path(O6_PROJECT, "00_raw_data/GSE10846/GSE10846_series_matrix.txt.gz"),
    annotation = file.path(O6_PROJECT, "04_tables/GSE10846/GSE10846_GPL570_probe_to_symbol_mapping.csv"),
    expected_n = 420L, collapse = "highest_mean", platform = "GPL570"),
  GSE181063 = list(path = file.path(O6_PROJECT, "00_raw_data/GSE181063/GSE181063_series_matrix.txt.gz"),
    annotation = file.path(O6_PROJECT, "results/GSE181063_extension/GSE181063_illuminaHumanv4_probe_to_symbol_mapping.csv"),
    expected_n = 1310L, collapse = "highest_mad", platform = "IlluminaHumanv4")
)

matching_rows <- list(); score_rows <- list(); summary_rows <- list(); geometry_rows <- list()
new_full <- list()

for (dataset_id in names(datasets)) {
  cfg <- datasets[[dataset_id]]
  message("Reading and preprocessing ", dataset_id)
  expr_probe <- o6_parse_series_matrix(cfg$path, cfg$expected_n)
  annotation <- o6_read_annotation(cfg$annotation, rownames(expr_probe))
  collapsed <- o6_collapse_highest(expr_probe, annotation, cfg$collapse)
  expr <- collapsed$expression
  if (ncol(expr) != cfg$expected_n || anyDuplicated(colnames(expr))) stop("Sample identity failure: ", dataset_id)

  full_detected <- lapply(sets$full, intersect, y = rownames(expr))
  core_detected <- lapply(sets$core, intersect, y = rownames(expr))
  if (any(lengths(full_detected) < 5L) || any(lengths(core_detected) < 5L)) {
    stop("Frozen minSize=5 gene coverage failure: ", dataset_id)
  }
  full_score <- o6_ssgsea(expr, sets$full)
  core_score <- o6_ssgsea(expr, sets$core)
  if (!identical(colnames(full_score), colnames(core_score)) ||
      !identical(colnames(full_score), colnames(expr))) stop("Paired sample identity failure: ", dataset_id)
  new_full[[dataset_id]] <- full_score

  for (set_type in c("FULL", "CORE")) {
    gs <- if (set_type == "FULL") sets$full else sets$core
    detected <- if (set_type == "FULL") full_detected else core_detected
    for (program_id in names(gs)) {
      missing <- setdiff(gs[[program_id]], detected[[program_id]])
      matching_rows[[length(matching_rows) + 1L]] <- data.frame(
        dataset_id = dataset_id, platform = cfg$platform,
        input_matrix = normalizePath(cfg$path, winslash = "/"),
        annotation_path = normalizePath(cfg$annotation, winslash = "/"),
        annotation_sha256 = o6_sha256(cfg$annotation), preprocessing_rule = toupper(cfg$collapse),
        program = o6_program_name(program_id), program_id = program_id, set_type = set_type,
        canonical_gene_count = length(gs[[program_id]]), matched_gene_count = length(detected[[program_id]]),
        missing_gene_count = length(missing), matched_genes = paste(detected[[program_id]], collapse = ";"),
        missing_genes = paste(missing, collapse = ";"), mapping_status = "EVALUABLE",
        exclusion_reason = "", sample_count = ncol(expr),
        sample_identifiers = paste(colnames(expr), collapse = ";"), stringsAsFactors = FALSE)
    }
  }

  for (set_type in c("FULL", "CORE")) {
    score <- if (set_type == "FULL") full_score else core_score
    score_rows[[length(score_rows) + 1L]] <- data.frame(dataset_id = dataset_id,
      sample_id = rep(colnames(score), each = nrow(score)), set_type = set_type,
      program = rep(vapply(rownames(score), o6_program_name, character(1)), times = ncol(score)),
      program_id = rep(rownames(score), times = ncol(score)), score = as.vector(score),
      GSVA_version = as.character(utils::packageVersion("GSVA")),
      ssgsea_minSize = 5L, ssgsea_maxSize = 500L, ssgsea_alpha = 0.25,
      ssgsea_normalize = TRUE, parallel_backend = "BiocParallel::SerialParam", stringsAsFactors = FALSE)
  }

  common_programs <- names(sets$core)
  for (program_id in common_programs) {
    f <- full_score[program_id, ]; c <- core_score[program_id, ]
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(dataset_id = dataset_id,
      program = o6_program_name(program_id), program_id = program_id, sample_count = length(f),
      full_gene_count = length(sets$full[[program_id]]), core_gene_count = length(sets$core[[program_id]]),
      full_genes_present = length(full_detected[[program_id]]), core_genes_present = length(core_detected[[program_id]]),
      full_mean = mean(f), full_sd = stats::sd(f), core_mean = mean(c), core_sd = stats::sd(c),
      full_core_Pearson = o6_safe_cor(f, c, "pearson"), full_core_Spearman = o6_safe_cor(f, c, "spearman"),
      evaluation_status = "EVALUABLE", stringsAsFactors = FALSE)
  }

  pairs <- utils::combn(common_programs, 2L, simplify = FALSE)
  for (method in c("pearson", "spearman")) {
    full_cor <- stats::cor(t(full_score[common_programs, , drop = FALSE]), method = method)
    core_cor <- stats::cor(t(core_score[common_programs, , drop = FALSE]), method = method)
    full_values <- core_values <- numeric(length(pairs))
    for (j in seq_along(pairs)) {
      p <- pairs[[j]]; full_values[j] <- full_cor[p[1L], p[2L]]; core_values[j] <- core_cor[p[1L], p[2L]]
      geometry_rows[[length(geometry_rows) + 1L]] <- data.frame(dataset_id = dataset_id,
        correlation_method = method, row_type = "PAIRWISE_ASSOCIATION", program_1 = o6_program_name(p[1L]),
        program_1_id = p[1L], program_2 = o6_program_name(p[2L]), program_2_id = p[2L],
        full_correlation = full_values[j], core_correlation = core_values[j],
        core_minus_full = core_values[j] - full_values[j],
        direction_concordant = sign(full_values[j]) == sign(core_values[j]),
        rank_concordance_ten_pairs = NA_real_, direction_concordance_fraction = NA_real_, stringsAsFactors = FALSE)
    }
    geometry_rows[[length(geometry_rows) + 1L]] <- data.frame(dataset_id = dataset_id,
      correlation_method = method, row_type = "TEN_PAIR_SUMMARY", program_1 = "ALL_FIVE_PROGRAM_PAIRS",
      program_1_id = "ALL", program_2 = "ALL_FIVE_PROGRAM_PAIRS", program_2_id = "ALL",
      full_correlation = NA_real_, core_correlation = NA_real_, core_minus_full = NA_real_,
      direction_concordant = NA, rank_concordance_ten_pairs = o6_safe_cor(full_values, core_values, "spearman"),
      direction_concordance_fraction = mean(sign(full_values) == sign(core_values)), stringsAsFactors = FALSE)
  }
  rm(expr_probe, expr, full_score, core_score, collapsed); gc(verbose = FALSE)
}

o6_write_csv_once(do.call(rbind, matching_rows), file.path(O6_OUTPUTS, "bulk/PAIRED_FULL_AND_CORE_GENE_MATCHING.csv"))
o6_write_csv_once(do.call(rbind, score_rows), file.path(O6_OUTPUTS, "bulk/PAIRED_FULL_AND_CORE_SCORES.csv"))
o6_write_csv_once(do.call(rbind, summary_rows), file.path(O6_OUTPUTS, "bulk/FULL_CORE_SCORE_SUMMARY.csv"))
o6_write_csv_once(do.call(rbind, geometry_rows), file.path(O6_OUTPUTS, "bulk/FULL_CORE_PROGRAM_CORRELATION_COMPARISON.csv"))

map_historical_labels <- function(x) {
  out <- rep(NA_character_, length(x)); low <- tolower(x)
  out[grepl("macrophage", low)] <- "macrophage_rich"
  out[grepl("t cell|t-cell", low)] <- "t_cell_inflamed"
  out[grepl("antigen|immune-inflamed", low)] <- "antigen_presentation"
  out[grepl("stromal|fibrotic", low)] <- "stromal_fibrotic"
  out[grepl("cold|exclusion", low)] <- "immune_cold_exclusion"
  out[grepl("proliferative|cycling", low)] <- "proliferative_cycling"
  out
}

historical_path <- file.path(O6_PROJECT, "02_processed_data/GSE31312_spatial_niche_ssGSEA_scores.csv")
historical <- utils::read.csv(historical_path, check.names = FALSE)
if (nrow(historical) != 6L) stop("Historical raw authority orientation is not 6 x samples")
historical_ids <- map_historical_labels(historical[[1L]])
if (anyNA(historical_ids) || anyDuplicated(historical_ids)) stop("Historical raw program labels cannot be mapped")
historical_mat <- as.matrix(historical[, -1L, drop = FALSE]); storage.mode(historical_mat) <- "double"
rownames(historical_mat) <- historical_ids

primary_path <- file.path(O6_REVISION,
  "05e_stage4_GSE31312_execution_attempt2/01_score_space_validation/GSE31312_primary_score_matrix_498x6.csv")
primary <- utils::read.csv(primary_path, check.names = FALSE)
primary_samples <- primary[[1L]]
primary_ids <- map_historical_labels(names(primary)[-1L])
primary_mat <- t(as.matrix(primary[, -1L, drop = FALSE])); storage.mode(primary_mat) <- "double"
rownames(primary_mat) <- primary_ids; colnames(primary_mat) <- primary_samples

comparison_rows <- list(); new <- new_full$GSE31312
for (authority_name in c("HISTORICAL_RAW_6x498", "FINAL_PRIMARY_Z_498x6")) {
  target <- if (authority_name == "HISTORICAL_RAW_6x498") historical_mat else primary_mat
  transformed_new <- if (authority_name == "HISTORICAL_RAW_6x498") new else t(scale(t(new)))
  samples_agree <- setequal(colnames(transformed_new), colnames(target))
  programs_agree <- setequal(rownames(transformed_new), rownames(target))
  for (program_id in unname(O6_PROGRAM_IDS)) {
    common <- intersect(colnames(transformed_new), colnames(target))
    x <- transformed_new[program_id, common]; y <- target[program_id, common]
    fit <- stats::lm(y ~ x)
    difference <- x - y
    comparison_rows[[length(comparison_rows) + 1L]] <- data.frame(
      authority = authority_name, authority_path = normalizePath(if (authority_name == "HISTORICAL_RAW_6x498") historical_path else primary_path, winslash = "/"),
      authority_sha256 = o6_sha256(if (authority_name == "HISTORICAL_RAW_6x498") historical_path else primary_path),
      new_score_scale = if (authority_name == "HISTORICAL_RAW_6x498") "RAW_SSGSEA" else "WITHIN_PROGRAM_Z_SCORE",
      authority_orientation = if (authority_name == "HISTORICAL_RAW_6x498") "PROGRAMS_BY_SAMPLES" else "SAMPLES_BY_PROGRAMS_TRANSPOSED_FOR_COMPARISON",
      sample_set_agreement = samples_agree, program_set_agreement = programs_agree,
      sample_order_agreement_before_matching = identical(colnames(transformed_new), colnames(target)),
      program = o6_program_name(program_id), program_id = program_id, matched_sample_count = length(common),
      Pearson = o6_safe_cor(x, y, "pearson"), Spearman = o6_safe_cor(x, y, "spearman"),
      mean_absolute_difference = mean(abs(difference)), root_mean_squared_difference = sqrt(mean(difference^2)),
      maximum_absolute_difference = max(abs(difference)), scaling_intercept_authority_on_new = coef(fit)[[1L]],
      scaling_slope_authority_on_new = coef(fit)[[2L]], comparison_status = "DESCRIPTIVE_NOT_FORCED_EQUAL",
      historical_GSVA_version = "UNKNOWN_MAY_2026", new_GSVA_version = as.character(utils::packageVersion("GSVA")),
      stringsAsFactors = FALSE)
  }
}
o6_write_csv_once(do.call(rbind, comparison_rows), file.path(O6_OUTPUTS, "bulk/FULL_RESCORE_VS_HISTORICAL_AUTHORITY.csv"))
cat("BULK_EXECUTION=PASS\n")
