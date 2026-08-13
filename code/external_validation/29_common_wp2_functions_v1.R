DLBCL_PROJECT_ROOT <- normalizePath(
  Sys.getenv("DLBCL_PROJECT_ROOT", unset = getwd()),
  winslash = "/", mustWork = FALSE
)
DLBCL_SUPPLEMENTARY_CODE_ROOT <- normalizePath(
  Sys.getenv("DLBCL_SUPPLEMENTARY_CODE_ROOT",
    unset = file.path(DLBCL_PROJECT_ROOT, "Supplementary_Code")),
  winslash = "/", mustWork = FALSE
)

# Shared functions for WP2 real external continuous validation.
# Scientific decisions are frozen in Amendment 023 and must not be changed here.

options(stringsAsFactors = FALSE)

WP2 <- list(
  root = DLBCL_PROJECT_ROOT,
  stage_rel = "revision_2026_reviewer_response/05z_wp2_real_external_continuous_validation",
  baseline = "0848dfab53b5e0ab36723f60f40ed03dff2bd226",
  branch = "revision/wp2-real-external-continuous-validation",
  seed = 20260730L,
  bootstrap_replicates = 2000L,
  token_name = "DLBCL_REVISION_ALLOW_WP2_EXTERNAL",
  token_value = "EXPLICITLY_APPROVED_WP2_EXTERNAL_CONTINUOUS_VALIDATION",
  programs = c(
    "Macrophage-rich program",
    "T cell-inflamed program",
    "Immune-inflamed / antigen-presentation program",
    "Stromal / fibrotic program",
    "Immune-cold / exclusion-associated program",
    "Proliferative / cycling program"
  )
)
WP2$stage <- file.path(WP2$root, WP2$stage_rel)
WP2$renv_library <- file.path(
  WP2$root,
  "revision_2026_reviewer_response/04b_stage4_environment_freeze/stage4_renv_project/renv/library/windows/R-4.5/x86_64-w64-mingw32"
)
WP2$paths <- list(
  GSE10846_expression = file.path(WP2$root, "00_raw_data/GSE10846/GSE10846_series_matrix.txt.gz"),
  GSE10846_metadata = file.path(WP2$root, "02_processed_data/GSE10846/GSE10846_phenotype_raw.csv"),
  GSE10846_annotation = file.path(WP2$root, "04_tables/GSE10846/GSE10846_GPL570_probe_to_symbol_mapping.csv"),
  GSE181063_expression = file.path(WP2$root, "00_raw_data/GSE181063/GSE181063_series_matrix.txt.gz"),
  GSE181063_metadata = file.path(WP2$root, "results/GSE181063_extension/GSE181063_sample_metadata_from_series_matrix.csv"),
  GSE181063_annotation = file.path(WP2$root, "results/GSE181063_extension/GSE181063_illuminaHumanv4_probe_to_symbol_mapping.csv"),
  canonical = file.path(WP2$root, "revision_2026_reviewer_response/05x_wp1_continuous_score_freeze/WP1_CANONICAL_PROGRAM_CONTRACT.csv"),
  projection_parameters = file.path(WP2$root, "revision_2026_reviewer_response/05x_wp1_continuous_score_freeze/WP1_CONTINUOUS_PROJECTION_PARAMETERS.csv"),
  pca_authority = file.path(WP2$root, "revision_2026_reviewer_response/05x_wp1_continuous_score_freeze/WP1_FROZEN_PCA_AUTHORITY.csv"),
  discovery_correlations = file.path(WP2$root, "revision_2026_reviewer_response/05x_wp1_continuous_score_freeze/WP1_PRIMARY_PROGRAM_PEARSON_CORRELATIONS_LONG.csv"),
  discovery_program_pc = file.path(WP2$root, "revision_2026_reviewer_response/05x_wp1_continuous_score_freeze/WP1_PROGRAM_PC_PEARSON_CORRELATIONS.csv"),
  discovery_raw_scores = file.path(WP2$root, "revision_2026_reviewer_response/05e_stage4_GSE31312_execution_attempt2/01_score_space_validation/GSE31312_historical_raw_score_matrix_6x498.csv"),
  scoring_authority = file.path(WP2$root, "revision_2026_reviewer_response/02_canonical_manifest_and_pipeline_rebuild/scripts/03_score_GSE31312_programs.R"),
  protected_manifest = file.path(WP2$root, "revision_2026_reviewer_response/05x_wp1_continuous_score_freeze/WP1_PROTECTED_PATH_BASELINE_MANIFEST.csv")
)

wp2_require_packages <- function() {
  required <- c("data.table", "digest", "matrixStats", "GSVA", "BiocParallel")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing frozen packages: ", paste(missing, collapse = ";"), call. = FALSE)
  versions <- vapply(required, function(x) as.character(utils::packageVersion(x)), character(1))
  if (versions[["GSVA"]] != "2.4.9") stop("Frozen GSVA 2.4.9 required", call. = FALSE)
  invisible(versions)
}

wp2_sha <- function(path) digest::digest(file = path, algo = "sha256")
wp2_norm <- function(path) gsub("\\\\", "/", normalizePath(path, winslash = "/", mustWork = FALSE))
wp2_rel <- function(path) {
  x <- wp2_norm(path)
  prefix <- paste0(wp2_norm(WP2$root), "/")
  ifelse(startsWith(x, prefix), substring(x, nchar(prefix) + 1L), x)
}
wp2_timestamp <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

wp2_assert_token <- function() {
  if (!identical(Sys.getenv(WP2$token_name, unset = ""), WP2$token_value)) {
    stop("WP2 real-execution token is absent or invalid", call. = FALSE)
  }
  invisible(TRUE)
}

wp2_assert_r451 <- function() {
  if (!(R.version$major == "4" && startsWith(R.version$minor, "5.1"))) {
    stop("R 4.5.1 is required", call. = FALSE)
  }
}

wp2_create_once <- function(path) {
  if (file.exists(path)) stop("Create-once output already exists: ", path, call. = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

wp2_write_csv <- function(x, path) {
  wp2_create_once(path)
  utils::write.csv(x, path, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
  invisible(path)
}

wp2_write_text <- function(x, path) {
  wp2_create_once(path)
  writeLines(x, path, useBytes = TRUE)
  invisible(path)
}

wp2_save_rds <- function(x, path) {
  wp2_create_once(path)
  saveRDS(x, path, compress = "gzip")
  invisible(path)
}

wp2_program_contract <- function() {
  x <- utils::read.csv(WP2$paths$canonical, check.names = FALSE, stringsAsFactors = FALSE)
  if (!identical(as.integer(x$program_order), 1:6)) stop("Frozen program order mismatch", call. = FALSE)
  if (!identical(x$program_name, WP2$programs)) stop("Frozen program names mismatch", call. = FALSE)
  if (any(x$canonical_gene_count != 22L)) stop("Each program must contain 22 genes", call. = FALSE)
  x$genes <- strsplit(x$canonical_gene_members, ";", fixed = TRUE)
  if (any(lengths(x$genes) != 22L)) stop("Canonical memberships failed parsing", call. = FALSE)
  x
}

wp2_fread_gzip <- function(path, ...) {
  tmp <- tempfile(fileext = ".txt")
  input <- gzfile(path, open = "rb")
  output <- file(tmp, open = "wb")
  on.exit({
    try(close(input), silent = TRUE)
    try(close(output), silent = TRUE)
    unlink(tmp)
  }, add = TRUE)
  repeat {
    chunk <- readBin(input, what = "raw", n = 8L * 1024L * 1024L)
    if (!length(chunk)) break
    writeBin(chunk, output)
  }
  close(input)
  close(output)
  size <- as.numeric(file.info(tmp)$size)
  tail_start <- max(0, size - 1024L * 1024L)
  tail_con <- file(tmp, open = "rb")
  seek(tail_con, where = tail_start, origin = "start")
  tail_raw <- readBin(tail_con, what = "raw", n = size - tail_start)
  close(tail_con)
  marker <- grepRaw("!series_matrix_table_end", tail_raw, fixed = TRUE)[[1L]]
  if (!length(marker) || is.na(marker)) stop("Series matrix end marker not found", call. = FALSE)
  truncate_at <- tail_start + marker - 1L
  truncate_con <- file(tmp, open = "r+b")
  seek(truncate_con, where = truncate_at, origin = "start")
  truncate(truncate_con)
  close(truncate_con)
  data.table::fread(tmp, ...)
}

wp2_read_series_matrix <- function(dataset_id, expression_path, metadata_path, expected_n) {
  wp2_assert_token()
  meta <- utils::read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"sample" %in% names(meta)) stop("Metadata sample column missing", call. = FALSE)
  if (anyDuplicated(meta$sample)) stop("Duplicated metadata GSM", call. = FALSE)

  tab <- wp2_fread_gzip(
    expression_path,
    skip = "!series_matrix_table_begin",
    header = TRUE,
    data.table = FALSE,
    check.names = FALSE,
    showProgress = TRUE,
    quote = "\""
  )
  if (ncol(tab) < 2L) stop("Expression table did not parse", call. = FALSE)
  id_col <- names(tab)[[1L]]
  probe <- gsub('^"|"$', "", as.character(tab[[id_col]]))
  keep <- !startsWith(probe, "!") & nzchar(probe)
  tab <- tab[keep, , drop = FALSE]
  probe <- probe[keep]
  if (anyDuplicated(probe)) stop("Duplicated expression probe IDs for ", dataset_id, call. = FALSE)
  sample_ids <- gsub('^"|"$', "", names(tab)[-1L])
  if (length(sample_ids) != expected_n) stop("Unexpected expression sample count for ", dataset_id, call. = FALSE)
  if (anyDuplicated(sample_ids)) stop("Duplicated expression GSM", call. = FALSE)
  if (!setequal(sample_ids, meta$sample)) stop("Expression/metadata GSM sets differ", call. = FALSE)
  join_index <- match(sample_ids, meta$sample)
  if (anyNA(join_index)) stop("Explicit metadata join failed", call. = FALSE)
  meta <- meta[join_index, , drop = FALSE]

  expr <- as.matrix(tab[, -1L, drop = FALSE])
  suppressWarnings(storage.mode(expr) <- "double")
  rownames(expr) <- probe
  colnames(expr) <- sample_ids
  if (any(!is.finite(expr))) stop("Non-finite or nonnumeric expression value", call. = FALSE)
  if (any(colSums(!is.na(expr)) == 0L)) stop("All-NA expression sample", call. = FALSE)
  q <- stats::quantile(expr, c(0, .25, .5, .75, .99, 1), names = FALSE, na.rm = TRUE)
  scale_status <- if (q[[5L]] <= 100 && q[[6L]] <= 1000) "PLATFORM_NORMALIZED_LOG_LIKE_SCALE_PASS" else "UNEXPECTED_SCALE_STOP"
  if (scale_status != "PLATFORM_NORMALIZED_LOG_LIKE_SCALE_PASS") stop("Expression scale gate failed; no transformation applied", call. = FALSE)
  list(
    dataset_id = dataset_id,
    expression = expr,
    metadata = meta,
    audit = data.frame(
      dataset_id = dataset_id,
      probe_id_column = id_col,
      feature_count = nrow(expr),
      duplicated_probe_id_count = anyDuplicated(probe),
      expression_sample_count = ncol(expr),
      metadata_sample_count = nrow(meta),
      duplicated_expression_GSM = anyDuplicated(sample_ids),
      duplicated_metadata_GSM = anyDuplicated(meta$sample),
      explicit_join_complete = identical(sample_ids, meta$sample),
      numeric_parse_complete = TRUE,
      nonfinite_expression_count = sum(!is.finite(expr)),
      all_NA_sample_count = sum(colSums(!is.na(expr)) == 0L),
      minimum = q[[1L]], Q1 = q[[2L]], median = q[[3L]], Q3 = q[[4L]], Q99 = q[[5L]], maximum = q[[6L]],
      scale_status = scale_status,
      transformation_applied = "NONE",
      clinical_fields_used = FALSE,
      stringsAsFactors = FALSE
    )
  )
}

wp2_clean_annotation <- function(annotation_path, expression_probe_ids) {
  ann <- utils::read.csv(annotation_path, check.names = FALSE, stringsAsFactors = FALSE)
  ann <- ann[, c("PROBEID", "SYMBOL"), drop = FALSE]
  ann$PROBEID <- trimws(as.character(ann$PROBEID))
  ann$SYMBOL <- toupper(trimws(as.character(ann$SYMBOL)))
  invalid <- is.na(ann$PROBEID) | !nzchar(ann$PROBEID) | is.na(ann$SYMBOL) | !nzchar(ann$SYMBOL)
  ann_valid <- unique(ann[!invalid & ann$PROBEID %in% expression_probe_ids, , drop = FALSE])
  probe_symbol_n <- table(ann_valid$PROBEID)
  ambiguous_probes <- names(probe_symbol_n)[probe_symbol_n > 1L]
  excluded_rows <- list(
    blank = unique(ann$PROBEID[invalid]),
    ambiguous = ambiguous_probes,
    absent = setdiff(unique(ann$PROBEID[!invalid]), expression_probe_ids)
  )
  excluded_reasons <- c(
    blank = "BLANK_OR_NA_IDENTIFIER_OR_SYMBOL",
    ambiguous = "ONE_PROBE_MAPS_TO_MULTIPLE_GENES",
    absent = "PROBE_NOT_PRESENT_IN_EXPRESSION_AUTHORITY"
  )
  excluded <- do.call(rbind, lapply(names(excluded_rows), function(name) data.frame(
    PROBEID = excluded_rows[[name]],
    reason = rep(excluded_reasons[[name]], length(excluded_rows[[name]])),
    stringsAsFactors = FALSE
  )))
  if (is.null(excluded)) excluded <- data.frame(PROBEID = character(), reason = character(), stringsAsFactors = FALSE)
  ann_valid <- ann_valid[!ann_valid$PROBEID %in% ambiguous_probes, , drop = FALSE]
  list(annotation = ann_valid, excluded = unique(excluded[nzchar(excluded$PROBEID), , drop = FALSE]))
}

wp2_collapse_expression <- function(expr, annotation) {
  keep <- rownames(expr) %in% annotation$PROBEID
  expr_use <- expr[keep, , drop = FALSE]
  ann_use <- annotation[match(rownames(expr_use), annotation$PROBEID), , drop = FALSE]
  groups <- split(seq_len(nrow(expr_use)), ann_use$SYMBOL)
  genes <- sort(names(groups))
  primary <- matrix(NA_real_, nrow = length(genes), ncol = ncol(expr_use), dimnames = list(genes, colnames(expr_use)))
  sensitivity <- primary
  row_mad <- matrixStats::rowMads(expr_use, na.rm = TRUE)
  mapping_rows <- vector("list", length(genes))

  for (i in seq_along(genes)) {
    idx <- groups[[genes[[i]]]]
    probes <- rownames(expr_use)[idx]
    if (length(idx) == 1L) {
      primary[i, ] <- expr_use[idx, ]
    } else {
      primary[i, ] <- matrixStats::colMedians(expr_use[idx, , drop = FALSE], na.rm = TRUE)
    }
    max_mad <- max(row_mad[idx], na.rm = TRUE)
    tied <- probes[row_mad[idx] == max_mad]
    selected <- sort(tied)[[1L]]
    selected_idx <- idx[match(selected, probes)]
    sensitivity[i, ] <- expr_use[selected_idx, ]
    mapping_rows[[i]] <- data.frame(
      gene_symbol = genes[[i]],
      probe_count = length(idx),
      probes = paste(sort(probes), collapse = ";"),
      primary_rule = "MEDIAN_ACROSS_PROBES_PER_GENE",
      sensitivity_rule = "HIGHEST_MAD_PROBE_LEXICAL_TIEBREAK",
      selected_sensitivity_probe = selected,
      selected_probe_MAD = max_mad,
      MAD_tie_count = length(tied),
      outcome_used = FALSE,
      program_membership_used_for_selection = FALSE,
      stringsAsFactors = FALSE
    )
  }
  list(primary = primary, sensitivity = sensitivity, mapping = do.call(rbind, mapping_rows))
}

wp2_coverage <- function(dataset_id, collapse_rule, gene_symbols) {
  contract <- wp2_program_contract()
  rows <- lapply(seq_len(nrow(contract)), function(i) {
    genes <- contract$genes[[i]]
    detected <- genes[genes %in% gene_symbols]
    missing <- setdiff(genes, detected)
    data.frame(
      dataset_id = dataset_id,
      collapse_rule = collapse_rule,
      program_order = contract$program_order[[i]],
      program = contract$program_name[[i]],
      canonical_gene_count = 22L,
      detected_canonical_gene_count = length(detected),
      coverage_fraction = length(detected) / 22,
      detected_genes = paste(detected, collapse = ";"),
      missing_genes = paste(missing, collapse = ";"),
      substituted_genes = "",
      coverage_gate = if (length(detected) >= 18L) "PASS" else "FAIL",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

wp2_score_programs <- function(gene_expr, coverage_rows) {
  wp2_assert_token()
  if (any(coverage_rows$coverage_gate != "PASS")) stop("Coverage gate blocks scoring", call. = FALSE)
  contract <- wp2_program_contract()
  gene_sets <- setNames(lapply(contract$genes, intersect, y = rownames(gene_expr)), contract$program_name)
  param <- GSVA::ssgseaParam(
    exprData = gene_expr,
    geneSets = gene_sets,
    minSize = 5L,
    maxSize = 500L,
    normalize = TRUE
  )
  scored <- GSVA::gsva(param, verbose = FALSE, BPPARAM = BiocParallel::SerialParam(progressbar = FALSE))
  scored <- scored[WP2$programs, , drop = FALSE]
  out <- data.frame(sample = colnames(scored), t(scored), check.names = FALSE)
  if (anyDuplicated(out$sample)) stop("Duplicated scored sample", call. = FALSE)
  if (any(!is.finite(as.matrix(out[, WP2$programs, drop = FALSE])))) stop("Non-finite score", call. = FALSE)
  if (any(vapply(out[, WP2$programs, drop = FALSE], stats::sd, numeric(1)) <= 0)) stop("Zero-SD score", call. = FALSE)
  out
}

wp2_skewness <- function(x) {
  s <- stats::sd(x)
  if (!is.finite(s) || s == 0) return(NA_real_)
  mean((x - mean(x))^3) / s^3
}

wp2_tied_fraction <- function(x) max(tabulate(match(x, unique(x)))) / length(x)

wp2_score_qc <- function(dataset_id, primary, sensitivity) {
  rows <- lapply(seq_along(WP2$programs), function(i) {
    p <- primary[[WP2$programs[[i]]]]
    s <- sensitivity[[WP2$programs[[i]]]]
    q <- stats::quantile(p, c(.25, .5, .75), names = FALSE)
    data.frame(
      dataset_id = dataset_id,
      program_order = i,
      program = WP2$programs[[i]],
      n = length(p), mean = mean(p), SD = stats::sd(p), min = min(p), Q1 = q[[1L]], median = q[[2L]], Q3 = q[[3L]], max = max(p), IQR = stats::IQR(p),
      skewness = wp2_skewness(p), tied_fraction = wp2_tied_fraction(p), nonfinite_count = sum(!is.finite(p)),
      primary_vs_highest_MAD_Pearson = stats::cor(p, s, method = "pearson"),
      primary_vs_highest_MAD_Spearman = stats::cor(p, s, method = "spearman"),
      absolute_mean_difference = mean(abs(p - s)),
      rank_concordance = stats::cor(rank(p), rank(s), method = "pearson"),
      mapping_sensitivity_status = if (stats::cor(p, s, method = "spearman") >= .85) "PASS" else "LOW_CONCORDANCE",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

wp2_correlation_long <- function(scores, dataset_id, collapse_rule) {
  m <- stats::cor(as.matrix(scores[, WP2$programs, drop = FALSE]), method = "pearson")
  idx <- which(upper.tri(m), arr.ind = TRUE)
  long <- data.frame(
    dataset_id = dataset_id,
    collapse_rule = collapse_rule,
    program_1 = rownames(m)[idx[, 1]],
    program_2 = colnames(m)[idx[, 2]],
    program_1_order = idx[, 1],
    program_2_order = idx[, 2],
    pearson_r = m[idx],
    sign = ifelse(m[idx] > 0, "positive", ifelse(m[idx] < 0, "negative", "zero")),
    stringsAsFactors = FALSE
  )
  list(matrix = m, long = long)
}

wp2_discovery_edges <- function() {
  x <- utils::read.csv(WP2$paths$discovery_correlations, check.names = FALSE, stringsAsFactors = FALSE)
  x <- x[x$program_1_order < x$program_2_order, , drop = FALSE]
  x <- x[order(x$program_1_order, x$program_2_order), , drop = FALSE]
  if (nrow(x) != 15L) stop("Discovery edge authority must contain 15 edges", call. = FALSE)
  x
}

wp2_structural_metrics <- function(external_long) {
  discovery <- wp2_discovery_edges()
  external_long <- external_long[order(external_long$program_1_order, external_long$program_2_order), , drop = FALSE]
  rho <- suppressWarnings(stats::cor(discovery$pearson_r, external_long$pearson_r, method = "spearman"))
  signs <- sum(sign(discovery$pearson_r) == sign(external_long$pearson_r))
  status <- if (!is.finite(rho) || rho < .25 || signs <= 8L) {
    "WEAK_OR_NONREPLICATION"
  } else if (rho >= .50 && signs >= 12L) {
    "SUPPORTIVE_STRUCTURAL_REPLICATION"
  } else {
    "PARTIAL_REPLICATION"
  }
  c(matrix_vector_spearman = rho, sign_concordance_n = signs, structural_replication_status = status)
}

wp2_bootstrap_structure <- function(scores, dataset_id, replicates = WP2$bootstrap_replicates, seed = WP2$seed) {
  set.seed(seed + ifelse(dataset_id == "GSE10846", 10846L, 181063L))
  x <- as.matrix(scores[, WP2$programs, drop = FALSE])
  discovery <- wp2_discovery_edges()$pearson_r
  rho <- rep(NA_real_, replicates)
  signs <- rep(NA_integer_, replicates)
  failure <- character(replicates)
  for (b in seq_len(replicates)) {
    idx <- sample.int(nrow(x), nrow(x), replace = TRUE)
    cm <- tryCatch(stats::cor(x[idx, , drop = FALSE], method = "pearson"), error = identity)
    if (inherits(cm, "error") || any(!is.finite(cm))) {
      failure[[b]] <- if (inherits(cm, "error")) conditionMessage(cm) else "NONFINITE_CORRELATION"
      next
    }
    edge <- cm[upper.tri(cm)]
    # upper.tri follows column-major order; reorder by explicit edge indices.
    ij <- which(upper.tri(cm), arr.ind = TRUE)
    ord <- order(ij[, 1], ij[, 2])
    edge <- edge[ord]
    rho[[b]] <- suppressWarnings(stats::cor(discovery, edge, method = "spearman"))
    signs[[b]] <- sum(sign(discovery) == sign(edge))
  }
  ok <- is.finite(rho) & is.finite(signs)
  if (!any(ok)) stop("All structural bootstrap replicates failed", call. = FALSE)
  list(
    replicates = data.frame(dataset_id = dataset_id, replicate = seq_len(replicates), matrix_vector_spearman = rho, sign_concordance_n = signs, success = ok, failure_reason = failure, stringsAsFactors = FALSE),
    summary = data.frame(
      dataset_id = dataset_id,
      requested_replicates = replicates,
      successful_replicates = sum(ok),
      failed_replicates = sum(!ok),
      matrix_rho_CI_lower = stats::quantile(rho[ok], .025, names = FALSE),
      matrix_rho_CI_upper = stats::quantile(rho[ok], .975, names = FALSE),
      sign_concordance_CI_lower = stats::quantile(signs[ok], .025, names = FALSE),
      sign_concordance_CI_upper = stats::quantile(signs[ok], .975, names = FALSE),
      seed = seed + ifelse(dataset_id == "GSE10846", 10846L, 181063L),
      stringsAsFactors = FALSE
    )
  )
}

wp2_matrix_data_frame <- function(m) {
  data.frame(program = rownames(m), m, check.names = FALSE, stringsAsFactors = FALSE)
}

wp2_verify_input_registry <- function() {
  registry_path <- file.path(WP2$stage, "WP2_REAL_INPUT_REGISTRY.csv")
  registry <- utils::read.csv(registry_path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("input_id", "absolute_path", "sha256", "size", "read_authorized", "required")
  if (!all(required %in% names(registry))) stop("Input registry schema mismatch", call. = FALSE)
  rows <- lapply(seq_len(nrow(registry)), function(i) {
    path <- registry$absolute_path[[i]]
    exists <- file.exists(path)
    actual_size <- if (exists) as.numeric(file.info(path)$size) else NA_real_
    actual_sha <- if (exists) wp2_sha(path) else NA_character_
    data.frame(
      input_id = registry$input_id[[i]],
      exists = exists,
      expected_size = registry$size[[i]],
      actual_size = actual_size,
      size_match = exists && identical(actual_size, as.numeric(registry$size[[i]])),
      expected_sha256 = registry$sha256[[i]],
      actual_sha256 = actual_sha,
      sha256_match = exists && identical(actual_sha, registry$sha256[[i]]),
      read_authorized = isTRUE(registry$read_authorized[[i]]),
      required = isTRUE(registry$required[[i]]),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (any(!out$exists | !out$size_match | !out$sha256_match | !out$read_authorized | !out$required)) {
    stop("WP2 input registry gate failed", call. = FALSE)
  }
  out
}

wp2_verify_protected_paths <- function() {
  manifest <- utils::read.csv(
    file.path(WP2$stage, "WP2_PROTECTED_PATH_BASELINE_MANIFEST.csv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  rows <- lapply(seq_len(nrow(manifest)), function(i) {
    path <- file.path(WP2$root, manifest$relative_path[[i]])
    exists <- file.exists(path)
    actual_size <- if (exists) as.numeric(file.info(path)$size) else NA_real_
    actual_sha <- if (exists) wp2_sha(path) else NA_character_
    data.frame(
      relative_path = manifest$relative_path[[i]],
      exists = exists,
      size_match = exists && identical(actual_size, as.numeric(manifest$file_size_bytes[[i]])),
      sha256_match = exists && identical(actual_sha, manifest$sha256[[i]]),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (any(!out$exists | !out$size_match | !out$sha256_match)) {
    stop("Protected-path integrity gate failed", call. = FALSE)
  }
  out
}

wp2_other_r_process_count <- function() {
  task_lines <- suppressWarnings(system2("tasklist", c("/FO", "CSV", "/NH"), stdout = TRUE, stderr = TRUE))
  if (!length(task_lines)) return(NA_integer_)
  is_r <- grepl('^"(R|Rscript)\\.exe"', task_lines, ignore.case = TRUE)
  pids <- suppressWarnings(as.integer(sub('^"[^"]+","([0-9]+)".*$', '\\1', task_lines[is_r])))
  # Windows Rscript uses an Rscript.exe wrapper plus an R.exe worker.
  max(0L, sum(!is.na(pids)) - 2L)
}

wp2_loadings <- function() {
  p <- utils::read.csv(WP2$paths$projection_parameters, check.names = FALSE, stringsAsFactors = FALSE)
  p <- p[order(p$program_order), , drop = FALSE]
  if (!identical(p$program, WP2$programs)) stop("Projection parameter program order mismatch", call. = FALSE)
  m <- as.matrix(p[, c("PC1_loading", "PC2_loading")])
  rownames(m) <- p$program
  colnames(m) <- c("PC1", "PC2")
  list(parameters = p, matrix = m)
}

wp2_strategy_b <- function(scores, dataset_id) {
  x <- as.matrix(scores[, WP2$programs, drop = FALSE])
  z <- scale(x)
  if (any(!is.finite(z))) stop("Strategy B standardization nonfinite", call. = FALSE)
  load <- wp2_loadings()$matrix
  coords <- z %*% load
  coord_df <- data.frame(sample = scores$sample, PC1_B = coords[, 1], PC2_B = coords[, 2], coordinate_scale = "COHORT_RELATIVE", loading_direction = "WP1_FROZEN_UNFLIPPED", stringsAsFactors = FALSE)
  cors <- do.call(rbind, lapply(seq_along(WP2$programs), function(i) do.call(rbind, lapply(1:2, function(k) data.frame(
    dataset_id = dataset_id, program_order = i, program = WP2$programs[[i]], component = paste0("PC", k), pearson_r = stats::cor(x[, i], coords[, k]), loading = load[i, k], loading_sign_consistent = sign(stats::cor(x[, i], coords[, k])) == sign(load[i, k]), stringsAsFactors = FALSE
  )))))
  variance <- data.frame(dataset_id = dataset_id, component = c("PC1_B", "PC2_B"), projected_variance = apply(coords, 2, stats::var), coordinate_scale = "COHORT_RELATIVE", stringsAsFactors = FALSE)
  list(coordinates = coord_df, correlations = cors, variance = variance)
}

wp2_read_discovery_raw_scores <- function() {
  d <- utils::read.csv(WP2$paths$discovery_raw_scores, row.names = 1L, check.names = FALSE, stringsAsFactors = FALSE)
  old_to_new <- c(
    "Macrophage-rich niche" = WP2$programs[[1]],
    "T cell-inflamed niche" = WP2$programs[[2]],
    "Immune-inflamed / antigen-presentation niche" = WP2$programs[[3]],
    "Stromal/fibrotic niche" = WP2$programs[[4]],
    "Immune-cold / excluded niche" = WP2$programs[[5]],
    "Proliferative malignant B-cell niche" = WP2$programs[[6]]
  )
  mapped <- unname(old_to_new[rownames(d)])
  if (anyNA(mapped)) stop("Discovery raw score names cannot be mapped", call. = FALSE)
  m <- as.matrix(d)
  storage.mode(m) <- "double"
  rownames(m) <- mapped
  t(m[WP2$programs, , drop = FALSE])
}

wp2_strategy_a_diagnostics <- function(primary, sensitivity, coverage_primary) {
  p <- wp2_loadings()$parameters
  x <- as.matrix(primary[, WP2$programs, drop = FALSE])
  s <- as.matrix(sensitivity[, WP2$programs, drop = FALSE])
  train <- wp2_read_discovery_raw_scores()
  z <- sweep(sweep(x, 2, p$training_mean, "-"), 2, p$training_standard_deviation, "/")
  rows <- lapply(seq_along(WP2$programs), function(i) {
    external_iqr <- stats::IQR(x[, i])
    training_iqr <- stats::IQR(train[, i])
    median_z <- stats::median(z[, i])
    extreme <- mean(abs(z[, i]) > 5)
    tied <- wp2_tied_fraction(x[, i])
    sens_rho <- stats::cor(x[, i], s[, i], method = "spearman")
    data.frame(
      program_order = i, program = WP2$programs[[i]], detected_genes = coverage_primary$detected_canonical_gene_count[[i]],
      all_scores_finite = all(is.finite(x[, i])), score_SD = stats::sd(x[, i]), absolute_median_z = abs(median_z),
      external_IQR = external_iqr, training_IQR = training_iqr, IQR_ratio = external_iqr / training_iqr,
      extreme_abs_z_gt5_fraction = extreme, maximum_tied_score_fraction = tied,
      primary_sensitivity_spearman = sens_rho,
      coverage_pass_22 = coverage_primary$detected_canonical_gene_count[[i]] == 22L,
      finite_pass = all(is.finite(x[, i])), SD_pass = stats::sd(x[, i]) > 0,
      median_z_pass = abs(median_z) <= 2,
      IQR_ratio_pass = external_iqr / training_iqr >= .5 && external_iqr / training_iqr <= 2,
      extreme_z_pass = extreme <= .05,
      tied_fraction_pass = tied <= .10,
      mapping_sensitivity_pass = sens_rho >= .85,
      stringsAsFactors = FALSE
    )
  })
  diag <- do.call(rbind, rows)
  matrix_pass <- all(is.finite(z))
  covariance_pass <- all(is.finite(stats::cov(z)))
  gate_columns <- c("coverage_pass_22", "finite_pass", "SD_pass", "median_z_pass", "IQR_ratio_pass", "extreme_z_pass", "tied_fraction_pass", "mapping_sensitivity_pass")
  all_program_pass <- all(diag[, gate_columns, drop = FALSE])
  status <- if (matrix_pass && covariance_pass && all_program_pass) "PASS" else "BLOCKED"
  list(diagnostics = diag, z = z, matrix_pass = matrix_pass, covariance_pass = covariance_pass, status = status)
}

wp2_assert_no_forbidden_columns <- function(files) {
  forbidden <- c("selected_k", "cluster", "cluster_label", "assigned_class", "subtype", "taxonomy", "nearest_centroid", "winner_program", "survival", "hazard_ratio", "logrank", "cutpoint")
  allowed <- c("final_k_status", "taxonomy_status")
  hits <- list()
  for (path in files) {
    if (tolower(tools::file_ext(path)) != "csv") next
    cols <- names(utils::read.csv(path, nrows = 1L, check.names = FALSE, stringsAsFactors = FALSE))
    bad <- cols[vapply(tolower(cols), function(col) any(vapply(tolower(forbidden), function(x) grepl(x, col, fixed = TRUE), logical(1))), logical(1))]
    bad <- setdiff(bad, allowed)
    if (length(bad)) hits[[wp2_rel(path)]] <- bad
  }
  if (length(hits)) stop("Forbidden scientific output columns: ", paste(names(hits), vapply(hits, paste, collapse = ";", FUN.VALUE = character(1)), collapse = " | "), call. = FALSE)
  invisible(TRUE)
}
