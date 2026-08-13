DLBCL_PROJECT_ROOT <- normalizePath(
  Sys.getenv("DLBCL_PROJECT_ROOT", unset = getwd()),
  winslash = "/", mustWork = FALSE
)
DLBCL_SUPPLEMENTARY_CODE_ROOT <- normalizePath(
  Sys.getenv("DLBCL_SUPPLEMENTARY_CODE_ROOT",
    unset = file.path(DLBCL_PROJECT_ROOT, "Supplementary_Code")),
  winslash = "/", mustWork = FALSE
)

source(file.path(DLBCL_SUPPLEMENTARY_CODE_ROOT, "orthogonal_TME_benchmark", "00_06r_common.R"))
ensure_dirs()
.libPaths(c(R6R$local_lib, .libPaths()))
require_pkg("MCPcounter")
require_pkg("data.table")
require_pkg("matrixStats")

message("[06r science] Starting at ", timestamp())

write_csv_science <- function(x, path) {
  if (file.exists(path)) {
    tryCatch(utils::read.csv(path, nrows = 1L, check.names = FALSE), error = function(e) stop("Existing CSV is not parseable: ", path, call. = FALSE))
    message("[06r science] Reusing existing create-once output: ", path)
    return(invisible(path))
  }
  write_csv_once(x, path)
}

write_text_science <- function(lines, path) {
  if (file.exists(path)) {
    message("[06r science] Reusing existing create-once text output: ", path)
    return(invisible(path))
  }
  write_text_once(lines, path)
}

if (!file.exists(file.path(R6R$protocol, "06R_PROTOCOL_VALIDATION.csv"))) stop("Protocol validation missing", call. = FALSE)
pv <- utils::read.csv(file.path(R6R$protocol, "06R_PROTOCOL_VALIDATION.csv"), check.names = FALSE, stringsAsFactors = FALSE)
if (any(pv$status != "PASS")) stop("Protocol validation did not pass", call. = FALSE)

programs <- load_program_contract()
program_ids <- programs$program_id
names(program_ids) <- programs$program_name
canonical_union <- canonical_gene_union()
mcp_genes <- read_mcp_genes()
mcp_genes_for_mcp <- mcp_genes
names(mcp_genes_for_mcp)[names(mcp_genes_for_mcp) == "gene_symbol"] <- "HUGO symbols"
names(mcp_genes_for_mcp)[names(mcp_genes_for_mcp) == "population"] <- "Cell population"
official_sets <- marker_sets_from_genes(mcp_genes)
disjoint_sets <- lapply(official_sets, setdiff, y = canonical_union)

write_csv_science(long_marker_df(disjoint_sets, "GLOBAL_DISJOINT_FROM_CANONICAL_121"),
               file.path(R6R$outputs, "marker_sets/06R_MCP_MARKER_SETS_GLOBAL_DISJOINT.csv"))

overlap_rows <- do.call(rbind, lapply(names(official_sets), function(pop) {
  mcp <- unique(official_sets[[pop]])
  by_program <- do.call(rbind, lapply(seq_len(nrow(programs)), function(i) {
    genes <- programs$gene_list[[i]]
    overlap <- intersect(mcp, genes)
    data.frame(
      comparison_scope = "PROGRAM",
      population = pop,
      program_id = programs$program_id[[i]],
      program_name = programs$program_name[[i]],
      mcp_marker_set_size = length(mcp),
      curated_program_size = 22L,
      overlap_gene_count = length(overlap),
      overlap_genes = paste(overlap, collapse = ";"),
      jaccard_index = length(overlap) / length(union(mcp, genes)),
      proportion_mcp_markers_overlapping_program = length(overlap) / length(mcp),
      proportion_program_overlapping_mcp_markers = length(overlap) / 22,
      stringsAsFactors = FALSE
    )
  }))
  union_overlap <- intersect(mcp, canonical_union)
  rbind(by_program, data.frame(
    comparison_scope = "CANONICAL_121_UNION",
    population = pop,
    program_id = "canonical_121_union",
    program_name = "Canonical 121-gene union",
    mcp_marker_set_size = length(mcp),
    curated_program_size = length(canonical_union),
    overlap_gene_count = length(union_overlap),
    overlap_genes = paste(union_overlap, collapse = ";"),
    jaccard_index = length(union_overlap) / length(union(mcp, canonical_union)),
    proportion_mcp_markers_overlapping_program = length(union_overlap) / length(mcp),
    proportion_program_overlapping_mcp_markers = length(union_overlap) / length(canonical_union),
    stringsAsFactors = FALSE
  ))
}))
write_csv_science(overlap_rows, file.path(R6R$outputs, "marker_sets/06R_MCP_PROGRAM_GENE_OVERLAP.csv"))

disjoint_summary <- do.call(rbind, lapply(names(official_sets), function(pop) {
  removed <- intersect(official_sets[[pop]], canonical_union)
  remain <- disjoint_sets[[pop]]
  data.frame(
    population = pop,
    original_marker_count = length(official_sets[[pop]]),
    canonical_overlap_count = length(removed),
    remaining_global_disjoint_marker_count = length(remain),
    removed_genes = paste(removed, collapse = ";"),
    remaining_genes = paste(remain, collapse = ";"),
    stringsAsFactors = FALSE
  )
}))

datasets_primary <- c("GSE31312", "GSE10846", "GSE181063")
standard_long <- list()
standard_wide <- list()
mcp_qc <- list()
disjoint_long <- list()
disjoint_wide <- list()
disjoint_qc <- list()
score_store <- list()

for (ds in datasets_primary) {
  message("[06r science] Standard and disjoint scoring: ", ds)
  aligned <- align_expression_scores(ds, primary = TRUE)
  expr <- aligned$expression
  scores <- aligned$scores
  estimate <- MCPcounter::MCPcounter.estimate(expr, featuresType = "HUGO_symbols", genes = mcp_genes_for_mcp)
  estimate <- as.matrix(estimate)
  standard_wide[[ds]] <- data.frame(dataset = ds, sample_id = colnames(estimate), t(estimate), check.names = FALSE)
  standard_long[[ds]] <- do.call(rbind, lapply(rownames(estimate), function(pop) {
    data.frame(dataset = ds, sample_id = colnames(estimate), population = pop, mcp_score = as.numeric(estimate[pop, ]), stringsAsFactors = FALSE)
  }))
  mcp_qc[[ds]] <- do.call(rbind, lapply(names(official_sets), function(pop) {
    official <- unique(official_sets[[pop]])
    present <- intersect(official, rownames(expr))
    vals <- if (pop %in% rownames(estimate)) estimate[pop, ] else rep(NA_real_, ncol(expr))
    data.frame(
      dataset = ds,
      population = pop,
      official_marker_count = length(official),
      genes_present_count = length(present),
      genes_missing_count = length(setdiff(official, present)),
      genes_present = paste(present, collapse = ";"),
      genes_missing = paste(setdiff(official, present), collapse = ";"),
      duplicate_gene_handling = "input matrices require unique gene symbols",
      platform_specific_mapping = "HUGO gene-symbol expression matrix",
      nonfinite_score_count = sum(!is.finite(vals)),
      score_sd = stats::sd(vals, na.rm = TRUE),
      status = if (all(is.finite(vals)) && stats::sd(vals) > 0) "PASS" else "FAIL",
      stringsAsFactors = FALSE
    )
  }))
  dj_scores <- list()
  dq <- list()
  for (pop in names(disjoint_sets)) {
    official <- unique(official_sets[[pop]])
    disjoint <- unique(disjoint_sets[[pop]])
    detected <- intersect(disjoint, rownames(expr))
    finite_detected <- detected[apply(expr[detected, , drop = FALSE], 1, function(x) all(is.finite(x)))]
    z <- row_z(expr[finite_detected, , drop = FALSE])
    eligible <- z$kept
    status <- if (length(eligible) >= 5L) "EVALUABLE" else "NOT_EVALUABLE_TOO_FEW_DISJOINT_MARKERS"
    vals <- if (status == "EVALUABLE") colMeans(z$z, na.rm = TRUE) else rep(NA_real_, ncol(expr))
    dj_scores[[pop]] <- vals
    dq[[pop]] <- data.frame(
      dataset = ds,
      population = pop,
      original_marker_count = length(official),
      canonical_overlap_count = length(intersect(official, canonical_union)),
      global_disjoint_marker_count = length(disjoint),
      detected_disjoint_marker_count = length(detected),
      eligible_detected_finite_nonzero_variance_disjoint_marker_count = length(eligible),
      eligible_genes = paste(eligible, collapse = ";"),
      excluded_detected_zero_variance_or_nonfinite_genes = paste(z$excluded, collapse = ";"),
      status = status,
      stringsAsFactors = FALSE
    )
  }
  dj_mat <- do.call(cbind, dj_scores)
  colnames(dj_mat) <- names(dj_scores)
  rownames(dj_mat) <- colnames(expr)
  disjoint_wide[[ds]] <- data.frame(dataset = ds, sample_id = rownames(dj_mat), dj_mat, check.names = FALSE)
  disjoint_long[[ds]] <- do.call(rbind, lapply(colnames(dj_mat), function(pop) {
    data.frame(dataset = ds, sample_id = rownames(dj_mat), population = pop, disjoint_score = as.numeric(dj_mat[, pop]), stringsAsFactors = FALSE)
  }))
  disjoint_qc[[ds]] <- do.call(rbind, dq)
  score_store[[ds]] <- list(program = scores, standard = standard_wide[[ds]], disjoint = disjoint_wide[[ds]])
}

write_csv_science(do.call(rbind, standard_long), file.path(R6R$outputs, "standard_mcp/06R_STANDARD_MCP_SCORES_LONG.csv"))
write_csv_science(do.call(rbind, standard_wide), file.path(R6R$outputs, "standard_mcp/06R_STANDARD_MCP_SCORES_WIDE.csv"))
write_csv_science(do.call(rbind, mcp_qc), file.path(R6R$outputs, "input_qc/06R_MCP_SCORE_QC.csv"))
write_csv_science(disjoint_summary, file.path(R6R$outputs, "marker_sets/06R_MCP_GLOBAL_DISJOINT_MARKER_SET_SUMMARY.csv"))
write_csv_science(do.call(rbind, disjoint_long), file.path(R6R$outputs, "disjoint_benchmark/06R_DISJOINT_MARKER_SCORES_LONG.csv"))
write_csv_science(do.call(rbind, disjoint_wide), file.path(R6R$outputs, "disjoint_benchmark/06R_DISJOINT_MARKER_SCORES_WIDE.csv"))
write_csv_science(do.call(rbind, disjoint_qc), file.path(R6R$outputs, "disjoint_benchmark/06R_DISJOINT_SCORE_QC.csv"))

correlation_rows <- function(ds, bench_df, bench_cols, value_type) {
  prog <- score_store[[ds]]$program
  bench <- bench_df[bench_df$dataset == ds, , drop = FALSE]
  bench <- bench[match(prog$sample_id, bench$sample_id), , drop = FALSE]
  do.call(rbind, lapply(programs$program_id, function(pid) do.call(rbind, lapply(bench_cols, function(pop) {
    v <- cor_pair(prog[[pid]], bench[[pop]])
    data.frame(dataset = ds, program_id = pid, program_name = programs$program_name[match(pid, programs$program_id)], population = pop,
               n_complete = v[["n_complete"]], pearson_r = v[["pearson_r"]], pearson_ci_low = v[["pearson_ci_low"]],
               pearson_ci_high = v[["pearson_ci_high"]], spearman_rho = v[["spearman_rho"]],
               benchmark_type = value_type,
               evaluation_status = if (is.finite(v[["pearson_r"]])) "PASS" else "NOT_EVALUABLE", stringsAsFactors = FALSE)
  }))))
}

std_all <- do.call(rbind, standard_wide)
dj_all <- do.call(rbind, disjoint_wide)
pop_cols <- names(official_sets)
std_cor <- do.call(rbind, lapply(datasets_primary, correlation_rows, bench_df = std_all, bench_cols = pop_cols, value_type = "STANDARD_MCP_COUNTER"))
dj_cor <- do.call(rbind, lapply(datasets_primary, correlation_rows, bench_df = dj_all, bench_cols = pop_cols, value_type = "GLOBAL_DISJOINT_MARKER_SCORE"))
write_csv_science(std_cor, file.path(R6R$outputs, "correlations/06R_STANDARD_MCP_PROGRAM_CORRELATIONS.csv"))
write_csv_science(dj_cor, file.path(R6R$outputs, "correlations/06R_DISJOINT_MCP_PROGRAM_CORRELATIONS.csv"))

delta <- merge(std_cor, dj_cor, by = c("dataset", "program_id", "program_name", "population"), suffixes = c("_standard", "_disjoint"), all = TRUE)
delta$delta_abs_pearson <- abs(delta$pearson_r_standard) - abs(delta$pearson_r_disjoint)
delta$delta_abs_spearman <- abs(delta$spearman_rho_standard) - abs(delta$spearman_rho_disjoint)
write_csv_science(delta, file.path(R6R$outputs, "correlations/06R_STANDARD_VS_DISJOINT_CORRELATION_CHANGE.csv"))

message("[06r science] GSE10846 all-420 sensitivity")
aligned420 <- align_expression_scores("GSE10846", primary = FALSE)
expr420 <- aligned420$expression
scores420 <- aligned420$scores
est420 <- MCPcounter::MCPcounter.estimate(expr420, featuresType = "HUGO_symbols", genes = mcp_genes_for_mcp)
est420 <- as.matrix(est420)
dj420 <- list()
for (pop in names(disjoint_sets)) {
  detected <- intersect(disjoint_sets[[pop]], rownames(expr420))
  z <- row_z(expr420[detected, , drop = FALSE])
  dj420[[pop]] <- if (length(z$kept) >= 5L) colMeans(z$z, na.rm = TRUE) else rep(NA_real_, ncol(expr420))
}
dj420 <- do.call(cbind, dj420)
colnames(dj420) <- names(disjoint_sets)
rownames(dj420) <- colnames(expr420)
sens_rows <- rbind(
  do.call(rbind, lapply(programs$program_id, function(pid) do.call(rbind, lapply(rownames(est420), function(pop) {
    v <- cor_pair(scores420[[pid]], as.numeric(est420[pop, scores420$sample_id]))
    data.frame(dataset = "GSE10846", cohort_definition = "ALL420_SENSITIVITY", excluded_from_primary_profiles = paste(R6R$excluded_gse10846, collapse = ";"), benchmark_type = "STANDARD_MCP_COUNTER", program_id = pid, program_name = programs$program_name[match(pid, programs$program_id)], population = pop, n_complete = v[["n_complete"]], pearson_r = v[["pearson_r"]], pearson_ci_low = v[["pearson_ci_low"]], pearson_ci_high = v[["pearson_ci_high"]], spearman_rho = v[["spearman_rho"]], stringsAsFactors = FALSE)
  })))),
  do.call(rbind, lapply(programs$program_id, function(pid) do.call(rbind, lapply(colnames(dj420), function(pop) {
    v <- cor_pair(scores420[[pid]], dj420[scores420$sample_id, pop])
    data.frame(dataset = "GSE10846", cohort_definition = "ALL420_SENSITIVITY", excluded_from_primary_profiles = paste(R6R$excluded_gse10846, collapse = ";"), benchmark_type = "GLOBAL_DISJOINT_MARKER_SCORE", program_id = pid, program_name = programs$program_name[match(pid, programs$program_id)], population = pop, n_complete = v[["n_complete"]], pearson_r = v[["pearson_r"]], pearson_ci_low = v[["pearson_ci_low"]], pearson_ci_high = v[["pearson_ci_high"]], spearman_rho = v[["spearman_rho"]], stringsAsFactors = FALSE)
  }))))
)
write_csv_science(sens_rows, file.path(R6R$outputs, "correlations/06R_GSE10846_ALL420_SENSITIVITY.csv"))

predict_rows <- list()
cv_rows <- list()
oof_rows <- list()
for (ds in datasets_primary) {
  prog <- score_store[[ds]]$program
  disj <- score_store[[ds]]$disjoint
  disj <- disj[match(prog$sample_id, disj$sample_id), , drop = FALSE]
  x <- as.matrix(disj[, pop_cols, drop = FALSE])
  mode(x) <- "numeric"
  finite_cols <- apply(x, 2, function(z) all(is.finite(z)) && stats::sd(z) > 0)
  x <- x[, finite_cols, drop = FALSE]
  folds <- make_folds(nrow(prog), 10L, R6R$seed + match(ds, datasets_primary))
  for (pid in programs$program_id) {
    y <- prog[[pid]]
    fit <- safe_lm_predictability(y, x)
    if (is.null(fit)) {
      predict_rows[[paste(ds, pid)]] <- data.frame(dataset = ds, program_id = pid, program_name = programs$program_name[match(pid, programs$program_id)], n = length(y), benchmark_predictors = ncol(x), design_rank = NA_integer_, condition_number = NA_real_, r_squared = NA_real_, adjusted_r_squared = NA_real_, residual_sd = NA_real_, residual_variance_fraction = NA_real_, observed_fitted_pearson = NA_real_, score_residual_pearson = NA_real_, max_abs_predictor_correlation = NA_real_, model_status = "NOT_EVALUABLE", stringsAsFactors = FALSE)
    } else {
      predict_rows[[paste(ds, pid)]] <- data.frame(dataset = ds, program_id = pid, program_name = programs$program_name[match(pid, programs$program_id)], n = fit$n, benchmark_predictors = fit$predictors, design_rank = fit$rank, condition_number = fit$condition, r_squared = fit$r2, adjusted_r_squared = fit$adj_r2, residual_sd = fit$residual_sd, residual_variance_fraction = fit$residual_variance_fraction, observed_fitted_pearson = fit$observed_fitted_pearson, score_residual_pearson = fit$score_residual_pearson, max_abs_predictor_correlation = fit$max_abs_predictor_correlation, model_status = "PASS", stringsAsFactors = FALSE)
    }
    cv <- cv_lm(y, x, folds)
    cv_rows[[paste(ds, pid)]] <- data.frame(dataset = ds, program_id = pid, program_name = programs$program_name[match(pid, programs$program_id)], n = length(y), folds = 10L, seed = R6R$seed, benchmark_predictors = ncol(x), CV_R2 = cv$cv_r2, OOF_Pearson = cv$pearson, OOF_Spearman = cv$spearman, RMSE = cv$rmse, cv_status = "PASS", stringsAsFactors = FALSE)
    oof_rows[[paste(ds, pid)]] <- data.frame(dataset = ds, sample_id = prog$sample_id, program_id = pid, observed_program_score = y, oof_prediction = cv$oof, fold = folds, stringsAsFactors = FALSE)
  }
}
pred <- do.call(rbind, predict_rows)
cv <- do.call(rbind, cv_rows)
write_csv_science(pred, file.path(R6R$outputs, "predictability/06R_PROGRAM_ABUNDANCE_PREDICTABILITY.csv"))
write_csv_science(cv, file.path(R6R$outputs, "predictability/06R_PROGRAM_ABUNDANCE_10FOLD_CV.csv"))
write_csv_science(do.call(rbind, oof_rows), file.path(R6R$outputs, "predictability/06R_PROGRAM_ABUNDANCE_OOF_PREDICTIONS.csv"))

strongest <- function(cor_df, type_label) {
  ok <- is.finite(cor_df$pearson_r)
  cor_df <- cor_df[ok, , drop = FALSE]
  do.call(rbind, lapply(split(cor_df, paste(cor_df$dataset, cor_df$program_id, sep = "|")), function(d) {
    d <- d[which.max(abs(d$pearson_r)), , drop = FALSE]
    data.frame(dataset = d$dataset, program_id = d$program_id, program_name = d$program_name, benchmark_type = type_label, strongest_population = d$population, strongest_pearson_r = d$pearson_r, strongest_spearman_rho = d$spearman_rho, stringsAsFactors = FALSE)
  }))
}
std_strong <- strongest(std_cor, "STANDARD_MCP_COUNTER")
dj_strong <- strongest(dj_cor, "GLOBAL_DISJOINT_MARKER_SCORE")
summary_rows <- merge(std_strong, dj_strong, by = c("dataset", "program_id", "program_name"), suffixes = c("_standard", "_disjoint"), all = TRUE)
summary_rows <- merge(summary_rows, pred[, c("dataset", "program_id", "adjusted_r_squared", "residual_variance_fraction")], by = c("dataset", "program_id"), all.x = TRUE)
summary_rows <- merge(summary_rows, cv[, c("dataset", "program_id", "CV_R2")], by = c("dataset", "program_id"), all.x = TRUE)
write_csv_science(summary_rows, file.path(R6R$outputs, "cross_cohort/06R_PROGRAM_LEVEL_BENCHMARK_SUMMARY.csv"))

cohort_pairs <- list(c("GSE31312", "GSE10846"), c("GSE31312", "GSE181063"), c("GSE10846", "GSE181063"))
structure_rows <- list()
bench_list <- list(STANDARD_MCP_COUNTER = std_cor, GLOBAL_DISJOINT_MARKER_SCORE = dj_cor)
for (nm in names(bench_list)) {
  df <- bench_list[[nm]]
  for (pair in cohort_pairs) {
    a <- df[df$dataset == pair[[1]], c("program_id", "population", "pearson_r")]
    b <- df[df$dataset == pair[[2]], c("program_id", "population", "pearson_r")]
    m <- merge(a, b, by = c("program_id", "population"), suffixes = c("_a", "_b"))
    ok <- is.finite(m$pearson_r_a) & is.finite(m$pearson_r_b)
    structure_rows[[paste(nm, pair, collapse = "|")]] <- data.frame(
      benchmark_type = nm,
      cohort_a = pair[[1]],
      cohort_b = pair[[2]],
      vector_length = sum(ok),
      pearson_correlation = suppressWarnings(stats::cor(m$pearson_r_a[ok], m$pearson_r_b[ok], method = "pearson")),
      spearman_correlation = suppressWarnings(stats::cor(m$pearson_r_a[ok], m$pearson_r_b[ok], method = "spearman")),
      stringsAsFactors = FALSE
    )
  }
}
write_csv_science(do.call(rbind, structure_rows), file.path(R6R$outputs, "cross_cohort/06R_CROSS_COHORT_CORRELATION_STRUCTURE.csv"))

estimate_status <- data.frame(method = "ESTIMATE", status = if (requireNamespace("estimate", quietly = TRUE)) "AVAILABLE_NOT_RUN_BY_06R_IMPLEMENTATION" else "OPTIONAL_ESTIMATE_NOT_RUN", notes = "06r primary module is complete with MCP-counter alone.", local_hits = "", stringsAsFactors = FALSE)
kotlov_hits <- list.files(R6R$revision, pattern = "Kotlov|LME|lymphoma.microenvironment|lymphoma_microenvironment", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
kotlov_rel_hits <- if (length(kotlov_hits)) vapply(kotlov_hits, rel_path, character(1)) else character()
kotlov_status <- data.frame(method = "Kotlov continuous LME signatures", status = if (length(kotlov_hits)) "KOTLOV_CANDIDATE_FILES_FOUND_NOT_RUN_SOURCE_NOT_LOCALLY_VERIFIED" else "KOTLOV_CONTINUOUS_COMPARISON_NOT_RUN_SOURCE_NOT_LOCALLY_VERIFIED", notes = "Exact published Kotlov gene-set source was not locally verified by 06r, so no comparison was run.", local_hits = paste(kotlov_rel_hits, collapse = ";"), stringsAsFactors = FALSE)
write_csv_science(rbind(estimate_status, kotlov_status), file.path(R6R$outputs, "optional_existing_frameworks/06R_OPTIONAL_METHOD_STATUS.csv"))

report <- c(
  "# 06R Final Report",
  "",
  "Final module status is assigned by the validation script.",
  "",
  "## Scientific Rationale",
  "This module benchmarks the six frozen continuous curated programs against broad MCP-counter immune and stromal abundance scores, then repeats key comparisons using MCP-counter-derived marker scores made globally disjoint from all 121 canonical genes.",
  "",
  "## Cohort Definitions",
  "Primary cohorts: GSE31312 n=498, GSE10846 clinical biopsy n=414 after excluding GSM361239-GSM361240-GSM361241-GSM361242-GSM361243-GSM361244, and GSE181063 n=1310. GSE10846 all 420 profiles are reported as sensitivity only.",
  "",
  "## MCP-counter Provenance",
  paste0("Official repository: ", R6R$expected_mcp_repo),
  paste0("Pinned commit: ", R6R$expected_mcp_commit),
  "Package version: 1.2.0",
  "",
  "## Bounded Interpretation",
  "The outputs quantify concordance between curated program scores and broad microenvironment abundance estimates. Persistence after exact gene-overlap removal is descriptive and does not prove causal independence, biological novelty, cell-of-origin, prognostic utility, or a discrete ecosystem taxonomy.",
  "",
  "## Output Locations",
  "Major CSV outputs are under 02_execution_outputs/standard_mcp, marker_sets, disjoint_benchmark, correlations, predictability, cross_cohort, optional_existing_frameworks, and validation."
)
write_text_science(report, file.path(R6R$outputs, "06R_FINAL_REPORT.md"))

message("[06r science] Completed at ", timestamp())
