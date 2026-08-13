DLBCL_PROJECT_ROOT <- normalizePath(
  Sys.getenv("DLBCL_PROJECT_ROOT", unset = getwd()),
  winslash = "/", mustWork = FALSE
)
DLBCL_SUPPLEMENTARY_CODE_ROOT <- normalizePath(
  Sys.getenv("DLBCL_SUPPLEMENTARY_CODE_ROOT",
    unset = file.path(DLBCL_PROJECT_ROOT, "Supplementary_Code")),
  winslash = "/", mustWork = FALSE
)

options(stringsAsFactors = FALSE)

R6R <- local({
  root <- DLBCL_PROJECT_ROOT
  revision <- file.path(root, "revision_2026_reviewer_response")
  module <- file.path(revision, "06r_bulk_orthogonal_tme_benchmark")
  list(
    root = root,
    revision = revision,
    module = module,
    token = "AUTHORIZE_06R_ORTHOGONAL_TME_BENCHMARK_SEED_20260807",
    human_authorization = "PROVIDED",
    seed = 20260807L,
    expected_mcp_repo = "https://github.com/ebecht/MCPcounter.git",
    expected_mcp_commit = "b6eac73e91c246fcff0bb1a5c68a816cd588fc48",
    local_lib = file.path(module, "00_environment/library"),
    mcp_src = file.path(module, "00_environment/src/MCPcounter_official_06r"),
    protocol = file.path(module, "00_protocol_freeze"),
    auth = file.path(module, "00_execution_authorization"),
    scripts = file.path(module, "01_execution_scripts"),
    outputs = file.path(module, "02_execution_outputs"),
    logs = file.path(module, "03_execution_logs"),
    program_contract = file.path(revision, "05x_wp1_continuous_score_freeze/WP1_CANONICAL_PROGRAM_CONTRACT.csv"),
    gse31312_expression = file.path(root, "02_processed_data/GSE31312_gene_expression_matrix.csv"),
    gse31312_scores = file.path(revision, "05e_stage4_GSE31312_execution_attempt2/01_score_space_validation/GSE31312_primary_score_matrix_498x6.csv"),
    gse10846_expression = file.path(revision, "05z_wp2_real_external_continuous_validation/local_only/GSE10846_primary_median_gene_expression.rds"),
    gse10846_scores = file.path(revision, "05z_wp2_real_external_continuous_validation/03_scores/GSE10846_PRIMARY_SIX_PROGRAM_SCORES.csv"),
    gse181063_expression = file.path(revision, "05z_wp2_real_external_continuous_validation/local_only/GSE181063_primary_median_gene_expression.rds"),
    gse181063_scores = file.path(revision, "05z_wp2_real_external_continuous_validation/03_scores/GSE181063_PRIMARY_SIX_PROGRAM_SCORES.csv"),
    wp2_input_registry = file.path(revision, "05z_wp2_real_external_continuous_validation/WP2_REAL_INPUT_REGISTRY.csv"),
    gse31312_input_registry = file.path(revision, "05e_stage4_GSE31312_execution_attempt2/00_run_control/STAGE4C1_ATTEMPT2_INPUT_HASH_VERIFICATION.csv"),
    excluded_gse10846 = c("GSM361239", "GSM361240", "GSM361241", "GSM361242", "GSM361243", "GSM361244"),
    protected_paths = c(
      "05x_wp1_continuous_score_freeze",
      "05z_wp2_real_external_continuous_validation",
      "06o_source_grounded_program_sensitivity",
      "06p_wp3_matched_null_and_depth_sensitivity",
      "06q_wp3_matched_null_balance_qc"
    )
  )
})

ensure_dirs <- function() {
  dirs <- c(
    R6R$auth, R6R$protocol, R6R$local_lib, dirname(R6R$mcp_src), R6R$scripts, R6R$logs,
    file.path(R6R$outputs, "input_qc"),
    file.path(R6R$outputs, "marker_sets"),
    file.path(R6R$outputs, "standard_mcp"),
    file.path(R6R$outputs, "disjoint_benchmark"),
    file.path(R6R$outputs, "correlations"),
    file.path(R6R$outputs, "predictability"),
    file.path(R6R$outputs, "cross_cohort"),
    file.path(R6R$outputs, "optional_existing_frameworks"),
    file.path(R6R$outputs, "validation")
  )
  invisible(vapply(dirs, dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))
}

require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Required package unavailable: ", pkg, call. = FALSE)
}

sha256 <- function(path) {
  require_pkg("digest")
  digest::digest(file = path, algo = "sha256")
}

norm_path <- function(path) gsub("\\\\", "/", normalizePath(path, winslash = "/", mustWork = FALSE))

rel_path <- function(path) {
  p <- norm_path(path)
  root <- paste0(norm_path(R6R$root), "/")
  if (startsWith(p, root)) substring(p, nchar(root) + 1L) else p
}

timestamp <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

write_csv_once <- function(x, path) {
  if (file.exists(path)) stop("Create-once output already exists: ", path, call. = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
  invisible(path)
}

write_text_once <- function(lines, path) {
  if (file.exists(path)) stop("Create-once output already exists: ", path, call. = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path, useBytes = TRUE)
  invisible(path)
}

file_record <- function(input_id, dataset_id, role, path, note = "") {
  exists <- file.exists(path)
  info <- if (exists) file.info(path) else NULL
  data.frame(
    input_id = input_id,
    dataset_id = dataset_id,
    role = role,
    absolute_path = norm_path(path),
    relative_path = rel_path(path),
    exists = exists,
    size_bytes = if (exists) as.numeric(info$size) else NA_real_,
    sha256 = if (exists) sha256(path) else NA_character_,
    mtime = if (exists) format(info$mtime, "%Y-%m-%dT%H:%M:%S%z") else NA_character_,
    notes = note,
    stringsAsFactors = FALSE
  )
}

git_status_record <- function() {
  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  setwd(R6R$revision)
  tracked <- system2("git", c("diff", "--name-only"), stdout = TRUE, stderr = TRUE)
  staged <- system2("git", c("diff", "--cached", "--name-only"), stdout = TRUE, stderr = TRUE)
  short <- system2("git", c("status", "--short"), stdout = TRUE, stderr = TRUE)
  data.frame(
    tracked_change_count = length(tracked[nzchar(tracked)]),
    staged_change_count = length(staged[nzchar(staged)]),
    status_short = paste(short, collapse = " | "),
    stringsAsFactors = FALSE
  )
}

load_program_contract <- function() {
  x <- utils::read.csv(R6R$program_contract, check.names = FALSE, stringsAsFactors = FALSE)
  x <- x[order(x$program_order), , drop = FALSE]
  x$gene_list <- strsplit(x$canonical_gene_members, ";", fixed = TRUE)
  if (!identical(as.integer(x$program_order), 1:6)) stop("Program order contract failed", call. = FALSE)
  if (any(lengths(x$gene_list) != 22L)) stop("Each canonical program must contain 22 memberships", call. = FALSE)
  x
}

program_id_map <- function() {
  pc <- load_program_contract()
  setNames(pc$program_id, pc$program_name)
}

canonical_gene_union <- function() {
  pc <- load_program_contract()
  sort(unique(unlist(pc$gene_list, use.names = FALSE)))
}

load_scores <- function(dataset_id, primary = TRUE) {
  pmap <- program_id_map()
  if (dataset_id == "GSE31312") {
    x <- utils::read.csv(R6R$gse31312_scores, check.names = FALSE, stringsAsFactors = FALSE)
    names(x)[1] <- "sample_id"
  } else if (dataset_id == "GSE10846") {
    x <- utils::read.csv(R6R$gse10846_scores, check.names = FALSE, stringsAsFactors = FALSE)
    names(x)[names(x) == "sample"] <- "sample_id"
  } else if (dataset_id == "GSE181063") {
    x <- utils::read.csv(R6R$gse181063_scores, check.names = FALSE, stringsAsFactors = FALSE)
    names(x)[names(x) == "sample"] <- "sample_id"
  } else {
    stop("Unknown dataset: ", dataset_id, call. = FALSE)
  }
  for (nm in names(pmap)) if (nm %in% names(x)) names(x)[names(x) == nm] <- pmap[[nm]]
  keep <- c("sample_id", unname(pmap))
  missing <- setdiff(keep, names(x))
  if (length(missing)) stop("Score table missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  x <- x[, keep, drop = FALSE]
  if (dataset_id == "GSE10846" && primary) x <- x[!x$sample_id %in% R6R$excluded_gse10846, , drop = FALSE]
  x
}

load_expression <- function(dataset_id, primary = TRUE) {
  require_pkg("data.table")
  if (dataset_id == "GSE31312") {
    dt <- data.table::fread(R6R$gse31312_expression, data.table = FALSE, check.names = FALSE, showProgress = FALSE)
    genes <- as.character(dt[[1]])
    m <- as.matrix(dt[, -1, drop = FALSE])
    storage.mode(m) <- "double"
    rownames(m) <- toupper(genes)
    colnames(m) <- names(dt)[-1]
  } else if (dataset_id == "GSE10846") {
    m <- readRDS(R6R$gse10846_expression)
    rownames(m) <- toupper(rownames(m))
    if (primary) m <- m[, !colnames(m) %in% R6R$excluded_gse10846, drop = FALSE]
  } else if (dataset_id == "GSE181063") {
    m <- readRDS(R6R$gse181063_expression)
    rownames(m) <- toupper(rownames(m))
  } else {
    stop("Unknown dataset: ", dataset_id, call. = FALSE)
  }
  if (anyDuplicated(rownames(m))) stop("Duplicated gene symbols in expression matrix for ", dataset_id, call. = FALSE)
  m
}

align_expression_scores <- function(dataset_id, primary = TRUE) {
  expr <- load_expression(dataset_id, primary = primary)
  scores <- load_scores(dataset_id, primary = primary)
  if (!setequal(colnames(expr), scores$sample_id)) {
    miss_e <- setdiff(scores$sample_id, colnames(expr))
    miss_s <- setdiff(colnames(expr), scores$sample_id)
    stop(dataset_id, " expression/score sample mismatch: missing_expr=", paste(miss_e, collapse = ";"),
         " missing_score=", paste(miss_s, collapse = ";"), call. = FALSE)
  }
  expr <- expr[, scores$sample_id, drop = FALSE]
  list(expression = expr, scores = scores)
}

read_mcp_genes <- function() {
  path <- file.path(R6R$mcp_src, "Signatures/genes.txt")
  if (!file.exists(path)) stop("Official MCP-counter genes.txt unavailable: ", path, call. = FALSE)
  x <- utils::read.table(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, quote = "\"")
  names(x)[names(x) == "HUGO symbols"] <- "gene_symbol"
  names(x)[names(x) == "Cell population"] <- "population"
  x$gene_symbol <- toupper(x$gene_symbol)
  x
}

marker_sets_from_genes <- function(mcp_genes) {
  split(mcp_genes$gene_symbol, mcp_genes$population)
}

long_marker_df <- function(marker_sets, set_type) {
  do.call(rbind, lapply(names(marker_sets), function(pop) {
    genes <- marker_sets[[pop]]
    if (!length(genes)) genes <- ""
    data.frame(set_type = set_type, population = pop, gene_symbol = genes, stringsAsFactors = FALSE)
  }))
}

cor_pair <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  if (n < 4L || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) {
    return(c(n_complete = n, pearson_r = NA_real_, pearson_ci_low = NA_real_, pearson_ci_high = NA_real_, spearman_rho = NA_real_))
  }
  r <- suppressWarnings(stats::cor(x[ok], y[ok], method = "pearson"))
  rho <- suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
  z <- atanh(max(min(r, 0.999999), -0.999999))
  se <- 1 / sqrt(n - 3)
  ci <- tanh(z + c(-1, 1) * 1.96 * se)
  c(n_complete = n, pearson_r = r, pearson_ci_low = ci[[1]], pearson_ci_high = ci[[2]], spearman_rho = rho)
}

row_z <- function(m) {
  sds <- apply(m, 1, stats::sd, na.rm = TRUE)
  means <- rowMeans(m, na.rm = TRUE)
  ok <- is.finite(sds) & sds > 0 & is.finite(means)
  z <- sweep(sweep(m[ok, , drop = FALSE], 1, means[ok], "-"), 1, sds[ok], "/")
  list(z = z, kept = rownames(m)[ok], excluded = rownames(m)[!ok])
}

condition_number <- function(x) {
  if (ncol(x) == 0L) return(NA_real_)
  out <- tryCatch(kappa(x, exact = TRUE), error = function(e) NA_real_)
  as.numeric(out)
}

design_rank <- function(x) qr(cbind(intercept = 1, x))$rank

safe_lm_predictability <- function(y, x) {
  ok <- is.finite(y) & apply(x, 1, function(z) all(is.finite(z)))
  y <- y[ok]
  x <- x[ok, , drop = FALSE]
  keep <- apply(x, 2, stats::sd) > 0
  x <- x[, keep, drop = FALSE]
  if (length(y) < 4L || ncol(x) < 1L) return(NULL)
  dat <- data.frame(y = y, x, check.names = FALSE)
  fit <- stats::lm(y ~ ., data = dat)
  pred <- stats::fitted(fit)
  r2 <- summary(fit)$r.squared
  adj <- summary(fit)$adj.r.squared
  paircor <- suppressWarnings(stats::cor(x, method = "pearson"))
  maxcor <- if (ncol(x) > 1L) max(abs(paircor[upper.tri(paircor)]), na.rm = TRUE) else NA_real_
  list(
    n = length(y), predictors = ncol(x), rank = design_rank(x), condition = condition_number(cbind(1, x)),
    r2 = r2, adj_r2 = adj, residual_sd = stats::sd(stats::residuals(fit)),
    residual_variance_fraction = 1 - r2,
    observed_fitted_pearson = suppressWarnings(stats::cor(y, pred, method = "pearson")),
    score_residual_pearson = suppressWarnings(stats::cor(y, stats::residuals(fit), method = "pearson")),
    max_abs_predictor_correlation = maxcor
  )
}

make_folds <- function(n, k, seed) {
  set.seed(seed)
  sample(rep(seq_len(k), length.out = n))
}

cv_lm <- function(y, x, folds) {
  oof <- rep(NA_real_, length(y))
  for (fold in sort(unique(folds))) {
    train <- folds != fold
    test <- folds == fold
    xtr <- x[train, , drop = FALSE]
    xte <- x[test, , drop = FALSE]
    keep <- apply(xtr, 2, stats::sd) > 0
    xtr <- xtr[, keep, drop = FALSE]
    xte <- xte[, keep, drop = FALSE]
    if (ncol(xtr) == 0L) {
      oof[test] <- mean(y[train])
      next
    }
    mu <- colMeans(xtr)
    sdv <- apply(xtr, 2, stats::sd)
    xtr <- sweep(sweep(xtr, 2, mu, "-"), 2, sdv, "/")
    xte <- sweep(sweep(xte, 2, mu, "-"), 2, sdv, "/")
    qr_x <- qr(xtr)
    piv <- qr_x$pivot[seq_len(qr_x$rank)]
    xtr <- xtr[, piv, drop = FALSE]
    xte <- xte[, piv, drop = FALSE]
    dat <- data.frame(y = y[train], xtr, check.names = FALSE)
    fit <- stats::lm(y ~ ., data = dat)
    oof[test] <- as.numeric(stats::predict(fit, newdata = data.frame(xte, check.names = FALSE)))
  }
  denom <- sum((y - mean(y))^2)
  cv_r2 <- 1 - sum((y - oof)^2) / denom
  list(
    oof = oof,
    cv_r2 = cv_r2,
    pearson = suppressWarnings(stats::cor(y, oof, method = "pearson")),
    spearman = suppressWarnings(stats::cor(y, oof, method = "spearman")),
    rmse = sqrt(mean((y - oof)^2))
  )
}

protected_snapshot <- function() {
  rows <- lapply(R6R$protected_paths, function(rel) {
    p <- file.path(R6R$revision, rel)
    files <- list.files(p, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
    files <- files[file.info(files)$isdir %in% FALSE]
    info <- file.info(files)
    data.frame(
      protected_path = rel,
      exists = dir.exists(p),
      file_count = length(files),
      total_size_bytes = sum(as.numeric(info$size), na.rm = TRUE),
      latest_mtime = if (length(files)) format(max(info$mtime, na.rm = TRUE), "%Y-%m-%dT%H:%M:%S%z") else NA_character_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
