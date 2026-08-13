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

O6_PROJECT <- DLBCL_PROJECT_ROOT
O6_REVISION <- file.path(O6_PROJECT, "revision_2026_reviewer_response")
O6_ROOT <- file.path(O6_REVISION, "06o_source_grounded_program_sensitivity")
O6_PROTOCOL <- file.path(O6_ROOT, "00_protocol_freeze")
O6_AMENDMENT <- file.path(O6_PROTOCOL, "amendments/AMENDMENT_001_06O_SEED_AND_FDR_FREEZE")
O6_SCRIPTS <- file.path(O6_ROOT, "01_execution_scripts")
O6_OUTPUTS <- file.path(O6_ROOT, "01_execution_outputs")
O6_LOGS <- file.path(O6_ROOT, "02_execution_logs")
O6_FROZEN_LIBRARY <- file.path(O6_REVISION,
  "04b_stage4_environment_freeze/stage4_renv_project/renv/library/windows/R-4.5/x86_64-w64-mingw32")
O6_WP3_LIBRARY <- file.path(O6_REVISION, "06b_wp3b_spatial_scope_method_resolution/.wp3_r_library")
O6_TOKEN <- "AUTHORIZE_06O_SOURCE_GROUNDED_CORE_EXECUTION_SEED_20260805"
O6_MASTER_SEED <- 20260805L
O6_PROGRAM_IDS <- c(
  "Macrophage-rich" = "macrophage_rich",
  "T cell-inflamed" = "t_cell_inflamed",
  "Antigen-presentation" = "antigen_presentation",
  "Stromal/fibrotic" = "stromal_fibrotic",
  "Immune-cold/exclusion" = "immune_cold_exclusion",
  "Proliferative/cycling" = "proliferative_cycling"
)

o6_set_library <- function() {
  stopifnot(dir.exists(O6_FROZEN_LIBRARY), dir.exists(O6_WP3_LIBRARY))
  .libPaths(unique(c(O6_WP3_LIBRARY, O6_FROZEN_LIBRARY, .Library)))
  invisible(.libPaths())
}

o6_require_token <- function() {
  path <- file.path(O6_ROOT, "00_execution_authorization/06O_EXECUTION_TOKEN.txt")
  if (!file.exists(path) || !identical(readChar(path, file.info(path)$size, useBytes = TRUE), O6_TOKEN)) {
    stop("Exact 06o execution token is absent or invalid", call. = FALSE)
  }
  invisible(TRUE)
}

o6_require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Required frozen packages are unavailable: ", paste(missing, collapse = ";"), call. = FALSE)
  invisible(TRUE)
}

o6_write_csv_once <- function(x, path) {
  if (file.exists(path)) stop("Create-once refusal: ", path, call. = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  invisible(path)
}

o6_write_text_once <- function(text, path) {
  if (file.exists(path)) stop("Create-once refusal: ", path, call. = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(text, path, useBytes = TRUE)
  invisible(path)
}

o6_sha256 <- function(path) {
  o6_require_packages("digest")
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

o6_gene_sets <- function() {
  authority <- utils::read.csv(file.path(O6_PROTOCOL, "06O_GENESET_AUTHORITY.csv"), check.names = FALSE)
  stopifnot(nrow(authority) == 132L, all(authority$program %in% names(O6_PROGRAM_IDS)))
  authority$gene <- toupper(authority$gene)
  full <- split(authority$gene, factor(authority$program, levels = names(O6_PROGRAM_IDS)))
  names(full) <- unname(O6_PROGRAM_IDS[names(full)])
  core_rows <- authority[authority$direct_source_core_membership == "TRUE", , drop = FALSE]
  core <- split(core_rows$gene, factor(core_rows$program,
    levels = setdiff(names(O6_PROGRAM_IDS), "Immune-cold/exclusion")))
  names(core) <- unname(O6_PROGRAM_IDS[names(core)])
  expected <- c(macrophage_rich = 17L, t_cell_inflamed = 15L,
    antigen_presentation = 8L, stromal_fibrotic = 16L, proliferative_cycling = 9L)
  stopifnot(identical(lengths(full), setNames(rep(22L, 6L), unname(O6_PROGRAM_IDS))),
    identical(lengths(core), expected))
  list(authority = authority, full = full, core = core)
}

o6_program_name <- function(program_id) names(O6_PROGRAM_IDS)[match(program_id, O6_PROGRAM_IDS)]

o6_parse_series_matrix <- function(path, expected_samples) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  tab <- utils::read.delim(con, header = TRUE, quote = "\"", comment.char = "!",
    check.names = FALSE, stringsAsFactors = FALSE)
  if (!nrow(tab) || ncol(tab) != expected_samples + 1L) stop("Series matrix dimensions are invalid: ", path)
  probes <- as.character(tab[[1L]])
  if (anyDuplicated(probes) || any(!nzchar(probes))) stop("Invalid probe identifiers in ", path)
  expr <- as.matrix(tab[, -1L, drop = FALSE])
  storage.mode(expr) <- "double"
  rownames(expr) <- probes
  if (any(!is.finite(expr))) stop("Non-finite expression value in ", path)
  expr
}

o6_read_annotation <- function(path, expression_probes) {
  ann <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("PROBEID", "SYMBOL") %in% names(ann))) stop("Annotation lacks PROBEID/SYMBOL: ", path)
  ann <- ann[, c("PROBEID", "SYMBOL"), drop = FALSE]
  ann$PROBEID <- trimws(as.character(ann$PROBEID))
  ann$SYMBOL <- toupper(trimws(as.character(ann$SYMBOL)))
  ann <- ann[!is.na(ann$PROBEID) & nzchar(ann$PROBEID) & !is.na(ann$SYMBOL) &
    nzchar(ann$SYMBOL) & ann$PROBEID %in% expression_probes, , drop = FALSE]
  unique(ann)
}

o6_collapse_highest <- function(expr_probe, annotation, rule = c("highest_mean", "highest_mad")) {
  rule <- match.arg(rule)
  metric <- if (rule == "highest_mean") {
    rowMeans(expr_probe[annotation$PROBEID, , drop = FALSE], na.rm = TRUE)
  } else {
    apply(expr_probe[annotation$PROBEID, , drop = FALSE], 1L, stats::mad, na.rm = TRUE)
  }
  annotation$selection_metric <- metric
  annotation <- annotation[order(annotation$SYMBOL, -annotation$selection_metric, annotation$PROBEID), , drop = FALSE]
  selected <- annotation[!duplicated(annotation$SYMBOL), , drop = FALSE]
  expr_gene <- expr_probe[selected$PROBEID, , drop = FALSE]
  rownames(expr_gene) <- selected$SYMBOL
  list(expression = expr_gene, selected = selected, rule = rule)
}

o6_ssgsea <- function(expression, gene_sets) {
  o6_require_packages(c("GSVA", "BiocParallel"))
  sets <- lapply(gene_sets, function(g) intersect(g, rownames(expression)))
  param <- GSVA::ssgseaParam(exprData = expression, geneSets = sets,
    minSize = 5L, maxSize = 500L, alpha = 0.25, normalize = TRUE)
  as.matrix(GSVA::gsva(param, BPPARAM = BiocParallel::SerialParam(), verbose = FALSE))
}

o6_safe_cor <- function(x, y, method) {
  if (length(x) != length(y) || length(x) < 3L || stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
  stats::cor(x, y, method = method, use = "complete.obs")
}

o6_read_tsv_gz <- function(path) utils::read.delim(gzfile(path), header = FALSE,
  quote = "", comment.char = "", check.names = FALSE)

o6_read_spatial_counts <- function(area) {
  base <- file.path(O6_PROJECT, "00_raw_data/GSE276542/standard_10x", area, "filtered_feature_bc_matrix")
  features <- o6_read_tsv_gz(file.path(base, "features.tsv.gz"))
  barcodes <- o6_read_tsv_gz(file.path(base, "barcodes.tsv.gz"))[[1L]]
  triplet <- Matrix::readMM(file.path(base, "matrix.mtx.gz"))
  counts <- methods::as(triplet, "CsparseMatrix")
  if (!methods::is(counts, "dgCMatrix")) stop("Unexpected sparse matrix class for ", area)
  symbols <- toupper(as.character(features[[2L]]))
  if (anyDuplicated(symbols)) {
    counts <- Matrix::sparseMatrix(i = match(symbols[counts@i + 1L], unique(symbols)),
      j = rep.int(seq_len(ncol(counts)), diff(counts@p)), x = counts@x,
      dims = c(length(unique(symbols)), ncol(counts)),
      dimnames = list(unique(symbols), barcodes))
  } else {
    dimnames(counts) <- list(symbols, barcodes)
  }
  counts
}

o6_ucell <- function(counts, gene_sets, area) {
  object <- SeuratObject::CreateSeuratObject(counts = counts, assay = "Spatial", project = area,
    min.cells = 0, min.features = 0)
  object <- Seurat::NormalizeData(object, assay = "Spatial", normalization.method = "LogNormalize",
    scale.factor = 10000, verbose = FALSE)
  object <- UCell::AddModuleScore_UCell(object, features = gene_sets, assay = "Spatial", slot = "data",
    BPPARAM = BiocParallel::SerialParam())
  meta <- object[[]]
  cols <- paste0(names(gene_sets), "_UCell")
  out <- data.frame(barcode = rownames(meta), meta[, cols, drop = FALSE], check.names = FALSE)
  names(out)[-1L] <- names(gene_sets)
  rm(object); gc(verbose = FALSE)
  out
}

o6_edge_index <- function(barcodes, edges) {
  from <- match(edges$from, barcodes); to <- match(edges$to, barcodes)
  if (anyNA(from) || anyNA(to) || any(from >= to) || anyDuplicated(paste(from, to))) {
    stop("Frozen adjacency is incompatible with authority barcode order")
  }
  list(from = from, to = to, weight = as.numeric(edges$weight), s0 = 2 * sum(as.numeric(edges$weight)))
}

o6_moran <- function(x, edge) {
  z <- x - mean(x); den <- sum(z^2)
  if (!is.finite(den) || den == 0) return(NA_real_)
  length(x) / edge$s0 * (2 * sum(edge$weight * z[edge$from] * z[edge$to])) / den
}

o6_geary <- function(x, edge) {
  z <- x - mean(x); den <- sum(z^2)
  if (!is.finite(den) || den == 0) return(NA_real_)
  delta <- x[edge$from] - x[edge$to]
  (length(x) - 1) / (2 * edge$s0) * (2 * sum(edge$weight * delta^2)) / den
}

o6_permutation <- function(x, edge, endpoint, seed, permutations = 9999L) {
  fun <- if (endpoint == "CORE_MORAN") o6_moran else o6_geary
  set.seed(as.integer(seed))
  observed <- fun(x, edge)
  null <- replicate(permutations, fun(sample(x, replace = FALSE), edge))
  center <- mean(null)
  p <- (1 + sum(abs(null - center) >= abs(observed - center))) / (permutations + 1)
  list(observed = observed, p = p, n = length(null), null_mean = center,
    null_sd = stats::sd(null), null_min = min(null), null_max = max(null), all_finite = all(is.finite(null)))
}

o6_authority_paths <- function() {
  reg <- utils::read.csv(file.path(O6_REVISION,
    "06n1_protocol_review_material_bundle/06_spatial_data_structure/SPATIAL_OBJECT_REGISTRY.csv"),
    check.names = FALSE)
  list(
    registry = reg,
    moran_geary = file.path(O6_REVISION,
      "06n1_protocol_review_material_bundle/04_wp3_spatial_authorities/WP3_FINAL_MORAN_GEARY_AUTHORITY.csv")
  )
}

o6_area_roles <- function() {
  x <- utils::read.csv(file.path(O6_REVISION,
    "06b_wp3b_spatial_scope_method_resolution/WP3_CAPTURE_AREA_ROLE_FREEZE.csv"), check.names = FALSE)
  stopifnot(nrow(x) == 9L, sum(x$analysis_role == "PRIMARY_DLBCL") == 5L,
    sum(x$analysis_role == "CONTEXT_ONLY") == 4L)
  x
}
