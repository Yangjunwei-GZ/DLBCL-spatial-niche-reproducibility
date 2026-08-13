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

WP3_ROOT <- DLBCL_PROJECT_ROOT
WP3_CONTRACT_DIR <- file.path(WP3_ROOT, "revision_2026_reviewer_response/06b_wp3b_spatial_scope_method_resolution")
WP3_FAILED_OUTPUT <- file.path(WP3_ROOT, "revision_2026_reviewer_response/06c_wp3_real_spatial_continuous_analysis")
WP3_INTERRUPTED_OUTPUT <- file.path(WP3_FAILED_OUTPUT, "continuation_v2")
WP3_FUTURE_OUTPUT <- file.path(WP3_FAILED_OUTPUT, "continuation_v3")
WP3_CONTINUATION_CONTRACT_DIR <- file.path(WP3_ROOT, "revision_2026_reviewer_response/06h_wp3_interrupted_continuation_sparse_resume")
WP3_TOKEN_NAME <- "DLBCL_REVISION_ALLOW_WP3_SPATIAL"
WP3_TOKEN_VALUE <- "EXPLICITLY_APPROVED_WP3_CANONICAL_CONTINUOUS_SPATIAL_ANALYSIS"
WP3_SEED <- 20260730L
WP3_FORBIDDEN_FIELDS <- c("spot_class", "highest_program_label", "patient_assignment",
  "specimen_assignment", "region_class", "taxonomy", "selected_k", "cluster")

wp3_require_r <- function() {
  if (paste(R.version$major, R.version$minor, sep = ".") != "4.5.1") {
    stop("WP3 requires R 4.5.1")
  }
  frozen <- file.path(WP3_ROOT,
    "revision_2026_reviewer_response/04b_stage4_environment_freeze/stage4_renv_project/renv/library/windows/R-4.5/x86_64-w64-mingw32")
  if (!dir.exists(frozen)) stop("Frozen R library is missing: ", frozen)
  isolated <- file.path(WP3_CONTRACT_DIR, ".wp3_r_library")
  .libPaths(c(isolated, frozen, .Library))
  invisible(frozen)
}

wp3_require_runtime_packages <- function() {
  required <- c("Matrix", "SeuratObject", "Seurat", "sctransform", "future", "glmGamPoi", "BiocParallel", "UCell", "GSVA")
  available <- vapply(required, requireNamespace, logical(1), quietly = TRUE)
  if (!all(available)) stop("WP3 runtime packages unavailable before expression read: ",
                            paste(required[!available], collapse = ","))
  resolution_path <- file.path(WP3_CONTRACT_DIR, "WP3_RUNTIME_PACKAGE_RESOLUTION.csv")
  if (!file.exists(resolution_path)) stop("WP3 package resolution contract missing")
  resolution <- read.csv(resolution_path, check.names = FALSE)
  if (!setequal(resolution$package, required) || !all(resolution$load_status == "PASS")) {
    stop("WP3 package resolution contract is incomplete")
  }
  invisible(resolution)
}

wp3_require_token <- function(node) {
  if (!node %in% c("36a", "36b", "36c", "36d", "36e", "36f", "36g", "36h", "36i", "36j", "36k", "36v", "36")) {
    stop("Unregistered WP3 node: ", node)
  }
  if (!identical(Sys.getenv(WP3_TOKEN_NAME, unset = ""), WP3_TOKEN_VALUE)) {
    stop("WP3 execution token absent or invalid")
  }
  invisible(TRUE)
}

wp3_assert_output_root <- function(path) {
  requested <- normalizePath(path, winslash = "/", mustWork = FALSE)
  allowed <- normalizePath(WP3_FUTURE_OUTPUT, winslash = "/", mustWork = FALSE)
  if (!identical(tolower(requested), tolower(allowed))) stop("Unauthorized WP3 output root: ", requested)
  invisible(requested)
}

wp3_mkdir_once <- function(path) {
  if (file.exists(path) || dir.exists(path)) stop("create-once refusal: ", path)
  if (!dir.create(path, recursive = TRUE, showWarnings = FALSE)) stop("Cannot create directory: ", path)
  invisible(path)
}

wp3_write_csv_once <- function(x, path) {
  if (file.exists(path)) stop("create-once refusal: ", path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  invisible(path)
}

wp3_save_rds_once <- function(x, path) {
  if (file.exists(path)) stop("create-once refusal: ", path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, path, compress = "xz")
  invisible(path)
}

wp3_read_program_contract <- function() {
  path <- file.path(WP3_ROOT,
    "revision_2026_reviewer_response/05x_wp1_continuous_score_freeze/WP1_CANONICAL_PROGRAM_CONTRACT.csv")
  x <- read.csv(path, check.names = FALSE)
  stopifnot(nrow(x) == 6L, all(x$canonical_gene_count == 22L), all(x$unique_gene_count == 22L))
  x$genes <- strsplit(x$canonical_gene_members, ";", fixed = TRUE)
  x
}

wp3_read_area_contract <- function() {
  x <- read.csv(file.path(WP3_CONTRACT_DIR, "WP3_CAPTURE_AREA_ROLE_FREEZE.csv"), check.names = FALSE)
  stopifnot(nrow(x) == 9L, sum(x$analysis_role == "PRIMARY_DLBCL") == 5L,
            sum(x$analysis_role == "CONTEXT_ONLY") == 4L)
  x
}

wp3_read_tsv_gz <- function(path, header = FALSE) {
  read.delim(gzfile(path), header = header, quote = "", comment.char = "", check.names = FALSE)
}

wp3_read_coordinates <- function(path) {
  first <- readLines(path, n = 1L, warn = FALSE)
  has_header <- grepl("barcode", first, ignore.case = TRUE)
  x <- read.csv(path, header = has_header, check.names = FALSE)
  if (!has_header) names(x)[seq_len(min(6L, ncol(x)))] <-
    c("barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col")[seq_len(min(6L, ncol(x)))]
  needed <- c("barcode", "array_row", "array_col")
  if (!all(needed %in% names(x))) stop("Coordinates lack required columns: ", path)
  x
}

wp3_area_input_paths <- function(capture_area_id) {
  base <- file.path(WP3_ROOT, "00_raw_data/GSE276542/standard_10x", capture_area_id)
  list(
    base = base,
    matrix = file.path(base, "filtered_feature_bc_matrix/matrix.mtx.gz"),
    features = file.path(base, "filtered_feature_bc_matrix/features.tsv.gz"),
    barcodes = file.path(base, "filtered_feature_bc_matrix/barcodes.tsv.gz"),
    coordinates = file.path(base, "spatial/tissue_positions_list.csv")
  )
}

wp3_validate_area_structure <- function(capture_area_id) {
  p <- wp3_area_input_paths(capture_area_id)
  required <- unlist(p[c("matrix", "features", "barcodes", "coordinates")])
  if (!all(file.exists(required))) stop("Missing real spatial input for ", capture_area_id)
  features <- wp3_read_tsv_gz(p$features)
  barcodes <- wp3_read_tsv_gz(p$barcodes)[[1L]]
  coords <- wp3_read_coordinates(p$coordinates)
  if (anyDuplicated(barcodes) || anyDuplicated(coords$barcode)) stop("Duplicate barcodes: ", capture_area_id)
  if (!all(barcodes %in% coords$barcode)) stop("Incomplete barcode-coordinate join: ", capture_area_id)
  list(paths = p, features = features, barcodes = barcodes,
       coordinates = coords[match(barcodes, coords$barcode), , drop = FALSE])
}

wp3_read_real_counts <- function(structure) {
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix package unavailable")
  triplet_counts <- Matrix::readMM(structure$paths$matrix)
  counts <- methods::as(triplet_counts, "CsparseMatrix")
  if (!methods::is(counts, "dgCMatrix")) {
    stop("Unexpected compressed sparse matrix class after Matrix::readMM: ",
         paste(class(counts), collapse = ";"))
  }
  if (nrow(counts) != nrow(structure$features) || ncol(counts) != length(structure$barcodes)) {
    stop("Matrix dimensions do not match features/barcodes")
  }
  symbols <- as.character(structure$features[[2L]])
  if (anyDuplicated(symbols)) {
    counts <- Matrix::sparseMatrix(
      i = match(symbols[counts@i + 1L], unique(symbols)),
      j = rep.int(seq_len(ncol(counts)), diff(counts@p)), x = counts@x,
      dims = c(length(unique(symbols)), ncol(counts)),
      dimnames = list(unique(symbols), structure$barcodes))
  } else {
    dimnames(counts) <- list(symbols, structure$barcodes)
  }
  counts
}

wp3_eligibility_label <- function(n) {
  ifelse(n >= 16L, "PRIMARY_ELIGIBLE", ifelse(n >= 12L, "EXPLORATORY_ONLY", "INELIGIBLE"))
}

wp3_compute_eligibility <- function(feature_symbols, capture_area_id, area_role) {
  pc <- wp3_read_program_contract()
  rows <- lapply(seq_len(nrow(pc)), function(i) {
    genes <- pc$genes[[i]]
    detected <- genes[genes %in% feature_symbols]
    label <- wp3_eligibility_label(length(detected))
    if (area_role == "CONTEXT_ONLY" && label == "PRIMARY_ELIGIBLE") label <- "CONTEXT_PRIMARY_COVERAGE"
    data.frame(capture_area_id = capture_area_id, area_role = area_role,
      program_order = pc$program_order[i], program_id = pc$program_id[i], program_name = pc$program_name[i],
      canonical_gene_count = 22L, detected_gene_count = length(detected),
      coverage_fraction = length(detected) / 22, detected_genes = paste(detected, collapse = ";"),
      missing_genes = paste(setdiff(genes, detected), collapse = ";"), eligibility = label,
      substituted_genes = "", stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

wp3_gene_sets_from_eligibility <- function(eligibility, include_exploratory = TRUE) {
  keep <- eligibility$eligibility %in% c("PRIMARY_ELIGIBLE", "CONTEXT_PRIMARY_COVERAGE")
  if (include_exploratory) keep <- keep | eligibility$eligibility == "EXPLORATORY_ONLY"
  out <- strsplit(eligibility$detected_genes[keep], ";", fixed = TRUE)
  names(out) <- eligibility$program_id[keep]
  out
}

wp3_make_seurat <- function(counts, capture_area_id) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) stop("SeuratObject unavailable")
  SeuratObject::CreateSeuratObject(counts = counts, assay = "Spatial", project = capture_area_id,
                                   min.cells = 0, min.features = 0)
}

wp3_extract_ucell <- function(object, program_ids) {
  meta <- object[[]]
  cols <- paste0(program_ids, "_UCell")
  if (!all(cols %in% names(meta))) stop("UCell columns missing")
  data.frame(barcode = rownames(meta), meta[, cols, drop = FALSE], check.names = FALSE)
}

wp3_primary_ucell <- function(counts, gene_sets, capture_area_id) {
  if (!requireNamespace("Seurat", quietly = TRUE) || !requireNamespace("UCell", quietly = TRUE)) {
    stop("Seurat and UCell are required")
  }
  object <- wp3_make_seurat(counts, capture_area_id)
  object <- Seurat::NormalizeData(object, assay = "Spatial", normalization.method = "LogNormalize",
                                  scale.factor = 10000, verbose = FALSE)
  object <- UCell::AddModuleScore_UCell(object, features = gene_sets, assay = "Spatial", slot = "data",
                                        BPPARAM = BiocParallel::SerialParam())
  list(object = object, scores = wp3_extract_ucell(object, names(gene_sets)))
}

wp3_sct_ucell <- function(counts, gene_sets, capture_area_id) {
  if (!requireNamespace("Seurat", quietly = TRUE) || !requireNamespace("UCell", quietly = TRUE)) {
    stop("Seurat and UCell are required")
  }
  if (!requireNamespace("future", quietly = TRUE)) stop("future is required for frozen sequential SCTransform")
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::sequential)
  object <- wp3_make_seurat(counts, capture_area_id)
  object <- Seurat::SCTransform(object, assay = "Spatial", new.assay.name = "SCT",
    vst.flavor = "v2", method = "glmGamPoi", vars.to.regress = NULL,
    return.only.var.genes = FALSE, verbose = FALSE)
  object <- UCell::AddModuleScore_UCell(object, features = gene_sets, assay = "SCT", slot = "data",
                                        BPPARAM = BiocParallel::SerialParam())
  list(object = object, scores = wp3_extract_ucell(object, names(gene_sets)))
}

wp3_ssgsea <- function(normalized_matrix, gene_sets) {
  if (!requireNamespace("GSVA", quietly = TRUE)) stop("GSVA unavailable")
  param <- GSVA::ssgseaParam(exprData = normalized_matrix, geneSets = gene_sets,
    minSize = 5L, maxSize = 500L, alpha = 0.25, normalize = TRUE)
  scores <- GSVA::gsva(param, verbose = FALSE)
  data.frame(barcode = colnames(scores), t(as.matrix(scores)), check.names = FALSE)
}

wp3_score_concordance <- function(primary, sensitivity, capture_area_id, comparison) {
  common <- intersect(primary$barcode, sensitivity$barcode)
  p <- primary[match(common, primary$barcode), , drop = FALSE]
  s <- sensitivity[match(common, sensitivity$barcode), , drop = FALSE]
  programs <- intersect(setdiff(names(p), "barcode"), setdiff(names(s), "barcode"))
  do.call(rbind, lapply(programs, function(program) data.frame(
    capture_area_id = capture_area_id, comparison = comparison, program_id = program,
    n_spots = length(common), spearman_rho = cor(p[[program]], s[[program]], method = "spearman"),
    pearson_r = cor(p[[program]], s[[program]], method = "pearson"), stringsAsFactors = FALSE)))
}

wp3_visium_edges <- function(coords, capture_area_id = "area") {
  needed <- c("barcode", "array_row", "array_col")
  if (!all(needed %in% names(coords))) stop("Synthetic/real coordinates lack required columns")
  if (anyDuplicated(coords$barcode)) stop("Duplicate coordinate barcode")
  key <- paste(coords$array_row, coords$array_col, sep = ":")
  lookup <- setNames(seq_len(nrow(coords)), key)
  offsets <- rbind(c(0L, 2L), c(0L, -2L), c(1L, 1L), c(1L, -1L), c(-1L, 1L), c(-1L, -1L))
  edges <- vector("list", nrow(coords))
  for (i in seq_len(nrow(coords))) {
    target <- paste(coords$array_row[i] + offsets[, 1L], coords$array_col[i] + offsets[, 2L], sep = ":")
    j <- unname(lookup[target]); j <- j[!is.na(j) & j > i]
    if (length(j)) edges[[i]] <- data.frame(capture_area_id = capture_area_id,
      from = coords$barcode[i], to = coords$barcode[j], weight = 1, stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, edges)
  if (is.null(out)) out <- data.frame(capture_area_id = character(), from = character(), to = character(), weight = numeric())
  out
}

wp3_edge_index <- function(barcodes, edges) {
  needed <- c("from", "to", "weight")
  if (!all(needed %in% names(edges))) stop("Edge list lacks required columns")
  from_index <- match(edges$from, barcodes)
  to_index <- match(edges$to, barcodes)
  if (anyNA(from_index) || anyNA(to_index)) stop("Edge barcode is absent from score order")
  if (any(from_index >= to_index)) stop("Edges must be unique, forward-oriented, and non-self")
  if (any(!is.finite(edges$weight)) || any(edges$weight <= 0)) stop("Edge weights must be finite and positive")
  key <- paste(from_index, to_index, sep = ":")
  if (anyDuplicated(key)) stop("Duplicate undirected edge")
  data.frame(from_index = from_index, to_index = to_index,
             weight = as.numeric(edges$weight), stringsAsFactors = FALSE)
}

wp3_edge_s0 <- function(edge_index) 2 * sum(edge_index$weight)

wp3_moran_edge <- function(x, edge_index) {
  z <- x - mean(x); denominator <- sum(z^2); s0 <- wp3_edge_s0(edge_index)
  if (!length(x) || s0 == 0 || denominator == 0) return(NA_real_)
  numerator <- 2 * sum(edge_index$weight * z[edge_index$from_index] * z[edge_index$to_index])
  length(x) / s0 * numerator / denominator
}

wp3_geary_edge <- function(x, edge_index) {
  z <- x - mean(x); denominator <- sum(z^2); s0 <- wp3_edge_s0(edge_index)
  if (!length(x) || s0 == 0 || denominator == 0) return(NA_real_)
  delta <- x[edge_index$from_index] - x[edge_index$to_index]
  numerator <- 2 * sum(edge_index$weight * delta^2)
  (length(x) - 1) / (2 * s0) * numerator / denominator
}

wp3_bivariate_moran_edge <- function(x, y, edge_index) {
  zx <- as.numeric(scale(x)); zy <- as.numeric(scale(y)); s0 <- wp3_edge_s0(edge_index)
  if (!length(x) || s0 == 0 || any(!is.finite(c(zx, zy)))) return(NA_real_)
  numerator <- sum(edge_index$weight *
    (zx[edge_index$from_index] * zy[edge_index$to_index] +
     zx[edge_index$to_index] * zy[edge_index$from_index]))
  length(x) / s0 * numerator / sum(zx^2)
}

wp3_permutation_test_edge <- function(x, edge_index, statistic = c("moran", "geary"), permutations = 9999L, seed = WP3_SEED) {
  statistic <- match.arg(statistic)
  fun <- if (statistic == "moran") wp3_moran_edge else wp3_geary_edge
  set.seed(seed); observed <- fun(x, edge_index)
  null <- replicate(permutations, fun(sample(x, replace = FALSE), edge_index))
  p <- (1 + sum(abs(null - mean(null)) >= abs(observed - mean(null)))) / (permutations + 1)
  list(observed = observed, permutation_p = p, permutations = permutations)
}

wp3_bivariate_permutation_edge <- function(x, y, edge_index, permutations = 9999L, seed = WP3_SEED) {
  set.seed(seed); observed <- wp3_bivariate_moran_edge(x, y, edge_index)
  null <- replicate(permutations, wp3_bivariate_moran_edge(x, sample(y, replace = FALSE), edge_index))
  p <- (1 + sum(abs(null - mean(null)) >= abs(observed - mean(null)))) / (permutations + 1)
  list(observed = observed, permutation_p = p, permutations = permutations)
}

wp3_assign_bh_family <- function(area_role, eligibility, endpoint, capture_area_id = NA_character_, contains_antigen = FALSE) {
  if (endpoint == "BIVARIATE") {
    if (area_role == "CONTEXT_ONLY") return(paste0(endpoint, "_CONTEXT_", capture_area_id))
    if (eligibility == "EXPLORATORY_ONLY" || contains_antigen) return(paste0(endpoint, "_EXPLORATORY_ANTIGEN_", capture_area_id))
    return(paste0(endpoint, "_PRIMARY_", capture_area_id))
  }
  if (area_role == "CONTEXT_ONLY") return(paste0(endpoint, "_CONTEXT"))
  if (eligibility == "EXPLORATORY_ONLY" || contains_antigen) return(paste0(endpoint, "_EXPLORATORY_ANTIGEN"))
  paste0(endpoint, "_PRIMARY_DLBCL")
}

wp3_adjust_families <- function(x, p_col = "permutation_p", family_col = "bh_family") {
  x$fdr <- NA_real_
  for (family in unique(x[[family_col]])) {
    idx <- which(x[[family_col]] == family & is.finite(x[[p_col]]))
    x$fdr[idx] <- p.adjust(x[[p_col]][idx], method = "BH")
  }
  x
}

wp3_project_pc <- function(scores, eligibility, area_role) {
  allowed <- eligibility$eligibility %in% c("PRIMARY_ELIGIBLE", "CONTEXT_PRIMARY_COVERAGE")
  if (sum(allowed) != 6L || nrow(eligibility) != 6L) stop("PC prohibited: all six programs must be primary-eligible")
  params <- read.csv(file.path(WP3_ROOT,
    "revision_2026_reviewer_response/05x_wp1_continuous_score_freeze/WP1_CONTINUOUS_PROJECTION_PARAMETERS.csv"),
    check.names = FALSE)
  ids <- wp3_read_program_contract()$program_id
  if (!all(ids %in% names(scores))) stop("PC prohibited: missing program score")
  z <- scale(as.matrix(scores[, ids, drop = FALSE]))
  if (any(!is.finite(z))) stop("PC prohibited: non-finite within-capture z score")
  data.frame(barcode = scores$barcode,
    PC1 = as.vector(z %*% params$PC1_loading),
    PC2 = as.vector(z %*% params$PC2_loading),
    pc_role = if (area_role == "CONTEXT_ONLY") "CONTEXT_ONLY_EXPLORATORY" else "DLBCL_EXPLORATORY",
    stringsAsFactors = FALSE)
}

wp3_forbidden_field_check <- function(x) {
  bad <- intersect(tolower(names(x)), WP3_FORBIDDEN_FIELDS)
  if (length(bad)) stop("Forbidden output fields: ", paste(bad, collapse = ","))
  TRUE
}

wp3_sha256_windows <- function(path) {
  old_target <- Sys.getenv("WP3_SHA256_TARGET", unset = NA_character_)
  on.exit(if (is.na(old_target)) Sys.unsetenv("WP3_SHA256_TARGET") else
    Sys.setenv(WP3_SHA256_TARGET = old_target), add = TRUE)
  Sys.setenv(WP3_SHA256_TARGET = normalizePath(path, winslash = "\\", mustWork = TRUE))
  command <- "(Get-FileHash -LiteralPath $env:WP3_SHA256_TARGET -Algorithm SHA256).Hash.ToLowerInvariant()"
  output <- system2("powershell.exe", c("-NoLogo", "-NoProfile", "-NonInteractive",
    "-Command", shQuote(command)), stdout = TRUE, stderr = TRUE)
  candidates <- tolower(trimws(output[grepl("^[0-9A-Fa-f]{64}$", trimws(output))]))
  if (length(candidates) != 1L) stop("Cannot obtain unique SHA-256 for: ", path)
  candidates[[1L]]
}

wp3_validate_interrupted_scene <- function(
    failed_root = WP3_FAILED_OUTPUT,
    manifest_path = file.path(WP3_CONTINUATION_CONTRACT_DIR,
                              "WP3_INTERRUPTED_CONTINUATION_V2_FROZEN_MANIFEST.csv")) {
  if (!dir.exists(failed_root)) stop("Frozen interrupted output root is missing")
  if (!file.exists(manifest_path)) stop("Frozen interrupted manifest is missing")
  manifest <- read.csv(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("relative_path", "file_role", "sha256", "size", "mtime",
                "complete_output", "partial_output", "reuse_authorized", "overwrite_allowed")
  if (!all(required %in% names(manifest)) || !nrow(manifest)) stop("Frozen interrupted manifest is invalid")
  file_rows <- manifest$file_role != "EMPTY_DIRECTORY"
  observed <- list.files(failed_root, recursive = TRUE, full.names = FALSE,
                         all.files = TRUE, no.. = TRUE)
  observed <- sort(gsub("\\\\", "/", observed[file.info(file.path(failed_root, observed))$isdir %in% FALSE]))
  observed <- observed[!startsWith(observed, "continuation_v3/")]
  if (!identical(observed, sort(manifest$relative_path[file_rows]))) stop("Frozen interrupted file set changed")
  for (i in which(file_rows)) {
    path <- file.path(failed_root, manifest$relative_path[[i]])
    info <- file.info(path)
    expected_mtime <- sub("([+-][0-9]{2}):([0-9]{2})$", "\\1\\2", manifest$mtime[[i]])
    expected_mtime <- as.POSIXct(expected_mtime, format = "%Y-%m-%dT%H:%M:%OS%z")
    mtime_matches <- !is.na(expected_mtime) &&
      abs(as.numeric(info$mtime) - as.numeric(expected_mtime)) <= 0.001
    if (!identical(unname(info$size), as.numeric(manifest$size[[i]])) ||
        !identical(wp3_sha256_windows(path), tolower(manifest$sha256[[i]])) ||
        !mtime_matches) {
      stop("Frozen interrupted hash, size, or mtime changed: ", manifest$relative_path[[i]])
    }
  }
  if (any(manifest$overwrite_allowed != "FALSE")) stop("Frozen interrupted overwrite permission is invalid")
  invisible(manifest)
}

wp3_validate_complete_v2_caparea1 <- function(root = WP3_INTERRUPTED_OUTPUT) {
  expected <- c(
    "01_input_qc" = 1L, "02_eligibility" = 6L, "03_primary_ucell" = 3384L,
    "04_sct_ucell" = 3384L, "05_ssgsea" = 3384L, "06_score_qc" = 18L,
    "07_method_concordance" = 12L, "08_adjacency" = 9866L, "09_moran" = 6L,
    "10_geary" = 6L, "11_bivariate_moran" = 15L, "12_exploratory_pc" = 3384L)
  area_files <- list.files(root, pattern = "Cap[.]area1", recursive = TRUE, full.names = TRUE)
  if (length(area_files) != 12L || any(file.info(area_files)$size <= 0)) {
    stop("Cap.area1 complete file set is missing, extra, or zero-byte")
  }
  if (length(list.files(root, pattern = "[.](tmp|partial|lock)$", recursive = TRUE))) {
    stop("Interrupted V2 contains a temporary area file")
  }
  loaded <- list()
  for (directory in names(expected)) {
    path <- list.files(file.path(root, directory), pattern = "Cap[.]area1.*[.]csv$", full.names = TRUE)
    if (length(path) != 1L) stop("Cap.area1 output missing or duplicated: ", directory)
    value <- tryCatch(read.csv(path, check.names = FALSE, stringsAsFactors = FALSE),
      error = function(e) stop("Cap.area1 CSV is truncated or unreadable: ", directory))
    if (nrow(value) != expected[[directory]]) stop("Cap.area1 row count mismatch: ", directory)
    wp3_forbidden_field_check(value); loaded[[directory]] <- value
  }
  score_names <- c("03_primary_ucell", "04_sct_ucell", "05_ssgsea")
  barcodes <- loaded[[score_names[[1L]]]]$barcode
  if (anyDuplicated(barcodes) || !all(vapply(loaded[score_names], function(x)
      identical(x$barcode, barcodes) && all(is.finite(as.matrix(x[, -1L, drop = FALSE]))), logical(1)))) {
    stop("Cap.area1 score barcode or finite-value validation failed")
  }
  invisible(loaded)
}

wp3_node_cli <- function(node, fun) {
  wp3_require_r(); wp3_require_token(node)
  on.exit(Sys.unsetenv(WP3_TOKEN_NAME), add = TRUE)
  wp3_require_runtime_packages()
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args)) stop("Node requires an input state RDS path")
  state <- readRDS(args[[1L]])
  result <- fun(state)
  output <- if (length(args) >= 2L) args[[2L]] else sub("[.]rds$", paste0("_", node, ".rds"), args[[1L]])
  wp3_save_rds_once(result, output)
  invisible(result)
}
