DLBCL_PROJECT_ROOT <- normalizePath(
  Sys.getenv("DLBCL_PROJECT_ROOT", unset = getwd()),
  winslash = "/", mustWork = FALSE
)
DLBCL_SUPPLEMENTARY_CODE_ROOT <- normalizePath(
  Sys.getenv("DLBCL_SUPPLEMENTARY_CODE_ROOT",
    unset = file.path(DLBCL_PROJECT_ROOT, "Supplementary_Code")),
  winslash = "/", mustWork = FALSE
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Expected one attempt-directory argument.")
attempt_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
options(stringsAsFactors = FALSE, warn = 1)
root <- DLBCL_PROJECT_ROOT
base <- file.path(root, "revision_2026_reviewer_response/06t_figure5_singlecell_contextualization")
a008 <- file.path(base, "00_protocol_freeze/amendments/AMENDMENT_008_FIGURE5_PHASEB_CELLSET_AND_EXECUTION_FREEZE")
manifest_path <- file.path(a008, "FIGURE5_PRIMARY_DLBCL_CELL_MANIFEST.csv.gz")
  program_path <- file.path(root, "revision_2026_reviewer_response/05x_wp1_continuous_score_freeze/WP1_CANONICAL_PROGRAM_CONTRACT.csv")
  counts_path <- file.path(root, "00_raw_data/GSE182434/GSE182434_raw_count_matrix.txt.gz")
  annotation_path <- file.path(root, "00_raw_data/GSE182434/GSE182434_cell_annotation.txt.gz")
  lib_path <- file.path(root, "revision_2026_reviewer_response/06b_wp3b_spatial_scope_method_resolution/.wp3_r_library")
  ucell_description_path <- file.path(lib_path, "UCell/DESCRIPTION")
.libPaths(c(lib_path, .libPaths()))
warnings_seen <- character()
log_path <- file.path(attempt_dir, "FIGURE5_PHASEC_EXECUTION_LOG.txt")
log_line <- function(...) cat(paste0(...), "\n", file = log_path, append = TRUE)
write.csv.safe <- function(x, name) data.table::fwrite(x, file.path(attempt_dir, name), quote = TRUE)
write.csv.gz <- function(x, name) data.table::fwrite(x, file.path(attempt_dir, name), compress = "gzip", quote = TRUE)
finite_matrix <- function(x) all(is.finite(as.matrix(x)))

main <- function() {
  suppressPackageStartupMessages({
    library(data.table); library(Matrix); library(UCell); library(BiocParallel)
    library(Seurat); library(SeuratObject); library(ggplot2)
  })
  BiocParallel::register(BiocParallel::SerialParam())
  log_line("start_time: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  log_line("scientific_scope: A009_FIGURE5_PHASEC")
  log_line("filtering_run: FALSE")
  log_line("clustering_run: FALSE")
  log_line("statistical_testing_run: FALSE")
  log_line("taxonomy_run: FALSE")
  log_line("final_figure_rendered: FALSE")

  manifest <- fread(manifest_path)
  stopifnot(nrow(manifest) == 14368L, uniqueN(manifest$ID) == 14368L,
            all(manifest$primary_retained), sum(manifest$low_feature_sensitivity_flag) == 84L,
            uniqueN(manifest$Patient) == 4L, uniqueN(manifest$Sample) == 8L)
  programs <- fread(program_path)
  stopifnot(nrow(programs) == 6L, all(programs$canonical_gene_count == 22L))
  features <- setNames(strsplit(programs$canonical_gene_members, ";", fixed = TRUE), programs$program_id)
  stopifnot(sum(lengths(features)) == 132L, length(unique(unlist(features))) == 121L)

  expected <- fread(file.path(a008, "FIGURE5_PHASEB_INPUT_AUTHORITY.csv"))
  hash_map <- setNames(expected$sha256, expected$authority_role)
  observed_hash <- c(
    raw_counts = digest::digest(file = counts_path, algo = "sha256"),
    official_annotation = digest::digest(file = annotation_path, algo = "sha256"),
    primary_manifest = digest::digest(file = manifest_path, algo = "sha256"),
    canonical_program_contract = digest::digest(file = program_path, algo = "sha256"),
    UCell_DESCRIPTION = digest::digest(file = ucell_description_path, algo = "sha256")
  )
  stopifnot(observed_hash[["raw_counts"]] == hash_map[["raw_counts"]],
            observed_hash[["official_annotation"]] == hash_map[["official_annotation"]],
            observed_hash[["canonical_program_contract"]] == hash_map[["canonical_program_contract"]],
            observed_hash[["UCell_DESCRIPTION"]] == hash_map[["UCell_DESCRIPTION"]])
  primary_manifest_expected <- fread(file.path(a008, "FIGURE5_PRIMARY_DLBCL_CELL_MANIFEST_SHA256.csv"))$sha256[[1L]]
  stopifnot(observed_hash[["primary_manifest"]] == primary_manifest_expected)

  formals_txt <- capture.output(print(formals(UCell::ScoreSignatures_UCell)))
  writeLines(formals_txt, file.path(attempt_dir, "FIGURE5_UCELL_FORMALS.txt"))
  stopifnot(as.character(getRversion()) == "4.5.1", as.character(packageVersion("UCell")) == "2.14.0")
  stopifnot(identical(formals(UCell::ScoreSignatures_UCell)$maxRank, 1500))
  synthetic <- Matrix(matrix(rep(c(5, 0, 1, 0), 50), nrow = 20L, ncol = 10L), sparse = TRUE)
  rownames(synthetic) <- paste0("G", 1:20); colnames(synthetic) <- paste0("C", 1:10)
  synthetic_features <- list(synthetic_signature = paste0("G", 1:10))
  smoke1 <- UCell::ScoreSignatures_UCell(synthetic, synthetic_features, maxRank = 15,
    ties.method = "average", missing_genes = "skip", BPPARAM = SerialParam(), ncores = 1,
    name = "_raw_ucell")
  smoke2 <- UCell::ScoreSignatures_UCell(synthetic, synthetic_features, maxRank = 15,
    ties.method = "average", missing_genes = "skip", BPPARAM = SerialParam(), ncores = 1,
    name = "_raw_ucell")
  smoke_pass <- identical(smoke1, smoke2) && finite_matrix(smoke1)
  if (!smoke_pass) stop("UCELL_RUNTIME_COMPATIBILITY_FAILED")
  built <- packageDescription("UCell")$Built
  runtime_lines <- c("# UCell runtime compatibility QC", "",
    paste0("- R.version.string: ", R.version.string),
    paste0("- UCell version: ", packageVersion("UCell")),
    paste0("- UCell Built: ", built),
    "- Synthetic signature genes: 10", "- Synthetic maxRank: 15 (explicit)",
    paste0("- Synthetic finite output: ", finite_matrix(smoke1)),
    paste0("- Synthetic exact repeatability: ", identical(smoke1, smoke2)),
    "- Compatibility status: UCELL_BUILD_WARNING_NONBLOCKING")
  writeLines(runtime_lines, file.path(attempt_dir, "FIGURE5_UCELL_RUNTIME_COMPATIBILITY_QC.md"))

  con <- gzfile(counts_path, "rt"); on.exit(try(close(con), silent = TRUE), add = TRUE)
  header <- strsplit(readLines(con, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1L]]
  matrix_ids <- header[-1L]
  stopifnot(length(matrix_ids) == 28416L, !anyDuplicated(matrix_ids))
  selected <- match(manifest$ID, matrix_ids)
  stopifnot(!anyNA(selected), !anyDuplicated(selected))
  chunk_rows <- 16L; i_parts <- list(); j_parts <- list(); x_parts <- list()
  gene_names <- character(); row_offset <- 0L; chunk_index <- 0L
  repeat {
    lines <- readLines(con, n = chunk_rows, warn = FALSE)
    if (!length(lines)) break
    parts <- strsplit(lines, "\t", fixed = TRUE)
    if (any(lengths(parts) != length(matrix_ids) + 1L)) stop("Malformed matrix row.")
    genes <- vapply(parts, `[[`, character(1), 1L)
    vals <- lapply(parts, function(z) as.numeric(z[selected + 1L]))
    if (any(vapply(vals, function(z) anyNA(z) || any(!is.finite(z)) || any(z < 0), logical(1)))) stop("Invalid count value.")
    nz <- lapply(vals, function(z) which(z != 0))
    lens <- lengths(nz)
    chunk_index <- chunk_index + 1L
    i_parts[[chunk_index]] <- rep.int(row_offset + seq_along(parts), lens)
    j_parts[[chunk_index]] <- unlist(nz, use.names = FALSE)
    x_parts[[chunk_index]] <- unlist(Map(function(z, q) z[q], vals, nz), use.names = FALSE)
    gene_names <- c(gene_names, genes); row_offset <- row_offset + length(parts)
    if (row_offset %% 2048L < chunk_rows) log_line("sparse_rows_processed: ", row_offset)
  }
  close(con)
  counts <- sparseMatrix(i = unlist(i_parts, use.names = FALSE),
    j = unlist(j_parts, use.names = FALSE), x = unlist(x_parts, use.names = FALSE),
    dims = c(row_offset, nrow(manifest)), dimnames = list(gene_names, manifest$ID),
    giveCsparse = TRUE)
  counts <- as(counts, "dgCMatrix")
  rm(i_parts, j_parts, x_parts); gc()
  stopifnot(nrow(counts) == 49632L, !anyDuplicated(gene_names), ncol(counts) == 14368L,
            identical(colnames(counts), manifest$ID))
  gene_gate <- data.table(program_id = programs$program_id,
    matched = vapply(features, function(g) sum(g %in% rownames(counts)), integer(1)))
  stopifnot(all(gene_gate$matched == 22L))
  log_line("raw_matrix_dimensions: ", nrow(counts), "x", length(matrix_ids), " (selected columns: ", ncol(counts), ")")
  log_line("raw_matrix_cell_ids_unique: TRUE")
  log_line("raw_matrix_gene_ids_unique: TRUE")
  log_line("sparse_dimensions: ", nrow(counts), "x", ncol(counts))
  log_line("sparse_nnzero: ", nnzero(counts))
  log_line("sparse_object_size_bytes: ", as.numeric(object.size(counts)))
  log_line("full_dense_coercion: FALSE")

  preflight <- data.table(
    gate = c("raw_count_hash", "annotation_hash", "primary_manifest_hash", "canonical_program_hash",
      "UCell_DESCRIPTION_hash", "raw_matrix_dimensions", "raw_matrix_cell_ids_unique",
      "raw_matrix_gene_ids_unique", "primary_cells",
      "primary_membership", "primary_order", "patients", "nested_samples", "low_feature_retained",
      "programs", "memberships", "unique_genes", "gene_gate", "matrix_class", "dense_coercion",
      "R_version", "UCell_version", "smoke_test", "maxRank_explicit", "UMAP_workers"),
    expected = c(hash_map[["raw_counts"]], hash_map[["official_annotation"]], primary_manifest_expected,
      hash_map[["canonical_program_contract"]], hash_map[["UCell_DESCRIPTION"]], "49632x28416", "TRUE", "TRUE",
      "14368", "TRUE", "TRUE", "4", "8", "84", "6", "132", "121", "6x22",
      "dgCMatrix", "FALSE", "4.5.1", "2.14.0", "PASS", "1500", "1"),
    observed = c(observed_hash[["raw_counts"]], observed_hash[["official_annotation"]], observed_hash[["primary_manifest"]],
      observed_hash[["canonical_program_contract"]], observed_hash[["UCell_DESCRIPTION"]],
      paste0(nrow(counts), "x", length(matrix_ids)), !anyDuplicated(matrix_ids), !anyDuplicated(gene_names),
      nrow(manifest), setequal(manifest$ID, colnames(counts)), identical(manifest$ID, colnames(counts)),
      uniqueN(manifest$Patient), uniqueN(manifest$Sample), sum(manifest$low_feature_sensitivity_flag),
      nrow(programs), sum(lengths(features)), length(unique(unlist(features))), "6x22",
      class(counts)[1L], "FALSE", as.character(getRversion()), as.character(packageVersion("UCell")),
      ifelse(smoke_pass, "PASS", "FAIL"), "1500", BiocParallel::bpnworkers(BiocParallel::bpparam())))
  preflight[, status := ifelse(as.character(expected) == as.character(observed), "PASS", "FAIL")]
  fwrite(preflight, file.path(attempt_dir, "FIGURE5_PHASEC_PREFLIGHT_QC.csv"))
  stopifnot(all(preflight$status == "PASS"))

  u <- UCell::ScoreSignatures_UCell(matrix = counts, features = features,
    precalc.ranks = NULL, maxRank = 1500, w_neg = 1, name = "_raw_ucell",
    chunk.size = 100, missing_genes = "skip", BPPARAM = SerialParam(), ncores = 1,
    ties.method = "average", force.gc = FALSE)
  score_df <- as.data.frame(u, check.names = FALSE)
  stopifnot(identical(rownames(score_df), manifest$ID))
  score_cols <- paste0(programs$program_id, "_raw_ucell")
  stopifnot(identical(colnames(score_df), score_cols), finite_matrix(score_df))
  meta_cols <- c("ID", "Patient", "Sample", "CellType", "TumorNormal", "COO", "Tissue", "low_feature_sensitivity_flag")
  raw_out <- cbind(as.data.frame(manifest)[, meta_cols], score_df)
  write.csv.gz(raw_out, "FIGURE5_RAW_UCELL_SCORES.csv.gz")
  qfun <- function(v) quantile(v, c(0, .01, .05, .25, .5, .75, .95, .99, 1), type = 7, names = FALSE)
  score_qc <- rbindlist(lapply(seq_along(score_cols), function(k) {
    v <- score_df[[score_cols[k]]]; q <- qfun(v)
    data.table(program = programs$program_id[k], n_cells = length(v), n_genes_expected = 22L,
      n_genes_matched = gene_gate$matched[k], minimum = q[1], p01 = q[2], p05 = q[3],
      p25 = q[4], median = q[5], mean = mean(v), p75 = q[6], p95 = q[7], p99 = q[8],
      maximum = q[9], sd = sd(v), n_NA = sum(is.na(v)), n_Inf = sum(is.infinite(v)),
      n_unique_values = uniqueN(v))
  }))
  fwrite(score_qc, file.path(attempt_dir, "FIGURE5_UCELL_SCORE_QC.csv"))
  stopifnot(all(score_qc$n_cells == 14368L), all(score_qc$n_genes_matched == 22L),
            all(score_qc$n_NA == 0L), all(score_qc$n_Inf == 0L))

  display_cols <- paste0(programs$program_id, "_display_z")
  display <- as.data.frame(lapply(score_df, function(v) (v - mean(v)) / sd(v)), check.names = FALSE)
  colnames(display) <- display_cols
  display_qc <- data.table(program = programs$program_id,
    n_cells = nrow(display), mean_display_z = vapply(display, mean, numeric(1)),
    sd_display_z = vapply(display, sd, numeric(1)),
    abs_mean = abs(vapply(display, mean, numeric(1))),
    abs_sd_minus_1 = abs(vapply(display, sd, numeric(1)) - 1), tolerance = 1e-10)
  display_qc[, status := ifelse(abs_mean <= tolerance & abs_sd_minus_1 <= tolerance, "PASS", "FAIL")]
  fwrite(display_qc, file.path(attempt_dir, "FIGURE5_DISPLAY_Z_QC.csv"))
  stopifnot(all(display_qc$status == "PASS"))
  score_display <- cbind(raw_out, display)
  write.csv.gz(score_display, "FIGURE5_UCELL_SCORES_WITH_DISPLAY_Z.csv.gz")

  metadata <- as.data.frame(manifest)
  rownames(metadata) <- metadata$ID
  seu <- CreateSeuratObject(counts = counts, meta.data = metadata, assay = "RNA", min.cells = 0, min.features = 0)
  stopifnot(inherits(LayerData(seu, assay = "RNA", layer = "counts"), "sparseMatrix"), ncol(seu) == 14368L)
  seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  variable_features <- VariableFeatures(seu)
  stopifnot(length(variable_features) == 3000L)
  fwrite(data.table(feature_order = seq_along(variable_features), gene = variable_features),
    file.path(attempt_dir, "FIGURE5_VARIABLE_FEATURES.csv"))
  seu <- ScaleData(seu, features = variable_features, verbose = FALSE)
  set.seed(20260807)
  seu <- RunPCA(seu, features = variable_features, npcs = 30, seed.use = 20260807,
    verbose = FALSE)
  pca <- Embeddings(seu, "pca")[, 1:30, drop = FALSE]
  pca_reordered <- !identical(rownames(pca), manifest$ID)
  pca <- pca[match(manifest$ID, rownames(pca)), , drop = FALSE]
  stopifnot(identical(rownames(pca), manifest$ID), finite_matrix(pca))
  pca_out <- cbind(data.frame(ID = manifest$ID), as.data.frame(pca))
  write.csv.gz(pca_out, "FIGURE5_PCA_COORDINATES.csv.gz")
  sdev <- Stdev(seu, "pca")[1:30]
  pca_var <- data.table(PC = paste0("PC", 1:30), component = 1:30,
    standard_deviation = sdev, variance = sdev^2,
    proportion_of_computed_variance = sdev^2 / sum(sdev^2),
    cumulative_proportion = cumsum(sdev^2 / sum(sdev^2)), seed = 20260807)
  fwrite(pca_var, file.path(attempt_dir, "FIGURE5_PCA_VARIANCE_EXPLAINED.csv"))

  stopifnot(BiocParallel::bpnworkers(BiocParallel::bpparam()) == 1L)
  log_line("umap_n_threads: ", BiocParallel::bpnworkers(BiocParallel::bpparam()))
  set.seed(20260807)
  seu <- RunUMAP(seu, reduction = "pca", dims = 1:30, reduction.name = "umap",
    reduction.key = "UMAP_", metric = "cosine", n.neighbors = 30, min.dist = 0.3,
    spread = 1.0, seed.use = 20260807, verbose = FALSE)
  primary_umap <- Embeddings(seu, "umap")
  umap_reordered <- !identical(rownames(primary_umap), manifest$ID)
  primary_umap <- primary_umap[match(manifest$ID, rownames(primary_umap)), , drop = FALSE]
  set.seed(20260807)
  seu <- RunUMAP(seu, reduction = "pca", dims = 1:30, reduction.name = "umap_replicate",
    reduction.key = "UMAPR_", metric = "cosine", n.neighbors = 30, min.dist = 0.3,
    spread = 1.0, seed.use = 20260807, verbose = FALSE)
  replicate_umap <- Embeddings(seu, "umap_replicate")
  replicate_umap <- replicate_umap[match(manifest$ID, rownames(replicate_umap)), , drop = FALSE]
  stopifnot(finite_matrix(primary_umap), finite_matrix(replicate_umap))
  dif <- abs(primary_umap - replicate_umap)
  repro <- data.table(max_abs_difference_UMAP1 = max(dif[,1]),
    max_abs_difference_UMAP2 = max(dif[,2]), mean_abs_difference = mean(dif),
    coordinate_correlation = cor(as.vector(primary_umap), as.vector(replicate_umap)),
    seed = 20260807, primary_run = "first", replicate_use = "QC_ONLY")
  fwrite(repro, file.path(attempt_dir, "FIGURE5_UMAP_REPRODUCIBILITY_QC.csv"))
  umap_out <- data.table(ID = manifest$ID, Patient = manifest$Patient, Sample = manifest$Sample,
    CellType = manifest$CellType, Tissue = manifest$Tissue,
    low_feature_sensitivity_flag = manifest$low_feature_sensitivity_flag,
    UMAP_1 = primary_umap[,1], UMAP_2 = primary_umap[,2])
  write.csv.gz(umap_out, "FIGURE5_UMAP_COORDINATES.csv.gz")

  panel_a <- umap_out[, .(ID, Patient, Sample, CellType, UMAP_1, UMAP_2)]
  write.csv.gz(panel_a, "FIGURE5_PANEL_A_SOURCE_DATA.csv.gz")
  panel_a_counts <- panel_a[, .(n_cells=.N, n_patients=uniqueN(Patient), n_samples=uniqueN(Sample)), by=CellType]
  fwrite(panel_a_counts, file.path(attempt_dir, "FIGURE5_PANEL_A_CELL_COUNTS.csv"))
  panel_b <- rbindlist(lapply(seq_len(nrow(programs)), function(k) data.table(
    ID = manifest$ID, Patient = manifest$Patient, Sample = manifest$Sample, CellType = manifest$CellType,
    program = programs$program_id[k], program_name = programs$program_name[k],
    raw_ucell = score_df[[score_cols[k]]], display_z = display[[display_cols[k]]],
    UMAP_1 = primary_umap[,1], UMAP_2 = primary_umap[,2])))
  stopifnot(nrow(panel_b) == 86208L, all(panel_b[, .N, by=ID]$N == 6L))
  write.csv.gz(panel_b, "FIGURE5_PANEL_B_SOURCE_DATA.csv.gz")
  celltypes <- c("B cells", "Plasma cells", "T cells CD4", "T cells CD8", "Tregs", "TFH",
    "NK cells", "Monocytes and Macrophages", "Others")
  panel_c <- rbindlist(lapply(seq_len(nrow(programs)), function(k) {
    q75 <- quantile(score_df[[score_cols[k]]], probs=.75, type=7, na.rm=FALSE)
    rbindlist(lapply(seq_along(celltypes), function(j) {
      idx <- which(manifest$CellType == celltypes[j])
      data.table(program=programs$program_id[k], program_name=programs$program_name[k],
        program_order=programs$program_order[k], CellType=celltypes[j], celltype_order=j,
        n_cells=length(idx), mean_raw_ucell=mean(score_df[[score_cols[k]]][idx]),
        mean_display_z=mean(display[[display_cols[k]]][idx]), pooled_program_q75=as.numeric(q75),
        n_above_q75=sum(score_df[[score_cols[k]]][idx] > q75),
        fraction_above_q75=mean(score_df[[score_cols[k]]][idx] > q75))
    }))
  }))
  stopifnot(nrow(panel_c) == 54L)
  fwrite(panel_c, file.path(attempt_dir, "FIGURE5_PANEL_C_SOURCE_DATA.csv"))
  patient_levels <- unique(manifest$Patient)
  panel_d <- rbindlist(lapply(seq_len(nrow(programs)), function(k) {
    grid <- CJ(Patient=patient_levels, CellType=celltypes, unique=TRUE)
    v <- data.table(Patient=manifest$Patient, CellType=manifest$CellType, value=score_df[[score_cols[k]]])
    agg <- v[, .(n_cells=.N, mean_raw_ucell=mean(value), median_raw_ucell=median(value)), by=.(Patient,CellType)]
    z <- merge(grid, agg, by=c("Patient","CellType"), all.x=TRUE, sort=FALSE)
    z[is.na(n_cells), n_cells := 0L]
    z[, `:=`(program=programs$program_id[k], program_name=programs$program_name[k], program_order=programs$program_order[k])]
    z
  }), use.names=TRUE)
  setcolorder(panel_d, c("program","program_name","program_order","Patient","CellType","n_cells","mean_raw_ucell","median_raw_ucell"))
  stopifnot(nrow(panel_d) == 216L, all(is.na(panel_d[n_cells==0]$mean_raw_ucell)), all(is.na(panel_d[n_cells==0]$median_raw_ucell)))
  fwrite(panel_d, file.path(attempt_dir, "FIGURE5_PANEL_D_SOURCE_DATA.csv"), na="NA")
  coverage <- CJ(Patient=patient_levels, CellType=celltypes, unique=TRUE)
  coverage <- merge(coverage, data.table(Patient=manifest$Patient,CellType=manifest$CellType)[,.N,by=.(Patient,CellType)], by=c("Patient","CellType"), all.x=TRUE)
  setnames(coverage,"N","n_cells"); coverage[is.na(n_cells),n_cells:=0L]; coverage[,present:=n_cells>0]
  fwrite(coverage, file.path(attempt_dir, "FIGURE5_PANEL_D_COVERAGE_QC.csv"))

  low_qc <- data.table(metric=c("flagged_expected","flagged_scored","flagged_in_PCA_input","flagged_in_UMAP_output","flagged_missing"),
    expected=c(84,84,84,84,0), observed=c(84,sum(raw_out$low_feature_sensitivity_flag),
      sum(manifest$low_feature_sensitivity_flag),sum(umap_out$low_feature_sensitivity_flag),
      sum(!manifest$ID[manifest$low_feature_sensitivity_flag] %in% umap_out$ID)))
  low_qc[,status:=ifelse(expected==observed,"PASS","FAIL")]
  fwrite(low_qc,file.path(attempt_dir,"FIGURE5_LOW_FEATURE_PRIMARY_RETENTION_QC.csv"))
  id_sources <- list(primary_manifest=manifest$ID, raw_UCell=raw_out$ID,
    display_z=score_display$ID, Seurat_metadata=colnames(seu), PCA=pca_out$ID,
    UMAP=umap_out$ID, Panel_A=panel_a$ID, Panel_B=unique(panel_b$ID),
    Panel_C_aggregation_input=raw_out$ID, Panel_D_aggregation_input=raw_out$ID)
  id_qc <- rbindlist(lapply(names(id_sources), function(nm) data.table(source=nm,
    n_rows=length(id_sources[[nm]]), n_unique=uniqueN(id_sources[[nm]]),
    same_cell_membership=setequal(id_sources[[nm]],manifest$ID),
    same_ID_order=identical(id_sources[[nm]],manifest$ID),
    reordered_for_authority_output=if(nm=="PCA")pca_reordered else if(nm=="UMAP")umap_reordered else FALSE)))
  fwrite(id_qc,file.path(attempt_dir,"FIGURE5_CELL_IDENTITY_CONSISTENCY_QC.csv"))
  stopifnot(all(id_qc$same_cell_membership), all(id_qc$same_ID_order))

  score_palette <- colorRampPalette(c("#2166AC","#F7F7F7","#B2182B"))(101)
  score_col <- function(z){q <- pmax(-2.5,pmin(2.5,z));score_palette[1L+round((q+2.5)/5*100)]}
  ct <- factor(manifest$CellType, levels=celltypes); ctcols <- setNames(rainbow(length(celltypes),s=.65,v=.75),celltypes)
  draw_page <- function(){par(mfrow=c(2,4),mar=c(2,2,3,1),oma=c(1,1,3,1));plot(primary_umap,pch=16,cex=.18,col=ctcols[as.character(ct)],axes=FALSE,main="Official CellType");for(k in seq_len(6))plot(primary_umap,pch=16,cex=.18,col=score_col(display[[k]]),axes=FALSE,main=programs$program_name[k]);plot.new();mtext("SCIENTIFIC QC ONLY - NOT FOR SUBMISSION",outer=TRUE,cex=1.2,font=2)}
  pdf(file.path(attempt_dir,"FIGURE5_PHASEC_SCIENTIFIC_QC_ONLY.pdf"),width=12,height=8,useDingbats=FALSE);draw_page();par(mfrow=c(2,3),mar=c(4,4,3,1),oma=c(1,1,3,1));for(k in seq_len(6))hist(score_df[[k]],breaks=50,col="#6BAED6",border=NA,main=programs$program_name[k],xlab="raw UCell");mtext("SCIENTIFIC QC ONLY - NOT FOR SUBMISSION",outer=TRUE,cex=1.2,font=2);dev.off()
  png(file.path(attempt_dir,"FIGURE5_PHASEC_SCIENTIFIC_QC_ONLY_preview.png"),width=2400,height=1600,res=200);draw_page();dev.off()

  gates <- data.table(gate=c("primary_cells","no_exclusion","low_feature_retained","gene_gate","UCell_runtime","smoke_test","raw_scores","display_z","variable_features","PCA","UMAP","Panel_A","Panel_B","Panel_C","Panel_D","patients","samples","no_clustering","no_statistics","no_taxonomy","no_malignant_inference"),
    expected=c(14368,0,84,"6x22","PASS","PASS","PASS","PASS",3000,"14368x30",14368,14368,86208,54,216,4,8,FALSE,FALSE,FALSE,FALSE),
    observed=c(nrow(manifest),sum(!manifest$primary_retained),sum(manifest$low_feature_sensitivity_flag),"6x22","PASS",ifelse(smoke_pass,"PASS","FAIL"),ifelse(all(score_qc$n_NA==0 & score_qc$n_Inf==0),"PASS","FAIL"),ifelse(all(display_qc$status=="PASS"),"PASS","FAIL"),length(variable_features),paste(dim(pca),collapse="x"),nrow(umap_out),nrow(panel_a),nrow(panel_b),nrow(panel_c),nrow(panel_d),uniqueN(manifest$Patient),uniqueN(manifest$Sample),FALSE,FALSE,FALSE,FALSE))
  gates[,status:=ifelse(as.character(expected)==as.character(observed),"PASS","FAIL")]
  fwrite(gates,file.path(attempt_dir,"FIGURE5_PHASEC_COMPLETION_GATE.csv"))
  stopifnot(all(gates$status=="PASS"))
  writeLines(capture.output(sessionInfo()),file.path(attempt_dir,"sessionInfo.txt"))
  writeLines(c(runtime_lines,"","## Captured warnings",if(length(warnings_seen))paste0("- ",unique(warnings_seen)) else "- NONE"),file.path(attempt_dir,"FIGURE5_UCELL_RUNTIME_COMPATIBILITY_QC.md"))
  log_line("end_time: ",format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"))
  log_line("R_scientific_status: COMPLETE")
}

tryCatch(withCallingHandlers(main(), warning=function(w){warnings_seen <<- c(warnings_seen,conditionMessage(w));invokeRestart("muffleWarning")}),
  error=function(e){log_line("R_scientific_status: FAILED");log_line("error: ",conditionMessage(e));writeLines(capture.output(sessionInfo()),file.path(attempt_dir,"sessionInfo.txt"));quit(status=1L)})
