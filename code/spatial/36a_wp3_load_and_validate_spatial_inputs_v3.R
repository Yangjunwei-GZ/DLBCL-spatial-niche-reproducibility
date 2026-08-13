DLBCL_PROJECT_ROOT <- normalizePath(
  Sys.getenv("DLBCL_PROJECT_ROOT", unset = getwd()),
  winslash = "/", mustWork = FALSE
)
DLBCL_SUPPLEMENTARY_CODE_ROOT <- normalizePath(
  Sys.getenv("DLBCL_SUPPLEMENTARY_CODE_ROOT",
    unset = file.path(DLBCL_PROJECT_ROOT, "Supplementary_Code")),
  winslash = "/", mustWork = FALSE
)

source(file.path(DLBCL_SUPPLEMENTARY_CODE_ROOT, "spatial", "wp3_common_v3.R"))

wp3_node_36a <- function(state) {
  stopifnot(length(state$capture_area_id) == 1L)
  structure <- wp3_validate_area_structure(state$capture_area_id)
  counts <- wp3_read_real_counts(structure)
  input_qc <- data.frame(capture_area_id = state$capture_area_id,
    feature_count = nrow(counts), spot_count = ncol(counts), coordinate_rows = nrow(structure$coordinates),
    barcode_coordinate_join_complete = all(structure$barcodes == structure$coordinates$barcode),
    duplicate_barcodes = anyDuplicated(structure$barcodes), expression_values_read = TRUE,
    input_matrix = normalizePath(structure$paths$matrix, winslash = "/"), stringsAsFactors = FALSE)
  list(capture_area_id = state$capture_area_id, counts = counts,
       features = structure$features, barcodes = structure$barcodes,
       coordinates = structure$coordinates, input_qc = input_qc, expression_values_read = TRUE,
       matrix_class = class(counts)[[1L]], sparse_coercion_warning_count = 0L,
       other_warning_count = 0L, stack_imbalance = FALSE)
}

if (sys.nframe() == 0L) wp3_node_cli("36a", wp3_node_36a)
