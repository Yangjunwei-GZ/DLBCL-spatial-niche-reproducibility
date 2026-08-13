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

wp3_node_36g <- function(state) {
  state$edges <- wp3_visium_edges(state$coordinates, state$capture_area_id)
  state$edge_index <- wp3_edge_index(state$barcodes, state$edges)
  degree <- tabulate(c(state$edge_index$from_index, state$edge_index$to_index),
                     nbins = length(state$barcodes))
  if (any(state$edge_index$from_index >= state$edge_index$to_index) || max(degree) > 6) {
    stop("Invalid first-order Visium graph")
  }
  state
}

if (sys.nframe() == 0L) wp3_node_cli("36g", wp3_node_36g)
