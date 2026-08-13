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

wp3_node_36b <- function(state) {
  areas <- wp3_read_area_contract()
  role <- areas$analysis_role[match(state$capture_area_id, areas$capture_area_id)]
  if (is.na(role)) stop("Capture area absent from frozen role contract")
  state$eligibility <- wp3_compute_eligibility(rownames(state$counts), state$capture_area_id, role)
  state$area_role <- role
  state
}

if (sys.nframe() == 0L) wp3_node_cli("36b", wp3_node_36b)
