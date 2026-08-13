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

wp3_node_36e <- function(state) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) stop("SeuratObject unavailable")
  sets <- wp3_gene_sets_from_eligibility(state$eligibility, include_exploratory = TRUE)
  normalized <- as.matrix(SeuratObject::LayerData(state$primary_object, assay = "Spatial", layer = "data"))
  state$ssgsea_scores <- wp3_ssgsea(normalized, sets)
  state$primary_object <- NULL
  state
}

if (sys.nframe() == 0L) wp3_node_cli("36e", wp3_node_36e)
