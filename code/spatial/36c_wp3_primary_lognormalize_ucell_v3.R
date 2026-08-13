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

wp3_node_36c <- function(state) {
  sets <- wp3_gene_sets_from_eligibility(state$eligibility, include_exploratory = TRUE)
  scored <- wp3_primary_ucell(state$counts, sets, state$capture_area_id)
  names(scored$scores) <- sub("_UCell$", "", names(scored$scores))
  state$primary_object <- scored$object
  state$primary_scores <- scored$scores
  state
}

if (sys.nframe() == 0L) wp3_node_cli("36c", wp3_node_36c)
