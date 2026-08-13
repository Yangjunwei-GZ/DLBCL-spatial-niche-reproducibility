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

wp3_node_36j <- function(state) {
  eligible <- state$eligibility$eligibility %in% c("PRIMARY_ELIGIBLE", "CONTEXT_PRIMARY_COVERAGE")
  state$pc_status <- if (all(eligible)) "AUTHORIZED_EXPLORATORY" else "PROHIBITED_INCOMPLETE_PRIMARY_COVERAGE"
  state$pc <- if (all(eligible)) wp3_project_pc(state$primary_scores, state$eligibility, state$area_role) else NULL
  state
}

if (sys.nframe() == 0L) wp3_node_cli("36j", wp3_node_36j)
