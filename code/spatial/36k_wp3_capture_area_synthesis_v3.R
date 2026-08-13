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

wp3_node_36k <- function(state) {
  if (is.null(state$area_results) || !length(state$area_results)) stop("No area results supplied")
  ac <- do.call(rbind, lapply(state$area_results, `[[`, "autocorrelation"))
  primary <- ac$area_role == "PRIMARY_DLBCL" & ac$eligibility == "PRIMARY_ELIGIBLE"
  x <- ac[primary, , drop = FALSE]
  state$across_area_summary <- do.call(rbind, lapply(split(x, x$program_id), function(z) data.frame(
    program_id = z$program_id[1L], eligible_area_count = nrow(z), median_Moran_I = median(z$Moran_I),
    minimum_Moran_I = min(z$Moran_I), maximum_Moran_I = max(z$Moran_I),
    positive_direction_count = sum(z$Moran_I > 0), synthesis = "DESCRIPTIVE_ONLY",
    patient_meta_analysis = FALSE, stringsAsFactors = FALSE)))
  state
}

if (sys.nframe() == 0L) wp3_node_cli("36k", wp3_node_36k)
