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

wp3_node_36h <- function(state) {
  scores <- state$primary_scores[match(state$barcodes, state$primary_scores$barcode), , drop = FALSE]
  sct_scores <- state$sct_scores[match(state$barcodes, state$sct_scores$barcode), , drop = FALSE]
  ssgsea_scores <- state$ssgsea_scores[match(state$barcodes, state$ssgsea_scores$barcode), , drop = FALSE]
  rows <- lapply(seq_len(nrow(state$eligibility)), function(i) {
    program <- state$eligibility$program_id[i]
    if (!program %in% names(scores) || state$eligibility$eligibility[i] == "INELIGIBLE") return(NULL)
    m <- wp3_permutation_test_edge(scores[[program]], state$edge_index, "moran", 9999L, WP3_SEED + i)
    g <- wp3_permutation_test_edge(scores[[program]], state$edge_index, "geary", 9999L, WP3_SEED + 100L + i)
    sct_moran <- wp3_moran_edge(sct_scores[[program]], state$edge_index)
    ssgsea_moran <- wp3_moran_edge(ssgsea_scores[[program]], state$edge_index)
    family <- wp3_assign_bh_family(state$area_role, state$eligibility$eligibility[i], "MORAN", state$capture_area_id)
    data.frame(capture_area_id = state$capture_area_id, area_role = state$area_role,
      program_id = program, eligibility = state$eligibility$eligibility[i],
      Moran_I = m$observed, Moran_p = m$permutation_p, Geary_C = g$observed,
      Geary_p = g$permutation_p, SCT_UCell_Moran_I = sct_moran,
      ssGSEA_Moran_I = ssgsea_moran,
      SCT_direction_concordant = sign(m$observed) == sign(sct_moran),
      ssGSEA_direction_concordant = sign(m$observed) == sign(ssgsea_moran),
      permutations = 9999L, bh_family = family)
  })
  state$autocorrelation <- do.call(rbind, rows)
  temp <- transform(state$autocorrelation, permutation_p = Moran_p)
  temp <- wp3_adjust_families(temp)
  state$autocorrelation$Moran_FDR <- temp$fdr
  temp$permutation_p <- temp$Geary_p
  state$autocorrelation$Geary_FDR <- wp3_adjust_families(temp)$fdr
  state
}

if (sys.nframe() == 0L) wp3_node_cli("36h", wp3_node_36h)
