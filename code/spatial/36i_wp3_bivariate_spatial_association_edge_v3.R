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

wp3_node_36i <- function(state) {
  eligible <- state$eligibility[state$eligibility$eligibility != "INELIGIBLE", , drop = FALSE]
  pairs <- combn(eligible$program_id, 2L, simplify = FALSE)
  scores <- state$primary_scores[match(state$barcodes, state$primary_scores$barcode), , drop = FALSE]
  state$bivariate <- do.call(rbind, lapply(seq_along(pairs), function(i) {
    pair <- pairs[[i]]; info <- eligible[match(pair, eligible$program_id), ]
    exploratory <- any(info$eligibility == "EXPLORATORY_ONLY")
    b <- wp3_bivariate_permutation_edge(scores[[pair[1L]]], scores[[pair[2L]]],
      state$edge_index, 9999L, WP3_SEED + 200L + i)
    data.frame(capture_area_id = state$capture_area_id, program_1 = pair[1L], program_2 = pair[2L],
      bivariate_Moran_I = b$observed, permutation_p = b$permutation_p, permutations = 9999L,
      exploratory = exploratory, bh_family = wp3_assign_bh_family(state$area_role,
        if (exploratory) "EXPLORATORY_ONLY" else "PRIMARY_ELIGIBLE", "BIVARIATE",
        state$capture_area_id, contains_antigen = exploratory), stringsAsFactors = FALSE)
  }))
  state$bivariate <- wp3_adjust_families(state$bivariate)
  state
}

if (sys.nframe() == 0L) wp3_node_cli("36i", wp3_node_36i)
