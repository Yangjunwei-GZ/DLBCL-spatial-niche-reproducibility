DLBCL_PROJECT_ROOT <- normalizePath(
  Sys.getenv("DLBCL_PROJECT_ROOT", unset = getwd()),
  winslash = "/", mustWork = FALSE
)
DLBCL_SUPPLEMENTARY_CODE_ROOT <- normalizePath(
  Sys.getenv("DLBCL_SUPPLEMENTARY_CODE_ROOT",
    unset = file.path(DLBCL_PROJECT_ROOT, "Supplementary_Code")),
  winslash = "/", mustWork = FALSE
)

# Node 29e: compare six-program correlation structure with the frozen GSE31312 authority.

if (!exists("WP2", inherits = TRUE)) source(file.path(DLBCL_SUPPLEMENTARY_CODE_ROOT, "external_validation", "29_common_wp2_functions_v1.R"))

wp2_node_29e <- function() {
  wp2_assert_token()
  discovery_edges <- wp2_discovery_edges()
  discovery_matrix <- matrix(NA_real_, 6L, 6L, dimnames = list(WP2$programs, WP2$programs))
  diag(discovery_matrix) <- 1
  for (i in seq_len(nrow(discovery_edges))) {
    a <- discovery_edges$program_1_order[[i]]
    b <- discovery_edges$program_2_order[[i]]
    discovery_matrix[a, b] <- discovery_edges$pearson_r[[i]]
    discovery_matrix[b, a] <- discovery_edges$pearson_r[[i]]
  }
  if (any(!is.finite(discovery_matrix))) stop("Discovery correlation matrix is incomplete", call. = FALSE)
  wp2_write_csv(wp2_matrix_data_frame(discovery_matrix), file.path(WP2$stage, "05_structure/GSE31312_PEARSON_CORRELATION_MATRIX.csv"))

  primary_edges <- list()
  all_edges <- list()
  point_rows <- list()
  bootstrap_rows <- list()
  bootstrap_summary <- list()

  for (dataset in c("GSE10846", "GSE181063")) {
    primary <- utils::read.csv(file.path(WP2$stage, "03_scores", paste0(dataset, "_PRIMARY_SIX_PROGRAM_SCORES.csv")), check.names = FALSE, stringsAsFactors = FALSE)
    sensitivity <- utils::read.csv(file.path(WP2$stage, "03_scores", paste0(dataset, "_SENSITIVITY_SIX_PROGRAM_SCORES.csv")), check.names = FALSE, stringsAsFactors = FALSE)
    if (!identical(primary$sample, sensitivity$sample)) stop("Primary/sensitivity sample mismatch for ", dataset, call. = FALSE)

    primary_cor <- wp2_correlation_long(primary, dataset, "MEDIAN_ACROSS_PROBES_PER_GENE")
    sensitivity_cor <- wp2_correlation_long(sensitivity, dataset, "HIGHEST_MAD_PROBE_PER_GENE")
    wp2_write_csv(
      wp2_matrix_data_frame(primary_cor$matrix),
      file.path(WP2$stage, "05_structure", paste0(dataset, "_PRIMARY_PEARSON_CORRELATION_MATRIX.csv"))
    )
    primary_edges[[dataset]] <- primary_cor$long
    all_edges[[paste0(dataset, "_primary")]] <- primary_cor$long
    all_edges[[paste0(dataset, "_sensitivity")]] <- sensitivity_cor$long

    metric <- wp2_structural_metrics(primary_cor$long)
    point_rows[[dataset]] <- data.frame(
      dataset_id = dataset,
      sample_count = nrow(primary),
      matrix_vector_spearman = as.numeric(metric[["matrix_vector_spearman"]]),
      sign_concordance_n = as.integer(metric[["sign_concordance_n"]]),
      sign_concordance_denominator = 15L,
      structural_replication_status = unname(metric[["structural_replication_status"]]),
      comparison_authority = "GSE31312_FROZEN_15_EDGE_PEARSON_VECTOR",
      stringsAsFactors = FALSE
    )

    boot <- wp2_bootstrap_structure(primary, dataset)
    bootstrap_rows[[dataset]] <- boot$replicates
    bootstrap_summary[[dataset]] <- boot$summary
  }

  wp2_write_csv(do.call(rbind, all_edges), file.path(WP2$stage, "05_structure/WP2_MAPPING_SENSITIVITY_CORRELATION_MATRICES_LONG.csv"))
  wp2_write_csv(do.call(rbind, primary_edges), file.path(WP2$stage, "05_structure/WP2_PRIMARY_CORRELATION_EDGES_LONG.csv"))
  wp2_write_csv(do.call(rbind, point_rows), file.path(WP2$stage, "05_structure/WP2_STRUCTURAL_REPLICATION_POINT_ESTIMATES.csv"))
  wp2_write_csv(do.call(rbind, bootstrap_rows), file.path(WP2$stage, "05_structure/WP2_STRUCTURAL_BOOTSTRAP_REPLICATES.csv"))
  wp2_write_csv(do.call(rbind, bootstrap_summary), file.path(WP2$stage, "05_structure/WP2_STRUCTURAL_BOOTSTRAP_SUMMARY.csv"))
  invisible(do.call(rbind, point_rows))
}

if (!identical(Sys.getenv("WP2_NODE_SOURCE_ONLY", unset = ""), "TRUE")) wp2_node_29e()
