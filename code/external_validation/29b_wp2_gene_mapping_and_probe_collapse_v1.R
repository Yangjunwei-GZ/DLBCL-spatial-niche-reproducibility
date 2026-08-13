DLBCL_PROJECT_ROOT <- normalizePath(
  Sys.getenv("DLBCL_PROJECT_ROOT", unset = getwd()),
  winslash = "/", mustWork = FALSE
)
DLBCL_SUPPLEMENTARY_CODE_ROOT <- normalizePath(
  Sys.getenv("DLBCL_SUPPLEMENTARY_CODE_ROOT",
    unset = file.path(DLBCL_PROJECT_ROOT, "Supplementary_Code")),
  winslash = "/", mustWork = FALSE
)

# Node 29b: clean probe mappings and create frozen primary/sensitivity gene matrices.

if (!exists("WP2", inherits = TRUE)) source(file.path(DLBCL_SUPPLEMENTARY_CODE_ROOT, "external_validation", "29_common_wp2_functions_v1.R"))

wp2_node_29b <- function() {
  wp2_assert_token()
  datasets <- c("GSE10846", "GSE181063")
  excluded_all <- list()
  coverage_all <- list()
  for (dataset in datasets) {
    parsed <- readRDS(file.path(WP2$stage, "local_only", paste0(dataset, "_parsed_probe_expression.rds")))
    annotation_path <- WP2$paths[[paste0(dataset, "_annotation")]]
    cleaned <- wp2_clean_annotation(annotation_path, rownames(parsed$expression))
    collapsed <- wp2_collapse_expression(parsed$expression, cleaned$annotation)
    wp2_save_rds(collapsed$primary, file.path(WP2$stage, "local_only", paste0(dataset, "_primary_median_gene_expression.rds")))
    wp2_save_rds(collapsed$sensitivity, file.path(WP2$stage, "local_only", paste0(dataset, "_sensitivity_highest_MAD_gene_expression.rds")))
    mapping <- collapsed$mapping
    mapping$dataset_id <- dataset
    mapping <- mapping[, c("dataset_id", setdiff(names(mapping), "dataset_id"))]
    wp2_write_csv(mapping, file.path(WP2$stage, "02_mapping", paste0(dataset, "_GENE_PROBE_COLLAPSE_AUDIT.csv")))
    excluded <- cleaned$excluded
    excluded$dataset_id <- dataset
    excluded_all[[dataset]] <- excluded[, c("dataset_id", "PROBEID", "reason")]
    coverage_all[[paste0(dataset, "_primary")]] <- wp2_coverage(dataset, "MEDIAN_ACROSS_PROBES_PER_GENE", rownames(collapsed$primary))
    coverage_all[[paste0(dataset, "_sensitivity")]] <- wp2_coverage(dataset, "HIGHEST_MAD_PROBE_PER_GENE", rownames(collapsed$sensitivity))
    rm(parsed, collapsed)
    gc()
  }
  coverage <- do.call(rbind, coverage_all)
  if (any(coverage$detected_canonical_gene_count < 18L)) stop("Primary or sensitivity coverage below 18/22", call. = FALSE)
  if (any(nzchar(coverage$substituted_genes))) stop("Gene substitution detected", call. = FALSE)
  wp2_write_csv(do.call(rbind, excluded_all), file.path(WP2$stage, "02_mapping/WP2_EXCLUDED_PROBES.csv"))
  wp2_write_csv(coverage, file.path(WP2$stage, "WP2_REAL_CANONICAL_COVERAGE.csv"))
  invisible(coverage)
}

if (!identical(Sys.getenv("WP2_NODE_SOURCE_ONLY", unset = ""), "TRUE")) wp2_node_29b()
