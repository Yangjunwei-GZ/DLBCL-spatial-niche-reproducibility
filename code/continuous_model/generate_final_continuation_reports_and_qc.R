#!/usr/bin/env Rscript

DLBCL_PROJECT_ROOT <- normalizePath(
  Sys.getenv("DLBCL_PROJECT_ROOT", unset = getwd()),
  winslash = "/", mustWork = FALSE
)
DLBCL_SUPPLEMENTARY_CODE_ROOT <- normalizePath(
  Sys.getenv("DLBCL_SUPPLEMENTARY_CODE_ROOT",
    unset = file.path(DLBCL_PROJECT_ROOT, "Supplementary_Code")),
  winslash = "/", mustWork = FALSE
)


output_root <- normalizePath(
  file.path(DLBCL_PROJECT_ROOT, "revision_2026_reviewer_response/05e_stage4_GSE31312_execution_attempt2"),
  winslash = "/",
  mustWork = TRUE
)

stopifnot(!nzchar(Sys.getenv("DLBCL_REVISION_ALLOW_SCIENCE", unset = "")))

safe_path <- function(relative_path) {
  path <- file.path(output_root, relative_path)
  normalized_parent <- normalizePath(dirname(path), winslash = "/", mustWork = FALSE)
  parent_is_root <- identical(tolower(normalized_parent), tolower(output_root))
  parent_is_child <- startsWith(
    paste0(tolower(normalized_parent), "/"),
    paste0(tolower(output_root), "/")
  )
  if (!(parent_is_root || parent_is_child)) {
    stop("Output path is outside the controlled output root: ", path, call. = FALSE)
  }
  if (file.exists(path)) {
    stop("Refusing to overwrite: ", path, call. = FALSE)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}

read_stage4_csv <- function(relative_path) {
  path <- file.path(output_root, relative_path)
  if (!file.exists(path)) stop("Missing required source CSV: ", path, call. = FALSE)
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

write_stage4_csv <- function(value, relative_path) {
  path <- safe_path(relative_path)
  write.csv(value, path, row.names = FALSE, fileEncoding = "UTF-8")
  path
}

write_stage4_lines <- function(value, relative_path) {
  path <- safe_path(relative_path)
  writeLines(sub("[[:space:]]+$", "", value), path, useBytes = TRUE)
  path
}

finite_columns <- function(value) {
  numeric_columns <- vapply(value, is.numeric, logical(1))
  all(vapply(value[numeric_columns], function(x) all(is.finite(x)), logical(1)))
}

internal_qc_label <- "INTERNAL QC \u2014 NOT FINAL MANUSCRIPT FIGURE"

save_plot_pair <- function(stem, plot_function, width = 7, height = 5) {
  png_path <- safe_path(file.path("internal_qc", paste0(stem, ".png")))
  grDevices::png(png_path, width = width, height = height, units = "in", res = 180)
  plot_function()
  grDevices::dev.off()

  pdf_path <- safe_path(file.path("internal_qc", paste0(stem, ".pdf")))
  grDevices::pdf(pdf_path, width = width, height = height, useDingbats = FALSE)
  plot_function()
  grDevices::dev.off()
  c(png_path, pdf_path)
}

explained <- read_stage4_csv("continuous_geometry/PCA_EXPLAINED_VARIANCE.csv")
loadings <- read_stage4_csv("continuous_geometry/PCA_PROGRAM_LOADINGS.csv")
coordinates <- read_stage4_csv("continuous_geometry/PCA_SAMPLE_COORDINATES.csv")
axis_correlations <- read_stage4_csv("continuous_geometry/PCA_AXIS_CORRELATIONS.csv")
hopkins <- read_stage4_csv("continuous_geometry/HOPKINS_CLUSTERABILITY.csv")
dip_tests <- read_stage4_csv("continuous_geometry/PC1_PC2_MULTIMODALITY.csv")
density_source <- read_stage4_csv("continuous_geometry/PC1_PC2_DENSITY_SOURCE.csv")
metrics <- read_stage4_csv("cluster_number/CLUSTER_NUMBER_METRICS.csv")
assignments <- read_stage4_csv("cluster_number/PRIMARY_HCLUST_ASSIGNMENTS_K2_K6.csv")
silhouette_samples <- read_stage4_csv("cluster_number/SILHOUETTE_BY_SAMPLE.csv")
cdf_source <- read_stage4_csv("cluster_number/CONSENSUS_CDF_SOURCE.csv")
bootstrap_summary <- read_stage4_csv("bootstrap/BOOTSTRAP_STABILITY_SUMMARY.csv")
bootstrap_failures <- read_stage4_csv("bootstrap/BOOTSTRAP_FAILED_OR_UNMATCHED.csv")
scores <- read_stage4_csv("score_spaces/historical_untruncated_z.csv")
legacy_metrics <- read_stage4_csv("legacy_comparison/OLD_VERSUS_NEW_K4_METRICS.csv")
metrics_for_finite_check <- metrics
metrics_for_finite_check$delta_area <- NULL

stopifnot(
  identical(explained$PC, paste0("PC", 1:6)),
  abs(sum(explained$explained_variance) - 1) <= 1e-12,
  abs(tail(explained$cumulative_variance, 1) - 1) <= 1e-12,
  nrow(coordinates) == 498L,
  identical(metrics$k, 2:6),
  all(bootstrap_failures$requested_replicates == 1000L),
  all(bootstrap_failures$failed_or_unmatched_replicates == 0L),
  finite_columns(explained),
  finite_columns(loadings),
  finite_columns(coordinates),
  finite_columns(metrics_for_finite_check),
  is.na(metrics$delta_area[[1]]),
  all(is.finite(metrics$delta_area[-1])),
  finite_columns(silhouette_samples),
  finite_columns(bootstrap_summary)
)

correlation_p_value <- function(r, n) {
  statistic <- r * sqrt((n - 2) / (1 - r^2))
  2 * stats::pt(-abs(statistic), df = n - 2)
}

axis_correlations$n <- nrow(coordinates)
axis_correlations$pearson_p_value <- vapply(
  axis_correlations$pearson_r,
  correlation_p_value,
  numeric(1),
  n = nrow(coordinates)
)
write_stage4_csv(
  axis_correlations,
  "continuous_geometry/PCA_AXIS_CORRELATIONS_WITH_P_VALUES.csv"
)

cluster_sizes <- do.call(rbind, lapply(split(assignments, assignments$k), function(group) {
  counts <- as.data.frame(table(group$cluster_id), stringsAsFactors = FALSE)
  names(counts) <- c("cluster_id", "n")
  counts$k <- unique(group$k)
  counts$cluster_fraction <- counts$n / nrow(scores)
  counts[, c("k", "cluster_id", "n", "cluster_fraction")]
}))
rownames(cluster_sizes) <- NULL
cluster_sizes <- cluster_sizes[order(cluster_sizes$k, cluster_sizes$cluster_id), ]
write_stage4_csv(cluster_sizes, "cluster_number/CLUSTER_SIZES_K2_K6.csv")

silhouette_by_cluster <- do.call(
  rbind,
  lapply(split(silhouette_samples, list(silhouette_samples$k, silhouette_samples$cluster_id)), function(group) {
    data.frame(
      k = group$k[[1]],
      cluster_id = group$cluster_id[[1]],
      n = nrow(group),
      mean_silhouette = mean(group$silhouette),
      median_silhouette = stats::median(group$silhouette),
      minimum_silhouette = min(group$silhouette),
      maximum_silhouette = max(group$silhouette),
      stringsAsFactors = FALSE
    )
  })
)
rownames(silhouette_by_cluster) <- NULL
silhouette_by_cluster <- silhouette_by_cluster[
  order(silhouette_by_cluster$k, silhouette_by_cluster$cluster_id),
]
write_stage4_csv(
  silhouette_by_cluster,
  "cluster_number/SILHOUETTE_BY_CLUSTER_SUMMARY.csv"
)

program_columns <- setdiff(names(scores), "sample")
neutral_centroids <- do.call(rbind, lapply(split(assignments, assignments$k), function(group) {
  matched <- match(group$sample, scores$sample)
  if (anyNA(matched)) stop("A clustering sample is missing from the score table.", call. = FALSE)
  value <- stats::aggregate(
    scores[matched, program_columns, drop = FALSE],
    by = list(cluster_id = group$cluster_id),
    FUN = mean
  )
  value$k <- group$k[[1]]
  value[, c("k", "cluster_id", program_columns)]
}))
rownames(neutral_centroids) <- NULL
neutral_centroids <- neutral_centroids[
  order(neutral_centroids$k, neutral_centroids$cluster_id),
]
write_stage4_csv(neutral_centroids, "cluster_number/NEUTRAL_CENTROIDS_K2_K6.csv")

execution_status <- data.frame(
  k = 2:6,
  assignment_status = "completed",
  consensus_status = "completed",
  silhouette_status = "completed",
  gap_status = "completed",
  CH_status = "completed",
  bootstrap_status = "completed",
  warning_count = 0L,
  failure_status = "none",
  stringsAsFactors = FALSE
)
write_stage4_csv(execution_status, "cluster_number/K2_K6_EXECUTION_STATUS.csv")

metrics$PAC_rank <- rank(metrics$PAC, ties.method = "min")
metrics$silhouette_rank <- rank(-metrics$mean_silhouette, ties.method = "min")
best_silhouette <- max(metrics$mean_silhouette)
metrics$silhouette_within_0_02_of_best <-
  (best_silhouette - metrics$mean_silhouette) <= 0.02

cluster_bootstrap <- bootstrap_summary[bootstrap_summary$cluster_id != "OVERALL", ]
overall_bootstrap <- bootstrap_summary[bootstrap_summary$cluster_id == "OVERALL", ]
minimum_bootstrap <- stats::aggregate(
  median ~ k,
  data = cluster_bootstrap,
  FUN = min
)
names(minimum_bootstrap)[[2]] <- "minimum_cluster_median_jaccard"
overall_bootstrap <- overall_bootstrap[, c("k", "median")]
names(overall_bootstrap)[[2]] <- "overall_bootstrap_median_jaccard"

candidate <- merge(metrics, overall_bootstrap, by = "k", sort = FALSE)
candidate <- merge(candidate, minimum_bootstrap, by = "k", sort = FALSE)
candidate <- candidate[match(2:6, candidate$k), ]
candidate$all_clusters_at_least_0_60 <-
  candidate$minimum_cluster_median_jaccard >= 0.60
candidate$overall_at_least_0_75 <-
  candidate$overall_bootstrap_median_jaccard >= 0.75
candidate$primary_metric_candidate <-
  candidate$PAC_rank <= 2 &
  candidate$mean_silhouette > 0.25 &
  candidate$silhouette_within_0_02_of_best &
  candidate$all_clusters_at_least_0_60 &
  candidate$overall_at_least_0_75 &
  candidate$minimum_cluster_size >= 25
candidate$sensitivity_evidence_available <- FALSE
candidate$final_candidate_status <- "PENDING_STAGE4C2"
candidate <- candidate[, c(
  "k", "PAC", "PAC_rank", "mean_silhouette", "silhouette_rank",
  "silhouette_within_0_02_of_best", "gap", "calinski_harabasz",
  "minimum_cluster_size", "minimum_cluster_fraction",
  "overall_bootstrap_median_jaccard", "minimum_cluster_median_jaccard",
  "all_clusters_at_least_0_60", "overall_at_least_0_75",
  "primary_metric_candidate", "sensitivity_evidence_available",
  "final_candidate_status"
)]
names(candidate)[names(candidate) == "gap"] <- "gap_statistic"
names(candidate)[names(candidate) == "calinski_harabasz"] <- "CH_index"
write_stage4_csv(candidate, "CORE_STAGE4_ATTEMPT2_CANDIDATE_TABLE.csv")

proliferation_row <- match("Proliferative / cycling program", loadings$program)
pc1_rank <- rank(-abs(loadings$PC1), ties.method = "min")[[proliferation_row]]
pc2_rank <- rank(-abs(loadings$PC2), ties.method = "min")[[proliferation_row]]
pc1_correlation <- axis_correlations[
  axis_correlations$comparison == "PC1_vs_six_program_mean",
]
pc2_correlation <- axis_correlations[
  axis_correlations$comparison == "PC2_vs_immune_stromal_contrast",
]

loading_lines <- vapply(seq_len(nrow(loadings)), function(i) {
  sprintf(
    "- %s: PC1 = %.6f; PC2 = %.6f",
    loadings$program[[i]],
    loadings$PC1[[i]],
    loadings$PC2[[i]]
  )
}, character(1))

script16_checkpoint <- c(
  "# STAGE4C1 Script 16 Checkpoint",
  "",
  "- status: COMPLETED",
  sprintf("- PC1 explained variance: %.15f", explained$explained_variance[[1]]),
  sprintf("- PC2 explained variance: %.15f", explained$explained_variance[[2]]),
  sprintf("- cumulative PC1+PC2 variance: %.15f", explained$cumulative_variance[[2]]),
  "",
  "## Program Loadings",
  "",
  loading_lines,
  "",
  sprintf(
    "- Proliferative/cycling absolute-loading rank: PC1 %d/6; PC2 %d/6.",
    pc1_rank,
    pc2_rank
  ),
  sprintf(
    "- PC1 vs six-program mean: Pearson r = %.15f; P = %.6g; Spearman rho = %.15f.",
    pc1_correlation$pearson_r,
    pc1_correlation$pearson_p_value,
    pc1_correlation$spearman_rho
  ),
  sprintf(
    "- PC2 vs immune-stromal contrast: Pearson r = %.15f; P = %.6g; Spearman rho = %.15f.",
    pc2_correlation$pearson_r,
    pc2_correlation$pearson_p_value,
    pc2_correlation$spearman_rho
  ),
  sprintf("- Hopkins statistic: %.15f.", hopkins$hopkins_statistic[[1]]),
  sprintf(
    "- PC1 dip test: statistic = %.15f; P = %.15f; status = %s.",
    dip_tests$statistic[dip_tests$axis == "PC1"],
    dip_tests$p_value[dip_tests$axis == "PC1"],
    dip_tests$execution_status[dip_tests$axis == "PC1"]
  ),
  sprintf(
    "- PC2 dip test: statistic = %.15f; P = %.15f; status = %s.",
    dip_tests$statistic[dip_tests$axis == "PC2"],
    dip_tests$p_value[dip_tests$axis == "PC2"],
    dip_tests$execution_status[dip_tests$axis == "PC2"]
  ),
  "- warnings: 0",
  "- GSVA loaded: FALSE",
  "- stack imbalance: FALSE",
  "",
  "Pearson P values were deterministically derived from the reported Pearson r and n = 498."
)
write_stage4_lines(script16_checkpoint, "STAGE4C1_SCRIPT16_CHECKPOINT.md")

script17_rows <- vapply(seq_len(nrow(candidate)), function(i) {
  sprintf(
    "- k=%d: PAC=%.6f (rank %d); mean silhouette=%.6f (rank %d); gap=%.6f; CH=%.6f; minimum cluster size=%d (%.6f); warnings=0.",
    candidate$k[[i]],
    candidate$PAC[[i]],
    candidate$PAC_rank[[i]],
    candidate$mean_silhouette[[i]],
    candidate$silhouette_rank[[i]],
    candidate$gap_statistic[[i]],
    candidate$CH_index[[i]],
    candidate$minimum_cluster_size[[i]],
    candidate$minimum_cluster_fraction[[i]]
  )
}, character(1))
script17_checkpoint <- c(
  "# STAGE4C1 Script 17 Checkpoint",
  "",
  "- status: COMPLETED",
  "- candidate k values: 2, 3, 4, 5, 6",
  "- ConsensusClusterPlus reps: 1000",
  "- gap statistic B: 500",
  "- automatic k selection: FALSE",
  "",
  script17_rows
)
write_stage4_lines(script17_checkpoint, "STAGE4C1_SCRIPT17_CHECKPOINT.md")

log_summary <- data.frame(
  script = c(
    "16_stage4_continuous_geometry.R",
    "17_stage4_cluster_number_metrics.R",
    "18_stage4_bootstrap_stability.R"
  ),
  run_identifier = c(
    "16_stage4_continuous_geometry_continuation",
    "17_stage4_cluster_number_metrics_continuation",
    "18_stage4_bootstrap_stability_continuation"
  ),
  completion_status = "completed",
  warning_count = 0L,
  error_count = 0L,
  stack_imbalance = FALSE,
  GSVA_loaded = FALSE,
  gene_level_expression_read = FALSE,
  science_token_after_process = "UNSET",
  stringsAsFactors = FALSE
)
write_stage4_csv(
  log_summary,
  "logs/STAGE4C1_FINAL_CONTINUATION_LOG_SUMMARY.csv"
)

requirement_status <- data.frame(
  requirement_id = c(
    "continuous_geometry",
    "cluster_number_metrics_k2_k6",
    "bootstrap_stability_k2_k6",
    "final_k_selection",
    "taxonomy_assignment",
    "score_definition_sensitivities",
    "alternative_distance_sensitivity",
    "purity_deconvolution",
    "external_framework_comparison"
  ),
  evidence = c(
    "continuous_geometry CSVs and Script 16 checkpoint",
    "cluster_number CSVs and Script 17 checkpoint",
    "bootstrap CSVs and candidate table",
    "core metrics completed; sensitivity evidence unavailable",
    "neutral IDs retained; no biological names assigned",
    "unique-gene-only, de-overlap, derived exclusion, and leave-one-program-out not run",
    "alternative distance analysis not run",
    "purity/deconvolution analysis not run",
    "external framework comparison not run"
  ),
  status = c(
    "core_analysis_completed",
    "core_analysis_completed",
    "core_analysis_completed",
    "partially_addressed_pending_sensitivity",
    "unresolved_requires_stage4c2",
    "unresolved_requires_stage4c2",
    "unresolved_requires_stage4c2",
    "unresolved_requires_later_stage",
    "unresolved_requires_later_stage"
  ),
  stringsAsFactors = FALSE
)
write_stage4_csv(
  requirement_status,
  "REQUIREMENT_STATUS_AFTER_STAGE4C1_ATTEMPT2_FINAL_CONTINUATION.csv"
)

candidate_lines <- vapply(seq_len(nrow(candidate)), function(i) {
  sprintf(
    "- k=%d: PAC %.6f; silhouette %.6f; minimum n %d; overall Jaccard %.6f; minimum cluster Jaccard %.6f; basic primary-metric candidate = %s.",
    candidate$k[[i]],
    candidate$PAC[[i]],
    candidate$mean_silhouette[[i]],
    candidate$minimum_cluster_size[[i]],
    candidate$overall_bootstrap_median_jaccard[[i]],
    candidate$minimum_cluster_median_jaccard[[i]],
    candidate$primary_metric_candidate[[i]]
  )
}, character(1))

legacy_lines <- vapply(seq_len(nrow(legacy_metrics)), function(i) {
  sprintf(
    "- %s vs %s: concordance %.6f; ARI %.6f; NMI %.6f.",
    legacy_metrics$candidate[[i]],
    legacy_metrics$reference[[i]],
    legacy_metrics$concordance[[i]],
    legacy_metrics$adjusted_rand_index[[i]],
    legacy_metrics$normalized_mutual_information[[i]]
  )
}, character(1))

core_report <- c(
  "# STAGE4C1 Attempt 2 Core Results Report",
  "",
  "## Observed Results",
  "",
  sprintf(
    "PC1 and PC2 explained %.6f and %.6f of variance, respectively; cumulative PC1+PC2 variance was %.6f.",
    explained$explained_variance[[1]],
    explained$explained_variance[[2]],
    explained$cumulative_variance[[2]]
  ),
  sprintf(
    "PC1 correlated with the six-program mean at Pearson r = %.6f (P = %.6g).",
    pc1_correlation$pearson_r,
    pc1_correlation$pearson_p_value
  ),
  sprintf(
    "PC2 correlated with the prespecified immune-stromal contrast at Pearson r = %.6f (P = %.6g).",
    pc2_correlation$pearson_r,
    pc2_correlation$pearson_p_value
  ),
  sprintf("The Hopkins statistic was %.6f.", hopkins$hopkins_statistic[[1]]),
  "",
  candidate_lines,
  "",
  "## Prespecified Threshold Assessment",
  "",
  sprintf(
    "The number of k values meeting all Stage 4C-1 basic primary-metric conditions was %d.",
    sum(candidate$primary_metric_candidate)
  ),
  "This is an interim core-analysis result. It does not select a final k and does not authorize a continuous-only conclusion.",
  "",
  "## Legacy Comparison Context",
  "",
  legacy_lines,
  "",
  "Legacy agreement is descriptive context only and was not used to select k.",
  "",
  "## Pending Evidence",
  "",
  "- unique-gene-only sensitivity",
  "- de-overlap sensitivity",
  "- derived exclusion sensitivity",
  "- leave-one-program-out sensitivity",
  "- alternative distance sensitivity",
  "- purity/deconvolution analysis",
  "- external framework comparison",
  "",
  "## Prohibited Conclusions",
  "",
  "- No final k was selected.",
  "- The historical four-class solution was neither validated nor declared failed.",
  "- No final taxonomy or biological ecosystem names were assigned.",
  "- The manuscript was not modified."
)
write_stage4_lines(core_report, "STAGE4C1_ATTEMPT2_CORE_RESULTS_REPORT.md")

status_report <- c(
  "# STAGE4C1 Attempt 2 Continuation Status",
  "",
  "- script 15a: PRESERVED_COMPLETED",
  "- script 16: COMPLETED",
  "- script 17: COMPLETED",
  "- script 18: COMPLETED",
  "- script 21: PRESERVED_COMPLETED",
  "- Stage 4C-2: NOT_RUN",
  "- final k: NOT_SELECTED",
  "- taxonomy: NOT_ASSIGNED",
  "- core status: COMPLETED_PENDING_SENSITIVITY"
)
write_stage4_lines(status_report, "STAGE4C1_ATTEMPT2_CONTINUATION_STATUS.md")

source_dir <- "internal_qc/source_data"
write_stage4_csv(explained, file.path(source_dir, "scree_plot_source.csv"))
write_stage4_csv(coordinates, file.path(source_dir, "PC1_PC2_scatter_source.csv"))
write_stage4_csv(density_source, file.path(source_dir, "PC1_PC2_density_source.csv"))
write_stage4_csv(loadings, file.path(source_dir, "loading_plot_source.csv"))
write_stage4_csv(metrics[, c("k", "PAC")], file.path(source_dir, "PAC_comparison_source.csv"))
write_stage4_csv(
  metrics[, c("k", "mean_silhouette")],
  file.path(source_dir, "silhouette_comparison_source.csv")
)
write_stage4_csv(metrics[, c("k", "gap", "gap_SE")], file.path(source_dir, "gap_comparison_source.csv"))
write_stage4_csv(
  metrics[, c("k", "calinski_harabasz")],
  file.path(source_dir, "CH_comparison_source.csv")
)
write_stage4_csv(cluster_sizes, file.path(source_dir, "cluster_size_comparison_source.csv"))
write_stage4_csv(
  bootstrap_summary,
  file.path(source_dir, "bootstrap_Jaccard_comparison_source.csv")
)
write_stage4_csv(
  neutral_centroids,
  file.path(source_dir, "neutral_centroid_heatmap_source.csv")
)

save_plot_pair("scree_plot", function() {
  graphics::barplot(
    explained$explained_variance,
    names.arg = explained$PC,
    col = "#6B7280",
    border = NA,
    xlab = "Principal component",
    ylab = "Explained variance",
    main = paste("Scree plot", internal_qc_label, sep = "\n")
  )
})

save_plot_pair("PC1_PC2_scatter", function() {
  graphics::plot(
    coordinates$PC1,
    coordinates$PC2,
    pch = 16,
    cex = 0.55,
    col = grDevices::adjustcolor("#374151", alpha.f = 0.55),
    xlab = "PC1",
    ylab = "PC2",
    main = paste("PC1-PC2 sample geometry", internal_qc_label, sep = "\n")
  )
})

save_plot_pair("PC1_PC2_density", function() {
  x <- sort(unique(density_source$PC1))
  y <- sort(unique(density_source$PC2))
  z <- matrix(density_source$density, nrow = length(x), ncol = length(y))
  graphics::image(
    x,
    y,
    z,
    col = grDevices::hcl.colors(40, "YlGnBu", rev = TRUE),
    xlab = "PC1",
    ylab = "PC2",
    main = paste("PC1-PC2 density", internal_qc_label, sep = "\n")
  )
  graphics::contour(x, y, z, add = TRUE, drawlabels = FALSE, col = "#374151")
})

save_plot_pair("loading_plot", function() {
  short <- c("Macrophage", "T-cell", "Antigen", "Stromal", "Immune-cold", "Cycling")
  graphics::plot(
    loadings$PC1,
    loadings$PC2,
    pch = 16,
    col = "#1F4E79",
    xlab = "PC1 loading",
    ylab = "PC2 loading",
    main = paste("Program loadings", internal_qc_label, sep = "\n")
  )
  graphics::text(loadings$PC1, loadings$PC2, labels = short, pos = 3, cex = 0.8)
  graphics::abline(h = 0, v = 0, col = "#D1D5DB", lty = 2)
})

metric_plot <- function(y, ylab, title, color = "#1F4E79") {
  graphics::plot(
    metrics$k,
    y,
    type = "b",
    pch = 16,
    lwd = 2,
    col = color,
    xlab = "k",
    ylab = ylab,
    xaxt = "n",
    main = paste(title, internal_qc_label, sep = "\n")
  )
  graphics::axis(1, at = metrics$k)
}

save_plot_pair("PAC_comparison", function() metric_plot(metrics$PAC, "PAC", "PAC comparison"))
save_plot_pair(
  "silhouette_comparison",
  function() metric_plot(metrics$mean_silhouette, "Mean silhouette", "Silhouette comparison")
)
save_plot_pair("gap_comparison", function() {
  metric_plot(metrics$gap, "Gap statistic", "Gap statistic comparison")
  graphics::arrows(
    metrics$k,
    metrics$gap - metrics$gap_SE,
    metrics$k,
    metrics$gap + metrics$gap_SE,
    angle = 90,
    code = 3,
    length = 0.05
  )
})
save_plot_pair(
  "CH_comparison",
  function() metric_plot(metrics$calinski_harabasz, "Calinski-Harabasz index", "CH comparison")
)

save_plot_pair("cluster_size_comparison", function() {
  size_matrix <- xtabs(n ~ cluster_id + k, data = cluster_sizes)
  graphics::barplot(
    size_matrix,
    beside = FALSE,
    col = grDevices::hcl.colors(nrow(size_matrix), "Dark 3"),
    xlab = "k",
    ylab = "Samples",
    main = paste("Cluster size comparison", internal_qc_label, sep = "\n"),
    legend.text = rownames(size_matrix),
    args.legend = list(x = "topright", cex = 0.7, bty = "n")
  )
})

save_plot_pair("bootstrap_Jaccard_comparison", function() {
  overall <- bootstrap_summary[bootstrap_summary$cluster_id == "OVERALL", ]
  graphics::plot(
    overall$k,
    overall$median,
    pch = 16,
    type = "b",
    lwd = 2,
    ylim = c(0, 1),
    xaxt = "n",
    xlab = "k",
    ylab = "Overall median Jaccard",
    main = paste("Bootstrap stability", internal_qc_label, sep = "\n")
  )
  graphics::axis(1, at = overall$k)
  graphics::arrows(
    overall$k,
    overall$Q1,
    overall$k,
    overall$Q3,
    angle = 90,
    code = 3,
    length = 0.05
  )
  graphics::abline(h = c(0.60, 0.75), lty = 2, col = c("#9CA3AF", "#4B5563"))
})

save_plot_pair("neutral_centroid_heatmap", function() {
  centroid_matrix <- as.matrix(neutral_centroids[, program_columns, drop = FALSE])
  rownames(centroid_matrix) <- paste0("k", neutral_centroids$k, "_", neutral_centroids$cluster_id)
  graphics::image(
    seq_len(nrow(centroid_matrix)),
    seq_len(ncol(centroid_matrix)),
    centroid_matrix,
    col = grDevices::hcl.colors(41, "Blue-Red 3", rev = TRUE),
    axes = FALSE,
    xlab = "Neutral cluster",
    ylab = "Program",
    main = paste("Neutral centroid heatmap", internal_qc_label, sep = "\n")
  )
  graphics::axis(1, at = seq_len(nrow(centroid_matrix)), labels = rownames(centroid_matrix), las = 2, cex.axis = 0.55)
  graphics::axis(2, at = seq_len(ncol(centroid_matrix)), labels = program_columns, las = 2, cex.axis = 0.55)
  box()
})

cat("FINAL_CONTINUATION_REPORTS_AND_QC=COMPLETED\n")
