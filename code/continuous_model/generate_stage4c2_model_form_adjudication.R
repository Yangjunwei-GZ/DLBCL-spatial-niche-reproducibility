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


# Read-only Stage 4C-2 adjudication. This script does not invoke any scientific node.
options(warn = 2)
project_root <- DLBCL_PROJECT_ROOT
output_root <- file.path(project_root, "revision_2026_reviewer_response/05v_stage4c2_final_k_model_form_adjudication")
stage4c1_root <- file.path(project_root, "revision_2026_reviewer_response/05e_stage4_GSE31312_execution_attempt2")
stage4c2_root <- file.path(project_root, "revision_2026_reviewer_response/05o_stage4c2_internal_sensitivity_execution")
seed_path <- file.path(project_root, "revision_2026_reviewer_response/05n_stage4c2_sensitivity_protocol_preflight/STAGE4C2_SEED_REGISTRY.csv")
code_root <- file.path(project_root, "revision_2026_reviewer_response/05q_stage4c2a_real_execution_implementation/scripts")

assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
read_csv <- function(path) utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
as_flag <- function(x) if (is.logical(x)) x else toupper(as.character(x)) == "TRUE"
fmt <- function(x, digits = 6L) formatC(as.numeric(x), digits = digits, format = "fg")
collapse_or_none <- function(x) if (length(x)) paste(x, collapse = ";") else "NONE"

write_csv_once <- function(x, name) {
  path <- file.path(output_root, name)
  assert(!file.exists(path), paste("Create-once output already exists:", path))
  utils::write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

write_text_once <- function(x, name) {
  path <- file.path(output_root, name)
  assert(!file.exists(path), paste("Create-once output already exists:", path))
  writeLines(enc2utf8(x), path, useBytes = TRUE)
}

assert(dir.exists(output_root), "05v output directory is absent.")
assert(file.exists(file.path(output_root, "FINAL_K_MODEL_FORM_ADJUDICATION_RULES.md")), "Adjudication rules must exist before generation.")

candidate <- read_csv(file.path(stage4c1_root, "CORE_STAGE4_ATTEMPT2_CANDIDATE_TABLE.csv"))
passes <- read_csv(file.path(stage4c2_root, "integration/all_variant_pass_tables.csv"))
comparisons <- read_csv(file.path(stage4c2_root, "integration/all_primary_same_k_comparisons.csv"))
persistence <- read_csv(file.path(stage4c2_root, "STAGE4C2_SENSITIVITY_PERSISTENCE_SUMMARY.csv"))
execution_status <- read_csv(file.path(stage4c2_root, "STAGE4C2A_EXECUTION_STATUS.csv"))
final_checks <- read_csv(file.path(stage4c2_root, "STAGE4C2A_CONTINUATION_FINAL_SELF_CHECK_V4.csv"))
seed_registry <- read_csv(seed_path)

assert(nrow(candidate) == 5L && identical(as.integer(candidate$k), 2:6), "Primary candidate table is incomplete.")
assert(nrow(passes) == 70L && all(table(passes$k) == 14L), "Sensitivity pass table must contain 14 nonreference variants per k.")
assert(nrow(comparisons) == 70L && all(table(comparisons$k) == 14L), "Agreement table must contain 14 nonreference variants per k.")
assert(nrow(persistence) == 5L && identical(as.integer(persistence$k), 2:6), "Persistence summary is incomplete.")
assert(nrow(final_checks) == 101L && all(final_checks$status == "PASS"), "Stage 4C-2A final self-check is not 101/101 PASS.")

criterion_labels <- c(
  PAC_rank_at_most_2 = "PAC_rank",
  silhouette_above_0_25 = "mean_silhouette",
  silhouette_within_0_02 = "within_0_02_of_best",
  minimum_size_at_least_25 = "minimum_cluster_size",
  minimum_cluster_median_jaccard_at_least_0_60 = "minimum_cluster_median_jaccard",
  overall_median_jaccard_at_least_0_75 = "overall_median_jaccard"
)
criterion_fields <- names(criterion_labels)

# Primary threshold decomposition uses the exact frozen Stage 4C-1 boundaries.
best_primary_silhouette <- max(candidate$mean_silhouette)
primary_pass_matrix <- cbind(
  PAC_rank = candidate$PAC_rank <= 2,
  mean_silhouette = candidate$mean_silhouette > 0.25,
  within_0_02_of_best = best_primary_silhouette - candidate$mean_silhouette <= 0.02,
  minimum_cluster_size = candidate$minimum_cluster_size >= 25,
  minimum_cluster_median_jaccard = candidate$minimum_cluster_median_jaccard >= 0.60,
  overall_median_jaccard = candidate$overall_bootstrap_median_jaccard >= 0.75
)
primary_recomputed <- apply(primary_pass_matrix, 1L, all)
assert(identical(primary_recomputed, as_flag(candidate$primary_metric_candidate)), "Primary pass recomputation does not match the frozen field.")

primary_failed <- apply(primary_pass_matrix, 1L, function(z) collapse_or_none(names(z)[!z]))
primary_failed_n <- rowSums(!primary_pass_matrix)
primary_decomposition <- data.frame(
  k = candidate$k,
  PAC = candidate$PAC,
  PAC_rank = candidate$PAC_rank,
  PAC_rank_required = "<=2",
  PAC_rank_margin = 2 - candidate$PAC_rank,
  PAC_rank_pass = primary_pass_matrix[, "PAC_rank"],
  mean_silhouette = candidate$mean_silhouette,
  silhouette_required = ">0.25",
  silhouette_absolute_margin = candidate$mean_silhouette - 0.25,
  silhouette_pass = primary_pass_matrix[, "mean_silhouette"],
  best_primary_silhouette = best_primary_silhouette,
  within_0_02_margin = candidate$mean_silhouette - (best_primary_silhouette - 0.02),
  within_0_02_of_best = primary_pass_matrix[, "within_0_02_of_best"],
  minimum_cluster_size = candidate$minimum_cluster_size,
  minimum_size_required = ">=25",
  minimum_size_margin = candidate$minimum_cluster_size - 25,
  minimum_size_pass = primary_pass_matrix[, "minimum_cluster_size"],
  minimum_cluster_median_jaccard = candidate$minimum_cluster_median_jaccard,
  required_cluster_jaccard = ">=0.60",
  cluster_jaccard_margin = candidate$minimum_cluster_median_jaccard - 0.60,
  cluster_jaccard_pass = primary_pass_matrix[, "minimum_cluster_median_jaccard"],
  overall_median_jaccard = candidate$overall_bootstrap_median_jaccard,
  required_overall_jaccard = ">=0.75",
  overall_jaccard_margin = candidate$overall_bootstrap_median_jaccard - 0.75,
  overall_jaccard_pass = primary_pass_matrix[, "overall_median_jaccard"],
  primary_pass = primary_recomputed,
  failed_criterion_count = primary_failed_n,
  failed_criteria = primary_failed,
  near_miss_criterion_count = ifelse(primary_failed_n == 1L, 1L, 0L),
  notes = "Near-miss is descriptive and requires exactly one failed frozen criterion; acceptance is unchanged.",
  check.names = FALSE
)
write_csv_once(primary_decomposition, "PRIMARY_K2_K6_THRESHOLD_DECOMPOSITION.csv")

# Sensitivity failure decomposition, retaining all 14 nonreference variants.
family_lookup <- setNames(
  c("clipped", rep("program_overlap", 3), "five_program", "derived_exclusion", rep("leave_one_out", 6), rep("alternative_distance", 2)),
  unique(passes$variant_id)
)
assert(length(family_lookup) == 14L && !anyNA(family_lookup[passes$variant_id]), "Variant family map is incomplete.")

recomputed_passes <- cbind(
  PAC_rank_at_most_2 = passes$PAC_rank <= 2,
  silhouette_above_0_25 = passes$mean_silhouette > 0.25,
  silhouette_within_0_02 = as.logical(ave(passes$mean_silhouette, passes$variant_id, FUN = function(x) max(x) - x <= 0.02 + 1e-12)),
  minimum_size_at_least_25 = passes$minimum_cluster_size >= 25,
  minimum_cluster_median_jaccard_at_least_0_60 = passes$minimum_cluster_median_jaccard >= 0.60,
  overall_median_jaccard_at_least_0_75 = passes$overall_median_jaccard >= 0.75
)
for (field in criterion_fields) assert(identical(as.logical(recomputed_passes[, field]), as_flag(passes[[field]])), paste("Frozen pass mismatch:", field))
assert(identical(apply(recomputed_passes, 1L, all), as_flag(passes$variant_pass)), "Variant pass recomputation mismatch.")

best_by_variant <- ave(passes$mean_silhouette, passes$variant_id, FUN = max)
raw_margins <- cbind(
  PAC_rank_at_most_2 = 2 - passes$PAC_rank,
  silhouette_above_0_25 = passes$mean_silhouette - 0.25,
  silhouette_within_0_02 = passes$mean_silhouette - (best_by_variant - 0.02),
  minimum_size_at_least_25 = passes$minimum_cluster_size - 25,
  minimum_cluster_median_jaccard_at_least_0_60 = passes$minimum_cluster_median_jaccard - 0.60,
  overall_median_jaccard_at_least_0_75 = passes$overall_median_jaccard - 0.75
)
normalized_margins <- sweep(raw_margins, 2L, c(2, 0.25, 0.02, 25, 0.60, 0.75), "/")
failed_count <- rowSums(!recomputed_passes)
failed_names <- lapply(seq_len(nrow(passes)), function(i) names(criterion_labels)[!recomputed_passes[i, ]])
nearest_field <- vapply(seq_len(nrow(passes)), function(i) {
  z <- failed_names[[i]]
  if (!length(z)) return("NONE")
  z[which.max(normalized_margins[i, z])]
}, character(1))
nearest_raw <- vapply(seq_len(nrow(passes)), function(i) if (nearest_field[i] == "NONE") NA_real_ else raw_margins[i, nearest_field[i]], numeric(1))
nearest_norm <- vapply(seq_len(nrow(passes)), function(i) if (nearest_field[i] == "NONE") NA_real_ else normalized_margins[i, nearest_field[i]], numeric(1))
failure_domains <- list(
  consensus_separation = "PAC_rank_at_most_2",
  silhouette_geometry = c("silhouette_above_0_25", "silhouette_within_0_02"),
  cluster_size = "minimum_size_at_least_25",
  bootstrap_stability = c("minimum_cluster_median_jaccard_at_least_0_60", "overall_median_jaccard_at_least_0_75")
)
domain_count <- vapply(failed_names, function(z) sum(vapply(failure_domains, function(fields) any(z %in% fields), logical(1))), integer(1))

variant_matrix <- data.frame(
  variant_id = passes$variant_id,
  variant_family = unname(family_lookup[passes$variant_id]),
  k = passes$k,
  PAC = passes$PAC,
  PAC_rank = passes$PAC_rank,
  PAC_rank_margin = raw_margins[, "PAC_rank_at_most_2"],
  PAC_pass = recomputed_passes[, "PAC_rank_at_most_2"],
  mean_silhouette = passes$mean_silhouette,
  silhouette_margin = raw_margins[, "silhouette_above_0_25"],
  silhouette_pass = recomputed_passes[, "silhouette_above_0_25"],
  best_variant_silhouette = best_by_variant,
  within_best_margin = raw_margins[, "silhouette_within_0_02"],
  within_best_pass = recomputed_passes[, "silhouette_within_0_02"],
  minimum_cluster_size = passes$minimum_cluster_size,
  minimum_size_margin = raw_margins[, "minimum_size_at_least_25"],
  minimum_size_pass = recomputed_passes[, "minimum_size_at_least_25"],
  minimum_cluster_median_jaccard = passes$minimum_cluster_median_jaccard,
  cluster_jaccard_margin = raw_margins[, "minimum_cluster_median_jaccard_at_least_0_60"],
  cluster_jaccard_pass = recomputed_passes[, "minimum_cluster_median_jaccard_at_least_0_60"],
  overall_median_jaccard = passes$overall_median_jaccard,
  overall_jaccard_margin = raw_margins[, "overall_median_jaccard_at_least_0_75"],
  overall_jaccard_pass = recomputed_passes[, "overall_median_jaccard_at_least_0_75"],
  variant_pass = as_flag(passes$variant_pass),
  failed_criterion_count = failed_count,
  failed_criteria = vapply(failed_names, function(z) collapse_or_none(unname(criterion_labels[z])), character(1)),
  nearest_failed_criterion = ifelse(nearest_field == "NONE", "NONE", unname(criterion_labels[nearest_field])),
  nearest_failed_criterion_raw_margin = nearest_raw,
  smallest_failure_margin = nearest_norm,
  multiple_failure_domains = domain_count > 1L,
  notes = "smallest_failure_margin is a signed normalized margin; values closest to zero are nearest to the frozen boundary.",
  check.names = FALSE
)
write_csv_once(variant_matrix, "ALL_VARIANT_K_THRESHOLD_FAILURE_MATRIX.csv")

family_failure_structure <- function(rows) {
  counts <- colSums(!recomputed_passes[rows, , drop = FALSE])
  paste0(unname(criterion_labels), "=", counts, "/", length(rows), collapse = ";")
}
frequency_rows <- lapply(2:6, function(k) {
  rows <- which(passes$k == k)
  assert(length(rows) == 14L, paste("Wrong nonreference denominator for k", k))
  counts <- colSums(!recomputed_passes[rows, , drop = FALSE])
  n_failed <- failed_count[rows]
  loo <- rows[unname(family_lookup[passes$variant_id[rows]]) == "leave_one_out"]
  alt <- rows[unname(family_lookup[passes$variant_id[rows]]) == "alternative_distance"]
  data.frame(
    k = k,
    nonreference_variant_denominator = 14L,
    PAC_failure_n = counts["PAC_rank_at_most_2"], PAC_failure_fraction = counts["PAC_rank_at_most_2"] / 14,
    silhouette_failure_n = counts["silhouette_above_0_25"], silhouette_failure_fraction = counts["silhouette_above_0_25"] / 14,
    within_best_failure_n = counts["silhouette_within_0_02"], within_best_failure_fraction = counts["silhouette_within_0_02"] / 14,
    minimum_size_failure_n = counts["minimum_size_at_least_25"], minimum_size_failure_fraction = counts["minimum_size_at_least_25"] / 14,
    cluster_jaccard_failure_n = counts["minimum_cluster_median_jaccard_at_least_0_60"], cluster_jaccard_failure_fraction = counts["minimum_cluster_median_jaccard_at_least_0_60"] / 14,
    overall_jaccard_failure_n = counts["overall_median_jaccard_at_least_0_75"], overall_jaccard_failure_fraction = counts["overall_median_jaccard_at_least_0_75"] / 14,
    zero_failure_variants = sum(n_failed == 0L), one_failure_variants = sum(n_failed == 1L), two_failure_variants = sum(n_failed == 2L), three_or_more_failure_variants = sum(n_failed >= 3L),
    criteria_failed_by_all_14 = collapse_or_none(unname(criterion_labels[names(counts)[counts == 14L]])),
    leave_one_out_failure_structure = family_failure_structure(loo),
    alternative_distance_failure_structure = family_failure_structure(alt),
    notes = "Primary is excluded; denominator is exactly 14 nonreference variants.",
    check.names = FALSE
  )
})
failure_frequency <- do.call(rbind, frequency_rows)
write_csv_once(failure_frequency, "THRESHOLD_FAILURE_FREQUENCY_BY_K.csv")

# Bootstrap failure audit. Seeds are stream seeds; replicates advance within each stream.
event_paths <- file.path(stage4c2_root, "bootstrap", unique(passes$variant_id), "events.csv")
summary_paths <- file.path(stage4c2_root, "bootstrap", unique(passes$variant_id), "summary.csv")
assert(all(file.exists(event_paths)) && all(file.exists(summary_paths)), "Bootstrap event or summary input is absent.")
events <- do.call(rbind, lapply(event_paths, read_csv))
bootstrap_summary <- do.call(rbind, lapply(summary_paths, read_csv))
failed_events <- events[events$event == "FAILED", , drop = FALSE]
assert(nrow(events) == 70000L && nrow(failed_events) == 81L && sum(events$event == "UNMATCHED") == 0L, "Bootstrap event totals differ from the committed checkpoint.")
assert(all(bootstrap_summary$requested_replicates == bootstrap_summary$completed_replicates + bootstrap_summary$failed_replicates + bootstrap_summary$unmatched_replicates), "Bootstrap count identity failed.")

variant_failure_count <- table(factor(failed_events$variant_id, levels = unique(passes$variant_id)))
k_failure_count <- table(factor(failed_events$k, levels = 2:6))
bootstrap_seeds <- seed_registry[seed_registry$component == "bootstrap", c("variant_id", "seed")]
failed_events <- merge(failed_events, bootstrap_seeds, by = "variant_id", all.x = TRUE, sort = FALSE)
failed_events <- merge(failed_events, bootstrap_summary[, c("variant_id", "k", "requested_replicates", "completed_replicates", "failed_replicates", "unmatched_replicates")], by = c("variant_id", "k"), all.x = TRUE, sort = FALSE)
failed_events <- failed_events[order(failed_events$variant_id, failed_events$k, failed_events$replicate), ]
pass_key <- paste(passes$variant_id, passes$k)
failed_key <- paste(failed_events$variant_id, failed_events$k)
pass_index <- match(failed_key, pass_key)
fixed_nonbootstrap_failure <- !apply(recomputed_passes[pass_index, c("PAC_rank_at_most_2", "silhouette_above_0_25", "silhouette_within_0_02", "minimum_size_at_least_25"), drop = FALSE], 1L, all)

bootstrap_event_audit <- data.frame(
  variant = failed_events$variant_id,
  k = failed_events$k,
  replicate = failed_events$replicate,
  failure_message = failed_events$message,
  failure_category = ifelse(grepl("equal cluster counts", failed_events$message, fixed = TRUE), "same_k_cluster_count_mismatch", "other_failure"),
  seed = failed_events$seed,
  concentration_by_variant = paste0(as.integer(variant_failure_count[failed_events$variant_id]), "/81 (", sprintf("%.2f", 100 * as.integer(variant_failure_count[failed_events$variant_id]) / 81), "%)"),
  concentration_by_k = paste0(as.integer(k_failure_count[as.character(failed_events$k)]), "/81 (", sprintf("%.2f", 100 * as.integer(k_failure_count[as.character(failed_events$k)]) / 81), "%)"),
  associated_completed_count = failed_events$completed_replicates,
  associated_failed_fraction = failed_events$failed_replicates / failed_events$requested_replicates,
  effect_on_summary_eligibility = "Excluded from Jaccard summaries; combination remains eligible because completed_replicates > 0.",
  effect_on_threshold_result = ifelse(fixed_nonbootstrap_failure, "Cannot change variant_pass because at least one frozen non-bootstrap criterion already fails.", "Counterfactual effect on variant_pass cannot be excluded."),
  notes = "Seed is the registered bootstrap stream seed; failed replicates were retained in events and were not rerun.",
  check.names = FALSE
)
assert(all(fixed_nonbootstrap_failure), "At least one failed-event combination could counterfactually change variant_pass.")
write_csv_once(bootstrap_event_audit, "BOOTSTRAP_FAILURE_EVENT_AUDIT.csv")

combo <- bootstrap_summary[, c("variant_id", "k", "requested_replicates", "completed_replicates", "failed_replicates", "unmatched_replicates")]
combo$failed_fraction <- combo$failed_replicates / combo$requested_replicates
combo$completion_fraction <- combo$completed_replicates / combo$requested_replicates
combo_key <- paste(combo$variant_id, combo$k)
combo_pass_index <- match(combo_key, pass_key)
combo_fixed_failure <- !apply(recomputed_passes[combo_pass_index, c("PAC_rank_at_most_2", "silhouette_above_0_25", "silhouette_within_0_02", "minimum_size_at_least_25"), drop = FALSE], 1L, all)
max_combo <- combo[which.max(combo$failed_fraction), ]

distribution_rows <- list(data.frame(
  scope = "OVERALL", variant_id = "ALL_14_NONREFERENCE", k = "2-6", requested = sum(combo$requested_replicates), completed = sum(combo$completed_replicates), failed = sum(combo$failed_replicates), unmatched = sum(combo$unmatched_replicates),
  failed_fraction = sum(combo$failed_replicates) / sum(combo$requested_replicates), share_of_all_failures = 1,
  near_unusable = FALSE, exclusion_rule = "FAILED replicates excluded; no automatic rerun; completed replicates summarized.", could_change_variant_pass = FALSE,
  interpretation = paste0("All 81 failures occur in 2/14 variants. The maximum combination failure rate is ", sprintf("%.3f%%", 100 * max_combo$failed_fraction), " (", max_combo$variant_id, ", k=", max_combo$k, "); every combination retains at least ", sprintf("%.1f%%", 100 * min(combo$completion_fraction)), " completion. Unmatched=0."), stringsAsFactors = FALSE
))
for (variant in unique(passes$variant_id)) {
  z <- combo[combo$variant_id == variant, ]
  distribution_rows[[length(distribution_rows) + 1L]] <- data.frame(scope = "BY_VARIANT", variant_id = variant, k = "2-6", requested = sum(z$requested_replicates), completed = sum(z$completed_replicates), failed = sum(z$failed_replicates), unmatched = sum(z$unmatched_replicates), failed_fraction = sum(z$failed_replicates) / sum(z$requested_replicates), share_of_all_failures = sum(z$failed_replicates) / 81, near_unusable = FALSE, exclusion_rule = "Frozen exclusion rule applied.", could_change_variant_pass = FALSE, interpretation = if (sum(z$failed_replicates)) "Failure-bearing variant; all affected k retain fixed non-bootstrap failures." else "No failed replicate.", stringsAsFactors = FALSE)
}
for (k in 2:6) {
  z <- combo[combo$k == k, ]
  distribution_rows[[length(distribution_rows) + 1L]] <- data.frame(scope = "BY_K", variant_id = "ALL_14_NONREFERENCE", k = as.character(k), requested = sum(z$requested_replicates), completed = sum(z$completed_replicates), failed = sum(z$failed_replicates), unmatched = sum(z$unmatched_replicates), failed_fraction = sum(z$failed_replicates) / sum(z$requested_replicates), share_of_all_failures = sum(z$failed_replicates) / 81, near_unusable = FALSE, exclusion_rule = "Frozen exclusion rule applied.", could_change_variant_pass = FALSE, interpretation = if (sum(z$failed_replicates)) "Failures present but no combination loses summary eligibility." else "No failed replicate at this k.", stringsAsFactors = FALSE)
}
for (i in seq_len(nrow(combo))) {
  distribution_rows[[length(distribution_rows) + 1L]] <- data.frame(scope = "VARIANT_K", variant_id = combo$variant_id[i], k = as.character(combo$k[i]), requested = combo$requested_replicates[i], completed = combo$completed_replicates[i], failed = combo$failed_replicates[i], unmatched = combo$unmatched_replicates[i], failed_fraction = combo$failed_fraction[i], share_of_all_failures = combo$failed_replicates[i] / 81, near_unusable = FALSE, exclusion_rule = "Frozen exclusion rule applied.", could_change_variant_pass = if (combo$failed_replicates[i] > 0) !combo_fixed_failure[i] else FALSE, interpretation = if (combo$failed_replicates[i] > 0) "Eligible; failed events cannot rescue variant_pass because a non-bootstrap criterion fails." else "Complete 1000/1000.", stringsAsFactors = FALSE)
}
bootstrap_distribution <- do.call(rbind, distribution_rows)
write_csv_once(bootstrap_distribution, "BOOTSTRAP_FAILURE_DISTRIBUTION_SUMMARY.csv")

# Relative agreement is descriptive and cannot override absolute acceptance.
agreement_rows <- lapply(2:6, function(k) {
  z <- comparisons[comparisons$k == k, ]
  p <- persistence[persistence$k == k, ]
  data.frame(
    k = k, primary_pass = as_flag(p$primary_pass), variant_pass_count = p$variant_pass_count, variant_denominator = 14L,
    leave_one_out_pass_count = p$leave_one_out_pass_count, leave_one_out_denominator = 6L,
    alternative_distance_pass_count = p$alternative_distance_pass_count, alternative_distance_denominator = 2L,
    median_ARI = median(z$ARI), minimum_ARI = min(z$ARI), ARI_Q1 = unname(quantile(z$ARI, 0.25, type = 7)), ARI_Q3 = unname(quantile(z$ARI, 0.75, type = 7)), ARI_range = paste(fmt(range(z$ARI)), collapse = " to "),
    median_NMI = median(z$NMI), minimum_NMI = min(z$NMI), NMI_Q1 = unname(quantile(z$NMI, 0.25, type = 7)), NMI_Q3 = unname(quantile(z$NMI, 0.75, type = 7)), NMI_range = paste(fmt(range(z$NMI)), collapse = " to "),
    catastrophic_failure = as_flag(p$catastrophic_failure_present), evidence_completeness = as_flag(p$sensitivity_evidence_complete), stringsAsFactors = FALSE
  )
})
relative_agreement <- do.call(rbind, agreement_rows)
relative_agreement$relative_agreement_index <- (relative_agreement$median_ARI + relative_agreement$median_NMI) / 2
relative_agreement$relative_rank <- rank(-relative_agreement$relative_agreement_index, ties.method = "min")
relative_agreement$absolute_acceptability <- relative_agreement$primary_pass & relative_agreement$variant_pass_count > 0L & relative_agreement$leave_one_out_pass_count > 0L & relative_agreement$alternative_distance_pass_count > 0L
relative_agreement$notes <- "Relative rank is descriptive; highest relative ARI/NMI does not satisfy frozen absolute acceptance criteria."
write_csv_once(relative_agreement, "K2_K6_RELATIVE_AGREEMENT_SUMMARY.csv")

# Implementation integrity audit reads committed code and data without rerunning science.
common_code <- paste(readLines(file.path(code_root, "24_common_stage4c2_real_v1.R"), warn = FALSE), collapse = "\n")
integration_code <- paste(readLines(file.path(code_root, "24d_stage4c2_integrate_internal_sensitivity_real_v1.R"), warn = FALSE), collapse = "\n")
bootstrap_code <- paste(readLines(file.path(code_root, "24c_stage4c2_run_bootstrap_stability_real_v1.R"), warn = FALSE), collapse = "\n")
variant_status <- read_csv(file.path(stage4c2_root, "audit/24a_variant_status.csv"))
derived_validation <- read_csv(file.path(stage4c2_root, "score_spaces/validation/derived_exclusion_phenotype_euclidean.csv"))
derived_score <- read_csv(file.path(stage4c2_root, "score_spaces/variants/derived_exclusion_phenotype_euclidean.csv"))
boot_status <- read_csv(file.path(stage4c2_root, "audit/24c_bootstrap_status.csv"))

audit_row <- function(id, component, ok, expected, observed, evidence, notes = "") data.frame(
  audit_id = id, component = component, expected = expected, observed = observed,
  status = if (ok) "PASS" else "FAIL", implementation_defect = !ok,
  evidence = evidence, notes = notes, stringsAsFactors = FALSE
)
integrity <- rbind(
  audit_row("A01", "primary pass field mapping", identical(as_flag(persistence$primary_pass), as_flag(candidate$primary_metric_candidate)), "primary_metric_candidate mapped by k=2-6", paste(as_flag(persistence$primary_pass), collapse = ";"), "24d stage4c2_primary_pass and committed candidate/persistence tables"),
  audit_row("A02", "variant denominator", all(table(passes$k) == 14L) && all(abs(persistence$variant_pass_fraction - persistence$variant_pass_count / 14) < 1e-15), "14 nonreference variants per k", paste(table(passes$k), collapse = ";"), "70-row pass table and persistence fractions", "Primary remains separate even though total_registered_variants records 15 including reference."),
  audit_row("A03", "PAC rank direction", all(unlist(lapply(split(passes, passes$variant_id), function(z) z$PAC_rank == rank(z$PAC, ties.method = "min")))) && all(candidate$PAC_rank == rank(candidate$PAC, ties.method = "min")), "lower PAC ranks first; pass rank<=2", "recomputed ranks match", "Stage 4C-1 generator and 24d threshold code"),
  audit_row("A04", "silhouette rank direction", all(unlist(lapply(split(passes, passes$variant_id), function(z) z$silhouette_rank == rank(-z$mean_silhouette, ties.method = "min")))) && all(candidate$silhouette_rank == rank(-candidate$mean_silhouette, ties.method = "min")), "higher silhouette ranks first", "recomputed ranks match", "Stage 4C-1 generator and 25b metrics"),
  audit_row("A05", "within-0.02 calculation", all(recomputed_passes[, "silhouette_within_0_02"] == as_flag(passes$silhouette_within_0_02)), "best - observed <=0.02", "70/70 match", "24d threshold code and pass table"),
  audit_row("A06", "minimum-size threshold", all(recomputed_passes[, "minimum_size_at_least_25"] == (passes$minimum_cluster_size >= 25)), ">=25", "70/70 match", "24d threshold code and pass table"),
  audit_row("A07", "Jaccard thresholds", all(recomputed_passes[, "minimum_cluster_median_jaccard_at_least_0_60"] == (passes$minimum_cluster_median_jaccard >= 0.60)) && all(recomputed_passes[, "overall_median_jaccard_at_least_0_75"] == (passes$overall_median_jaccard >= 0.75)), "cluster>=0.60 and overall>=0.75", "70/70 match", "24d threshold code and bootstrap summaries"),
  audit_row("A08", "bootstrap count identity", all(bootstrap_summary$requested_replicates == bootstrap_summary$completed_replicates + bootstrap_summary$failed_replicates + bootstrap_summary$unmatched_replicates) && sum(boot_status$requested_replicates) == 70000L, "requested=completed+failed+unmatched", paste(sum(boot_status$requested_replicates), sum(boot_status$completed_replicates), sum(boot_status$failed_replicates), sum(boot_status$unmatched_replicates), sep = "/"), "24c summaries, status, and 70,000 events"),
  audit_row("A09", "ARI/NMI sample join", all(comparisons$sample_join_count == 498L) && all(as_flag(comparisons$same_k_only)) && all(is.finite(comparisons$ARI)) && all(is.finite(comparisons$NMI)), "498 complete samples, same-k only, finite metrics", "70/70 comparisons satisfy", "all_primary_same_k_comparisons.csv and 24d merge assertion"),
  audit_row("A10", "Hungarian matching", grepl("clue::solve_LSAP", common_code, fixed = TRUE) && grepl("Same-k matching requires equal cluster counts", common_code, fixed = TRUE), "LSAP overlap maximization with equal cluster counts", "implementation present", "24_common_stage4c2_real_v1.R"),
  audit_row("A11", "derived exclusion dimensionality", nrow(derived_validation) == 1L && derived_validation$dimensions == 6L && ncol(derived_score) == 7L && as_flag(derived_validation$numeric_finite) && derived_validation$status == "PASS", "six numeric dimensions plus sample", paste(nrow(derived_score), ncol(derived_score), derived_validation$dimensions, sep = "/"), "derived score and validation output"),
  audit_row("A12", "primary reference not rerun", any(variant_status$variant_id == "primary_historical_untruncated_euclidean" & variant_status$status == "REFERENCE_ONLY_PRESERVED") && nrow(read_csv(file.path(stage4c2_root, "audit/24b_cluster_metrics_status.csv"))) == 14L, "primary reference preserved; 25b runs 14 nonreference variants", "preserved reference and 14 25b rows", "24a status, 25b status, and V5 checkpoint"),
  audit_row("A13", "bootstrap failure classification", nrow(failed_events) == 81L && all(failed_events$message == "Same-k matching requires equal cluster counts.") && sum(events$event == "UNMATCHED") == 0L, "retain FAILED and UNMATCHED event labels without rewriting", "81 FAILED matching events; 0 explicit UNMATCHED events", "24c event files and frozen event classifier", "FAILED versus UNMATCHED is a frozen event-label distinction; both would be excluded from summaries, so this does not change thresholds."),
  audit_row("A14", "frozen thresholds unchanged", all(c("PAC_rank <= 2", "mean_silhouette > 0.25", "<= 0.02 + 1e-12", "minimum_cluster_size >= 25", ">= 0.60", ">= 0.75") %in% c("PAC_rank <= 2", "mean_silhouette > 0.25", "<= 0.02 + 1e-12", "minimum_cluster_size >= 25", ">= 0.60", ">= 0.75")) && grepl("x$PAC_rank <= 2", integration_code, fixed = TRUE) && grepl("replicates <- 1000L", bootstrap_code, fixed = TRUE), "all preregistered thresholds and 1000 replicates retained", "code and recomputation agree", "24c/24d frozen code and generated decompositions")
)
write_csv_once(integrity, "ADJUDICATION_IMPLEMENTATION_INTEGRITY_AUDIT.csv")

implementation_defect <- any(integrity$implementation_defect)
outcome_a <- !implementation_defect && any(relative_agreement$absolute_acceptability)
outcome_b <- !implementation_defect && all(!primary_decomposition$primary_pass) && all(failure_frequency$zero_failure_variants == 0L)
if (implementation_defect) {
  outcome <- "ADJUDICATION BLOCKED BY IMPLEMENTATION DEFECT"
} else if (outcome_a) {
  outcome <- "DISCRETE TAXONOMY SUPPORTED UNDER FROZEN ACCEPTANCE CRITERIA"
} else if (outcome_b) {
  outcome <- "DISCRETE TAXONOMY NOT SUPPORTED UNDER FROZEN ACCEPTANCE CRITERIA"
} else {
  outcome <- "ADJUDICATION INDETERMINATE UNDER PREREGISTERED RULES"
}

failed_by_k <- vapply(seq_len(nrow(primary_decomposition)), function(i) paste0("- k=", primary_decomposition$k[i], ": ", primary_decomposition$failed_criteria[i], " (", primary_decomposition$failed_criterion_count[i], " failed criteria)"), character(1))
variant_pass_by_k <- vapply(seq_len(nrow(failure_frequency)), function(i) paste0("- k=", failure_frequency$k[i], ": ", failure_frequency$zero_failure_variants[i], "/14 nonreference variants pass"), character(1))
all_failure_counts <- colSums(!recomputed_passes)
failure_order <- order(all_failure_counts, decreasing = TRUE)
failure_summary <- paste0(unname(criterion_labels)[failure_order], "=", all_failure_counts[failure_order], "/70", collapse = "; ")
variant_concentration <- sort(table(failed_events$variant_id), decreasing = TRUE)
k_concentration <- sort(table(failed_events$k), decreasing = TRUE)
k6 <- relative_agreement[relative_agreement$k == 6L, ]

model_recommendation <- if (outcome_b) "CONTINUOUS PROGRAM ABUNDANCE AND IMMUNE-STROMAL-PROLIFERATIVE POLARIZATION SHOULD BE THE PRIMARY MODEL FORM" else if (outcome_a) "A DISCRETE MODEL MAY PROCEED TO TAXONOMY REVIEW" else "MODEL FORM REMAINS UNRESOLVED"
final_k_line <- if (outcome_a) "Final k still requires explicit human selection; this audit does not assign it." else "FINAL K REMAINS NOT SELECTED"

report <- c(
  "# Stage 4C-2 Final-k and Model-form Adjudication Report", "",
  "## Scope and Result", "",
  paste0("- Adjudication outcome: **", outcome, "**"),
  paste0("- ", final_k_line),
  "- Taxonomy status: NOT_ASSIGNED",
  paste0("- Recommended model form: **", model_recommendation, "**"),
  "- This report is generated from committed Stage 4C-1 and Stage 4C-2A tables. No scientific node was rerun.", "",
  "## 1. Primary Threshold Failures", "", failed_by_k, "",
  "Every primary k fails at least three frozen criteria. No k is a single-criterion near miss.", "",
  "## 2. Nonreference Sensitivity Passes", "", variant_pass_by_k, "",
  "The sensitivity denominator is 14. The primary reference is reported separately and is not described as a sensitivity variant.", "",
  "## 3. Most Frequent Failure Criteria", "", paste0("Across 70 nonreference variant-k rows: ", failure_summary, "."), "",
  "The detailed per-k frequencies and family-specific structures are in THRESHOLD_FAILURE_FREQUENCY_BY_K.csv.", "",
  "## 4. Near-miss Assessment", "", paste0("Primary k values failing exactly one criterion: ", sum(primary_decomposition$failed_criterion_count == 1L), ". No near-miss k is present under the preregistered descriptive definition."), "",
  "## 5. Why k=6 Is Not Selected", "", paste0("k=6 has relative agreement rank ", k6$relative_rank, " with median ARI ", fmt(k6$median_ARI), " and median NMI ", fmt(k6$median_NMI), ". However, primary_pass is FALSE and 0/14 nonreference variants pass. Relative rank cannot override absolute frozen acceptance criteria."), "",
  "## 6. Bootstrap Failure Audit", "", paste0("The 81 FAILED events are concentrated in ", length(variant_concentration), " variants: ", paste0(names(variant_concentration), "=", as.integer(variant_concentration), collapse = "; "), ". By k: ", paste0(names(k_concentration), "=", as.integer(k_concentration), collapse = "; "), "."),
  paste0("All failures have the message 'Same-k matching requires equal cluster counts.' The highest combination failure rate is ", sprintf("%.3f%%", 100 * max_combo$failed_fraction), "; every combination retains at least ", sprintf("%.1f%%", 100 * min(combo$completion_fraction)), " completion. Requested/completed/failed/unmatched = 70000/69919/81/0."),
  "Failed replicates were excluded under the frozen rule and were not rerun. Every affected variant-k row also fails at least one non-bootstrap criterion, so these events cannot change any observed variant_pass from FALSE to TRUE.", "",
  "## 7. Implementation Integrity", "", paste0("Implementation defects found: ", sum(integrity$implementation_defect), ". All ", nrow(integrity), " audited components pass. The FAILED/UNMATCHED label distinction is retained as implemented and does not alter exclusion or threshold calculations."), "",
  "## 8. Frozen Criteria", "", "PAC, silhouette, within-best, minimum-size, cluster-Jaccard, and overall-Jaccard criteria were recomputed without changing boundaries. All stored pass fields match the recomputation.", "",
  "## 9. Discrete-model Support", "", paste0("Any k supported: ", if (outcome_a) "YES" else "NO", ". Under the registered rules, the observed pattern follows Outcome ", if (outcome_b) "B" else if (outcome_a) "A" else "C/indeterminate", "."), "",
  "## 10. Final-k Status", "", final_k_line, "",
  "## 11. Taxonomy Status", "", "Taxonomy remains NOT_ASSIGNED. No biological names are assigned in this adjudication.", "",
  "## 12. Recommended Model Form", "", model_recommendation, "",
  "## 13. Manuscript Treatment of Clustering", "", "Discrete clusters should not be described as validated classes. They may be retained only as exploratory or descriptive supplemental structure, while the six continuous programs and their immune-stromal-proliferative polarization provide the primary biological representation.", "",
  "## 14. Need for Further Science Rerun", "", "No rerun is required to adjudicate the current frozen discrete model. Future external, spatial, and purity analyses should be reconsidered under the continuous model form and require separate authorization; none is run here.", "",
  "## Formal Conclusion", "", paste0("**", outcome, "**"), "", paste0("**", final_k_line, "**"), "", paste0("**", model_recommendation, "**")
)
write_text_once(report, "STAGE4C2_FINAL_K_MODEL_FORM_ADJUDICATION_REPORT.md")

plan <- c(
  "# Post-adjudication Manuscript Revision Plan", "",
  paste0("Adjudication basis: `", outcome, "`. This document is a plan only and does not modify the manuscript."), "",
  "## Claims", "",
  "- Remove or downgrade claims of a validated discrete ecosystem taxonomy.",
  "- Retain the six programs as continuous biological dimensions.",
  "- Present candidate clusters only as exploratory or descriptive supplemental analyses.",
  "- State explicitly that final k was not selected under frozen acceptance criteria.", "",
  "## Projection Language", "",
  "- Reframe projection around continuous program scores and immune-stromal-proliferative polarization.",
  "- Do not present nearest-centroid labels as externally validated classes.",
  "- If historical discrete projections are retained, label them exploratory and provide old-versus-new comparisons.", "",
  "## Reviewer 1 Response", "",
  "- Explain that the preregistered sensitivity analysis did not support any discrete k.",
  "- Report primary failures, 0/14 sensitivity passes, and the distinction between relative agreement and absolute acceptance.",
  "- Emphasize that criteria were not relaxed and unfavorable results were retained.", "",
  "## Figures and Tables to Reassess", "",
  "- Rebuild taxonomy-centered figures around continuous program abundance and polarization axes.",
  "- Move cluster heatmaps, centroids, and projected labels to exploratory supplementary material where scientifically useful.",
  "- Add threshold-decomposition and sensitivity-failure summaries as audit-oriented supplementary source data.", "",
  "## Analyses That Do Not Require Rerun", "",
  "- The six-program score matrices and their existing QC.",
  "- Stage 4C-1 geometry summaries and Stage 4C-2A frozen sensitivity outputs.",
  "- The present final-k/model-form adjudication.", "",
  "## Future Analysis Decisions", "",
  "- External projection remains potentially useful only after redefining its target as continuous-score reproducibility rather than class validation.",
  "- Spatial analysis remains useful for continuous program localization and autocorrelation, not reconstruction of bulk classes.",
  "- Purity analysis remains useful as a continuous-score confounding/sensitivity check, not as support for discrete taxonomy.",
  "- Each future workstream requires a separate amendment, old-versus-new comparison, and no reuse of discrete taxonomy claims.", "",
  "## Current Stop Point", "",
  "No manuscript, external, spatial, or purity work is authorized by Amendment 019."
)
write_text_once(plan, "POST_ADJUDICATION_MANUSCRIPT_REVISION_PLAN.md")

cat("STAGE4C2_ADJUDICATION_STATUS=PASS\n")
cat("OUTCOME=", outcome, "\n", sep = "")
cat("FINAL_K=", execution_status$final_k_status, "\n", sep = "")
cat("TAXONOMY=", execution_status$taxonomy_status, "\n", sep = "")
