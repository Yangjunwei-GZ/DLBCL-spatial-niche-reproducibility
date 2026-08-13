DLBCL_PROJECT_ROOT <- normalizePath(
  Sys.getenv("DLBCL_PROJECT_ROOT", unset = getwd()),
  winslash = "/", mustWork = FALSE
)
DLBCL_SUPPLEMENTARY_CODE_ROOT <- normalizePath(
  Sys.getenv("DLBCL_SUPPLEMENTARY_CODE_ROOT",
    unset = file.path(DLBCL_PROJECT_ROOT, "Supplementary_Code")),
  winslash = "/", mustWork = FALSE
)

source(file.path(DLBCL_SUPPLEMENTARY_CODE_ROOT, "source_grounded_sensitivity", "01_common.R"))
o6_require_token(); o6_set_library(); o6_require_packages("digest")

read_o <- function(path) utils::read.csv(file.path(O6_OUTPUTS, path), check.names = FALSE, stringsAsFactors = FALSE)
input_hash <- read_o("run_control/06O_INPUT_HASH_REGISTRY.csv")
script_manifest <- read_o("run_control/06O_SCRIPT_MANIFEST.csv")
seed <- read_o("run_control/06O_SEED_REGISTRY.csv")
matching <- read_o("bulk/PAIRED_FULL_AND_CORE_GENE_MATCHING.csv")
scores <- read_o("bulk/PAIRED_FULL_AND_CORE_SCORES.csv")
historical <- read_o("bulk/FULL_RESCORE_VS_HISTORICAL_AUTHORITY.csv")
bulk_summary <- read_o("bulk/FULL_CORE_SCORE_SUMMARY.csv")
geometry <- read_o("bulk/FULL_CORE_PROGRAM_CORRELATION_COMPARISON.csv")
gate <- read_o("spatial/FULL_REGENERATION_AUTHORITY_CHECK.csv")
detection <- read_o("spatial/CORE_GENE_DETECTION_BY_AREA.csv")
concordance <- read_o("spatial/FULL_CORE_SPOT_SCORE_CONCORDANCE.csv")
spatial <- read_o("spatial/CORE_MORAN_GEARY.csv")
permutations <- read_o("spatial/CORE_SCORE_LABEL_PERMUTATIONS.csv")
fdr <- read_o("spatial/06O_SENSITIVITY_FDR.csv")
protected_baseline <- read_o("run_control/06O_PROTECTED_PATH_HASH_BASELINE.csv")

checks <- list()
add <- function(id, check, pass, observed, expected, mandatory = TRUE) {
  checks[[length(checks) + 1L]] <<- data.frame(check_id = id, check = check,
    status = if (isTRUE(pass)) "PASS" else "FAIL", observed = as.character(observed),
    expected = as.character(expected), mandatory = mandatory, stringsAsFactors = FALSE)
}

add("06O-V001", "Exact token", identical(readChar(file.path(O6_ROOT,
  "00_execution_authorization/06O_EXECUTION_TOKEN.txt"), nchar(O6_TOKEN), useBytes = TRUE), O6_TOKEN), "validated", O6_TOKEN)
add("06O-V002", "Preflight input hashes", all(input_hash$verification_status == "PASS"),
  sum(input_hash$verification_status == "PASS"), nrow(input_hash))

file_hash_rows <- input_hash[input_hash$verification_method == "FILE_SHA256", , drop = FALSE]
rehash_ok <- vapply(seq_len(nrow(file_hash_rows)), function(i) file.exists(file_hash_rows$path[i]) &&
  identical(o6_sha256(file_hash_rows$path[i]), file_hash_rows$actual_sha256[i]) &&
  identical(as.numeric(file.info(file_hash_rows$path[i])$size), as.numeric(file_hash_rows$actual_size_bytes[i])), logical(1))
add("06O-V003", "All file inputs unchanged after execution", all(rehash_ok), sum(rehash_ok), length(rehash_ok))

runtime_path <- file.path(O6_OUTPUTS, "run_control/06O_RUNTIME_AND_PACKAGE_CAPTURE.txt")
runtime_text <- paste(readLines(runtime_path, warn = FALSE), collapse = "\n")
add("06O-V004", "Runtime capture and GSVA 2.4.9", file.exists(runtime_path) &&
  grepl("GSVA\\s+2.4.9", runtime_text) && grepl("BiocParallel::SerialParam", runtime_text),
  "runtime captured", "GSVA 2.4.9; SerialParam")

script_now <- sort(list.files(O6_SCRIPTS, pattern = "\\.(R|ps1)$", full.names = TRUE))
script_ok <- setequal(normalizePath(script_now, winslash = "/"), script_manifest$absolute_path) &&
  all(vapply(script_manifest$absolute_path, o6_sha256, character(1)) == script_manifest$sha256)
add("06O-V005", "Exact script hashes", script_ok, nrow(script_manifest), length(script_now))

sets <- o6_gene_sets()
add("06O-V006", "Canonical and core memberships", nrow(sets$authority) == 132L &&
  identical(lengths(sets$core), c(macrophage_rich = 17L, t_cell_inflamed = 15L,
    antigen_presentation = 8L, stromal_fibrotic = 16L, proliferative_cycling = 9L)),
  paste(nrow(sets$authority), paste(lengths(sets$core), collapse = ";"), sep = "/"), "132/17;15;8;16;9")
add("06O-V007", "No immune-cold direct-source core", !any(matching$set_type == "CORE" &
  matching$program_id == "immune_cold_exclusion") && all(detection$core_gene_count[detection$program_id ==
  "immune_cold_exclusion"] == 0L), "excluded", "NOT_EVALUABLE")

sample_expected <- c(GSE31312 = 498L, GSE10846 = 420L, GSE181063 = 1310L)
sample_observed <- vapply(names(sample_expected), function(d) length(unique(scores$sample_id[scores$dataset_id == d])), integer(1))
paired_ok <- all(sample_observed == sample_expected) && all(vapply(names(sample_expected), function(d) {
  f <- scores$sample_id[scores$dataset_id == d & scores$set_type == "FULL"]
  c <- scores$sample_id[scores$dataset_id == d & scores$set_type == "CORE"]
  setequal(f, c)
}, logical(1)))
add("06O-V008", "Bulk paired sample identity", paired_ok,
  paste(names(sample_observed), sample_observed, sep = "=", collapse = ";"), "498;420;1310 and paired")
add("06O-V009", "Bulk outputs complete and finite", nrow(matching) == 33L && nrow(bulk_summary) == 15L &&
  nrow(historical) == 12L && nrow(geometry) == 66L && all(is.finite(scores$score)),
  paste(nrow(matching), nrow(bulk_summary), nrow(historical), nrow(geometry), sep = ";"), "33;15;12;66")

add("06O-V010", "All nine full UCell authority gates", nrow(gate) == 9L &&
  all(gate$authority_gate_status == "PASS") && all(gate$maximum_absolute_difference <= 1e-12),
  paste(table(gate$authority_gate_status), collapse = ";"), "9 PASS at <=1e-12")
add("06O-V011", "Spatial evaluable combinations", nrow(spatial) == 39L &&
  sum(detection$core_eligibility == "EVALUABLE") == 39L && sum(concordance$evaluation_status == "EVALUABLE") == 39L,
  paste(nrow(spatial), sum(detection$core_eligibility == "EVALUABLE"),
    sum(concordance$evaluation_status == "EVALUABLE"), sep = ";"), "39;39;39")
add("06O-V012", "Spatial NOT_EVALUABLE combinations", sum(detection$core_eligibility == "NOT_EVALUABLE") == 15L &&
  sum(detection$program_id == "antigen_presentation" & detection$core_eligibility == "NOT_EVALUABLE") == 6L &&
  sum(detection$program_id == "immune_cold_exclusion" & detection$core_eligibility == "NOT_EVALUABLE") == 9L,
  paste(sum(detection$core_eligibility == "NOT_EVALUABLE"),
    sum(detection$program_id == "antigen_presentation" & detection$core_eligibility == "NOT_EVALUABLE"),
    sum(detection$program_id == "immune_cold_exclusion" & detection$core_eligibility == "NOT_EVALUABLE"), sep = ";"), "15;6;9")

seed_key <- paste(seed$capture_area_id, seed$program, sub("_PERMUTATION_ENDPOINT$", "", seed$analysis_type))
perm_key <- paste(permutations$capture_area_id, permutations$program, permutations$endpoint)
seed_match <- match(perm_key, seed_key)
add("06O-V013", "Exact 78 endpoint seeds", nrow(permutations) == 78L && !anyNA(seed_match) &&
  all(permutations$derived_seed == seed$derived_seed[seed_match]) && length(unique(permutations$derived_seed)) == 78L,
  paste(nrow(permutations), length(unique(permutations$derived_seed)), sep = ";"), "78;78 exact")
add("06O-V014", "Exactly 9999 completed permutations per endpoint", all(permutations$requested_permutations == 9999L &
  permutations$completed_permutations == 9999L & permutations$all_null_statistics_finite),
  paste(range(permutations$completed_permutations), collapse = ";"), "9999;9999 and finite")

family_summary <- fdr[fdr$row_type == "FAMILY_SUMMARY", , drop = FALSE]
active_tests <- fdr[fdr$row_type == "TEST", , drop = FALSE]
add("06O-V015", "Six separate FDR families", nrow(family_summary) == 6L &&
  sum(active_tests$endpoint == "CORE_MORAN") == 39L && sum(active_tests$endpoint == "CORE_GEARY") == 39L &&
  all(family_summary$observed_test_count == family_summary$expected_test_count),
  paste(nrow(family_summary), sum(active_tests$endpoint == "CORE_MORAN"),
    sum(active_tests$endpoint == "CORE_GEARY"), sep = ";"), "6;39;39")

current_protected <- do.call(rbind, lapply(unique(protected_baseline$protected_root), function(root_name) {
  root <- file.path(O6_REVISION, root_name)
  files <- sort(list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE))
  files <- files[!file.info(files)$isdir]
  data.frame(protected_root = root_name, path = normalizePath(files, winslash = "/"),
    size_bytes = file.info(files)$size, sha256 = vapply(files, o6_sha256, character(1)), stringsAsFactors = FALSE)
}))
base_key <- paste(protected_baseline$protected_root, tolower(protected_baseline$path))
current_key <- paste(current_protected$protected_root, tolower(current_protected$path))
ord <- match(base_key, current_key)
protected_ok <- nrow(current_protected) == nrow(protected_baseline) && !anyNA(ord) &&
  all(as.numeric(protected_baseline$size_bytes) == current_protected$size_bytes[ord]) &&
  all(protected_baseline$sha256 == current_protected$sha256[ord])
add("06O-V016", "Protected trees unchanged including 06p", protected_ok,
  paste(nrow(protected_baseline), nrow(current_protected), sep = ";"), "identical file sets, sizes, SHA-256")
add("06O-V017", "06p not executed", !dir.exists(file.path(O6_REVISION,
  "06p_wp3_matched_null_and_depth_sensitivity/01_execution_outputs")), "no execution directory", "06p unauthorized")

tracked <- system2("git", c("-C", shQuote(O6_PROJECT), "diff", "--name-only"), stdout = TRUE)
staged <- system2("git", c("-C", shQuote(O6_PROJECT), "diff", "--cached", "--name-only"), stdout = TRUE)
add("06O-V018", "No tracked or staged Git changes", !length(tracked) && !length(staged),
  paste(length(tracked), length(staged), sep = ";"), "0;0")

csv_files <- list.files(O6_OUTPUTS, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
parse_ok <- vapply(csv_files, function(p) !inherits(try(utils::read.csv(p, check.names = FALSE), silent = TRUE), "try-error"), logical(1))
add("06O-V019", "All current CSV outputs parse", all(parse_ok), sum(parse_ok), length(parse_ok))
forbidden <- c("cluster", "selected_k", "taxonomy", "centroid", "survival", "pathway", "ligand_receptor",
  "patient_assignment", "specimen_assignment", "region_assignment")
output_columns <- unique(tolower(unlist(lapply(csv_files, function(p) names(utils::read.csv(p, nrows = 1L, check.names = FALSE))))))
add("06O-V020", "No prohibited analytical fields", !any(output_columns %in% forbidden),
  paste(intersect(output_columns, forbidden), collapse = ";"), "none")
add("06O-V021", "Negative and heterogeneous results retained", nrow(scores) > 0L && nrow(spatial) == 39L &&
  !anyNA(bulk_summary$full_core_Spearman), "all complete rows retained", "no selective omission")

validation <- do.call(rbind, checks)
mandatory_pass <- all(validation$status[validation$mandatory] == "PASS")
o6_write_csv_once(validation, file.path(O6_OUTPUTS, "validation/06O_EXECUTION_VALIDATION.csv"))

fmt <- function(x, digits = 3L) ifelse(is.finite(x), formatC(x, digits = digits, format = "f"), "NA")
bulk_lines <- unlist(lapply(split(bulk_summary, bulk_summary$dataset_id), function(x) c(
  paste0("### ", unique(x$dataset_id)),
  paste0("- ", x$program, ": Pearson ", fmt(x$full_core_Pearson),
    "; Spearman ", fmt(x$full_core_Spearman), "; full/core genes present ",
    x$full_genes_present, "/", x$core_genes_present)
)))
gate_lines <- paste0("- ", gate$capture_area_id, ": ", gate$authority_gate_status,
  "; max absolute difference = ", format(gate$maximum_absolute_difference, scientific = TRUE))
spatial_summary <- aggregate(cbind(Spearman, Pearson) ~ area_role + program, data = concordance[concordance$evaluation_status == "EVALUABLE", ],
  FUN = function(x) stats::median(x, na.rm = TRUE))
spatial_lines <- paste0("- ", spatial_summary$area_role, " / ", spatial_summary$program,
  ": median Pearson ", fmt(spatial_summary$Pearson), "; median Spearman ", fmt(spatial_summary$Spearman))

bulk_min <- min(bulk_summary$full_core_Spearman, na.rm = TRUE)
spatial_min <- min(concordance$Spearman[concordance$evaluation_status == "EVALUABLE"], na.rm = TRUE)
retained_direction <- mean(c(spatial$Moran_direction_concordant, spatial$Geary_departure_direction_concordant))
interpretation <- if (bulk_min >= 0.8 && spatial_min >= 0.8 && retained_direction >= 0.9) {
  "ROBUST_WITHIN_EVALUABLE_DIRECT_SOURCE_CORES"
} else if (stats::median(bulk_summary$full_core_Spearman) >= 0.5 &&
           stats::median(concordance$Spearman[concordance$evaluation_status == "EVALUABLE"]) >= 0.5) {
  "PARTIAL_ROBUSTNESS_WITH_HETEROGENEITY"
} else "MATERIAL_SENSITIVITY_TO_AUTHOR_CURATED_EXTENSIONS"

report <- c(
  "# 06o Source-Grounded Program Sensitivity: Final Report",
  "",
  paste0("- execution_status: ", if (mandatory_pass) "PASSED_FINAL_VALIDATION" else "FAILED_FINAL_VALIDATION"),
  paste0("- authority_status: ", if (mandatory_pass) "FINAL_06O_AUTHORITY" else "NOT_AUTHORITY"),
  paste0("- interpretation_class: ", interpretation),
  "- final_k: NOT_SELECTED",
  "- taxonomy: NOT_ASSIGNED",
  "- 06p executed: FALSE",
  "",
  "## Bulk Paired Rescoring",
  "All three bulk datasets were rescored as paired full/core analyses from the same dataset-specific preprocessed expression matrix using GSVA 2.4.9, ssgseaParam(minSize=5, maxSize=500, alpha=0.25, normalize=TRUE), and SerialParam.",
  bulk_lines,
  "",
  "The newly generated GSE31312 full scores were compared descriptively with both the historical raw score authority and the final within-program z-score authority. The complete per-program correlations, absolute differences, RMSE, maximum differences, and fitted scaling relationships are retained in FULL_RESCORE_VS_HISTORICAL_AUTHORITY.csv. Equality was not required because the May 2026 historical GSVA version is unknown.",
  "",
  "Five-program full/core correlation geometry, including all ten pairwise associations, core-minus-full differences, rank concordance, and direction concordance, is retained without filtering in FULL_CORE_PROGRAM_CORRELATION_COMPARISON.csv.",
  "",
  "## Spatial Full-Regeneration Authority Gate",
  gate_lines,
  "",
  paste0("Evaluable core program-area combinations: ", nrow(spatial), ". NOT_EVALUABLE combinations: ",
    sum(detection$core_eligibility == "NOT_EVALUABLE"), "."),
  "Immune-cold/exclusion has no direct-source core and is NOT_EVALUABLE in all nine capture areas. Antigen-presentation is NOT_EVALUABLE in Cap.area4-Cap.area9 because fewer than five direct-source core genes are detected; no substitution or imputation was performed.",
  "",
  "## Full-versus-Core Spot Scores",
  spatial_lines,
  "",
  "## Core Moran and Geary",
  paste0("All ", nrow(spatial), " evaluable combinations completed Moran and Geary statistics with exactly 9,999 score-label permutations per endpoint. Core-minus-full differences, effect-retention ratios, NOT_APPLICABLE denominators, and direction concordance are reported in CORE_MORAN_GEARY.csv."),
  paste0("Moran BH families: PRIMARY_DLBCL=", sum(fdr$row_type == "TEST" & fdr$endpoint == "CORE_MORAN" & fdr$role_family == "PRIMARY_DLBCL"),
    ", CONTEXT_ONLY=", sum(fdr$row_type == "TEST" & fdr$endpoint == "CORE_MORAN" & fdr$role_family == "CONTEXT_ONLY"),
    ", EXPLORATORY_ANTIGEN=0 (NOT_RUN)."),
  paste0("Geary BH families: PRIMARY_DLBCL=", sum(fdr$row_type == "TEST" & fdr$endpoint == "CORE_GEARY" & fdr$role_family == "PRIMARY_DLBCL"),
    ", CONTEXT_ONLY=", sum(fdr$row_type == "TEST" & fdr$endpoint == "CORE_GEARY" & fdr$role_family == "CONTEXT_ONLY"),
    ", EXPLORATORY_ANTIGEN=0 (NOT_RUN)."),
  "",
  "## Interpretation and Limitations",
  paste0("The complete pattern is classified as ", interpretation, ". This is a sensitivity assessment of dependence on author-curated extensions, not a new taxonomy or a validation of discrete classes."),
  "The analyses remain capture-area and spot-level. They do not establish patient-level replication, prognosis, pathway activity, or cell-cell communication. Direct-source cores are smaller than the canonical programs, low spatial gene detection limits inference, and the historical bulk comparison is descriptive because its package version is unknown. Negative and heterogeneous findings are retained.",
  "",
  paste0("Outputs may be used in manuscript/reviewer-response drafting only if authority_status is FINAL_06O_AUTHORITY. Current status: ",
    if (mandatory_pass) "ELIGIBLE_FOR_CAUTIOUS_USE" else "NOT_ELIGIBLE"),
  ""
)
o6_write_text_once(report, file.path(O6_OUTPUTS, "06O_FINAL_REPORT.md"))

decision <- c("# 06o Authority Decision", "",
  paste0("- validation: ", if (mandatory_pass) "PASS" else "FAIL"),
  paste0("- authority: ", if (mandatory_pass) "FINAL_06O_AUTHORITY_ASSIGNED" else "06O_AUTHORITY_NOT_ASSIGNED"),
  "- scope: paired full/core bulk and spatial direct-source-core sensitivity only",
  "- final_k: NOT_SELECTED", "- taxonomy: NOT_ASSIGNED", "- 06p: NOT_EXECUTED",
  "- manuscript_use: cautious sensitivity reporting only; no discrete class, prognostic, pathway, LR, or patient-level claim")
o6_write_text_once(decision, file.path(O6_OUTPUTS, "validation/06O_AUTHORITY_DECISION.md"))

manifest_files <- sort(list.files(O6_OUTPUTS, recursive = TRUE, full.names = TRUE))
manifest_files <- manifest_files[basename(manifest_files) != "06O_OUTPUT_MANIFEST.csv"]
manifest <- data.frame(relative_path = substring(normalizePath(manifest_files, winslash = "/"),
    nchar(normalizePath(O6_OUTPUTS, winslash = "/")) + 2L),
  absolute_path = normalizePath(manifest_files, winslash = "/"),
  size_bytes = file.info(manifest_files)$size,
  sha256 = vapply(manifest_files, o6_sha256, character(1)),
  manifest_self_excluded = TRUE, stringsAsFactors = FALSE)
o6_write_csv_once(manifest, file.path(O6_OUTPUTS, "validation/06O_OUTPUT_MANIFEST.csv"))

if (!mandatory_pass) stop("06o final validator failed; authority not assigned")
cat("FINAL_VALIDATION=PASS\nAUTHORITY=FINAL_06O_AUTHORITY_ASSIGNED\n")
