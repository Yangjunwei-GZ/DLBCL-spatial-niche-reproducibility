DLBCL_PROJECT_ROOT <- normalizePath(
  Sys.getenv("DLBCL_PROJECT_ROOT", unset = getwd()),
  winslash = "/", mustWork = FALSE
)
DLBCL_SUPPLEMENTARY_CODE_ROOT <- normalizePath(
  Sys.getenv("DLBCL_SUPPLEMENTARY_CODE_ROOT",
    unset = file.path(DLBCL_PROJECT_ROOT, "Supplementary_Code")),
  winslash = "/", mustWork = FALSE
)

# WP2 real external continuous validation orchestrator.
# This is the only authorized entry point for real GSE10846/GSE181063 execution.

source(file.path(DLBCL_SUPPLEMENTARY_CODE_ROOT, "external_validation", "29_common_wp2_functions_v1.R"))
.libPaths(c(WP2$renv_library, .Library))
wp2_assert_r451()
wp2_require_packages()
wp2_assert_token()
on.exit(Sys.unsetenv(WP2$token_name), add = TRUE)

if (!identical(trimws(system("git branch --show-current", intern = TRUE)), WP2$branch)) stop("Wrong WP2 execution branch", call. = FALSE)
if (!identical(trimws(system("git rev-parse HEAD", intern = TRUE)), WP2$baseline)) stop("Wrong WP2 baseline HEAD", call. = FALSE)
if (wp2_other_r_process_count() != 0L) stop("Another R/Rscript process is active", call. = FALSE)
wp2_verify_input_registry()
wp2_verify_protected_paths()

output_manifest <- utils::read.csv(file.path(WP2$stage, "WP2_REAL_OUTPUT_MANIFEST.csv"), check.names = FALSE, stringsAsFactors = FALSE)
preexisting_real <- vapply(output_manifest$relative_path, function(path) file.exists(file.path(WP2$stage, path)), logical(1))
if (any(preexisting_real)) stop("Create-once WP2 real output already exists: ", paste(output_manifest$relative_path[preexisting_real], collapse = ";"), call. = FALSE)

node_scripts <- c(
  `29a` = "29a_wp2_parse_and_map_external_expression_v1.R",
  `29b` = "29b_wp2_gene_mapping_and_probe_collapse_v1.R",
  `29c` = "29c_wp2_external_program_scoring_v1.R",
  `29d` = "29d_wp2_score_qc_and_mapping_sensitivity_v1.R",
  `29e` = "29e_wp2_program_structure_validation_v1.R",
  `29f` = "29f_wp2_strategy_b_frozen_loading_projection_v1.R",
  `29g` = "29g_wp2_strategy_a_scale_diagnostics_and_projection_v1.R",
  `29h` = "29h_wp2_cross_cohort_integration_v1.R"
)
node_functions <- c(
  `29a` = "wp2_node_29a", `29b` = "wp2_node_29b", `29c` = "wp2_node_29c", `29d` = "wp2_node_29d",
  `29e` = "wp2_node_29e", `29f` = "wp2_node_29f", `29g` = "wp2_node_29g", `29h` = "wp2_node_29h"
)

old_source_only <- Sys.getenv("WP2_NODE_SOURCE_ONLY", unset = "")
Sys.setenv(WP2_NODE_SOURCE_ONLY = "TRUE")
on.exit({
  if (nzchar(old_source_only)) Sys.setenv(WP2_NODE_SOURCE_ONLY = old_source_only) else Sys.unsetenv("WP2_NODE_SOURCE_ONLY")
}, add = TRUE)
for (script in node_scripts) source(file.path(DLBCL_SUPPLEMENTARY_CODE_ROOT, "external_validation", script), local = .GlobalEnv)

runtime_rows <- list()
node_rows <- list()
for (node in names(node_functions)) {
  warnings <- character()
  started <- Sys.time()
  status <- "COMPLETED"
  error_message <- ""
  stack_before <- sys.nframe()
  tryCatch(
    withCallingHandlers(
      do.call(get(node_functions[[node]], envir = .GlobalEnv), list()),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      status <<- "FAILED"
      error_message <<- conditionMessage(e)
    }
  )
  finished <- Sys.time()
  stack_after <- sys.nframe()
  runtime_rows[[node]] <- data.frame(
    node = node,
    started_at = format(started, "%Y-%m-%dT%H:%M:%S%z"),
    finished_at = format(finished, "%Y-%m-%dT%H:%M:%S%z"),
    elapsed_seconds = as.numeric(difftime(finished, started, units = "secs")),
    status = status,
    warning_count = length(warnings),
    warning_messages = paste(warnings, collapse = " | "),
    stack_before = stack_before,
    stack_after = stack_after,
    stack_imbalance = stack_before != stack_after,
    error_message = error_message,
    stringsAsFactors = FALSE
  )
  for (dataset in c("GSE10846", "GSE181063")) {
    node_rows[[paste(node, dataset, sep = "_")]] <- data.frame(
      node = node,
      dataset_id = dataset,
      status = status,
      failure_isolated = ifelse(status == "FAILED", "COHORT_OUTPUTS_RETAINED_IF_ALREADY_CREATED", "NOT_APPLICABLE"),
      error_message = error_message,
      stringsAsFactors = FALSE
    )
  }
  if (status == "FAILED" || length(warnings) > 0L || stack_before != stack_after) break
}

completed_nodes <- names(runtime_rows)
missing_nodes <- setdiff(names(node_functions), completed_nodes)
if (length(missing_nodes)) {
  for (node in missing_nodes) {
    runtime_rows[[node]] <- data.frame(
      node = node, started_at = "", finished_at = "", elapsed_seconds = NA_real_, status = "NOT_RUN_AFTER_CONTROLLED_FAILURE",
      warning_count = 0L, warning_messages = "", stack_before = NA_integer_, stack_after = NA_integer_, stack_imbalance = FALSE,
      error_message = "Prior node failed", stringsAsFactors = FALSE
    )
    for (dataset in c("GSE10846", "GSE181063")) {
      node_rows[[paste(node, dataset, sep = "_")]] <- data.frame(
        node = node, dataset_id = dataset, status = "NOT_RUN_AFTER_CONTROLLED_FAILURE", failure_isolated = "TRUE",
        error_message = "Prior node failed", stringsAsFactors = FALSE
      )
    }
  }
}

runtime <- do.call(rbind, runtime_rows)
nodes <- do.call(rbind, node_rows)
wp2_write_csv(runtime, file.path(WP2$stage, "logs/WP2_RUNTIME_AUDIT.csv"))
wp2_write_csv(nodes, file.path(WP2$stage, "checkpoints/WP2_NODE_STATUS.csv"))

if (any(runtime$status != "COMPLETED") || any(runtime$warning_count != 0L) || any(runtime$stack_imbalance)) {
  stop("WP2 orchestrator stopped after a controlled node failure", call. = FALSE)
}

# The final validator must inspect the required post-run UNSET token state.
Sys.unsetenv(WP2$token_name)
source(file.path(DLBCL_SUPPLEMENTARY_CODE_ROOT, "external_validation", "29v_wp2_external_continuous_validation_validator_v1.R"), local = .GlobalEnv)
wp2_node_29v()
cat("WP2 EXTERNAL CONTINUOUS VALIDATION COMPLETED\n")
