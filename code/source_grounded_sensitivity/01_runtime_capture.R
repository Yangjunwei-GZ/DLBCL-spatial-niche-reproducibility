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
o6_require_token(); o6_set_library()

packages <- c("GSVA", "BiocParallel", "UCell", "Seurat", "SeuratObject", "Matrix",
  "BiocManager", "renv", "digest", "Biobase", "AnnotationDbi")
o6_require_packages(packages)
invisible(lapply(packages, requireNamespace, quietly = TRUE))

script_files <- sort(list.files(O6_SCRIPTS, pattern = "\\.(R|ps1)$", full.names = TRUE))
manifest <- data.frame(
  script_name = basename(script_files),
  absolute_path = normalizePath(script_files, winslash = "/"),
  size_bytes = file.info(script_files)$size,
  sha256 = vapply(script_files, o6_sha256, character(1)),
  stringsAsFactors = FALSE
)
o6_write_csv_once(manifest, file.path(O6_OUTPUTS, "run_control/06O_SCRIPT_MANIFEST.csv"))

seed <- utils::read.csv(file.path(O6_AMENDMENT, "06O_MASTER_SEED_REGISTRY.csv"), check.names = FALSE)
stopifnot(nrow(seed) == 78L, length(unique(seed$derived_seed)) == 78L)
seed$execution_status <- "AUTHORIZED_PENDING_ENDPOINT_EXECUTION"
o6_write_csv_once(seed, file.path(O6_OUTPUTS, "run_control/06O_SEED_REGISTRY.csv"))

annotation_paths <- c(
  GPL570 = file.path(O6_PROJECT, "04_tables/GSE10846/GSE10846_GPL570_probe_to_symbol_mapping.csv"),
  IlluminaHumanv4 = file.path(O6_PROJECT, "results/GSE181063_extension/GSE181063_illuminaHumanv4_probe_to_symbol_mapping.csv")
)
stopifnot(all(file.exists(annotation_paths)))
pkg <- data.frame(package = packages,
  version = vapply(packages, function(x) as.character(utils::packageVersion(x)), character(1)),
  namespace_path = vapply(packages, function(x) normalizePath(getNamespaceInfo(asNamespace(x), "path"), winslash = "/"), character(1)))
parallel_vars <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
  "RCPP_PARALLEL_NUM_THREADS", "MC_CORES", "R_FUTURE_PLAN")
capture <- c(
  "06o runtime and package capture",
  paste("captured_at", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), sep = "\t"),
  paste("R_version", R.version.string, sep = "\t"),
  paste("platform", R.version$platform, sep = "\t"),
  paste("Bioconductor_version", as.character(BiocManager::version()), sep = "\t"),
  paste("master_seed", O6_MASTER_SEED, sep = "\t"),
  paste("parallel_backend", "BiocParallel::SerialParam", sep = "\t"),
  paste("library_paths", paste(.libPaths(), collapse = ";"), sep = "\t"),
  paste("locale", paste(names(Sys.getlocale()), Sys.getlocale(), sep = "=", collapse = ";"), sep = "\t"),
  paste("BLAS", extSoftVersion()[["BLAS"]], sep = "\t"),
  paste("LAPACK", extSoftVersion()[["LAPACK"]], sep = "\t"),
  paste("parallel_environment", paste(parallel_vars, Sys.getenv(parallel_vars, unset = "UNSET"), sep = "=", collapse = ";"), sep = "\t"),
  "",
  "PACKAGE_VERSIONS",
  capture.output(print(pkg, row.names = FALSE)),
  "",
  "LOCAL_ANNOTATION_DEPENDENCIES",
  paste(names(annotation_paths), normalizePath(annotation_paths, winslash = "/"), file.info(annotation_paths)$size,
    vapply(annotation_paths, o6_sha256, character(1)), sep = "\t"),
  "",
  "LOADED_NAMESPACES",
  sort(loadedNamespaces()),
  "",
  "RENV_STATUS",
  capture.output(try(renv::status(project = file.path(O6_REVISION,
    "04b_stage4_environment_freeze/stage4_renv_project")), silent = TRUE)),
  "",
  "SESSION_INFO",
  capture.output(sessionInfo())
)
o6_write_text_once(capture, file.path(O6_OUTPUTS, "run_control/06O_RUNTIME_AND_PACKAGE_CAPTURE.txt"))
cat("RUNTIME_CAPTURE=PASS\n")
