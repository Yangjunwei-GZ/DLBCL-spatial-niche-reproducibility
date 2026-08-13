DLBCL_PROJECT_ROOT <- normalizePath(
  Sys.getenv("DLBCL_PROJECT_ROOT", unset = getwd()),
  winslash = "/", mustWork = FALSE
)
DLBCL_SUPPLEMENTARY_CODE_ROOT <- normalizePath(
  Sys.getenv("DLBCL_SUPPLEMENTARY_CODE_ROOT",
    unset = file.path(DLBCL_PROJECT_ROOT, "Supplementary_Code")),
  winslash = "/", mustWork = FALSE
)

stage4b_environment_root <- file.path(DLBCL_SUPPLEMENTARY_CODE_ROOT, "environment")
stage4b_renv_project <- stage4b_environment_root
stage4b_lock_path <- file.path(stage4b_renv_project, "renv.lock")
stage4b_expected_lock_sha256 <- "abf1607763905c0afbbafd75834d28ac2865781064601cb374ea02e7425a736e"
stage4b_expected_R_major_minor <- "4.5"
stage4b_required_token <- "EXPLICITLY_APPROVED_FUTURE_STAGE"
stage4c1_required_packages <- c("ConsensusClusterPlus", "cluster", "digest", "MASS")
stage4c1_optional_packages <- "diptest"
stage4c2_gsva_required_packages <- c(
  "ConsensusClusterPlus", "cluster", "GSVA", "Biobase", "data.table",
  "digest", "MASS", "Matrix", "ggplot2", "clue"
)
stage4c2_gsva_optional_packages <- c("corpcor", "diptest")
stage4b_required_packages <- stage4c1_required_packages
stage4b_optional_packages <- stage4c1_optional_packages
stage4b_R_distribution_packages <- c("cluster", "MASS", "Matrix")

stage4b_assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

stage4b_normalize <- function(path, must_work = FALSE) {
  normalizePath(path, winslash = "/", mustWork = must_work)
}

stage4b_science_gate <- function() {
  token <- Sys.getenv("DLBCL_REVISION_ALLOW_SCIENCE", unset = "")
  stage4b_assert(
    identical(token, stage4b_required_token),
    paste0(
      "Scientific execution is disabled. Set DLBCL_REVISION_ALLOW_SCIENCE=",
      stage4b_required_token,
      " only in a separately authorized future stage."
    )
  )
}

stage4b_environment_activate <- function(
  require_science_token = TRUE,
  verbose = TRUE,
  profile = c("stage4c1", "stage4c2_gsva")
) {
  profile <- match.arg(profile)
  if (isTRUE(require_science_token)) stage4b_science_gate()
  required_packages <- if (identical(profile, "stage4c1")) {
    stage4c1_required_packages
  } else {
    stage4c2_gsva_required_packages
  }
  optional_packages <- if (identical(profile, "stage4c1")) {
    stage4c1_optional_packages
  } else {
    stage4c2_gsva_optional_packages
  }

  state_root <- file.path(stage4b_environment_root, "renv_state")
  Sys.setenv(
    RENV_PROJECT = stage4b_renv_project,
    RENV_PATHS_ROOT = state_root,
    RENV_PATHS_CACHE = file.path(state_root, "cache"),
    RENV_PATHS_BINARY = file.path(state_root, "binary"),
    RENV_PATHS_SOURCE = file.path(state_root, "source"),
    RENV_PATHS_SANDBOX = file.path(state_root, "sandbox"),
    RENV_DOWNLOAD_METHOD = "libcurl",
    RENV_CONFIG_CONSENT = "TRUE",
    RENV_CONFIG_AUTO_SNAPSHOT = "FALSE"
  )

  activate_path <- file.path(stage4b_renv_project, "renv", "activate.R")
  stage4b_assert(file.exists(activate_path), paste("Missing renv activation script:", activate_path))
  source(activate_path, local = FALSE)

  observed_R <- paste(R.version$major, strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1], sep = ".")
  stage4b_assert(
    identical(observed_R, stage4b_expected_R_major_minor),
    paste("Stage 4B requires R", stage4b_expected_R_major_minor, "but found", observed_R)
  )
  stage4b_assert(file.exists(stage4b_lock_path), paste("Missing renv.lock:", stage4b_lock_path))

  digest_path <- find.package("digest", quiet = TRUE)
  stage4b_assert(nzchar(digest_path), "Required package digest is unavailable after renv activation.")
  observed_lock_sha256 <- digest::digest(file = stage4b_lock_path, algo = "sha256")
  stage4b_assert(
    identical(observed_lock_sha256, stage4b_expected_lock_sha256),
    paste("renv.lock SHA-256 mismatch. Expected", stage4b_expected_lock_sha256, "observed", observed_lock_sha256)
  )

  active_libs <- stage4b_normalize(.libPaths(), must_work = TRUE)
  approved_root <- paste0(tolower(stage4b_normalize(stage4b_environment_root, must_work = TRUE)), "/")
  R_distribution_root <- paste0(tolower(stage4b_normalize(R.home("library"), must_work = TRUE)), "/")
  stage4b_assert(
    all(startsWith(paste0(tolower(active_libs), "/"), approved_root)),
    paste("Unapproved library root after renv activation:", paste(active_libs, collapse = " | "))
  )

  package_rows <- lapply(c(required_packages, optional_packages), function(package) {
    candidates <- stage4b_normalize(
      file.path(active_libs, package)[file.exists(file.path(active_libs, package))],
      must_work = TRUE
    )
    available <- length(candidates) > 0L
    required <- package %in% required_packages
    if (required) {
      stage4b_assert(available, paste("Required Stage 4 package is unavailable:", package))
      stage4b_assert(length(candidates) == 1L, paste("Required package resolves ambiguously:", package))
    }
    resolved <- if (available) stage4b_normalize(find.package(package), must_work = TRUE) else ""
    if (available) {
      in_stage4b <- startsWith(paste0(tolower(resolved), "/"), approved_root)
      allowed_R_distribution <- package %in% stage4b_R_distribution_packages &&
        startsWith(paste0(tolower(resolved), "/"), R_distribution_root)
      stage4b_assert(
        in_stage4b || allowed_R_distribution,
        paste("Package resolved outside Stage 4B environment:", package, resolved)
      )
    }
    data.frame(
      package = package,
      required = required,
      available = available,
      version = if (available) as.character(packageVersion(package)) else "",
      resolved_path = resolved,
      candidate_count = length(candidates),
      candidate_paths = paste(candidates, collapse = " | "),
      stringsAsFactors = FALSE
    )
  })
  package_resolution <- do.call(rbind, package_rows)

  context <- list(
    environment_root = stage4b_environment_root,
    renv_project = stage4b_renv_project,
    lock_path = stage4b_lock_path,
    lock_sha256 = observed_lock_sha256,
    R_major_minor = observed_R,
    library_paths = active_libs,
    package_resolution = package_resolution,
    optional_available = setNames(
      package_resolution$available[package_resolution$package %in% optional_packages],
      package_resolution$package[package_resolution$package %in% optional_packages]
    ),
    profile = profile,
    required_packages = required_packages,
    optional_packages = optional_packages
  )
  options(stage4b.environment.context = context)

  if (identical(profile, "stage4c1")) {
    stage4b_assert(
      !("GSVA" %in% loadedNamespaces()),
      "GSVA was loaded during Stage 4C-1 environment activation."
    )
  }

  if (isTRUE(verbose)) {
    message("Stage 4B renv active.")
    message(".libPaths(): ", paste(active_libs, collapse = " | "))
    message("renv.lock SHA-256: ", observed_lock_sha256)
  }
  invisible(context)
}
