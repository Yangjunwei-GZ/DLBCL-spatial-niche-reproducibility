DLBCL_PROJECT_ROOT <- normalizePath(
  Sys.getenv("DLBCL_PROJECT_ROOT", unset = getwd()),
  winslash = "/", mustWork = FALSE
)
DLBCL_SUPPLEMENTARY_CODE_ROOT <- normalizePath(
  Sys.getenv("DLBCL_SUPPLEMENTARY_CODE_ROOT",
    unset = file.path(DLBCL_PROJECT_ROOT, "Supplementary_Code")),
  winslash = "/", mustWork = FALSE
)
DLBCL_SUPPLEMENTARY_DATA_ROOT <- normalizePath(
  Sys.getenv("DLBCL_SUPPLEMENTARY_DATA_ROOT",
    unset = file.path(DLBCL_PROJECT_ROOT, "Supplementary_Data")),
  winslash = "/", mustWork = FALSE
)
DLBCL_SUPPLEMENTARY_TABLES_ROOT <- normalizePath(
  Sys.getenv("DLBCL_SUPPLEMENTARY_TABLES_ROOT",
    unset = file.path(DLBCL_PROJECT_ROOT, "Supplementary_Tables")),
  winslash = "/", mustWork = FALSE
)

# Figure 1 Panel B data-integrity repair
#
# Core conclusion: the manuscript uses six continuous canonical programs whose
# exact 132 memberships span 121 unique genes; no discrete taxonomy is implied.
# Evidence chain: Panel A states the continuous workflow, Panel B shows the exact
# binary membership matrix, and Panel C limits the cohort scope to four datasets.
# Figure archetype: schematic-led composite.
# Export contract: 183 mm wide, one-page vector PDF, true 600 dpi LZW RGB TIFF,
# and an R-rendered preview PNG.

options(stringsAsFactors = FALSE, warn = 1)

# R on Windows must use a UTF-8 native locale to open the authority workbook's
# Unicode directory path. This changes path handling only, never data values.
if (.Platform$OS.type == "windows") {
  locale_result <- suppressWarnings(
    Sys.setlocale("LC_ALL", "English_United States.utf8")
  )
  if (is.na(locale_result)) {
    stop("A Windows UTF-8 locale is required to read the authority workbook path.")
  }
}

required_packages <- c(
  "readxl", "ggplot2", "patchwork", "dplyr", "tidyr",
  "ragg", "digest", "zip"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
})

# -----------------------------------------------------------------------------
# 1. Frozen inputs and new output directory
# -----------------------------------------------------------------------------

project_root <- normalizePath(
  DLBCL_PROJECT_ROOT,
  winslash = "/",
  mustWork = TRUE
)

output_dir <- file.path(
  project_root,
  "reproduced_figures/Figure_1_Panel_B"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

authority_workbook <- normalizePath(
  file.path(
    DLBCL_SUPPLEMENTARY_TABLES_ROOT,
    "DLBCL_continuous_model_Supplementary_Tables_FINAL_SUBMISSION_PUBLICATION_READY.xlsx"
  ),
  winslash = "/",
  mustWork = TRUE
)
authority_sheet <- "S2A_Program_genes"
multiply_sign <- intToUtf8(215L)
figure_input_dir <- file.path(
  DLBCL_SUPPLEMENTARY_DATA_ROOT, "figure_rendering", "Figure_1"
)
canonical_manifest <- file.path(figure_input_dir, "canonical_programs_v2.csv")
frozen_s2a_export <- file.path(
  figure_input_dir, "source_S2A_memberships_export_2026-07-30.csv"
)
independent_signature_workbook <- file.path(
  figure_input_dir, "spatial_niche_signatures.xlsx"
)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) {
  stop("Run this workflow with Rscript so its exact plotting script can be archived.")
}
plot_script <- normalizePath(
  sub("^--file=", "", script_arg),
  winslash = "/",
  mustWork = TRUE
)

input_paths <- c(
  authority_workbook,
  canonical_manifest,
  frozen_s2a_export,
  independent_signature_workbook
)
if (any(!file.exists(input_paths))) {
  stop("Missing required frozen input(s): ", paste(input_paths[!file.exists(input_paths)], collapse = "; "))
}

program_order <- c(
  "Macrophage-rich program",
  "T cell-inflamed program",
  "Immune-inflamed / antigen-presentation program",
  "Stromal / fibrotic program",
  "Immune-cold / exclusion-associated program",
  "Proliferative / cycling program"
)

display_name <- c(
  "Macrophage-rich program" = "Macrophage-rich",
  "T cell-inflamed program" = "T cell-inflamed",
  "Immune-inflamed / antigen-presentation program" =
    "Antigen-presentation",
  "Stromal / fibrotic program" = "Stromal/fibrotic",
  "Immune-cold / exclusion-associated program" =
    "Immune-cold/exclusion",
  "Proliferative / cycling program" = "Proliferative/cycling"
)

# -----------------------------------------------------------------------------
# 2. Read the final workbook first, then cross-check every frozen authority
# -----------------------------------------------------------------------------

workbook_raw <- readxl::read_excel(
  authority_workbook,
  sheet = authority_sheet,
  skip = 3,
  .name_repair = "minimal"
)
required_workbook_columns <- c("Program name", "Gene symbol")
if (!all(required_workbook_columns %in% names(workbook_raw))) {
  stop("Final workbook S2A lacks required columns: ",
       paste(setdiff(required_workbook_columns, names(workbook_raw)), collapse = ", "))
}

workbook_pairs <- workbook_raw |>
  transmute(
    program = as.character(.data[["Program name"]]),
    gene = as.character(.data[["Gene symbol"]]),
    source_excel_row = row_number() + 4L
  ) |>
  filter(!is.na(program) | !is.na(gene))

canonical_raw <- read.csv(
  canonical_manifest,
  check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)
canonical_pairs <- canonical_raw |>
  filter(toupper(as.character(canonical_membership)) == "TRUE") |>
  transmute(
    program = as.character(normalized_program_name),
    gene = as.character(gene_symbol),
    membership_index = as.integer(membership_index)
  )

source_export_raw <- read.csv(
  frozen_s2a_export,
  check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)
source_export_pairs <- source_export_raw |>
  transmute(
    program = as.character(normalized_program_name),
    gene = as.character(gene_symbol),
    membership_index = as.integer(membership_index)
  )

independent_raw <- readxl::read_excel(
  independent_signature_workbook,
  sheet = "spatial_niche_signatures",
  .name_repair = "minimal"
)
required_independent_columns <- c("niche_name", "gene_symbol")
if (!all(required_independent_columns %in% names(independent_raw))) {
  stop("Independent signature workbook lacks required columns: ",
       paste(setdiff(required_independent_columns, names(independent_raw)), collapse = ", "))
}

independent_observed <- independent_raw |>
  transmute(
    independent_program_name = as.character(niche_name),
    independent_gene = as.character(gene_symbol),
    source_excel_row = row_number() + 1L
  )

# The frozen canonical manifest links every membership to an exact row in the
# independent workbook. Use that explicit row-level mapping; never fuzzy-match
# or silently translate the independent workbook's historical program labels.
independent_links <- canonical_raw |>
  filter(toupper(as.character(canonical_membership)) == "TRUE") |>
  transmute(
    program = as.character(normalized_program_name),
    gene = as.character(gene_symbol),
    source_excel_row = as.integer(sub(
      ".*::row=", "", as.character(source_independent_workbook)
    ))
  )
independent_joined <- independent_links |>
  left_join(independent_observed, by = "source_excel_row")
independent_link_error <- any(
  is.na(independent_joined$independent_gene) |
    independent_joined$gene != independent_joined$independent_gene
)
independent_pairs <- independent_joined |>
  transmute(program, gene = independent_gene, source_excel_row)

pair_key <- function(x) paste(x$program, x$gene, sep = "\r")
pair_set_difference <- function(a, b) {
  a[!(pair_key(a) %in% pair_key(b)), c("program", "gene"), drop = FALSE]
}

source_checks <- list(
  final_workbook = workbook_pairs,
  canonical_manifest = canonical_pairs,
  frozen_s2a_export = source_export_pairs,
  independent_signature_workbook = independent_pairs
)

source_anomalies <- lapply(names(source_checks), function(source_name) {
  x <- source_checks[[source_name]]
  data.frame(
    source = source_name,
    empty_program = sum(is.na(x$program) | x$program == ""),
    empty_gene = sum(is.na(x$gene) | x$gene == ""),
    program_whitespace = sum(!is.na(x$program) & x$program != trimws(x$program)),
    gene_whitespace = sum(!is.na(x$gene) & x$gene != trimws(x$gene)),
    duplicated_pairs = sum(duplicated(pair_key(x))),
    unknown_programs = sum(!is.na(x$program) & !(x$program %in% program_order)),
    stringsAsFactors = FALSE
  )
}) |>
  bind_rows()

case_inconsistency <- function(x) {
  z <- split(x$gene, tolower(x$gene))
  sum(vapply(z, function(v) length(unique(v[!is.na(v)])) > 1L, logical(1)))
}
source_anomalies$case_inconsistent_gene_groups <- vapply(
  source_checks,
  case_inconsistency,
  integer(1)
)

wb_vs_canonical_a <- pair_set_difference(workbook_pairs, canonical_pairs)
wb_vs_canonical_b <- pair_set_difference(canonical_pairs, workbook_pairs)
wb_vs_export_a <- pair_set_difference(workbook_pairs, source_export_pairs)
wb_vs_export_b <- pair_set_difference(source_export_pairs, workbook_pairs)
wb_vs_independent_a <- pair_set_difference(workbook_pairs, independent_pairs)
wb_vs_independent_b <- pair_set_difference(independent_pairs, workbook_pairs)

ordered_workbook <- workbook_pairs |>
  mutate(program_order = match(program, program_order)) |>
  arrange(program_order, source_excel_row)
ordered_canonical <- canonical_pairs |>
  mutate(program_order = match(program, program_order)) |>
  arrange(program_order, membership_index)
within_program_order_agrees <- identical(
  pair_key(ordered_workbook),
  pair_key(ordered_canonical)
)

definition_conflict <- any(source_anomalies[, -1] != 0L) ||
  any(vapply(
    list(
      wb_vs_canonical_a, wb_vs_canonical_b,
      wb_vs_export_a, wb_vs_export_b,
      wb_vs_independent_a, wb_vs_independent_b
    ),
    nrow,
    integer(1)
  ) > 0L) ||
  independent_link_error ||
  !within_program_order_agrees

if (definition_conflict) {
  conflict_path <- file.path(output_dir, "DATA_SOURCE_CONFLICT_REPORT.md")
  conflict_lines <- c(
    "# Figure 1 Panel B data-source conflict",
    "",
    "Final figure generation was stopped because frozen sources did not pass the exact agreement gate.",
    "",
    paste0("- Final workbook: `", authority_workbook, "`"),
    paste0("- Sheet: `", authority_sheet, "`"),
    paste0("- Canonical manifest: `", canonical_manifest, "`"),
    paste0("- Frozen S2A export: `", frozen_s2a_export, "`"),
    paste0("- Independent signature workbook: `", independent_signature_workbook, "`"),
    "",
    "## Source anomaly counts",
    "",
    paste(capture.output(print(source_anomalies, row.names = FALSE)), collapse = "\n"),
    "",
    paste0("- Independent-workbook frozen row links agree: ", !independent_link_error),
    paste0("- Within-program frozen order agrees: ", within_program_order_agrees),
    paste0("- Workbook-only versus canonical: ", nrow(wb_vs_canonical_a)),
    paste0("- Canonical-only versus workbook: ", nrow(wb_vs_canonical_b)),
    paste0("- Workbook-only versus frozen export: ", nrow(wb_vs_export_a)),
    paste0("- Frozen-export-only versus workbook: ", nrow(wb_vs_export_b)),
    paste0("- Workbook-only versus independent workbook: ", nrow(wb_vs_independent_a)),
    paste0("- Independent-workbook-only versus workbook: ", nrow(wb_vs_independent_b))
  )
  writeLines(conflict_lines, conflict_path, useBytes = TRUE)
  stop("Frozen program-gene authorities disagree. See ", conflict_path)
}

# Remove only a transient conflict report created by an earlier failed run of
# this same new script after the exact frozen-source gate now passes.
transient_conflict_report <- file.path(output_dir, "DATA_SOURCE_CONFLICT_REPORT.md")
if (file.exists(transient_conflict_report)) unlink(transient_conflict_report)

# -----------------------------------------------------------------------------
# 3. Build the exact 6 x 121 binary membership matrix
# -----------------------------------------------------------------------------

manifest <- ordered_workbook |>
  select(program, gene, source_excel_row) |>
  mutate(
    program_order = match(program, program_order),
    membership_index = ave(
      source_excel_row,
      program,
      FUN = function(z) seq_along(z)
    )
  ) |>
  arrange(program_order, membership_index)

# No explicit global gene display order exists in S2A. Following the frozen rule,
# genes are added on first appearance while scanning the six fixed programs.
unique_gene_order <- unique(manifest$gene)
gene_order_lookup <- setNames(seq_along(unique_gene_order), unique_gene_order)

membership_matrix <- matrix(
  0L,
  nrow = length(program_order),
  ncol = length(unique_gene_order),
  dimnames = list(program_order, unique_gene_order)
)
membership_matrix[cbind(
  match(manifest$program, program_order),
  unname(gene_order_lookup[manifest$gene])
)] <- 1L

row_sums <- rowSums(membership_matrix)
total_present <- sum(membership_matrix)
unique_gene_n <- length(unique_gene_order)
unique_pair_n <- nrow(unique(manifest[c("program", "gene")]))
extra_memberships <- total_present - unique_gene_n

assertions <- c(
  matrix_dimensions = identical(dim(membership_matrix), c(6L, 121L)),
  total_present = identical(as.integer(total_present), 132L),
  all_row_sums_22 = all(row_sums == 22L),
  unique_gene_count = identical(as.integer(unique_gene_n), 121L),
  unique_program_gene_pairs = identical(as.integer(unique_pair_n), 132L),
  binary_cells_only = all(membership_matrix %in% c(0L, 1L)),
  no_pair_duplicates = !any(duplicated(pair_key(manifest))),
  no_empty_program = !any(is.na(manifest$program) | manifest$program == ""),
  no_empty_gene = !any(is.na(manifest$gene) | manifest$gene == ""),
  no_unknown_program = all(manifest$program %in% program_order),
  no_whitespace_anomaly = all(manifest$program == trimws(manifest$program)) &&
    all(manifest$gene == trimws(manifest$gene)),
  no_case_inconsistency = case_inconsistency(manifest) == 0L,
  all_frozen_sources_agree = TRUE,
  frozen_within_program_order_agrees = within_program_order_agrees
)

if (!all(assertions)) {
  stop(
    "Panel B assertions failed: ",
    paste(names(assertions)[!assertions], collapse = ", "),
    ". FINAL outputs were not generated."
  )
}

# -----------------------------------------------------------------------------
# 4. Write machine-readable QC and source data
# -----------------------------------------------------------------------------

matrix_qc <- data.frame(
  program = unname(display_name[program_order]),
  expected_memberships = 22L,
  observed_memberships = as.integer(row_sums[program_order]),
  pass = row_sums[program_order] == 22L,
  stringsAsFactors = FALSE
)

shared_genes <- manifest |>
  group_by(gene) |>
  summarise(
    n_programs = n_distinct(program),
    programs = paste(
      unname(
        .env$display_name[
          .env$program_order[.env$program_order %in% unique(.data$program)]
        ]
      ),
      collapse = "; "
    ),
    extra_memberships = n_programs - 1L,
    .groups = "drop"
  ) |>
  filter(n_programs > 1L) |>
  arrange(unname(gene_order_lookup[gene]))

matrix_long <- manifest |>
  transmute(
    program = unname(display_name[program]),
    gene,
    program_order,
    gene_order = unname(gene_order_lookup[gene]),
    present = 1L
  ) |>
  arrange(program_order, gene_order)

matrix_wide <- as.data.frame(membership_matrix, check.names = FALSE) |>
  tibble::rownames_to_column("program") |>
  mutate(program = unname(display_name[program]))

write.csv(
  matrix_qc,
  file.path(output_dir, "PANEL_B_MATRIX_QC.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  shared_genes,
  file.path(output_dir, "PANEL_B_SHARED_GENES.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  matrix_long,
  file.path(output_dir, "PANEL_B_MATRIX_LONG.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  matrix_wide,
  file.path(output_dir, "PANEL_B_MATRIX_WIDE.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

input_manifest <- data.frame(
  role = c(
    "primary_authority_workbook",
    "frozen_canonical_manifest",
    "frozen_S2A_export",
    "independent_signature_cross_check"
  ),
  file_name = basename(input_paths),
  absolute_path = normalizePath(input_paths, winslash = "/", mustWork = TRUE),
  sheet = c(authority_sheet, NA, NA, "spatial_niche_signatures"),
  file_size_bytes = unname(file.info(input_paths)$size),
  sha256 = vapply(
    input_paths,
    digest::digest,
    character(1),
    algo = "sha256",
    file = TRUE
  ),
  stringsAsFactors = FALSE
)
write.csv(
  input_manifest,
  file.path(output_dir, "INPUT_MANIFEST_SHA256.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# -----------------------------------------------------------------------------
# 5. Generate the QC report and final figure legend
# -----------------------------------------------------------------------------

shared_lines <- paste0(
  "- `", shared_genes$gene, "`: ", shared_genes$n_programs,
  " programs (", shared_genes$programs, ")"
)
assertion_lines <- paste0(
  "- ", names(assertions), ": ", ifelse(assertions, "PASS", "FAIL")
)

qc_report <- c(
  "# Figure 1 Panel B QC report",
  "",
  "## Authority",
  "",
  paste0("- Authoritative input: `", authority_workbook, "`"),
  paste0("- Sheet: `", authority_sheet, "`"),
  "- Program column: `Program name`",
  "- Gene column: `Gene symbol`",
  paste0("- Workbook manifest rows: ", nrow(workbook_pairs)),
  paste0("- Unique program-gene pairs: ", unique_pair_n),
  paste0("- Input workbook SHA-256: `", input_manifest$sha256[1], "`"),
  "- Frozen cross-checks: canonical manifest, S2A export, and independent signature workbook all agree exactly on the 132 memberships.",
  "- Ordering rule: fixed program order; within-program S2A row order (confirmed identical to frozen `membership_index`); global genes ordered by first appearance.",
  "",
  "## Matrix",
  "",
  paste0("Matrix dimensions: 6 ", multiply_sign, " 121"),
  "Total present cells: 132",
  "Row sums: 22, 22, 22, 22, 22, 22",
  paste0("Unique genes: ", unique_gene_n),
  paste0("Shared genes: ", nrow(shared_genes)),
  paste0("Extra cross-program memberships: 132 - 121 = ", extra_memberships),
  "The extra-membership count is distinct from the shared-gene count; each shared gene's program count is reported below.",
  "",
  "## Programs",
  "",
  paste0("- ", unname(display_name[program_order]), ": 22 memberships"),
  "",
  "## Shared genes",
  "",
  shared_lines,
  "",
  "## Assertions",
  "",
  assertion_lines,
  "",
  "## Integrity statement",
  "",
  "No gene symbol or scientific program definition was added, deleted, replaced, corrected, or case-normalized. Only the user-specified display labels were applied.",
  "No existing Figure 1, workbook, or manifest was overwritten.",
  "",
  "FINAL_PANEL_B_STATUS = PASS"
)
writeLines(
  qc_report,
  file.path(output_dir, "FIGURE1_PANEL_B_QC_REPORT.md"),
  useBytes = TRUE
)

figure_legend <- c(
  "# Figure 1. Canonical continuous-program framework and study design",
  "",
  "**A, Continuous-model study design.** Six canonical programs are scored in the GSE31312 discovery cohort, followed by continuous-architecture assessment and prespecified model adjudication. No value of k was selected and no taxonomy was assigned; the retained representation comprises six continuous scores plus PC1/PC2. External structural replication uses GSE10846 and GSE181063, and continuous spatial analysis uses GSE276542.",
  "",
  "**B, Canonical program membership.** Exact membership matrix for six prespecified 22-gene programs across 121 unique genes. Coral cells indicate membership and light-gray cells indicate absence. Shared genes retain membership in every prespecified program to which they belong. The 132 memberships therefore comprise 121 unique genes and 11 additional cross-program memberships.",
  "",
  "**C, Cohort scope.** GSE31312 is the discovery cohort; GSE10846 and GSE181063 provide external structural replication; GSE276542 provides continuous spatial analysis. No class, centroid, survival, spatial-taxonomy, or ligand-receptor inference is shown."
)
writeLines(
  figure_legend,
  file.path(output_dir, "Figure_1_LEGEND_FINAL.md"),
  useBytes = TRUE
)

# -----------------------------------------------------------------------------
# 6. Draw Panels A-C in R; Panel B is a true 6 x 121 vector tile matrix
# -----------------------------------------------------------------------------

theme_set(
  theme_void(base_family = "Arial") +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      plot.title = element_text(
        family = "Arial", size = 8.8, face = "bold",
        hjust = 0, margin = margin(b = 4)
      )
    )
)

workflow_boxes <- data.frame(
  xmin = c(0.15, 3.50, 6.85, 3.50, 0.15, 6.85),
  xmax = c(3.15, 6.50, 9.85, 6.50, 3.15, 9.85),
  ymin = c(1.35, 1.35, 1.35, 0.15, 0.15, 0.15),
  ymax = c(2.15, 2.15, 2.15, 0.95, 0.95, 0.95),
  label = c(
    "Six canonical programs\nGSE31312 discovery cohort",
    "Continuous program scoring",
    "Continuous architecture\nPrespecified model adjudication",
    "No k selected\nTaxonomy not assigned\nSix scores + PC1/PC2",
    "GSE10846/GSE181063\nExternal structural replication",
    "GSE276542\nContinuous spatial analysis"
  ),
  fill = c(
    "#F2F2F2", "#F2F2F2", "#F2F2F2",
    "#FCE8E5", "#E8EEF2", "#E8EEF2"
  )
)
workflow_arrows <- data.frame(
  x = c(3.17, 6.52, 8.35, 3.48, 6.52),
  xend = c(3.46, 6.81, 5.00, 3.19, 6.81),
  y = c(1.75, 1.75, 1.33, 0.55, 0.55),
  yend = c(1.75, 1.75, 0.97, 0.55, 0.55)
)

p_a <- ggplot() +
  geom_rect(
    data = workflow_boxes,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
    colour = "#666666", linewidth = 0.28
  ) +
  scale_fill_identity() +
  geom_text(
    data = workflow_boxes,
    aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
    family = "Arial", size = 2.65, lineheight = 0.94, colour = "#202020"
  ) +
  geom_segment(
    data = workflow_arrows,
    aes(x = x, xend = xend, y = y, yend = yend),
    linewidth = 0.35, colour = "#555555",
    arrow = grid::arrow(length = grid::unit(1.3, "mm"), type = "closed")
  ) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 2.35), clip = "off", expand = FALSE) +
  labs(title = "A  Continuous-model study design") +
  theme(plot.margin = margin(3, 3, 4, 3))

heat_df <- expand.grid(
  program = program_order,
  gene = unique_gene_order,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
) |>
  mutate(
    present = as.integer(membership_matrix[cbind(
      match(program, program_order),
      match(gene, unique_gene_order)
    )]),
    program_display = factor(
      unname(display_name[program]),
      levels = rev(unname(display_name[program_order]))
    ),
    gene = factor(gene, levels = unique_gene_order)
  )

p_b <- ggplot(heat_df, aes(x = gene, y = program_display, fill = factor(present))) +
  geom_tile(colour = "white", linewidth = 0.075) +
  scale_fill_manual(
    values = c("0" = "#F1F1F1", "1" = "#D8756A"),
    breaks = c("0", "1"),
    labels = c("Absent", "Present"),
    name = NULL,
    drop = FALSE
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(add = 0.08)) +
  labs(
    title = "B  Canonical program membership",
    x = "Canonical genes ordered by first program membership",
    y = NULL,
    caption = paste0(
      "6 ", multiply_sign, " 22 memberships = 132; 121 unique genes.\n",
      "Shared genes retain membership in each prespecified program."
    )
  ) +
  theme_minimal(base_family = "Arial", base_size = 7.5) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 7.2, colour = "#202020", hjust = 1),
    axis.title.x = element_text(size = 7.2, colour = "#303030", margin = margin(t = 5)),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.justification = "right",
    legend.direction = "horizontal",
    legend.text = element_text(size = 7.0),
    legend.key.width = grid::unit(4.0, "mm"),
    legend.key.height = grid::unit(2.5, "mm"),
    legend.margin = margin(0, 0, 1, 0),
    plot.title = element_text(size = 8.8, face = "bold", margin = margin(b = 1)),
    plot.caption = element_text(size = 7.0, colour = "#404040", hjust = 0, lineheight = 1.05, margin = margin(t = 5)),
    plot.margin = margin(2, 3, 4, 3)
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

cohort_df <- data.frame(
  x = 1:4,
  cohort = c("GSE31312", "GSE10846", "GSE181063", "GSE276542"),
  role = c(
    "Discovery cohort",
    "External structural\nreplication",
    "External structural\nreplication",
    "Continuous spatial\nanalysis"
  ),
  stringsAsFactors = FALSE
)

p_c <- ggplot(cohort_df, aes(x = x, y = 1)) +
  geom_tile(
    width = 0.88, height = 0.62,
    fill = "#F2F2F2", colour = "#777777", linewidth = 0.28
  ) +
  geom_text(
    aes(label = paste0(cohort, "\n", role)),
    family = "Arial", size = 2.7, lineheight = 0.95, colour = "#202020"
  ) +
  coord_cartesian(xlim = c(0.45, 4.55), ylim = c(0.62, 1.38), expand = FALSE) +
  labs(title = "C  Cohort scope") +
  theme(plot.margin = margin(3, 3, 3, 3))

figure_1 <- p_a / p_b / p_c +
  patchwork::plot_layout(heights = c(1.22, 1.52, 0.66)) +
  patchwork::plot_annotation(
    theme = theme(plot.background = element_rect(fill = "white", colour = NA))
  )

# -----------------------------------------------------------------------------
# 7. Export submission files from the selected R backend only
# -----------------------------------------------------------------------------

width_in <- 183 / 25.4
height_in <- 215 / 25.4
pdf_path <- file.path(output_dir, "Figure_1_CONTINUOUS_MODEL_FINAL.pdf")
tiff_path <- file.path(output_dir, "Figure_1_CONTINUOUS_MODEL_FINAL_600dpi.tiff")
preview_path <- file.path(output_dir, "Figure_1_CONTINUOUS_MODEL_FINAL_preview.png")

grDevices::cairo_pdf(
  pdf_path,
  width = width_in,
  height = height_in,
  family = "Arial",
  onefile = FALSE,
  bg = "white"
)
print(figure_1)
grDevices::dev.off()

ragg::agg_tiff(
  tiff_path,
  width = width_in,
  height = height_in,
  units = "in",
  res = 600,
  compression = "lzw",
  background = "white"
)
print(figure_1)
grDevices::dev.off()

ragg::agg_png(
  preview_path,
  width = width_in,
  height = height_in,
  units = "in",
  res = 300,
  background = "white"
)
print(figure_1)
grDevices::dev.off()

# -----------------------------------------------------------------------------
# 8. Record runtime provenance and build the submission package
# -----------------------------------------------------------------------------

runtime_lines <- c(
  paste0("R version: ", R.version.string),
  paste0("Platform: ", R.version$platform),
  paste0("Figure backend: R only"),
  paste0("Figure dimensions inches: ", width_in, " x ", height_in),
  "TIFF resolution: 600 dpi",
  "TIFF compression: LZW",
  paste0(
    required_packages, ": ",
    vapply(required_packages, function(x) as.character(utils::packageVersion(x)), character(1))
  )
)
writeLines(
  runtime_lines,
  file.path(output_dir, "FIGURE1_RUNTIME_PROVENANCE.txt"),
  useBytes = TRUE
)

package_files <- file.path(
  output_dir,
  c(
    "Figure_1_CONTINUOUS_MODEL_FINAL.pdf",
    "Figure_1_CONTINUOUS_MODEL_FINAL_600dpi.tiff",
    "Figure_1_CONTINUOUS_MODEL_FINAL_preview.png",
    "Figure_1_LEGEND_FINAL.md",
    "FIGURE1_PANEL_B_QC_REPORT.md",
    "PANEL_B_MATRIX_QC.csv",
    "PANEL_B_MATRIX_WIDE.csv",
    "PANEL_B_MATRIX_LONG.csv",
    "PANEL_B_SHARED_GENES.csv",
    "INPUT_MANIFEST_SHA256.csv",
    "FIGURE1_RUNTIME_PROVENANCE.txt",
    basename(plot_script)
  )
)
if (any(!file.exists(package_files))) {
  stop("Submission package member missing: ",
       paste(package_files[!file.exists(package_files)], collapse = "; "))
}

zip_path <- file.path(
  output_dir,
  "Figure_1_CONTINUOUS_MODEL_FINAL_SUBMISSION_PACKAGE.zip"
)
if (file.exists(zip_path)) file.remove(zip_path)
old_wd <- setwd(output_dir)
on.exit(setwd(old_wd), add = TRUE)
zip::zipr(
  zipfile = basename(zip_path),
  files = basename(package_files),
  recurse = FALSE,
  include_directories = FALSE
)
setwd(old_wd)

cat("FINAL_PANEL_B_STATUS = PASS\n")
cat("Matrix dimensions: 6 x 121\n")
cat("Total present cells: 132\n")
cat("Row sums: ", paste(row_sums, collapse = ", "), "\n", sep = "")
cat("Unique genes: ", unique_gene_n, "\n", sep = "")
cat("Shared genes: ", nrow(shared_genes), "\n", sep = "")
cat("Extra memberships: ", extra_memberships, "\n", sep = "")
cat("Output directory: ", output_dir, "\n", sep = "")
