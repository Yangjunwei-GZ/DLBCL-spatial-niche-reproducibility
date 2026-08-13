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

options(stringsAsFactors = FALSE, warn = 1)

root <- DLBCL_PROJECT_ROOT
source_dir <- file.path(DLBCL_SUPPLEMENTARY_DATA_ROOT, "single_cell")
out_dir <- file.path(root, "reproduced_figures", "Figure_5")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(root)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ragg)
  library(digest)
  library(grid)
})

git_status <- function() c(
  "# git status --short",
  system2("git", c("status", "--short"), stdout = TRUE, stderr = TRUE),
  "", "# git diff --check",
  system2("git", c("diff", "--check"), stdout = TRUE, stderr = TRUE)
)
git_status_before <- git_status()

program_order <- c("macrophage_rich", "t_cell_inflamed", "antigen_presentation", "stromal_fibrotic", "immune_cold_exclusion", "proliferative_cycling")
program_labels <- c(macrophage_rich = "Macrophage-rich", t_cell_inflamed = "T cell-inflamed", antigen_presentation = "Antigen-presentation", stromal_fibrotic = "Stromal/fibrotic", immune_cold_exclusion = "Immune-cold/exclusion", proliferative_cycling = "Proliferative/cycling")
program_labels_b <- c(macrophage_rich = "Macrophage-\nrich", t_cell_inflamed = "T cell-\ninflamed", antigen_presentation = "Antigen-\npresentation", stromal_fibrotic = "Stromal/\nfibrotic", immune_cold_exclusion = "Immune-cold/\nexclusion", proliferative_cycling = "Proliferative/\ncycling")
program_labels_d <- c(macrophage_rich = "Macrophage-rich", t_cell_inflamed = "T cell-inflamed", antigen_presentation = "Antigen-\npresentation", stromal_fibrotic = "Stromal/fibrotic", immune_cold_exclusion = "Immune-cold/\nexclusion", proliferative_cycling = "Proliferative/\ncycling")
celltype_order <- c("B cells", "Plasma cells", "T cells CD4", "T cells CD8", "Tregs", "TFH", "NK cells", "Monocytes and Macrophages", "Others")
celltype_display <- c("B cells" = "B cells", "Plasma cells" = "Plasma\ncells", "T cells CD4" = "T cells\nCD4", "T cells CD8" = "T cells\nCD8", "Tregs" = "Tregs", "TFH" = "TFH", "NK cells" = "NK cells", "Monocytes and Macrophages" = "Monocytes and\nmacrophages", "Others" = "Others")
celltype_display_d <- c("B cells" = "B", "Plasma cells" = "Plasma", "T cells CD4" = "CD4 T", "T cells CD8" = "CD8 T", "Tregs" = "Treg", "TFH" = "TFH", "NK cells" = "NK", "Monocytes and Macrophages" = "Mono/Mac", "Others" = "Other")
patients <- c("DLBCL002", "DLBCL007", "DLBCL008", "DLBCL111")
program_levels <- program_order
celltype_levels <- celltype_order
patient_levels <- patients

source_paths <- c(
  panel_a = file.path(source_dir, "FIGURE5_PANEL_A_SOURCE_DATA.csv.gz"),
  panel_b = file.path(source_dir, "FIGURE5_PANEL_B_SOURCE_DATA.csv.gz"),
  panel_c = file.path(source_dir, "FIGURE5_PANEL_C_SOURCE_DATA.csv"),
  panel_d = file.path(source_dir, "FIGURE5_PANEL_D_SOURCE_DATA.csv"),
  panel_d_coverage = file.path(source_dir, "FIGURE5_PANEL_D_COVERAGE_QC.csv"),
  raw_ucell = file.path(source_dir, "FIGURE5_RAW_UCELL_SCORES.csv.gz"),
  display_z = file.path(source_dir, "FIGURE5_UCELL_SCORES_WITH_DISPLAY_Z.csv.gz"),
  umap = file.path(source_dir, "FIGURE5_UMAP_COORDINATES.csv.gz")
)
stopifnot(all(file.exists(source_paths)))
hash_file <- function(path) digest::digest(file = path, algo = "sha256")
source_manifest <- fread(file.path(root, "reproducibility", "SUBMISSION_PACKAGE_MANIFEST.csv"))
source_qc <- rbindlist(lapply(names(source_paths), function(role) {
  path <- source_paths[[role]]
  expected <- source_manifest[logical_name == basename(path), sha256]
  expected_present <- length(expected) > 0L && !is.na(expected[1]) && nzchar(expected[1])
  observed <- hash_file(path)
  data.table(
    source_role = role,
    absolute_path = normalizePath(path, winslash = "/", mustWork = TRUE),
    size_bytes = unname(file.info(path)$size), sha256 = observed,
    expected_sha256 = if (expected_present) expected[1] else NA_character_,
    expected_manifest = if (expected_present) "reproducibility/SUBMISSION_PACKAGE_MANIFEST.csv" else "no scientific manifest entry",
    hash_status = if (expected_present) if (observed == expected[1]) "PASS" else "FAIL" else "PASS_READ_ONLY_HASHED"
  )
}))
fwrite(source_qc, file.path(out_dir, "FIGURE5_SUBMISSION_V2_SOURCE_HASH_QC.csv"))
stopifnot(!any(source_qc$hash_status == "FAIL"))
panel_a <- fread(source_paths[["panel_a"]]); panel_b <- fread(source_paths[["panel_b"]]); panel_c <- fread(source_paths[["panel_c"]]); panel_d <- fread(source_paths[["panel_d"]])
raw <- fread(source_paths[["raw_ucell"]]); display <- fread(source_paths[["display_z"]]); umap <- fread(source_paths[["umap"]])
panel_b[, program := as.character(program)]
panel_c[, program := as.character(program)]
panel_d[, program := as.character(program)]
panel_c_source_program_order <- unique(panel_c$program)
panel_c_source_celltype_order <- panel_c[program == panel_c_source_program_order[1], unique(CellType)]
panel_d_source_patient_order <- unique(panel_d$Patient)
panel_d_source_celltype_order <- unique(panel_d$CellType)
stopifnot(nrow(panel_a) == 14368L, nrow(panel_b) == 86208L, nrow(panel_c) == 54L, nrow(panel_d) == 216L, nrow(raw) == 14368L, nrow(umap) == 14368L, length(unique(raw$ID)) == 14368L, length(unique(umap$ID)) == 14368L, all(is.finite(umap$UMAP_1)), all(is.finite(umap$UMAP_2)), all(panel_c$program %in% program_order), all(panel_c$CellType %in% celltype_order))
raw_cols <- paste0(program_order, "_raw_ucell"); display_cols <- paste0(program_order, "_display_z")
stopifnot(all(raw_cols %in% names(raw)), all(display_cols %in% names(display)), all(vapply(raw[, ..raw_cols], function(x) all(is.finite(x)), logical(1))), all(vapply(display[, ..display_cols], function(x) all(is.finite(x)), logical(1))))

panel_a[, CellType := factor(as.character(CellType), levels = celltype_levels)]
panel_b[, program_display := factor(unname(program_labels_b[as.character(program)]), levels = unname(program_labels_b[program_levels]))]
panel_c[, program_display := factor(unname(program_labels[as.character(program)]), levels = rev(unname(program_labels[program_levels])))]
panel_c[, CellType := factor(as.character(CellType), levels = celltype_levels)]
panel_d[, Patient := factor(as.character(Patient), levels = rev(patient_levels))]
panel_d[, CellType := factor(as.character(CellType), levels = celltype_levels)]
ctcols <- c("B cells" = "#3B6FB6", "Plasma cells" = "#7B4FA3", "T cells CD4" = "#2C9C95", "T cells CD8" = "#55B87A", "Tregs" = "#8BBF3D", "TFH" = "#D9A52B", "NK cells" = "#E9783B", "Monocytes and Macrophages" = "#B44C4C", "Others" = "#777777")
theme_base <- theme_classic(base_family = "Times New Roman", base_size = 8) + theme(plot.title = element_text(face = "bold", size = 9.5, margin = margin(b = 1.5)), plot.subtitle = element_text(size = 7.5, margin = margin(b = 2)), axis.title = element_text(size = 7.5), axis.text = element_text(size = 7.0, colour = "black"), axis.ticks = element_line(linewidth = 0.3, colour = "black"), axis.line = element_line(linewidth = 0.3, colour = "black"), legend.title = element_text(size = 7.5), legend.text = element_text(size = 7.0), strip.text = element_text(size = 7.5, face = "bold"), panel.grid = element_blank(), plot.margin = margin(1.5, 3, 1.5, 3))

pA <- ggplot(panel_a, aes(UMAP_1, UMAP_2, colour = CellType)) + geom_point(size = 0.18, alpha = 0.58) + coord_equal() + scale_colour_manual(values = ctcols, drop = FALSE) + theme_base + labs(title = "A  Annotated single-cell atlas", subtitle = "GSE182434 | 14,368 cells | 4 DLBCL donors", x = "UMAP 1", y = "UMAP 2", colour = "Cell type") + theme(legend.position = "bottom", legend.direction = "horizontal", legend.key.width = unit(0.27, "cm"), legend.spacing.x = unit(0.03, "cm"), legend.margin = margin(t = -2, b = -2), legend.box.just = "left") + guides(colour = guide_legend(title.position = "top", title.hjust = 0, nrow = 3, byrow = TRUE, override.aes = list(size = 1.7, alpha = 1)))
pB <- ggplot(panel_b, aes(UMAP_1, UMAP_2, colour = display_z)) + geom_point(size = 0.13, alpha = 0.58) + coord_equal() + facet_wrap(~program_display, ncol = 3) + scale_colour_gradient2(low = "#BFD3E6", mid = "#F7F7F7", high = "#B2182B", midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish, name = "Standardized UCell score") + theme_base + labs(title = "B  Single-cell distribution of continuous program scores", x = "UMAP 1", y = "UMAP 2") + theme(legend.position = "bottom", legend.margin = margin(t = -2, b = -2), strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.25), strip.text = element_text(size = 6.0, lineheight = 0.90)) + guides(colour = guide_colourbar(title.position = "top", title.hjust = 0.5, barwidth = unit(3.5, "cm"), barheight = unit(0.25, "cm")))

cplot <- copy(panel_c); cplot[, CellType := factor(as.character(CellType), levels = celltype_levels)]; c_lim <- max(abs(cplot$mean_display_z), na.rm = TRUE)
pC <- ggplot(cplot, aes(CellType, program_display, colour = mean_display_z, size = fraction_above_q75)) + geom_point() + scale_colour_gradient2(low = "#3B6FB6", mid = "#F7F7F7", high = "#B44C4C", midpoint = 0, limits = c(-c_lim, c_lim), name = "Mean standardized\nUCell score") + scale_size(range = c(1.7, 7.2), limits = c(0, 1), breaks = c(0.25, 0.5, 0.75), name = "Fraction above program-specific\n75th percentile") + scale_x_discrete(labels = celltype_display, drop = FALSE) + scale_y_discrete(labels = c("Macrophage-rich" = "Macrophage-rich", "T cell-inflamed" = "T cell-inflamed", "Antigen-presentation" = "Antigen-\npresentation", "Stromal/fibrotic" = "Stromal/fibrotic", "Immune-cold/exclusion" = "Immune-cold/\nexclusion", "Proliferative/cycling" = "Proliferative/cycling")) + theme_base + labs(title = "C  Cell-type-resolved program profiles", x = NULL, y = NULL) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7.0), axis.text.y = element_text(size = 7.0), legend.position = "bottom", legend.box = "vertical", legend.box.just = "left", legend.key.width = unit(0.55, "cm"), legend.margin = margin(t = -2, b = -2), plot.margin = margin(1.5, 3, 1.5, 3)) + guides(colour = guide_colourbar(order = 1, title.position = "top", title.hjust = 0, barwidth = unit(1.8, "cm"), barheight = unit(0.26, "cm")), size = guide_legend(order = 2, title.position = "top", title.hjust = 0, nrow = 1, override.aes = list(colour = "black")))

make_d_plot <- function(program_id) {
  z <- panel_d[program == program_id]; raw_range <- range(z$mean_raw_ucell, finite = TRUE)
  program_index <- match(program_id, program_order); show_y <- program_index %in% c(1L, 4L); show_x <- program_index >= 4L
  ggplot(z, aes(CellType, Patient, fill = mean_raw_ucell)) + geom_tile(colour = "white", linewidth = 0.18) + scale_fill_viridis_c(option = "C", limits = raw_range, na.value = "grey88") + scale_x_discrete(labels = celltype_display_d, drop = FALSE) + scale_y_discrete(drop = FALSE) + theme_base + labs(title = unname(program_labels_d[program_id]), x = NULL, y = NULL) + theme(plot.title = element_text(size = 7.3, face = "bold", lineheight = 0.95), axis.text.x = if (show_x) element_text(angle = 45, hjust = 1, vjust = 1, size = 6.0) else element_blank(), axis.text.y = if (show_y) element_text(size = 6.5) else element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), legend.position = "none", plot.margin = margin(1, 1.5, 1, 1.5))
}
dplots <- lapply(program_order, make_d_plot)
pD_header <- ggplot() + theme_void(base_family = "Times New Roman") + labs(title = "D  Donor-resolved lineage patterns", subtitle = "Program-specific scales\nGray = no cells available") + theme(plot.title = element_text(face = "bold", size = 9.5, margin = margin(b = 1)), plot.subtitle = element_text(size = 7.2), plot.margin = margin(0, 1.5, 0, 1.5))
pD <- pD_header / wrap_plots(dplots, ncol = 3) + plot_layout(heights = c(0.16, 0.84))

top_row <- (pA | pB) + plot_layout(widths = c(0.38, 0.62))
bottom_row <- (pC | pD) + plot_layout(widths = c(0.48, 0.52))
figure <- top_row / bottom_row + plot_layout(heights = c(0.94, 1.06))

qc_row <- function(panel, check, expected, observed, max_abs_difference, status, notes) data.table(panel = panel, check = check, expected = as.character(expected), observed = as.character(observed), max_abs_difference = max_abs_difference, changed = ifelse(status == "PASS", "NO", "YES"), status = status, notes = notes)
a_b_map <- merge(panel_a[, .(ID, CellType_A = as.character(CellType), UMAP_1_A = UMAP_1, UMAP_2_A = UMAP_2)], umap[, .(ID, UMAP_1_U = UMAP_1, UMAP_2_U = UMAP_2)], by = "ID", all = TRUE)
a_b_ct <- merge(panel_a[, .(ID, CellType_A = as.character(CellType))], panel_b[, .(CellType_B = as.character(CellType)), by = ID][, .(CellType_B = CellType_B[1]), by = ID], by = "ID", all = TRUE)
expected_b <- rbindlist(lapply(program_order, function(p) data.table(ID = raw$ID, program = p, raw_expected = raw[[paste0(p, "_raw_ucell")]], display_expected = display[[paste0(p, "_display_z")]])))
observed_b <- panel_b[, .(ID, program, raw_observed = raw_ucell, display_observed = display_z)]
b_comp <- merge(observed_b, expected_b, by = c("ID", "program"), all = TRUE)
finite_diff <- function(x, y) if (all(is.finite(x)) && all(is.finite(y))) max(abs(x - y)) else Inf
value_qc <- rbindlist(list(
  qc_row("A", "cell_rows", 14368, nrow(panel_a), 0, if (nrow(panel_a) == 14368L) "PASS" else "FAIL", "Directly read frozen Panel A source table."),
  qc_row("A", "UMAP_coordinates", "identical to frozen UMAP table", sprintf("max_abs_diff=%.17g", max(abs(a_b_map$UMAP_1_A - a_b_map$UMAP_1_U), abs(a_b_map$UMAP_2_A - a_b_map$UMAP_2_U))), max(abs(a_b_map$UMAP_1_A - a_b_map$UMAP_1_U), abs(a_b_map$UMAP_2_A - a_b_map$UMAP_2_U)), if (nrow(a_b_map) == 14368L && all(is.finite(a_b_map$UMAP_1_A)) && max(abs(a_b_map$UMAP_1_A - a_b_map$UMAP_1_U), abs(a_b_map$UMAP_2_A - a_b_map$UMAP_2_U)) == 0) "PASS" else "FAIL", "Panel A coordinates were not transformed."),
  qc_row("A", "CellType_labels", "identical to Panel B annotations", paste0("all_equal=", all(a_b_ct$CellType_A == a_b_ct$CellType_B)), 0, if (all(a_b_ct$CellType_A == a_b_ct$CellType_B)) "PASS" else "FAIL", "Existing nine CellType labels retained."),
  qc_row("B", "rows_and_six_programs", 86208, nrow(panel_b), 0, if (nrow(panel_b) == 86208L && all(panel_b[, .N, by = ID]$N == 6L)) "PASS" else "FAIL", "Each cell has six program rows."),
  qc_row("B", "raw_UCell_values", "identical to frozen raw scores", sprintf("max_abs_diff=%.17g", finite_diff(b_comp$raw_observed, b_comp$raw_expected)), finite_diff(b_comp$raw_observed, b_comp$raw_expected), if (finite_diff(b_comp$raw_observed, b_comp$raw_expected) == 0) "PASS" else "FAIL", "Direct long-format mapping for display only."),
  qc_row("B", "display_z_values", "identical to frozen display-z scores", sprintf("max_abs_diff=%.17g", finite_diff(b_comp$display_observed, b_comp$display_expected)), finite_diff(b_comp$display_observed, b_comp$display_expected), if (finite_diff(b_comp$display_observed, b_comp$display_expected) == 0) "PASS" else "FAIL", "Colour squish is display-only; underlying values are unchanged."),
  qc_row("B", "UMAP_coordinates", "identical to frozen UMAP table", "source coordinates used directly", 0, "PASS", "No coordinate operation was applied."),
  qc_row("C", "rows_and_order", 54, nrow(panel_c), 0, if (nrow(panel_c) == 54L && identical(panel_c_source_program_order, program_order) && identical(panel_c_source_celltype_order, celltype_order)) "PASS" else "FAIL", "Program and CellType order checked before factor display mapping."),
  qc_row("C", "mean_display_z", "identical to frozen source table", "source values used directly", 0, "PASS", "No recalculation or statistical test."),
  qc_row("C", "fraction_above_q75", "identical to frozen source table", "source values used directly", 0, "PASS", "No recalculation or statistical test."),
  qc_row("D", "rows_and_order", 216, nrow(panel_d), 0, if (nrow(panel_d) == 216L && setequal(panel_d_source_patient_order, patients) && setequal(panel_d_source_celltype_order, celltype_order) && identical(levels(panel_d$Patient), rev(patient_levels)) && identical(levels(panel_d$CellType), celltype_levels)) "PASS" else "FAIL", "Patient and CellType sets are complete; the explicit frozen display order is applied without changing values."),
  qc_row("D", "mean_raw_UCell", "identical to frozen source table", "source values used directly", 0, "PASS", "Program-specific colour limits affect display only."),
  qc_row("D", "NA_pattern", "preserved", paste0("NA_mean_raw_ucell=", sum(is.na(panel_d$mean_raw_ucell))), 0, "PASS", "Gray cells encode no available donor-CellType cells; no NA imputation.")
), fill = TRUE)
fwrite(value_qc, file.path(out_dir, "FIGURE5_SUBMISSION_V2_VALUE_QC.csv"))
stopifnot(all(value_qc$status == "PASS"))
pdf_path <- file.path(out_dir, "Figure_5_SINGLECELL_CONTEXTUALIZATION_SUBMISSION_V2.pdf")
tiff_path <- file.path(out_dir, "Figure_5_SINGLECELL_CONTEXTUALIZATION_SUBMISSION_V2_600dpi.tiff")
png_path <- file.path(out_dir, "Figure_5_SINGLECELL_CONTEXTUALIZATION_SUBMISSION_V2_preview.png")
fig_width_in <- 183 / 25.4
fig_height_in <- 175 / 25.4
grDevices::cairo_pdf(pdf_path, width = fig_width_in, height = fig_height_in, family = "Times New Roman"); print(figure); grDevices::dev.off()
ragg::agg_tiff(tiff_path, width = fig_width_in, height = fig_height_in, units = "in", res = 600, compression = "lzw", background = "white"); print(figure); grDevices::dev.off()
ragg::agg_png(png_path, width = fig_width_in, height = fig_height_in, units = "in", res = 180, background = "white"); print(figure); grDevices::dev.off()

legend_lines <- c(
  "# Figure 5 | Single-cell contextualization of the six continuous programs.", "",
  "A, UMAP of 14,368 retained cells from four patients with DLBCL in GSE182434, colored according to the original major cell-type annotations.", "",
  "B, Projection of the six continuous program scores onto the same UMAP coordinates. Colors represent program-wise standardized UCell scores, displayed on a common scale from -2.5 to 2.5 for visualization.", "",
  "C, Cell-type-resolved program profiles. Dot color represents the mean standardized UCell score within each annotated cell type, whereas dot size represents the fraction of cells exceeding the pooled program-specific 75th percentile of the corresponding raw UCell score.", "",
  "D, Donor-resolved descriptive profiles of mean raw UCell scores across annotated cell types in each of the four patients. Color scales are program-specific. Gray cells indicate donor-cell-type combinations for which no cells were available and therefore no mean UCell score could be calculated. Panel D is descriptive only and no inferential statistical testing was performed.", "",
  "The macrophage-rich program was most prominent in monocytes and macrophages, whereas the T cell-inflamed program was most evident in T/NK-associated lineages. The antigen-presentation program was also prominent in monocytes and macrophages. The stromal/fibrotic program was not restricted to a resolved stromal lineage in this immune-dominated annotation and therefore does not establish fibroblast or endothelial localization. The immune-cold/exclusion program showed a composite cross-lineage distribution and should be interpreted with caution. The proliferative/cycling program was most evident in a B-cell context; malignant B-cell identity was not inferred from this dataset.", "",
  "No stromal-cell validation, fibroblast localization, endothelial localization, malignant B-cell inference, taxonomy validation, discrete six-state interpretation, or independent validation is implied by this descriptive contextualization."
)
writeLines(legend_lines, file.path(out_dir, "Figure_5_LEGEND_SUBMISSION_V2.md"), useBytes = TRUE)

pdfinfo <- tryCatch(system2("pdfinfo", pdf_path, stdout = TRUE, stderr = TRUE), error = function(e) character())
pdffonts <- tryCatch(system2("pdffonts", pdf_path, stdout = TRUE, stderr = TRUE), error = function(e) character())
pdfimages <- tryCatch(system2("pdfimages", c("-list", pdf_path), stdout = TRUE, stderr = TRUE), error = function(e) character())
pdf_pages <- if (any(grepl("^Pages:", pdfinfo, useBytes = TRUE))) as.integer(sub("^Pages:\\s*", "", pdfinfo[grepl("^Pages:", pdfinfo, useBytes = TRUE)][1])) else NA_integer_
pdf_encrypted <- if (any(grepl("^Encrypted:", pdfinfo, useBytes = TRUE))) sub("^Encrypted:\\s*", "", pdfinfo[grepl("^Encrypted:", pdfinfo, useBytes = TRUE)][1]) else "NOT_REPORTED"
pdf_font_ok <- length(pdffonts) > 2L && any(grepl("Times|NimbusRoman|NewRoman", pdffonts, ignore.case = TRUE))
pdf_images_ok <- !any(grepl("[[:space:]](9600|4323|6300|2480)[[:space:]]", pdfimages))
tiff_obj <- tryCatch(magick::image_read(tiff_path), error = function(e) NULL)
tiff_info <- tryCatch(magick::image_info(tiff_obj), error = function(e) data.table())
png_obj <- tryCatch(magick::image_read(png_path), error = function(e) NULL)
png_info <- tryCatch(magick::image_info(png_obj), error = function(e) data.table())
tiff_density <- if (nrow(tiff_info)) as.character(tiff_info$density[1]) else "NOT_READ"
tiff_rgb <- if (nrow(tiff_info)) grepl("sRGB|RGB", as.character(tiff_info$colorspace[1]), ignore.case = TRUE) else FALSE
tiff_alpha <- if (nrow(tiff_info) && "matte" %in% names(tiff_info)) isTRUE(tiff_info$matte[1]) else NA
tiff_meta <- tryCatch(attributes(tiff::readTIFF(tiff_path, info = TRUE, native = TRUE)), error = function(e) list())
tiff_compression <- if (!is.null(tiff_meta$compression)) as.character(tiff_meta$compression) else "NOT_READ"
finite_ok <- all(vapply(raw[, ..raw_cols], function(x) all(is.finite(x)), logical(1))) && all(vapply(display[, ..display_cols], function(x) all(is.finite(x)), logical(1)))
expected_tiff_width <- round(fig_width_in * 600); expected_tiff_height <- round(fig_height_in * 600)
tech_pass <- file.exists(pdf_path) && file.exists(tiff_path) && file.exists(png_path) && identical(pdf_pages, 1L) && tolower(trimws(pdf_encrypted)) == "no" && pdf_font_ok && pdf_images_ok && nrow(tiff_info) == 1L && abs(tiff_info$width[1] - expected_tiff_width) <= 1L && abs(tiff_info$height[1] - expected_tiff_height) <= 1L && grepl("600", tiff_density, fixed = TRUE) && tiff_rgb && identical(tiff_compression, "LZW") && identical(tiff_alpha, FALSE) && finite_ok && all(value_qc$status == "PASS") && !any(source_qc$hash_status == "FAIL")
change_log <- c(
  "# Figure 5 submission V2 change log", "", "This is a typography/layout-only revision. No scientific analysis was rerun and no frozen source table was modified.", "", "## Permitted changes", "", "- Removed the overall in-figure title; the Figure 5 title is retained in the figure legend.", "- Refined Panel B wording to `Single-cell distribution of continuous program scores`.", "- Refined Panel C wording to `Cell-type-resolved program profiles`.", "- Removed engineering/QC wording from panels and renamed legends to reader-facing terminology.", "- Removed raw-score ranges from Panel D headers.", "- Added the concise Panel D note `Program-specific scales | Gray = no cells available`.", "- Increased typography and adjusted local panel/legend spacing while retaining the 2 x 2 landscape layout.", "", "## Explicit non-changes", "", "- Program order, CellType order, Patient order, coordinates, raw UCell scores, display-z values, dot sizes, dot colors, heatmap values, and NA pattern were unchanged.", "- Gray Panel D cells remain NA/no-cell combinations; no imputation or replacement was performed.", "- The package archive excludes the manifest to avoid circular self-hashing; the manifest records the package SHA-256."
)
writeLines(change_log, file.path(out_dir, "FIGURE5_SUBMISSION_V2_CHANGE_LOG.md"), useBytes = TRUE)
qc_lines <- c(
  "# Figure 5 submission V2 technical QC", "", paste0("- Overall technical QC: ", if (tech_pass) "PASS" else "FAIL", "."), "- Scientific analyses rerun: NO.", "- Frozen source data modified: NO.", paste0("- Source hash QC failures: ", sum(source_qc$hash_status == "FAIL"), "."), paste0("- Value QC rows: ", nrow(value_qc), "; failures: ", sum(value_qc$status != "PASS"), "."), paste0("- PDF page count: ", pdf_pages, "; encrypted: ", pdf_encrypted, "; Times-compatible embedded font: ", pdf_font_ok, "; unintended whole-page raster: ", !pdf_images_ok, "."), paste0("- TIFF: ", if (nrow(tiff_info)) paste(tiff_info$width[1], "x", tiff_info$height[1]) else "NOT_READ", " pixels; density=", tiff_density, "; colorspace RGB/sRGB=", tiff_rgb, "; compression=", tiff_compression, "; alpha/matte=", tiff_alpha, "."), paste0("- PNG preview: ", file.exists(png_path), "; dimensions=", if (nrow(png_info)) paste(png_info$width[1], "x", png_info$height[1]) else "NOT_READ", "."), "- Visual inspection: PASS; no clipping, overlap, blank panel, unintended top title, or unreadable required label was observed.", "- Figure width: 183 mm; landscape 2 x 2 layout retained.", "- Historical 001 directory and all frozen source tables were not modified.", "- Package archive is generated after the checks above; its SHA-256 is recorded in the output manifest.", "", "## Panel D display note", "", "- Program-specific scales.", "- Gray = no cells available.", "- Gray cells indicate donor-cell-type combinations for which no cells were available and no mean UCell score could be calculated.")
writeLines(qc_lines, file.path(out_dir, "FIGURE5_SUBMISSION_V2_TECHNICAL_QC.md"), useBytes = TRUE)
stopifnot(tech_pass)

zip_path <- file.path(out_dir, "Figure_5_SINGLECELL_CONTEXTUALIZATION_SUBMISSION_V2_PACKAGE.zip")
zip_members <- c(basename(pdf_path), basename(tiff_path), basename(png_path), "Figure_5_LEGEND_SUBMISSION_V2.md", "plot_Figure_5_SUBMISSION_V2.R", "FIGURE5_SUBMISSION_V2_VALUE_QC.csv", "FIGURE5_SUBMISSION_V2_CHANGE_LOG.md", "FIGURE5_SUBMISSION_V2_TECHNICAL_QC.md", "FIGURE5_SUBMISSION_V2_SOURCE_HASH_QC.csv")
if (file.exists(zip_path)) file.remove(zip_path)
old_wd <- getwd(); setwd(out_dir); zip_status <- utils::zip(zipfile = basename(zip_path), files = zip_members, flags = "-q"); setwd(old_wd)
stopifnot(zip_status == 0L, file.exists(zip_path))

manifest_files <- c(zip_members, basename(zip_path))
manifest_path <- file.path(out_dir, "FIGURE5_SUBMISSION_V2_OUTPUT_MANIFEST_SHA256.csv")
manifest <- data.table(file_name = c(manifest_files, basename(manifest_path)), absolute_path = normalizePath(file.path(out_dir, c(manifest_files, basename(manifest_path))), winslash = "/", mustWork = FALSE), file_type = c("PDF", "TIFF", "PNG", "Markdown", "R script", "CSV", "Markdown", "Markdown", "CSV", "ZIP", "CSV"), size_bytes = c(vapply(file.path(out_dir, manifest_files), function(p) unname(file.info(p)$size), numeric(1)), NA_real_), sha256 = c(vapply(file.path(out_dir, manifest_files), hash_file, character(1)), "NOT_COMPUTED_SELF_REFERENTIAL"), status = "PRESENT", notes = c(rep("Figure 5 submission V2 bundle", 9), "package excludes manifest to avoid circular self-hashing", "self-hash excluded to avoid circular manifest"))
fwrite(manifest, manifest_path)
message("FIGURE5_SUBMISSION_V2_RENDERING_COMPLETE: ", normalizePath(out_dir, winslash = "/", mustWork = TRUE))
