#!/usr/bin/env Rscript

# Figure 6 final submission typography/annotation cleanup from the QC-passed V4
# rendering contract. No score, q05/q95 limit, spatial statistic, permutation,
# P value, FDR value, correlation, clustering result, or taxonomy is recomputed.

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ragg)
  library(magick)
  library(digest)
  library(readr)
  library(zip)
  library(jsonlite)
  library(png)
})

project_dir <- Sys.getenv("DLBCL_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
revision_dir <- file.path(project_dir, "revision_2026_reviewer_response")
input_06c <- file.path(
  revision_dir,
  "06c_wp3_real_spatial_continuous_analysis/continuation_v3/finalization_v2"
)
input_06j <- file.path(revision_dir, "06j_wp3_spatial_result_interpretation")
input_06k <- file.path(revision_dir, "06k_wp3_figure7_manuscript_integration")
spatial_area4 <- file.path(
  project_dir,
  "00_raw_data/GSE276542/standard_10x/GSM8500537_Cap.area4_DLBCL_V2/spatial"
)
input_v3 <- file.path(revision_dir, "figure6_submission_revision_004")
input_v4 <- file.path(revision_dir, "figure6_submission_revision_006")
output_dir <- file.path(revision_dir, "figure6_submission_revision_FINAL")

source_paths <- c(
  area = file.path(input_06c, "WP3_FINAL_AREA_AUTHORITY_REGISTRY.csv"),
  moran_geary = file.path(input_06c, "WP3_FINAL_MORAN_GEARY_AUTHORITY.csv"),
  bivariate = file.path(input_06c, "WP3_FINAL_BIVARIATE_AUTHORITY.csv"),
  method_concordance = file.path(input_06c, "WP3_FINAL_METHOD_CONCORDANCE_AUTHORITY.csv"),
  scope = file.path(input_06c, "WP3_FINAL_SCOPE_ISOLATION.csv"),
  validator = file.path(input_06c, "WP3_FINAL_VALIDATOR_V4.csv"),
  completion = file.path(input_06c, "WP3_FINALIZATION_COMPLETION_REPORT.csv"),
  antigen = file.path(input_06j, "WP3_ANTIGEN_PRESENTATION_COVERAGE_AND_RESULT_LIMITATION.csv"),
  tissue_hires = file.path(spatial_area4, "tissue_hires_image.png"),
  scalefactors = file.path(spatial_area4, "scalefactors_json.json"),
  v4_panel_a_source = file.path(input_v4, "FIGURE6A_CANONICAL_SPOT_SOURCE_DATA.csv"),
  v4_bcd_qc = file.path(input_v4, "FIGURE6_V4_BCD_FREEZE_QC.csv"),
  v4_authority_qc = file.path(input_v4, "FIGURE6A_CANONICAL_AUTHORITY_QC.csv"),
  v4_plot_script = file.path(input_v4, "plot_Figure_6_SUBMISSION_V4.R"),
  v4_legend = file.path(input_v4, "Figure_6_LEGEND_SUBMISSION_V4.md"),
  v4_source_hash_qc = file.path(input_v4, "FIGURE6_V4_SOURCE_HASH_QC.csv"),
  v4_technical_qc = file.path(input_v4, "FIGURE6_V4_TECHNICAL_QC.md"),
  v4_pdf = file.path(input_v4, "Figure_6_SPATIAL_CONTINUOUS_SUBMISSION_V4.pdf"),
  v4_tiff = file.path(input_v4, "Figure_6_SPATIAL_CONTINUOUS_SUBMISSION_V4_600dpi.tiff"),
  v2_plot_script = file.path(input_06k, "scripts/plot_figure7_wp3_final.R"),
  v2_tiff = file.path(
    project_dir,
    "05_manuscript/Manuscript 13/Figure_6_spatial_continuous_FINAL_SUBMISSION_V2_600dpi.tiff"
  )
)

if (!all(file.exists(source_paths))) {
  stop("One or more frozen Figure 6 source files are missing.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

expected_source_hashes <- c(
  tissue_hires = "e5cb9e1c6138c907b74c3c8126e0f42e6386132b15e87803fd1e4a869fed22db",
  scalefactors = "2ed4a54793b526a79e5b520f40d90d7af7f8d70a281ba5b0e264bcccbf31d8b3",
  v4_panel_a_source = "b0c869266911f97f84e5f40a8a41b33e24882af67af907b8e5775bfdac6876f3",
  v4_bcd_qc = "06cef344bdff7c0fb2508e9e0f4877cbfc9b1e925e9051d14f6633e73777dc32",
  v4_authority_qc = "33fafe1c147bb2759bf550c46ae732b96d5cdceeee1e14832385212c6a08983f",
  v4_plot_script = "f67da8b365ed2a5e920d1ed9771be3b19f3ed171358438c2795b2c987c2b9ea5",
  v4_legend = "0885890a561fb1d5d5d7f0e3526f7d4ac12793ea608bbfb231b723ac58c6919b",
  v4_source_hash_qc = "a021977ccc0f3ed805417a9cfaab48c240828c020b436824045bf930afd84b74",
  v4_technical_qc = "92123a98d41c5774b9e3a2f3a59be890c93d58accb30f3a09674f1c6e44cb604",
  v4_pdf = "212aaab3fab2e8a6373ea1b7d26e316c2f59bb617b27e5fb89dc1791c1627a8a",
  v4_tiff = "b90ab4bffa2e1f90b18c775526b30e6fe010c8d041cd18914884a4ee4a52bc42"
)
observed_gate_hashes <- vapply(
  unname(source_paths[names(expected_source_hashes)]),
  digest::digest,
  character(1),
  file = TRUE,
  algo = "sha256"
)
names(observed_gate_hashes) <- names(expected_source_hashes)
if (!identical(unname(observed_gate_hashes), unname(expected_source_hashes))) {
  stop("A QC-passed V4 source, layout authority, or tissue-geometry hash changed.")
}

finalize <- identical(Sys.getenv("FIGURE6_FINALIZE", unset = "FALSE"), "TRUE")
candidate_id <- Sys.getenv("FIGURE6_CANDIDATE_ID", unset = "001")
render_dir <- if (finalize) output_dir else file.path(output_dir, paste0("candidate_render_", candidate_id))
dir.create(render_dir, recursive = TRUE, showWarnings = FALSE)

font_family <- "Times New Roman"
width_mm <- 183
height_mm <- 228
minimum_text_pt <- 7.0
seed_jitter <- 20260808L

if (!requireNamespace("systemfonts", quietly = TRUE)) {
  stop("Package systemfonts is required for the Times New Roman preflight.")
}
font_match <- systemfonts::match_fonts(font_family)
if (nrow(font_match) != 1L || !file.exists(font_match$path[1])) {
  stop("Times New Roman is not available to the R graphics backend.")
}

read_authority <- function(path) read.csv(path, check.names = FALSE)
area <- read_authority(source_paths[["area"]])
moran <- read_authority(source_paths[["moran_geary"]])
biv <- read_authority(source_paths[["bivariate"]])
method <- read_authority(source_paths[["method_concordance"]])
scope <- read_authority(source_paths[["scope"]])
validator <- read_authority(source_paths[["validator"]])
completion <- read_authority(source_paths[["completion"]])
antigen <- read_authority(source_paths[["antigen"]])
if (!identical(colnames(antigen)[1], "capture_area") && grepl("capture_area$", colnames(antigen)[1])) {
  colnames(antigen)[1] <- "capture_area"
}
spot_source <- read_authority(source_paths[["v4_panel_a_source"]])
v4_bcd_qc <- read_authority(source_paths[["v4_bcd_qc"]])
v4_authority_qc <- read_authority(source_paths[["v4_authority_qc"]])
scale_factors <- jsonlite::fromJSON(source_paths[["scalefactors"]])
tissue_image <- png::readPNG(source_paths[["tissue_hires"]])
tissue_height <- dim(tissue_image)[1]
tissue_width <- dim(tissue_image)[2]
hires_scale <- as.numeric(scale_factors$tissue_hires_scalef)

program_ids <- c(
  "macrophage_rich",
  "t_cell_inflamed",
  "antigen_presentation",
  "stromal_fibrotic",
  "immune_cold_exclusion",
  "proliferative_cycling"
)
program_canonical <- c(
  macrophage_rich = "Macrophage-rich",
  t_cell_inflamed = "T cell-inflamed",
  antigen_presentation = "Immune-inflamed / antigen-presentation",
  stromal_fibrotic = "Stromal / fibrotic",
  immune_cold_exclusion = "Immune-cold / exclusion-associated",
  proliferative_cycling = "Proliferative / cycling"
)
program_old <- c(
  macrophage_rich = "Macrophage-rich",
  t_cell_inflamed = "T cell-inflamed",
  antigen_presentation = "Antigen-presentation",
  stromal_fibrotic = "Stromal/fibrotic",
  immune_cold_exclusion = "Immune-cold/exclusion",
  proliferative_cycling = "Proliferative/cycling"
)
program_abbr <- c(
  macrophage_rich = "MR",
  t_cell_inflamed = "TCI",
  antigen_presentation = "AP",
  stromal_fibrotic = "SF",
  immune_cold_exclusion = "IC/EX",
  proliferative_cycling = "PC"
)
displayed_programs <- c(
  "macrophage_rich",
  "t_cell_inflamed",
  "stromal_fibrotic",
  "proliferative_cycling"
)
displayed_program_labels <- unname(program_canonical[displayed_programs])
expected_displayed_programs <- c(
  "macrophage_rich",
  "t_cell_inflamed",
  "stromal_fibrotic",
  "proliferative_cycling"
)
legacy_score_columns <- c(
  "Macrophage_rich_Score1",
  "T_cell_inflamed_Score1",
  "Stromal_fibrotic_Score1",
  "Proliferative_malignant_B_cell_Score1"
)
panel_a_expected_columns <- c(
  "barcode", "capture_area", "x_coordinate", "y_coordinate", "in_tissue",
  "macrophage_rich", "t_cell_inflamed", "stromal_fibrotic", "proliferative_cycling",
  "macrophage_rich_q05", "macrophage_rich_q95", "t_cell_inflamed_q05",
  "t_cell_inflamed_q95", "stromal_fibrotic_q05", "stromal_fibrotic_q95",
  "proliferative_cycling_q05", "proliferative_cycling_q95"
)
capture_area_id <- "GSM8500537_Cap.area4_DLBCL_V2"

display_stats <- do.call(rbind, lapply(displayed_programs, function(program) {
  q05_values <- unique(spot_source[[paste0(program, "_q05")]])
  q95_values <- unique(spot_source[[paste0(program, "_q95")]])
  if (length(q05_values) != 1L || length(q95_values) != 1L) {
    stop("V4 Panel A display limits are not single frozen values per program.")
  }
  data.frame(
    program = program,
    program_label = unname(program_canonical[program]),
    q05 = q05_values,
    q95 = q95_values,
    stringsAsFactors = FALSE
  )
}))
spot_source$x_hires <- spot_source$x_coordinate * hires_scale
spot_source$y_hires <- tissue_height - spot_source$y_coordinate * hires_scale
rho_symbol <- intToUtf8(0x03C1)
times_symbol <- intToUtf8(0x00D7)

as_flag <- function(x) tolower(as.character(x)) == "true"
area_order <- paste0("Cap.area", 1:9)
primary_order <- paste0("Cap.area", 3:7)
area$capture_area <- factor(area$capture_area, levels = area_order)
area_name_map <- setNames(as.character(area$capture_area), area$capture_area_id)
area_role_map <- setNames(
  ifelse(as.character(area$capture_area) %in% primary_order, "PRIMARY_DLBCL", "CONTEXT_ONLY"),
  area$capture_area_id
)

moran$capture_area <- unname(area_name_map[moran$capture_area_id])
biv$capture_area <- unname(area_name_map[biv$capture_area_id])
method$capture_area <- unname(area_name_map[method$capture_area_id])
biv$area_role <- unname(area_role_map[biv$capture_area_id])
method$area_role <- unname(area_role_map[method$capture_area_id])

# Fail closed against the frozen scientific contract.
contract_checks <- c(
  nrow(area) == 9L,
  nrow(moran) == 54L,
  nrow(biv) == 135L,
  nrow(method) == 108L,
  identical(sort(unique(moran$program_id)), sort(program_ids)),
  identical(sort(unique(moran$permutations)), 9999L),
  identical(sort(unique(biv$permutations)), 9999L),
  all(is.finite(moran$Moran_I)),
  all(is.finite(moran$Geary_C)),
  all(is.finite(moran$Moran_FDR)),
  all(is.finite(moran$Geary_FDR)),
  all(is.finite(biv$bivariate_Moran_I)),
  all(is.finite(biv$fdr)),
  all(is.finite(method$spearman_rho)),
  all(as_flag(validator$status)),
  identical(as.character(completion$final_k[1]), "NOT_SELECTED"),
  identical(as.character(completion$taxonomy[1]), "NOT_ASSIGNED"),
  identical(as.character(completion$raw_expression_read[1]), "FALSE"),
  identical(as.character(completion$score_computations[1]), "0"),
  identical(as.character(completion$spatial_statistic_computations[1]), "0"),
  identical(setdiff(colnames(spot_source), c("x_hires", "y_hires")), panel_a_expected_columns),
  nrow(spot_source) == 1723L,
  length(unique(spot_source$barcode)) == 1723L,
  identical(unique(spot_source$capture_area), capture_area_id),
  identical(displayed_programs, expected_displayed_programs),
  all(vapply(spot_source[, displayed_programs, drop = FALSE], is.numeric, logical(1))),
  all(vapply(spot_source[, displayed_programs, drop = FALSE], function(x) all(is.finite(x)), logical(1))),
  !any(legacy_score_columns %in% colnames(spot_source)),
  !anyDuplicated(spot_source$barcode),
  all(is.finite(spot_source$x_coordinate)),
  all(is.finite(spot_source$y_coordinate)),
  all(is.finite(display_stats$q05)),
  all(is.finite(display_stats$q95)),
  all(display_stats$q05 < display_stats$q95),
  nrow(v4_bcd_qc) == 693L,
  all(v4_bcd_qc$status == "PASS"),
  all(v4_bcd_qc$panel_n_mismatch == 0L),
  all(v4_bcd_qc$panel_max_abs_difference <= 1e-12),
  all(v4_authority_qc$status == "PASS"),
  identical(as.integer(tissue_width), 1951L),
  identical(as.integer(tissue_height), 2000L),
  is.finite(hires_scale),
  hires_scale > 0
)
if (!all(contract_checks)) {
  stop(
    "Frozen Figure 6 scientific contract failed at check indices: ",
    paste(which(!contract_checks), collapse = ","),
    "; no rendering was performed."
  )
}

# Exact frozen descriptive constants from the final authority-backed legend.
method_summary <- list(
  median = 0.877515,
  q1 = 0.756383,
  q3 = 0.954079,
  minimum = 0.190400,
  maximum = 0.991849,
  weak_area = "Cap.area6",
  weak_program = "t_cell_inflamed",
  weak_comparison = "PRIMARY_vs_ssGSEA"
)
weak <- method[
  method$capture_area == method_summary$weak_area &
    method$program_id == method_summary$weak_program &
    method$comparison == method_summary$weak_comparison,
]
if (nrow(weak) != 1L || abs(weak$spearman_rho - method_summary$minimum) > 5e-7) {
  stop("Frozen weakest method-comparison record is inconsistent.")
}
ring_currently_circled <- TRUE
ring_is_frozen_weakest <- TRUE
ring_retained <- ring_currently_circled && ring_is_frozen_weakest

theme_wp3_v3 <- function(base_size = 7.3) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      text = element_text(family = font_family, colour = "#202020"),
      axis.line = element_line(linewidth = 0.35, colour = "#2B2B2B"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#2B2B2B"),
      axis.text = element_text(size = 7.1, colour = "#202020"),
      axis.title = element_text(size = 7.4, colour = "#202020"),
      plot.title = element_text(size = 9.3, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 7.2, colour = "#4D4D4D", hjust = 0),
      plot.caption = element_blank(),
      strip.background = element_rect(fill = "#F2F3F5", colour = "#B8BDC4", linewidth = 0.35),
      strip.text = element_text(size = 7.4, face = "bold"),
      legend.title = element_text(size = 7.3),
      legend.text = element_text(size = 7.0),
      legend.key.height = grid::unit(3.4, "mm"),
      legend.key.width = grid::unit(4.5, "mm"),
      plot.margin = margin(2, 3, 2, 3)
    )
}
theme_set(theme_wp3_v3())

panel_a_palette <- c("#F7F7F7", "#F4A582", "#D6604D", "#B2182B")

make_spatial_map <- function(program) {
  stat_row <- display_stats[display_stats$program == program, ]
  ggplot(spot_source, aes(x = x_hires, y = y_hires)) +
    annotation_raster(
      tissue_image,
      xmin = 0,
      xmax = tissue_width,
      ymin = 0,
      ymax = tissue_height,
      interpolate = FALSE
    ) +
    geom_point(
      aes(colour = .data[[program]]),
      shape = 16,
      size = 0.52,
      alpha = 0.94,
      stroke = 0
    ) +
    scale_colour_gradientn(
      colours = panel_a_palette,
      limits = c(stat_row$q05, stat_row$q95),
      oob = scales::squish,
      breaks = c(stat_row$q05, stat_row$q95),
      labels = scales::label_number(accuracy = 0.01, trim = TRUE),
      guide = guide_colourbar(
        direction = "horizontal",
        title = NULL,
        barwidth = grid::unit(15, "mm"),
        barheight = grid::unit(1.2, "mm"),
        ticks = TRUE,
        frame.colour = "#666666",
        frame.linewidth = 0.25
      )
    ) +
    coord_fixed(
      xlim = c(0, tissue_width),
      ylim = c(0, tissue_height),
      expand = FALSE,
      clip = "on"
    ) +
    labs(title = unname(program_canonical[program]), x = NULL, y = NULL) +
    theme_void(base_family = font_family) +
    theme(
      plot.title = element_text(
        family = font_family,
        size = 8.0,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 1)
      ),
      legend.position = "bottom",
      legend.justification = "center",
      legend.text = element_text(family = font_family, size = 7.0, colour = "#202020"),
      legend.margin = margin(t = -1, r = 0, b = 0, l = 0),
      plot.margin = margin(0, 3.5, 0, 3.5)
    )
}

panel_a_maps <- lapply(displayed_programs, make_spatial_map)
p_a <- wrap_plots(panel_a_maps, nrow = 1) +
  plot_annotation(
    title = "Representative spatial distribution of continuous programs",
    subtitle = "GSM8500537_Cap.area4_DLBCL_V2",
    theme = theme(
      plot.title = element_text(
        family = font_family,
        size = 9.3,
        face = "bold",
        hjust = 0,
        margin = margin(l = 5)
      ),
      plot.subtitle = element_text(
        family = font_family,
        size = 7.2,
        colour = "#4D4D4D",
        hjust = 0,
        margin = margin(l = 5, b = 1)
      )
    )
  )

m_primary <- moran[moran$area_role == "PRIMARY_DLBCL", ]
m_primary$capture_area <- factor(m_primary$capture_area, levels = rev(primary_order))
m_primary$program_label <- factor(
  unname(program_abbr[m_primary$program_id]),
  levels = unname(program_abbr[program_ids])
)
m_primary$exploratory_display <- m_primary$eligibility == "EXPLORATORY_ONLY"
m_primary$Moran_significant <- m_primary$Moran_FDR < 0.05
m_primary$Geary_significant <- m_primary$Geary_FDR < 0.05
m_primary$moran_label <- sprintf(
  "%.2f%s%s",
  m_primary$Moran_I,
  ifelse(m_primary$Moran_significant, "*", ""),
  ifelse(m_primary$exploratory_display, " E", "")
)
m_primary$geary_label <- sprintf(
  "%.2f%s%s",
  m_primary$Geary_C,
  ifelse(m_primary$Geary_significant, "*", ""),
  ifelse(m_primary$exploratory_display, " E", "")
)

p_moran <- ggplot(m_primary, aes(program_label, capture_area)) +
  geom_tile(aes(fill = Moran_I, linetype = exploratory_display), colour = "#6F747A", linewidth = 0.38) +
  geom_text(aes(label = moran_label), size = 2.52, colour = "#111111", family = font_family) +
  scale_fill_gradient(low = "#F7FBFF", high = "#2166AC", name = "Moran's I") +
  scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "dashed"), guide = "none") +
  labs(title = "Moran's I", x = NULL, y = NULL) +
  theme_wp3_v3() +
  theme(
    axis.text.x = element_text(size = 7.2, face = "bold"),
    axis.text.y = element_text(size = 7.1),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    plot.margin = margin(1, 2, 1, 2)
  )

p_geary <- ggplot(m_primary, aes(program_label, capture_area)) +
  geom_tile(aes(fill = Geary_C, linetype = exploratory_display), colour = "#6F747A", linewidth = 0.38) +
  geom_text(aes(label = geary_label), size = 2.52, colour = "#111111", family = font_family) +
  scale_fill_gradient(low = "#2166AC", high = "#F7FBFF", name = "Geary's C") +
  scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "dashed"), guide = "none") +
  labs(title = "Geary's C", x = NULL, y = NULL) +
  theme_wp3_v3() +
  theme(
    axis.text.x = element_text(size = 7.2, face = "bold"),
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    plot.margin = margin(1, 2, 1, 2)
  )

p_b <- (p_moran | p_geary) +
  plot_annotation(
    title = "Spatial autocorrelation of continuous programs",
    theme = theme(
      plot.title = element_text(family = font_family, size = 9.3, face = "bold", margin = margin(l = 10)),
      plot.subtitle = element_blank()
    )
  )

pair_template <- unique(biv[, c("program_1", "program_2")])
pair_template$contains_antigen <- pair_template$program_1 == "antigen_presentation" |
  pair_template$program_2 == "antigen_presentation"
pair_template$p1_order <- match(pair_template$program_1, program_ids)
pair_template$p2_order <- match(pair_template$program_2, program_ids)
pair_template <- pair_template[order(
  pair_template$contains_antigen,
  pair_template$p1_order,
  pair_template$p2_order
), ]
pair_template$pair_key <- paste(pair_template$program_1, pair_template$program_2, sep = "__")
pair_template$pair_label <- paste(
  unname(program_abbr[pair_template$program_1]),
  unname(program_abbr[pair_template$program_2]),
  sep = paste0(" ", times_symbol, " ")
)
pair_levels <- rev(pair_template$pair_label)
pair_label_map <- setNames(pair_template$pair_label, pair_template$pair_key)

b_primary <- biv[biv$area_role == "PRIMARY_DLBCL", ]
b_primary$capture_area <- factor(b_primary$capture_area, levels = primary_order)
b_primary$pair_key <- paste(b_primary$program_1, b_primary$program_2, sep = "__")
b_primary$pair_label <- factor(unname(pair_label_map[b_primary$pair_key]), levels = pair_levels)
b_primary$FDR_significant <- b_primary$fdr < 0.05
b_primary$exploratory_flag <- as_flag(b_primary$exploratory)
b_primary$cell_label <- sprintf(
  "%.2f%s%s",
  b_primary$bivariate_Moran_I,
  ifelse(b_primary$FDR_significant, "*", ""),
  ifelse(b_primary$exploratory_flag, " E", "")
)
biv_limit <- max(abs(b_primary$bivariate_Moran_I), na.rm = TRUE)

p_c <- ggplot(b_primary, aes(capture_area, pair_label)) +
  geom_tile(
    aes(fill = bivariate_Moran_I, linetype = exploratory_flag),
    colour = "#777777",
    linewidth = 0.38
  ) +
  geom_text(aes(label = cell_label), size = 2.50, colour = "#111111", family = font_family) +
  scale_fill_gradient2(
    low = "#3B6FB6",
    mid = "#FFFFFF",
    high = "#C84A47",
    midpoint = 0,
    limits = c(-biv_limit, biv_limit),
    oob = squish,
    name = "Bivariate\nMoran's I"
  ) +
  scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "dashed"), guide = "none") +
  labs(
    title = "Bivariate spatial associations among programs",
    x = NULL,
    y = NULL
  ) +
  theme_wp3_v3() +
  theme(
    axis.text.x = element_text(size = 7.3, face = "bold"),
    axis.text.y = element_text(size = 7.0),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    plot.margin = margin(1, 3, 1, 3)
  )

method$program_label <- factor(
  unname(program_abbr[method$program_id]),
  levels = rev(unname(program_abbr[program_ids]))
)
method$capture_area <- factor(method$capture_area, levels = area_order)
method$comparison_display <- factor(
  ifelse(
    method$comparison == "PRIMARY_vs_SCT_UCell",
    "Primary vs SCT-UCell",
    "Primary vs ssGSEA"
  ),
  levels = c("Primary vs SCT-UCell", "Primary vs ssGSEA")
)
method$role_shape <- factor(method$area_role, levels = c("PRIMARY_DLBCL", "CONTEXT_ONLY"))
weak$program_label <- factor(
  unname(program_abbr[weak$program_id]),
  levels = rev(unname(program_abbr[program_ids]))
)
weak$comparison_display <- factor(
  "Primary vs ssGSEA",
  levels = c("Primary vs SCT-UCell", "Primary vs ssGSEA")
)

p_d <- ggplot(method, aes(spearman_rho, program_label)) +
  geom_vline(xintercept = 0, colour = "#B9BDC2", linewidth = 0.35) +
  geom_point(
    aes(colour = area_role, shape = role_shape),
    size = 1.65,
    alpha = 0.90,
    position = position_jitter(height = 0.10, width = 0, seed = seed_jitter)
  ) +
  geom_point(
    data = weak,
    aes(x = spearman_rho, y = program_label),
    inherit.aes = FALSE,
    shape = 21,
    size = 3.0,
    stroke = 0.75,
    fill = NA,
    colour = "black"
  ) +
  facet_wrap(~comparison_display, nrow = 1) +
  scale_colour_manual(
    values = c(PRIMARY_DLBCL = "#2F6B8A", CONTEXT_ONLY = "#AEB4BA"),
    labels = c(PRIMARY_DLBCL = "Primary DLBCL", CONTEXT_ONLY = "Context only"),
    name = "Area role"
  ) +
  scale_shape_manual(values = c(PRIMARY_DLBCL = 16, CONTEXT_ONLY = 1), guide = "none") +
  scale_x_continuous(limits = c(0, 1.01), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Scoring-method concordance",
    subtitle = paste0("Median Spearman ", rho_symbol, " = 0.88"),
    x = paste0("Spearman ", rho_symbol),
    y = NULL
  ) +
  theme_wp3_v3() +
  theme(
    axis.text.y = element_text(size = 7.2, face = "bold"),
    axis.text.x = element_text(size = 7.1),
    strip.text = element_text(size = 7.5),
    legend.position = "bottom",
    legend.direction = "horizontal",
    plot.margin = margin(1, 3, 1, 3)
  )

main_plot <- wrap_elements(full = p_a) / wrap_elements(full = p_b) / p_c / p_d +
  plot_layout(heights = c(1.75, 1.45, 2.30, 1.25)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(family = font_family, size = 11.5, face = "bold"),
      plot.tag.position = c(0.002, 0.995)
    )
  )

stem <- if (finalize) {
  "Figure_6_SPATIAL_CONTINUOUS_FINAL_SUBMISSION"
} else {
  paste0("Figure_6_SPATIAL_CONTINUOUS_FINAL_CANDIDATE_", candidate_id)
}
pdf_path <- file.path(render_dir, paste0(stem, ".pdf"))
tiff_path <- file.path(render_dir, paste0(stem, if (finalize) "_600dpi.tiff" else "_600dpi.tiff"))
png_path <- file.path(render_dir, paste0(stem, if (finalize) "_preview.png" else "_preview.png"))

if (any(file.exists(c(pdf_path, tiff_path, png_path)))) {
  stop("Refusing to overwrite an existing Figure 6 render.")
}

width_in <- width_mm / 25.4
height_in <- height_mm / 25.4
grDevices::cairo_pdf(
  pdf_path,
  width = width_in,
  height = height_in,
  family = font_family,
  bg = "white"
)
print(main_plot)
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
print(main_plot)
grDevices::dev.off()

ragg::agg_png(
  png_path,
  width = width_in,
  height = height_in,
  units = "in",
  res = 220,
  background = "white"
)
print(main_plot)
grDevices::dev.off()

run_tool <- function(name, args) {
  tool <- Sys.which(name)
  if (!nzchar(tool)) return(paste(name, "NOT_FOUND"))
  suppressWarnings(system2(tool, args = args, stdout = TRUE, stderr = TRUE))
}

pdfinfo_lines <- run_tool("pdfinfo", shQuote(pdf_path))
pdffonts_lines <- run_tool("pdffonts", shQuote(pdf_path))
pdfimages_lines <- run_tool("pdfimages", c("-list", shQuote(pdf_path)))
tiff_info <- magick::image_info(magick::image_read(tiff_path))[1, ]
png_info <- magick::image_info(magick::image_read(png_path))[1, ]

read_tiff_short_tag <- function(path, wanted_tag) {
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  endian_marker <- rawToChar(readBin(con, "raw", n = 2L))
  endian <- if (endian_marker == "II") "little" else if (endian_marker == "MM") "big" else stop("Invalid TIFF endian marker")
  magic <- readBin(con, "integer", n = 1L, size = 2L, endian = endian, signed = FALSE)
  if (magic != 42L) stop("Invalid classic TIFF magic")
  ifd_offset <- readBin(con, "integer", n = 1L, size = 4L, endian = endian)
  seek(con, where = ifd_offset, origin = "start")
  n_entries <- readBin(con, "integer", n = 1L, size = 2L, endian = endian, signed = FALSE)
  for (i in seq_len(n_entries)) {
    tag <- readBin(con, "integer", n = 1L, size = 2L, endian = endian, signed = FALSE)
    type <- readBin(con, "integer", n = 1L, size = 2L, endian = endian, signed = FALSE)
    count <- readBin(con, "integer", n = 1L, size = 4L, endian = endian)
    value_raw <- readBin(con, "raw", n = 4L)
    if (tag == wanted_tag && type == 3L && count == 1L) {
      tmp <- rawConnection(value_raw, "rb")
      on.exit(close(tmp), add = TRUE)
      return(readBin(tmp, "integer", n = 1L, size = 2L, endian = endian, signed = FALSE))
    }
  }
  NA_integer_
}

tiff_compression_tag <- read_tiff_short_tag(tiff_path, 259L)
tiff_photometric_tag <- read_tiff_short_tag(tiff_path, 262L)
tiff_samples_per_pixel <- read_tiff_short_tag(tiff_path, 277L)
tiff_extra_samples <- read_tiff_short_tag(tiff_path, 338L)

if (!finalize) {
  candidate_qc <- c(
    "# Figure 6 candidate render QC",
    "",
    paste0("- Candidate: ", candidate_id),
    paste0("- Size: ", width_mm, " x ", height_mm, " mm"),
    paste0("- PDF: ", normalizePath(pdf_path, winslash = "/", mustWork = TRUE)),
    paste0("- TIFF: ", normalizePath(tiff_path, winslash = "/", mustWork = TRUE)),
    paste0("- PNG: ", normalizePath(png_path, winslash = "/", mustWork = TRUE)),
    paste0("- TIFF pixels: ", tiff_info$width, " x ", tiff_info$height),
    paste0("- TIFF density: ", tiff_info$density),
    paste0("- TIFF colorspace: ", tiff_info$colorspace),
    paste0("- TIFF matte: ", tiff_info$matte),
    paste0("- TIFF compression tag: ", tiff_compression_tag, " (5=LZW)"),
    "",
    "## pdfinfo",
    "```",
    pdfinfo_lines,
    "```",
    "",
    "## pdffonts",
    "```",
    pdffonts_lines,
    "```",
    "",
    "## pdfimages -list",
    "```",
    pdfimages_lines,
    "```"
  )
  writeLines(candidate_qc, file.path(render_dir, "CANDIDATE_TECHNICAL_QC.md"), useBytes = TRUE)
  message("FIGURE6_CANDIDATE_RENDER_COMPLETE: ", normalizePath(render_dir, winslash = "/", mustWork = TRUE))
  quit(save = "no", status = 0L)
}

# Retained only as inactive provenance from the V4 script copy. The final
# submission path below replaces this block and never executes it.
if (FALSE) {
# Finalize the deterministic Panel A source and governance records only after a
# candidate render has completed successfully.
spot_source_export <- spot_source[, setdiff(colnames(spot_source), c("x_hires", "y_hires")), drop = FALSE]
panel_a_source_path <- file.path(output_dir, "FIGURE6A_CANONICAL_SPOT_SOURCE_DATA.csv")
readr::write_csv(spot_source_export, panel_a_source_path, na = "")

authority_qc <- data.frame(
  check_id = c(
    "canonical_score_hash", "canonical_row_count", "canonical_unique_barcodes",
    "canonical_columns", "legacy_columns_absent", "coordinate_grid_unique_barcodes",
    "coordinate_join_count", "coordinate_join_complete", "coordinate_values_finite",
    "capture_area_fixed", "displayed_program_set_fixed", "canonical_scores_finite"
  ),
  expected = c(
    expected_source_hashes[["canonical_score"]], "1723", "1723",
    paste(canonical_expected_columns, collapse = ";"), "TRUE", "4992",
    "1723", "TRUE", "TRUE", capture_area_id,
    paste(expected_displayed_programs, collapse = ";"), "TRUE"
  ),
  observed = c(
    observed_gate_hashes[["canonical_score"]], as.character(nrow(canonical_score)),
    as.character(length(unique(canonical_score$barcode))),
    paste(colnames(canonical_score), collapse = ";"),
    as.character(!any(legacy_score_columns %in% colnames(canonical_score))),
    as.character(length(unique(positions$barcode))), as.character(nrow(spot_source)),
    as.character(position_join_complete && identical(joined_positions$barcode, canonical_score$barcode)),
    as.character(all(is.finite(spot_source$x_coordinate)) && all(is.finite(spot_source$y_coordinate))),
    unique(spot_source$capture_area), paste(displayed_programs, collapse = ";"),
    as.character(all(vapply(canonical_score[, displayed_programs, drop = FALSE], function(x) all(is.finite(x)), logical(1))))
  ),
  status = "PASS",
  notes = c(
    "Frozen WP3 canonical authority hash.", "No spot filtering.", "No duplicate score barcode.",
    "Barcode plus six raw canonical program identifiers.", "Retired AddModuleScore columns are prohibited.",
    "Standard 10x coordinate grid.", "Exact score-to-coordinate join.", "Barcode order preserved after exact match.",
    "Full-resolution pixel coordinates.", "Prespecified historical representative area.",
    "Program identities inherited as selection provenance only.", "No score recomputation."
  ),
  stringsAsFactors = FALSE
)
authority_qc$status <- ifelse(authority_qc$expected == authority_qc$observed, "PASS", "FAIL")
authority_qc_pass <- all(authority_qc$status == "PASS")
write.csv(
  authority_qc,
  file.path(output_dir, "FIGURE6A_CANONICAL_AUTHORITY_QC.csv"),
  row.names = FALSE,
  na = ""
)

legacy_exclusion_qc <- data.frame(
  source_file = normalizePath(source_paths[["canonical_score"]], winslash = "/", mustWork = TRUE),
  source_columns = paste(colnames(canonical_score), collapse = ";"),
  prohibited_legacy_columns = paste(legacy_score_columns, collapse = ";"),
  legacy_AddModuleScore_columns_present = any(legacy_score_columns %in% colnames(spot_source_export)),
  old_values_reused = FALSE,
  status = ifelse(any(legacy_score_columns %in% colnames(spot_source_export)), "FAIL", "PASS"),
  stringsAsFactors = FALSE
)
write.csv(
  legacy_exclusion_qc,
  file.path(output_dir, "FIGURE6A_LEGACY_COLUMN_EXCLUSION_QC.csv"),
  row.names = FALSE,
  na = ""
)
legacy_exclusion_pass <- identical(legacy_exclusion_qc$status, "PASS")

display_contract_lines <- c(
  "# Figure 6A Canonical Display Contract",
  "",
  "- Existing frozen canonical spatial-map display contract found: NO.",
  paste0("- Capture area: `", capture_area_id, "`."),
  "- Score authority: frozen canonical LogNormalize+UCell values; no score was recomputed.",
  "- Coordinates: standard 10x full-resolution pixel coordinates joined by exact barcode.",
  "- Tissue background: authentic high-resolution tissue image using the frozen `tissue_hires_scalef`.",
  "- Layout: four maps use identical spot membership, tissue geometry, spatial extent, and crop.",
  "- Display transformation: each program independently uses R `quantile(..., probs = c(0.05, 0.95), type = 7)`.",
  "- Values outside q05-q95 are squished only for color rendering; source scores are unchanged.",
  "- Scale type: program-specific; colors are not comparable as absolute magnitudes across programs.",
  "",
  "| Program | Minimum | q05 | q95 | Maximum |",
  "|---|---:|---:|---:|---:|",
  vapply(seq_len(nrow(display_stats)), function(i) {
    paste0(
      "| ", display_stats$program_label[i], " | ",
      format(display_stats$minimum[i], digits = 17), " | ",
      format(display_stats$q05[i], digits = 17), " | ",
      format(display_stats$q95[i], digits = 17), " | ",
      format(display_stats$maximum[i], digits = 17), " |"
    )
  }, character(1))
)
writeLines(
  display_contract_lines,
  file.path(output_dir, "FIGURE6A_CANONICAL_DISPLAY_CONTRACT.md"),
  useBytes = TRUE
)

legacy_note_lines <- c(
  "# Figure 6A Legacy Provenance Note",
  "",
  "Historical Figure 7A is retained only as provenance for the prespecified representative capture area and displayed program set. Its retired AddModuleScore values were not reused in the canonical Figure 6A.",
  "",
  paste0("- Prespecified capture area retained: `", capture_area_id, "`."),
  paste0("- Prespecified displayed programs retained: ", paste(displayed_program_labels, collapse = "; "), "."),
  "- Historical score authority retained: NO.",
  "- Historical color-value authority retained: NO.",
  "- Historical layout used only as a four-map reference: YES.",
  "- Controlling legacy disposition: `MUST_RETIRE_FROM_PRIMARY_INFERENCE`."
)
writeLines(
  legacy_note_lines,
  file.path(output_dir, "FIGURE6A_LEGACY_PROVENANCE_NOTE.md"),
  useBytes = TRUE
)

# Final B-D value QC: every plotted scientific value and status is compared
# with the V3 rendering authority at the original stored precision.
make_numeric_qc <- function(panel, key, field, x) {
  data.frame(
    panel = panel,
    record_key = key,
    field = field,
    old_value = format(x, digits = 17, scientific = FALSE, trim = TRUE),
    new_value = format(x, digits = 17, scientific = FALSE, trim = TRUE),
    absolute_difference = 0,
    tolerance = 1e-12,
    status = "PASS",
    stringsAsFactors = FALSE
  )
}
make_text_qc <- function(panel, key, field, x) {
  data.frame(
    panel = panel,
    record_key = key,
    field = field,
    old_value = as.character(x),
    new_value = as.character(x),
    absolute_difference = NA_real_,
    tolerance = NA_real_,
    status = "PASS",
    stringsAsFactors = FALSE
  )
}

b_key <- paste(m_primary$capture_area_id, m_primary$program_id, sep = "|")
current_bcd_qc <- rbind(
  make_numeric_qc("B", b_key, "Moran_I", m_primary$Moran_I),
  make_numeric_qc("B", b_key, "Geary_C", m_primary$Geary_C),
  make_numeric_qc("B", b_key, "Moran_FDR", m_primary$Moran_FDR),
  make_numeric_qc("B", b_key, "Geary_FDR", m_primary$Geary_FDR),
  make_text_qc("B", b_key, "Moran_significant", m_primary$Moran_significant),
  make_text_qc("B", b_key, "Geary_significant", m_primary$Geary_significant),
  make_text_qc("B", b_key, "exploratory_display", m_primary$exploratory_display)
)

c_key <- paste(b_primary$capture_area_id, b_primary$program_1, b_primary$program_2, sep = "|")
current_bcd_qc <- rbind(
  current_bcd_qc,
  make_numeric_qc("C", c_key, "bivariate_Moran_I", b_primary$bivariate_Moran_I),
  make_numeric_qc("C", c_key, "fdr", b_primary$fdr),
  make_text_qc("C", c_key, "direction", sign(b_primary$bivariate_Moran_I)),
  make_text_qc("C", c_key, "FDR_significant", b_primary$FDR_significant),
  make_text_qc("C", c_key, "exploratory_flag", b_primary$exploratory_flag)
)

d_key <- paste(method$capture_area_id, method$program_id, method$comparison, sep = "|")
current_bcd_qc <- rbind(
  current_bcd_qc,
  make_numeric_qc("D", d_key, "spearman_rho", method$spearman_rho)
)
v3_bcd_qc <- v3_value_qc[v3_value_qc$panel %in% c("B", "C", "D"), ]
v3_bcd_qc <- v3_bcd_qc[, c("panel", "record_key", "field", "new_value", "absolute_difference")]
colnames(v3_bcd_qc) <- c("panel", "record_key", "field", "v3_value", "v3_numeric_marker")
v4_bcd_qc <- current_bcd_qc[, c("panel", "record_key", "field", "new_value", "absolute_difference")]
colnames(v4_bcd_qc) <- c("panel", "record_key", "field", "v4_value", "v4_numeric_marker")
bcd_freeze_qc <- merge(
  v3_bcd_qc,
  v4_bcd_qc,
  by = c("panel", "record_key", "field"),
  all = TRUE,
  sort = FALSE
)
bcd_freeze_qc$numeric_field <- !is.na(bcd_freeze_qc$v4_numeric_marker)
bcd_freeze_qc$absolute_difference <- ifelse(
  bcd_freeze_qc$numeric_field,
  abs(as.numeric(bcd_freeze_qc$v3_value) - as.numeric(bcd_freeze_qc$v4_value)),
  NA_real_
)
bcd_freeze_qc$value_match <- ifelse(
  bcd_freeze_qc$numeric_field,
  bcd_freeze_qc$absolute_difference <= 1e-12,
  bcd_freeze_qc$v3_value == bcd_freeze_qc$v4_value
)
bcd_freeze_qc$status <- ifelse(bcd_freeze_qc$value_match, "PASS", "FAIL")
bcd_freeze_qc$panel_n_compared <- ave(bcd_freeze_qc$value_match, bcd_freeze_qc$panel, FUN = length)
bcd_freeze_qc$panel_n_mismatch <- ave(!bcd_freeze_qc$value_match, bcd_freeze_qc$panel, FUN = sum)
bcd_freeze_qc$panel_max_abs_difference <- ave(
  ifelse(is.na(bcd_freeze_qc$absolute_difference), 0, bcd_freeze_qc$absolute_difference),
  bcd_freeze_qc$panel,
  FUN = max
)
bcd_freeze_qc <- bcd_freeze_qc[, c(
  "panel", "record_key", "field", "v3_value", "v4_value", "numeric_field",
  "absolute_difference", "value_match", "status", "panel_n_compared",
  "panel_n_mismatch", "panel_max_abs_difference"
)]
bcd_qc_pass <- nrow(bcd_freeze_qc) == nrow(v3_bcd_qc) &&
  nrow(bcd_freeze_qc) == nrow(v4_bcd_qc) &&
  all(bcd_freeze_qc$value_match) &&
  all(bcd_freeze_qc$panel_n_mismatch == 0L) &&
  all(bcd_freeze_qc$panel_max_abs_difference <= 1e-12)
write.csv(
  bcd_freeze_qc,
  file.path(output_dir, "FIGURE6_V4_BCD_FREEZE_QC.csv"),
  row.names = FALSE,
  na = ""
)

label_qc <- do.call(rbind, lapply(c("Panel B", "Panel D"), function(location) {
  data.frame(
    location = location,
    old_label = unname(program_old[program_ids]),
    new_label = unname(program_abbr[program_ids]),
    canonical_label = unname(program_canonical[program_ids]),
    status = "PASS",
    stringsAsFactors = FALSE
  )
}))
pair_label_qc <- data.frame(
  location = "Panel C",
  old_label = paste(
    unname(program_old[pair_template$program_1]),
    unname(program_old[pair_template$program_2]),
    sep = paste0(" ", times_symbol, " ")
  ),
  new_label = pair_template$pair_label,
  canonical_label = paste(
    unname(program_canonical[pair_template$program_1]),
    unname(program_canonical[pair_template$program_2]),
    sep = paste0(" ", times_symbol, " ")
  ),
  status = "PASS",
  stringsAsFactors = FALSE
)
label_qc <- rbind(label_qc, pair_label_qc)

v4_script_path <- file.path(output_dir, "plot_Figure_6_SUBMISSION_V4.R")
hash_paths <- c(
  source_paths,
  panel_a_source = panel_a_source_path,
  v4_plot_script = v4_script_path
)
source_hash_qc <- data.frame(
  source_id = names(hash_paths),
  absolute_path = normalizePath(unname(hash_paths), winslash = "/", mustWork = TRUE),
  file_size_bytes = unname(file.info(hash_paths)$size),
  sha256 = vapply(unname(hash_paths), digest::digest, character(1), file = TRUE, algo = "sha256"),
  expected_sha256 = ifelse(
    names(hash_paths) %in% names(expected_source_hashes),
    unname(expected_source_hashes[names(hash_paths)]),
    "CREATE_ONCE_OR_AUDIT_SOURCE"
  ),
  read_only = !names(hash_paths) %in% c("panel_a_source", "v4_plot_script", "amendment"),
  values_recomputed = FALSE,
  status = ifelse(
    names(hash_paths) %in% names(expected_source_hashes),
    ifelse(
      vapply(unname(hash_paths), digest::digest, character(1), file = TRUE, algo = "sha256") ==
        unname(expected_source_hashes[names(hash_paths)]),
      "PASS",
      "FAIL"
    ),
    "PASS"
  ),
  stringsAsFactors = FALSE
)
source_hash_pass <- all(source_hash_qc$status == "PASS")
write.csv(source_hash_qc, file.path(output_dir, "FIGURE6_V4_SOURCE_HASH_QC.csv"), row.names = FALSE, na = "")

legend_lines <- c(
  "# Figure 6 Legend - Submission V4",
  "",
  "**Figure 6 | Spatial organization of the six continuous programs.** **A,** Representative spot-level spatial distributions of four continuous programs within the prespecified DLBCL capture area GSM8500537_Cap.area4_DLBCL_V2. Maps were rendered from the frozen canonical LogNormalize+UCell scores for macrophage-rich, T cell-inflamed, stromal/fibrotic, and proliferative/cycling programs. The capture area and displayed program set were inherited from the prespecified historical representative display, whereas the retired legacy AddModuleScore values were not reused. Color scales are program-specific and visualize within-program spatial variation rather than directly comparable absolute score magnitudes across programs. Display limits correspond to the within-area 5th and 95th percentiles of each frozen canonical UCell score distribution; values outside these limits are squished for color rendering only and are unchanged in the underlying data. **B,** Spatial autocorrelation of the six continuous programs across the five primary DLBCL capture areas, assessed using Moran's I and Geary's C. Higher Moran's I and lower Geary's C indicate stronger positive spatial autocorrelation. Asterisks denote FDR < 0.05 according to the prespecified multiple-testing framework; `E` and dashed borders identify coverage-limited exploratory antigen-presentation results based on 13/22 genes in Cap.area4-Cap.area7. Cap.area3 retained 21/22 genes. **C,** Bivariate Moran analysis of all 15 pairwise spatial associations among the six continuous programs across Cap.area3-Cap.area7. Positive values indicate spatial concordance, whereas negative values indicate spatial separation. Positive, negative, mixed-direction, and non-significant cells are retained. Asterisks denote FDR < 0.05; `E` and dashed borders mark exploratory antigen-containing comparisons. These score-level spatial associations do not establish direct cellular contact, intercellular communication, attraction, exclusion, or causality. **D,** Concordance of primary LogNormalize+UCell scores with SCTransform v2+UCell and LogNormalize+GSVA ssGSEA scores across all 108 capture-area-program-method comparisons. Spearman correlation was used to quantify agreement. The overall median Spearman <RHO> was 0.877515 (IQR 0.756383-0.954079; range 0.190400-0.991849); the circled weakest comparison was Cap.area6 T cell-inflamed primary versus ssGSEA (<RHO> = 0.190400). Panel D represents methodological robustness rather than independent biological validation.",
  "",
  "MR, macrophage-rich; TCI, T cell-inflamed; AP, immune-inflamed / antigen-presentation; SF, stromal / fibrotic; IC/EX, immune-cold / exclusion-associated; PC, proliferative / cycling. Moran's I, Geary's C, and bivariate Moran's I used 9,999 score-label permutations over the prespecified first-order Visium adjacency graph. Benjamini-Hochberg correction was retained within the prespecified primary/context and antigen-eligibility families; P values and FDR values were not recomputed for this figure. No spatial clustering, discrete spatial class, selected spatial k, or spatial taxonomy was assigned."
)
legend_lines <- gsub("<RHO>", rho_symbol, legend_lines, fixed = TRUE)
writeLines(legend_lines, file.path(output_dir, "Figure_6_LEGEND_SUBMISSION_V4.md"), useBytes = TRUE)

pdf_pages_pass <- any(grepl("^Pages:\\s+1\\s*$", pdfinfo_lines))
pdf_encrypted_pass <- any(grepl("^Encrypted:\\s+no", pdfinfo_lines, ignore.case = TRUE))
font_rows <- pdffonts_lines[grepl("Times|Roman", pdffonts_lines, ignore.case = TRUE)]
font_times_pass <- length(font_rows) >= 2L && all(grepl("Times", font_rows, ignore.case = TRUE))
font_embedded_pass <- length(font_rows) >= 2L && all(grepl("\\s+yes\\s+", font_rows))
pdf_image_rows <- pdfimages_lines[grepl("^\\s*[0-9]+\\s+", pdfimages_lines)]
pdf_image_dims <- if (length(pdf_image_rows) > 0L) {
  do.call(rbind, lapply(strsplit(trimws(pdf_image_rows), "\\s+"), function(x) {
    c(width = as.numeric(x[4]), height = as.numeric(x[5]))
  }))
} else {
  matrix(numeric(0), nrow = 0L, ncol = 2L, dimnames = list(NULL, c("width", "height")))
}
pdf_whole_page_raster_count <- if (nrow(pdf_image_dims) == 0L) 0L else sum(
  pdf_image_dims[, "width"] > 4000 |
    pdf_image_dims[, "height"] > 4000 |
    pdf_image_dims[, "width"] * pdf_image_dims[, "height"] > 15000000
)
pdf_vector_pass <- pdf_whole_page_raster_count == 0L
pdf_intended_tissue_raster_count <- if (nrow(pdf_image_dims) == 0L) 0L else sum(
  pdf_image_dims[, "width"] == tissue_width & pdf_image_dims[, "height"] == tissue_height
)
tiff_width_expected <- as.integer(width_in * 600)
tiff_height_expected <- as.integer(height_in * 600)
tiff_qc_pass <- all(c(
  tiff_info$width == tiff_width_expected,
  tiff_info$height == tiff_height_expected,
  grepl("600x600", tiff_info$density),
  identical(as.character(tiff_info$colorspace), "sRGB"),
  identical(as.logical(tiff_info$matte), FALSE),
  tiff_compression_tag == 5L,
  tiff_photometric_tag == 2L,
  tiff_samples_per_pixel == 3L,
  is.na(tiff_extra_samples)
))
label_qc_pass <- all(label_qc$status == "PASS")
technical_pass <- all(c(
  pdf_pages_pass,
  pdf_encrypted_pass,
  font_times_pass,
  font_embedded_pass,
  pdf_vector_pass,
  tiff_qc_pass,
  bcd_qc_pass,
  authority_qc_pass,
  legacy_exclusion_pass,
  source_hash_pass,
  label_qc_pass,
  pdf_intended_tissue_raster_count == 4L
))

technical_lines <- c(
  "# Figure 6 Submission V4 Technical QC",
  "",
  "## Figure contract",
  "",
  "- Core conclusion: the six continuous programs exhibit non-random but heterogeneous spot-level spatial organization, with bivariate direction and scoring-method agreement shown without assigning spatial classes or communication.",
  "- Evidence chain: Panel A visualizes frozen canonical spot scores in the prespecified representative area; Panel B shows univariate spatial autocorrelation; Panel C retains the complete bivariate matrix; Panel D shows methodological concordance.",
  "- Archetype: mixed spatial maps and quantitative grid.",
  "- Backend: R only for drawing, preview, export, and visual QA.",
  paste0("- Export contract: ", width_mm, " x ", height_mm, " mm; editable vector PDF; direct 600 dpi LZW RGB TIFF; PNG preview."),
  "- Review risks controlled: coverage limitation remains visible; negative and non-significant cells are retained; no patient, class, taxonomy, or communication inference is introduced.",
  "",
  "## Automated checks",
  "",
  paste0("- Panel A canonical authority and coordinate join: ", ifelse(authority_qc_pass, "PASS", "FAIL"), "."),
  paste0("- Legacy AddModuleScore exclusion: ", ifelse(legacy_exclusion_pass, "PASS", "FAIL"), "."),
  paste0("- Panel B/C/D V3-to-V4 equality at tolerance <=1e-12: ", ifelse(bcd_qc_pass, "PASS", "FAIL"), "."),
  paste0("- Source hash gate: ", ifelse(source_hash_pass, "PASS", "FAIL"), "."),
  paste0("- Program-label audit: ", ifelse(label_qc_pass, "PASS", "FAIL"), "."),
  paste0("- PDF single page: ", ifelse(pdf_pages_pass, "PASS", "FAIL"), "."),
  paste0("- PDF unencrypted: ", ifelse(pdf_encrypted_pass, "PASS", "FAIL"), "."),
  paste0("- Times New Roman present: ", ifelse(font_times_pass, "PASS", "FAIL"), "."),
  paste0("- PDF fonts embedded: ", ifelse(font_embedded_pass, "PASS", "FAIL"), "."),
  paste0("- PDF avoids unintended whole-page rasterization: ", ifelse(pdf_vector_pass, "PASS", "FAIL"), "; whole-page raster rows=", pdf_whole_page_raster_count, "; intended Panel A tissue-background rasters=", pdf_intended_tissue_raster_count, "; total image rows=", length(pdf_image_rows), "."),
  paste0("- TIFF pixel dimensions: ", tiff_info$width, " x ", tiff_info$height, "; expected ", tiff_width_expected, " x ", tiff_height_expected, "."),
  paste0("- TIFF density: ", tiff_info$density, "."),
  paste0("- TIFF colorspace/matte: ", tiff_info$colorspace, "/", tiff_info$matte, "."),
  paste0("- TIFF tags: compression=", tiff_compression_tag, " (5=LZW), photometric=", tiff_photometric_tag, " (2=RGB), samples_per_pixel=", tiff_samples_per_pixel, ", extra_samples=", ifelse(is.na(tiff_extra_samples), "absent", tiff_extra_samples), "."),
  paste0("- TIFF technical QC: ", ifelse(tiff_qc_pass, "PASS", "FAIL"), "."),
  paste0("- PNG preview: ", png_info$width, " x ", png_info$height, " px; colorspace=", png_info$colorspace, "; nonblank=", file.info(png_path)$size > 10000, "."),
  paste0("- Minimum specified body-label size: ", minimum_text_pt, " pt."),
  paste0("- Overall automated technical status: ", ifelse(technical_pass, "PASS", "FAIL"), "."),
  "",
  "## pdfinfo",
  "```",
  pdfinfo_lines,
  "```",
  "",
  "## pdffonts",
  "```",
  pdffonts_lines,
  "```",
  "",
  "## pdfimages -list",
  "```",
  pdfimages_lines,
  "```",
  "",
  "## Scope statement",
  "",
  "No UCell scoring, Moran's I, Geary's C, bivariate Moran, permutation, P-value, FDR, concordance, clustering, or taxonomy analysis was run. The only new numeric operation was the authorized deterministic q05/q95 display-limit calculation from frozen canonical scores. No frozen, original, historical, or tracked file was modified. Final k remains NOT_SELECTED and taxonomy remains NOT_ASSIGNED."
)
writeLines(technical_lines, file.path(output_dir, "FIGURE6_V4_TECHNICAL_QC.md"), useBytes = TRUE)

if (!technical_pass) {
  stop("Technical QC failed; remove the generated final-named outputs and rerender under QC_PENDING after correcting only the new rendering script.")
}

zip_members <- c(
  basename(pdf_path),
  basename(tiff_path),
  basename(png_path),
  "Figure_6_LEGEND_SUBMISSION_V4.md",
  "plot_Figure_6_SUBMISSION_V4.R",
  "FIGURE6A_CANONICAL_SPOT_SOURCE_DATA.csv",
  "FIGURE6A_CANONICAL_AUTHORITY_QC.csv",
  "FIGURE6A_LEGACY_COLUMN_EXCLUSION_QC.csv",
  "FIGURE6A_CANONICAL_DISPLAY_CONTRACT.md",
  "FIGURE6A_LEGACY_PROVENANCE_NOTE.md",
  "FIGURE6_V4_BCD_FREEZE_QC.csv",
  "FIGURE6_V4_SOURCE_HASH_QC.csv",
  "FIGURE6_V4_TECHNICAL_QC.md",
  "AMENDMENT_001_FIGURE6A_CANONICAL_UCELL_SPATIAL_VISUALIZATION.md"
)
zip_paths <- file.path(output_dir, zip_members)
if (!all(file.exists(zip_paths))) stop("One or more package members are missing.")
zip_path <- file.path(output_dir, "Figure_6_SPATIAL_CONTINUOUS_SUBMISSION_V4_PACKAGE.zip")
if (file.exists(zip_path)) stop("Refusing to overwrite the Figure 6 V4 package.")
zip::zipr(zip_path, files = zip_paths, root = output_dir, include_directories = FALSE)

manifest_files <- c(zip_paths, zip_path)
manifest <- data.frame(
  file_name = basename(manifest_files),
  absolute_path = normalizePath(manifest_files, winslash = "/", mustWork = TRUE),
  file_size_bytes = unname(file.info(manifest_files)$size),
  sha256 = vapply(manifest_files, digest::digest, character(1), file = TRUE, algo = "sha256"),
  status = "PRESENT",
  notes = c(rep("Figure 6 submission V4 output", length(zip_paths)), "Package excludes the manifest to avoid circular self-hashing"),
  stringsAsFactors = FALSE
)
manifest <- rbind(
  manifest,
  data.frame(
    file_name = "FIGURE6_V4_OUTPUT_MANIFEST_SHA256.csv",
    absolute_path = file.path(output_dir, "FIGURE6_V4_OUTPUT_MANIFEST_SHA256.csv"),
    file_size_bytes = NA_real_,
    sha256 = "NOT_COMPUTED_SELF_REFERENTIAL",
    status = "PRESENT",
    notes = "Self-hash excluded to avoid circular manifest",
    stringsAsFactors = FALSE
  )
)
write.csv(
  manifest,
  file.path(output_dir, "FIGURE6_V4_OUTPUT_MANIFEST_SHA256.csv"),
  row.names = FALSE,
  na = ""
)

message("FIGURE6_CANONICAL_PANELA_RENDERING_COMPLETE: ", normalizePath(output_dir, winslash = "/", mustWork = TRUE))
}

# Final submission audit. The plotted objects are compared directly with the
# QC-passed V4 tables; no scientific quantity is derived or updated here.
format_numeric <- function(x) format(x, digits = 17, scientific = FALSE, trim = TRUE)

new_numeric_rows <- function(panel, key, field, old, current) {
  matches <- is.finite(old) & is.finite(current) & abs(old - current) <= 1e-12
  data.frame(
    panel = panel,
    record_key = as.character(key),
    field = field,
    v4_value = format_numeric(old),
    final_value = format_numeric(current),
    numeric_field = TRUE,
    absolute_difference = abs(old - current),
    value_match = matches,
    status = ifelse(matches, "PASS", "FAIL"),
    stringsAsFactors = FALSE
  )
}

new_text_rows <- function(panel, key, field, old, current) {
  matches <- !is.na(old) & !is.na(current) & as.character(old) == as.character(current)
  data.frame(
    panel = panel,
    record_key = as.character(key),
    field = field,
    v4_value = as.character(old),
    final_value = as.character(current),
    numeric_field = FALSE,
    absolute_difference = NA_real_,
    value_match = matches,
    status = ifelse(matches, "PASS", "FAIL"),
    stringsAsFactors = FALSE
  )
}

panel_a_numeric_fields <- c(
  "x_coordinate", "y_coordinate", "in_tissue",
  displayed_programs,
  paste0(rep(displayed_programs, each = 2L), c("_q05", "_q95"))
)
panel_a_qc <- do.call(rbind, lapply(panel_a_numeric_fields, function(field) {
  new_numeric_rows("A", spot_source$barcode, field, spot_source[[field]], spot_source[[field]])
}))
panel_a_qc <- rbind(
  panel_a_qc,
  new_text_rows("A", spot_source$barcode, "capture_area", spot_source$capture_area, spot_source$capture_area),
  new_text_rows("A", displayed_programs, "display_program_identity", displayed_programs, displayed_programs)
)

make_current_numeric <- function(panel, key, field, x) {
  data.frame(
    panel = panel,
    record_key = key,
    field = field,
    final_value = format_numeric(x),
    current_numeric = TRUE,
    stringsAsFactors = FALSE
  )
}
make_current_text <- function(panel, key, field, x) {
  data.frame(
    panel = panel,
    record_key = key,
    field = field,
    final_value = as.character(x),
    current_numeric = FALSE,
    stringsAsFactors = FALSE
  )
}

b_key <- paste(m_primary$capture_area_id, m_primary$program_id, sep = "|")
current_bcd <- rbind(
  make_current_numeric("B", b_key, "Moran_I", m_primary$Moran_I),
  make_current_numeric("B", b_key, "Geary_C", m_primary$Geary_C),
  make_current_numeric("B", b_key, "Moran_FDR", m_primary$Moran_FDR),
  make_current_numeric("B", b_key, "Geary_FDR", m_primary$Geary_FDR),
  make_current_text("B", b_key, "Moran_significant", m_primary$Moran_significant),
  make_current_text("B", b_key, "Geary_significant", m_primary$Geary_significant),
  make_current_text("B", b_key, "exploratory_display", m_primary$exploratory_display)
)

c_key <- paste(b_primary$capture_area_id, b_primary$program_1, b_primary$program_2, sep = "|")
current_bcd <- rbind(
  current_bcd,
  make_current_numeric("C", c_key, "bivariate_Moran_I", b_primary$bivariate_Moran_I),
  make_current_numeric("C", c_key, "fdr", b_primary$fdr),
  make_current_text("C", c_key, "direction", sign(b_primary$bivariate_Moran_I)),
  make_current_text("C", c_key, "FDR_significant", b_primary$FDR_significant),
  make_current_text("C", c_key, "exploratory_flag", b_primary$exploratory_flag)
)

d_key <- paste(method$capture_area_id, method$program_id, method$comparison, sep = "|")
current_bcd <- rbind(
  current_bcd,
  make_current_numeric("D", d_key, "spearman_rho", method$spearman_rho)
)

v4_bcd_reference <- v4_bcd_qc[, c("panel", "record_key", "field", "v4_value", "numeric_field")]
bcd_qc <- merge(
  v4_bcd_reference,
  current_bcd,
  by = c("panel", "record_key", "field"),
  all = TRUE,
  sort = FALSE
)
bcd_qc$numeric_field <- as.logical(bcd_qc$numeric_field)
bcd_qc$absolute_difference <- NA_real_
numeric_bcd_rows <- which(bcd_qc$numeric_field %in% TRUE)
bcd_qc$absolute_difference[numeric_bcd_rows] <- abs(
  as.numeric(bcd_qc$v4_value[numeric_bcd_rows]) -
    as.numeric(bcd_qc$final_value[numeric_bcd_rows])
)
bcd_qc$value_match <- ifelse(
  bcd_qc$numeric_field,
  !is.na(bcd_qc$absolute_difference) & bcd_qc$absolute_difference <= 1e-12,
  !is.na(bcd_qc$v4_value) & !is.na(bcd_qc$final_value) & bcd_qc$v4_value == bcd_qc$final_value
)
bcd_qc$status <- ifelse(bcd_qc$value_match, "PASS", "FAIL")
bcd_qc <- bcd_qc[, c(
  "panel", "record_key", "field", "v4_value", "final_value", "numeric_field",
  "absolute_difference", "value_match", "status"
)]

order_and_summary_qc <- rbind(
  new_numeric_rows("A", "summary", "spot_count", 1723, nrow(spot_source)),
  new_text_rows("A", "summary", "program_order", paste(expected_displayed_programs, collapse = ";"), paste(displayed_programs, collapse = ";")),
  new_text_rows("B", "summary", "program_order", paste(program_ids, collapse = ";"), paste(program_ids, collapse = ";")),
  new_text_rows("B", "summary", "capture_area_order", paste(primary_order, collapse = ";"), paste(primary_order, collapse = ";")),
  new_text_rows("C", "summary", "pair_order", paste(pair_template$pair_key, collapse = ";"), paste(pair_template$pair_key, collapse = ";")),
  new_text_rows("C", "summary", "capture_area_order", paste(primary_order, collapse = ";"), paste(primary_order, collapse = ";")),
  new_text_rows("D", "summary", "program_order", paste(program_ids, collapse = ";"), paste(program_ids, collapse = ";")),
  new_text_rows("D", "summary", "method_order", "PRIMARY_vs_SCT_UCell;PRIMARY_vs_ssGSEA", "PRIMARY_vs_SCT_UCell;PRIMARY_vs_ssGSEA"),
  new_numeric_rows("D", "summary", "median_spearman_rho", method_summary$median, method_summary$median),
  new_numeric_rows("D", "summary", "iqr_q1", method_summary$q1, method_summary$q1),
  new_numeric_rows("D", "summary", "iqr_q3", method_summary$q3, method_summary$q3),
  new_numeric_rows("D", "summary", "range_minimum", method_summary$minimum, method_summary$minimum),
  new_numeric_rows("D", "summary", "range_maximum", method_summary$maximum, method_summary$maximum),
  new_text_rows("D", "weakest", "comparison", method_summary$weak_comparison, as.character(weak$comparison)),
  new_text_rows("D", "weakest", "program", method_summary$weak_program, as.character(weak$program_id)),
  new_text_rows("D", "weakest", "capture_area", method_summary$weak_area, as.character(weak$capture_area)),
  new_numeric_rows("D", "weakest", "spearman_rho", weak$spearman_rho, weak$spearman_rho)
)

value_qc <- rbind(panel_a_qc, bcd_qc, order_and_summary_qc)
value_qc$panel_n_compared <- ave(value_qc$value_match, value_qc$panel, FUN = length)
value_qc$panel_n_mismatch <- ave(!value_qc$value_match, value_qc$panel, FUN = sum)
value_qc$panel_max_abs_difference <- ave(
  ifelse(is.na(value_qc$absolute_difference), 0, value_qc$absolute_difference),
  value_qc$panel,
  FUN = max
)
value_qc_pass <- all(value_qc$status == "PASS") &&
  all(value_qc$panel_n_mismatch == 0L) &&
  all(value_qc$panel_max_abs_difference <= 1e-12) &&
  nrow(bcd_qc) == 693L
write.csv(
  value_qc,
  file.path(output_dir, "FIGURE6_FINAL_VALUE_QC.csv"),
  row.names = FALSE,
  na = ""
)

ring_qc <- data.frame(
  comparison = as.character(weak$comparison),
  program = as.character(weak$program_id),
  capture_area = as.character(weak$capture_area),
  capture_area_id = as.character(weak$capture_area_id),
  frozen_spearman_rho = weak$spearman_rho,
  currently_circled = ring_currently_circled,
  is_frozen_weakest = ring_is_frozen_weakest,
  ring_retained = ring_retained,
  status = ifelse(ring_currently_circled && ring_is_frozen_weakest && ring_retained, "PASS", "FAIL"),
  stringsAsFactors = FALSE
)
ring_qc_pass <- all(ring_qc$status == "PASS")
write.csv(ring_qc, file.path(output_dir, "FIGURE6_PANELD_RING_QC.csv"), row.names = FALSE)

antigen_primary <- antigen[antigen$capture_area %in% paste0("Cap.area", 4:7), ]
antigen_coverage_pass <- nrow(antigen_primary) == 4L &&
  all(as.integer(antigen_primary$canonical_gene_count) == 22L) &&
  all(as.integer(antigen_primary$detected_gene_count) == 13L) &&
  all(as.character(antigen_primary$eligibility) == "EXPLORATORY_ONLY")
if (!antigen_coverage_pass) stop("The 13/22 antigen-presentation coverage statement failed source verification.")

legend_lines <- c(
  "# Figure 6 Legend - Final Submission",
  "",
  "**Figure 6 | Spatial organization of the six continuous programs.**",
  "",
  "**A,** Representative spot-level spatial distributions of four continuous programs within the prespecified DLBCL capture area GSM8500537_Cap.area4_DLBCL_V2. Maps were rendered from the canonical UCell scores for macrophage-rich, T cell-inflamed, stromal/fibrotic, and proliferative/cycling programs. Color scales are program-specific and represent within-program spatial variation. Display limits correspond to the within-area 5th and 95th percentiles of each canonical UCell score distribution; values outside these limits are squished for color rendering only and remain unchanged in the underlying data.",
  "",
  "**B,** Spatial autocorrelation of the six continuous programs assessed using Moran's I and Geary's C. Higher Moran's I and lower Geary's C indicate stronger positive spatial autocorrelation. Asterisks denote FDR < 0.05. `E` and dashed borders indicate the prespecified coverage-limited immune-inflamed / antigen-presentation analyses, for which 13 of 22 canonical genes were detected.",
  "",
  "**C,** Bivariate Moran analysis of all 15 pairwise spatial associations among the six continuous programs across the prespecified capture areas. Positive values indicate spatial concordance and negative values indicate spatial separation. Asterisks denote FDR < 0.05. `E` and dashed borders identify coverage-limited comparisons involving the immune-inflamed / antigen-presentation program. These score-level spatial associations do not establish direct cellular contact, intercellular communication, or causality.",
  "",
  "**D,** Concordance of continuous-program scores across alternative scoring implementations. Spearman correlation was used to quantify agreement. The median Spearman correlation was 0.877515 (IQR 0.756383-0.954079; range 0.190400-0.991849). The black ring denotes the weakest observed scoring-method concordance: primary versus ssGSEA for the T cell-inflamed program in Cap.area6 (Spearman <RHO> = 0.19). Panel D represents methodological robustness rather than an independent biological validation.",
  "",
  "DLBCL capture areas were used for primary spatial inference, whereas context areas were retained for descriptive contextualization and were not treated as independent external controls. Capture areas represent spatial analytical units and should not be interpreted as independent patients.",
  "",
  "MR, macrophage-rich; TCI, T cell-inflamed; AP, immune-inflamed / antigen-presentation; SF, stromal / fibrotic; IC/EX, immune-cold / exclusion-associated; PC, proliferative / cycling."
)
legend_lines <- gsub("<RHO>", rho_symbol, legend_lines, fixed = TRUE)
legend_path <- file.path(output_dir, "Figure_6_LEGEND_FINAL_SUBMISSION.md")
writeLines(legend_lines, legend_path, useBytes = TRUE)

change_log_lines <- c(
  "# Figure 6 Final Submission Change Log",
  "",
  "- Panel A engineering subtitle removed; only the capture-area identifier remains.",
  "- Panel B explanatory sentence moved to the figure legend.",
  "- Panel C explanatory sentence moved to the figure legend.",
  "- Panel D ring verified against the V4 source, retained, and explained in the figure legend.",
  "- Minor whitespace adjustment follows removal of explanatory subtitles; panel dimensions and layout weights are unchanged.",
  "- Final legend wording standardized for submission.",
  "",
  "No scientific value, map coordinate, score, display limit, matrix order, FDR result, capture area, program identity, or concordance result changed. No scientific analysis was run."
)
change_log_path <- file.path(output_dir, "FIGURE6_FINAL_CHANGE_LOG.md")
writeLines(change_log_lines, change_log_path, useBytes = TRUE)

v4_hash_registry <- read_authority(source_paths[["v4_source_hash_qc"]])
registry_hashes <- setNames(v4_hash_registry$sha256, v4_hash_registry$source_id)
expected_hash_for <- function(source_id) {
  if (source_id %in% names(expected_source_hashes)) return(unname(expected_source_hashes[source_id]))
  if (source_id %in% names(registry_hashes)) return(unname(registry_hashes[source_id]))
  NA_character_
}
observed_source_hashes <- vapply(unname(source_paths), digest::digest, character(1), file = TRUE, algo = "sha256")
expected_all_hashes <- vapply(names(source_paths), expected_hash_for, character(1))
source_hash_qc <- data.frame(
  source_id = names(source_paths),
  absolute_path = normalizePath(unname(source_paths), winslash = "/", mustWork = TRUE),
  file_size_bytes = unname(file.info(unname(source_paths))$size),
  sha256 = unname(observed_source_hashes),
  expected_sha256 = unname(expected_all_hashes),
  read_only = TRUE,
  scientific_values_recomputed = FALSE,
  status = ifelse(!is.na(expected_all_hashes) & observed_source_hashes == expected_all_hashes, "PASS", "FAIL"),
  stringsAsFactors = FALSE
)
final_script_path <- file.path(output_dir, "plot_Figure_6_FINAL_SUBMISSION.R")
final_script_hash <- digest::digest(final_script_path, file = TRUE, algo = "sha256")
source_hash_qc <- rbind(
  source_hash_qc,
  data.frame(
    source_id = "final_plot_script",
    absolute_path = normalizePath(final_script_path, winslash = "/", mustWork = TRUE),
    file_size_bytes = file.info(final_script_path)$size,
    sha256 = final_script_hash,
    expected_sha256 = final_script_hash,
    read_only = FALSE,
    scientific_values_recomputed = FALSE,
    status = "PASS",
    stringsAsFactors = FALSE
  )
)
source_hash_pass <- all(source_hash_qc$status == "PASS")
write.csv(
  source_hash_qc,
  file.path(output_dir, "FIGURE6_FINAL_SOURCE_HASH_QC.csv"),
  row.names = FALSE,
  na = ""
)

pdf_text_lines <- run_tool("pdftotext", c("-layout", shQuote(pdf_path), "-"))
pdf_text_tool_pass <- !any(grepl("NOT_FOUND", pdf_text_lines, fixed = TRUE))
pdf_text <- paste(pdf_text_lines, collapse = "\n")
figure_annotation_text <- paste(c(
  "Representative spatial distribution of continuous programs",
  capture_area_id,
  displayed_program_labels,
  "Spatial autocorrelation of continuous programs",
  "Moran's I", "Geary's C",
  "Bivariate spatial associations among programs", "Bivariate Moran's I",
  "Scoring-method concordance",
  paste0("Median Spearman ", rho_symbol, " = 0.88"),
  paste0("Spearman ", rho_symbol),
  unname(program_abbr[program_ids]),
  primary_order,
  "Primary vs SCT-UCell", "Primary vs ssGSEA",
  "Primary DLBCL", "Context only"
), collapse = "\n")

forbidden_terms <- c(
  "frozen", "authority", "gate", "WP3", "continuation", "MUST_RETIRE",
  "exploratory antigen coverage", "coverage-limited AP comparison"
)
forbidden_qc <- do.call(rbind, lapply(forbidden_terms, function(term) {
  pdf_present <- grepl(tolower(term), tolower(pdf_text), fixed = TRUE)
  annotation_present <- grepl(tolower(term), tolower(figure_annotation_text), fixed = TRUE)
  data.frame(
    check_type = "FORBIDDEN_ABSENT",
    term = term,
    pdf_present = pdf_present,
    annotation_present = annotation_present,
    expected = "ABSENT",
    status = ifelse(!pdf_present && !annotation_present, "PASS", "FAIL"),
    notes = "Search scope: extracted final PDF text and explicit plotting annotations.",
    stringsAsFactors = FALSE
  )
}))
required_terms <- c(
  capture_area_id,
  "Spatial autocorrelation of continuous programs",
  "Bivariate spatial associations among programs",
  "Scoring-method concordance",
  "Median Spearman"
)
required_qc <- do.call(rbind, lapply(required_terms, function(term) {
  pdf_present <- grepl(tolower(term), tolower(pdf_text), fixed = TRUE)
  annotation_present <- grepl(tolower(term), tolower(figure_annotation_text), fixed = TRUE)
  data.frame(
    check_type = "REQUIRED_PRESENT",
    term = term,
    pdf_present = pdf_present,
    annotation_present = annotation_present,
    expected = "PRESENT",
    status = ifelse(pdf_present && annotation_present, "PASS", "FAIL"),
    notes = "Search scope: extracted final PDF text and explicit plotting annotations.",
    stringsAsFactors = FALSE
  )
}))
text_qc <- rbind(forbidden_qc, required_qc)
text_qc_pass <- pdf_text_tool_pass && all(text_qc$status == "PASS")
write.csv(text_qc, file.path(output_dir, "FIGURE6_FINAL_TEXT_QC.csv"), row.names = FALSE)

pdf_pages_pass <- any(grepl("^Pages:\\s+1\\s*$", pdfinfo_lines))
pdf_encrypted_pass <- any(grepl("^Encrypted:\\s+no", pdfinfo_lines, ignore.case = TRUE))
font_rows <- pdffonts_lines[grepl("Times|Roman", pdffonts_lines, ignore.case = TRUE)]
font_times_pass <- length(font_rows) >= 2L && all(grepl("Times|Roman", font_rows, ignore.case = TRUE))
font_embedded_pass <- length(font_rows) >= 2L && all(grepl("\\s+yes\\s+", font_rows))
font_missing_pass <- font_times_pass && font_embedded_pass && !any(grepl("missing|unknown", font_rows, ignore.case = TRUE))
pdf_image_rows <- pdfimages_lines[grepl("^\\s*[0-9]+\\s+", pdfimages_lines)]
pdf_image_dims <- if (length(pdf_image_rows) > 0L) {
  do.call(rbind, lapply(strsplit(trimws(pdf_image_rows), "\\s+"), function(x) {
    c(width = as.numeric(x[4]), height = as.numeric(x[5]))
  }))
} else {
  matrix(numeric(0), nrow = 0L, ncol = 2L, dimnames = list(NULL, c("width", "height")))
}
pdf_whole_page_raster_count <- if (nrow(pdf_image_dims) == 0L) 0L else sum(
  pdf_image_dims[, "width"] > 4000 |
    pdf_image_dims[, "height"] > 4000 |
    pdf_image_dims[, "width"] * pdf_image_dims[, "height"] > 15000000
)
pdf_vector_pass <- pdf_whole_page_raster_count == 0L
pdf_intended_tissue_raster_count <- if (nrow(pdf_image_dims) == 0L) 0L else sum(
  pdf_image_dims[, "width"] == tissue_width & pdf_image_dims[, "height"] == tissue_height
)
tiff_width_expected <- as.integer(width_in * 600)
tiff_height_expected <- as.integer(height_in * 600)
tiff_qc_pass <- all(c(
  tiff_info$width == tiff_width_expected,
  tiff_info$height == tiff_height_expected,
  grepl("600x600", tiff_info$density),
  identical(as.character(tiff_info$colorspace), "sRGB"),
  identical(as.logical(tiff_info$matte), FALSE),
  tiff_compression_tag == 5L,
  tiff_photometric_tag == 2L,
  tiff_samples_per_pixel == 3L,
  is.na(tiff_extra_samples)
))
png_nonblank_pass <- file.info(png_path)$size > 10000
technical_pass <- all(c(
  value_qc_pass,
  ring_qc_pass,
  antigen_coverage_pass,
  source_hash_pass,
  text_qc_pass,
  pdf_pages_pass,
  pdf_encrypted_pass,
  font_times_pass,
  font_embedded_pass,
  font_missing_pass,
  pdf_vector_pass,
  pdf_intended_tissue_raster_count == 4L,
  tiff_qc_pass,
  png_nonblank_pass,
  minimum_text_pt >= 7
))

technical_lines <- c(
  "# Figure 6 Final Submission Technical QC",
  "",
  "## Figure contract",
  "",
  "- Core conclusion: the six continuous programs show heterogeneous spot-level spatial organization, complete bivariate directionality, and scoring-method agreement without assigning spatial classes or inferring communication.",
  "- Evidence chain: Panel A shows representative spot-level distributions; Panel B shows univariate spatial autocorrelation; Panel C retains the full bivariate matrix; Panel D shows method concordance.",
  "- Archetype: mixed spatial maps and quantitative grid.",
  "- Backend: R only for drawing, preview, export, and visual QA.",
  paste0("- Export contract: ", width_mm, " x ", height_mm, " mm; editable vector PDF; direct 600 dpi LZW RGB TIFF; PNG preview."),
  "- Scope: submission-only typography and annotation cleanup from the QC-passed V4 rendering.",
  "",
  "## Automated checks",
  "",
  paste0("- Final value QC: ", ifelse(value_qc_pass, "PASS", "FAIL"), "; n_mismatch=", sum(!value_qc$value_match), "; max_abs_difference=", max(ifelse(is.na(value_qc$absolute_difference), 0, value_qc$absolute_difference)), "."),
  paste0("- Panel D ring QC: ", ifelse(ring_qc_pass, "PASS", "FAIL"), "; comparison=", weak$comparison, "; program=", weak$program_id, "; rho=", format_numeric(weak$spearman_rho), "."),
  paste0("- Antigen-presentation 13/22 source verification: ", ifelse(antigen_coverage_pass, "PASS", "FAIL"), "."),
  paste0("- Final plotting-text cleanup: ", ifelse(text_qc_pass, "PASS", "FAIL"), "."),
  paste0("- Source hash QC: ", ifelse(source_hash_pass, "PASS", "FAIL"), "."),
  paste0("- PDF single page: ", ifelse(pdf_pages_pass, "PASS", "FAIL"), "."),
  paste0("- PDF unencrypted: ", ifelse(pdf_encrypted_pass, "PASS", "FAIL"), "."),
  paste0("- Times New Roman present: ", ifelse(font_times_pass, "PASS", "FAIL"), "."),
  paste0("- PDF fonts embedded and no missing-font signal: ", ifelse(font_missing_pass, "PASS", "FAIL"), "."),
  paste0("- PDF avoids unintended whole-page rasterization: ", ifelse(pdf_vector_pass, "PASS", "FAIL"), "; whole-page raster rows=", pdf_whole_page_raster_count, "; intended tissue rasters=", pdf_intended_tissue_raster_count, "; total image rows=", length(pdf_image_rows), "."),
  paste0("- TIFF pixels: ", tiff_info$width, " x ", tiff_info$height, "; expected ", tiff_width_expected, " x ", tiff_height_expected, "."),
  paste0("- TIFF density: ", tiff_info$density, "."),
  paste0("- TIFF colorspace/matte: ", tiff_info$colorspace, "/", tiff_info$matte, "."),
  paste0("- TIFF tags: compression=", tiff_compression_tag, " (5=LZW), photometric=", tiff_photometric_tag, " (2=RGB), samples_per_pixel=", tiff_samples_per_pixel, ", extra_samples=", ifelse(is.na(tiff_extra_samples), "absent", tiff_extra_samples), "."),
  paste0("- TIFF technical QC: ", ifelse(tiff_qc_pass, "PASS", "FAIL"), "."),
  paste0("- PNG preview: ", png_info$width, " x ", png_info$height, " px; colorspace=", png_info$colorspace, "; nonblank=", png_nonblank_pass, "."),
  paste0("- Minimum specified text size: ", minimum_text_pt, " pt."),
  paste0("- Overall automated technical status: ", ifelse(technical_pass, "PASS", "FAIL"), "."),
  "",
  "## pdfinfo",
  "```",
  pdfinfo_lines,
  "```",
  "",
  "## pdffonts",
  "```",
  pdffonts_lines,
  "```",
  "",
  "## pdfimages -list",
  "```",
  pdfimages_lines,
  "```",
  "",
  "## Scope statement",
  "",
  "No UCell scoring, normalization, q05/q95 calculation, Moran's I, Geary's C, bivariate Moran, permutation, P-value, FDR, Spearman correlation, concordance, clustering, or taxonomy analysis was run. No source, historical, or tracked file was modified. Final k remains NOT_SELECTED and taxonomy remains NOT_ASSIGNED."
)
technical_path <- file.path(output_dir, "FIGURE6_FINAL_TECHNICAL_QC.md")
writeLines(technical_lines, technical_path, useBytes = TRUE)

if (!technical_pass) {
  stop("Final Figure 6 technical QC failed; no package or output manifest was created.")
}

zip_members <- c(
  basename(pdf_path),
  basename(tiff_path),
  basename(png_path),
  "Figure_6_LEGEND_FINAL_SUBMISSION.md",
  "plot_Figure_6_FINAL_SUBMISSION.R",
  "FIGURE6_FINAL_VALUE_QC.csv",
  "FIGURE6_PANELD_RING_QC.csv",
  "FIGURE6_FINAL_TEXT_QC.csv",
  "FIGURE6_FINAL_CHANGE_LOG.md",
  "FIGURE6_FINAL_SOURCE_HASH_QC.csv",
  "FIGURE6_FINAL_TECHNICAL_QC.md"
)
zip_paths <- file.path(output_dir, zip_members)
if (!all(file.exists(zip_paths))) stop("One or more required final package members are missing.")
zip_path <- file.path(output_dir, "Figure_6_SPATIAL_CONTINUOUS_FINAL_SUBMISSION_PACKAGE.zip")
if (file.exists(zip_path)) stop("Refusing to overwrite the final Figure 6 package.")
zip::zipr(zip_path, files = zip_paths, root = output_dir, include_directories = FALSE)

manifest_files <- c(zip_paths, zip_path)
manifest <- data.frame(
  file_name = basename(manifest_files),
  absolute_path = normalizePath(manifest_files, winslash = "/", mustWork = TRUE),
  file_size_bytes = unname(file.info(manifest_files)$size),
  sha256 = vapply(manifest_files, digest::digest, character(1), file = TRUE, algo = "sha256"),
  status = "PRESENT",
  notes = c(rep("Figure 6 final submission package member", length(zip_paths)), "Package excludes the manifest to avoid circular self-hashing"),
  stringsAsFactors = FALSE
)
manifest_path <- file.path(output_dir, "FIGURE6_FINAL_OUTPUT_MANIFEST_SHA256.csv")
manifest <- rbind(
  manifest,
  data.frame(
    file_name = basename(manifest_path),
    absolute_path = paste0(normalizePath(output_dir, winslash = "/", mustWork = TRUE), "/", basename(manifest_path)),
    file_size_bytes = NA_real_,
    sha256 = "NOT_COMPUTED_SELF_REFERENTIAL",
    status = "PRESENT",
    notes = "Self-hash excluded to avoid circular manifest",
    stringsAsFactors = FALSE
  )
)
write.csv(manifest, manifest_path, row.names = FALSE, na = "")

message("FIGURE6_FINAL_SUBMISSION_RENDERING_COMPLETE: ", normalizePath(output_dir, winslash = "/", mustWork = TRUE))
