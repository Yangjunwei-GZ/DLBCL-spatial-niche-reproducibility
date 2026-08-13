from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
PROJECT_ROOT = Path(os.environ.get("DLBCL_PROJECT_ROOT", REPOSITORY_ROOT))
ROOT = PROJECT_ROOT / "revision_2026_reviewer_response/07b_suppfig_s4_primary_spatial_ucell_final"
LOCAL_LIB = PROJECT_ROOT / "revision_2026_reviewer_response/06z_suppfig_s2_final_submission/python_lib"
sys.path.insert(0, str(LOCAL_LIB))

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle
from PIL import Image


mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 6.4,
    "axes.titlesize": 7.1,
    "axes.titleweight": "bold",
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "savefig.facecolor": "white",
    "figure.facecolor": "white",
})

AREAS = ["Cap.area3", "Cap.area4", "Cap.area5", "Cap.area6", "Cap.area7"]
PROGRAMS = [
    ("macrophage_rich_UCell", "MR", "Macrophage-rich"),
    ("t_cell_inflamed_UCell", "TCI", "T cell-inflamed"),
    ("antigen_presentation_UCell", "AP", "Antigen-presentation"),
    ("stromal_fibrotic_UCell", "SF", "Stromal/fibrotic"),
    ("immune_cold_exclusion_UCell", "IC/EX", "Immune-cold/exclusion"),
    ("proliferative_cycling_UCell", "PC", "Proliferative/cycling"),
]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


rows = read_csv(ROOT / "S4_PRIMARY_SPATIAL_UCELL_SPOT_LEVEL_SOURCE.csv")
coverage = read_csv(ROOT / "S4_PRIMARY_SPATIAL_PROGRAM_COVERAGE_AUTHORITY.csv")
coverage_status = {(r["Area"], r["Program_ID"] + "_UCell"): r["Coverage_status"] for r in coverage}
area_rows = {area: [r for r in rows if r["Area"] == area] for area in AREAS}

scales = {}
for program_id, _, _ in PROGRAMS:
    values = np.array([float(row[program_id]) for row in rows], dtype=float)
    low, high = np.quantile(values, [0.01, 0.99])
    if not np.isfinite(low) or not np.isfinite(high) or high <= low:
        raise ValueError(f"Invalid display scale for {program_id}: {low}, {high}")
    scales[program_id] = (float(low), float(high))

cmap = LinearSegmentedColormap.from_list(
    "s4_blue",
    ["#F4F7F8", "#D8E7ED", "#A8C7D5", "#6F9FB5", "#376F8A", "#173F5F"],
)

width_mm = 183.0
height_mm = 202.0
fig = plt.figure(figsize=(width_mm / 25.4, height_mm / 25.4), constrained_layout=False)
grid = fig.add_gridspec(
    5,
    6,
    left=0.115,
    right=0.965,
    bottom=0.105,
    top=0.885,
    wspace=0.11,
    hspace=0.16,
)

fig.suptitle(
    "Supplementary Fig. S4 | Canonical UCell spatial distributions\nacross five primary DLBCL capture areas",
    x=0.055,
    y=0.975,
    ha="left",
    va="top",
    fontsize=10.0,
    fontweight="bold",
    color="#20282E",
    linespacing=1.05,
)

axes = {}
for i, area in enumerate(AREAS):
    x = np.array([float(r["y"]) for r in area_rows[area]])
    y = np.array([float(r["x"]) for r in area_rows[area]])
    xpad = max((x.max() - x.min()) * 0.025, 1.0)
    ypad = max((y.max() - y.min()) * 0.025, 1.0)
    for j, (program_id, abbreviation, _) in enumerate(PROGRAMS):
        ax = fig.add_subplot(grid[i, j])
        axes[(area, program_id)] = ax
        values = np.array([float(r[program_id]) for r in area_rows[area]])
        low, high = scales[program_id]
        ax.scatter(
            x,
            y,
            c=values,
            cmap=cmap,
            vmin=low,
            vmax=high,
            s=2.15 if len(values) < 2500 else 1.35,
            marker="o",
            linewidths=0,
            rasterized=False,
        )
        ax.set_xlim(x.min() - xpad, x.max() + xpad)
        ax.set_ylim(y.max() + ypad, y.min() - ypad)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        if i == 0:
            ax.set_title(abbreviation, fontsize=7.5, pad=4.0, color="#20282E")
        if j == 0:
            ax.text(-0.16, 0.5, area, transform=ax.transAxes, rotation=90, ha="center", va="center", fontsize=7.0, fontweight="bold", color="#20282E")
        if coverage_status[(area, program_id)] == "EXPLORATORY_ONLY":
            ax.text(
                0.92,
                0.92,
                "E",
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=6.8,
                fontweight="bold",
                color="#47545C",
                bbox={"boxstyle": "round,pad=0.15", "facecolor": "white", "edgecolor": "#69757D", "linewidth": 0.6},
            )
            ax.add_patch(Rectangle((0.01, 0.01), 0.98, 0.98, transform=ax.transAxes, fill=False, edgecolor="#69757D", linewidth=0.65, linestyle=(0, (2.2, 2.2)), clip_on=False))

color_y = 0.075
segment_width = 0.090
start_x = 0.145
for j, (program_id, abbreviation, _) in enumerate(PROGRAMS):
    cax = fig.add_axes([start_x + j * 0.140, color_y, segment_width, 0.010])
    boundaries = np.linspace(scales[program_id][0], scales[program_id][1], 65, dtype=float)
    gradient = ((boundaries[:-1] + boundaries[1:]) / 2.0).reshape(1, -1)
    cax.pcolormesh(boundaries, [0.0, 1.0], gradient, cmap=cmap, vmin=scales[program_id][0], vmax=scales[program_id][1], shading="flat", edgecolors="face", linewidth=0, rasterized=False)
    cax.set_xlim(*scales[program_id])
    cax.set_ylim(0, 1)
    cax.set_yticks([])
    cax.set_xticks(list(scales[program_id]))
    cax.set_xticklabels([f"{scales[program_id][0]:.2f}", f"{scales[program_id][1]:.2f}"], fontsize=4.8)
    cax.tick_params(axis="x", length=1.5, pad=1.5)
    cax.set_title(abbreviation, fontsize=5.3, pad=1.0, color="#47545C")
    for spine in cax.spines.values():
        spine.set_visible(False)

fig.text(0.055, 0.044, "E, exploratory-only AP display (13/22 detected canonical genes).", ha="left", va="center", fontsize=5.7, color="#47545C")
fig.text(0.055, 0.027, "All maps show frozen canonical UCell scores; 1st-99th percentile display limits are shared across areas within each program.", ha="left", va="center", fontsize=5.7, color="#47545C")

visual = ROOT / "visual_qc"
visual.mkdir(exist_ok=True)
mode = sys.argv[1] if len(sys.argv) > 1 else "draft1"
draft_name = "Supplementary_Fig_S4_DRAFT2.png" if mode == "draft2" else "Supplementary_Fig_S4_DRAFT1.png"
fig.savefig(visual / draft_name, dpi=300, facecolor="white")

if mode == "final":
    pdf = ROOT / "Supplementary_Fig_S4_FINAL_SUBMISSION.pdf"
    svg = ROOT / "Supplementary_Fig_S4_FINAL_SUBMISSION.svg"
    png = ROOT / "Supplementary_Fig_S4_FINAL_SUBMISSION.png"
    tiff = ROOT / "Supplementary_Fig_S4_FINAL_SUBMISSION.tiff"
    fig.savefig(pdf, facecolor="white")
    fig.savefig(svg, facecolor="white")
    fig.savefig(png, dpi=600, facecolor="white")
    image = Image.open(png).convert("RGB")
    image.save(tiff, format="TIFF", compression="tiff_lzw", dpi=(600, 600))
    metadata = {
        "status": "PASS",
        "backend": "Python/matplotlib",
        "matplotlib_version": mpl.__version__,
        "width_mm": width_mm,
        "height_mm": height_mm,
        "png_dimensions": list(image.size),
        "tiff_mode": image.mode,
        "tiff_compression": "LZW",
        "tiff_dpi": [600, 600],
        "display_scaling": "Per-program shared 1st-99th percentile limits across all five areas; source values unchanged.",
        "display_scales": {program_id: {"lower": scales[program_id][0], "upper": scales[program_id][1]} for program_id, _, _ in PROGRAMS},
        "subpanels": 30,
    }
    (ROOT / "technical_qc/S4_FINAL_EXPORT_METADATA.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metadata, indent=2))
else:
    print(json.dumps({"status": mode.upper() + "_CREATED", "path": str(visual / draft_name), "display_scales": scales}, indent=2))

plt.close(fig)
