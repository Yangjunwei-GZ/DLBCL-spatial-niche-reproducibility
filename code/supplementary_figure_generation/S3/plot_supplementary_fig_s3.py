from __future__ import annotations

import csv
import json
import math
import os
import sys
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
PROJECT_ROOT = Path(os.environ.get("DLBCL_PROJECT_ROOT", REPOSITORY_ROOT))
ROOT = PROJECT_ROOT / "revision_2026_reviewer_response/07a_suppfig_s3_single_cell_metadata_final"
LOCAL_LIB = PROJECT_ROOT / "revision_2026_reviewer_response/06z_suppfig_s2_final_submission/python_lib"
sys.path.insert(0, str(LOCAL_LIB))

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import FancyBboxPatch, Rectangle
from PIL import Image


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 7.2,
        "axes.titlesize": 9.0,
        "axes.titleweight": "bold",
        "axes.labelsize": 8.0,
        "axes.linewidth": 0.75,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "xtick.labelsize": 7.0,
        "ytick.labelsize": 7.0,
        "legend.fontsize": 6.7,
        "legend.frameon": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
    }
)

NAVY = "#173F5F"
BLUE = "#2F6B8A"
MID_BLUE = "#7EA9BE"
PALE_BLUE = "#DCEAF1"
LIGHT_BLUE = "#EFF5F8"
DARK = "#20282E"
MID_GRAY = "#69757D"
LIGHT_GRAY = "#E9EDF0"
PALE_GRAY = "#F4F6F7"
GRID = "#D9E0E4"

PATIENTS = ["DLBCL002", "DLBCL007", "DLBCL008", "DLBCL111"]
CELLTYPES = [
    "B cells",
    "T cells CD8",
    "T cells CD4",
    "Tregs",
    "TFH",
    "Monocytes and Macrophages",
    "NK cells",
    "Plasma cells",
    "Others",
]


def read_csv(name: str) -> list[dict[str, str]]:
    with (ROOT / name).open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def add_panel_label(ax, letter: str, title: str, title_size: float = 8.2) -> None:
    ax.text(-0.09, 1.105, letter, transform=ax.transAxes, fontsize=10.5, fontweight="bold", color=DARK, va="bottom")
    ax.text(0.0, 1.105, title, transform=ax.transAxes, fontsize=title_size, fontweight="bold", color=DARK, va="bottom", linespacing=1.05)


cohort = read_csv("GSE182434_S3_COHORT_SUMMARY.csv")
patients = read_csv("GSE182434_S3_PATIENT_COUNTS.csv")
samples = read_csv("GSE182434_S3_PATIENT_SAMPLE_COUNTS.csv")
celltypes = read_csv("GSE182434_S3_CELLTYPE_COUNTS.csv")
patient_celltypes = read_csv("GSE182434_S3_PATIENT_CELLTYPE_COUNTS.csv")

patient_counts = {row["Patient"]: int(row["N_cells"]) for row in patients}
sample_counts = {(row["Patient"], row["Sample"]): int(row["N_cells"]) for row in samples}
celltype_counts = {row["CellType"]: int(row["N_cells"]) for row in celltypes}
celltype_pct = {row["CellType"]: float(row["Percent_of_DLBCL_subset"]) for row in celltypes}
pt_counts = {(row["Patient"], row["CellType"]): int(row["N_cells"]) for row in patient_celltypes}
pt_pct = {(row["Patient"], row["CellType"]): float(row["Within_patient_percent"]) for row in patient_celltypes}

width_mm = 183.0
height_mm = 174.0
fig = plt.figure(figsize=(width_mm / 25.4, height_mm / 25.4), constrained_layout=False)
outer = fig.add_gridspec(
    2,
    2,
    left=0.13,
    right=0.965,
    bottom=0.075,
    top=0.885,
    width_ratios=[0.96, 1.24],
    height_ratios=[0.80, 1.18],
    wspace=0.46,
    hspace=0.46,
)

fig.suptitle(
    "Supplementary Fig. S3 |\nSingle-cell cohort composition and frozen annotation structure in GSE182434",
    x=0.055,
    y=0.982,
    ha="left",
    va="top",
    fontsize=9.6,
    fontweight="bold",
    color=DARK,
)

# Panel A
ax_a = fig.add_subplot(outer[0, 0])
ax_a.set_xlim(0, 1)
ax_a.set_ylim(0, 1)
ax_a.axis("off")
add_panel_label(ax_a, "A", "GSE182434 cohort composition\nand primary-analysis boundary", 7.4)

def block(ax, x, y, w, h, face, edge, title, value, subtitle="", title_color=DARK, value_color=NAVY):
    patch = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.015", facecolor=face, edgecolor=edge, linewidth=0.9)
    ax.add_patch(patch)
    ax.text(x + 0.025, y + h - 0.055, title, ha="left", va="top", fontsize=7.4, fontweight="bold", color=title_color)
    ax.text(x + 0.025, y + h * 0.47, value, ha="left", va="center", fontsize=10.2, fontweight="bold", color=value_color)
    if subtitle:
        ax.text(x + 0.025, y + 0.035, subtitle, ha="left", va="bottom", fontsize=6.5, color=MID_GRAY)

block(ax_a, 0.01, 0.40, 0.25, 0.36, LIGHT_BLUE, MID_BLUE, "Full GSE182434", "28,416", "cells | frozen annotation")
ax_a.annotate("", xy=(0.34, 0.58), xytext=(0.26, 0.58), arrowprops=dict(arrowstyle="-|>", lw=1.1, color=MID_GRAY))

components = [
    ("DLBCL", "14,368", PALE_BLUE, BLUE, NAVY),
    ("FL", "10,219", PALE_GRAY, GRID, DARK),
    ("Tonsil", "3,829", PALE_GRAY, GRID, DARK),
]
ys = [0.66, 0.45, 0.24]
for (label, value, face, edge, value_color), y in zip(components, ys):
    patch = FancyBboxPatch((0.34, y), 0.22, 0.14, boxstyle="round,pad=0.01,rounding_size=0.014", facecolor=face, edgecolor=edge, linewidth=0.9)
    ax_a.add_patch(patch)
    ax_a.text(0.45, y + 0.098, label, ha="center", va="center", fontsize=6.4, fontweight="bold", color=DARK)
    ax_a.text(0.45, y + 0.043, value, ha="center", va="center", fontsize=7.8, fontweight="bold", color=value_color)

ax_a.annotate("", xy=(0.67, 0.73), xytext=(0.56, 0.73), arrowprops=dict(arrowstyle="-|>", lw=1.2, color=BLUE))
block(
    ax_a,
    0.67,
    0.54,
    0.32,
    0.34,
    "#E4F0F5",
    BLUE,
    "Primary DLBCL subset",
    "14,368",
    "cells | 4 patients | 8 samples",
    value_color=NAVY,
)
ax_a.text(0.01, 0.08, "Primary single-cell analyses were restricted to the DLBCL subset.", fontsize=6.2, color=NAVY, fontweight="bold", ha="left")

# Panel B
ax_b = fig.add_subplot(outer[0, 1])
add_panel_label(ax_b, "B", "Cell representation across DLBCL\npatients and samples")
patient_y = np.arange(len(PATIENTS))[::-1]
sample_colors = {"B": "#4E8098", "NB": "#A7C5D3"}
for y, patient in zip(patient_y, PATIENTS):
    left = 0
    patient_samples = [key[1] for key in sample_counts if key[0] == patient]
    for sample in sorted(patient_samples, key=lambda value: value.endswith("NB")):
        n = sample_counts[(patient, sample)]
        suffix = "NB" if sample.endswith("NB") else "B"
        ax_b.barh(y, n, left=left, height=0.58, color=sample_colors[suffix], edgecolor="white", linewidth=0.8)
        if n >= 750:
            ax_b.text(left + n / 2, y, suffix, ha="center", va="center", fontsize=6.2, color="white" if suffix == "B" else DARK, fontweight="bold")
        left += n
    ax_b.text(left + 90, y, f"{left:,}", va="center", ha="left", fontsize=7.0, color=DARK, fontweight="bold")
ax_b.set_yticks(patient_y, PATIENTS)
ax_b.set_xlabel("Cells")
ax_b.set_xlim(0, 6100)
ax_b.grid(axis="x", color=GRID, linewidth=0.6)
ax_b.set_axisbelow(True)
ax_b.legend(
    handles=[
        Rectangle((0, 0), 1, 1, facecolor=sample_colors["B"], label="B-cell sorted sample"),
        Rectangle((0, 0), 1, 1, facecolor=sample_colors["NB"], label="Non-B-cell sorted sample"),
    ],
    loc="lower right",
    bbox_to_anchor=(1.0, 1.005),
    ncol=2,
    handlelength=1.0,
    columnspacing=1.0,
)
ax_b.text(0.0, -0.23, "Cells are nested within patients and samples; counts are descriptive.", transform=ax_b.transAxes, fontsize=6.5, color=MID_GRAY)

# Panel C
ax_c = fig.add_subplot(outer[1, 0])
add_panel_label(ax_c, "C", "Cell-type composition of the primary DLBCL subset")
order_c = sorted(CELLTYPES, key=lambda value: celltype_counts[value], reverse=True)
y_c = np.arange(len(order_c))[::-1]
values_c = [celltype_counts[value] for value in order_c]
ax_c.barh(y_c, values_c, height=0.64, color="#6F9FB5", edgecolor="white", linewidth=0.7)
labels_c = ["Monocytes and\nMacrophages" if value == "Monocytes and Macrophages" else value for value in order_c]
ax_c.set_yticks(y_c, labels_c)
ax_c.set_xlabel("Cells")
ax_c.set_xlim(0, 6700)
ax_c.grid(axis="x", color=GRID, linewidth=0.6)
ax_c.set_axisbelow(True)
for y, celltype in zip(y_c, order_c):
    n = celltype_counts[celltype]
    ax_c.text(n + 90, y, f"{n:,} ({celltype_pct[celltype]:.1f}%)", va="center", ha="left", fontsize=6.7, color=DARK)

# Panel D
ax_d = fig.add_subplot(outer[1, 1])
add_panel_label(ax_d, "D", "Frozen cell-type composition\nacross DLBCL patients")
matrix_pct = np.array([[pt_pct[(patient, celltype)] for patient in PATIENTS] for celltype in CELLTYPES])
matrix_n = np.array([[pt_counts[(patient, celltype)] for patient in PATIENTS] for celltype in CELLTYPES])
cmap = LinearSegmentedColormap.from_list("s3blue", ["#F5F8FA", "#D5E6EE", "#85B2C6", "#2F6B8A", "#173F5F"])
vmax = 60.0
im = ax_d.pcolormesh(
    np.arange(len(PATIENTS) + 1) - 0.5,
    np.arange(len(CELLTYPES) + 1) - 0.5,
    matrix_pct,
    cmap=cmap,
    vmin=0,
    vmax=vmax,
    shading="flat",
    edgecolors="white",
    linewidth=1.0,
)
ax_d.set_xlim(-0.5, len(PATIENTS) - 0.5)
ax_d.set_ylim(len(CELLTYPES) - 0.5, -0.5)
ax_d.set_xticks(np.arange(len(PATIENTS)), PATIENTS)
labels_d = ["Monocytes and\nMacrophages" if value == "Monocytes and Macrophages" else value for value in CELLTYPES]
ax_d.set_yticks(np.arange(len(CELLTYPES)), labels_d)
ax_d.tick_params(axis="x", top=True, bottom=False, labeltop=True, labelbottom=False, length=0, pad=5)
ax_d.tick_params(axis="y", length=0)
for i in range(len(CELLTYPES)):
    for j in range(len(PATIENTS)):
        value = matrix_pct[i, j]
        text_color = "white" if value >= 28 else DARK
        ax_d.text(j, i, f"{matrix_n[i, j]:,}", ha="center", va="center", fontsize=6.8, color=text_color, fontweight="bold" if value >= 28 else "normal")
cbar = fig.colorbar(im, ax=ax_d, fraction=0.042, pad=0.035)
cbar.solids.set_rasterized(False)
cbar.solids.set_edgecolor("face")
cbar.solids.set_linewidth(0)
cbar.set_label("Within-patient cells (%)", rotation=270, labelpad=13)
cbar.set_ticks([0, 15, 30, 45, 60])
cbar.ax.set_yticklabels(["0", "15", "30", "45", "60"])
ax_d.text(0.0, -0.105, "Color: within-patient proportion  |  Labels: raw cell counts", transform=ax_d.transAxes, fontsize=6.7, color=MID_GRAY)
for spine in ax_d.spines.values():
    spine.set_visible(False)

visual = ROOT / "visual_qc"
visual.mkdir(parents=True, exist_ok=True)
export = sys.argv[1] if len(sys.argv) > 1 else "draft"
draft_name = "Supplementary_Fig_S3_DRAFT2.png" if export == "draft2" else "Supplementary_Fig_S3_DRAFT1.png"
draft = visual / draft_name
fig.savefig(draft, dpi=300, facecolor="white")

if export == "final":
    pdf = ROOT / "Supplementary_Fig_S3_FINAL_SUBMISSION.pdf"
    svg = ROOT / "Supplementary_Fig_S3_FINAL_SUBMISSION.svg"
    png = ROOT / "Supplementary_Fig_S3_FINAL_SUBMISSION.png"
    tiff = ROOT / "Supplementary_Fig_S3_FINAL_SUBMISSION.tiff"
    fig.savefig(pdf, facecolor="white")
    fig.savefig(svg, facecolor="white")
    fig.savefig(png, dpi=600, facecolor="white")
    rgba = Image.open(png).convert("RGB")
    rgba.save(tiff, format="TIFF", compression="tiff_lzw", dpi=(600, 600))
    metadata = {
        "status": "PASS",
        "backend": "Python/matplotlib",
        "matplotlib_version": mpl.__version__,
        "width_mm": width_mm,
        "height_mm": height_mm,
        "png_dimensions": [rgba.width, rgba.height],
        "tiff_mode": "RGB",
        "tiff_compression": "LZW",
        "tiff_dpi": [600, 600],
    }
    (ROOT / "technical_qc/FINAL_EXPORT_METADATA.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metadata, indent=2))
else:
    print(json.dumps({"status": f"{draft.stem}_CREATED", "path": str(draft), "width_mm": width_mm, "height_mm": height_mm}, indent=2))

plt.close(fig)
