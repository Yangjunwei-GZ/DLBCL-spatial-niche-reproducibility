from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_ROOT = SCRIPT_DIR.parent
LOCAL_LIB = OUTPUT_ROOT / "python_lib"
sys.path.insert(0, str(LOCAL_LIB))

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
import numpy as np
from PIL import Image


REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
PROJECT_ROOT = Path(os.environ.get("DLBCL_PROJECT_ROOT", REPOSITORY_ROOT))
AUTHORITY_ROOT = PROJECT_ROOT / "revision_2026_reviewer_response/06y_suppfig_s2_discrete_taxonomy_nonretention_authority"
DATA_ROOT = AUTHORITY_ROOT / "data_authority"
MATRIX_ROOT = DATA_ROOT / "consensus_matrices"
AUTHORITY_QC = AUTHORITY_ROOT / "qc"

SUMMARY_PATH = DATA_ROOT / "SUPPLEMENTARY_FIG_S2_K2_K6_MODEL_ADJUDICATION_AUTHORITY.csv"
CDF_PATH = DATA_ROOT / "SUPPLEMENTARY_FIG_S2_CONSENSUS_CDF_DELTA_AREA_PLOTTING.csv"
THRESHOLD_PATH = DATA_ROOT / "SUPPLEMENTARY_FIG_S2_PRESPECIFIED_ACCEPTANCE_THRESHOLDS.csv"
AUTHORITY_REPORT = AUTHORITY_ROOT / "SUPPLEMENTARY_FIG_S2_DATA_AUTHORITY_REPORT.md"
AUTHORITY_SELF_CHECK = AUTHORITY_QC / "SUPPLEMENTARY_FIG_S2_AUTHORITY_SELF_CHECK.csv"
MATRIX_QC_PATH = AUTHORITY_QC / "SUPPLEMENTARY_FIG_S2_CONSENSUS_MATRIX_QC.csv"
PYTHON_AUTHORITY_QC = AUTHORITY_QC / "SUPPLEMENTARY_FIG_S2_INDEPENDENT_PYTHON_VALIDATION.json"

EXPECTED_HASHES = {
    SUMMARY_PATH: "e06a95bc745c8d5d63234bd131f29b8c1ef20049f29096e77cdf5f066e42dbb0",
    CDF_PATH: "ba37b2eb5227c963c6934fb120e528464f2153b2d2dce59b8de2e0d458011f5c",
    THRESHOLD_PATH: "9ae40a4200a048450eafc1c8e30aa27ef719a04b1a015d7f02581d05aaa2fba1",
}

CONCLUSION_LINE_1 = "No candidate k satisfied all prespecified retention criteria."
CONCLUSION_LINE_2 = "No discrete ecosystem taxonomy was assigned."
GOVERNING_CONCLUSION = f"{CONCLUSION_LINE_1[:-1]}; {CONCLUSION_LINE_2[0].lower() + CONCLUSION_LINE_2[1:]}"

BLUE_SCALE = ["#D6E6EF", "#B9D4E2", "#8CB8CD", "#568FAE", "#1E5A78"]
LINE_COLORS = ["#244C66", "#416E88", "#638EA5", "#87AFC1", "#A8C8D5"]
NAVY = "#173B56"
TEXT = "#22272B"
MID = "#687078"
LIGHT = "#EDF1F3"
GRID = "#D5DBDF"
PASS_FILL = "#9FC3D5"
FAIL_FILL = "#ECEFF1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def as_bool(value: str) -> bool:
    return str(value).upper() == "TRUE"


def verify_authority() -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]], dict[int, np.ndarray], dict[int, list[dict[str, str]]], dict[str, object]]:
    for path, expected in EXPECTED_HASHES.items():
        if not path.is_file():
            raise RuntimeError(f"Missing frozen authority: {path}")
        observed = sha256(path)
        if observed != expected:
            raise RuntimeError(f"Frozen authority hash mismatch: {path}; expected={expected}; observed={observed}")

    summary = read_csv(SUMMARY_PATH)
    cdf = read_csv(CDF_PATH)
    thresholds = read_csv(THRESHOLD_PATH)
    self_check = read_csv(AUTHORITY_SELF_CHECK)
    matrix_qc = read_csv(MATRIX_QC_PATH)
    independent_qc = json.loads(PYTHON_AUTHORITY_QC.read_text(encoding="utf-8"))
    authority_report = AUTHORITY_REPORT.read_text(encoding="utf-8")

    if len(summary) != 5 or [int(row["k"]) for row in summary] != [2, 3, 4, 5, 6]:
        raise RuntimeError("Summary authority does not contain exactly k=2-6")
    if any(as_bool(row["all_criteria_met"]) for row in summary):
        raise RuntimeError("A candidate k unexpectedly passed the frozen criteria")
    if [int(row["failed_criterion_count"]) for row in summary] != [3, 5, 6, 5, 5]:
        raise RuntimeError("Failed-criterion counts differ from frozen authority")
    if any("NOT_SELECTED" not in row["final_decision"] or "NOT_ASSIGNED" not in row["final_decision"] for row in summary):
        raise RuntimeError("Final decision fields do not retain NOT_SELECTED/NOT_ASSIGNED")
    if len(cdf) != 505 or any(sum(int(row["k"]) == k for row in cdf) != 101 for k in range(2, 7)):
        raise RuntimeError("CDF authority must contain 101 points per k")
    if len(thresholds) != 6 or any(as_bool(row["criteria_redefined"]) for row in thresholds):
        raise RuntimeError("Threshold authority is incomplete or marked redefined")
    if any(row["status"] != "PASS" for row in self_check) or len(self_check) != 18:
        raise RuntimeError("06y authority self-check is not 18/18 PASS")
    if len(matrix_qc) != 5 or any(row["assignments_match"] != "TRUE" for row in matrix_qc):
        raise RuntimeError("06y matrix QC is incomplete")
    if independent_qc.get("status") != "PASS" or independent_qc.get("scientific_reanalysis_performed") is not False:
        raise RuntimeError("06y independent validation is not PASS")
    if GOVERNING_CONCLUSION not in authority_report:
        raise RuntimeError("Governing conclusion is absent from the authority report")

    matrices: dict[int, np.ndarray] = {}
    orders: dict[int, list[dict[str, str]]] = {}
    matrix_hashes: dict[str, str] = {}
    for k in range(2, 7):
        matrix_path = MATRIX_ROOT / f"SUPPLEMENTARY_FIG_S2_CONSENSUS_MATRIX_K{k}_PLOTTING.csv.gz"
        order_path = MATRIX_ROOT / f"SUPPLEMENTARY_FIG_S2_CONSENSUS_ORDER_K{k}.csv"
        order = read_csv(order_path)
        with gzip.open(matrix_path, "rt", encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle)
            header = next(reader)
            rows = list(reader)
        if len(rows) != 498 or len(header) != 499 or any(len(row) != 499 for row in rows):
            raise RuntimeError(f"Consensus matrix k={k} is not 498 x 498")
        if len(order) != 498 or [row["sample_id"] for row in order] != header[1:] or [row[0] for row in rows] != header[1:]:
            raise RuntimeError(f"Consensus matrix k={k} sample order differs from frozen order")
        matrix = np.asarray([[float(value) for value in row[1:]] for row in rows], dtype=float)
        if not np.isfinite(matrix).all() or matrix.min() < 0 or matrix.max() > 1:
            raise RuntimeError(f"Consensus matrix k={k} contains invalid values")
        if not np.allclose(matrix, matrix.T, rtol=0, atol=1e-14) or not np.array_equal(np.diag(matrix), np.ones(498)):
            raise RuntimeError(f"Consensus matrix k={k} symmetry or diagonal check failed")
        matrices[k] = matrix
        orders[k] = order
        matrix_hashes[str(matrix_path)] = sha256(matrix_path)
        matrix_hashes[str(order_path)] = sha256(order_path)

    authority_qc = {
        "status": "PASS",
        "summary_rows": len(summary),
        "cdf_rows": len(cdf),
        "threshold_rows": len(thresholds),
        "consensus_matrices": 5,
        "matrix_dimensions": "498 x 498 each",
        "candidate_k": [2, 3, 4, 5, 6],
        "all_candidate_k_fail": True,
        "failed_criteria_counts": [3, 5, 6, 5, 5],
        "final_k": "NOT_SELECTED",
        "taxonomy": "NOT_ASSIGNED",
        "scientific_reanalysis_performed": False,
        "input_hashes": {str(path): sha256(path) for path in EXPECTED_HASHES} | matrix_hashes,
    }
    return summary, cdf, thresholds, matrices, orders, authority_qc


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "font.size": 6.4,
            "axes.titlesize": 7.2,
            "axes.labelsize": 6.5,
            "xtick.labelsize": 5.8,
            "ytick.labelsize": 5.8,
            "legend.fontsize": 5.8,
            "axes.linewidth": 0.65,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "text.color": TEXT,
            "axes.labelcolor": TEXT,
            "axes.edgecolor": TEXT,
            "xtick.color": TEXT,
            "ytick.color": TEXT,
        }
    )


def panel_label(ax: mpl.axes.Axes, letter: str, title: str, x: float = -0.08) -> None:
    ax.text(x, 1.09, letter, transform=ax.transAxes, fontsize=8.3, fontweight="bold", va="bottom", ha="left", clip_on=False)
    ax.text(0.0, 1.09, title, transform=ax.transAxes, fontsize=7.3, fontweight="bold", va="bottom", ha="left", clip_on=False)


def plot_panel_a(fig: mpl.figure.Figure, sub_spec: mpl.gridspec.SubplotSpec, matrices: dict[int, np.ndarray], orders: dict[int, list[dict[str, str]]]) -> list[mpl.axes.Axes]:
    cmap = LinearSegmentedColormap.from_list("consensus_blue", ["#FFFFFF", "#DCEAF1", "#9DC3D4", "#4F87A4", "#173F5F"])
    grid = sub_spec.subgridspec(1, 6, width_ratios=[1, 1, 1, 1, 1, 0.055], wspace=0.08)
    axes: list[mpl.axes.Axes] = []
    image = None
    for index, k in enumerate(range(2, 7)):
        ax = fig.add_subplot(grid[0, index])
        axes.append(ax)
        image = ax.imshow(matrices[k], cmap=cmap, vmin=0, vmax=1, interpolation="nearest", rasterized=True, aspect="equal")
        ax.set_title(f"k = {k}", fontsize=6.8, fontweight="normal", pad=2.2)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(0.45)
            spine.set_color("#7B8388")
        partition = [row["neutral_partition_id"] for row in orders[k]]
        transitions = [i for i in range(1, len(partition)) if partition[i] != partition[i - 1]]
        for boundary in transitions:
            ax.axhline(boundary - 0.5, color="#374047", linewidth=0.28, alpha=0.75)
            ax.axvline(boundary - 0.5, color="#374047", linewidth=0.28, alpha=0.75)
    cax = fig.add_subplot(grid[0, 5])
    colorbar = fig.colorbar(image, cax=cax, ticks=[0, 0.5, 1])
    colorbar.set_label("Consensus", fontsize=5.7, labelpad=2)
    colorbar.ax.tick_params(labelsize=5.3, length=2, width=0.5)
    colorbar.outline.set_linewidth(0.45)
    panel_label(axes[0], "A", "Consensus matrices, candidate k = 2-6", x=-0.23)
    return axes


def plot_panel_b(ax: mpl.axes.Axes, cdf: list[dict[str, str]]) -> None:
    panel_label(ax, "B", "Consensus CDF")
    for color, k in zip(LINE_COLORS, range(2, 7)):
        rows = [row for row in cdf if int(row["k"]) == k]
        x = np.asarray([float(row["consensus"]) for row in rows])
        y = np.asarray([float(row["CDF"]) for row in rows])
        ax.plot(x, y, color=color, linewidth=1.15, label=f"k={k}")
    ax.set_xlabel("Consensus index")
    ax.set_ylabel("CDF")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.grid(color=GRID, linewidth=0.45, alpha=0.7)
    ax.legend(ncol=5, loc="lower right", columnspacing=0.8, handlelength=1.5, handletextpad=0.35, borderaxespad=0.2)


def plot_panel_c(fig: mpl.figure.Figure, sub_spec: mpl.gridspec.SubplotSpec, summary: list[dict[str, str]]) -> list[mpl.axes.Axes]:
    grid = sub_spec.subgridspec(1, 2, wspace=0.40)
    axes = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])]
    panel_label(axes[0], "C", "CDF area and delta area", x=-0.30)
    k = np.asarray([int(row["k"]) for row in summary])
    area = np.asarray([float(row["consensus_CDF_area"]) for row in summary])
    delta = np.asarray([np.nan if row["delta_area"] == "" else float(row["delta_area"]) for row in summary])
    for ax, values, ylabel, subtitle in [
        (axes[0], area, "CDF area", "CDF area"),
        (axes[1], delta, "Delta area", "Delta area"),
    ]:
        ax.plot(k, values, color=NAVY, linewidth=0.9, zorder=1)
        ax.scatter(k, values, s=17, c=LINE_COLORS, edgecolors="white", linewidths=0.45, zorder=2)
        ax.set_title(subtitle, fontsize=6.4, fontweight="normal", pad=2)
        ax.set_xlabel("Candidate k")
        ax.set_ylabel(ylabel)
        ax.set_xticks(k)
        ax.grid(axis="y", color=GRID, linewidth=0.45, alpha=0.7)
    axes[0].set_ylim(0.30, 0.84)
    axes[1].set_ylim(0, 0.25)
    axes[1].text(2, 0.012, "NA", color=MID, fontsize=6.0, ha="center", va="bottom", fontweight="bold")
    return axes


def draw_cell(ax: mpl.axes.Axes, x: float, y: float, width: float, height: float, face: str, edge: str = "white", linewidth: float = 0.7) -> None:
    ax.add_patch(Rectangle((x, y), width, height, facecolor=face, edgecolor=edge, linewidth=linewidth))


def plot_panel_d(ax: mpl.axes.Axes, summary: list[dict[str, str]]) -> None:
    panel_label(ax, "D", "Quantitative model diagnostics")
    metrics = [
        ("PAC", "PAC", lambda value: f"{value:.3f}"),
        ("Mean silhouette", "silhouette", lambda value: f"{value:.3f}"),
        ("Minimum cluster size", "minimum_cluster_size", lambda value: f"{value:.0f}"),
        ("Minimum cluster\nmedian Jaccard", "minimum_cluster_median_jaccard", lambda value: f"{value:.3f}"),
        ("Overall median\nJaccard", "overall_median_jaccard", lambda value: f"{value:.3f}"),
    ]
    n_rows = len(metrics)
    n_cols = len(summary)
    ax.set_xlim(-1.85, n_cols)
    ax.set_ylim(n_rows, 0)
    for row_index, (label, field, formatter) in enumerate(metrics):
        ax.text(-1.77, row_index + 0.5, label, ha="left", va="center", fontsize=5.8)
        for col_index, row in enumerate(summary):
            draw_cell(ax, col_index, row_index, 1, 1, "#F7F8F9" if row_index % 2 == 0 else "#EDF1F3", edge="#D8DDE0", linewidth=0.55)
            ax.text(col_index + 0.5, row_index + 0.5, formatter(float(row[field])), ha="center", va="center", fontsize=5.8, fontweight="normal")
    for col_index, row in enumerate(summary):
        ax.text(col_index + 0.5, -0.14, f"k={row['k']}", ha="center", va="bottom", fontsize=5.9, fontweight="bold")
    ax.text(-1.77, -0.14, "Frozen metric", ha="left", va="bottom", fontsize=5.7, color=MID)
    ax.set_axis_off()


def plot_panel_e(ax: mpl.axes.Axes, summary: list[dict[str, str]], thresholds: list[dict[str, str]]) -> None:
    panel_label(ax, "E", "Prespecified retention criteria adjudication")
    criteria = [
        ("PAC rank <=2", "PAC_criterion_met"),
        ("Mean silhouette >0.25", "silhouette_absolute_criterion_met"),
        ("Within-best silhouette\n(within 0.02)", "silhouette_within_best_criterion_met"),
        ("Minimum cluster size >=25", "minimum_cluster_size_criterion_met"),
        ("Minimum cluster median\nJaccard >=0.60", "minimum_cluster_jaccard_criterion_met"),
        ("Overall median\nJaccard >=0.75", "overall_jaccard_criterion_met"),
    ]
    if len(thresholds) != len(criteria):
        raise RuntimeError("Panel E criterion count differs from authority")
    n_rows = len(criteria)
    n_cols = len(summary)
    ax.set_xlim(-3.75, n_cols + 0.15)
    ax.set_ylim(n_rows + 2.0, -0.55)
    for row_index, (label, field) in enumerate(criteria):
        ax.text(-3.65, row_index + 0.5, label, ha="left", va="center", fontsize=5.15, linespacing=0.92)
        for col_index, row in enumerate(summary):
            passed = as_bool(row[field])
            draw_cell(ax, col_index, row_index, 1, 1, PASS_FILL if passed else FAIL_FILL, edge="white", linewidth=0.75)
            ax.text(col_index + 0.5, row_index + 0.5, "PASS" if passed else "FAIL", ha="center", va="center", fontsize=5.1, fontweight="bold", color=NAVY if passed else "#596168")
    for col_index, row in enumerate(summary):
        ax.text(col_index + 0.5, -0.20, f"k={row['k']}", ha="center", va="bottom", fontsize=5.9, fontweight="bold")
        ax.text(col_index + 0.5, n_rows + 0.27, row["failed_criterion_count"], ha="center", va="center", fontsize=6.0, fontweight="bold")
    ax.text(-3.65, n_rows + 0.27, "Failed criteria", ha="left", va="center", fontsize=5.6, fontweight="bold")

    box_y = n_rows + 0.78
    ax.add_patch(Rectangle((-3.65, box_y), n_cols + 3.70, 1.0, facecolor="#F4F7F9", edgecolor=NAVY, linewidth=0.75))
    ax.text(-3.43, box_y + 0.35, CONCLUSION_LINE_1, ha="left", va="center", fontsize=5.45, fontweight="bold", color=NAVY)
    ax.text(-3.43, box_y + 0.70, CONCLUSION_LINE_2, ha="left", va="center", fontsize=5.45, color=NAVY)
    ax.set_axis_off()


def build_figure(summary: list[dict[str, str]], cdf: list[dict[str, str]], thresholds: list[dict[str, str]], matrices: dict[int, np.ndarray], orders: dict[int, list[dict[str, str]]], version: str) -> mpl.figure.Figure:
    width_mm = 183.0
    height_mm = 205.0 if version == "draft1" else 190.0
    fig = plt.figure(figsize=(width_mm / 25.4, height_mm / 25.4), facecolor="white")
    fig.suptitle("Supplementary Fig. S2 | Extended diagnostics for candidate-k model adjudication", x=0.07, y=0.986, ha="left", va="top", fontsize=9.2, fontweight="bold", color=TEXT)
    grid = fig.add_gridspec(
        3,
        2,
        left=0.075,
        right=0.965,
        top=0.925,
        bottom=0.055,
        height_ratios=[1.30, 1.72, 2.70],
        width_ratios=[1.0, 1.13],
        hspace=0.62 if version == "draft1" else 0.48,
        wspace=0.28 if version == "draft1" else 0.25,
    )
    plot_panel_a(fig, grid[0, :], matrices, orders)
    ax_b = fig.add_subplot(grid[1, 0])
    plot_panel_b(ax_b, cdf)
    plot_panel_c(fig, grid[1, 1], summary)
    ax_d = fig.add_subplot(grid[2, 0])
    plot_panel_d(ax_d, summary)
    ax_e = fig.add_subplot(grid[2, 1])
    plot_panel_e(ax_e, summary, thresholds)
    return fig


def save_draft(fig: mpl.figure.Figure, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, facecolor="white")


def save_final(fig: mpl.figure.Figure) -> dict[str, object]:
    base = OUTPUT_ROOT / "Supplementary_Fig_S2_FINAL_SUBMISSION"
    pdf_path = base.with_suffix(".pdf")
    svg_path = base.with_suffix(".svg")
    png_path = base.with_suffix(".png")
    tiff_path = base.with_suffix(".tiff")
    fig.savefig(pdf_path, facecolor="white")
    fig.savefig(svg_path, facecolor="white")
    fig.savefig(png_path, dpi=600, facecolor="white")
    temporary_tiff = OUTPUT_ROOT / "technical_qc/Supplementary_Fig_S2_FINAL_SUBMISSION_RGBA_TEMP.tiff"
    fig.savefig(temporary_tiff, dpi=600, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    with Image.open(temporary_tiff) as rgba_image:
        rgb_image = Image.new("RGB", rgba_image.size, "white")
        if rgba_image.mode == "RGBA":
            rgb_image.paste(rgba_image, mask=rgba_image.getchannel("A"))
        else:
            rgb_image.paste(rgba_image.convert("RGB"))
        rgb_image.save(tiff_path, format="TIFF", compression="tiff_lzw", dpi=(600, 600))
    temporary_tiff.unlink()
    with Image.open(tiff_path) as image:
        tiff_qc = {
            "width_px": image.width,
            "height_px": image.height,
            "mode": image.mode,
            "format": image.format,
            "compression": image.info.get("compression"),
            "dpi": tuple(float(value) for value in image.info.get("dpi", ())),
            "has_alpha": image.mode in {"RGBA", "LA"},
        }
    return {
        "width_mm": 183.0,
        "height_mm": 190.0,
        "pdf": str(pdf_path),
        "svg": str(svg_path),
        "png": str(png_path),
        "tiff": str(tiff_path),
        "tiff_qc": tiff_qc,
    }


def write_legend(path: Path) -> None:
    legend = (
        "Supplementary Fig. S2 | Extended diagnostics for candidate-k model adjudication. "
        "(A) Frozen consensus matrices for candidate partitions k=2-6, shown in the corresponding frozen consensus order with a common 0-1 color scale. Neutral identifiers C1-Ck denote plotting partitions only and are not biological classes. "
        "(B) Frozen consensus cumulative distribution function (CDF) curves for k=2-6. "
        "(C) Frozen CDF area and delta-area diagnostics. Delta area for k=2 is a prespecified structural missing value and is displayed as NA without imputation. "
        "(D) Frozen PAC, mean silhouette, minimum cluster size, minimum cluster median Jaccard, and overall median Jaccard values. "
        "(E) Adjudication against the six prespecified retention criteria: PAC rank <=2, mean silhouette >0.25, silhouette within 0.02 of the best value in the same score space, minimum cluster size >=25, minimum cluster median Jaccard >=0.60, and overall median Jaccard >=0.75. "
        "Candidate k=2, 3, 4, 5, and 6 failed 3, 5, 6, 5, and 5 criteria, respectively. Lower PAC or apparent block structure does not override the prespecified retention criteria. "
        "No candidate k satisfied all prespecified retention criteria; therefore, no discrete ecosystem taxonomy was assigned. Final k remained NOT_SELECTED and taxonomy remained NOT_ASSIGNED."
    )
    path.write_text(legend + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["draft1", "draft2", "final"], required=True)
    args = parser.parse_args()
    configure_matplotlib()
    summary, cdf, thresholds, matrices, orders, authority_qc = verify_authority()
    (OUTPUT_ROOT / "technical_qc").mkdir(parents=True, exist_ok=True)
    (OUTPUT_ROOT / "technical_qc/AUTHORITY_PREFLIGHT.json").write_text(json.dumps(authority_qc, indent=2) + "\n", encoding="utf-8")
    fig = build_figure(summary, cdf, thresholds, matrices, orders, args.mode)
    if args.mode in {"draft1", "draft2"}:
        output = OUTPUT_ROOT / f"visual_qc/Supplementary_Fig_S2_{args.mode.upper()}.png"
        save_draft(fig, output)
        payload = {"status": f"{args.mode.upper()}_RENDERED", "output": str(output), "scientific_reanalysis_performed": False}
    else:
        exports = save_final(fig)
        write_legend(OUTPUT_ROOT / "Supplementary_Fig_S2_FINAL_LEGEND.txt")
        payload = {"status": "FINAL_RENDERED_PENDING_QC", "exports": exports, "scientific_reanalysis_performed": False}
        (OUTPUT_ROOT / "technical_qc/FINAL_EXPORT_METADATA.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    plt.close(fig)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
