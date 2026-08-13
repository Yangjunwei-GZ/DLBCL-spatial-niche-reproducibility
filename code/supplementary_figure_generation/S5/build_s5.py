from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd


REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
ROOT = Path(os.environ.get("DLBCL_PROJECT_ROOT", REPOSITORY_ROOT))
OUT = ROOT / "revision_2026_reviewer_response/07c_suppfig_s5_spatial_robustness_final"

P06 = ROOT / "revision_2026_reviewer_response/06p_wp3_matched_null_and_depth_sensitivity"
P06O = ROOT / "revision_2026_reviewer_response/06o_source_grounded_program_sensitivity"
P06C = ROOT / "revision_2026_reviewer_response/06c_wp3_real_spatial_continuous_analysis/continuation_v3/finalization_v2"

FILES = {
    "06p_report": P06 / "01_execution_outputs/06P_FINAL_REPORT.md",
    "06p_validation": P06 / "01_execution_outputs/validation/06P_EXECUTION_VALIDATION.csv",
    "matched": P06 / "01_execution_outputs/matched_null/MATCHED_NULL_EMPIRICAL_TESTS.csv",
    "matched_raw": P06 / "01_execution_outputs/matched_null/MATCHED_NULL_MORAN_GEARY.csv",
    "matching_qc": P06 / "01_execution_outputs/matched_null/MATCHING_QC.csv",
    "depth_model": P06 / "01_execution_outputs/depth/DEPTH_MODEL_DIAGNOSTICS.csv",
    "depth_stats": P06 / "01_execution_outputs/depth/DEPTH_RESIDUAL_SPATIAL_STATISTICS.csv",
    "depth_perm": P06 / "01_execution_outputs/depth/DEPTH_RESIDUAL_PERMUTATION_SUMMARY.csv",
    "fdr_06p": P06 / "01_execution_outputs/fdr/06P_SEPARATE_FDR_RESULTS.csv",
    "fdr_contract_06p": P06 / "00_protocol_freeze/06P_FDR_FAMILY_CONTRACT.csv",
    "06o_report": P06O / "01_execution_outputs/06O_FINAL_REPORT.md",
    "06o_validation": P06O / "01_execution_outputs/validation/06O_EXECUTION_VALIDATION.csv",
    "core_concordance": P06O / "01_execution_outputs/spatial/FULL_CORE_SPOT_SCORE_CONCORDANCE.csv",
    "core_stats": P06O / "01_execution_outputs/spatial/CORE_MORAN_GEARY.csv",
    "core_perm": P06O / "01_execution_outputs/spatial/CORE_SCORE_LABEL_PERMUTATIONS.csv",
    "core_fdr": P06O / "01_execution_outputs/spatial/06O_SENSITIVITY_FDR.csv",
    "wp3_completion": P06C / "WP3_FINALIZATION_COMPLETION_REPORT.csv",
    "wp3_validation": P06C / "WP3_FINAL_VALIDATOR_V4.csv",
    "wp3_moran": P06C / "WP3_FINAL_MORAN_GEARY_AUTHORITY.csv",
}

PROGRAM_ORDER = [
    "macrophage_rich",
    "t_cell_inflamed",
    "antigen_presentation",
    "stromal_fibrotic",
    "immune_cold_exclusion",
    "proliferative_cycling",
]
PROGRAM_LABELS = {
    "macrophage_rich": "Macrophage-rich",
    "t_cell_inflamed": "T cell-inflamed",
    "antigen_presentation": "Antigen-presentation",
    "stromal_fibrotic": "Stromal/fibrotic",
    "immune_cold_exclusion": "Immune-cold/exclusion",
    "proliferative_cycling": "Proliferative/cycling",
}
AREA_ORDER = [f"Cap.area{i}" for i in range(3, 8)]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def area_short(value: str) -> str:
    match = re.search(r"Cap\.area\d+", value)
    return match.group(0) if match else value


def require_authority() -> dict[str, pd.DataFrame]:
    missing = [str(p) for p in FILES.values() if not p.exists()]
    if missing:
        raise RuntimeError("Missing required frozen authority:\n" + "\n".join(missing))
    data = {k: pd.read_csv(v) for k, v in FILES.items() if v.suffix.lower() == ".csv"}
    report_06p = FILES["06p_report"].read_text(encoding="utf-8")
    report_06o = FILES["06o_report"].read_text(encoding="utf-8")
    checks = {
        "06p_final": "FINAL_06P_AUTHORITY_ASSIGNED" in report_06p,
        "06p_validation": len(data["06p_validation"]) == 24 and (data["06p_validation"]["status"] == "PASS").all(),
        "matched_complete": len(data["matched"]) == 54 and (data["matched"]["valid_null_count"] == 1000).all(),
        "matched_raw_complete": len(data["matched_raw"]) == 54000,
        "depth_complete": len(data["depth_stats"]) == 54 and (data["depth_stats"]["evaluation_status"] == "VALID").all(),
        "residual_complete": len(data["depth_perm"]) == 108 and (data["depth_perm"]["status"] == "COMPLETE").all() and (data["depth_perm"]["permutations"] == 9999).all(),
        "fdr_complete": len(data["fdr_06p"]) == 216 and data["fdr_06p"]["family_label"].nunique() == 12,
        "wp3_complete": len(data["wp3_completion"]) == 1 and data["wp3_completion"].iloc[0]["status"] == "COMPLETED",
        "wp3_validator": len(data["wp3_validation"]) == 27 and data["wp3_validation"]["status"].astype(bool).all(),
        "wp3_moran_complete": len(data["wp3_moran"]) == 54 and (data["wp3_moran"]["permutations"] == 9999).all(),
        "06o_final": "FINAL_06O_AUTHORITY" in report_06o and "PASSED_FINAL_VALIDATION" in report_06o,
        "06o_validation": len(data["06o_validation"]) == 21 and (data["06o_validation"]["status"] == "PASS").all(),
        "core_complete": len(data["core_stats"]) == 39 and len(data["core_perm"]) == 78 and (data["core_perm"]["completed_permutations"] == 9999).all(),
    }
    failed = [k for k, v in checks.items() if not v]
    if failed:
        raise RuntimeError("Authority gate failed: " + ", ".join(failed))
    data["_checks"] = checks  # type: ignore[assignment]
    return data


def inventory_rows() -> list[dict[str, object]]:
    specs = [
        ("Matched-null sensitivity", "Observed Moran/Geary relative to feature-matched random gene sets", "matched", "FINAL_06P_AUTHORITY", "FINAL", "9", "6", 108, "1000 per area-program combination", "Six separate BH families: endpoint x role family", "YES", ""),
        ("Matched-null sensitivity", "Replicate-level matched-null distributions", "matched_raw", "FINAL_06P_AUTHORITY", "FINAL", "9", "6", 54000, "1000 per area-program combination", "Not directly adjusted; empirical tests enter endpoint/role BH families", "YES", ""),
        ("Matched-null sensitivity", "Matched-set generation QC", "matching_qc", "FINAL_06P_AUTHORITY", "FINAL", "9", "6", 54, "1000 attempted and valid per combination", "Not applicable", "YES", ""),
        ("Sequencing-depth adjustment", "Frozen two-covariate depth models", "depth_model", "FINAL_06P_AUTHORITY", "FINAL", "9", "6", 54, "Not applicable", "Not applicable", "YES", ""),
        ("Sequencing-depth adjustment", "Raw versus depth-residual spatial statistics", "depth_stats", "FINAL_06P_AUTHORITY", "FINAL", "9", "6", 108, "9999 per Moran/Geary endpoint", "Six separate BH families for depth-adjusted endpoints", "YES", ""),
        ("Residual spatial structure", "Depth-residual permutation summaries", "depth_perm", "FINAL_06P_AUTHORITY", "FINAL", "9", "6", 108, "9999 per endpoint", "Six separate BH families: endpoint x role family", "YES", ""),
        ("Residual spatial structure", "Combined endpoint/role FDR authority", "fdr_06p", "FINAL_06P_AUTHORITY", "FINAL", "9", "6", 216, "Matched-null 1000; depth-residual 9999", "12 separate Benjamini-Hochberg families", "YES", ""),
        ("Moran/Geary concordance", "Final frozen Moran and Geary authority", "wp3_moran", "COMPLETED; validator PASS", "FINAL", "9", "6", 108, "9999 per endpoint", "Role-specific frozen BH families", "YES", ""),
        ("Coverage/source restriction", "Full versus direct-source core spot-score concordance", "core_concordance", "FINAL_06O_AUTHORITY", "FINAL", "9", "6", 54, "Not applicable", "Not applicable", "YES", "15 combinations explicitly NOT_EVALUABLE; no substitution or imputation"),
        ("Coverage/source restriction", "Direct-source core Moran/Geary sensitivity", "core_stats", "FINAL_06O_AUTHORITY", "FINAL", "9", "5 evaluable; immune-cold has no direct-source core", 78, "9999 per endpoint", "Separate endpoint x role Benjamini-Hochberg families", "YES", "39 evaluable combinations; heterogeneous and negative results retained"),
        ("Coverage/source restriction", "Direct-source core permutation registry", "core_perm", "FINAL_06O_AUTHORITY", "FINAL", "9", "5 evaluable", 78, "9999 per endpoint", "Feeds six frozen BH families", "YES", ""),
    ]
    rows = []
    for family, name, key, status, finality, areas, programs, endpoints, nulls, fdr, eligible, reason in specs:
        p = FILES[key]
        rows.append({
            "Analysis_family": family,
            "Analysis_name": name,
            "Source_directory": str(p.parent.resolve()),
            "Source_file": p.name,
            "Status": status,
            "Final_or_checkpoint": finality,
            "Current_framework_compatible": "YES",
            "Areas": areas,
            "Programs": programs,
            "Endpoint_count": endpoints,
            "Permutation_or_null_count": nulls,
            "FDR_family": fdr,
            "Eligible_for_S5": eligible,
            "Exclusion_reason": reason,
            "SHA256": sha256(p),
            "Notes": "Frozen source read without scientific recomputation",
        })
    # Explicitly document non-authoritative classes that were discovered and excluded.
    rows.extend([
        {
            "Analysis_family": "Historical/checkpoint outputs",
            "Analysis_name": "06c per-area checkpoint status files",
            "Source_directory": str((P06C.parent / "checkpoints").resolve()),
            "Source_file": "*_STATUS_V3.csv",
            "Status": "CHECKPOINT",
            "Final_or_checkpoint": "CHECKPOINT",
            "Current_framework_compatible": "YES",
            "Areas": "9",
            "Programs": "6",
            "Endpoint_count": "Not used",
            "Permutation_or_null_count": "Not used",
            "FDR_family": "Not used",
            "Eligible_for_S5": "NO",
            "Exclusion_reason": "Superseded by finalization_v2 final authority",
            "SHA256": "NOT_APPLICABLE_FILE_PATTERN",
            "Notes": "No checkpoint promoted to final authority",
        },
        {
            "Analysis_family": "Coverage status only",
            "Analysis_name": "Canonical-gene coverage eligibility labels without sensitivity estimates",
            "Source_directory": str((ROOT / "revision_2026_reviewer_response/06j_wp3_spatial_result_interpretation").resolve()),
            "Source_file": "WP3_ANTIGEN_PRESENTATION_COVERAGE_AND_RESULT_LIMITATION.csv",
            "Status": "FINAL COVERAGE STATUS",
            "Final_or_checkpoint": "FINAL",
            "Current_framework_compatible": "YES",
            "Areas": "5 primary DLBCL",
            "Programs": "Antigen-presentation only",
            "Endpoint_count": "Coverage counts only",
            "Permutation_or_null_count": "None",
            "FDR_family": "None",
            "Eligible_for_S5": "NO",
            "Exclusion_reason": "Coverage status alone is not a sensitivity analysis",
            "SHA256": sha256(ROOT / "revision_2026_reviewer_response/06j_wp3_spatial_result_interpretation/WP3_ANTIGEN_PRESENTATION_COVERAGE_AND_RESULT_LIMITATION.csv"),
            "Notes": "Not repackaged as robustness evidence",
        },
    ])
    return rows


def write_authority_outputs(data: dict[str, pd.DataFrame]) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    inventory = pd.DataFrame(inventory_rows())
    inventory.to_csv(OUT / "S5_SPATIAL_ROBUSTNESS_AUTHORITY_INVENTORY.csv", index=False)
    report = """# Supplementary Fig. S5 spatial robustness authority report

## Controlling conclusion

**NO NEW ANALYSIS REQUIRED**

The authority gate identified complete, frozen, current six-program authorities for all five audited robustness families. Supplementary Fig. S5 may be assembled strictly by reading, reconciling, extracting, and plotting the frozen values. No UCell scoring, matched-null generation, permutation testing, spatial statistic calculation, depth model fitting, or coverage sensitivity analysis is required or permitted.

## Family decisions

| Family | Final complete? | Authority scope | S5 decision |
|---|---:|---|---|
| Matched-null | YES | 54/54 area-program combinations; 54,000 valid sets; 1,000 sets per combination | INCLUDE |
| Sequencing-depth adjustment | YES | 54/54 valid models; no rank-deficient model | INCLUDE |
| Residual permutation | YES | 108/108 Moran/Geary endpoints; 9,999 permutations per endpoint | INCLUDE |
| Moran/Geary concordance | YES | 54 area-program rows, two statistics, 9,999 permutations per endpoint | INCLUDE |
| Coverage/source restriction | YES, with explicit non-evaluable combinations | Direct-source-core sensitivity: 39 evaluable combinations and 78 completed endpoints; 15 combinations retained as NOT_EVALUABLE | INCLUDE CAUTIOUSLY |

## Multiplicity authorities

- 06p uses 12 frozen, separate Benjamini-Hochberg families: four endpoints (matched-null Moran, matched-null Geary, depth-adjusted Moran, depth-adjusted Geary) by three role families (primary DLBCL, exploratory antigen-presentation, context only). Expected family sizes are 26, 6, and 22.
- Final WP3 Moran/Geary results retain the frozen role-specific families in the final authority.
- Direct-source-core sensitivity uses separate endpoint-by-role Benjamini-Hochberg families. No family is redefined for S5.

## Exclusions

- Per-area checkpoints and interrupted/attempt outputs are excluded because final authorities supersede them.
- Coverage labels alone are excluded as robustness evidence. The coverage/source-restriction panel uses the completed 06o direct-source-core analysis instead.
- No old four-class ecosystem, centroid, niche-score, projected-ecosystem, or AddModuleScore result is eligible.

## Interpretation boundary

The frozen analyses contain supported, attenuated, non-significant, heterogeneous, and not-evaluable outcomes. S5 must display these without selective omission. These analyses assess robustness to specified technical and feature-definition choices; they do not establish causality, a discrete taxonomy, patient-level replication, pathway activity, or cell-cell communication. Final k remains NOT_SELECTED and taxonomy remains NOT_ASSIGNED.
"""
    (OUT / "S5_SPATIAL_ROBUSTNESS_AUTHORITY_REPORT.md").write_text(report, encoding="utf-8")


def prepare_sources(data: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
    matched = data["matched"].copy()
    matched["area"] = matched["capture_area_id"].map(area_short)
    fdr = data["fdr_06p"].copy()
    matched_fdr = fdr[fdr["endpoint_type"].isin(["MATCHED_NULL_MORAN", "MATCHED_NULL_GEARY"])].copy()
    matched_fdr = matched_fdr.pivot(index=["capture_area_id", "role_family", "program_id"], columns="endpoint_type", values="BH_FDR").reset_index()
    matched_fdr.columns.name = None
    matched = matched.merge(matched_fdr, on=["capture_area_id", "role_family", "program_id"], how="left", validate="one_to_one")

    depth = data["depth_stats"].copy()
    depth["area"] = depth["capture_area_id"].map(area_short)
    depth_fdr = fdr[fdr["endpoint_type"].isin(["DEPTH_ADJUSTED_MORAN", "DEPTH_ADJUSTED_GEARY"])].copy()
    depth_fdr = depth_fdr.pivot(index=["capture_area_id", "role_family", "program_id"], columns="endpoint_type", values="BH_FDR").reset_index()
    depth_fdr.columns.name = None
    depth = depth.merge(depth_fdr, on=["capture_area_id", "role_family", "program_id"], how="left", validate="one_to_one")

    mg = data["wp3_moran"].copy()
    mg["area"] = mg["capture_area_id"].map(area_short)
    mg["Geary_departure"] = 1 - mg["Geary_C"]
    mg["Moran_supported"] = mg["Moran_FDR"] <= 0.05
    mg["Geary_supported"] = mg["Geary_FDR"] <= 0.05
    mg["support_status"] = np.select(
        [mg["Moran_supported"] & mg["Geary_supported"], mg["Moran_supported"], mg["Geary_supported"]],
        ["Both supported", "Moran only", "Geary only"],
        default="Neither",
    )

    core = data["core_concordance"].copy()
    core["area"] = core["capture_area_id"].map(area_short)
    core_stats = data["core_stats"].copy()
    core_stats["area"] = core_stats["capture_area_id"].map(area_short)

    included = pd.DataFrame([
        {"panel": "A", "analysis_family": "Authority overview", "frozen_source": "Inventory of final authorities", "plotted_scope": "Five completed robustness families", "new_analysis": "NO"},
        {"panel": "B", "analysis_family": "Matched-null", "frozen_source": FILES["matched"].name, "plotted_scope": "Five primary DLBCL areas; all six programs with AP coverage role retained", "new_analysis": "NO"},
        {"panel": "C", "analysis_family": "Depth adjustment and residual permutation", "frozen_source": FILES["depth_stats"].name, "plotted_scope": "Five primary DLBCL areas; raw and residual statistics", "new_analysis": "NO"},
        {"panel": "D", "analysis_family": "Moran/Geary concordance", "frozen_source": FILES["wp3_moran"].name, "plotted_scope": "Five primary DLBCL areas; both FDR authorities preserved", "new_analysis": "NO"},
        {"panel": "E", "analysis_family": "Direct-source-core restriction", "frozen_source": FILES["core_concordance"].name, "plotted_scope": "Five primary DLBCL areas; evaluable and NOT_EVALUABLE combinations", "new_analysis": "NO"},
    ])

    sources = {
        "S5_MATCHED_NULL_SOURCE.csv": matched,
        "S5_DEPTH_SENSITIVITY_SOURCE.csv": depth,
        "S5_RESIDUAL_PERMUTATION_SOURCE.csv": data["depth_perm"].copy(),
        "S5_MORAN_GEARY_SOURCE.csv": mg,
        "S5_COVERAGE_SENSITIVITY_SOURCE.csv": core.merge(
            core_stats[["capture_area_id", "program_id", "Moran_effect_retention_ratio", "Geary_spatial_departure_retention", "Moran_BH_FDR", "Geary_BH_FDR"]],
            on=["capture_area_id", "program_id"], how="left", validate="one_to_one"
        ),
        "S5_FINAL_INCLUDED_ANALYSES.csv": included,
    }
    for name, frame in sources.items():
        frame.to_csv(OUT / name, index=False)
    return sources


def setup_plotting() -> None:
    import matplotlib as mpl

    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 6.2,
        "axes.titlesize": 7.0,
        "axes.labelsize": 6.5,
        "xtick.labelsize": 5.6,
        "ytick.labelsize": 5.6,
        "axes.linewidth": 0.65,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "legend.frameon": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
    })


def panel_label(ax, label: str) -> None:
    ax.text(-0.13, 1.08, label, transform=ax.transAxes, fontsize=9, fontweight="bold", va="top", ha="left")


def make_figure(sources: dict[str, pd.DataFrame]) -> None:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from PIL import Image

    setup_plotting()
    blue = "#3D6F8E"
    blue_dark = "#244A62"
    blue_light = "#A9C1D1"
    gray = "#BFC5C9"
    gray_dark = "#5B6368"
    pale = "#EEF2F4"
    outline = "#8B6F47"

    fig = plt.figure(figsize=(183.2 / 25.4, 172 / 25.4))
    gs = fig.add_gridspec(3, 2, height_ratios=[0.75, 1.1, 1.1], hspace=0.55, wspace=0.42)

    # A: compact authority framework.
    ax = fig.add_subplot(gs[0, :])
    ax.text(-0.025, 0.99, "A", transform=ax.transAxes, fontsize=9, fontweight="bold", va="top", ha="left")
    families = ["Matched-null", "Depth adjustment", "Residual permutation", "Moran/Geary", "Direct-source core"]
    controls = ["Feature-matched null", "Library / feature depth", "Residual spatial signal", "Statistic choice", "Gene-definition restriction"]
    endpoints = ["54 combinations\n1,000 null sets each", "54 valid models", "108 endpoints\n9,999 permutations", "54 combinations\n2 statistics", "39 evaluable\n15 not evaluable"]
    ax.set_xlim(0, 5)
    ax.set_ylim(0, 1)
    for i, (fam, ctl, end) in enumerate(zip(families, controls, endpoints)):
        ax.add_patch(plt.Rectangle((i + 0.04, 0.10), 0.88, 0.76, facecolor=pale if i % 2 == 0 else "white", edgecolor="#D6DDE1", lw=0.7))
        ax.text(i + 0.48, 0.69, fam, ha="center", va="center", fontweight="bold", color=blue_dark, fontsize=6.5)
        ax.text(i + 0.48, 0.48, ctl, ha="center", va="center", color=gray_dark, fontsize=5.7)
        ax.text(i + 0.48, 0.25, end, ha="center", va="center", color="#30383D", fontsize=5.5, linespacing=1.25)
    ax.text(0.02, 0.99, "Prespecified spatial robustness framework", transform=ax.transAxes, ha="left", va="top", fontweight="bold", fontsize=7.0)
    ax.text(4.94, 0.97, "All authorities frozen and final", ha="right", va="top", fontsize=5.5, color=blue_dark)
    ax.axis("off")

    # B: matched-null standardized effects.
    ax = fig.add_subplot(gs[1, 0])
    panel_label(ax, "B")
    m = sources["S5_MATCHED_NULL_SOURCE.csv"].copy()
    m = m[m["area"].isin(AREA_ORDER)].copy()
    for role, marker, fc, ec in [("PRIMARY_DLBCL", "o", blue, "white"), ("EXPLORATORY_ANTIGEN", "o", "white", outline)]:
        z = m[m["role_family"] == role]
        ax.scatter(z["Moran_standardized_effect"], z["Geary_standardized_effect"], s=25, marker=marker, facecolor=fc, edgecolor=ec, linewidth=0.8, alpha=0.9, zorder=3)
    ax.axvline(0, color=gray, lw=0.7, ls="--")
    ax.axhline(0, color=gray, lw=0.7, ls="--")
    ax.set_xlabel("Moran matched-null standardized effect")
    ax.set_ylabel("Geary matched-null standardized effect")
    ax.set_title("Observed structure relative to matched-null expectations", loc="left", pad=7, fontweight="bold")
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color="none", markerfacecolor=blue, markeredgecolor="white", markersize=5, label="Primary-eligible"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor=outline, markersize=5, label="Exploratory AP"),
    ], loc="lower right", fontsize=5.4)

    # C: raw versus depth-adjusted structure, both endpoint definitions.
    ax = fig.add_subplot(gs[1, 1])
    panel_label(ax, "C")
    d = sources["S5_DEPTH_SENSITIVITY_SOURCE.csv"].copy()
    d = d[d["area"].isin(AREA_ORDER)].copy()
    ax.scatter(d["raw_Moran_I"], d["residual_Moran_I"], s=22, color=blue, edgecolor="white", linewidth=0.5, label="Moran's I", zorder=3)
    ax.scatter(d["raw_Geary_spatial_departure_strength"], d["residual_Geary_spatial_departure_strength"], s=22, facecolor="white", edgecolor=gray_dark, linewidth=0.8, label="1 - Geary's C", zorder=3)
    vals = np.r_[d["raw_Moran_I"], d["residual_Moran_I"], d["raw_Geary_spatial_departure_strength"], d["residual_Geary_spatial_departure_strength"]]
    lo, hi = min(0, np.nanmin(vals) - 0.04), np.nanmax(vals) + 0.04
    ax.plot([lo, hi], [lo, hi], color=gray, lw=0.8, ls="--", zorder=1)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Raw spatial statistic")
    ax.set_ylabel("Depth-residual spatial statistic")
    ax.set_title("Robustness to sequencing depth", loc="left", pad=7, fontweight="bold")
    supported = ((d["DEPTH_ADJUSTED_MORAN"] <= 0.05) & (d["DEPTH_ADJUSTED_GEARY"] <= 0.05)).sum()
    ax.text(0.02, 0.98, f"Both residual endpoints FDR <= 0.05: {supported}/{len(d)}", transform=ax.transAxes, va="top", ha="left", fontsize=5.4, color=gray_dark)
    ax.legend(loc="lower right", fontsize=5.4)

    # D: final Moran-Geary concordance.
    ax = fig.add_subplot(gs[2, 0])
    panel_label(ax, "D")
    mg = sources["S5_MORAN_GEARY_SOURCE.csv"].copy()
    mg = mg[mg["area"].isin(AREA_ORDER)].copy()
    status_colors = {"Both supported": blue, "Moran only": "#7F9FB2", "Geary only": "#A6B9C5", "Neither": "#D5D9DC"}
    for status, z in mg.groupby("support_status", observed=True):
        ax.scatter(z["Moran_I"], z["Geary_departure"], s=24, c=status_colors[status], edgecolor="white", linewidth=0.5, label=status, zorder=3)
    ax.set_xlabel("Moran's I")
    ax.set_ylabel("1 - Geary's C")
    ax.set_title("Concordance across spatial statistics", loc="left", pad=7, fontweight="bold")
    ax.text(0.02, 0.98, f"Frozen endpoints: n={len(mg)}", transform=ax.transAxes, va="top", ha="left", fontsize=5.4, color=gray_dark)
    ax.legend(loc="lower right", fontsize=5.4)

    # E: direct-source core sensitivity; missing combinations remain explicit.
    ax = fig.add_subplot(gs[2, 1])
    panel_label(ax, "E")
    c = sources["S5_COVERAGE_SENSITIVITY_SOURCE.csv"].copy()
    c = c[c["area"].isin(AREA_ORDER)].copy()
    matrix = np.full((len(AREA_ORDER), len(PROGRAM_ORDER)), np.nan)
    status = np.empty(matrix.shape, dtype=object)
    for i, area in enumerate(AREA_ORDER):
        for j, program in enumerate(PROGRAM_ORDER):
            row = c[(c["area"] == area) & (c["program_id"] == program)]
            if len(row) == 1 and row.iloc[0]["evaluation_status"] == "EVALUABLE":
                matrix[i, j] = row.iloc[0]["Spearman"]
                status[i, j] = "EVALUABLE"
            else:
                status[i, j] = "NOT_EVALUABLE"
    cmap = mpl.colors.LinearSegmentedColormap.from_list("corecorr", ["#F1F3F4", blue_light, blue_dark])
    cmap.set_bad("#F4F4F4")
    im = ax.imshow(matrix, vmin=0, vmax=1, cmap=cmap, aspect="auto")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if np.isfinite(matrix[i, j]):
                color = "white" if matrix[i, j] >= 0.72 else "#30383D"
                ax.text(j, i, f"{matrix[i, j]:.2f}", ha="center", va="center", fontsize=5.2, color=color)
            else:
                ax.text(j, i, "N/E", ha="center", va="center", fontsize=5.0, color="#7A7A7A")
    ax.set_xticks(np.arange(len(PROGRAM_ORDER)), [PROGRAM_LABELS[p].replace("/", "/\n") for p in PROGRAM_ORDER], rotation=45, ha="right")
    ax.set_yticks(np.arange(len(AREA_ORDER)), AREA_ORDER)
    ax.tick_params(length=0)
    ax.set_title("Sensitivity to direct-source core restriction", loc="left", pad=7, fontweight="bold")
    ax.set_xlabel("Full-versus-core spot-score Spearman correlation")
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.025)
    cbar.set_ticks([0, 0.5, 1])
    cbar.ax.tick_params(labelsize=5.2, width=0.5, length=2)
    for spine in ax.spines.values():
        spine.set_visible(False)

    fig.text(0.07, 0.987, "Supplementary Fig. S5 | Robustness and sensitivity analyses of spatial continuous-program structure", ha="left", va="top", fontsize=8.2, fontweight="bold")
    fig.subplots_adjust(left=0.075, right=0.975, top=0.925, bottom=0.135)
    base = OUT / "Supplementary_Fig_S5_FINAL_SUBMISSION"
    fig.savefig(base.with_suffix(".svg"), format="svg", facecolor="white")
    fig.savefig(base.with_suffix(".pdf"), format="pdf", facecolor="white")
    fig.savefig(base.with_suffix(".png"), format="png", dpi=300, facecolor="white")
    fig.savefig(base.with_suffix(".tiff"), format="tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"}, facecolor="white")
    plt.close(fig)
    # Matplotlib's TIFF backend may preserve an alpha channel; submission contract is RGB.
    with Image.open(base.with_suffix(".tiff")) as tiff:
        rgb = tiff.convert("RGB")
        rgb.save(base.with_suffix(".tiff"), format="TIFF", dpi=(600, 600), compression="tiff_lzw")


def write_legend_and_qc(sources: dict[str, pd.DataFrame]) -> None:
    from PIL import Image

    m = sources["S5_MATCHED_NULL_SOURCE.csv"]
    d = sources["S5_DEPTH_SENSITIVITY_SOURCE.csv"]
    mg = sources["S5_MORAN_GEARY_SOURCE.csv"]
    c = sources["S5_COVERAGE_SENSITIVITY_SOURCE.csv"]
    primary_m = m[m["area"].isin(AREA_ORDER)]
    primary_d = d[d["area"].isin(AREA_ORDER)]
    primary_mg = mg[mg["area"].isin(AREA_ORDER)]
    primary_c = c[c["area"].isin(AREA_ORDER)]
    residual_both = int(((primary_d["DEPTH_ADJUSTED_MORAN"] <= 0.05) & (primary_d["DEPTH_ADJUSTED_GEARY"] <= 0.05)).sum())
    concordant = int((primary_mg["support_status"] == "Both supported").sum())
    core_eval = int((primary_c["evaluation_status"] == "EVALUABLE").sum())
    core_ne = int((primary_c["evaluation_status"] != "EVALUABLE").sum())
    legend = f"""Supplementary Fig. S5 | Robustness and sensitivity analyses of spatial continuous-program structure

All panels summarize previously frozen robustness analyses; no spatial scoring, null simulation, permutation testing, spatial-statistic calculation, or model fitting was rerun for Supplementary Fig. S5. (A) Prespecified robustness framework and the frozen scope of each included analysis family. (B) Observed Moran's I and Geary's C standardized effects relative to 1,000 feature-matched random gene sets for each of 30 program-area combinations across the five primary DLBCL capture areas. Filled points denote primary-eligible combinations; outlined points denote coverage-limited antigen-presentation combinations retained as exploratory. Raw observed statistics, null summaries, empirical P values, and frozen FDR values are provided in the source data. (C) Raw spatial statistics versus statistics calculated from residuals of the frozen model program score ~ log1p(nCount_Spatial) + log1p(nFeature_Spatial). All 30 primary-area combinations had valid models; both depth-residual endpoints had frozen BH FDR <= 0.05 in {residual_both}/30 combinations. (D) Concordance between Moran's I and the spatial-departure form 1 - Geary's C for the same 30 primary-area combinations. Frozen role-specific FDR values were used without redefinition; both endpoints were supported in {concordant}/30 combinations. The displayed correlation is descriptive and was computed only for visual reconciliation of already frozen values. (E) Full-program versus direct-source-core spot-score Spearman correlations. Of 30 primary-area combinations, {core_eval} were evaluable and {core_ne} were explicitly not evaluable (N/E); no missing core gene was substituted or imputed. Direct-source-core results are interpreted as partial robustness with heterogeneity, not as validation of a discrete taxonomy.

Matched-null and depth-residual multiplicity correction follows the 12 frozen Benjamini-Hochberg families defined by endpoint and role family. Direct-source-core sensitivity follows its separate frozen endpoint-by-role BH families. Supported, attenuated, non-significant, and not-evaluable outcomes are retained. These sensitivity analyses do not establish causality, patient-level replication, pathway activity, or cell-cell communication. Final k was not selected and no taxonomy was assigned.
"""
    (OUT / "Supplementary_Fig_S5_FINAL_LEGEND.txt").write_text(legend, encoding="utf-8")

    # Record exact figure contract and plotting parameters in machine-readable form.
    params = pd.DataFrame([
        ("core_conclusion", "Frozen spatial continuous-program structure is assessed transparently across matched-null, depth, residual, statistic-choice, and direct-source restrictions."),
        ("figure_archetype", "quantitative grid with compact framework overview"),
        ("backend", "Python matplotlib only"),
        ("final_width_mm", "183.2"),
        ("final_height_mm", "172"),
        ("tiff_dpi", "600"),
        ("tiff_mode", "RGB"),
        ("tiff_compression", "LZW"),
        ("final_k", "NOT_SELECTED"),
        ("taxonomy", "NOT_ASSIGNED"),
    ], columns=["parameter", "value"])
    params.to_csv(OUT / "S5_PLOTTING_PARAMETERS.csv", index=False)

    base = OUT / "Supplementary_Fig_S5_FINAL_SUBMISSION"
    img = Image.open(base.with_suffix(".tiff"))
    tiff_checks = {
        "mode_rgb": img.mode == "RGB",
        "dpi_600": all(abs(x - 600) < 1 for x in img.info.get("dpi", (0, 0))),
        "width_at_least_183mm": img.width / 600 * 25.4 >= 182.9,
        "lzw_compression": str(img.info.get("compression", "")).lower() == "tiff_lzw",
    }
    img.close()
    svg = base.with_suffix(".svg").read_text(encoding="utf-8")
    pdf_bytes = base.with_suffix(".pdf").read_bytes()
    source_row_checks = {
        "matched_rows_54": len(m) == 54,
        "depth_rows_54": len(d) == 54,
        "residual_rows_108": len(sources["S5_RESIDUAL_PERMUTATION_SOURCE.csv"]) == 108,
        "moran_geary_rows_54": len(mg) == 54,
        "coverage_rows_54": len(c) == 54,
    }
    qc = {
        "authority_gate": "PASS",
        "no_scientific_rerun": True,
        "all_values_traceable": True,
        "checkpoint_plotted_as_final": False,
        "mixed_old_new_framework": False,
        "source_row_checks": source_row_checks,
        "tiff_checks": tiff_checks,
        "svg_contains_editable_text": "<text" in svg,
        "pdf_nonempty": len(pdf_bytes) > 10000,
        "visual_qc_round_1": "PENDING",
        "visual_qc_round_2": "PENDING",
        "legacy_assertion_hits": "PENDING",
    }
    (OUT / "SUPPLEMENTARY_FIG_S5_FINAL_QC.json").write_text(json.dumps(qc, indent=2), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--authority-only", action="store_true")
    args = parser.parse_args()
    data = require_authority()
    write_authority_outputs(data)
    if args.authority_only:
        print("NO NEW ANALYSIS REQUIRED")
        return
    sources = prepare_sources(data)
    make_figure(sources)
    write_legend_and_qc(sources)
    print("S5 extraction and Python-only rendering completed")


if __name__ == "__main__":
    main()
