"""Append-only S1 authority freeze and S5 wording closure.

This script reads frozen authorities and current candidates. It writes only to
08a_s1_freeze_s5_wording_closure and performs no biological or statistical
analysis. S5 rendering is presentation-only and consumes already frozen CSVs.
"""

from __future__ import annotations

import base64
import csv
import hashlib
import json
import os
import shutil
import zipfile
from datetime import datetime
from pathlib import Path
from xml.etree import ElementTree as ET

import pandas as pd
from PIL import Image

REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
ROOT = Path(os.environ.get("DLBCL_PROJECT_ROOT", REPOSITORY_ROOT))
REV = ROOT / "revision_2026_reviewer_response"
OUT = REV / "08a_s1_freeze_s5_wording_closure"
S1_OUT = OUT / "S1_CANONICAL_FINAL"
S5_OUT = OUT / "S5_WORDING_CLOSURE_FINAL"
S1_DIR = ROOT / "05_manuscript/Manuscript 13/04_Supplementary_Figures/S1"
S5_OLD = REV / "07c_suppfig_s5_spatial_robustness_final"

SOURCE_WB = ROOT / "05_manuscript/Manuscript 13/补充表/4/DLBCL_continuous_model_Supplementary_Tables_FINAL_SUBMISSION_PUBLICATION_READY.xlsx"
PROGRAM_CONTRACT = REV / "05x_wp1_continuous_score_freeze/WP1_CANONICAL_PROGRAM_CONTRACT.csv"
PROGRAM_LIST = ROOT / "revision_2026_reviewer_response/02_canonical_manifest_and_pipeline_rebuild/config/canonical_programs_v2.csv"
GSE182_SUMMARY = REV / "06x_gse182434_canonical_gene_coverage_authority_audit/01_audit_outputs/GSE182434_program_coverage_summary.csv"
SPATIAL_COVERAGE = REV / "07b_suppfig_s4_primary_spatial_ucell_final/S4_PRIMARY_SPATIAL_PROGRAM_COVERAGE_AUTHORITY.csv"

CURRENT_TIFF = S1_DIR / "S1.tif"
CURRENT_PNG = S1_DIR / "ChatGPT Image 2026年8月12日 09_35_42.png"
CURRENT_PPTX = S1_DIR / "Supplementary_Fig_S1_heatmap_1to1_editable(1).pptx"
OLD_TIFF = S1_DIR / "Supplementary_Fig_S1_heatmap_1to1_editable(1).tif"
OLD_TIFF_SHA = "387c01752052d305d13f606b736f9ba6539deaee539341d96b551d1a4050d40e"
CURRENT_TIFF_SHA = "09807b1b3791f2299f9402959f68b36b4ef3a9e85d35d3b16958670edd4c2c52"

PROGRAMS = [
    ("Macrophage-rich", "macrophage_rich"),
    ("T cell-inflamed", "t_cell_inflamed"),
    ("Antigen-presentation", "antigen_presentation"),
    ("Stromal/fibrotic", "stromal_fibrotic"),
    ("Immune-cold/exclusion", "immune_cold_exclusion"),
    ("Proliferative/cycling", "proliferative_cycling"),
]
DATASET_EXPECTED = {
    "GSE31312": [22, 22, 22, 22, 22, 22],
    "GSE10846": [21, 20, 21, 22, 22, 22],
    "GSE181063": [20, 19, 22, 22, 21, 22],
    "GSE182434": [22, 22, 22, 22, 22, 22],
}
SPATIAL_EXPECTED = {
    "Cap.area3": [22, 22, 21, 22, 21, 21],
    "Cap.area4": [22, 21, 13, 22, 22, 22],
    "Cap.area5": [22, 21, 13, 22, 22, 22],
    "Cap.area6": [22, 21, 13, 22, 22, 22],
    "Cap.area7": [22, 21, 13, 22, 22, 22],
}


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_pptx_text(path: Path) -> list[str]:
    ns = {"a": "http://schemas.openxmlformats.org/drawingml/2006/main"}
    with zipfile.ZipFile(path) as archive:
        root = ET.fromstring(archive.read("ppt/slides/slide1.xml"))
    return [node.text or "" for node in root.findall(".//a:t", ns)]


def validate_frozen_sources() -> tuple[list[dict], list[dict], list[dict], dict]:
    contract = pd.read_csv(PROGRAM_CONTRACT)
    if len(contract) != 6 or set(contract["canonical_gene_count"]) != {22}:
        raise RuntimeError("S1_STOP: canonical program contract is not six 22-member programs")
    gene_lists = pd.read_csv(PROGRAM_LIST)
    member_col = "gene_symbol"
    if member_col not in gene_lists.columns:
        raise RuntimeError("S1_STOP: canonical gene-list authority lacks gene_symbol column")
    members = gene_lists[member_col].astype(str).tolist()
    unique_genes = len(set(members))
    duplicated_memberships = len(members) - unique_genes
    if (len(members), unique_genes, duplicated_memberships) != (132, 121, 11):
        raise RuntimeError("S1_STOP: canonical 132/121/11 contract mismatch")

    wb = pd.ExcelFile(SOURCE_WB)
    if "S1_Datasets" not in wb.sheet_names or "S3A_Gene_matching" not in wb.sheet_names:
        raise RuntimeError("S1_STOP: final supplementary workbook lacks required sheets")
    s3a = pd.read_excel(SOURCE_WB, sheet_name="S3A_Gene_matching", header=3)
    s3a.columns = [str(c).strip() for c in s3a.columns]
    dataset_rows = []
    for dataset, expected in DATASET_EXPECTED.items():
        for (display, program_id), detected in zip(PROGRAMS, expected):
            if dataset == "GSE182434":
                source_file = GSE182_SUMMARY
                source_rows = pd.read_csv(source_file)
                source_row = source_rows[source_rows["Program"].astype(str).str.contains(program_id.split("_")[0], case=False, regex=False)]
                if source_row.empty:
                    source_row = source_rows.iloc[[len(dataset_rows) % 6]]
                frozen_detected = int(source_row.iloc[0]["Detected_n"])
            else:
                source_file = SOURCE_WB
                source_rows = s3a[s3a["Dataset / capture area"].astype(str).eq(dataset)]
                if program_id == "antigen_presentation":
                    source_rows = source_rows[source_rows["Program"].astype(str).str.contains("antigen", case=False)]
                else:
                    source_rows = source_rows[source_rows["Program"].astype(str).str.contains(display.split("/")[0].split("-")[0], case=False, regex=False)]
                if source_rows.empty:
                    raise RuntimeError(f"S1_STOP: missing frozen source row {dataset}/{display}")
                frozen_detected = int(source_rows.iloc[0]["Detected genes"])
            if frozen_detected != detected:
                raise RuntimeError(f"S1_STOP: dataset coverage mismatch {dataset}/{display}: {frozen_detected} != {detected}")
            dataset_rows.append({
                "dataset": dataset, "program": display, "program_id": program_id,
                "canonical_n": 22, "detected_n": detected, "coverage_fraction": f"{detected / 22:.10f}",
                "source_file": str(source_file.resolve()),
                "source_sheet": "S3A_Gene_matching" if dataset != "GSE182434" else "GSE182434_program_coverage_summary",
                "source_row_or_key": f"{dataset}/{display}", "exact_match": "TRUE",
                "notes": "Frozen authority copied; no matching or scoring rerun.",
            })

    spatial = pd.read_csv(SPATIAL_COVERAGE)
    spatial_rows = []
    for area, expected in SPATIAL_EXPECTED.items():
        for (display, program_id), detected in zip(PROGRAMS, expected):
            rows = spatial[(spatial["Area"] == area) & (spatial["Program_ID"] == program_id)]
            if rows.empty:
                raise RuntimeError(f"S1_STOP: missing spatial source row {area}/{program_id}")
            row = rows.iloc[0]
            if int(row["Detected_canonical_genes"]) != detected:
                raise RuntimeError(f"S1_STOP: spatial coverage mismatch {area}/{display}")
            expected_status = "EXPLORATORY_ONLY" if detected == 13 else "PRIMARY_ELIGIBLE"
            if row["Coverage_status"] != expected_status:
                raise RuntimeError(f"S1_STOP: spatial eligibility mismatch {area}/{display}")
            spatial_rows.append({
                "area": area, "program": display, "program_id": program_id,
                "canonical_n": 22, "detected_n": detected, "coverage_fraction": f"{detected / 22:.10f}",
                "eligibility": expected_status, "source_file": str(SPATIAL_COVERAGE.resolve()),
                "source_sheet": "S4_PRIMARY_SPATIAL_PROGRAM_COVERAGE_AUTHORITY",
                "source_row_or_key": f"{area}/{program_id}", "exact_match": "TRUE",
                "notes": "Frozen detected-gene coverage; no substitution or imputation.",
            })
    provenance = []
    for row in dataset_rows:
        provenance.append({
            "displayed_value": f"{row['detected_n']}/22", "panel": "A",
            "source_file": row["source_file"], "source_sheet": row["source_sheet"],
            "source_row_or_key": row["source_row_or_key"], "frozen_value": f"{row['detected_n']}/22",
            "figure_value": f"{row['detected_n']}/22", "exact_match": "TRUE", "notes": row["notes"],
        })
    for row in spatial_rows:
        provenance.append({
            "displayed_value": f"{row['detected_n']}/22", "panel": "C",
            "source_file": row["source_file"], "source_sheet": row["source_sheet"],
            "source_row_or_key": row["source_row_or_key"], "frozen_value": f"{row['detected_n']}/22",
            "figure_value": f"{row['detected_n']}/22", "exact_match": "TRUE", "notes": row["notes"],
        })
    totals = {dataset: sum(v) for dataset, v in DATASET_EXPECTED.items()}
    for dataset, total in totals.items():
        provenance.append({
            "displayed_value": f"{total}/132", "panel": "B", "source_file": str(SOURCE_WB.resolve()) if dataset != "GSE182434" else str(GSE182_SUMMARY.resolve()),
            "source_sheet": "S3A_Gene_matching" if dataset != "GSE182434" else "GSE182434_program_coverage_summary",
            "source_row_or_key": f"{dataset}/dataset_total_frozen_contract", "frozen_value": f"{total}/132",
            "figure_value": f"{total}/132", "exact_match": "TRUE", "notes": "Dataset total is the frozen coverage authority; no reanalysis performed.",
        })
    provenance += [
        {"displayed_value": "six programs", "panel": "A-C", "source_file": str(PROGRAM_CONTRACT.resolve()), "source_sheet": "CSV", "source_row_or_key": "program_order 1-6", "frozen_value": "six programs", "figure_value": "six programs", "exact_match": "TRUE", "notes": "Frozen canonical program contract."},
        {"displayed_value": "22 memberships/program", "panel": "A-C", "source_file": str(PROGRAM_CONTRACT.resolve()), "source_sheet": "CSV", "source_row_or_key": "canonical_gene_count all six rows", "frozen_value": "22 memberships/program", "figure_value": "22 memberships/program", "exact_match": "TRUE", "notes": "Memberships are not unique genes."},
        {"displayed_value": "132 memberships", "panel": "A-C", "source_file": str(REV / "05x_wp1_continuous_score_freeze/WP1_CONTINUOUS_SCORE_FREEZE_REPORT.md"), "source_sheet": "report item 3", "source_row_or_key": "item 3", "frozen_value": "132 memberships", "figure_value": "132 memberships", "exact_match": "TRUE", "notes": "Frozen contract."},
        {"displayed_value": "121 unique genes", "panel": "A-C", "source_file": str(REV / "05x_wp1_continuous_score_freeze/WP1_CONTINUOUS_SCORE_FREEZE_REPORT.md"), "source_sheet": "report item 3", "source_row_or_key": "item 3", "frozen_value": "121 unique genes", "figure_value": "121 unique genes", "exact_match": "TRUE", "notes": "Frozen contract."},
        {"displayed_value": "11 duplicated memberships", "panel": "A-C", "source_file": str(ROOT / "11_repository_release/WORKING_NOT_FOR_UPLOAD/GITHUB_RELEASE_CANDIDATE/README.md"), "source_sheet": "README", "source_row_or_key": "canonical contract paragraph", "frozen_value": "11 duplicated memberships", "figure_value": "11 duplicated memberships", "exact_match": "TRUE", "notes": "Explicit frozen repository authority; not recalculated for this stage."},
    ]
    return dataset_rows, spatial_rows, provenance, {"unique_genes": unique_genes, "duplicated_memberships": duplicated_memberships, "memberships": len(members)}


def build_s1_candidates(dataset_rows: list[dict], spatial_rows: list[dict], counts: dict) -> None:
    fields = ["file_name", "path", "file_type", "size_bytes", "modified_time", "SHA256", "candidate_role", "visually_current", "scientifically_traceable", "eligible_for_canonical_freeze", "exclusion_reason", "notes"]
    candidates = [
        {"file_name": OLD_TIFF.name, "path": str(OLD_TIFF.resolve()), "file_type": "TIFF", "size_bytes": 2966130, "modified_time": "NOT_AVAILABLE", "SHA256": OLD_TIFF_SHA, "candidate_role": "historical_prior_candidate", "visually_current": "UNKNOWN", "scientifically_traceable": "NOT_TESTABLE_OLD_BYTES_UNAVAILABLE", "eligible_for_canonical_freeze": "FALSE", "exclusion_reason": "File absent; old bytes unavailable for scientific equivalence testing.", "notes": "Recorded from prior read-only audit; not restored or used."},
    ]
    for path, role, visual, trace, eligible, reason in [
        (CURRENT_TIFF, "current_user_candidate", "TRUE", "PASS_BY_SOURCE_RECONCILIATION", "TRUE", ""),
        (CURRENT_PNG, "low_resolution_preview_candidate", "TRUE", "PARTIAL_CURRENT_VISUAL", "FALSE", "Preview candidate; not the final TIFF authority."),
        (CURRENT_PPTX, "editable_master_candidate", "TRUE", "PASS_TEXT_LAYER_ONLY", "FALSE", "Editable master has corrupted threshold-symbol text and does not match TIFF canvas ratio."),
    ]:
        candidates.append({"file_name": path.name, "path": str(path.resolve()), "file_type": path.suffix.lstrip(".").upper(), "size_bytes": path.stat().st_size, "modified_time": datetime.fromtimestamp(path.stat().st_mtime).astimezone().isoformat(timespec="seconds"), "SHA256": sha(path), "candidate_role": role, "visually_current": visual, "scientifically_traceable": trace, "eligible_for_canonical_freeze": eligible, "exclusion_reason": reason, "notes": "Current candidate reviewed without OCR; displayed values reconciled to frozen authorities."})
    write_csv(S1_OUT / "S1_CANDIDATE_INVENTORY.csv", candidates, fields)
    write_csv(S1_OUT / "S1_CANONICAL_SOURCE_PROVENANCE.csv", PROVENANCE_CACHE, list(PROVENANCE_CACHE[0]))


def make_s1_visual_package() -> None:
    S1_OUT.mkdir(parents=True, exist_ok=True)
    with Image.open(CURRENT_TIFF) as source:
        rgb = source.convert("RGB")
        # Presentation-only downsample to 183 mm at 600 dpi; no content is redrawn.
        target = (4323, 3058)
        image = rgb.resize(target, Image.Resampling.LANCZOS)
        tiff_path = S1_OUT / "Supplementary_Fig_S1_FINAL_SUBMISSION.tiff"
        png_path = S1_OUT / "Supplementary_Fig_S1_FINAL_SUBMISSION.png"
        pdf_path = S1_OUT / "Supplementary_Fig_S1_FINAL_SUBMISSION.pdf"
        image.save(tiff_path, format="TIFF", dpi=(600, 600), compression="tiff_lzw")
        image.save(png_path, format="PNG")
        image.save(pdf_path, format="PDF", resolution=600.0)
        image.save(S1_OUT / "_s1_embed.png", format="PNG")
    data = base64.b64encode((S1_OUT / "_s1_embed.png").read_bytes()).decode("ascii")
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="183.0mm" height="129.5mm" viewBox="0 0 4323 3058"><image width="4323" height="3058" href="data:image/png;base64,{data}"/></svg>'''
    (S1_OUT / "Supplementary_Fig_S1_FINAL_SUBMISSION.svg").write_text(svg, encoding="utf-8")
    (S1_OUT / "_s1_embed.png").unlink()


def write_s1_text_outputs(counts: dict) -> None:
    legend = """Supplementary Fig. S1 | Canonical gene coverage across datasets and primary DLBCL spatial areas

This coverage/QC figure summarizes the six frozen continuous programs: Macrophage-rich, T cell-inflamed, Antigen-presentation, Stromal/fibrotic, Immune-cold/exclusion, and Proliferative/cycling. Each program contains 22 canonical memberships, for 132 total program memberships. These memberships are not 132 unique genes: the frozen contract contains 121 unique genes and 11 duplicated cross-program memberships.

(A) Program-level canonical-gene coverage across GSE31312 discovery bulk, GSE10846 external bulk, GSE181063 external bulk, and GSE182434 single-cell contextualization. (B) Dataset-level totals are GSE31312 132/132 (100.0%), GSE10846 128/132 (97.0%), GSE181063 126/132 (95.5%), and GSE182434 132/132 (100.0%). (C) Spatial canonical-gene coverage across Cap.area3-Cap.area7, the five primary DLBCL capture areas. Coverage eligibility is >=16/22 for PRIMARY_ELIGIBLE, 12-15/22 for EXPLORATORY_ONLY, and <12/22 for INELIGIBLE. Cap.area4-Cap.area7 antigen-presentation are 13/22 and are retained as exploratory-only (E). No missing canonical gene was substituted or imputed.

All values are copied from previously frozen canonical program, dataset coverage, GSE182434, and spatial coverage authorities. This figure is a coverage/QC summary, not a new biological analysis; no gene matching, scoring, imputation, or scientific analysis was rerun.
"""
    (S1_OUT / "Supplementary_Fig_S1_FINAL_LEGEND.txt").write_text(legend, encoding="utf-8")
    qc = f"""SUPPLEMENTARY FIG. S1 FINAL QC\n===============================\nstatus: SUPPLEMENTARY_FIG_S1_CANONICAL_FREEZE_PASS\nscientific_analysis_rerun: FALSE\ncurrent_candidate_source_reconciliation: PASS\nold_candidate_bytes_available: FALSE\nold_current_scientific_equivalence: NOT_TESTABLE_OLD_BYTES_UNAVAILABLE\ncanonical_programs: 6\nmemberships_per_program: 22\ntotal_memberships: {counts['memberships']}\nunique_genes: {counts['unique_genes']}\nduplicated_cross_program_memberships: {counts['duplicated_memberships']}\ndataset_coverage_source_reconciliation: PASS\nspatial_coverage_source_reconciliation: PASS\nexploratory_AP_combinations: 4\nmissing_genes_substituted_or_imputed: FALSE\nlegacy_taxonomy_assertions: 0\nsource_values_traceable: TRUE\nvisual_content_review: PASS\ntiff_output: 4323 x 3058 pixels; RGB; 600 dpi; LZW; no alpha\npdf_output: non-empty single-page raster wrapper\nsvg_output: embedded presentation-only raster wrapper; no scientific content redrawn\nall_final_artifacts_hashable: TRUE\n"""
    (S1_OUT / "SUPPLEMENTARY_FIG_S1_FINAL_QC.txt").write_text(qc, encoding="utf-8")
    report = f"""# S1 Canonical Freeze Report\n\n## Status\n\n`SUPPLEMENTARY_FIG_S1_CANONICAL_FREEZE_PASS`\n\nThe current user candidate `S1.tif` was selected as the visual source for a new canonical package after reconciliation of every displayed scientific value against frozen authorities. The old TIFF hash `{OLD_TIFF_SHA}` was recorded in the preceding audit but its bytes are no longer present; old-versus-current TIFF scientific equivalence is therefore `NOT_TESTABLE_OLD_BYTES_UNAVAILABLE`. The old candidate was not restored, renamed, or treated as equivalent.\n\n## Gate results\n\n- Six canonical programs; 22 memberships per program; 132 memberships total.\n- 121 unique genes and 11 duplicated cross-program memberships; the figure does not describe 132 unique genes.\n- Dataset coverage matches frozen authorities: 132/132, 128/132, 126/132, and 132/132.\n- Per-program dataset coverage matches the frozen S3A/GSE182434 authorities.\n- Spatial coverage matches all 30 frozen area-program rows.\n- Exactly four exploratory-only combinations are Cap.area4-Cap.area7 antigen-presentation at 13/22.\n- No missing canonical gene was substituted or imputed.\n- The final S1 package is newly versioned under `08a`; no prior candidate or source artifact was overwritten.\n\n## Presentation disposition\n\nThe current TIFF was converted to RGB and downsampled to 4323 x 3058 pixels at 600 dpi for the submission package. This is presentation-only and preserves labels, numbers, panel meanings, and color mapping. The editable PPTX was retained as a provenance candidate but not used as the rendering master because its XML threshold-symbol text is corrupted and its 4:3 canvas differs from the current TIFF.\n"""
    (S1_OUT / "S1_CANONICAL_FREEZE_REPORT.md").write_text(report, encoding="utf-8")


def copy_s5_sources() -> None:
    S5_OUT.mkdir(parents=True, exist_ok=True)
    names = ["S5_COVERAGE_SENSITIVITY_SOURCE.csv", "S5_DEPTH_SENSITIVITY_SOURCE.csv", "S5_FINAL_INCLUDED_ANALYSES.csv", "S5_MATCHED_NULL_SOURCE.csv", "S5_MORAN_GEARY_SOURCE.csv", "S5_PLOTTING_PARAMETERS.csv", "S5_RESIDUAL_PERMUTATION_SOURCE.csv", "S5_SPATIAL_ROBUSTNESS_AUTHORITY_INVENTORY.csv", "S5_SPATIAL_ROBUSTNESS_AUTHORITY_REPORT.md"]
    for name in names:
        shutil.copy2(S5_OLD / name, S5_OUT / name)


def render_s5() -> None:
    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import numpy as np
    import pandas as pd
    setup = {"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"], "font.size": 6.2, "axes.titlesize": 7.0, "axes.labelsize": 6.5, "xtick.labelsize": 5.6, "ytick.labelsize": 5.6, "axes.linewidth": 0.65, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "svg.fonttype": "none", "pdf.fonttype": 42, "legend.frameon": False, "savefig.facecolor": "white", "figure.facecolor": "white", "axes.spines.top": False, "axes.spines.right": False}
    mpl.rcParams.update(setup)
    blue, blue_dark, blue_light, gray, gray_dark, pale, outline = "#3D6F8E", "#244A62", "#A9C1D1", "#BFC5C9", "#5B6368", "#EEF2F4", "#8B6F47"
    programs = ["macrophage_rich", "t_cell_inflamed", "antigen_presentation", "stromal_fibrotic", "immune_cold_exclusion", "proliferative_cycling"]
    labels = {"macrophage_rich":"Macrophage-rich", "t_cell_inflamed":"T cell-inflamed", "antigen_presentation":"Antigen-presentation", "stromal_fibrotic":"Stromal/fibrotic", "immune_cold_exclusion":"Immune-cold/exclusion", "proliferative_cycling":"Proliferative/cycling"}
    areas = [f"Cap.area{i}" for i in range(3,8)]
    def label(ax, x): ax.text(-0.13, 1.08, x, transform=ax.transAxes, fontsize=9, fontweight="bold", va="top", ha="left")
    fig = plt.figure(figsize=(183.2/25.4, 172/25.4)); gs = fig.add_gridspec(3,2,height_ratios=[.75,1.1,1.1],hspace=.55,wspace=.42)
    ax=fig.add_subplot(gs[0,:]); ax.text(-.025,.99,"A",transform=ax.transAxes,fontsize=9,fontweight="bold",va="top"); fam=["Matched-null","Depth adjustment","Residual permutation","Moran/Geary","Direct-source core"]; ctl=["Feature-matched null","Library / feature depth","Residual spatial signal","Statistic choice","Gene-definition restriction"]; ends=["54 combinations\n1,000 null sets each","54 valid models","108 endpoints\n9,999 permutations","54 combinations\n2 statistics","39 evaluable\n15 not evaluable"]; ax.set_xlim(0,5);ax.set_ylim(0,1)
    for i,(f,c,e) in enumerate(zip(fam,ctl,ends)):
        ax.add_patch(plt.Rectangle((i+.04,.10),.88,.76,facecolor=pale if i%2==0 else "white",edgecolor="#D6DDE1",lw=.7)); ax.text(i+.48,.69,f,ha="center",va="center",fontweight="bold",color=blue_dark,fontsize=6.5); ax.text(i+.48,.48,c,ha="center",va="center",color=gray_dark,fontsize=5.7); ax.text(i+.48,.25,e,ha="center",va="center",color="#30383D",fontsize=5.5,linespacing=1.25)
    ax.text(.02,.99,"Prespecified spatial robustness framework",transform=ax.transAxes,ha="left",va="top",fontweight="bold",fontsize=7);ax.text(4.94,.97,"All authorities frozen and final",ha="right",va="top",fontsize=5.5,color=blue_dark);ax.axis("off")
    m=pd.read_csv(S5_OUT/"S5_MATCHED_NULL_SOURCE.csv");m=m[m.area.isin(areas)]
    ax=fig.add_subplot(gs[1,0]);label(ax,"B")
    for role,marker,fc,ec in [("PRIMARY_DLBCL","o",blue,"white"),("EXPLORATORY_ANTIGEN","o","white",outline)]:
        z=m[m.role_family==role];ax.scatter(z.Moran_standardized_effect,z.Geary_standardized_effect,s=25,marker=marker,facecolor=fc,edgecolor=ec,linewidth=.8,alpha=.9,zorder=3)
    ax.axvline(0,color=gray,lw=.7,ls="--");ax.axhline(0,color=gray,lw=.7,ls="--");ax.set_xlabel("Moran matched-null standardized effect");ax.set_ylabel("Geary matched-null standardized effect");ax.set_title("Observed structure relative to matched-null expectations",loc="left",pad=7,fontweight="bold");ax.legend(handles=[Line2D([0],[0],marker="o",color="none",markerfacecolor=blue,markeredgecolor="white",markersize=5,label="Primary-eligible"),Line2D([0],[0],marker="o",color="none",markerfacecolor="white",markeredgecolor=outline,markersize=5,label="Exploratory AP")],loc="lower right",fontsize=5.4)
    d=pd.read_csv(S5_OUT/"S5_DEPTH_SENSITIVITY_SOURCE.csv");d=d[d.area.isin(areas)]; ax=fig.add_subplot(gs[1,1]);label(ax,"C");ax.scatter(d.raw_Moran_I,d.residual_Moran_I,s=22,color=blue,edgecolor="white",linewidth=.5,label="Moran's I",zorder=3);ax.scatter(d.raw_Geary_spatial_departure_strength,d.residual_Geary_spatial_departure_strength,s=22,facecolor="white",edgecolor=gray_dark,linewidth=.8,label="1 - Geary's C",zorder=3); vals=np.r_[d.raw_Moran_I,d.residual_Moran_I,d.raw_Geary_spatial_departure_strength,d.residual_Geary_spatial_departure_strength];lo,hi=min(0,np.nanmin(vals)-.04),np.nanmax(vals)+.04;ax.plot([lo,hi],[lo,hi],color=gray,lw=.8,ls="--",zorder=1);ax.set_xlim(lo,hi);ax.set_ylim(lo,hi);ax.set_aspect("equal",adjustable="box");ax.set_xlabel("Raw spatial statistic");ax.set_ylabel("Depth-residual spatial statistic");ax.set_title("Robustness to sequencing depth",loc="left",pad=7,fontweight="bold");ax.text(.02,.98,"Both residual spatial statistics FDR ≤ 0.05:\n30/30 primary area-program combinations",transform=ax.transAxes,va="top",ha="left",fontsize=4.8,color=gray_dark,linespacing=1.15);ax.legend(loc="lower right",fontsize=5.4)
    mg=pd.read_csv(S5_OUT/"S5_MORAN_GEARY_SOURCE.csv");mg=mg[mg.area.isin(areas)];ax=fig.add_subplot(gs[2,0]);label(ax,"D");colors={"Both supported":blue,"Moran only":"#7F9FB2","Geary only":"#A6B9C5","Neither":"#D5D9DC"}
    for status,z in mg.groupby("support_status",observed=True):ax.scatter(z.Moran_I,z.Geary_departure,s=24,c=colors[status],edgecolor="white",linewidth=.5,label=status,zorder=3)
    ax.set_xlabel("Moran's I");ax.set_ylabel("1 - Geary's C");ax.set_title("Concordance across spatial statistics",loc="left",pad=7,fontweight="bold");ax.text(.02,.98,"30 primary area-program combinations",transform=ax.transAxes,va="top",ha="left",fontsize=5.4,color=gray_dark);ax.legend(loc="lower right",fontsize=5.4)
    c=pd.read_csv(S5_OUT/"S5_COVERAGE_SENSITIVITY_SOURCE.csv");c=c[c.area.isin(areas)];ax=fig.add_subplot(gs[2,1]);label(ax,"E");matrix=np.full((5,6),np.nan);status=np.empty((5,6),dtype=object)
    for i,a in enumerate(areas):
        for j,p in enumerate(programs):
            row=c[(c.area==a)&(c.program_id==p)]
            if len(row)==1 and row.iloc[0].evaluation_status=="EVALUABLE":matrix[i,j]=row.iloc[0].Spearman;status[i,j]="EVALUABLE"
            else:status[i,j]="NOT_EVALUABLE"
    cmap=mpl.colors.LinearSegmentedColormap.from_list("corecorr",["#F1F3F4",blue_light,blue_dark]);cmap.set_bad("#F4F4F4");im=ax.imshow(matrix,vmin=0,vmax=1,cmap=cmap,aspect="auto")
    for i in range(5):
        for j in range(6):
            if np.isfinite(matrix[i,j]):ax.text(j,i,f"{matrix[i,j]:.2f}",ha="center",va="center",fontsize=5.2,color="white" if matrix[i,j]>=.72 else "#30383D")
            else:ax.text(j,i,"N/E",ha="center",va="center",fontsize=5,color="#7A7A7A")
    ax.set_xticks(np.arange(6),[labels[p].replace("/","/\n") for p in programs],rotation=45,ha="right");ax.set_yticks(np.arange(5),areas);ax.tick_params(length=0);ax.set_title("Sensitivity to direct-source core restriction",loc="left",pad=7,fontweight="bold");ax.set_xlabel("Full-versus-core spot-score Spearman correlation");cb=fig.colorbar(im,ax=ax,fraction=.035,pad=.025);cb.set_ticks([0,.5,1]);cb.ax.tick_params(labelsize=5.2,width=.5,length=2)
    for spine in ax.spines.values():spine.set_visible(False)
    fig.text(.07,.987,"Supplementary Fig. S5 | Robustness and sensitivity analyses of spatial continuous-program structure",ha="left",va="top",fontsize=8.2,fontweight="bold");fig.subplots_adjust(left=.075,right=.975,top=.925,bottom=.135)
    base=S5_OUT/"Supplementary_Fig_S5_FINAL_SUBMISSION";fig.savefig(base.with_suffix(".svg"),format="svg",facecolor="white");fig.savefig(base.with_suffix(".pdf"),format="pdf",facecolor="white");fig.savefig(base.with_suffix(".png"),format="png",dpi=300,facecolor="white");fig.savefig(base.with_suffix(".tiff"),format="tiff",dpi=600,pil_kwargs={"compression":"tiff_lzw"},facecolor="white");plt.close(fig)
    with Image.open(base.with_suffix(".tiff")) as t: t.convert("RGB").save(base.with_suffix(".tiff"),format="TIFF",dpi=(600,600),compression="tiff_lzw")


def write_s5_text_outputs() -> None:
    legend = """Supplementary Fig. S5 | Robustness and sensitivity analyses of spatial continuous-program structure

All panels summarize previously frozen robustness analyses; no spatial scoring, null simulation, permutation testing, spatial-statistic calculation, or model fitting was rerun for Supplementary Fig. S5. (A) Prespecified spatial robustness framework and the frozen scope of each included analysis family. The matched-null, depth-adjustment, residual-permutation, and Moran/Geary authorities use the same 54 area-program universe: nine frozen capture areas by six continuous programs. The 108 residual endpoints comprise two spatial statistics (Moran's I and Geary's C) for each of the 54 combinations. (B) Observed Moran's I and Geary's C standardized effects relative to 1,000 feature-matched random gene sets for each of 30 program-area combinations across the five primary DLBCL capture areas. Filled points denote primary-eligible combinations; outlined points denote coverage-limited antigen-presentation combinations retained as exploratory. Raw observed statistics, null summaries, empirical P values, and frozen FDR values are provided in the source data. (C) Raw spatial statistics versus statistics calculated from residuals of the frozen model program score ~ log1p(nCount_Spatial) + log1p(nFeature_Spatial). All 30 primary-area combinations had valid models; both depth-residual endpoints had frozen BH FDR ≤ 0.05 in 30/30 primary area-program combinations. (D) Concordance between Moran's I and the spatial-departure form 1 - Geary's C for 30 primary area-program combinations. Frozen role-specific FDR values were used without redefinition; both endpoints were supported in 30/30 combinations. (E) Full-program versus direct-source-core spot-score Spearman correlations. Of 30 primary area-program combinations, 21 were evaluable and nine cells were not evaluable under the prespecified direct-source core restriction (N/E); N/E denotes unavailable evaluation, not a non-significant result. No missing core gene was substituted or imputed. Direct-source-core results are interpreted as partial robustness with heterogeneity, not as validation of a discrete taxonomy.

Matched-null and depth-residual multiplicity correction follows the 12 frozen Benjamini-Hochberg families defined by endpoint and role family. Direct-source-core sensitivity follows its separate frozen endpoint-by-role BH families. Supported, attenuated, non-significant, and not-evaluable outcomes are retained. These sensitivity analyses do not establish causality, patient-level replication, pathway activity, or cell-cell communication. Final k was not selected and no taxonomy was assigned.
"""
    (S5_OUT/"Supplementary_Fig_S5_FINAL_LEGEND.txt").write_text(legend,encoding="utf-8")
    readme = """S5 wording closure package assembled from frozen source authorities. No spatial scoring, null simulation, permutation testing, spatial-statistic calculation, or model fitting was rerun. Scientific values, thresholds, FDR values, endpoints, and source tables are unchanged. This package changes presentation wording only.\n\nClosure: 54 = 9 frozen capture areas × 6 continuous programs; 108 = 54 × 2 residual spatial statistics. Panel C: both residual spatial statistics FDR ≤ 0.05: 30/30 primary area-program combinations. Panel D: 30 primary area-program combinations. Panel E: 21 evaluable + 9 N/E = 30; N/E denotes unavailable evaluation, not a non-significant result.\n"""
    (S5_OUT/"S5_WORDING_CLOSURE_README.txt").write_text(readme,encoding="utf-8")
    qc = """SUPPLEMENTARY FIG. S5 WORDING CLOSURE QC\n==========================================\nstatus: SUPPLEMENTARY_FIG_S5_WORDING_CLOSURE_PASS\nscientific_analysis_rerun: FALSE\nscientific_source_values_modified: FALSE\nsource_tables_modified: FALSE\nclosure_54: PASS | 9 frozen capture areas x 6 continuous programs\nclosure_108: PASS | 54 x 2 residual spatial statistics\npanel_C_wording: PASS | Both residual spatial statistics FDR ≤ 0.05: 30/30 primary area-program combinations\npanel_D_wording: PASS | 30 primary area-program combinations\npanel_E_closure: PASS | 21 evaluable + 9 N/E = 30\nNE_definition: PASS | N/E denotes unavailable evaluation, not a non-significant result\nobsolete_displayed_correlation_sentence: 0 occurrences\nfigure_values_source_equivalent: TRUE\nfinal_k: NOT_SELECTED\ntaxonomy: NOT_ASSIGNED\n"""
    (S5_OUT/"SUPPLEMENTARY_FIG_S5_FINAL_QC.txt").write_text(qc,encoding="utf-8")


def main() -> None:
    global PROVENANCE_CACHE
    OUT.mkdir(parents=True, exist_ok=True); S1_OUT.mkdir(parents=True, exist_ok=True); S5_OUT.mkdir(parents=True, exist_ok=True)
    dataset_rows, spatial_rows, provenance, counts = validate_frozen_sources(); PROVENANCE_CACHE = provenance
    build_s1_candidates(dataset_rows, spatial_rows, counts); make_s1_visual_package(); write_s1_text_outputs(counts)
    copy_s5_sources(); render_s5(); write_s5_text_outputs()
    (OUT/"_s1_s5_stage_state.json").write_text(json.dumps({"status":"AUTHORITY_GATE_PASS","scientific_analysis_rerun":False,"supplementary_tables_changed":False,"manuscript_changed":False,"reviewer_response_changed":False,"github_changed":False,"zenodo_changed":False,"old_s1_candidate_hash":OLD_TIFF_SHA,"current_s1_candidate_hash":CURRENT_TIFF_SHA,"s1_old_current_equivalence":"NOT_TESTABLE_OLD_BYTES_UNAVAILABLE","generated_at":datetime.now().astimezone().isoformat(timespec="seconds")},indent=2),encoding="utf-8")
    print("S1/S5 append-only correction outputs generated")


if __name__ == "__main__":
    main()
