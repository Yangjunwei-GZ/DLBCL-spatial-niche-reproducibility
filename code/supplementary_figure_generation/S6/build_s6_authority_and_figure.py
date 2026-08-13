from __future__ import annotations

import csv
import hashlib
import math
import os
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from PIL import Image
from docx import Document


REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
ROOT = Path(os.environ.get("DLBCL_PROJECT_ROOT", REPOSITORY_ROOT))
OUT = ROOT / "revision_2026_reviewer_response" / "07d_suppfig_s6_lr_expression_support_final"
SRC = ROOT / "results" / "spatial_LR_support"
SCRIPT = ROOT / "revision_2026_reviewer_response" / "source_snapshot" / "scripts" / "09_spatial_LR_coexpression_colocalization_GSE276542.R"
PREFLIGHT = SRC / "Supplementary_Fig_S9_Table_S22_preflight_report.txt"
MANUSCRIPT = ROOT / "revision_2026_reviewer_response" / "06l_wp3_manuscript_insertion" / "DLBCL_manuscript_WP3_integrated_clean.docx"
TRANSITION = ROOT / "revision_2026_reviewer_response" / "06m_global_continuous_model_transition"

AXES = [
    ("MIF-CD74", "MIF", "CD74"),
    ("CXCL9-CXCR3", "CXCL9", "CXCR3"),
    ("CXCL10-CXCR3", "CXCL10", "CXCR3"),
    ("HLA-DRA-CD4", "HLA-DRA", "CD4"),
    ("CD40LG-CD40", "CD40LG", "CD40"),
    ("LGALS9-HAVCR2", "LGALS9", "HAVCR2"),
    ("ICAM1-ITGAL", "ICAM1", "ITGAL"),
    ("IFNG-IFNGR1", "IFNG", "IFNGR1"),
    ("CD70-CD27", "CD70", "CD27"),
    ("CD86-CD28", "CD86", "CD28"),
]
AREAS = [
    "GSM8500536_Cap.area3_DLBCL_V1",
    "GSM8500537_Cap.area4_DLBCL_V2",
    "GSM8500538_Cap.area5_DLBCL_V2",
    "GSM8500539_Cap.area6_DLBCL_V2",
    "GSM8500540_Cap.area7_DLBCL_V2",
]
AREA_LABELS = ["Cap.area3", "Cap.area4", "Cap.area5", "Cap.area6", "Cap.area7"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def clean_area(value: str) -> str:
    match = re.search(r"Cap\.area(\d+)", str(value))
    return f"Cap.area{match.group(1)}" if match else str(value)


def load_sources() -> dict[str, pd.DataFrame]:
    files = {
        "detection": "S22B_gene_detection.csv",
        "correlation": "S22C_spot_level_correlation.csv",
        "neighborhood": "S22D_neighborhood_enrichment.csv",
        "old_context": "S22E_niche_context.csv",
        "maps": "Supplementary_Fig_S9D_spatial_map_source.csv",
        "areas": "GSE276542_spatial_LR_capture_area_summary.csv",
        "missing": "GSE276542_spatial_LR_missing_genes_by_capture_area.csv",
    }
    out = {}
    for key, name in files.items():
        path = SRC / name
        if not path.exists() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
        out[key] = pd.read_csv(path)
    return out


def validate_sources(data: dict[str, pd.DataFrame]) -> None:
    expected_pairs = [x[0] for x in AXES]
    for key in ("detection", "correlation", "neighborhood"):
        df = data[key]
        if len(df) != 50:
            raise ValueError(f"{key}: expected 50 rows, found {len(df)}")
        if sorted(df["pair_id"].unique()) != sorted(expected_pairs):
            raise ValueError(f"{key}: LR-axis mismatch")
        if sorted(df["capture_area"].unique()) != sorted(AREAS):
            raise ValueError(f"{key}: capture-area mismatch")
        if df.duplicated(["capture_area", "pair_id"]).any():
            raise ValueError(f"{key}: duplicate area-axis rows")
    det = data["detection"]
    if int((det["status"] == "tested").sum()) != 45:
        raise ValueError("Detection tested-row count is not 45")
    if int((det["status"] != "tested").sum()) != 5:
        raise ValueError("Detection not-tested-row count is not 5")
    expected_missing = {
        (AREAS[0], "LGALS9-HAVCR2", "LGALS9"),
        (AREAS[1], "HLA-DRA-CD4", "HLA-DRA"),
        (AREAS[2], "HLA-DRA-CD4", "HLA-DRA"),
        (AREAS[3], "HLA-DRA-CD4", "HLA-DRA"),
        (AREAS[4], "HLA-DRA-CD4", "HLA-DRA"),
    }
    observed_missing = set(
        det.loc[det["status"] != "tested", ["capture_area", "pair_id", "missing_gene"]]
        .itertuples(index=False, name=None)
    )
    if observed_missing != expected_missing:
        raise ValueError(f"Unexpected missing-gene rows: {observed_missing}")
    for key in ("correlation", "neighborhood"):
        df = data[key]
        tested = df["status"] == "tested"
        if int(tested.sum()) != 45 or int(df.loc[tested, "BH_FDR"].notna().sum()) != 45:
            raise ValueError(f"{key}: tested/FDR completeness failure")
    maps = data["maps"]
    if len(maps) != 3446 or set(maps["pair_id"]) != {"MIF-CD74", "CXCL9-CXCR3"}:
        raise ValueError("Representative map source is incomplete")
    if maps.groupby("pair_id")["barcode"].nunique().to_dict() != {
        "CXCL9-CXCR3": 1723,
        "MIF-CD74": 1723,
    }:
        raise ValueError("Representative map barcodes are incomplete")
    for col in ("x_coord", "y_coord", "ligand_expression", "receptor_expression"):
        if not np.isfinite(pd.to_numeric(maps[col], errors="coerce")).all():
            raise ValueError(f"Map source has non-finite {col}")


def export_frozen_source_tables(data: dict[str, pd.DataFrame]) -> None:
    names = {
        "detection": "S6_SOURCE_GENE_DETECTION.csv",
        "correlation": "S6_SOURCE_SAME_SPOT_CORRELATION.csv",
        "neighborhood": "S6_SOURCE_NEIGHBORHOOD_ENRICHMENT.csv",
        "maps": "S6_SOURCE_REPRESENTATIVE_GENE_EXPRESSION_MAPS.csv",
        "areas": "S6_SOURCE_CAPTURE_AREA_SUMMARY.csv",
    }
    for key, name in names.items():
        data[key].to_csv(OUT / name, index=False, encoding="utf-8", lineterminator="\n")


def build_historical_audit(data: dict[str, pd.DataFrame]) -> None:
    rows = []
    definitions = {
        "S9B": ("same-spot correlation", "Spearman rho; BH across 45 tested area-axis rows", True, ""),
        "S9C": ("neighborhood support", "k=6; 1,000 permutations; empirical P; BH across 45 tested area-axis rows", True, ""),
    }
    for panel, (program, method, eligible, reason) in definitions.items():
        frame = data["correlation"] if panel == "S9B" else data["neighborhood"]
        source = SRC / ("S22C_spot_level_correlation.csv" if panel == "S9B" else "S22D_neighborhood_enrichment.csv")
        for row in frame.itertuples(index=False):
            rows.append({
                "Historical_file": source.name,
                "Historical_panel": panel,
                "LR_axis": row.pair_id,
                "Ligand": row.ligand,
                "Receptor": row.receptor,
                "Spatial_area": clean_area(row.capture_area),
                "Program_or_score": program,
                "Score_method": "not applicable; gene expression support",
                "Statistical_method": method,
                "Current_framework_compatible": "TRUE",
                "Eligible_for_final_S6": "TRUE" if eligible else "FALSE",
                "Exclusion_reason": reason,
                "Source_path": str(source.resolve()),
                "SHA256": sha256(source),
                "Notes": "N/E retained when a prespecified gene was unavailable; no imputation.",
            })
    context_source = SRC / "S22E_niche_context.csv"
    for row in data["old_context"].itertuples(index=False):
        rows.append({
            "Historical_file": context_source.name,
            "Historical_panel": "S9 historical program-association component",
            "LR_axis": row.pair_id,
            "Ligand": row.ligand,
            "Receptor": row.receptor,
            "Spatial_area": clean_area(row.capture_area),
            "Program_or_score": row.niche_score,
            "Score_method": "historical niche score; not demonstrated as canonical UCell authority",
            "Statistical_method": "Spearman correlation of historical LR score with old niche score; BH FDR",
            "Current_framework_compatible": "FALSE",
            "Eligible_for_final_S6": "FALSE",
            "Exclusion_reason": "Old score framework is incompatible or uncertain under the final canonical UCell model.",
            "Source_path": str(context_source.resolve()),
            "SHA256": sha256(context_source),
            "Notes": "Excluded without recomputation.",
        })
    map_source = SRC / "Supplementary_Fig_S9D_spatial_map_source.csv"
    for axis, ligand, receptor in AXES[:2]:
        rows.append({
            "Historical_file": map_source.name,
            "Historical_panel": "S9D",
            "LR_axis": axis,
            "Ligand": ligand,
            "Receptor": receptor,
            "Spatial_area": "Cap.area4",
            "Program_or_score": "separate ligand and receptor expression",
            "Score_method": "frozen normalized spatial expression",
            "Statistical_method": "descriptive map; no new test",
            "Current_framework_compatible": "TRUE",
            "Eligible_for_final_S6": "TRUE",
            "Exclusion_reason": "",
            "Source_path": str(map_source.resolve()),
            "SHA256": sha256(map_source),
            "Notes": "Final S6 plots ligand and receptor separately and does not rely on the composite historical LR score.",
        })
    write_csv(
        OUT / "S6_HISTORICAL_S9_AUTHORITY_AUDIT.csv",
        rows,
        [
            "Historical_file", "Historical_panel", "LR_axis", "Ligand", "Receptor",
            "Spatial_area", "Program_or_score", "Score_method", "Statistical_method",
            "Current_framework_compatible", "Eligible_for_final_S6", "Exclusion_reason",
            "Source_path", "SHA256", "Notes",
        ],
    )


def build_inventory(data: dict[str, pd.DataFrame]) -> None:
    rows = []
    source_specs = [
        ("detection", "S22B_gene_detection.csv", "gene detection", "detection fractions", None, None),
        ("correlation", "S22C_spot_level_correlation.csv", "same-spot expression concordance", "Spearman rho", "p_value", "BH_FDR"),
        ("neighborhood", "S22D_neighborhood_enrichment.csv", "kNN neighborhood expression support", "log2 enrichment; k=6; 1,000 permutations", "empirical_p_value", "BH_FDR"),
    ]
    for key, name, analysis, statistic, pcol, fcol in source_specs:
        source = SRC / name
        for row in data[key].to_dict("records"):
            rows.append({
                "LR_axis": row["pair_id"],
                "Ligand": row["ligand"],
                "Receptor": row["receptor"],
                "Source_file": str(source.resolve()),
                "Source_analysis": analysis,
                "Spatial_scope": clean_area(row["capture_area"]),
                "Program_context": "none; expression-only final S6",
                "Program_score_method": "not used",
                "Statistic": statistic + (f"={row.get('spearman_rho')}" if key == "correlation" and pd.notna(row.get("spearman_rho")) else "") + (f"={row.get('log2_enrichment_ratio')}" if key == "neighborhood" and pd.notna(row.get("log2_enrichment_ratio")) else ""),
                "P_value": row.get(pcol, "") if pcol else "",
                "FDR": row.get(fcol, "") if fcol else "",
                "Frozen_status": "PASS preflight; tested" if row["status"] == "tested" else "PASS preflight; not evaluable due to missing gene",
                "Current_framework_compatible": "TRUE",
                "Eligible_for_S6": "TRUE",
                "Notes": f"status={row['status']}; missing_gene={row.get('missing_gene', '')}",
                "SHA256": sha256(source),
            })
    map_source = SRC / "Supplementary_Fig_S9D_spatial_map_source.csv"
    for axis, ligand, receptor in AXES[:2]:
        rows.append({
            "LR_axis": axis,
            "Ligand": ligand,
            "Receptor": receptor,
            "Source_file": str(map_source.resolve()),
            "Source_analysis": "representative separate-gene expression maps",
            "Spatial_scope": "Cap.area4",
            "Program_context": "none; expression-only final S6",
            "Program_score_method": "not used",
            "Statistic": "descriptive normalized ligand and receptor expression",
            "P_value": "",
            "FDR": "",
            "Frozen_status": "PASS preflight",
            "Current_framework_compatible": "TRUE",
            "Eligible_for_S6": "TRUE",
            "Notes": "Historical composite LR score excluded from final map display.",
            "SHA256": sha256(map_source),
        })
    write_csv(
        OUT / "S6_FINAL_LR_AUTHORITY_INVENTORY.csv",
        rows,
        [
            "LR_axis", "Ligand", "Receptor", "Source_file", "Source_analysis",
            "Spatial_scope", "Program_context", "Program_score_method", "Statistic",
            "P_value", "FDR", "Frozen_status", "Current_framework_compatible",
            "Eligible_for_S6", "Notes", "SHA256",
        ],
    )


def build_exclusion_log() -> None:
    rows = [
        {
            "Historical_panel": "S9 historical niche-context panel / S22E",
            "Historical_analysis": "LR score versus six historical niche scores",
            "Reason_excluded": "Program-score provenance is historical and not demonstrated as final canonical UCell authority.",
            "Old_framework": "niche score / historical program columns",
            "Current_replacement": "No program-association panel; expression-only S6.",
            "Scientific_result_rerun": "NO",
            "Notes": "All 300 rows retained in the historical audit but excluded from final scientific outputs.",
        },
        {
            "Historical_panel": "Historical S9 workflow endpoint",
            "Historical_analysis": "LR score versus niche scores",
            "Reason_excluded": "Workflow endpoint depends on the excluded historical score framework.",
            "Old_framework": "niche score",
            "Current_replacement": "Workflow is replaced by an auditable scope table without program associations.",
            "Scientific_result_rerun": "NO",
            "Notes": "Detection, same-spot correlation, and neighborhood support remain eligible.",
        },
        {
            "Historical_panel": "Historical S9D composite maps",
            "Historical_analysis": "sqrt(ligand expression x receptor expression) LR score",
            "Reason_excluded": "The composite score is not needed to show frozen gene-expression support and could be overinterpreted.",
            "Old_framework": "derived LR score",
            "Current_replacement": "Separate ligand and receptor expression maps from the same frozen source rows.",
            "Scientific_result_rerun": "NO",
            "Notes": "No new statistic; only frozen columns are visualized.",
        },
        {
            "Historical_panel": "Historical Figure 6C-D mechanism panels",
            "Historical_analysis": "Class-prioritized LR axes and expression dot plot",
            "Reason_excluded": "Class-dependent mechanism claim was retired when no discrete taxonomy was retained.",
            "Old_framework": "four-class ecosystem / class contrast",
            "Current_replacement": "Supplementary expression-only spatial support, without class labels or mechanism inference.",
            "Scientific_result_rerun": "NO",
            "Notes": "Historical axis provenance is retained; class-dependent conclusions are not.",
        },
    ]
    write_csv(
        OUT / "S6_HISTORICAL_RESULT_EXCLUSION_LOG.csv",
        rows,
        ["Historical_panel", "Historical_analysis", "Reason_excluded", "Old_framework", "Current_replacement", "Scientific_result_rerun", "Notes"],
    )


def build_compatibility_report() -> None:
    text = f"""# S6 Current Framework Compatibility Report

## Controlling status

`CURRENT_UCELL_COMPATIBLE`

Final disposition: `SUPPLEMENTARY_FIG_S6_EXPRESSION_ONLY_PASS`

This compatibility statement is qualified: the final S6 contains no program-association panel and therefore does not use historical niche scores or any replacement UCell-program association. The frozen expression-support components are independent of the retired score framework and are compatible with the current manuscript when interpreted descriptively.

## Authority gate

- The historical script explicitly registers 10 curated Figure 6 ligand-receptor axes.
- The frozen preflight status is PASS.
- The final primary spatial scope is the five DLBCL capture areas Cap.area3-Cap.area7.
- S22B contains 50 area-axis detection rows: 45 tested and 5 not evaluable because a prespecified gene was unavailable.
- S22C contains the frozen same-spot Spearman results and BH FDR values for all 45 tested combinations.
- S22D contains the frozen k=6, 1,000-permutation neighborhood results and BH FDR values for all 45 tested combinations.
- The representative map source contains 1,723 spots for each of MIF-CD74 and CXCL9-CXCR3 in Cap.area4.
- S22E uses historical niche-score columns and is excluded in full.
- No current canonical UCell LR-program association is asserted or created.

## Missingness

- Cap.area3 LGALS9-HAVCR2: LGALS9 unavailable.
- Cap.area4-Cap.area7 HLA-DRA-CD4: HLA-DRA unavailable.
- These five combinations remain N/E. No zero replacement, imputation, synonym substitution, or alternative gene was used.

## Interpretation boundary

Expression-based spatial support does not establish direct ligand-receptor signaling, physical cell-cell contact, receptor activation, functional communication, or causality.

## Provenance

- Historical generating script: `{SCRIPT.resolve()}`
- Frozen preflight: `{PREFLIGHT.resolve()}`
- Frozen result directory: `{SRC.resolve()}`
- Script SHA-256: `{sha256(SCRIPT)}`
- Preflight SHA-256: `{sha256(PREFLIGHT)}`

No LR discovery, correlation test, neighborhood test, permutation, FDR calculation, UCell scoring, pathway analysis, differential expression, clustering, or taxonomy analysis was run for S6.
"""
    (OUT / "S6_CURRENT_FRAMEWORK_COMPATIBILITY_REPORT.md").write_text(text, encoding="utf-8")


def collect_wording_hits() -> list[dict]:
    hits = []
    terms = re.compile(r"ligand|receptor|communication|signaling|interaction|activation|mediates", re.I)
    if MANUSCRIPT.exists():
        doc = Document(MANUSCRIPT)
        for idx, paragraph in enumerate(doc.paragraphs, start=1):
            text = paragraph.text.strip()
            if text and terms.search(text):
                evidence = "Expression-only spatial support" if re.search(r"ligand|receptor", text, re.I) else "No direct LR evidence"
                strong = bool(re.search(r"active signaling|confirmed communication|receptor activation|mediates|demonstrat(?:e|es|ed) communication", text, re.I))
                hits.append({
                    "location": f"{MANUSCRIPT.resolve()} paragraph {idx}",
                    "existing": text,
                    "evidence": evidence,
                    "recommendation": "Restrict to expression-based spatial support; do not infer communication, contact, activation, or causality." if strong else "Retain boundary language; verify supplement citation at manuscript integration.",
                    "change": "YES" if strong else "NO",
                })
    inventory = TRANSITION / "MANUSCRIPT_DISCRETE_MODEL_CLAIM_INVENTORY.csv"
    if inventory.exists():
        frame = pd.read_csv(inventory)
        for row in frame.to_dict("records"):
            claim = str(row.get("claim_text", row.get("Claim_text", "")))
            if not claim:
                for value in row.values():
                    if isinstance(value, str) and terms.search(value):
                        claim = value
                        break
            if claim and terms.search(claim):
                hits.append({
                    "location": f"{inventory.resolve()} record {row.get('claim_id', row.get('Claim_ID', 'historical'))}",
                    "existing": claim,
                    "evidence": "Historical claim inventory; class-dependent LR inference is retired.",
                    "recommendation": "Do not restore class-dependent LR or communication wording. Cite S6 only as expression support.",
                    "change": "YES" if re.search(r"ecosystem|communication|signaling|interaction", claim, re.I) else "NO",
                })
    return hits


def build_wording_audit() -> None:
    hits = collect_wording_hits()
    lines = [
        "# S6 Main-text LR Wording Audit",
        "",
        "This is a read-only wording audit. No manuscript file was modified.",
        "",
        "## Required interpretation",
        "",
        "Expression-based spatial support does not establish direct ligand-receptor signaling, physical cell-cell contact, receptor activation, functional communication, or causality.",
        "",
        "## Findings",
        "",
        "| Location | Existing wording | Evidence actually available | Recommended restrained wording | Change required |",
        "|---|---|---|---|---|",
    ]
    for hit in hits:
        esc = lambda value: str(value).replace("|", "\\|").replace("\n", " ")
        lines.append(f"| {esc(hit['location'])} | {esc(hit['existing'])} | {esc(hit['evidence'])} | {esc(hit['recommendation'])} | {hit['change']} |")
    lines.extend([
        "",
        "## S6 insertion recommendation",
        "",
        "Use: `spatial expression support`, `same-spot expression concordance`, and `neighborhood expression support`.",
        "",
        "Do not use: `active signaling`, `confirmed communication`, `receptor activation`, `physical contact`, `mediates`, or causal wording.",
        "",
        "The current final Figure 6 legend already states that score-level spatial associations do not establish direct cellular contact, intercellular communication, or causality. S6 should preserve the same boundary.",
    ])
    (OUT / "S6_MAIN_TEXT_LR_WORDING_AUDIT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def make_figure(data: dict[str, pd.DataFrame]) -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 6.5,
        "axes.titlesize": 7.2,
        "axes.labelsize": 6.5,
        "xtick.labelsize": 5.8,
        "ytick.labelsize": 5.8,
        "axes.linewidth": 0.6,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
    })
    fig = plt.figure(figsize=(183 / 25.4, 216 / 25.4), constrained_layout=False)
    outer = GridSpec(3, 1, figure=fig, height_ratios=[1.05, 2.15, 2.45], hspace=0.46, top=0.965, bottom=0.055, left=0.095, right=0.945)
    ax_a = fig.add_subplot(outer[0])
    mid = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[1], width_ratios=[1, 1], wspace=0.42)
    ax_b = fig.add_subplot(mid[0])
    ax_c = fig.add_subplot(mid[1])
    bot = GridSpecFromSubplotSpec(1, 4, subplot_spec=outer[2], wspace=0.12)
    map_axes = [fig.add_subplot(bot[i]) for i in range(4)]

    fig.suptitle("Supplementary Fig. S6 | Spatial expression support for curated ligand-receptor axes", x=0.095, y=0.992, ha="left", fontsize=9, fontweight="bold")

    # Panel A: compact authority scope table.
    ax_a.set_axis_off()
    ax_a.text(-0.035, 1.04, "a", transform=ax_a.transAxes, fontsize=9, fontweight="bold", va="top")
    ax_a.text(0.0, 1.04, "Curated axes and frozen primary-area scope", transform=ax_a.transAxes, fontsize=7.5, fontweight="bold", va="top")
    headers = ["Axis", "Ligand", "Receptor", "Tested", "N/E"]
    offsets = [0.00, 0.18, 0.30, 0.40, 0.475]
    for block in (0.01, 0.51):
        for offset, header in zip(offsets, headers):
            ax_a.text(block + offset, 0.84, header, transform=ax_a.transAxes, fontsize=5.9, fontweight="bold", color="#263238")
    for i, (axis, ligand, receptor) in enumerate(AXES):
        col = 0 if i < 5 else 1
        row = i if i < 5 else i - 5
        x0 = 0.01 + col * 0.50
        y = 0.68 - row * 0.135
        subset = data["detection"].query("pair_id == @axis")
        tested = int((subset["status"] == "tested").sum())
        ne = 5 - tested
        ax_a.text(x0, y, axis, transform=ax_a.transAxes, fontsize=5.9)
        ax_a.text(x0 + 0.18, y, ligand, transform=ax_a.transAxes, fontsize=5.9)
        ax_a.text(x0 + 0.30, y, receptor, transform=ax_a.transAxes, fontsize=5.9)
        ax_a.text(x0 + 0.40, y, f"{tested}/5", transform=ax_a.transAxes, fontsize=5.9)
        ax_a.text(x0 + 0.475, y, str(ne), transform=ax_a.transAxes, fontsize=5.9, color="#5F6B73")
    ax_a.text(0.01, -0.04, "Current canonical UCell LR-program association: not available; no association panel was generated.", transform=ax_a.transAxes, fontsize=5.7, color="#4D5A61")
    ax_a.add_patch(Rectangle((0, -0.10), 1, 1.01, transform=ax_a.transAxes, facecolor="none", edgecolor="#C7D0D5", linewidth=0.6))

    # Panel B: frozen same-spot correlations.
    cor = data["correlation"].copy()
    pivot = cor.pivot(index="pair_id", columns="capture_area", values="spearman_rho").reindex(index=[x[0] for x in AXES], columns=AREAS)
    fdr = cor.pivot(index="pair_id", columns="capture_area", values="BH_FDR").reindex(index=pivot.index, columns=pivot.columns)
    status = cor.pivot(index="pair_id", columns="capture_area", values="status").reindex(index=pivot.index, columns=pivot.columns)
    div_cmap = LinearSegmentedColormap.from_list("s6div", ["#315A7D", "#F5F6F6", "#C45A3C"])
    vmax = max(abs(np.nanmin(pivot.values)), abs(np.nanmax(pivot.values)))
    im = ax_b.imshow(pivot.values, cmap=div_cmap, norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax), aspect="auto")
    ax_b.set_xticks(range(5), AREA_LABELS, rotation=45, ha="right")
    ax_b.set_yticks(range(10), [x[0] for x in AXES])
    ax_b.set_title("b  Same-spot ligand-receptor expression concordance", loc="left", fontweight="bold", pad=7)
    for i in range(10):
        for j in range(5):
            if status.iloc[i, j] != "tested":
                ax_b.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, facecolor="#D9DEE1", edgecolor="white", linewidth=0.8))
                ax_b.text(j, i, "N/E", ha="center", va="center", fontsize=5.5, color="#4D565C")
            elif fdr.iloc[i, j] < 0.05:
                ax_b.text(j, i, "*", ha="center", va="center", fontsize=7.2, fontweight="bold", color="#111111")
    ax_b.set_xticks(np.arange(-0.5, 5, 1), minor=True)
    ax_b.set_yticks(np.arange(-0.5, 10, 1), minor=True)
    ax_b.grid(which="minor", color="white", linewidth=0.8)
    ax_b.tick_params(which="minor", bottom=False, left=False)
    for spine in ax_b.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(im, ax=ax_b, fraction=0.042, pad=0.03)
    cbar.set_label("Spearman rho", fontsize=5.8)
    cbar.ax.tick_params(labelsize=5.3, length=2)

    # Panel C: frozen neighborhood enrichment.
    neigh = data["neighborhood"].copy()
    x_lookup = {area: i for i, area in enumerate(AREAS)}
    y_lookup = {axis[0]: i for i, axis in enumerate(AXES)}
    tested = neigh[neigh["status"] == "tested"].copy()
    xs = tested["capture_area"].map(x_lookup).to_numpy()
    ys = tested["pair_id"].map(y_lookup).to_numpy()
    colors = tested["log2_enrichment_ratio"].to_numpy(float)
    neglog = -np.log10(np.clip(tested["BH_FDR"].to_numpy(float), 1e-300, 1))
    sizes = 7 + 18 * np.clip(neglog, 0, 2.5) / 2.5
    cmax = max(abs(np.nanmin(colors)), abs(np.nanmax(colors)))
    sc = ax_c.scatter(xs, ys, c=colors, s=sizes, cmap=div_cmap, norm=TwoSlopeNorm(vmin=-cmax, vcenter=0, vmax=cmax), edgecolor="#FFFFFF", linewidth=0.25)
    for row in neigh[neigh["status"] != "tested"].itertuples(index=False):
        x = x_lookup[row.capture_area]
        y = y_lookup[row.pair_id]
        ax_c.add_patch(Rectangle((x - 0.38, y - 0.38), 0.76, 0.76, facecolor="#D9DEE1", edgecolor="none"))
        ax_c.text(x, y, "N/E", ha="center", va="center", fontsize=5.1, color="#4D565C")
    for row in tested[tested["BH_FDR"] < 0.05].itertuples(index=False):
        ax_c.text(x_lookup[row.capture_area], y_lookup[row.pair_id], "*", ha="center", va="center", fontsize=6.8, fontweight="bold", color="#111111")
    ax_c.set_xlim(-0.6, 4.6)
    ax_c.set_ylim(9.6, -0.6)
    ax_c.set_xticks(range(5), AREA_LABELS, rotation=45, ha="right")
    ax_c.set_yticks(range(10), [x[0] for x in AXES])
    ax_c.set_title("c  Frozen kNN neighborhood expression support", loc="left", fontweight="bold", pad=7)
    ax_c.grid(color="#E8ECEE", linewidth=0.45, zorder=0)
    for spine in ax_c.spines.values():
        spine.set_visible(False)
    cbar2 = fig.colorbar(sc, ax=ax_c, fraction=0.042, pad=0.03)
    cbar2.set_label("log2 enrichment", fontsize=5.8)
    cbar2.ax.tick_params(labelsize=5.3, length=2)
    ax_c.text(0.98, -0.20, "Dot size: -log10 BH FDR | * BH FDR < 0.05", transform=ax_c.transAxes, fontsize=5.1, va="top", ha="right")

    # Panel D: separate ligand and receptor maps from frozen columns.
    maps = data["maps"].copy()
    specs = [
        ("MIF-CD74", "MIF", "ligand_expression"),
        ("MIF-CD74", "CD74", "receptor_expression"),
        ("CXCL9-CXCR3", "CXCL9", "ligand_expression"),
        ("CXCL9-CXCR3", "CXCR3", "receptor_expression"),
    ]
    seq_cmap = LinearSegmentedColormap.from_list("s6seq", ["#EDF3F4", "#73A9A7", "#175B63"])
    for idx, (axis, gene, value_col) in enumerate(specs):
        ax = map_axes[idx]
        sub = maps[maps["pair_id"] == axis].copy()
        values = sub[value_col].to_numpy(float)
        vmax_gene = float(np.nanpercentile(values, 99))
        artist = ax.scatter(sub["x_coord"], -sub["y_coord"], c=np.clip(values, 0, vmax_gene), s=2.0, cmap=seq_cmap, vmin=0, vmax=vmax_gene, linewidths=0, rasterized=True)
        ax.set_title(f"{gene}\n{axis}", fontsize=6.2, fontweight="bold", pad=2, y=0.98)
        ax.set_aspect("equal")
        ax.set_axis_off()
        cb = fig.colorbar(artist, ax=ax, orientation="horizontal", fraction=0.04, pad=0.035, aspect=18)
        cb.ax.tick_params(labelsize=4.8, length=1.5, pad=1)
        cb.set_label("Frozen normalized expression", fontsize=5.0, labelpad=1)
    fig.text(0.055, 0.348, "d", fontsize=9, fontweight="bold", va="top")
    fig.text(0.078, 0.348, "Representative separate-gene expression maps | Cap.area4", fontsize=7.5, fontweight="bold", va="top")
    fig.text(0.095, 0.018, "N/E: prespecified metric not evaluable from frozen source data; not a zero or nonsignificant result.\nExpression support does not establish signaling, contact, activation, communication, or causality.", fontsize=5.3, color="#3F4A50", va="bottom")

    base = OUT / "Supplementary_Fig_S6_FINAL_SUBMISSION"
    fig.savefig(base.with_suffix(".pdf"), bbox_inches=None)
    fig.savefig(base.with_suffix(".svg"), bbox_inches=None)
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches=None)
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches=None, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    with Image.open(base.with_suffix(".tiff")) as image:
        rgb = image.convert("RGB")
        if rgb.mode != image.mode:
            rgb.save(base.with_suffix(".tiff"), dpi=(600, 600), compression="tiff_lzw")


def write_legend() -> None:
    legend = """Supplementary Fig. S6 | Spatial expression support for curated ligand-receptor axes. a, Audit scope for 10 pre-existing curated axes across the five primary DLBCL capture areas in GSE276542. Areas tested reports the number of capture areas in which both prespecified genes were available. No current canonical UCell LR-program association was available or generated. b, Frozen same-spot Spearman correlations between ligand and receptor expression across spots. The Benjamini-Hochberg FDR values were retained from the original family of 45 tested capture-area-by-axis combinations; asterisks denote BH FDR < 0.05. Positive, negative, and nonsignificant results are all shown. c, Frozen k-nearest-neighbor neighborhood expression support using k=6 and 1,000 permutations. Color shows the frozen log2 enrichment ratio, dot size shows -log10 BH FDR, and asterisks denote BH FDR < 0.05 across the original 45 tested combinations. d, Separate frozen ligand and receptor expression maps for MIF-CD74 and CXCL9-CXCR3 in the prespecified representative primary DLBCL area GSM8500537_Cap.area4_DLBCL_V2. Maps use the frozen normalized gene-expression columns and do not display or recompute a composite LR score. N/E indicates that the prespecified expression-support metric could not be evaluated from the frozen source data because a ligand or receptor gene was unavailable; it does not indicate absence of signaling or a non-significant result. Same-spot expression compatibility is not evidence of a same-cell interaction, and neighborhood support reflects spatial proximity rather than physical contact. These analyses provide expression-based spatial support for the curated ligand-receptor axes and do not establish direct ligand-receptor signaling, physical cell-cell contact, receptor activation, functional communication, or causality.
"""
    (OUT / "Supplementary_Fig_S6_FINAL_LEGEND.txt").write_text(legend, encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    data = load_sources()
    validate_sources(data)
    export_frozen_source_tables(data)
    build_historical_audit(data)
    build_inventory(data)
    build_exclusion_log()
    build_compatibility_report()
    build_wording_audit()
    make_figure(data)
    write_legend()
    print("S6 authority gate and Python figure generation complete")


if __name__ == "__main__":
    main()
