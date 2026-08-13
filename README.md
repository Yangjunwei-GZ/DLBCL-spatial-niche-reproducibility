# DLBCL spatial niche reproducibility

## Overview

This repository supports a continuous-model analysis of the diffuse large B-cell lymphoma (DLBCL) microenvironment. It provides current analysis code, released derived source data, display figures, supplementary material, and environment records.

## Programs

1. Macrophage-rich
2. T cell-inflamed
3. Antigen-presentation
4. Stromal/fibrotic
5. Immune-cold/exclusion
6. Proliferative/cycling

Each program contains 22 canonical memberships: 132 memberships in total, 121 unique genes, and 11 cross-program duplicate memberships.

## Model-form conclusion

No candidate k from 2-6 satisfied all prespecified retention criteria. Final k is **NOT_SELECTED**, and no discrete ecosystem taxonomy was assigned.

## Repository structure

- `code/`: current component workflows and final figure-generation scripts.
- `environment/`: the retained lock, package notes, seeds, and session records.
- `source_data/`: machine-readable figure inputs and compact analysis summaries.
- `figures/`: one final submission TIFF and legend for each main Figure 1-6.
- `supplementary_figures/`: one PDF, legend, and source workbook for each of S1-S6.
- `supplementary_tables/`: the final Frontiers submission-ready Supplementary Tables workbook, [`DLBCL_continuous_model_Supplementary_Tables_FRONTIERS_FINAL_SUBMISSION_READY.xlsx`](supplementary_tables/DLBCL_continuous_model_Supplementary_Tables_FRONTIERS_FINAL_SUBMISSION_READY.xlsx).
- `manifests/`: public release hashes and their independent self-check.

## Public datasets

The study reuses GEO accessions GSE31312, GSE10846, GSE181063, GSE182434, and GSE276542. Raw GEO data are not redistributed. See [DATA_SOURCES.md](DATA_SOURCES.md).

## Reproduction

See [REPRODUCE.md](REPRODUCE.md). Some figures are reproduced from released derived authorities because a complete raw-to-figure one-command environment is not claimed.

## Figures and source data

Use [source_data/SOURCE_DATA_INDEX.csv](source_data/SOURCE_DATA_INDEX.csv) to map each retained source table to its figure or analysis. High-resolution submission formats and large detailed intermediates are reserved for the associated Zenodo archival release; no DOI is claimed here.

The released S1-S6 source workbooks and final Supplementary Tables are retained byte-for-byte. Some workbook cells preserve original provenance paths from the controlled analysis environment; these are metadata only and are not runtime dependencies. Public Markdown, CSV registries, and configurable script roots do not contain user-specific local paths.

## License

Project-authored code and original repository materials are released under the MIT License. This does not relicense third-party data, publications, software, or assets. See [THIRD_PARTY_DATA_AND_SOFTWARE_NOTICE.md](THIRD_PARTY_DATA_AND_SOFTWARE_NOTICE.md).

## Citation

Citation metadata are provided in [CITATION.cff](CITATION.cff). No Zenodo DOI has been assigned or claimed.
