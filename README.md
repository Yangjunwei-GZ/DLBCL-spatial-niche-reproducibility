# DLBCL spatial niche reproducibility

## Overview

This repository supports a continuous-model analysis of the diffuse large B-cell lymphoma (DLBCL) microenvironment. It provides current analysis code, released derived/source data, reproducibility documentation, environment records, and integrity/provenance manifests.

## Programs

1. Macrophage-rich
2. T cell-inflamed
3. Antigen-presentation
4. Stromal/fibrotic
5. Immune-cold/exclusion
6. Proliferative/cycling

These literature-informed curated composite programs each contain 22 recorded memberships: 132 memberships in total, 121 unique genes, and 11 shared cross-program memberships. They are analysis constructs rather than validated biological classes.

## Model-form conclusion

The candidate-k retention rules were prespecified before model-form adjudication. No candidate k from 2-6 was sufficiently stable under the full rule set. A possible k=2 structure remained stability-limited; final k is **NOT_SELECTED**, and no discrete ecosystem taxonomy was assigned.

## Repository structure

- `code/`: current component workflows and final figure-generation scripts.
- `environment/`: the retained lock, package notes, seeds, and session records.
- `source_data/`: machine-readable figure inputs and compact final sensitivity summaries.
- `manifests/`: current-tree hashes, provenance records, and repository quality-control audits.

## Public datasets

The study reuses GEO accessions GSE31312, GSE10846, GSE181063, GSE182434, and GSE276542. Raw GEO data are not redistributed. See [DATA_SOURCES.md](DATA_SOURCES.md).

## Reproduction

See [REPRODUCE.md](REPRODUCE.md). Some figures are reproduced from released derived inputs because a complete raw-to-figure one-command environment is not claimed.

## Figures and source data

Use [source_data/SOURCE_DATA_INDEX.csv](source_data/SOURCE_DATA_INDEX.csv) to map each retained source table to its figure or analysis. Summaries of provenance, candidate-k adjudication, harmonized bulk scoring, exploratory continuous-program survival analyses, orthogonal benchmarking, and spatial sensitivities are under [`source_data/final_sensitivity_summaries/`](source_data/final_sensitivity_summaries/). Large detailed intermediates remain reserved for the associated archival release.

Publication-ready main figures, supplementary figures, and the compiled Supplementary Tables workbook are distributed with the journal submission/publication and are not duplicated in this repository. The retained scripts and machine-readable source authorities document the reported analyses; publication artwork is not required for integrity verification of the released data.

## License

Project-authored code and original repository materials are released under the MIT License. This does not relicense third-party data, publications, software, or assets. See [THIRD_PARTY_DATA_AND_SOFTWARE_NOTICE.md](THIRD_PARTY_DATA_AND_SOFTWARE_NOTICE.md).

## Citation

Citation metadata are provided in [CITATION.cff](CITATION.cff), and the release scope is summarized in [RELEASE_NOTES.md](RELEASE_NOTES.md). The current v1.0.2 Zenodo DOI is [10.5281/zenodo.22008715](https://doi.org/10.5281/zenodo.22008715); the all-versions/latest-version concept DOI is [10.5281/zenodo.21960758](https://doi.org/10.5281/zenodo.21960758).
