# DLBCL continuous spatially informed microenvironment programs

This repository is the public reproducibility compendium for a continuous-model analysis of the diffuse large B-cell lymphoma (DLBCL) microenvironment.

## Current final model

The model contains six continuous gene-expression programs:

1. Macrophage-rich
2. T cell-inflamed
3. Antigen-presentation
4. Stromal/fibrotic
5. Immune-cold/exclusion
6. Proliferative/cycling

The canonical definition contains 22 memberships per program, 132 memberships in total, 121 unique genes, and 11 duplicated cross-program memberships.

## Model-form result

Candidate values k = 2-6 were evaluated against prespecified retention criteria. No candidate k satisfied all criteria. Therefore, final k = **NOT_SELECTED** and taxonomy = **NOT_ASSIGNED**. The repository does not present a discrete ecosystem classification as the current model.

## Contents

- `code/`: current component workflows and final figure-generation scripts.
- `environment/`: frozen lock, package records, seeds, and session information.
- `source_data/`: frozen machine-readable figure inputs and analysis summaries.
- `figures/`: final main Figures 1-6.
- `supplementary_figures/`: one canonical final package for each of S1-S6.
- `supplementary_tables/`: the publication-ready Supplementary Tables workbook.
- `manifests/`: file hashes and release registries.
- `docs/`: public audits and release notes.

## Public datasets

The study reuses GEO accessions GSE31312, GSE10846, GSE181063, GSE182434, and GSE276542. Raw GEO data are not redistributed. See [DATA_SOURCES.md](DATA_SOURCES.md).

## Reproducibility

See [REPRODUCE.md](REPRODUCE.md). The repository preserves validated component workflows and frozen machine-readable inputs; it does not claim a single-command reconstruction of every historical object from raw data.

## License

Project-authored code and original repository materials are released under the MIT License. This does not relicense third-party data, publications, software, or assets. See [THIRD_PARTY_DATA_AND_SOFTWARE_NOTICE.md](THIRD_PARTY_DATA_AND_SOFTWARE_NOTICE.md).

## Citation

Citation metadata are provided in [CITATION.cff](CITATION.cff). No Zenodo DOI has been assigned or claimed.
