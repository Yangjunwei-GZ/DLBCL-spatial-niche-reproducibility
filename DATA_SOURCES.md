# Data sources

Raw and repository-provided processed GEO data retain their original repository terms. This repository does not redistribute raw GEO expression matrices, raw sequencing reads, or source histology images.

| Accession | Study role | Primary publication already verified in the project | Public source | First processing component | Redistribution policy |
|---|---|---|---|---|---|
| GSE31312 | Discovery bulk cohort (498 profiles) | Visco et al. (2012) | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31312 | `code/continuous_model/28_wp1_continuous_score_freeze_v1.R` begins from the frozen score authority; historical raw scoring provenance is documented in `docs/` | Raw/processed GEO source files are not redistributed |
| GSE10846 | External bulk structural validation | Lenz et al. (2008) | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10846 | `code/external_validation/29_wp2_external_continuous_validation_orchestrator_v1.R` | Raw/processed GEO source files are not redistributed |
| GSE181063 | External bulk structural validation | Lacy et al. (2020); Painter et al. (2019) data resource | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181063 | `code/external_validation/29_wp2_external_continuous_validation_orchestrator_v1.R` | Raw/processed GEO source files are not redistributed |
| GSE182434 | Descriptive single-cell contextualization; full dataset 28,416 cells, primary DLBCL subset 14,368 cells from four patients | Steen et al. (2021) | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182434 | `code/single_cell/16_FIGURE5_PHASEC_SCIENTIFIC_EXECUTION.R` | Raw matrices and the full single-cell object are not redistributed |
| GSE276542 | Spatial continuous contextualization and validation across nine capture areas (five primary DLBCL areas) | Diaz-Herrero et al. (2026) | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276542 | `code/spatial/36_wp3_resume_orchestrator_v3.R` | Raw Visium matrices and source images are not redistributed |

Repository-derived source tables are project-generated outputs and are mapped to figures in `manifests/SOURCE_DATA_RELEASE_INVENTORY.csv`. They must not be used to infer patient, specimen, section, or region identities beyond explicit source metadata.
