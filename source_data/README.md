# Source data

This directory contains frozen machine-readable inputs and result summaries used by the final continuous-program analyses and figures.

- `figure_rendering/` contains portable inputs for main Figures 1-6.
- `single_cell/` contains frozen Figure 5 and GSE182434 contextualization sources.
- `SOURCE_DATA_INDEX.csv` is the public registry for all retained source-data files.

Raw expression matrices, sequencing reads, single-cell objects, raw Visium matrices, and source histology images are not redistributed. Obtain those data from the GEO accessions documented in `../DATA_SOURCES.md`.

The excluded `all_9_spot_signature_scores_long.csv` is a long-format 186,180-row table of the nine-area by six-program primary LogNormalize+UCell scores. It can be regenerated from GSE276542 using the retained `code/spatial/` workflow and its frozen program definitions. Current Figure 6 uses its compact released source authority, and S5 uses its released source workbook. The complete spot-level source table will be included in the associated Zenodo archival release.

Large detailed matched-null, depth, matching-balance, and orthogonal-benchmark tables are likewise reserved for the archival release when their display-ready/current values are already present in the S5 source workbook or final Supplementary Tables.

The final framework is six continuous programs plus PC1/PC2, with final k NOT_SELECTED and taxonomy NOT_ASSIGNED. No scientific values were regenerated or changed while assembling this repository.
