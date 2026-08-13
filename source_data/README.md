# Source data

This directory contains frozen machine-readable inputs and result summaries used by the final continuous-program analyses and figures.

- `figure_rendering/` contains portable inputs for main Figures 1-6.
- `single_cell/` contains frozen Figure 5 and GSE182434 contextualization sources.
- the numbered directories contain spatial scores, matched-null summaries, depth-adjustment outputs, matching-balance summaries, and the orthogonal TME benchmark.

Replicate-level matched-set membership partitions are intentionally not duplicated here; their partition manifest, seed registry, input registry, and final summary tests are retained. Raw expression matrices, sequencing reads, single-cell objects, raw Visium matrices, and source histology images are not redistributed. Obtain those data from the GEO accessions documented in `../DATA_SOURCES.md`.

The final framework is six continuous programs plus PC1/PC2, with final k NOT_SELECTED and taxonomy NOT_ASSIGNED. No scientific values were regenerated or changed while assembling this repository.
