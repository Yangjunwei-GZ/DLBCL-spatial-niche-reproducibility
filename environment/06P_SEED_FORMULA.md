# 06P deterministic seed reconstruction

The full 54,108-row registry is reserved for the associated Zenodo archival release. The GitHub package retains the deterministic seed contract below.

## Formula

For each capture area order `a` (1-9), program order `p` (1-6), analysis code `c`, and replicate or endpoint index `r`:

```text
derived_seed = 20,000,000 + a*1,000,000 + p*10,000 + c*100 + r
```

Analysis codes:

- `1`: matched-gene-set replicate, `r = 1..1000`
- `2`: depth-residual Moran permutation endpoint, `r = 1`
- `3`: depth-residual Geary permutation endpoint, `r = 1`

Master seed: `20260804`. Capture-area and role-family ordering are defined by the retained current scripts and frozen program order. No random set or scientific value was regenerated during repository slimming.
