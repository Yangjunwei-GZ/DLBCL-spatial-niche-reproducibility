# WP2 External Continuous Validation Report

## Execution status

WP2 EXTERNAL CONTINUOUS VALIDATION COMPLETED

- GSE10846 samples: 420
- GSE181063 samples: 1310
- Primary gene aggregation: sample-wise median across valid probes per uppercase symbol
- Sensitivity gene aggregation: highest cross-sample MAD probe with lexical probe-ID tiebreak
- Program scores: GSVA ssGSEA, six frozen canonical programs
- Cohorts were not merged or batch corrected

## Canonical coverage

- GSE10846: 1=21/22; 2=20/22; 3=21/22; 4=22/22; 5=22/22; 6=22/22
- GSE181063: 1=20/22; 2=19/22; 3=22/22; 4=22/22; 5=21/22; 6=22/22

No genes were substituted, borrowed, imputed, or restored from a historical shortlist.

## Structural replication

- GSE10846: rho=0.93571, signs=12/15, bootstrap rho 95% CI [0.88929, 0.96786], status=SUPPORTIVE_STRUCTURAL_REPLICATION
- GSE181063: rho=0.98214, signs=13/15, bootstrap rho 95% CI [0.97143, 0.99286], status=SUPPORTIVE_STRUCTURAL_REPLICATION

- Cross-cohort 15-edge Spearman: 0.97143; sign concordance: 14/15

Structural labels follow the preregistered HD6 thresholds and are not selected by P values.

## Continuous PC representations

- Strategy B: completed for both cohorts using cohort-relative z scores and frozen, unflipped WP1 PC1/PC2 loadings.
- GSE10846 Strategy A: BLOCKED.
- GSE181063 Strategy A: not authorized and not attempted.

## Scope controls

- Clinical outcome analysis: NOT RUN
- Spatial analysis: NOT RUN
- Purity analysis: NOT RUN
- Cluster or class assignment: NOT PRODUCED
- final k: NOT_SELECTED
- taxonomy: NOT_ASSIGNED

This WP2 analysis evaluates continuous program-level structural robustness. It does not validate a discrete ecosystem taxonomy or clinical outcome model.
