# Stage 4C-2 Final-k and Model-form Adjudication Report

## Scope and Result

- Adjudication outcome: **DISCRETE TAXONOMY NOT SUPPORTED UNDER FROZEN ACCEPTANCE CRITERIA**
- FINAL K REMAINS NOT SELECTED
- Taxonomy status: NOT_ASSIGNED
- Recommended model form: **CONTINUOUS PROGRAM ABUNDANCE AND IMMUNE-STROMAL-PROLIFERATIVE POLARIZATION SHOULD BE THE PRIMARY MODEL FORM**
- This report is generated from committed Stage 4C-1 and Stage 4C-2A tables. No scientific node was rerun.

## 1. Primary Threshold Failures

- k=2: PAC_rank;minimum_cluster_median_jaccard;overall_median_jaccard (3 failed criteria)
- k=3: PAC_rank;mean_silhouette;minimum_cluster_size;minimum_cluster_median_jaccard;overall_median_jaccard (5 failed criteria)
- k=4: PAC_rank;mean_silhouette;within_0_02_of_best;minimum_cluster_size;minimum_cluster_median_jaccard;overall_median_jaccard (6 failed criteria)
- k=5: mean_silhouette;within_0_02_of_best;minimum_cluster_size;minimum_cluster_median_jaccard;overall_median_jaccard (5 failed criteria)
- k=6: mean_silhouette;within_0_02_of_best;minimum_cluster_size;minimum_cluster_median_jaccard;overall_median_jaccard (5 failed criteria)

Every primary k fails at least three frozen criteria. No k is a single-criterion near miss.

## 2. Nonreference Sensitivity Passes

- k=2: 0/14 nonreference variants pass
- k=3: 0/14 nonreference variants pass
- k=4: 0/14 nonreference variants pass
- k=5: 0/14 nonreference variants pass
- k=6: 0/14 nonreference variants pass

The sensitivity denominator is 14. The primary reference is reported separately and is not described as a sensitivity variant.

## 3. Most Frequent Failure Criteria

Across 70 nonreference variant-k rows: overall_median_jaccard=70/70; minimum_cluster_median_jaccard=66/70; mean_silhouette=65/70; within_0_02_of_best=49/70; PAC_rank=42/70; minimum_cluster_size=30/70.

The detailed per-k frequencies and family-specific structures are in THRESHOLD_FAILURE_FREQUENCY_BY_K.csv.

## 4. Near-miss Assessment

Primary k values failing exactly one criterion: 0. No near-miss k is present under the preregistered descriptive definition.

## 5. Why k=6 Is Not Selected

k=6 has relative agreement rank 1 with median ARI 0.294313 and median NMI 0.446557. However, primary_pass is FALSE and 0/14 nonreference variants pass. Relative rank cannot override absolute frozen acceptance criteria.

## 6. Bootstrap Failure Audit

The 81 FAILED events are concentrated in 2 variants: unique_gene_only_euclidean=62; primary_shrinkage_mahalanobis_whitened=19. By k: 5=33; 6=26; 4=22.
All failures have the message 'Same-k matching requires equal cluster counts.' The highest combination failure rate is 2.200%; every combination retains at least 97.8% completion. Requested/completed/failed/unmatched = 70000/69919/81/0.
Failed replicates were excluded under the frozen rule and were not rerun. Every affected variant-k row also fails at least one non-bootstrap criterion, so these events cannot change any observed variant_pass from FALSE to TRUE.

## 7. Implementation Integrity

Implementation defects found: 0. All 14 audited components pass. The FAILED/UNMATCHED label distinction is retained as implemented and does not alter exclusion or threshold calculations.

## 8. Frozen Criteria

PAC, silhouette, within-best, minimum-size, cluster-Jaccard, and overall-Jaccard criteria were recomputed without changing boundaries. All stored pass fields match the recomputation.

## 9. Discrete-model Support

Any k supported: NO. Under the registered rules, the observed pattern follows Outcome B.

## 10. Final-k Status

FINAL K REMAINS NOT SELECTED

## 11. Taxonomy Status

Taxonomy remains NOT_ASSIGNED. No biological names are assigned in this adjudication.

## 12. Recommended Model Form

CONTINUOUS PROGRAM ABUNDANCE AND IMMUNE-STROMAL-PROLIFERATIVE POLARIZATION SHOULD BE THE PRIMARY MODEL FORM

## 13. Manuscript Treatment of Clustering

Discrete clusters should not be described as validated classes. They may be retained only as exploratory or descriptive supplemental structure, while the six continuous programs and their immune-stromal-proliferative polarization provide the primary biological representation.

## 14. Need for Further Science Rerun

No rerun is required to adjudicate the current frozen discrete model. Future external, spatial, and purity analyses should be reconsidered under the continuous model form and require separate authorization; none is run here.

## Formal Conclusion

**DISCRETE TAXONOMY NOT SUPPORTED UNDER FROZEN ACCEPTANCE CRITERIA**

**FINAL K REMAINS NOT SELECTED**

**CONTINUOUS PROGRAM ABUNDANCE AND IMMUNE-STROMAL-PROLIFERATIVE POLARIZATION SHOULD BE THE PRIMARY MODEL FORM**
