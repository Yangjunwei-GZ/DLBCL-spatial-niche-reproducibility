# Final-k and Model-form Adjudication Rules

## Scope

This adjudication is read-only. It uses the committed Stage 4C-1 primary results and committed Stage 4C-2A sensitivity results. It does not rerun clustering, bootstrap resampling, score construction, external projection, spatial analysis, purity analysis, or manuscript analysis.

The frozen acceptance criteria remain unchanged:

- PAC rank <= 2 (lower PAC is better).
- Mean silhouette > 0.25.
- Mean silhouette within 0.02 of the best value in the same score space.
- Minimum cluster size >= 25.
- Minimum cluster median Jaccard >= 0.60.
- Overall median Jaccard >= 0.75.

Primary results and the 14 nonreference sensitivity variants are reported separately. The sensitivity denominator is 14, not 15. The primary reference is not counted as a sensitivity variant.

## Descriptive Rules

- Signed margins are calculated as observed minus the acceptance boundary, except PAC rank margin, which is `2 - PAC_rank`. Positive margins pass; negative margins fail. The strict silhouette rule still requires `> 0.25` when its numeric margin is zero.
- A primary near-miss k is defined descriptively as a k failing exactly one frozen criterion. No additional tolerance band is introduced and near-miss status does not alter acceptance.
- For sensitivity rows, the nearest failed criterion is the failed criterion with the signed normalized margin closest to zero. Normalization uses the frozen boundary for scale comparability and does not alter any pass field.
- Relative agreement rank is based on the descending mean of median ARI and median NMI. It is descriptive only. Relative rank does not contribute to absolute acceptance.
- Bootstrap failures are retained as observed. Failed replicates are not rerun and do not enter Jaccard summaries. A combination remains summary-eligible when at least one replicate completed, exactly as implemented in the frozen pipeline.

## Outcome A: Discrete Taxonomy Supported

This outcome requires at least one k with all of the following:

- primary pass is TRUE;
- at least one nonreference sensitivity variant passes;
- at least one leave-one-out variant passes;
- at least one alternative-distance variant passes;
- no implementation defect invalidates adjudication.

The above family requirements prevent a k with systematic leave-one-out or alternative-distance failure from being accepted. Thresholds may not be relaxed to create this outcome.

## Outcome B: Discrete Taxonomy Not Supported

This outcome applies when:

- primary pass is FALSE for every k=2-6;
- sensitivity pass count is 0/14 for every k=2-6; and
- no implementation defect invalidates adjudication.

The formal conclusion is:

`DISCRETE TAXONOMY NOT SUPPORTED UNDER FROZEN ACCEPTANCE CRITERIA`

Under this outcome, final k remains `NOT_SELECTED`, taxonomy remains `NOT_ASSIGNED`, continuous program abundance and immune-stromal-proliferative polarization become the primary model form, and clustering may be retained only as exploratory or descriptive supplemental analysis.

## Outcome C: Implementation Defect

This outcome is used only if the integrity audit identifies a genuine code or data implementation defect that prevents a valid read-only adjudication:

`ADJUDICATION BLOCKED BY IMPLEMENTATION DEFECT`

Unfavorable results, low agreement, failed bootstrap replicates retained under the frozen rules, or failure to meet acceptance thresholds are not implementation defects.

## Otherwise

If none of Outcomes A-C is satisfied, the result is `ADJUDICATION INDETERMINATE UNDER PREREGISTERED RULES` and requires a new, explicitly approved amendment. No final k or taxonomy may be assigned automatically.
