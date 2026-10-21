# Decision 04: Disclosure of PCA Generalization Collapse

**Status**: ADOPTED AND ENFORCED  
**Date**: Phase 7 / Step 12  
**Decision**: Transparently report the out-of-fold generalization collapse of low-rank PCA classifiers ($0.4388$ Macro-F1) rather than suppressing the negative finding or reporting in-sample training accuracy ($89.10\%$).

## Context
In Step 12, an exploratory logistic regression model fit on an 18-component PCA projection achieved 89.10% training accuracy. However, evaluating this same PCA pipeline under 5-fold family-held-out cross-validation caused family Macro-F1 to collapse to $0.4388$ (worse than majority baseline).

## Decision Rationale
1. **Epistemological Honesty**: Reporting $89.10\%$ training accuracy as a "linear low-rank compression" would be a fraudulent misrepresentation of model generalization.
2. **Scientific Insight**: The held-out collapse provides valuable biological evidence: viral replication timing coordinates cannot be compressed into linear top-variance components. The discriminative information is distributed across non-linear high-dimensional embedding space.
