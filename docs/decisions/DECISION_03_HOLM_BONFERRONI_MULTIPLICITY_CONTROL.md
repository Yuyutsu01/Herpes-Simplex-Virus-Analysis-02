# Decision 03: Holm-Bonferroni Multiplicity Control

**Status**: ADOPTED AND ENFORCED  
**Date**: Phase 6 / Step 11  
**Decision**: Enforce Holm-Bonferroni step-down Family-Wise Error Rate (FWER) correction across all paired model comparisons at $lpha = 0.05$.

## Context
When comparing ESM-2 against classical baselines (Majority, Physicochemical Random Forest, Physicochemical Logistic Regression, Combined), nominal unadjusted p-values showed $p = 0.0201$ vs RF and $p = 0.0294$ vs LogReg. In informal machine learning literature, these are frequently reported with asterisks ($*$) claiming "statistical superiority".

## Decision Rationale
1. **Multiplicity Discipline**: Testing 4 paired hypotheses on the same 63 gene families inflates false positive risk.
2. **Holm Adjustment**: Ordering the p-values and adjusting by rank yields $p_{	ext{adj}} = 0.060185$ for both physicochemical models.
3. **Scientific Boundary**: Because $0.0602 > 0.05$, formal superiority cannot be claimed at $lpha=0.05$. It must be described as a strong empirical trend with positive bootstrap intervals ($[+0.0915, +0.4626]$ vs RF; $[+0.0540, +0.3976]$ vs LogReg).
