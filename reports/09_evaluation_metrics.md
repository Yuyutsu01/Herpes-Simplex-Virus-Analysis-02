# Step 9 Protocol: Prespecified Statistical Metrics and Uncertainty Protocol

**Project**: *Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes*  
**Date**: 2026-10-21 14:14:59 UTC  
**Deliverable**: `reports/09_evaluation_metrics.md`  

---

## 1. Statistical Principles

1. **Sequence-Weighted Metrics Prohibited as Primary**:
   Evaluating accuracy across 687 sequences assigns 30.8% of total score weight to a single gene (`RS1`). A model predicting ICP4 correctly on 212 sequences would achieve high sequence accuracy even if it failed on all other 62 viral genes.
2. **Gene-Family Grouped Averaging (Macro-Evaluation)**:
   All validation metrics must compute prediction error at the **gene-family level** ($N=63$).
3. **Confidence Interval Estimation**:
   Bootstrap resampling across gene families ($B = 1,000$ iterations) to construct 95% confidence intervals on Macro-F1 and Balanced Accuracy.
