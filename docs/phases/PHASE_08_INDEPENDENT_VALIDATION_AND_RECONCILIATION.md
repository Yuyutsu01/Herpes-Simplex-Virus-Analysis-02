# Phase 8: Independent Biological Validation & Reconciliation

**Consolidated Steps**: Step 13 (Biological Validation & Zero-Hardcoding), Step 13.1 (Discrepancy Reconciliation)  
**Primary Deliverables**: `reports/13_1_discrepancy_reconciliation.md`, `reports/13_1_independently_recomputed_metrics.csv`  

---

## 1. Scientific Objective
Audit all code, inputs, and outputs across the research pipeline to enforce zero-hardcoding, confirm deterministic reproducibility, and resolve historical narrative discrepancies that arose across earlier project steps.

## 2. Reconciled Discrepancies
1. **RS1 Sequence Count & Fold 0 Isolation**: Resolved ambiguity between 212 direct benchmark sequences vs. 129 cluster metadata entries. Proven that exactly 212 valid full-length sequences belong to `RS1` and all are strictly isolated in Fold 0.
2. **Model Hyperparameters**: Confirmed $C=1.0$, `penalty='l2'`, `solver='lbfgs'`, `max_iter=1000`, `random_state=42` across all Logistic Regression baselines.
3. **Feature Dimensionality**: Confirmed $D=8$ for physical, $D=480$ for ESM-2, and $D=488$ for combined models.
4. **Zero-Hardcoding**: All metrics derived by dynamically executing model training and prediction evaluation.

## 3. Validation & Test Outcomes
- Validated under `tests/test_13_independent_validation.py` and `test_13_1_discrepancy_reconciliation.py` (16 passing tests).
- All 687 out-of-fold prediction rows reproduced with absolute floating-point agreement ($< 1.11 \times 10^{-16}$).
