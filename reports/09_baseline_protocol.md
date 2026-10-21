# Step 9 Protocol: Prespecified Baseline Hierarchy

**Project**: *Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes*  
**Date**: 2026-10-21 14:14:59 UTC  
**Deliverable**: `reports/09_baseline_protocol.md`  

---

## 1. Baseline Model Hierarchy

Before fitting any advanced neural architectures or tuning classifiers, future modeling must evaluate four rigorous baselines:

### 1. Majority-Class Dummy Classifier
- Predicts the most frequent class in the training split.
- **Purpose**: Establishes the empirical performance floor under class imbalance.

### 2. Gene-Family Label Lookup (Diagnostic Diagnostic)
- Predicts the historical label of the gene family using exact symbol or homology group match.
- **Purpose**: Used **strictly as a diagnostic for data leakage**.
- If a model achieves near-100% test performance when family lookup is impossible under grouped holdouts, the model is exploiting unexpected leakage.

### 3. Physicochemical Property Baseline
- Linear logistic regression / Random Forest trained on 8 classical physical features:
  * Molecular weight, isoelectric point, GRAVY (hydropathy), aliphatic index, instability index, charged fraction, polar fraction, hydrophobic fraction.
- **Purpose**: Tests whether temporal kinetics correlates simply with bulk physical chemistry (e.g. acidic transcription factors vs basic capsid proteins).

### 4. Ridge-Regularized ESM-2 Linear Classifier
- Multinomial logistic regression with L2 regularization directly on the 480-dimensional mean-pooled ESM-2 embeddings.
- **Purpose**: Tests whether protein language representations encode temporal expression timing.

---

## 2. Prespecified Evaluation Metrics

Because of severe class imbalance at the gene-family level (N = 5 for IE vs N = 24 for L2):

1. **Macro-Averaged F1-Score (Macro-F1)**:
   Arithmetic mean of per-class F1-scores across all K classes: `Macro-F1 = (1/K) * sum(F1_k)`.
   *Gives equal weight to all classes regardless of sequence count.*
2. **Balanced Accuracy**:
   Arithmetic mean of per-class recalls: `Balanced Accuracy = (1/K) * sum(Recall_k)`.
3. **Per-Class Precision, Recall, and Confusion Matrices**.
4. **Gene-Family-Level Uncertainty Estimation**:
   * Metrics must be computed by averaging performance across the **63 independent gene families**, rather than sequence-weighted averages.
