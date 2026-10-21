# Phase 7: Representation Interpretation & Embedding Geometry

**Consolidated Steps**: Step 12 (Representation Interpretation), Step 12.1 (Geometry & PCA Audit)  
**Primary Deliverables**: `reports/12_1_embedding_geometry.csv`, `reports/12_1_pca_heldout_metrics.csv`  

---

## 1. Scientific Objective
Investigate the intrinsic dimensionality of the 480-dimensional ESM-2 embedding space, test whether low-dimensional linear projections generalize out-of-fold, and determine whether embeddings reflect causal kinetic mechanisms or distributed structural signatures.

## 2. Embedding Geometry & Effective Rank
- **Roy-Vetterli Effective Rank**:
  - *Standardized Correlation Matrix*: **13.63** (captures distributed degrees of freedom).
  - *Raw Unstandardized Covariance*: **5.93** (dominated by top 1-5 large variance directions).
- **Variance Distribution**: The first 18 principal components capture 90.0% of variance; 30 components capture 95.1%.

## 3. Low-Rank PCA Generalization Collapse
- In Step 12, an 18-component PCA classifier achieved 89.1% training accuracy.
- In Step 12.1, held-out 5-fold cross-validation revealed that out-of-fold Macro-F1 collapsed to **0.4388**.
- **Scientific Conclusion**: Unsupervised linear dimensionality reduction discards subtle discriminative coordinates necessary for out-of-fold viral phase discrimination. Viral replication timing is non-linear and distributed across high-dimensional latent space.

## 4. Non-Causal Boundary
ESM-2 embeddings correlate with globular catalytic domains in regulatory enzymes and transmembrane/repeat regions in structural proteins. The model leverages structural topology correlates, **not** causal transcriptional cascade mechanisms.

## 5. Validation & Test Outcomes
- Validated under `tests/test_12_representation_interpretation.py` and `test_12_1_interpretation_audit.py` (14 passing tests).
