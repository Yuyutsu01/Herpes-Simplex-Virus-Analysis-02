# Step 12: Interpretable ESM-2 Representation Analysis for HSV Temporal Classes

**Project**: *Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes*  
**Date**: 2026-10-21  
**Deliverable**: [`reports/12_representation_interpretation.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/12_representation_interpretation.md)  
**Associated Tables**:
- [`reports/12_embedding_coefficient_stability.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/12_embedding_coefficient_stability.csv)
- [`reports/12_feature_association.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/12_feature_association.csv)
- [`reports/12_sequence_level_diagnostics.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/12_sequence_level_diagnostics.csv)
- [`reports/12_validation_and_limitations.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/12_validation_and_limitations.md)
- [`reports/12_reproducibility_manifest.json`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/12_reproducibility_manifest.json)

---

## 1. Executive Summary & Epistemological Stance

Step 12 investigates which sequence-derived representations contribute to ESM-2's ability to predict binary replication timing (`PRE_REPLICATIVE` = IE + E vs `POST_REPLICATIVE` = L1 + L2) across 63 directly annotated HSV gene families.

### Crucial Methodological Boundaries:
1. **Predictive Association vs. Biological Mechanism**: An embedding dimension, amino-acid frequency, or sequence descriptor that correlates with temporal classes does **not** control replication timing. Viral temporal cascades are governed in trans by viral transcriptional transactivators (e.g. VP16/UL48, ICP4/RS1, ICP27/UL54), host RNA polymerase II recruitment, and viral DNA replication, not by intrinsic physical timing properties encoded within structural proteins.
2. **ESM-2 Multicollinearity Prevents Single-Dimension Reification**: The 480-dimensional ESM-2 embedding space is intrinsically low-dimensional across the viral proteome. Just **18 principal components capture 95% of the total variance**, and the covariance matrix has a condition number of **3.50e+15**. Individual linear weights are non-identifiable and can rotate substantially across folds without changing prediction accuracy.
3. **Absence of Univariate Compositional Drivers**: None of the 20 canonical amino acid frequencies or 8 classical bulk physicochemical descriptors remain statistically significant after multiple testing correction across the 28 tests (minimum Benjamini-Hochberg FDR $q = 0.4965 > 0.05$). ESM-2 does not rely on a simple compositional heuristic (e.g. "early proteins are basic" or "late proteins are hydrophobic").
4. **Token-Level Artifact Boundary**: Only 480-dimensional mean-pooled representations were computed and saved in Step 5. No token-level hidden states or multi-head attention maps are present on disk. Token-level attributions are **not** fabricated.

---

## 2. Embedding-Level Attribution & Coefficient Stability

Linear classifiers (Ridge logistic regression, $C=1.0$) were evaluated across the 5 cross-validation folds. Positive coefficients correspond to predicting `PRE_REPLICATIVE`, while negative coefficients correspond to predicting `POST_REPLICATIVE`.

### Fold-to-Fold Coefficient Consistency:
- **Mean Pairwise Pearson Correlation**: **0.7547** (range: 0.6433 to 0.8481)
- **Mean Pairwise Spearman Rank Correlation**: **0.7333**

### Top Influential Embedding Dimensions:

#### Top 5 Dimensions Predictive of PRE_REPLICATIVE:
| Dimension | Mean Coef | Coef Std | SNR (|Mean|/Std) | Mean Rank (1-480) | Range of Ranks |
| :---: | :---: | :---: | :---: | :---: | :---: |
| **Dim 469** | +0.4200 | 0.1122 | 3.74 | 10.2 | 2-23 |
| **Dim 56** | +0.3519 | 0.0364 | 9.66 | 15.6 | 10-29 |
| **Dim 155** | +0.3643 | 0.0993 | 3.67 | 18.2 | 2-38 |
| **Dim 419** | +0.3522 | 0.0853 | 4.13 | 20.4 | 7-49 |
| **Dim 17** | +0.2920 | 0.0455 | 6.41 | 33.6 | 15-62 |

#### Top 5 Dimensions Predictive of POST_REPLICATIVE:
| Dimension | Mean Coef | Coef Std | SNR (|Mean|/Std) | Mean Rank (1-480) | Range of Ranks |
| :---: | :---: | :---: | :---: | :---: | :---: |
| **Dim 295** | -0.5557 | 0.1229 | 4.52 | 2.0 | 1-4 |
| **Dim 126** | -0.4202 | 0.0872 | 4.82 | 8.2 | 3-21 |
| **Dim 411** | -0.4084 | 0.0597 | 6.84 | 9.4 | 1-20 |
| **Dim 279** | -0.4170 | 0.1140 | 3.66 | 11.2 | 2-25 |
| **Dim 304** | -0.3744 | 0.0312 | 12.00 | 11.6 | 7-20 |

### Multicollinearity and Low-Rank Geometry:
- **Covariance Condition Number**: **3.50e+15**
- **Dimensionality Spectrum**:
  - 50% variance captured by **1 PCs**
  - 80% variance captured by **5 PCs**
  - 90% variance captured by **10 PCs**
  - 95% variance captured by **18 PCs**
- **Pairwise Correlation Distribution**:
  - Mean absolute off-diagonal correlation: **0.3545** (Max: **0.9523**)
  - Pairs with $|r| > 0.50$: **28.4%**
  - Pairs with $|r| > 0.70$: **8.7%**
- **Sparsity / Regularization Sensitivity**:
  - Under $L_1$ regularization (Lasso), only **30 out of 480 dimensions** receive non-zero weights.
  - A compressed model using only the first **18 PCA components** achieves **93.6%** training accuracy.
  - *Epistemological Consequence*: Because dimensions are highly correlated and low-rank, individual dimensions cannot be attributed to isolated biological functions. They act as distributed coordinates in a shared latent space.

---

## 3. Sequence-Level Composition & Confounder Analysis

To determine whether ESM-2's predictions merely reflect trivial sequence differences, 28 candidate features (20 amino acid frequencies + 8 classical physicochemical descriptors) were evaluated across the 63 families. Family-level mean aggregation was strictly enforced to prevent large sequence clusters (e.g. RS1 with 133 sequences) from dominating.

### Candidate Feature Statistical Summary:

| Feature | Category | Mean PRE ($N=23$) | Mean POST ($N=40$) | Difference | Raw Mann-Whitney $p$ | BH FDR $q$ | Length-Adjusted Wald $p$ |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **`aa_K`** | AminoAcid_Composition | 0.0184 | 0.0136 | +0.0048 | 2.7402e-02 | 0.4965 | 0.1233 |
| **`aa_E`** | AminoAcid_Composition | 0.0485 | 0.0418 | +0.0067 | 4.9652e-02 | 0.4965 | 0.0663 |
| **`aa_D`** | AminoAcid_Composition | 0.0533 | 0.0472 | +0.0061 | 7.0932e-02 | 0.4965 | 0.1314 |
| **`isoelectric_point`** | Physicochemical | 7.3131 | 8.1342 | -0.8211 | 7.0932e-02 | 0.4965 | 0.1249 |
| **`sequence_length`** | Physicochemical | 610.0731 | 505.4920 | +104.5811 | 1.0213e-01 | 0.5200 | 0.3693 |
| **`molecular_weight`** | Physicochemical | 63027.4534 | 52879.4832 | +10147.9702 | 1.1144e-01 | 0.5200 | 0.3436 |
| **`aa_P`** | AminoAcid_Composition | 0.0813 | 0.0958 | -0.0145 | 1.7277e-01 | 0.6911 | 0.1386 |
| **`aa_I`** | AminoAcid_Composition | 0.0323 | 0.0269 | +0.0055 | 2.1688e-01 | 0.6994 | 0.0901 |

### Multiple Testing & Confounder Findings:
1. **Zero Significant Compositional Features**: While Lysine (`aa_K`, raw $p = 0.0274$) and Glutamate (`aa_E`, raw $p = 0.0496$) show unadjusted trends toward enrichment in pre-replicative proteins, **neither survives Benjamini-Hochberg FDR correction ($q = 0.4965$) or Holm FWER correction**.
2. **Sequence Length Confounding**: Pre-replicative proteins are modestly longer on average ($610.1$ vs $505.5$ amino acids; raw $p = 0.1021$). When controlling for $\log(\text{sequence\_length})$ via logistic regression, `aa_K` ($p = 0.1232$) and `aa_E` ($p = 0.0662$) are entirely non-significant at $\alpha = 0.05$.
3. **Exploratory vs. Confirmatory Integrity**: No individual amino acid or bulk descriptor accounts for the timing distinction. ESM-2's predictive capacity stems from distributed domain and structural contextualization rather than primary amino acid abundance.

---

## 4. Attention and Token-Level Artifact Boundary

### Architectural Documentation:
- **Model**: `facebook/esm2_t12_35M_UR50D` (12 transformer layers, 480 hidden dimensions).
- **Available Artifact**: Mean-pooled embedding array `esm2_t12_representatives_embeddings.npy` of shape `(3111, 480)`.
- **Missing Artifacts**: Per-token hidden states, query-key-value projection matrices, and multi-head attention maps were **not** preserved to disk during Step 5 embedding generation.
- **Scientific Integrity Decision**: Mean-pooled vectors cannot be mathematically inverted to reconstruct per-token attention maps. Rather than synthesizing or hallucinating token-level "heatmaps", we explicitly document that token-level attribution requires a dedicated forward-pass caching pipeline.
- **Biological Non-Sequiturs**: Embedding representations cannot infer promoter binding motifs, transcription initiation complexes, or chromatin state. Those biological mechanisms operate on viral genomic DNA and regulatory nucleoprotein complexes, not on mature protein amino acid embeddings.

---

## 5. Stratified Robustness by Evidence Tier

Model performance was recomputed across the independent biological evidence tiers:

| Stratification Cohort | Total Families | PRE Families | POST Families | Family Accuracy | Balanced Accuracy | Macro-F1 |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **All Direct Families** | 63 | 23 | 40 | **80.95% (51/63)** | **79.46%** | **0.7946** |
| **Protein-Validated (Transcript + Protein)** | 52 | 19 | 33 | **82.69% (43/52)** | **81.90%** | **0.8154** |
| **Transcript-Only Evidence** | 11 | 4 | 7 | **72.73% (8/11)** | **67.86%** | **0.6857** |
| **HSV-1 $\to$ HSV-2 Ortholog Diagnostic** | 7 | 3 | 4 | **100.0% (7/7)** | **100.0%** | **1.0000** |

### Stratification Analysis:
- The classifier performs **higher on protein-validated families** (Macro-F1 = 0.8154) than on transcript-only families (Macro-F1 = 0.6857).
- This divergence aligns with biological reality: transcript kinetics occasionally diverge from functional protein accumulation due to post-transcriptional delays, leaky late synthesis, or delayed structural assembly.
- The 12 misclassified families out of 63 include:
  - `HOM_RL1` (ICP34.5 / gamma1 34.5): Complex transcriptional and translational timing.
  - `HOM_UL13` (Protein kinase): Leaky late expression with early transcription kinetics.
  - `HOM_UL36` (Large tegument protein): Giant structural protein (~3,100 aa) with multi-phase kinetics.

---

## 6. Synthesis and Reproducibility

The findings in Step 12 confirm that ESM-2 captures robust predictive structure that correlates with viral replication stages, but that this signal:
1. Is **distributed across a low-rank manifold** (18 PCs = 95% variance), ruling out single-dimension functional attribution.
2. Is **independent of trivial bulk physicochemical or amino acid composition biases**.
3. Represents **phenomenological statistical association**, not molecular causality.

All underlying tables, diagnostics, and metrics are reproducible via [`scripts/12_representation_interpretation.py`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/scripts/12_representation_interpretation.py).
