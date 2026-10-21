# Step 12.1: Validation of Representation Interpretation and Feature Analysis

**Project**: *Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes*  
**Date**: 2026-10-21  
**Deliverable**: [`reports/12_1_interpretation_validation_audit.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/12_1_interpretation_validation_audit.md)  
**Associated Tables**:
- [`reports/12_1_data_alignment_checks.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/12_1_data_alignment_checks.csv)
- [`reports/12_1_pca_heldout_metrics.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/12_1_pca_heldout_metrics.csv)
- [`reports/12_1_embedding_geometry.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/12_1_embedding_geometry.csv)
- [`reports/12_1_composition_inference_audit.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/12_1_composition_inference_audit.csv)
- [`reports/12_1_coefficient_stability_audit.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/12_1_coefficient_stability_audit.csv)
- [`reports/12_1_reproducibility_manifest.json`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/12_1_reproducibility_manifest.json)


---

## 1. Executive Summary & Audit Conclusions

Step 12.1 independently audits all scripts, tables, and interpretation claims made in Step 12. Every major claim is evaluated against explicit data joins, out-of-fold generalization, and small-sample statistical theory.

### Principal Audit Findings:
1. **Data Alignment Verified with Zero Drift**: Explicit ID-based joins confirm exact 1-to-1 row alignment across FASTA record IDs, metadata, split manifest, and embedding rows ($N = 3,111$ total; $N = 687$ direct records across 63 unique homology groups).
2. **Correction of the PCA Interpretation Claim**: The 89.1% accuracy reported for PCA in Step 12 was strictly **in-sample training accuracy**. When evaluated under proper 5-fold cross-validation with PCA and scaling fitted fold-locally, PCA-compressed logistic regression ($k=18$) collapses to **$44.4\%$ Family Accuracy and $50.7\%$ Balanced Accuracy** (chance level). Unsupervised PCA discards the subtle, distributed directions required for temporal timing classification.
3. **Condition Number Demystified via Numerical Precision**: The previously reported condition number of $\sim 3.5 	imes 10^{15}$ is an artifact of dividing by near-zero numerical noise floor eigenvalues ($\lambda_{480} pprox 5.4 	imes 10^{-15}$). The true information-theoretic **Roy-Vetterli effective rank** of the standardized correlation matrix is **$	ext{erank} pprox 13.63$** degrees of freedom.
4. **Sequence Composition Non-Significance Reclassified as Inconclusive**: After Benjamini-Hochberg FDR correction across 28 features, zero features achieve $q < 0.05$. Cliff's Delta effect sizes (e.g., Lysine $d = +0.337$, Glutamate $d = +0.300$) show medium trends that become non-significant after sequence length adjustment ($p = 0.120$ and $p = 0.069$). This failure to reject the null is an **inconclusive sample-size limitation** ($N=63$), not definitive proof that amino-acid composition has no biological relationship to viral timing.
5. **Collinear Dimension Rank Swapping Documented**: Highly correlated embedding dimensions (e.g., Dim 304 and Dim 422, $r = +0.671$) fluctuate widely in magnitude across folds (Dim 422 ranks #1 in Fold 2, but #104 in Fold 0). Linear weights are model-dependent mathematical coordinates, not biological markers.

---

## 2. Section 1: Verification of Data Alignment and IDs

An exhaustive audit of all identifiers was performed using explicit ID-based merges ([`12_1_data_alignment_checks.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/12_1_data_alignment_checks.csv)):

| Audit Check | Source File | Reference File | Observed Value | Verdict |
| :--- | :--- | :--- | :--- | :---: |
| **Total Representatives** | `08_1_strict_temporal_benchmark.csv` | `esm2_t12_representatives_metadata.csv` | 3,111 rows in all files | `PASSED` |
| **Embedding Shape** | `esm2_t12_representatives_embeddings.npy` | Manifest | `(3111, 480)` float32 | `PASSED` |
| **Index Equality** | Manifest `embedding_row_index` | `numpy.arange(3111)` | Identical $0 \dots 3110$ | `PASSED` |
| **Record ID Match** | Benchmark `record_id` | Manifest `record_id` | Identical order (0 mismatches) | `PASSED` |
| **Direct Cohort** | Manifest `adjudicated_cohort` | Benchmark `adjudicated_cohort` | Exactly 687 records | `PASSED` |
| **Homology Groups** | Manifest `homology_group_id` | 63 Families | 23 PRE vs 40 POST (0 overlap) | `PASSED` |

---

## 3. Section 2: Validation of Family-Level Evaluation Closure

We audited whether multi-column grouping could inadvertently fragment families:
- **Homology Group Closure**: Every homology group maps to exactly 1 gene symbol, 1 temporal class, 1 measured molecule, and 1 cross-validation fold. Grouping by `homology_group_id` produces strictly 63 rows.
- **Aggregation Decision Rule**: Primary family predictions are recomputed by taking the **mean out-of-fold predicted probability** across all sequences in each family ($\hat{y} = 	ext{PRE}$ if $ar{p} \ge 0.5$).
- **Voting Comparison**: Sequence majority voting differs from probability averaging on **0 families** in this benchmark.
- **Stratified Metrics**:
  - **Overall Family Metrics ($N=63$)**: Accuracy = **$80.95\%$ (51/63)**, Balanced Accuracy = **$79.46\%$**, Macro-F1 = **$0.7946$**.
  - **Sequence-Level Metrics ($N=687$)**: Accuracy = **$81.66\%$**, Balanced Accuracy = **$81.83\%$**, Macro-F1 = **$0.8164$**.
  - **Protein-Validated Families ($N=52$)**: Macro-F1 = **$0.8154$**.
  - **Transcript-Only Families ($N=11$)**: Macro-F1 = **$0.6857$**.

---

## 4. Section 3: Correction of the PCA Generalization Claim

In Step 12, an in-sample accuracy of $89.1\%$ was reported for an 18-component PCA model. That analysis fitted and evaluated PCA on the same data. We re-evaluated PCA using **strict fold-local fitting** across the 5 family-held-out folds ([`12_1_pca_heldout_metrics.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/12_1_pca_heldout_metrics.csv)):

| Model Configuration | Components ($k$) | Protocol | Seq Accuracy | Fam Accuracy | Fam Balanced Acc | Fam Macro-F1 |
| :--- | :---: | :--- | :---: | :---: | :---: | :---: |
| `pca_18_in_sample` | 18 | **IN-SAMPLE (Step 12)** | **89.1%** | N/A | N/A | N/A |
| `pca_5_held_out` | 5 | 5-Fold Family Held-Out | 49.2% | 55.6% | 53.9% | 0.5358 |
| `pca_10_held_out` | 10 | 5-Fold Family Held-Out | 49.2% | 46.0% | 52.9% | 0.4536 |
| **`pca_18_held_out`** | **18** | **5-Fold Family Held-Out** | **47.3%** | **44.4%** | **50.7%** | **0.4394** |
| `pca_30_held_out` | 30 | 5-Fold Family Held-Out | 51.8% | 63.5% | 52.8% | 0.4849 |
| **`ridge_esm2_full_480`** | **480** | **5-Fold Family Held-Out** | **81.7%** | **81.0%** | **79.5%** | **0.7946** |

### Methodological Verdict:
> [!WARNING]
> **REVISION OF STEP 12 CLAIM**: The 89.1% PCA accuracy was an in-sample training metric. Under rigorous family-held-out evaluation, **PCA-compressed models collapse to chance level (Macro-F1 = 0.4394)**. Unsupervised PCA captures global evolutionary divergence and sequence length variance, but discards the specific, fine-grained coordinates that encode viral replication kinetics.

---

## 5. Section 4: Reassessment of Embedding Geometry

We recomputed the eigenvalue spectrum and investigated the condition number ([`12_1_embedding_geometry.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/12_1_embedding_geometry.csv)):

| Matrix Representation | Condition Number (Raw) | Truncated Condition Number ($\lambda_1 / \lambda_{479}$) | Roy-Vetterli Effective Rank | PCs for 95% Variance |
| :--- | :---: | :---: | :---: | :---: |
| **Standardized Correlation Matrix** | $3.43 	imes 10^{16}$ | **$1.72 	imes 10^{7}$** | **13.63** | 18 |
| **Raw Unstandardized Covariance** | $4.51 	imes 10^{16}$ | **$2.76 	imes 10^{7}$** | **5.93** | 18 |

### Explanation of the $10^{15}$ Condition Number:
The condition number $\lambda_{\max} / \lambda_{\min}$ is inflated because the smallest eigenvalue ($\lambda_{480} pprox 5.42 	imes 10^{-15}$) sits at the IEEE 754 floating-point machine precision noise floor. Truncating the single noise-floor eigenvalue yields a condition number of $\sim 1.7 	imes 10^7$. More informatively, the **Roy-Vetterli effective rank** is **$pprox 13.63$**, proving that the 480 embedding dimensions span approximately 14 independent degrees of freedom.

---

## 6. Section 5: Audit of Sequence Composition Inference

We evaluated 28 candidate features across the 63 families with effect sizes (Cliff's Delta with 95% bootstrap CIs) and length-adjusted ANCOVA ([`12_1_composition_inference_audit.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/12_1_composition_inference_audit.csv)):

| Feature | Mean PRE | Mean POST | Cliff's Delta [95% CI] | Raw MW $p$ | BH FDR $q$ | Length-Adjusted ANCOVA $p$ |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **`aa_K`** (Lysine) | $0.0184$ | $0.0136$ | $+0.337$ [$+0.05, +0.58$] | $0.0274$ | **$0.4965$** | $0.1204$ |
| **`aa_E`** (Glutamate) | $0.0485$ | $0.0418$ | $+0.300$ [$+0.03, +0.55$] | $0.0497$ | **$0.4965$** | $0.0688$ |
| **`aa_D`** (Aspartate) | $0.0533$ | $0.0472$ | $+0.276$ [$-0.01, +0.54$] | $0.0709$ | **$0.4965$** | $0.1478$ |
| **`isoelectric_point`** | $7.3131$ | $8.1342$ | $-0.276$ [$-0.54, +0.01$] | $0.0709$ | **$0.4965$** | $0.1310$ |
| **`sequence_length`** | $610.1$ | $505.5$ | $+0.250$ [$-0.05, +0.52$] | $0.1021$ | **$0.5200$** | N/A (Covariate) |
| **`molecular_weight`** | $63,027$ | $52,879$ | $+0.243$ [$-0.05, +0.51$] | $0.1114$ | **$0.5200$** | $0.3436$ |

### Scientific Verdict on Composition:
- **Status: INCONCLUSIVE DUE TO SAMPLE SIZE**.
- Zero features reach significance after FDR correction ($q \ge 0.4965$).
- When controlling for sequence length, Lysine ($p = 0.120$) and Glutamate ($p = 0.069$) fail to reach $lpha = 0.05$.
- *Epistemological Warning*: Failing to reject the null hypothesis in a sample of $N=63$ families does **not** prove that no composition differences exist. It indicates that the available sample size lacks statistical power to differentiate true subtle compositional bias from stochastic variation.

---

## 7. Section 6: Audit of Coefficient Stability and Collinear Rank Swapping

Inspection of linear model coefficients across the 5 cross-validation folds revealed multiple pairs of highly correlated dimensions that exchange importance across training subsets ([`12_1_coefficient_stability_audit.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/12_1_coefficient_stability_audit.csv)):

| Dimension | Mean Coef | Fold Rank Range (Min - Max) | Correlated Partner Dimension |
| :---: | :---: | :---: | :---: |
| **Dim 422** | $-0.3651$ | **Rank 1 to 104** (delta = 103) | Dim 304 ($r = +0.6705$) |
| **Dim 117** | $+0.3323$ | **Rank 2 to 164** (delta = 162) | Dim 464 ($r = +0.4789$) |
| **Dim 85** | $-0.3289$ | **Rank 4 to 163** (delta = 159) | Dim 115 ($r = +0.5120$) |
| **Dim 464** | $+0.3137$ | **Rank 3 to 100** (delta = 97) | Dim 117 ($r = +0.4789$) |
| **Dim 387** | $-0.3159$ | **Rank 7 to 80** (delta = 73) | Dim 279 ($r = +0.5641$) |

### Interpretation:
In Fold 2, Dim 422 is the #1 most important negative dimension, but in Fold 0 it drops to rank 104 as its correlated partner Dim 304 absorbs the weight. This confirms that individual linear coefficients are non-identifiable and reflect mathematical coordinate partitioning rather than distinct biological entities.

---

## 8. Required Classification of All Project Claims

| Claim / Finding | Epistemological Classification | Scientific Assessment |
| :--- | :--- | :--- |
| **Ridge ESM-2 achieves ~79.5% Balanced Accuracy on unseen families** | **Reproduced and supported** | Statistically validated under 5-fold cross-validation with fold-local scaling. |
| **18 PCs explain 95% of total embedding variance** | **Reproduced but descriptive only** | Accurately describes dataset variance, but these top components do not capture the temporal replication signal. |
| **89.1% accuracy for 18-component PCA model** | **Methodologically questionable (In-sample only)** | Overfitted training accuracy. Held-out accuracy collapses to 44.4% (chance level). |
| **Condition number is ~3.5 × 10^15** | **Methodologically questionable (Numerical artifact)** | Inflated by dividing by machine-precision noise floor eigenvalue (~ 10^-15). Effective rank is ~13.6. |
| **Amino acid composition has no association with replication class** | **Inconclusive due to sample size** | Moderate effect sizes observed ($d pprox 0.30$–$0.34$), but underpowered to survive FDR control at $N=63$. |
| **Individual embedding dimensions represent distinct biological switches** | **Not reproduced / Biologically unproven** | Disproven by collinear dimension rank swapping across folds. |
| **ESM-2 embeddings detect transcriptional promoter mechanisms** | **Biologically unfounded** | Mature protein sequences do not contain DNA promoter motifs or chromatin regulators. |
