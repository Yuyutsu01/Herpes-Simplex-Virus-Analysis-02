# Step 13.1: Resolution of Benchmark Contradictions and Independent Result Verification

**Project**: *Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes*  
**Date**: 2026-10-21  
**Lead Auditor**: Senior Computational Virologist, Machine-Learning Auditor, Biostatistician  
**Deliverable Report**: [`reports/13_1_discrepancy_reconciliation.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/13_1_discrepancy_reconciliation.md)  
**Primary Deliverable Tables**:
- [`reports/13_1_record_and_fold_reconciliation.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/13_1_record_and_fold_reconciliation.csv)
- [`reports/13_1_model_configuration_audit.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/13_1_model_configuration_audit.csv)
- [`reports/13_1_feature_dimension_audit.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/13_1_feature_dimension_audit.csv)
- [`reports/13_1_hardcoding_audit.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/13_1_hardcoding_audit.csv)
- [`reports/13_1_independently_recomputed_metrics.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/13_1_independently_recomputed_metrics.csv)
- [`reports/13_1_prediction_reproduction_comparison.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/13_1_prediction_reproduction_comparison.csv)
- [`reports/13_1_statistical_reconciliation.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/13_1_statistical_reconciliation.csv)
- [`reports/13_1_scientific_claims_audit.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/13_1_scientific_claims_audit.csv)
- [`reports/13_1_reproducibility_manifest.json`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/13_1_reproducibility_manifest.json)
- [`tests/test_13_1_discrepancy_reconciliation.py`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/tests/test_13_1_discrepancy_reconciliation.py)

---

## 1. Executive Summary & Resolution of Apparent Contradictions

Step 13.1 resolves three critical contradictions identified between the textual summaries and underlying code outputs:

1. **Resolution of the RS1 Count and Fold Contradiction**:
   - **Contradiction**: Earlier project records documented $212$ `RS1` sequences assigned to Fold 0 validation. The textual narrative in the previous Step 13 report stated that `RS1` contains $447$ sequences and described Fold 3 validation as swelling to $471$ sequences.
   - **Resolution**: Direct inspection of the source data files ([`09_group_split_manifest.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/data/processed/09_group_split_manifest.csv), [`10_out_of_fold_predictions.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/10_out_of_fold_predictions.csv), [`13_retrained_oof_predictions.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/13_retrained_oof_predictions.csv), and [`13_family_split_audit.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/13_family_split_audit.csv)) proves that `HOM_RS1` contains **exactly 212 records**, and $100\%$ of them are allocated to **Fold 0 validation**. In Fold 0, validation sequences swell to $290$ ($212$ from `RS1` + $78$ from $18$ other families). In Fold 3, validation sequences total $89$ ($0$ from `RS1`). The numbers $447$ and $471$ were narrative textual errors (conflating Fold 3 training totals of $598$ with validation), while all underlying computational code and CSV files were consistently accurate.

2. **Resolution of Model Hyperparameter Discrepancy (Regularization $C$)**:
   - **Contradiction**: The Step 13 report narrative stated that ESM-2 and combined models were trained with $C=0.1$, while code excerpts and [`10_model_configuration.json`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/10_model_configuration.json) indicated $C=1.0$.
   - **Resolution**: Both [`scripts/10_temporal_baseline_models.py`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/scripts/10_temporal_baseline_models.py) and [`scripts/13_independent_biological_validation.py`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/scripts/13_independent_biological_validation.py) executed with **$C = 1.0$**. Re-retraining with $C=1.0$ matches the saved Step 10 and Step 13 predictions to machine precision (max probability difference $\le 1.11 \times 10^{-16}$, $0$ prediction mismatches). The mention of $C=0.1$ was a typographical error in the narrative description.

3. **Resolution of Physicochemical Feature Dimensions ($8$ vs. $28$ features)**:
   - **Contradiction**: Step 13 code excerpts listed $8$ physicochemical descriptors, while the report narrative claimed $28$ features and a $508$-dimensional combined representation.
   - **Resolution**: The supervised predictive baselines in Step 10 and Step 13 strictly use **8 classical bulk physicochemical descriptors** (`sequence_length`, `molecular_weight`, `isoelectric_point`, `instability_index`, `aliphatic_index`, `gravy`, `charged_fraction`, `aromatic_fraction`). The combined representation is **488 dimensions** ($8 + 480$). The number $28$ comes from Step 12's univariate sequence composition study ([`12_1_composition_inference_audit.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/12_1_composition_inference_audit.csv)), which evaluated the $20$ individual amino acid frequencies (`aa_A` through `aa_Y`) alongside the $8$ bulk descriptors ($20 + 8 = 28$). Those $28$ features were exploratory diagnostic associations, never inputs to the supervised classifiers.

---

## 2. Independent Model Reproduction & Record-by-Record Comparison

All $5$ models were retrained from scratch across all $5$ folds using verified parameters ($C=1.0$, $8$ physicochemical descriptors) and compared against saved Step 10 and Step 13 predictions by `record_id`:

| Model Identifier | Evaluated Records | Max $\Delta P$ vs Step 10 | Class Mismatches vs Step 10 | Max $\Delta P$ vs Step 13 | Class Mismatches vs Step 13 | Reproduction Verdict |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **`dummy_majority`** | 687 | 0.00e+00 | 0 | 0.00e+00 | 0 | **`EXACT_REPRODUCIBILITY`** |
| **`logreg_physicochemical`** | 687 | 1.11e-16 | 0 | 1.11e-16 | 0 | **`EXACT_REPRODUCIBILITY`** |
| **`rf_physicochemical`** | 687 | 1.11e-16 | 0 | 1.11e-16 | 0 | **`EXACT_REPRODUCIBILITY`** |
| **`ridge_esm2`** | 687 | 1.11e-16 | 0 | 1.11e-16 | 0 | **`EXACT_REPRODUCIBILITY`** |
| **`combined_phys_esm2`** | 687 | 1.11e-16 | 0 | 1.11e-16 | 0 | **`EXACT_REPRODUCIBILITY`** |

**Key Audit Finding**: Max probability divergence across all $3,435$ prediction records is $\le 1.11 \times 10^{-16}$ (floating-point precision limit), with **zero** classification mismatches.

---

## 3. Independently Recomputed Benchmark Metrics

Metrics were recomputed independently from out-of-fold predictions. Evaluation is reported at both sequence level ($N=687$) and family level ($N=63$), alongside evidence-tier subgroup breakdowns:

| Model Identifier | Evaluation Unit | Accuracy | Balanced Accuracy | Macro-F1 | Precision (PRE) | Recall (PRE) | Precision (POST) | Recall (POST) | Protein Subgroup Macro-F1 ($N=52$) | Transcript Subgroup Macro-F1 ($N=11$) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`dummy_majority`** | Sequence Representative (N=687) | 0.2707 | 0.2692 | 0.2651 | 0.3098 | 0.3514 | 0.2172 | 0.1869 | N/A | N/A |
| **`dummy_majority`** | Homology Gene Family (N=63) | 0.4762 | 0.5321 | 0.4741 | 0.3864 | 0.7391 | 0.6842 | 0.3250 | 0.4759 | 0.4500 |
| **`logreg_physicochemical`** | Sequence Representative (N=687) | 0.5793 | 0.5815 | 0.5747 | 0.6151 | 0.4657 | 0.5569 | 0.6973 | N/A | N/A |
| **`logreg_physicochemical`** | Homology Gene Family (N=63) | 0.6032 | 0.5674 | 0.5679 | 0.4545 | 0.4348 | 0.6829 | 0.7000 | 0.6103 | 0.3125 |
| **`rf_physicochemical`** | Sequence Representative (N=687) | 0.5444 | 0.5465 | 0.5399 | 0.5688 | 0.4371 | 0.5287 | 0.6558 | N/A | N/A |
| **`rf_physicochemical`** | Homology Gene Family (N=63) | 0.5873 | 0.5179 | 0.5102 | 0.4000 | 0.2609 | 0.6458 | 0.7750 | 0.4895 | 0.6118 |
| **`ridge_esm2`** | Sequence Representative (N=687) | 0.9112 | 0.9111 | 0.9112 | 0.9117 | 0.9143 | 0.9107 | 0.9080 | N/A | N/A |
| **`ridge_esm2`** | Homology Gene Family (N=63) | 0.8095 | 0.7946 | 0.7946 | 0.7391 | 0.7391 | 0.8500 | 0.8500 | 0.8154 | 0.6857 |
| **`combined_phys_esm2`** | Sequence Representative (N=687) | 0.8967 | 0.8964 | 0.8966 | 0.8908 | 0.9086 | 0.9030 | 0.8843 | N/A | N/A |
| **`combined_phys_esm2`** | Homology Gene Family (N=63) | 0.7778 | 0.7696 | 0.7644 | 0.6800 | 0.7391 | 0.8421 | 0.8000 | 0.7786 | 0.6857 |

---

## 4. Statistical Inference & Multiplicity Reconciliation

Paired model comparisons against `ridge_esm2` across the $63$ independent families under exact two-tailed McNemar tests, $10,000$ paired bootstrap resamples, and Holm-Bonferroni step-down correction:

| Paired Comparison | Correct ($b / c$) | Discordant Pairs ($b+c$) | Raw McNemar $p$ | Holm-Adjusted $p$ | FWER Significant ($lpha=0.05$) | $\Delta \text{Macro-F1}$ [95% CI] | $\Delta \text{Balanced Accuracy}$ [95% CI] |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`ridge_esm2 vs dummy_majority`** | 26 / 5 | 31 | 0.000192 | **0.000769** | **Yes** | +0.3205 [0.1758, 0.4731] | +0.2625 [0.1163, 0.4157] |
| **`ridge_esm2 vs rf_physicochemical`** | 23 / 9 | 32 | 0.020062 | **0.060185** | **No** | +0.2844 [0.0915, 0.4626] | +0.2766 [0.0929, 0.4464] |
| **`ridge_esm2 vs logreg_physicochemical`** | 22 / 9 | 31 | 0.029449 | **0.060185** | **No** | +0.2267 [0.0540, 0.3976] | +0.2272 [0.0517, 0.3974] |
| **`ridge_esm2 vs combined_phys_esm2`** | 2 / 0 | 2 | 0.500000 | **0.500000** | **No** | +0.0301 [0.0000, 0.0769] | +0.0250 [0.0000, 0.0638] |

**Statistical Conclusion**:
- Superiority over the majority baseline is statistically conclusive ($p_{\text{adj}} = 0.000769$).
- Against physicochemical models, ESM-2 demonstrates positive point gains ($+0.227$ to $+0.284$) and strictly positive $95\%$ bootstrap confidence intervals.
- However, after Holm-Bonferroni correction across the $4$ declared comparisons, raw p-values ($p = 0.0201$ and $p = 0.0294$) adjust to **$p_{\text{adj}} = 0.0602 > 0.05$**.
- This result is formally classified as **inconclusive under strict FWER control**, proving that claims of statistically proven superiority over physicochemical models are premature.

---

## 5. Epistemological Classification of Scientific Claims

| Research Claim | Epistemological Status | Justification & Quantitative Evidence |
| :--- | :---: | :--- |
| **ESM-2 Embedding Superiority Over Majority Baseline** | **`Reproduced`** | Family Macro-F1 = 0.7946 vs 0.4741; McNemar b=26, c=5; Holm p_adj = 0.000769 < 0.05. Directly demonstrated on 63 held-out HSV-1/HSV-2 gene families under 5-fold CV. |
| **ESM-2 Superiority Over Classical Physicochemical Features** | **`Supported within scope`** | Delta Macro-F1 = +0.227 to +0.284; 95% CIs strictly positive; raw p = 0.0201 and 0.0294; Holm p_adj = 0.0602. Does not achieve formal statistical significance at alpha=0.05 under strict 4-test FWER control. |
| **Zero Leakage Family-Held-Out Protocol** | **`Reproduced`** | All 5 cross-validation folds have train-val family overlap = 0 (confirmed across all steps). Restricted to the declared 63 homology groups in HSV-1/HSV-2 direct benchmark. |
| **Protein-Validated vs Transcript-Only Evidence Stratification** | **`Supported within scope`** | ESM-2 Macro-F1 is 0.8154 on protein-validated families (N=52) vs 0.6857 on transcript-only (N=11). Subgroup sample size for transcript-only is small (N=11). |
| **Low-Dimensional PCA Embedding Generalization** | **`Invalid (as claimed in Step 12)`** | Held-out 5-fold CV Macro-F1 collapses to 0.4388, whereas in-sample training was 0.8910. Evaluated on k=18 PCA components. |
| **Univariate Physicochemical Feature Associations** | **`Inconclusive`** | All 28 tested features fail Benjamini-Hochberg FDR correction at alpha=0.05 (Step 12.1 audit). 28 sequence features (20 AA frequencies + 8 bulk properties). |
| **Cross-Species Transfer to Unseen Herpesviruses** | **`Invalid (as framed)`** | HSV-2 N=7 test families are core orthologs of genes seen during HSV-1 training. Evaluates sister-species sequence divergence, not generalization to novel viral clades. |
| **External Biological Validation** | **`Not assessed`** | Zero independent external proteomic datasets exist in the repository; no synthetic data fabricated. Requires prospective mass-spectrometry time-course data under metabolic arrest. |

---

## 6. Final Scientific Verdict & Completion Criteria

1. **RS1 Discrepancy Status**: **RESOLVED**. `HOM_RS1` contains $212$ records assigned to Fold 0 validation. The numbers $447$ and $471$ in prior text were narrative errors; all underlying data files are consistent.
2. **Hyperparameter Status**: **RESOLVED**. All linear models execute with $C=1.0$. Exact reproduction verified.
3. **Feature Dimension Status**: **RESOLVED**. Supervised baseline models use $8$ bulk descriptors ($488$ combined). Univariate composition analysis used $28$ descriptors ($508$ combined).
4. **Data Lineage Status**: **CONFIRMED**. Complete trace verified from FASTA records to manifest, predictions, metrics, and report tables.
5. **Statistical Status**: **CONFIRMED**. Holm correction enforced. Physicochemical comparisons remain suggestive but non-significant under FWER ($lpha=0.05$).
6. **Publication Readiness**: The technical implementation and reproducibility are publication-grade. Manuscript submission requires updating narrative text to reflect the resolved $212$ RS1 count, $C=1.0$, $8$ physicochemical descriptors, and Holm-adjusted p-values ($p_{\text{adj}} = 0.0602$).
