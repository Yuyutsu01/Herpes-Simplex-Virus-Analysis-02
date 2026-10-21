# Step 11: Robustness, Generalization, and Scientific Validity Audit

**Project**: *Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes*  
**Date**: 2026-10-21  
**Deliverable**: [`reports/11_robustness_generalization_audit.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/11_robustness_generalization_audit.md)  
**Manifest**: [`reports/11_reproducibility_manifest.json`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/11_reproducibility_manifest.json)  

---

## 1. Executive Summary and Audit Scope

Step 11 executes a rigorous robustness, stability, and generalization audit of the Step 10/10.1 baseline benchmark. Rather than seeking higher predictive accuracy, this audit systematically stress-tests the benchmark under leave-one-family-out perturbations, fold-level variation, evidence-tier stratification, paired statistical hypothesis testing, and cross-species analysis.

### Principal Findings:
1. **Integrity and Provenance Verified**: All authoritative artifacts, row alignments, and split assignments were verified with zero data leakage.
2. **ESM-2 Stability Under LOFO Perturbation**: In leave-one-family-out sensitivity analysis across all 63 families, ESM-2 Family Macro-F1 remained tightly bounded within **[0.7886, 0.8082]** (median: **0.7926**, range: **0.0195**). No single family artificially creates the ESM-2 advantage.
3. **Paired Statistical Superiority Confirmed**: In paired family-level hypothesis tests, ESM-2 is statistically significantly superior to the fixed family-majority baseline ($p = 0.0347$), the fold sequence-majority baseline ($p = 0.00019$), linear physicochemical regression ($p = 0.0294$), and random forest physicochemical classification ($p = 0.0201$).
4. **Combined Model Dilution Explained**: Combining 8 bulk physical features with 480 ESM-2 features under isotropic L2 regularization causes slight performance degradation (0.7644 vs 0.7946 Macro-F1). The 8 non-informative features add parameter noise that misclassifies two borderline families (`UL28` and `UL4`).
5. **Cross-Species Unseen Benchmark Infeasible**: Due to severe experimental label sparsity in HSV-2 direct literature (only 7 families directly annotated in HSV-2 vs 63 in HSV-1), a cross-species unseen-family benchmark is **statistically infeasible ($N=7$)**. We explicitly refrain from manufacturing an invalid split.

---

## 2. Artifact Integrity and Provenance Verification

All authoritative files were audited for row count, checksum, and group disjointness:

| Artifact Label | Relative Path | SHA-256 (Truncated) | Rows | Status |
| :--- | :--- | :---: | :---: | :---: |
| **Step 8.1 Strict Temporal Benchmark** | `data\processed\08_1_strict_temporal_benchmark.csv` | `60ab4442c1ddc676...` | 3111 | `VERIFIED_PRESENT` |
| **Step 9 Group Split Manifest** | `data\processed\09_group_split_manifest.csv` | `44c46b4bdd13022b...` | 3111 | `VERIFIED_PRESENT` |
| **ESM-2 Representatives Metadata** | `data\processed\esm2_t12_representatives_metadata.csv` | `3b7f4f7fa543f24e...` | 3111 | `VERIFIED_PRESENT` |
| **Step 10 Out-of-Fold Predictions** | `reports\10_out_of_fold_predictions.csv` | `c6a2b24952656cb0...` | 3435 | `VERIFIED_PRESENT` |
| **Step 10 Baseline Metrics Summary** | `reports\10_baseline_metrics.csv` | `613aa1dac00a13fa...` | 5 | `VERIFIED_PRESENT` |
| **Step 10 Family-Level Predictions** | `reports\10_family_level_metrics.csv` | `ec6439932e23a55f...` | 315 | `VERIFIED_PRESENT` |
| **Step 10.1 Recomputed Metrics** | `reports\10_1_recomputed_metrics.csv` | `bbe5fad71807759f...` | 5 | `VERIFIED_PRESENT` |
| **Step 10.1 Label Sensitivity** | `reports\10_1_label_evidence_sensitivity.csv` | `c13c017cbc89f19a...` | 5 | `VERIFIED_PRESENT` |
| **Step 10.1 Reproducibility Decision** | `reports\10_1_reproducibility_decision.md` | `742dfda3a40b8970...` | -1 | `VERIFIED_PRESENT` |

---

## 3. Recomputed Primary Benchmark with Dual Majority Baselines

The benchmark recomputed on identical evaluation data ($N=63$ families, $N=687$ sequences) compares the models against both a fold-specific sequence-majority baseline and a constant family-majority baseline (always predicting POST):

| Model ID | Baseline Archetype | Family Macro-F1 | Family Bal Acc | Family Acc | Seq Macro-F1 | Seq Bal Acc |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| **`dummy_majority`** | Fold-Specific Sequence Majority | **0.4741** | 0.5321 | 0.4762 | 0.2651 | 0.2692 |
| **`logreg_physicochemical`** | Learned Model | **0.5679** | 0.5674 | 0.6032 | 0.5747 | 0.5815 |
| **`rf_physicochemical`** | Learned Model | **0.5102** | 0.5179 | 0.5873 | 0.5399 | 0.5465 |
| **`ridge_esm2`** | Learned Model | **0.7946** | 0.7946 | 0.8095 | 0.9112 | 0.9111 |
| **`combined_phys_esm2`** | Learned Model | **0.7644** | 0.7696 | 0.7778 | 0.8966 | 0.8964 |
| **`fixed_family_majority`** | Constant Family Majority (POST) | **0.3883** | 0.5000 | 0.6349 | 0.3375 | 0.5000 |

*Interpretation*: The fixed family-majority baseline achieves Balanced Accuracy of **0.5000** and Macro-F1 of **0.3883**, confirming that ESM-2 (0.7946 Macro-F1) provides substantial discriminative power beyond class imbalance.

---

## 4. Leave-One-Family-Out (LOFO) Sensitivity

To test whether any single viral locus disproportionately dictates performance, each of the 63 families was iteratively omitted:

| Model ID | Base Family MF1 | LOFO Min MF1 | LOFO Median MF1 | LOFO Max MF1 | Range (Max - Min) |
| :--- | :---: | :---: | :---: | :---: | :---: |
| `dummy_majority` | 0.4741 | 0.4643 | 0.4817 | 0.4817 | **0.0175** |
| `logreg_physicochemical` | 0.5679 | 0.5550 | 0.5640 | 0.5773 | **0.0223** |
| `rf_physicochemical` | 0.5102 | 0.4912 | 0.5067 | 0.5185 | **0.0273** |
| `ridge_esm2` | 0.7946 | 0.7886 | 0.7926 | 0.8082 | **0.0195** |
| `combined_phys_esm2` | 0.7644 | 0.7581 | 0.7621 | 0.7773 | **0.0192** |

### Most Influential Families for ESM-2:
- **Omission Most Increases Score (Misclassified Loci)**:
  * **`UL14`** (`HOM_UL14`, true PRE_REPLICATIVE): LOFO MF1 = 0.8082 (+0.0136)
  * **`UL20`** (`HOM_UL20`, true POST_REPLICATIVE): LOFO MF1 = 0.8082 (+0.0136)
  * **`UL38`** (`HOM_UL38`, true POST_REPLICATIVE): LOFO MF1 = 0.8082 (+0.0136)
  * **`UL33`** (`HOM_UL33`, true POST_REPLICATIVE): LOFO MF1 = 0.8082 (+0.0136)
  * **`UL34`** (`HOM_UL34`, true POST_REPLICATIVE): LOFO MF1 = 0.8082 (+0.0136)
- **Omission Most Decreases Score (Correctly Predicted Loci)**:
  * **`RL2`** (`HOM_RL2`, true PRE_REPLICATIVE): LOFO MF1 = 0.7886 (-0.0059)
  * **`RS1`** (`HOM_RS1`, true PRE_REPLICATIVE): LOFO MF1 = 0.7886 (-0.0059)
  * **`UL12`** (`HOM_UL12`, true PRE_REPLICATIVE): LOFO MF1 = 0.7886 (-0.0059)
  * **`UL2`** (`HOM_UL2`, true PRE_REPLICATIVE): LOFO MF1 = 0.7886 (-0.0059)
  * **`UL23`** (`HOM_UL23`, true PRE_REPLICATIVE): LOFO MF1 = 0.7886 (-0.0059)

*Methodological Note*: LOFO on fixed predictions evaluates evaluation-set sensitivity; it does not measure retraining influence.

---

## 5. Fold-Level Stability & Ascertainment Distortion

Cross-validation fold performance metrics across the 5 held-out folds:

| Model ID | Mean Fold MF1 (± Std) | Median Fold MF1 | Min / Max Fold MF1 | Fold 0 Seqs | Fold 0 MF1 (Fam / Seq) |
| :--- | :---: | :---: | :---: | :---: | :---: |
| `dummy_majority` | 0.2999 ± 0.0965 | 0.3158 | 0.1538 / 0.4062 | 290 | 0.4062 / 0.1785 |
| `logreg_physicochemical` | 0.5976 ± 0.1839 | 0.5357 | 0.3942 / 0.8070 | 290 | 0.3942 / 0.4832 |
| `rf_physicochemical` | 0.4981 ± 0.1692 | 0.4904 | 0.2500 / 0.6944 | 290 | 0.4904 / 0.6044 |
| `ridge_esm2` | 0.7768 ± 0.1131 | 0.7750 | 0.6071 / 0.9231 | 290 | 0.8081 / 0.9273 |
| `combined_phys_esm2` | 0.7509 ± 0.0874 | 0.7708 | 0.6071 / 0.8452 | 290 | 0.7564 / 0.9112 |

*Effect of RS1 in Fold 0*: `RS1` (ICP4) contributes 212 of the 290 validation sequences in Fold 0. Its correct classification elevates sequence-level accuracy in Fold 0 to **0.9273**, whereas family-level accuracy remains **0.8081**, confirming the necessity of grouping to prevent clinical ascertainment bias from masquerading as biological generalization.

---

## 6. Label Evidence Stratification Audit

Evaluating performance across evidence confidence tiers ([`11_label_evidence_sensitivity.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/11_label_evidence_sensitivity.csv)):

| Model ID | All Direct (N=63) MF1 | Protein-Validated (N=52) MF1 | Transcript-Only (N=11) MF1 | Delta (Prot - Tx) |
| :--- | :---: | :---: | :---: | :---: |
| `dummy_majority` | **0.4741** | **0.4759** | 0.4500 | +0.0259 |
| `logreg_physicochemical` | **0.5679** | **0.6103** | 0.3125 | +0.2978 |
| `rf_physicochemical` | **0.5102** | **0.4895** | 0.6118 | -0.1222 |
| `ridge_esm2` | **0.7946** | **0.8154** | 0.6857 | +0.1297 |
| `combined_phys_esm2` | **0.7644** | **0.7786** | 0.6857 | +0.0928 |

*Observation*: ESM-2 performance increases to **0.8154 Macro-F1** on the 52 mass-spectrometry protein-validated loci, whereas transcript-only loci show lower performance (0.6857). Transcript timing reflects nascent transcription rates rather than protein accumulation.

---

## 7. Paired Model Comparison & Statistical Significance

Paired hypothesis tests comparing ESM-2 against every baseline across the exact same 63 families ([`11_paired_model_comparison.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/11_paired_model_comparison.csv)):

| Comparison | ESM-2 Better / Worse | Exact McNemar p-value | Significant (p < 0.05) | Paired Delta Macro-F1 (95% CI) |
| :--- | :---: | :---: | :---: | :---: |
| **`ridge_esm2 vs fixed_family_majority`** | 17 / 6 | **3.4690e-02** | `True` | [+0.275, +0.515] |
| **`ridge_esm2 vs dummy_majority`** | 26 / 5 | **1.9220e-04** | `True` | [+0.173, +0.474] |
| **`ridge_esm2 vs logreg_physicochemical`** | 22 / 9 | **2.9449e-02** | `True` | [+0.054, +0.396] |
| **`ridge_esm2 vs rf_physicochemical`** | 23 / 9 | **2.0062e-02** | `True` | [+0.104, +0.477] |
| **`ridge_esm2 vs combined_phys_esm2`** | 2 / 0 | **5.0000e-01** | `False` | [+0.000, +0.077] |

### Why Does the Combined Model Underperform ESM-2 Alone?
The combined model adds 8 bulk physicochemical features to the 480 ESM-2 features under uniform L2 regularization ($C=1.0$). Because the 8 physical features have near-chance discriminative power (0.51-0.56 Macro-F1), they introduce non-informative variance that degrades decision boundary precision, causing misclassification of two borderline loci (`UL28` and `UL4`).

---

## 8. Species Representation and Cross-Species Generalization Audit

Audit of species composition across the direct cohort:
- Total direct sequences: 687 (422 HSV-1, 265 HSV-2).
- Total direct families: 63.
- **HSV-1 Direct Families**: **63** (100% of direct families have HSV-1 experimental evidence).
- **HSV-2 Direct Families**: **7** (Only 7 families have direct experimental kinetics in HSV-2 literature: `RL2`, `RS1`, `UL19`, `UL23`, `UL27`, `US4`, `US6`).
- **Shared Families**: 7.
- **HSV-2 Exclusive Families**: 0.

### Official Feasibility Determination:
CRITICAL FEASIBILITY AUDIT: A formal cross-species generalization benchmark (training exclusively on HSV-1 and testing on unseen HSV-2 families) is STATISTICALLY INFEASIBLE. In authoritative literature, only 7 gene families have direct experimental kinetic annotations in HSV-2. All 7 families are shared orthologs with HSV-1. Evaluating an unseen-family split across species would leave at most N=7 validation families (or fewer if partitioned), which is severely underpowered and scientifically indefensible. Cross-species transfer on the 7 shared families evaluates ortholog sequence conservation, NOT generalization to unseen viral gene families.

### Ortholog Transfer Performance on the 7 Shared Families:
| Homology Group | Gene Symbol | True Class | ESM-2 Prediction | Correct? | HSV-1 Seqs | HSV-2 Seqs |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| `HOM_RL2` | **RL2** | PRE_REPLICATIVE | PRE_REPLICATIVE | `True` | 12 | 9 |
| `HOM_RS1` | **RS1** | PRE_REPLICATIVE | PRE_REPLICATIVE | `True` | 79 | 133 |
| `HOM_UL19` | **UL19** | POST_REPLICATIVE | POST_REPLICATIVE | `True` | 7 | 18 |
| `HOM_UL23` | **UL23** | PRE_REPLICATIVE | PRE_REPLICATIVE | `True` | 11 | 28 |
| `HOM_UL27` | **UL27** | POST_REPLICATIVE | POST_REPLICATIVE | `True` | 9 | 23 |
| `HOM_US4` | **US4** | POST_REPLICATIVE | POST_REPLICATIVE | `True` | 10 | 39 |
| `HOM_US6` | **US6** | POST_REPLICATIVE | POST_REPLICATIVE | `True` | 2 | 15 |

*Diagnostic Result*: ESM-2 correctly predicts **7 out of 7 (100.0%)** of the shared cross-species ortholog families.

---

## 9. Defensible Scientific Boundaries

1. **No Causal Mechanism**: ESM-2 captures statistical evolutionary constraints, not the transcription factors, promoter motifs, or chromatin regulators that mechanistically direct kinetic cascades.
2. **No Generalization to Distant Herpesviruses**: The evaluation is restricted to *Simplexvirus* (HSV-1/HSV-2). Predictions cannot be extrapolated to Cytomegalovirus (Beta) or Epstein-Barr virus (Gamma).
3. **Sample Size Constraints ($N=63$)**: Generalization uncertainty spans [0.683, 0.890] for Macro-F1 due to the evolutionary repertoire size of 63 HSV genes.
