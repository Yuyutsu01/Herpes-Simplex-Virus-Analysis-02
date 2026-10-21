# Step 11.1: Statistical Inference and Cross-Species Claim Audit

**Project**: *Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes*  
**Date**: 2026-10-21  
**Deliverable**: [`reports/11_1_statistical_and_species_audit.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/11_1_statistical_and_species_audit.md)  
**Associated Tables**: [`reports/11_1_paired_hypothesis_tests.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/11_1_paired_hypothesis_tests.csv), [`reports/11_1_paired_bootstrap_differences.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/11_1_paired_bootstrap_differences.csv), [`reports/11_1_cross_species_transfer_experiment.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/11_1_cross_species_transfer_experiment.csv)  

---

## 1. Executive Summary and Audit Scope

Step 11.1 provides an independent statistical and methodological audit of the Step 11 robustness report. It rigorously resolves multiple hypothesis testing, estimates uncertainty via 10,000 paired bootstrap resamples, verifies the family-held-out protocol, and reclassifies the cross-species claims with a true isolated transfer experiment.

### Principal Audit Outcomes:
1. **Holm-Bonferroni Correction Applied**: While unadjusted McNemar tests suggested nominal superiority of ESM-2 over physicochemical models ($p = 0.020$–$0.029$) and the fixed family majority ($p = 0.035$), **only the comparison against the dummy sequence-majority baseline survives FWER control** ($p_{\text{adj}} = 0.00096$). Comparisons against physicochemical models yield adjusted p-values of **$p_{\text{adj}} \approx 0.080$–$0.088$**, which do NOT meet the strict family-wise significance threshold ($\alpha = 0.05$) across 5 comparisons.
2. **10,000 Paired Bootstrap Differences Computed**: The 95% bootstrap intervals for paired differences in Macro-F1 strictly exclude zero for all baseline models (ESM-2 vs LR Phys: **[+0.0543, +0.4041]**; ESM-2 vs RF Phys: **[+0.0994, +0.4646]**; ESM-2 vs Fixed Maj: **[+0.2815, +0.5195]**). These intervals quantify evaluation-sample uncertainty conditional on fixed predictions.
3. **Protocol Integrity Verified**: Preprocessing was confirmed to be fold-local (`StandardScaler` fitted exclusively on training sets). Zero validation leakage was detected.
4. **Cross-Species Claim Reclassified & True Transfer Tested**: The Step 11 7/7 result was reclassified from 'cross-species transfer' to a 'family-held-out ortholog diagnostic'. A true isolated HSV-1 $\to$ HSV-2 experiment was trained on 422 HSV-1 sequences and evaluated on 265 HSV-2 sequences, achieving **7/7 (100%) family accuracy** and **97.7% sequence accuracy**. However, because all 7 test families are shared orthologs seen during training, this demonstrates **ortholog sequence transferability**, NOT generalization to unseen viral gene families across species.

---

## 2. Section A: Paired Model Comparisons and Multiple Testing Correction

The McNemar test evaluates whether the discordant prediction counts between two models on the 63 families are symmetric ($H_0: b = c$). To control the Family-Wise Error Rate across the 5 prespecified comparisons, the **Holm-Bonferroni step-down procedure** was applied:

| Comparison | ESM-2 Only (b) | Baseline Only (c) | Raw McNemar p | Holm Rank | Holm-Adjusted p | Significant (FWER < 0.05) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **`ridge_esm2 vs fixed_family_majority`** | 17 | 6 | 3.4690e-02 | - | **8.8348e-02** | `False` |
| **`ridge_esm2 vs dummy_majority`** | 26 | 5 | 1.9220e-04 | - | **9.6098e-04** | `True` |
| **`ridge_esm2 vs logreg_physicochemical`** | 22 | 9 | 2.9449e-02 | - | **8.8348e-02** | `False` |
| **`ridge_esm2 vs rf_physicochemical`** | 23 | 9 | 2.0062e-02 | - | **8.0246e-02** | `False` |
| **`ridge_esm2 vs combined_phys_esm2`** | 2 | 0 | 5.0000e-01 | - | **5.0000e-01** | `False` |

### Statistical Interpretation of Adjusted Evidence:
- **Dummy Baseline**: Strongly rejected ($p_{\text{adj}} = 0.00096 < 0.001$). ESM-2 is decisively superior to sequence-majority guessing.
- **Physicochemical Models**: Raw tests showed nominal significance ($p = 0.0201$ for RF, $p = 0.0294$ for LogReg). However, after Holm adjustment, their p-values rise to **$0.0802$** and **$0.0883$**. Therefore, we **refrain from claiming statistically significant superiority under strict FWER control** at $\alpha = 0.05$.
- **Fixed Family Majority**: Similarly rises from raw $p = 0.0347$ to adjusted **$p_{\text{adj}} = 0.0883$**.

---

## 3. Section B: High-Precision Paired Bootstrap (B = 10,000)

Paired difference resampling was conducted over $B = 10,000$ iterations by sampling the 63 gene families jointly for both models ([`11_1_paired_bootstrap_differences.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/11_1_paired_bootstrap_differences.csv)):

| Comparison | Obs $\Delta$ Macro-F1 (95% CI) | Obs $\Delta$ Bal Acc (95% CI) | Obs $\Delta$ Accuracy (95% CI) |
| :--- | :---: | :---: | :---: |
| **`ridge_esm2 vs fixed_family_majority`** | **+0.4062** [+0.2815, +0.5195] | +0.2946 [+0.1872, +0.3973] | +0.1746 [+0.0317, +0.3175] |
| **`ridge_esm2 vs dummy_majority`** | **+0.3205** [+0.1758, +0.4731] | +0.2625 [+0.1163, +0.4157] | +0.3333 [+0.1905, +0.4921] |
| **`ridge_esm2 vs logreg_physicochemical`** | **+0.2267** [+0.0543, +0.4041] | +0.2272 [+0.0517, +0.4047] | +0.2063 [+0.0476, +0.3651] |
| **`ridge_esm2 vs rf_physicochemical`** | **+0.2844** [+0.0994, +0.4646] | +0.2766 [+0.0994, +0.4480] | +0.2222 [+0.0635, +0.3813] |
| **`ridge_esm2 vs combined_phys_esm2`** | **+0.0301** [+0.0000, +0.0765] | +0.0250 [+0.0000, +0.0638] | +0.0317 [+0.0000, +0.0794] |

### Handling of Class Representation in Resampling:
- Number of resamples with missing class: **0 / 10,000**.
- *Mathematical Rationale*: With $N_{\text{pre}}=23$ and $N_{\text{post}}=40$ in a sample of size 63, the binomial probability of drawing zero PRE instances is $(40/63)^{63} \approx 7.7 \times 10^{-13}$. Zero resamples were degenerate; no ad-hoc metric substitutions were required.
- *Scope of Uncertainty*: These bootstrap intervals quantify **evaluation-sample uncertainty** conditional on the saved predictions. They do NOT reflect model-retraining variability, label annotation noise, or sampling from external viral families.

---

## 4. Section C: Family-Held-Out Protocol Verification

The complete modeling pipeline and manifest were re-inspected:
- **no_validation_family_in_training**: `VERIFIED`
- **fold_local_preprocessing_only**: `VERIFIED`
- **prespecified_hyperparameters_no_cv_tuning**: `VERIFIED`
- **model_selection_prior_to_fitting**: `VERIFIED`
- **all_five_folds_contributed**: `VERIFIED`

*Limitations of Single 5-Fold Split*: The benchmark evaluates a single fixed 5-fold Stratified Group K-Fold split. It cannot estimate split-to-split sampling variability (variance across different random seeds). Uncertainty intervals derived from these predictions reflect evaluation-family sampling variability only.

---

## 5. Section D: Cross-Species Claim Audit and True Transfer Experiment

### Claim Reclassification:
RECLASSIFICATION REQUIRED: The 7/7 result reported in Step 11 was derived from the 5-fold cross-validation out-of-fold predictions on the 7 shared ortholog families. In that 5-fold CV, both HSV-1 and HSV-2 sequences of each family were held out together in validation. It was NOT an isolated model trained on HSV-1 and tested on HSV-2. The phrase 'cross-species transfer performance' is REMOVED and officially reclassified as a 'family-held-out ortholog diagnostic'.

### True Isolated HSV-1 $\to$ HSV-2 Cross-Species Transfer Experiment:
- **Training**: 422 HSV-1 sequences across 63 gene families. `StandardScaler` fitted strictly on HSV-1.
- **Testing**: 265 HSV-2 sequences across the 7 directly annotated families (`RL2`, `RS1`, `UL19`, `UL23`, `UL27`, `US4`, `US6`). Zero HSV-2 data in training.
- **Family-Level Performance ($N=7$)**: **Accuracy = 100.0%**, **Balanced Accuracy = 100.0%**, **Macro-F1 = 1.0000**.
- **Sequence-Level Performance ($N=265$)**: **Accuracy = 97.74% (259/265)**, **Balanced Accuracy = 98.24%**, **Macro-F1 = 0.9757**.

### Per-Family Test Results on HSV-2:

| Homology Group | Gene Symbol | True Class | Pred Class | Mean Prob (Pre) | HSV-2 Seqs | Correct? |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| `HOM_RL2` | **RL2** | PRE_REPLICATIVE | PRE_REPLICATIVE | 0.9995 | 9 | `True` |
| `HOM_RS1` | **RS1** | PRE_REPLICATIVE | PRE_REPLICATIVE | 0.9981 | 133 | `True` |
| `HOM_UL19` | **UL19** | POST_REPLICATIVE | POST_REPLICATIVE | 0.0229 | 18 | `True` |
| `HOM_UL23` | **UL23** | PRE_REPLICATIVE | PRE_REPLICATIVE | 0.7930 | 28 | `True` |
| `HOM_UL27` | **UL27** | POST_REPLICATIVE | POST_REPLICATIVE | 0.0210 | 23 | `True` |
| `HOM_US4` | **US4** | POST_REPLICATIVE | POST_REPLICATIVE | 0.0009 | 39 | `True` |
| `HOM_US6` | **US6** | POST_REPLICATIVE | POST_REPLICATIVE | 0.0251 | 15 | `True` |

### Defensible Scientific Boundary for Cross-Species Generalization:
Because the 7 HSV-2 test families are shared orthologs seen in HSV-1 during training, this experiment measures ortholog sequence transferability between sister species for KNOWN genes. It does NOT demonstrate generalization to novel, unseen viral gene families across species. Furthermore, because authoritative HSV-2 literature directly annotates only N=7 families, broad claims of cross-species generalization are scientifically indefensible.
