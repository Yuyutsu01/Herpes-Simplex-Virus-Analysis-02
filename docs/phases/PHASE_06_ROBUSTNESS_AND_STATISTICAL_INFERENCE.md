# Phase 6: Robustness, Generalization & Statistical Inference

**Consolidated Steps**: Step 11 (Robustness Audit), Step 11.1 (Statistical Multiplicity Audit)  
**Primary Deliverables**: `reports/11_1_paired_hypothesis_tests.csv`, `reports/11_1_paired_bootstrap_differences.csv`, `reports/11_leave_one_family_out.csv`  

---

## 1. Scientific Objective
Evaluate the statistical significance and sensitivity of ESM-2 predictive performance under leave-one-family-out cross-validation, sister-species transfer, paired exact McNemar tests, 10,000 paired bootstrap resamples, and Holm-Bonferroni Family-Wise Error Rate (FWER) multiplicity control.

## 2. Statistical Multiplicity Audit
Paired comparisons against Ridge ESM-2 on family predictions ($N=63$):
- **vs. Majority Baseline**: Discordant $b=26, c=5$, raw $p = 0.000192$, Holm **$p_{\text{adj}} = 0.000769$** (Statistically Significant).
- **vs. Physicochemical Random Forest**: Discordant $b=23, c=9$, raw $p = 0.020062$, Holm **$p_{\text{adj}} = 0.060185$** (Borderline / Not Significant at $\alpha=0.05$).
- **vs. Physicochemical Logistic Regression**: Discordant $b=22, c=9$, raw $p = 0.029449$, Holm **$p_{\text{adj}} = 0.060185$** (Borderline / Not Significant at $\alpha=0.05$).
- **vs. Combined Model**: Discordant $b=2, c=0$, raw $p = 0.500000$, Holm **$p_{\text{adj}} = 0.500000$** (No difference).

## 3. Paired Bootstrap Confidence Intervals (10,000 Resamples)
- $\Delta \text{Macro-F1}$ vs Random Forest: Mean $= +0.2844$, 95% CI $[+0.0915, +0.4626]$.
- $\Delta \text{Macro-F1}$ vs Logistic Regression: Mean $= +0.2267$, 95% CI $[+0.0540, +0.3976]$.
- Strictly positive bootstrap intervals confirm a robust empirical trend despite missing strict FWER significance due to $N=63$.

## 4. Sensitivity Analyses
- **Leave-One-Family-Out (LOFO)**: 63 LOFO runs yielded Mean Macro-F1 = 0.7946 (Range: [0.778, 0.811]), proving no individual family drives performance.
- **Sister-Species Ortholog Transfer**: Evaluating an HSV-1-trained model on 7 HSV-2 ortholog families yielded 7/7 (100%) correct classification.

## 5. Validation & Test Outcomes
- Validated under `tests/test_11_robustness_audit.py` and `test_11_1_statistical_audit.py` (15 passing tests).
