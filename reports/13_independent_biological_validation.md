# Step 13: Independent Biological Validation, Data-Derived Metrics, and Zero-Hardcoding Audit

**Project**: *Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes*  
**Date**: 2026-10-21  
**Lead Auditor**: Senior Computational Virologist & Machine-Learning Auditor  
**Deliverable Report**: [`reports/13_independent_biological_validation.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/13_independent_biological_validation.md)  
**Associated Tables & Manifest**:
- [`reports/13_hardcoding_and_provenance_audit.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/13_hardcoding_and_provenance_audit.csv)
- [`reports/13_data_integrity_audit.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/13_data_integrity_audit.csv)
- [`reports/13_temporal_evidence_adjudication.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/13_temporal_evidence_adjudication.csv)
- [`reports/13_family_split_audit.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/13_family_split_audit.csv)
- [`reports/13_retrained_oof_predictions.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/13_retrained_oof_predictions.csv)
- [`reports/13_recomputed_model_metrics.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/13_recomputed_model_metrics.csv)
- [`reports/13_family_level_metrics.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/13_family_level_metrics.csv)
- [`reports/13_statistical_inference.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/13_statistical_inference.csv)
- [`reports/13_external_validation_feasibility.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/13_external_validation_feasibility.md)
- [`reports/13_reproducibility_manifest.json`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/13_reproducibility_manifest.json)

---

## 1. Executive Summary & Epistemological Stance

Step 13 conducts a comprehensive, data-derived biological validation and zero-hardcoding audit of the entire computational research pipeline. In accordance with strict reproducibility standards:
1. **Zero Hardcoded Numbers**: Every reported metric, p-value, sample size, and confidence interval is computed directly from raw and processed artifacts during runtime execution.
2. **Scratch Model Retraining**: All 5 baseline models (`dummy_majority`, `logreg_physicochemical`, `rf_physicochemical`, `ridge_esm2`, and `combined_phys_esm2`) were retrained from scratch across 5 cross-validation folds. Retrained probabilities matched the Step 10 saved predictions down to floating-point machine precision (maximum absolute difference $\le 1.11 \times 10^{-16}$), proving **100% deterministic reproducibility with zero seed drift**.
3. **Biological Independence Verified**: The 5-fold Stratified Group K-Fold split enforces **zero group leakage** ($\max |G_{\text{train}} \cap G_{\text{val}}| = 0$). All preprocessing (`StandardScaler`) was strictly fold-local.
4. **Correction of Overstated Claims**:
   - The Step 12 reported 89.1% PCA accuracy was an in-sample training metric. Under fold-local held-out cross-validation, PCA-compressed models collapse to chance level ($44.4\%$ Family Accuracy, $50.7\%$ Balanced Accuracy).
   - After Holm-Bonferroni correction across the 4 paired comparisons, only the dummy majority baseline is statistically significant ($p_{\text{adj}} = 0.00077$). Differences vs. physicochemical models ($p_{\text{adj}} \approx 0.0602$) do **not** survive strict FWER control at $\alpha = 0.05$.
   - The 7/7 result on HSV-2 represents an **ortholog sequence transfer diagnostic** for known genes, not novel cross-species family generalization.

---

## 2. Static No-Hardcoding and Provenance Audit

A systematic scan of scripts across Steps 8.1 through 12.1 categorized all findings into formal governance classes:

| Location | Classification | Description / Issue | Corrective Action / Resolution |
| :--- | :--- | :--- | :--- |
| `scripts/10_temporal_baseline_models.py:115` | `legitimate_configuration` | RANDOM_SEED = 42 | Retained as documented prespecified configuration parameter. |
| `scripts/10_temporal_baseline_models.py:169` | `legitimate_configuration` | binary_mapping = {'IE_ALPHA': 'PRE_REPLICATIVE', ...} | Retained as formal domain specification; validated against literature evidence sources. |
| `scripts/10_temporal_baseline_models.py:144-165` | `legitimate_configuration` | assert len(df_direct) == 687; assert num_fams == 63 | Retained as integrity assertions; dynamic calculations derive all subsequent metrics. |
| `scripts/10_temporal_baseline_models.py:220` | `legitimate_configuration` | LogisticRegression(C=1.0, penalty='l2') | Retained as prespecified baseline specification. |
| `scripts/12_representation_interpretation.py:347` | `suspicious` | pca_18_acc = accuracy_score(y, clf.predict(X_pca)) | Audited in Step 12.1 and Step 13: strictly re-evaluated using fold-local PCA within 5-fold CV, exposing held-out collapse to 44.4% accuracy. |
| `scripts/12_representation_interpretation.py:230` | `suspicious` | Condition number = 3.50e+15 reported in text | Audited in Step 12.1 and Step 13: computed Roy-Vetterli effective rank (13.6) and truncated condition number (1.7e7) to accurately represent geometry. |
| `reports/11_paired_model_comparison.csv:2` | `suspicious` | Raw McNemar p-values reported without FWER adjustment | Audited in Step 11.1 and Step 13: applied Holm-Bonferroni correction, revealing that only the dummy baseline comparison survives at alpha=0.05. |
| `reports/11_robustness_generalization_audit.md:76` | `suspicious` | Phrase 'cross-species transfer performance' used for 7/7 result | Audited in Step 11.1 and Step 13: phrase officially removed and reclassified as 'family-held-out ortholog diagnostic'; isolated transfer experiment trained strictly on HSV-1. |

*Table artifact*: [`reports/13_hardcoding_and_provenance_audit.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/13_hardcoding_and_provenance_audit.csv)

---

## 3. Data Identity, Embedding Alignment, and Negative Tests

The mapping between FASTA record IDs, metadata, split manifest, and embedding matrix rows was validated:

| Test Name | Expected Condition | Observed Value | Verdict |
| :--- | :--- | :--- | :---: |
| `direct_cohort_record_id_uniqueness` | 687 unique record IDs | total=687, unique=687 | `PASSED` |
| `embedding_matrix_integrity` | Shape (3111, 480), dtype float32, all values finite | shape=(3111, 480), dtype=float32, all_finite=True | `PASSED` |
| `embedding_row_index_validity` | All embedding_row_index values in [0, 3110] | min=15, max=3108 | `PASSED` |
| `negative_test_out_of_bounds_index_detection` | IndexError raised when accessing invalid embedding row index | IndexError raised successfully | `PASSED` |
| `negative_test_shuffled_row_mismatch_detection` | Exact record_id alignment fails upon deliberate row permutation | Mismatch detected upon row shuffling | `PASSED` |

*Table artifact*: [`reports/13_data_integrity_audit.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/13_data_integrity_audit.csv)

### Negative Testing Results:
- Deliberate injection of an out-of-bounds embedding row index (index 9999) was immediately caught with an `IndexError`.
- Row shuffling and index misalignment were detected with zero reliance on incidental array ordering.

---

## 4. Re-Adjudication of Temporal Evidence & Provenance

Temporal labels across 68 canonical genes were re-adjudicated against primary experimental literature:
- **Primary Direct Benchmark**: Exactly 63 gene families (52 dual protein-and-transcript validated, 11 transcript-only kinetics).
- **Experimental Sources**: Foundational metabolic arrest studies (`SRC_HONESS_1974`, `SRC_HONESS_1975`), modern 4sU-seq high-throughput transcriptomics (`SRC_RUTKOWSKI_2015`), SILAC mass spectrometry proteome time-courses (`SRC_KULEJ_2017`), ChIP-seq Pol II dynamics (`SRC_DREMEL_2019`), and authoritative consensus catalogs (`SRC_ROIZMAN_KNIPE_2013`).
- **Endpoint Integrity**: Dual-validated genes (`TRANSCRIPT_AND_PROTEIN`) establish that protein accumulation directly tracks transcript kinetics under uninhibited physiological conditions.

*Table artifact*: [`reports/13_temporal_evidence_adjudication.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/13_temporal_evidence_adjudication.csv)

---

## 5. Biological Independence & 5-Fold Group Splitting

The 5-fold Stratified Group K-Fold split was reconstructed from identifiers and clustering outputs:

| Fold | Train Seqs | Val Seqs | Train Fams | Val Fams (PRE / POST) | Leakage Status | RS1 Seqs in Val |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| Fold 0 | 397 | 290 | 44 | 19 (6 PRE / 13 POST) | `ZERO_LEAKAGE` | 212 |
| Fold 1 | 586 | 101 | 54 | 9 (5 PRE / 4 POST) | `ZERO_LEAKAGE` | 0 |
| Fold 2 | 580 | 107 | 52 | 11 (2 PRE / 9 POST) | `ZERO_LEAKAGE` | 0 |
| Fold 3 | 598 | 89 | 52 | 11 (4 PRE / 7 POST) | `ZERO_LEAKAGE` | 0 |
| Fold 4 | 587 | 100 | 50 | 13 (6 PRE / 7 POST) | `ZERO_LEAKAGE` | 0 |

*Table artifact*: [`reports/13_family_split_audit.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/13_family_split_audit.csv)

### Pseudoreplication Analysis:
- In Fold 0, the single `HOM_RS1` family contributes 212 of the 289 sequences ($73.4\%$).
- This structural asymmetry demonstrates why **family-level evaluation must be the primary metric**: sequence-level pooling is heavily distorted by single gene variants.

---

## 6. Baseline Model Retraining from Scratch

All 5 models were re-trained from scratch across the 5 folds with fold-local preprocessing:

| Model ID | Family Macro-F1 | Family Bal Acc | Family Acc | Prot Val F1 (N=52) | Tx Only F1 (N=11) | Seq Macro-F1 | Seq Bal Acc |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `dummy_majority` | **0.4741** | 0.5321 | 47.6% | 0.4759 | 0.4500 | 0.2651 | 0.2692 |
| `logreg_physicochemical` | **0.5679** | 0.5674 | 60.3% | 0.6103 | 0.3125 | 0.5747 | 0.5815 |
| `rf_physicochemical` | **0.5102** | 0.5179 | 58.7% | 0.4895 | 0.6118 | 0.5399 | 0.5465 |
| `ridge_esm2` | **0.7946** | 0.7946 | 81.0% | 0.8154 | 0.6857 | 0.9112 | 0.9111 |
| `combined_phys_esm2` | **0.7644** | 0.7696 | 77.8% | 0.7786 | 0.6857 | 0.8966 | 0.8964 |

*Table artifacts*: [`reports/13_family_level_metrics.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/13_family_level_metrics.csv), [`reports/13_recomputed_model_metrics.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/13_recomputed_model_metrics.csv)

### Reproducibility Finding:
Comparing retrained predictions against the saved Step 10 predictions row-by-row confirmed:
- `dummy_majority`: $0$ mismatches, max probability difference = $0.00$.
- `logreg_physicochemical`: $0$ mismatches, max probability difference = $1.11 \times 10^{-16}$.
- `rf_physicochemical`: $0$ mismatches, max probability difference = $1.11 \times 10^{-16}$.
- `ridge_esm2`: $0$ mismatches, max probability difference = $1.11 \times 10^{-16}$.
- `combined_phys_esm2`: $0$ mismatches, max probability difference = $1.11 \times 10^{-16}$.

---

## 7. Statistical Multiplicity and Paired Bootstrap (B = 10,000)

Exact McNemar tests and paired bootstrap differences ($B = 10,000$, seed=42) were recomputed on retrained family predictions:

| Comparison | Discordant (ESM2 vs Base) | Raw McNemar p | Holm-Adjusted p | FWER Status (alpha=0.05) | Obs Delta Macro-F1 (95% CI) | Obs Delta Bal Acc (95% CI) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| `ridge_esm2 vs dummy_majority` | 26 vs 5 | 1.9220e-04 | **7.6878e-04** | **Significant (p_adj < 0.05)** | +0.3205 [+0.1758, +0.4731] | +0.2625 [+0.1163, +0.4157] |
| `ridge_esm2 vs rf_physicochemical` | 23 vs 9 | 2.0062e-02 | **6.0185e-02** | Not Significant | +0.2844 [+0.0915, +0.4626] | +0.2766 [+0.0929, +0.4464] |
| `ridge_esm2 vs logreg_physicochemical` | 22 vs 9 | 2.9449e-02 | **6.0185e-02** | Not Significant | +0.2267 [+0.0540, +0.3976] | +0.2272 [+0.0517, +0.3974] |
| `ridge_esm2 vs combined_phys_esm2` | 2 vs 0 | 5.0000e-01 | **5.0000e-01** | Not Significant | +0.0301 [+0.0000, +0.0769] | +0.0250 [+0.0000, +0.0638] |

*Table artifact*: [`reports/13_statistical_inference.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/13_statistical_inference.csv)

### Multiplicity Finding:
- **Only the dummy baseline is rejected under strict FWER control** ($p_{\text{adj}} = 0.00077$).
- Comparisons against physicochemical models yield adjusted p-values of **$p_{\text{adj}} \approx 0.0602$**, which do **not** survive Holm-Bonferroni correction at $\alpha = 0.05$.
- However, 95% paired bootstrap intervals strictly exclude zero:
  - vs RF Phys: $\Delta \text{Macro-F1} = +0.2844$ [$+0.0915, +0.4626$]
  - vs LR Phys: $\Delta \text{Macro-F1} = +0.2267$ [$+0.0540, +0.3976$]

---

## 8. External Validation Feasibility Summary

As documented in [`reports/13_external_validation_feasibility.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/13_external_validation_feasibility.md):
- **Zero Genuine External Datasets in Repo**: All sequences derive from historical HSV-1/HSV-2 FASTA records.
- **Sister-Species Transfer vs Generalization**: The 7 HSV-2 test families are shared orthologs seen in HSV-1 during training. This evaluates ortholog transferability, not novel unseen-family generalization.
- **Future Protocol**: True external validation requires prospective mass spectrometry time-courses on independent primary clinical isolates or homologous profiling in other alphaherpesviruses (e.g. VZV / PRV).

---

## 9. Final Scientific Verdict and Classification of Project Claims

Every project claim is explicitly classified into one of the seven formal epistemological statuses:

### 1. Reproduced:
- **Ridge ESM-2 Baseline Performance**: Retrained from scratch to yield exact out-of-fold metrics: Family Balanced Accuracy = **$79.46\%$**, Macro-F1 = **$0.7946$**, Accuracy = **$80.95\%$ (51/63)**.
- **Deterministic Pipeline Execution**: Zero seed drift across all 5 models (max probability difference $\le 1.11 \times 10^{-16}$).
- **Group Disjointness**: Zero family leakage across all 5 cross-validation folds.

### 2. Supported Within Scope:
- **ESM-2 Superiority over Bulk Physicochemical Features**: Supported by consistent positive bootstrap margins ($\Delta \text{Macro-F1} > 0$), though marginally non-significant under strict Holm FWER control ($p_{\text{adj}} \approx 0.0602$).
- **Evidence Tier Consistency**: Higher performance on 52 protein-validated families (Macro-F1 = $0.8154$) than on 11 transcript-only families (Macro-F1 = $0.6857$).

### 3. Descriptive Only:
- **Low-Rank Intrinsic Geometry**: 18 principal components explain $95\%$ of embedding variance; Roy-Vetterli effective rank is $\approx 13.63$.
- **Bulk Amino Acid Trends**: Lysine ($d = +0.337$) and Glutamate ($d = +0.300$) show exploratory positive trends in pre-replicative proteins.

### 4. Inconclusive:
- **Sequence Composition Associations**: Zero features survive Benjamini-Hochberg FDR correction ($q \ge 0.4965$). Failing to reject the null is an inconclusive sample size constraint ($N=63$ families), not proof of no biological role.
- **Novel Viral Generalization**: Cannot be evaluated with available data ($N=7$ direct families in HSV-2).

### 5. Not Reproduced / Corrected:
- **89.1% PCA Model Accuracy**: Corrected from an in-sample training metric to a fold-local held-out cross-validation accuracy of **$44.4\%$** (Macro-F1 = $0.4394$).
- **$3.5 \times 10^{15}$ Condition Number**: Reclassified as a floating-point machine precision artifact on the smallest noise eigenvalue.
- **Raw McNemar Significance**: Unadjusted $p \approx 0.02$–$0.03$ against physicochemical models do **not** survive Holm-Bonferroni correction ($p_{\text{adj}} \approx 0.0602$).

### 6. Invalid:
- Attributing specific biological mechanisms (promoter activation, chromatin opening) directly to mature protein embedding dimensions.
- Calling the 7/7 ortholog diagnostic "cross-species transfer generalization".

### 7. Not Assessed:
- Token-level transformer attention attribution (multi-head attention tensors were not saved during Step 5 embedding extraction).

---

## 10. Required Synthesis & Manuscript Readiness Assessment

1. **Strongest Defensible Biological Conclusion**: Pre-trained protein language model embeddings (ESM-2) capture distributed structural and domain organization that distinguishes pre-replicative enzymes/regulators from post-replicative structural virion proteins among unseen HSV gene families (Macro-F1 $\approx 0.795$).
2. **Largest Remaining Threats to Validity**: Small family sample size ($N=63$), reliance on a single fixed 5-fold split, transcript-only kinetics for 11 families, and lack of external clinical validation.
3. **Reported Metrics That Changed**: PCA model held-out Macro-F1 revised from $89.1\%$ (in-sample) to $43.9\%$ (held-out); McNemar p-values adjusted for multiplicity ($p_{\text{adj}} \approx 0.0602$).
4. **Analyses That Could Not Be Completed**: Token-level attention attribution could not be performed because per-token tensors were not saved to disk.
5. **Publication Readiness Assessment**: The computational methodology, leakage controls, and code reproducibility are in a defensible state. However, manuscript submission requires **internal correction** to formalize the reclassified PCA and cross-species claims, followed by prospective independent experimental validation.
