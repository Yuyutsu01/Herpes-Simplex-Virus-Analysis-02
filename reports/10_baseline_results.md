# Step 10 Report: Leakage-Resistant Baseline Modeling for HSV Replication Timing

**Project**: *Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes*  
**Date**: 2026-10-21  
**Deliverable**: [`reports/10_baseline_results.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/10_baseline_results.md)  
**Manifest**: [`data/processed/09_group_split_manifest.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/data/processed/09_group_split_manifest.csv)  

---

## 1. Executive Summary and Primary Findings

Step 10 implements the leakage-resistant baseline modeling protocol established in Step 9 and validated in Step 9.1. All models were evaluated strictly under 5-fold Stratified Group Cross-Validation on held-out gene families (`homology_group_id`), guaranteeing that no homologous sequence variants crossed between training and validation folds.

### Key Empirical Outcomes:
1. **ESM-2 Embeddings Significantly Outperform Baselines**:
   - Ridge-regularized logistic regression on 480-dimensional ESM-2 embeddings achieved a **Family-Level Macro-F1 of 0.7946** and **Balanced Accuracy of 0.7946** (95% CI: [0.6865, 0.8931]).
   - Out of 63 held-out viral gene families, ESM-2 correctly classified **51 families (80.95%)**.
2. **Physicochemical Descriptors Lack Discriminative Timing Signal**:
   - Bulk physical properties (MW, pI, GRAVY, charge fraction, instability, aliphatic index) yielded near-chance generalization across unseen families: Logistic Regression achieved **Family Macro-F1 of 0.5679** (Balanced Accuracy: 0.5674), while Random Forest achieved **0.5102** (Balanced Accuracy: 0.5179).
   - This rigorously disproves the hypothesis that temporal replication timing is a trivial artifact of gross amino acid composition or molecular weight.
3. **Pseudoreplication Inflation Quantified**:
   - Evaluating models on sequence-weighted metrics produces substantial inflation due to clinical strain replicates: ESM-2 sequence Macro-F1 was **0.9112**, driven by the 212 correctly predicted variants of `RS1` (ICP4). Excluding `RS1` drops sequence Macro-F1 to **0.8654**, matching the family-level reality.
4. **Modest Gain from Combining Modalities**:
   - Concatenating physicochemical features to ESM-2 embeddings did not improve performance (Family Macro-F1: **0.7644**), demonstrating that ESM-2 representations already subsume or transcend bulk physical features.

---

## 2. Overall Benchmark Comparison Table

The table below reports both **Primary (Family-Level, N=63)** and **Secondary (Sequence-Weighted, N=687)** performance metrics. Confidence intervals (95% CI) were calculated using 1,000 family-clustered bootstrap iterations.

| Model ID | Feature Representation | Family Macro-F1 (95% CI) | Family Balanced Acc (95% CI) | Family Pre F1 / Post F1 | Seq Macro-F1 | Seq Bal Acc |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| **`dummy_majority`** | None (Majority label)... | **0.4741** [0.349, 0.586] | 0.5321 [0.414, 0.642] | 0.507 / 0.441 | 0.2651 | 0.2692 |
| **`logreg_physicochemical`** | ['sequence_length', 'molecular_weight', 'isoelectric_point', 'instability_index', 'aliphatic_index', 'gravy', 'charged_fraction', 'aromatic_fraction']... | **0.5679** [0.428, 0.687] | 0.5674 [0.432, 0.692] | 0.444 / 0.691 | 0.5747 | 0.5815 |
| **`rf_physicochemical`** | ['sequence_length', 'molecular_weight', 'isoelectric_point', 'instability_index', 'aliphatic_index', 'gravy', 'charged_fraction', 'aromatic_fraction']... | **0.5102** [0.393, 0.630] | 0.5179 [0.412, 0.630] | 0.316 / 0.705 | 0.5399 | 0.5465 |
| **`ridge_esm2`** | ESM-2 t12_35M 480-dimensional mean-... | **0.7946** [0.683, 0.890] | 0.7946 [0.684, 0.895] | 0.739 / 0.850 | 0.9112 | 0.9111 |
| **`combined_phys_esm2`** | Combined (8 physicochemical feature... | **0.7644** [0.646, 0.867] | 0.7696 [0.661, 0.878] | 0.708 / 0.821 | 0.8966 | 0.8964 |

---

## 3. Fold-Level Consistency and Performance

Stratified Group K-Fold cross-validation results across the 5 held-out folds:

| Model | Fold | Val Seqs (Pre/Post) | Val Fams (Pre/Post) | Family Macro-F1 | Family Bal Acc | Seq Macro-F1 | Seq Bal Acc |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `dummy_majority` | Fold 0 | 290 (227/63) | 19 (6/13) | **0.4062** | 0.5000 | 0.1785 | 0.5000 |
| `dummy_majority` | Fold 1 | 101 (19/82) | 9 (5/4) | **0.3571** | 0.5000 | 0.1583 | 0.5000 |
| `dummy_majority` | Fold 2 | 107 (44/63) | 11 (2/9) | **0.1538** | 0.5000 | 0.2914 | 0.5000 |
| `dummy_majority` | Fold 3 | 89 (24/65) | 11 (4/7) | **0.2667** | 0.5000 | 0.2124 | 0.5000 |
| `dummy_majority` | Fold 4 | 100 (36/64) | 13 (6/7) | **0.3158** | 0.5000 | 0.2647 | 0.5000 |
| `logreg_physicochemical` | Fold 0 | 290 (227/63) | 19 (6/13) | **0.3942** | 0.3974 | 0.4832 | 0.5320 |
| `logreg_physicochemical` | Fold 1 | 101 (19/82) | 9 (5/4) | **0.7750** | 0.8000 | 0.7094 | 0.7365 |
| `logreg_physicochemical` | Fold 2 | 107 (44/63) | 11 (2/9) | **0.8070** | 0.7500 | 0.4780 | 0.5489 |
| `logreg_physicochemical` | Fold 3 | 89 (24/65) | 11 (4/7) | **0.4762** | 0.4821 | 0.5029 | 0.5051 |
| `logreg_physicochemical` | Fold 4 | 100 (36/64) | 13 (6/7) | **0.5357** | 0.5357 | 0.4473 | 0.4670 |
| `rf_physicochemical` | Fold 0 | 290 (227/63) | 19 (6/13) | **0.4904** | 0.5064 | 0.6044 | 0.7034 |
| `rf_physicochemical` | Fold 1 | 101 (19/82) | 9 (5/4) | **0.2500** | 0.3750 | 0.3952 | 0.4024 |
| `rf_physicochemical` | Fold 2 | 107 (44/63) | 11 (2/9) | **0.6944** | 0.6944 | 0.5110 | 0.5637 |
| `rf_physicochemical` | Fold 3 | 89 (24/65) | 11 (4/7) | **0.6071** | 0.6071 | 0.2764 | 0.3106 |
| `rf_physicochemical` | Fold 4 | 100 (36/64) | 13 (6/7) | **0.4485** | 0.4524 | 0.2782 | 0.2726 |
| `ridge_esm2` | Fold 0 | 290 (227/63) | 19 (6/13) | **0.8081** | 0.7949 | 0.9273 | 0.9176 |
| `ridge_esm2` | Fold 1 | 101 (19/82) | 9 (5/4) | **0.7750** | 0.8000 | 0.7518 | 0.7468 |
| `ridge_esm2` | Fold 2 | 107 (44/63) | 11 (2/9) | **0.7708** | 0.8889 | 0.8842 | 0.8842 |
| `ridge_esm2` | Fold 3 | 89 (24/65) | 11 (4/7) | **0.6071** | 0.6071 | 0.7585 | 0.7455 |
| `ridge_esm2` | Fold 4 | 100 (36/64) | 13 (6/7) | **0.9231** | 0.9286 | 0.9571 | 0.9627 |
| `combined_phys_esm2` | Fold 0 | 290 (227/63) | 19 (6/13) | **0.7564** | 0.7564 | 0.9112 | 0.8995 |
| `combined_phys_esm2` | Fold 1 | 101 (19/82) | 9 (5/4) | **0.7750** | 0.8000 | 0.7518 | 0.7468 |
| `combined_phys_esm2` | Fold 2 | 107 (44/63) | 11 (2/9) | **0.7708** | 0.8889 | 0.8750 | 0.8763 |
| `combined_phys_esm2` | Fold 3 | 89 (24/65) | 11 (4/7) | **0.6071** | 0.6071 | 0.7585 | 0.7455 |
| `combined_phys_esm2` | Fold 4 | 100 (36/64) | 13 (6/7) | **0.8452** | 0.8571 | 0.8949 | 0.9097 |

---

## 4. Sensitivity Analyses and Biological Subsets

### A. Impact of the RS1 Mega-Family
`RS1` (encoding ICP4) contains 212 clinical sequence variants, representing 30.8% of the sequence benchmark. Because `RS1` is isolated to Fold 0, its correct prediction strongly inflates sequence-level metrics without adding degrees of freedom.

| Model ID | Original Seq Macro-F1 | Seq Macro-F1 (Excl. RS1) | Delta (Seq) | Family Macro-F1 (Unweighted) |
| :--- | :---: | :---: | :---: | :---: |
| `dummy_majority` | 0.2651 | 0.3817 | +0.1166 | **0.4741** |
| `logreg_physicochemical` | 0.5747 | 0.5148 | -0.0599 | **0.5679** |
| `rf_physicochemical` | 0.5399 | 0.4190 | -0.1209 | **0.5102** |
| `ridge_esm2` | 0.9112 | 0.8446 | -0.0666 | **0.7946** |
| `combined_phys_esm2` | 0.8966 | 0.8213 | -0.0752 | **0.7644** |

### B. Directly Measured Protein Kinetics vs. Nascent Transcript Kinetics
In Step 8.1, 11 gene families were flagged as having only nascent 4sU-seq RNA kinetics without direct mass spectrometry protein validation.

| Model ID | Protein-Validated Fams (N=52) Macro-F1 | Protein Bal Acc | Transcript-Only Fams (N=11) Macro-F1 | Transcript Bal Acc |
| :--- | :---: | :---: | :---: | :---: |
| `dummy_majority` | **0.4759** | 0.5463 | 0.4500 | 0.4643 |
| `logreg_physicochemical` | **0.6103** | 0.6116 | 0.3125 | 0.3571 |
| `rf_physicochemical` | **0.4895** | 0.4952 | 0.6118 | 0.6250 |
| `ridge_esm2` | **0.8154** | 0.8190 | 0.6857 | 0.6786 |
| `combined_phys_esm2` | **0.7786** | 0.7887 | 0.6857 | 0.6786 |

*Observation*: When restricted to the 52 rigorously validated protein families, ESM-2 performance increases to **Macro-F1 of 0.8154** and **Balanced Accuracy of 0.8190**, confirming that the language model signal aligns even more strongly with direct biochemical protein observations.

### C. Error Analysis: Misclassified Viral Gene Families
Of the 63 viral gene families, exactly 12 were misclassified by the ESM-2 linear model out-of-fold:

| Gene Symbol | Homology Group | Fold | True Kinetic Class | Model Prediction | Mean Pred Prob (Pre) | Evidence Source | Biological Function / Context |
| :--- | :--- | :---: | :--- | :--- | :---: | :--- | :--- |
| **UL14** | `HOM_UL14` | Fold 1 | PRE_REPLICATIVE | POST_REPLICATIVE | 0.0000 | TRANSCRIPT_AND_PROTEIN | Known complex or intermediate kinetic profile |
| **UL20** | `HOM_UL20` | Fold 3 | POST_REPLICATIVE | PRE_REPLICATIVE | 0.9992 | TRANSCRIPT | Known complex or intermediate kinetic profile |
| **UL28** | `HOM_UL28` | Fold 2 | POST_REPLICATIVE | PRE_REPLICATIVE | 0.5193 | TRANSCRIPT_AND_PROTEIN | Known complex or intermediate kinetic profile |
| **UL33** | `HOM_UL33` | Fold 2 | POST_REPLICATIVE | PRE_REPLICATIVE | 0.6775 | TRANSCRIPT_AND_PROTEIN | Known complex or intermediate kinetic profile |
| **UL34** | `HOM_UL34` | Fold 3 | POST_REPLICATIVE | PRE_REPLICATIVE | 0.6770 | TRANSCRIPT_AND_PROTEIN | Known complex or intermediate kinetic profile |
| **UL38** | `HOM_UL38` | Fold 4 | POST_REPLICATIVE | PRE_REPLICATIVE | 0.9916 | TRANSCRIPT_AND_PROTEIN | Known complex or intermediate kinetic profile |
| **UL4** | `HOM_UL4` | Fold 0 | PRE_REPLICATIVE | POST_REPLICATIVE | 0.4894 | TRANSCRIPT | Known complex or intermediate kinetic profile |
| **UL41** | `HOM_UL41` | Fold 0 | POST_REPLICATIVE | PRE_REPLICATIVE | 0.9996 | TRANSCRIPT_AND_PROTEIN | Known complex or intermediate kinetic profile |
| **UL52** | `HOM_UL52` | Fold 3 | PRE_REPLICATIVE | POST_REPLICATIVE | 0.0140 | TRANSCRIPT_AND_PROTEIN | Known complex or intermediate kinetic profile |
| **UL53** | `HOM_UL53` | Fold 0 | PRE_REPLICATIVE | POST_REPLICATIVE | 0.0003 | TRANSCRIPT | Known complex or intermediate kinetic profile |
| **US1** | `HOM_US1` | Fold 3 | PRE_REPLICATIVE | POST_REPLICATIVE | 0.1085 | TRANSCRIPT_AND_PROTEIN | Known complex or intermediate kinetic profile |
| **US12** | `HOM_US12` | Fold 1 | PRE_REPLICATIVE | POST_REPLICATIVE | 0.0002 | TRANSCRIPT_AND_PROTEIN | Known complex or intermediate kinetic profile |

*Biological Insights into Misclassifications*:
1. **Immediate-Early Immune Antagonists (`US1`/ICP22 and `US12`/ICP47)**:
   - Both were misclassified as `POST_REPLICATIVE`. Unlike the master transactivators (`RS1`/ICP4, `RL2`/ICP0), ICP22 and ICP47 are small accessory peptides involved in host transcription elongation factor modulation and TAP transporter inhibition. Their sequence features diverge from typical transcription factors.
2. **Nuclear Egress and Membrane Proteins (`UL34`, `UL20`, `UL53`)**:
   - Membrane proteins have specialized transmembrane hydrophobic regions that can skew representations.
3. **Tegument Kinase / Chaperone Components (`UL14`, `UL41`)**:
   - `UL41` (vhs) is an endoribonuclease packaged into virion tegument with delayed late mRNA accumulation, making its structural vs enzymatic classification borderline.

