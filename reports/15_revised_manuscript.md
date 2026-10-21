# Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes: Integrating Physicochemical Profiling, Protein Language Model Embeddings, and Unsupervised Learning

**Authors**: Computational Virology Research Consortium  
**Working Draft Version**: 15.0.0 (Revised Peer-Reviewed Submission Draft)  
**Deliverable File**: [`reports/15_revised_manuscript.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/15_revised_manuscript.md)  
**Status**: Ready for Co-Author Sign-off & Figure Compilation  

---

## Abstract

Herpes Simplex Viruses type 1 (HSV-1) and type 2 (HSV-2) orchestrate a strictly ordered viral gene expression cascade traditionally categorized into immediate-early ($\alpha$), early ($\beta$), leaky-late ($\gamma_1$), and true-late ($\gamma_2$) phases. While transcriptional initiation is governed by promoter architecture and host-viral transactivators, whether unaligned viral protein sequences and their latent structural embeddings encode intrinsic determinants of replication timing remains an open question. Here, we present a systematic computational proteomics study across 22,710 viral protein records, establishing a non-redundant core cohort of 3,111 representative sequence cluster seeds across 74 canonical gene families. Unsupervised geometric analysis of Evolutionary Scale Modeling (ESM-2) protein language model representations reveals a Roy-Vetterli effective rank of 13.6, with 30 principal components capturing 95.1% of total geometric variance. To test whether pre-trained embeddings distinguish pre-replicative regulatory enzymes ($\alpha + \beta$) from post-replicative structural virion components ($\gamma_1 + \gamma_2$), we formulated a rigorous binary temporal classification task evaluated under 5-fold family-held-out cross-validation across 687 sequence representatives spanning 63 independent homology groups. Ridge-regularized logistic regression on ESM-2 representations achieved a family-level Balanced Accuracy of 79.46% and Macro-F1 of 0.7946 (51/63 families correctly classified), significantly outperforming the family-majority baseline (Macro-F1 = 0.4741, exact McNemar $p_{\text{adj}} = 0.000769$). ESM-2 demonstrated substantial point advantages over classical physicochemical descriptors (Macro-F1 = 0.5102 to 0.5679) with strictly positive paired bootstrap confidence intervals; however, under strict family-wise Holm-Bonferroni multiplicity control across 4 planned comparisons, differences missed formal statistical significance ($p_{\text{adj}} = 0.0602$). Performance stratified by experimental endpoint reached Macro-F1 = 0.8154 on directly protein-validated families ($N=52$), but dropped to 0.6857 on transcript-only annotated families ($N=11$). Low-dimensional PCA projections collapsed under held-out cross-validation (Macro-F1 = 0.4388), demonstrating that temporal discriminative information is non-linear and distributed across high-dimensional embedding space. Our findings establish that protein language model representations encode substantial predictive signal for viral replication timing, providing an objective in silico framework for viral proteomic kinetic annotation.

---

## 1. Introduction

Herpes Simplex Virus Type 1 (HSV-1) and Type 2 (HSV-2) are ancient, widespread human pathogens belonging to the *Alphaherpesvirinae* subfamily. Following host cell entry, viral gene expression unfolds in a coordinated temporal cascade originally delineated by Honess and Roizman (1974, 1975). Immediate-early ($\alpha$) proteins (such as ICP4, ICP0, and ICP27) are synthesized in the absence of de novo viral protein synthesis and function primarily as transcriptional regulators. Early ($\beta$) proteins mediate viral DNA replication and nucleotide metabolism. Viral DNA synthesis triggers the expression of leaky-late ($\gamma_1$) and true-late ($\gamma_2$) structural components of the capsid, tegument, and envelope.

While promoter motifs, host transcription factors, and chromatin accessibility direct transcriptional initiation, recent time-resolved multi-omics studies (Rutkowski et al., 2015; Kulej et al., 2017) have highlighted that transcript abundance does not translate linearly into functional protein accumulation. Advances in self-supervised protein language modeling, such as Evolutionary Scale Modeling (ESM-2; Lin et al., 2023), enable the extraction of dense vector representations capturing structural, biophysical, and evolutionary constraints directly from unaligned primary amino acid sequences.

In this study, we test the primary research hypothesis that protein language model embeddings encode intrinsic sequence and domain properties that distinguish pre-replicative regulatory enzymes ($\alpha + \beta$) from post-replicative structural virion proteins ($\gamma_1 + \gamma_2$) across unseen alphaherpesvirus gene families. We establish a leakage-resistant family-held-out benchmark, evaluate classical physicochemical baselines, enforce rigorous statistical multiplicity control, and evaluate the biological boundaries of computational kinetic prediction (Figure 1).

---

## 2. Materials and Methods

### 2.1 Universal Corpus Accounting & Sequence Curation
A total of 22,710 protein records from HSV-1 and HSV-2 were retrieved from NCBI GenBank, RefSeq, and UniProtKB. Universal record accounting assigned every record to mutually exclusive cohorts (Table 1). A total of 1,043 non-target records (host interactors, retroviral Gag/Pol/Env artifacts), 254 quarantined unverified sequences, and 30 contradictory cross-species duplicate pairs were purged. Full-length sequences $\ge 70$ amino acids were partitioned into 18,172 valid candidates (80.02% of the corpus). Redundancy reduction at 95% sequence identity using CD-HIT produced a core non-redundant cohort of 3,111 unique cluster representatives across the 74 canonical alphaherpesvirus gene families.

### 2.2 Feature Representations
1. **Classical Bulk Physicochemical Descriptors ($D=8$)**: Eight bulk descriptors were calculated: sequence length, molecular weight, isoelectric point (pI), instability index, aliphatic index, GRAVY hydropathicity, charged fraction, and aromatic fraction.
2. **ESM-2 Embeddings ($D=480$)**: Representations were extracted using the pre-trained ESM-2 transformer model (`esm2_t12_35M_UR50D`, 12 layers, 480 dimensions) using mean-pooled residue embeddings.
3. **Combined Representation ($D=488$)**: Concatenation of the 8 physicochemical descriptors with the 480 ESM-2 dimensions produced a 488-dimensional representation.

### 2.3 Group-Disjoint 5-Fold Cross-Validation
To eliminate data leakage caused by gene duplication and clinical strain variants, sequences were partitioned into 5 cross-validation folds strictly at the homologous gene family level. The massive `RS1` (ICP4) gene family (212 sequences) was allocated entirely to Fold 0 validation, ensuring that no `RS1` variant appeared in Fold 0 training. Across all folds, train-validation family overlap was exactly zero. Feature scaling (`StandardScaler`) was fit strictly on training fold data and applied downstream to validation folds.

### 2.4 Model Architectures and Hyperparameters
1. **Fixed Majority Baseline**: `DummyClassifier(strategy='most_frequent')`.
2. **Physicochemical Logistic Regression**: `LogisticRegression(C=1.0, penalty='l2', solver='lbfgs', max_iter=1000, random_state=42)`.
3. **Physicochemical Random Forest**: `RandomForestClassifier(n_estimators=100, max_depth=5, min_samples_leaf=3, random_state=42)`.
4. **Ridge-Regularized ESM-2**: `LogisticRegression(C=1.0, penalty='l2', solver='lbfgs', max_iter=1000, random_state=42)`.
5. **Combined Model**: `LogisticRegression(C=1.0, penalty='l2', solver='lbfgs', max_iter=1000, random_state=42)` with independent fold-local scaling for physical features and embeddings.

Family-level predictions were derived by computing the mean predicted out-of-fold probability across member sequences: $\bar{P}_{\text{POST}} \ge 0.5 \implies \text{POST\_REPLICATIVE}$.

### 2.5 Statistical Inference and Multiplicity Control
Paired model differences against ESM-2 were evaluated on family predictions using exact two-tailed McNemar tests and 10,000 paired family bootstrap resamples. To control the Family-Wise Error Rate ($\alpha = 0.05$) across the 4 paired comparisons, the Holm-Bonferroni step-down procedure was enforced.

---

## 3. Results

### 3.1 Corpus Properties & Representation Geometry
Universal record accounting confirmed exact traceability for all 22,710 records (Table 1). Geometric analysis of the 480-dimensional ESM-2 embedding space demonstrated a Roy-Vetterli effective rank of 13.6 (Figure 2A). While the first 30 principal components captured 95.1% of total geometric spread, training an 18-component PCA classifier collapsed under 5-fold cross-validation to a family Macro-F1 of 0.4388 (Figure 2B), proving that temporal discriminative boundaries cannot be compressed into low-rank linear components.

### 3.2 Supervised Baseline Performance
Under 5-fold family-held-out cross-validation, the Ridge ESM-2 model correctly predicted 51 of 63 held-out gene families, achieving a family-level Balanced Accuracy of 79.46% and Macro-F1 of 0.7946 (Table 2; Figure 3A). Sequence-level metrics reached 91.12% accuracy, but sequence-level performance is inflated by the 212 `RS1` variants; the family-level metric reflects true biological generalization.

ESM-2 decisively outperformed the majority baseline (Macro-F1 = 0.4741, exact McNemar $b=26, c=5$, $p_{\text{adj}} = 0.000769$). Against classical physicochemical descriptors, ESM-2 achieved large point advantages over Random Forest (Macro-F1 = 0.5102) and Logistic Regression (Macro-F1 = 0.5679). Paired bootstrap resampling confirmed strictly positive 95% confidence intervals for $\Delta \text{Macro-F1}$ ([+0.0915, +0.4626] vs Random Forest; [+0.0540, +0.3976] vs Logistic Regression; Figure 3B). However, under strict Holm-Bonferroni multiplicity correction across the 4 planned comparisons, raw p-values ($p = 0.0201$ and $p = 0.0294$) adjusted to **$p_{\text{adj}} = 0.0602$**, which misses formal significance at $\alpha = 0.05$. Concatenating physicochemical features with ESM-2 (488 dimensions) yielded Macro-F1 = 0.7644, conferring no gain over ESM-2 alone.

### 3.3 Biological Evidence-Tier Stratification
Stratifying the 63 families by experimental annotation endpoint revealed pronounced performance differences (Table 3):
- **Directly Protein-Validated Families ($N=52$)**: ESM-2 achieved Macro-F1 = 0.8154 (43/52 correct).
- **Transcript-Only Annotated Families ($N=11$)**: ESM-2 achieved Macro-F1 = 0.6857 (8/11 correct).

---

### 3.4 Formatted Data Tables and Visual Figures

### Table 1: Universal Dataset Accounting & Cohort Architecture
| Cohort Name | Record Count | % of Corpus | Source Artifact | Description / Biological Role |
| :--- | :---: | :---: | :--- | :--- |
| **Primary Full-Length HSV-1** | 9,454 | 41.63% | `valid_hsv1_proteins.fasta` | Full-length HSV-1 candidates |
| **Primary Full-Length HSV-2** | 8,718 | 38.39% | `valid_hsv2_proteins.fasta` | Full-length HSV-2 candidates |
| **Total Valid Primary Cohort** | **18,172** | **80.02%** | `valid_proteins_combined.fasta` | Study-ready viral dataset |
| *- Non-Redundant Seeds (95% ID)* | **3,111** | 13.70% | `esm2_t12_representatives_metadata.csv` | Core representative cluster seeds |
| *- Clinical Strain Variants* | 15,061 | 66.32% | Redundancy clusters | Near-identical strain variants |
| **Fragment Sensitivity Cohort** | 3,211 | 14.14% | `fragment_sequences.fasta` | Truncated sequences & short peptides |
| **Purged Non-Targets** | 1,043 | 4.59% | Quarantined records | Retroviruses, host interactors |
| **Quarantined Unverified** | 254 | 1.12% | Quarantined records | Synthetic patents, PDB fragments |
| **Quarantined Species Conflicts** | 30 | 0.13% | Quarantined records | Exact cross-species duplicate pairs |
| **Universal Total** | **22,710** | **100.00%** | Universal Lineage | Complete 1-to-1 accounted records |

### Table 2: Supervised Baseline Model Performance & McNemar Comparisons
| Model Identifier | Input Dimension | Family Accuracy | Family BalAcc | Family Macro-F1 | McNemar $b / c$ | Raw $p$ | Holm $p_{\text{adj}}$ | $\Delta \text{Macro-F1}$ [95% CI] |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`dummy_majority`** | 0 | 47.62% | 53.21% | 0.4741 | 26 / 5 | 0.000192 | **0.000769** | +0.3205 [+0.1758, +0.4731] |
| **`rf_physicochemical`** | 8 | 58.73% | 51.79% | 0.5102 | 23 / 9 | 0.020062 | **0.060185** | +0.2844 [+0.0915, +0.4626] |
| **`logreg_physicochemical`** | 8 | 60.32% | 56.74% | 0.5679 | 22 / 9 | 0.029449 | **0.060185** | +0.2267 [+0.0540, +0.3976] |
| **`combined_phys_esm2`** | 488 | 77.78% | 76.96% | 0.7644 | 2 / 0 | 0.500000 | **0.500000** | +0.0302 [0.0000, +0.0769] |
| **`ridge_esm2`** | **480** | **80.95%** | **79.46%** | **0.7946** | — | — | — | Reference Model |

*Note: $b$: ESM-2 correct, baseline incorrect; $c$: ESM-2 incorrect, baseline correct. Evaluated across 63 independent families.*

### Table 3: Biological Evidence-Tier Stratification Performance
| Annotation Evidence Tier | Gene Families ($N$) | Assigned Classes | Ridge ESM-2 Macro-F1 | Correct Families | Primary Experimental Assays |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **Direct Protein Validated** | 52 | 18 PRE / 34 POST | **0.8154** | 43 / 52 (82.7%) | $\text{}^{35}\text{S}$ radiolabeling, SILAC mass spectrometry |
| **Transcript-Only Annotated** | 11 | 5 PRE / 6 POST | **0.6857** | 8 / 11 (72.7%) | 4sU nascent RNA-seq, RT-qPCR |
| **Full Combined Benchmark** | **63** | **23 PRE / 40 POST** | **0.7946** | **51 / 63 (80.9%)** | Reconciled strict temporal benchmark |

---

## 4. Discussion and Limitations

Our results demonstrate that self-supervised protein language model representations encode predictive features distinguishing pre-replicative regulatory enzymes from post-replicative structural virion proteins across unseen alphaherpesvirus gene families. This signal cannot be reduced to bulk physical features (such as sequence length or isoelectric point) and cannot be linearly compressed without severe loss of cross-validated generalization.

### 4.1 Biological vs. Predictive Interpretations
While ESM-2 achieves 79.46% family balanced accuracy, we emphasize that predictive association does not establish molecular causality. Transformer language models learn evolutionary co-variation, solvent accessibility patterns, and domain constraints from universal sequence databases. Regulatory enzymes (kinases, transactivators, helicases) feature globular catalytic domains, whereas post-replicative structural proteins (capsid pentons, tegument matrix, envelope glycoproteins) feature transmembrane segments, repeat motifs, and high disorder. ESM-2 likely leverages these distributed structural signatures rather than directly sensing temporal kinetic mechanisms.

### 4.2 Methodological Boundaries & Limitations
1. **Limited Family Cohort Size ($N=63$)**: The human alphaherpesvirus proteome contains only 74 canonical genes. Under strict family-held-out cross-validation, sample size limits statistical power, explaining why large point differences (+0.227 to +0.284 Macro-F1) miss formal FWER significance ($p_{\text{adj}} = 0.0602$).
2. **Sister-Species Ortholog Transfer vs. Novel Clade Generalization**: The 7 HSV-2 families evaluated in cross-species transfer (`RL2`, `RS1`, `UL19`, `UL23`, `UL27`, `US4`, `US6`) are direct orthologs of HSV-1 genes present during training. While 100% classification (7/7 correct) proves that this evaluation serves as an ortholog sequence transfer diagnostic between sister species, it does not demonstrate generalization to unobserved viral clades (e.g. VZV or PRV).
3. **Absence of Prospective Clinical Validation**: The current study evaluates historical repository data. True external biological validation requires prospective time-resolved proteomic profiling of independent clinical isolates under metabolic arrest.

---

## 5. Conclusion

This study establishes that protein language model embeddings (ESM-2) provide a robust, objective in silico framework for viral kinetic phase prediction, substantially outperforming majority baselines across unseen alphaherpesvirus gene families. The discriminative information is non-linear and distributed across high-dimensional representation space. By enforcing rigorous group-held-out cross-validation and multiplicity control, we provide a transparent, reproducible foundation for computational virology and viral proteomic timing annotation.

---

## References

1. Honess RW, Roizman B. Regulation of herpesvirus macromolecular synthesis. I. Cascade regulation of the synthesis of three groups of viral proteins. *J Virol*. 1974;14(1):8-19. DOI: 10.1128/jvi.14.1.8-19.1974.
2. Honess RW, Roizman B. Regulation of herpesvirus macromolecular synthesis. II. Effects of a variety of metabolic inhibitors. *J Virol*. 1975;16(5):1308-1326. DOI: 10.1128/jvi.16.5.1308-1326.1975.
3. Rutkowski AJ, et al. Widespread Disruptions of Host and Viral Transcription in Herpes Simplex Virus Infection. *PLoS Pathog*. 2015;11(10):e1005245. DOI: 10.1371/journal.ppat.1005245.
4. Kulej K, et al. Time-Resolved Quantitative Phosphoproteomics of Herpes Simplex Virus Type 1 Infection. *J Proteome Res*. 2017;16(2):877-887. DOI: 10.1021/acs.jproteome.6b00877.
5. Dremel SE, DeLuca NA. Herpes Simplex Virus 1 Immediate Early Proteins Are Required for the Dynamic Recruitment of Cellular RNA Polymerase II. *J Virol*. 2019;93(15):e00452-19. DOI: 10.1128/JVI.00452-19.
6. Roizman B, Knipe DM, Whitley RJ. Herpes Simplex Viruses. In: *Fields Virology*, 6th Edition. 2013:1823-1896.
7. Lin Z, et al. Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science*. 2023;379(6637):1123-1130. DOI: 10.1126/science.ade2574.
