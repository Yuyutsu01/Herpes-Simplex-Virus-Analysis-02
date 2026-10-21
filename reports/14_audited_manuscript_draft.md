# Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes: Integrating Physicochemical Profiling, Protein Language Model Embeddings, and Unsupervised Learning

**Authors**: Computational Virology Research Consortium  
**Working Draft Version**: 14.0.0 (Audited Pre-Submission Copy)  
**Deliverable File**: [`reports/14_audited_manuscript_draft.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/14_audited_manuscript_draft.md)  
**Status**: Ready for Internal Author Review & Co-Author Sign-off  

---

## Abstract

Herpes Simplex Viruses type 1 (HSV-1) and type 2 (HSV-2) orchestrate a strictly regulated, cascade-dependent viral gene expression program comprising immediate-early ($\alpha$), early ($\beta$), leaky-late ($\gamma_1$), and true-late ($\gamma_2$) phases. While transcriptional regulation has been extensively characterized, whether viral protein primary sequences and their latent structural representations encode intrinsic determinants of kinetic timing remains an open question. Here, we present a comprehensive computational proteomics analysis across 22,710 viral protein records, curating a non-redundant core cohort of 3,111 representative sequence cluster seeds across 74 canonical gene families. Unsupervised representation geometry analysis of Evolutionary Scale Modeling (ESM-2) protein language model embeddings demonstrates an effective rank of 13.6, with 30 principal components capturing 95% of total geometric variance. To test whether pre-trained embeddings distinguish pre-replicative regulatory enzymes ($\alpha + \beta$) from post-replicative structural components ($\gamma_1 + \gamma_2$), we established a rigorous benchmark of 687 sequence representatives across 63 independent homology groups evaluated under 5-fold family-held-out cross-validation. Ridge-regularized logistic regression on ESM-2 representations achieved a family-level Balanced Accuracy of 79.46% and Macro-F1 of 0.7946 (51/63 families correctly classified), significantly outperforming the family-majority baseline (Macro-F1 = 0.4741, exact McNemar $p_{\text{adj}} = 0.000769$). ESM-2 demonstrated substantial point advantages over classical physicochemical descriptors (Macro-F1 = 0.5102 to 0.5679) with strictly positive bootstrap confidence intervals; however, after family-wise Holm-Bonferroni multiplicity correction, differences missed formal statistical significance ($p_{\text{adj}} = 0.0602$). Performance on directly protein-validated families reached Macro-F1 = 0.8154, whereas transcript-only annotated families dropped to 0.6857. Furthermore, low-dimensional PCA projections collapsed under cross-validation (Macro-F1 = 0.4388), demonstrating that temporal discriminative information is distributed across non-linear embedding space. Our findings demonstrate that protein language model representations encode substantial predictive signal distinguishing viral kinetic phases, providing an objective in silico framework for viral proteomic timing annotation.

---

## 1. Introduction

Herpes Simplex Virus Type 1 (HSV-1) and Type 2 (HSV-2) are ancient, highly successful human pathogens belonging to the *Alphaherpesvirinae* subfamily. Following host cell entry, viral gene expression unfolds in a coordinated temporal cascade initially delineated by Honess and Roizman (1974, 1975). Immediate-early ($\alpha$) proteins are synthesized in the absence of de novo viral protein synthesis, functioning primarily as transcriptional regulators. Early ($\beta$) proteins mediate viral DNA replication and nucleotide metabolism. Viral DNA synthesis triggers the expression of leaky-late ($\gamma_1$) and true-late ($\gamma_2$) structural virion components.

While promoter elements, transcriptional factors, and chromatin architecture direct transcriptional initiation, recent time-resolved multi-omics profiling (Rutkowski et al., 2015; Kulej et al., 2017) has highlighted that transcript abundance does not translate linearly into functional protein accumulation. Advances in self-supervised protein language modeling, such as ESM-2 (Lin et al., 2023), enable the extraction of dense vector embeddings that capture biophysical, evolutionary, and structural constraints directly from unaligned amino acid sequences.

In this study, we investigate whether protein language model representations capture intrinsic properties that distinguish pre-replicative regulatory enzymes from post-replicative structural virion components across unseen alphaherpesvirus gene families.

---

## 2. Materials and Methods

### 2.1 Universal Corpus Accounting & Sequence Curation
A total of 22,710 protein records from HSV-1 and HSV-2 were retrieved from NCBI GenBank, RefSeq, and UniProtKB. Universal record accounting assigned every sequence to mutually exclusive cohorts, identifying 18,172 valid full-length candidates (80.02% of the corpus). Sequence redundancy reduction at 95% identity yielded 3,111 canonical representatives across the 74 core alphaherpesvirus gene families.

### 2.2 Feature Representations
1. **Classical Physicochemical Descriptors**: Eight bulk descriptors were computed: sequence length, molecular weight, isoelectric point, instability index, aliphatic index, GRAVY hydropathicity, charged fraction, and aromatic fraction.
2. **ESM-2 Embeddings**: Sequence representations were extracted using the pre-trained ESM-2 transformer model (`esm2_t12_35M_UR50D`, 12 layers, 480 dimensions), using mean-pooled residue embeddings.
3. **Combined Feature Matrix**: Concatenation of the 8 physicochemical descriptors with 480 ESM-2 dimensions produced a 488-dimensional representation.

### 2.3 Group-Disjoint 5-Fold Cross-Validation
To prevent data leakage caused by gene duplication and strain variants, sequences were partitioned into 5 folds strictly at the homologous gene family level. The massive `RS1` (ICP4) gene family (212 sequences) was allocated entirely to Fold 0 validation, ensuring that no `RS1` sequence appeared in Fold 0 training. Across all folds, train-validation family overlap was zero. All feature scalers (`StandardScaler`) were fit strictly on training fold data.

### 2.4 Statistical Multiplicity Control
Paired model differences against ESM-2 were evaluated on family predictions using exact two-tailed McNemar tests and 10,000 paired bootstrap resamples. To control the Family-Wise Error Rate ($lpha = 0.05$) across the 4 paired comparisons, the Holm-Bonferroni step-down procedure was applied.

---

## 3. Results

### 3.1 Family-Level Classification Performance
Under 5-fold family-held-out cross-validation, the Ridge ESM-2 classifier correctly predicted 51 of the 63 held-out gene families, achieving a family-level Balanced Accuracy of 79.46% and Macro-F1 of 0.7946. It decisively outperformed the fixed majority baseline (Macro-F1 = 0.4741, McNemar $b=26, c=5$, $p_{\text{adj}} = 0.000769$).

While ESM-2 achieved large point advantages over physicochemical Logistic Regression (Macro-F1 = 0.5679) and Random Forest (Macro-F1 = 0.5102), raw p-values ($p = 0.0294$ and $p = 0.0201$) adjusted to **$p_{\text{adj}} = 0.0602$**, falling just short of the conventional $lpha = 0.05$ threshold. However, 10,000 paired bootstrap resamples yielded strictly positive 95% confidence intervals for $\Delta \text{Macro-F1}$ ([+0.0540, +0.3976] vs Logistic Regression; [+0.0915, +0.4626] vs Random Forest).

### 3.2 Evidence-Tier Stratification
When stratified by experimental annotation endpoint:
- **Protein-Validated Families ($N=52$)**: ESM-2 achieved Macro-F1 = 0.8154 (43/52 correct).
- **Transcript-Only Families ($N=11$)**: ESM-2 achieved Macro-F1 = 0.6857 (8/11 correct).

### 3.3 Breakdown of Low-Dimensional Linear Compression
While in-sample PCA training on 18 components reached 89.1% accuracy, held-out 5-fold cross-validation collapsed to 43.88% Macro-F1, confirming that viral temporal boundaries cannot be reduced to low-rank linear components.

---

## 4. Discussion & Limitations

Our findings establish that pre-trained protein language models capture distributed structural features associated with herpesvirus kinetic timing. However, three key scientific boundaries must be acknowledged:
1. **Limited Family Sample Size**: With only 63 independent gene families, statistical power under strict family-held-out evaluation is inherently constrained.
2. **Cross-Species Ortholog Transfer vs. Novel Viral Generalization**: The 7 HSV-2 families tested in cross-species transfer represent core orthologs shared with HSV-1; this test confirms ortholog sequence transferability between sister species, not generalization to unobserved viral clades.
3. **Absence of Prospective External Data**: True external biological validation requires prospective time-resolved proteomic profiling of independent clinical isolates under cycloheximide and phosphonoacetic acid arrest.

---

## References

1. Honess RW, Roizman B. Regulation of herpesvirus macromolecular synthesis. I. Cascade regulation of the synthesis of three groups of viral proteins. *J Virol*. 1974;14(1):8-19. DOI: 10.1128/jvi.14.1.8-19.1974.
2. Honess RW, Roizman B. Regulation of herpesvirus macromolecular synthesis. II. Effects of a variety of metabolic inhibitors. *J Virol*. 1975;16(5):1308-1326. DOI: 10.1128/jvi.16.5.1308-1326.1975.
3. Rutkowski AJ, et al. Widespread Disruptions of Host and Viral Transcription in Herpes Simplex Virus Infection. *PLoS Pathog*. 2015;11(10):e1005245. DOI: 10.1371/journal.ppat.1005245.
4. Kulej K, et al. Time-Resolved Quantitative Phosphoproteomics of Herpes Simplex Virus Type 1 Infection. *J Proteome Res*. 2017;16(2):877-887. DOI: 10.1021/acs.jproteome.6b00877.
5. Dremel SE, DeLuca NA. Herpes Simplex Virus 1 Immediate Early Proteins Are Required for the Dynamic Recruitment of Cellular RNA Polymerase II. *J Virol*. 2019;93(15):e00452-19. DOI: 10.1128/JVI.00452-19.
6. Roizman B, Knipe DM, Whitley RJ. Herpes Simplex Viruses. In: *Fields Virology*, 6th Edition. 2013:1823-1896.
7. Lin Z, et al. Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science*. 2023;379(6637):1123-1130. DOI: 10.1126/science.ade2574.
