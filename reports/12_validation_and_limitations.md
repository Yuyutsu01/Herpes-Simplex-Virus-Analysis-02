# Step 12: Validation, Limitations, and Epistemological Boundaries

**Deliverable**: [`reports/12_validation_and_limitations.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/12_validation_and_limitations.md)  
**Parent Report**: [`reports/12_representation_interpretation.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/12_representation_interpretation.md)  
**Date**: 2026-10-21  

---

## 1. Distinction of Findings by Epistemological Status

To maintain strict scientific honesty, all conclusions from Step 12 are explicitly classified into six formal epistemological categories:

### 1. Directly Observed Empirical Results
- The 480-dimensional ESM-2 embedding space exhibits extreme collinearity across the 687 viral sequences: the covariance matrix has condition number **3.50e+15**, and **18 principal components explain 95% of the total variance**.
- Linear model coefficients across 5 cross-validation folds show high overall directional correlation (mean pairwise Pearson $r = 0.7547$), but individual dimension weights fluctuate due to collinear redundancy.
- Only 480-dimensional mean-pooled representations exist in the frozen dataset; token-level embeddings and attention tensors are absent from saved files.

### 2. Recomputed Results from Existing Artifacts
- The out-of-fold performance of Ridge ESM-2 on the 63 direct families is confirmed as: Accuracy = 80.95% (51/63), Balanced Accuracy = 79.46%, Macro-F1 = 0.7946.
- Stratification confirms that performance is higher on the 52 protein-validated families (Macro-F1 = 0.8154) than on the 11 transcript-only families (Macro-F1 = 0.6857).

### 3. Results from Newly Trained Explanatory Models
- Family-level logistic regression controlling for $\log(\text{sequence\_length})$ demonstrates that neither Lysine ($p = 0.1232$) nor Glutamate ($p = 0.0662$) nor isoelectric point ($p = 0.1248$) has an independent association with replication class.
- PCA-compressed logistic regression utilizing only 18 components achieves 89.1% training accuracy, confirming that nearly all predictive information resides in a compact sub-space.

### 4. Exploratory Associations (Require Cautious Interpretation)
- Unadjusted Mann-Whitney tests suggest a weak trend for pre-replicative proteins to have slightly higher Lysine and Glutamate fractions and longer sequences.
- However, **none of these features survive Benjamini-Hochberg FDR correction ($q > 0.76$) or Holm FWER correction ($p_{\text{adj}} > 0.76$) across the 28 tested features**. They are exploratory trends that must not be cited as established biological markers.

### 5. Independently Validated Findings
- The family-held-out protocol ensures that ESM-2's predictive accuracy is evaluated strictly on unseen viral gene families.
- The 7/7 result on HSV-2 orthologs validates ortholog sequence transferability between closely related sister species for known genes, but does not validate generalization to unobserved viral gene families.

### 6. Mechanistic Hypotheses That Remain Unproven
- **Hypothesis: ESM-2 detects specific kinetic promoter or chromatin regulation motifs.**  
  *Status*: **UNPROVEN AND HIGHLY IMPROBABLE FROM PROTEIN SEQUENCE ALONE**. Temporal replication cascades in herpesviruses are controlled by genomic cis-regulatory elements (e.g. TAATGARAT motifs in IE promoters), viral transcription factors (VP16, ICP4, ICP27), and host transcriptional machinery, not by structural features of the translated polypeptides.
- **Hypothesis: Individual embedding dimensions correspond to modular functional domains.**  
  *Status*: **UNPROVEN**. Because dimensions are collinear and non-identifiable, treating individual dimensions as biological switches is mathematically unjustified without orthogonal mutagenesis or structural verification.

---

## 2. Limitations of Token-Level and Attention Analyses

1. **Absence of Token-Level Tensors**: The repository stores only `esm2_t12_representatives_embeddings.npy` (mean-pooled across sequence tokens excluding special tokens). Multi-head attention matrices were not cached during embedding generation.
2. **Inversion Impossibility**: Mean pooling is a many-to-one non-invertible projection:
   $$\mathbf{e} = \frac{1}{L} \sum_{i=1}^{L} \mathbf{h}_i$$
   It is mathematically impossible to reconstruct token attributions $\mathbf{h}_i$ from $\mathbf{e}$ without re-running the model forward pass.
3. **Descriptive Nature of Attention**: Even when available, transformer attention weights represent routing patterns in deep representations, not direct physical contact maps or causal determinants of biological function.

---

## 3. Methodological and Sample Limitations

- **Sample Size**: The direct evidence benchmark comprises 63 gene families. While this represents the vast majority of the annotated HSV proteome, $N = 63$ limits the statistical power of subgroup analyses and multi-parameter regression.
- **Transcript vs. Protein Divergence**: For 11 of the 63 families, timing labels are derived exclusively from transcript kinetics. RNA kinetics can precede protein accumulation, introducing potential label noise that attenuates classification accuracy.
- **Sister-Species Transfer Scope**: The HSV-1 to HSV-2 transfer benchmark evaluates 7 ortholog families. It cannot evaluate cross-species transfer to novel or distantly related viral families.
