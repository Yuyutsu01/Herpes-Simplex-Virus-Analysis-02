# Scientific Scope and Limitations

This document explicitly defines the biological scope, analytical boundaries, statistical limitations, and interpretation constraints of the HSV-1 and HSV-2 temporal proteomics study.

---

## 1. Biological and Experimental Scope

### 1.1 Predictive Modeling vs. Biological Causality
The machine learning models developed in this project operate strictly on primary amino acid sequences. They predict whether a viral protein belongs to the **PRE-replicative** (immediate-early [IE] / early [E]) or **POST-replicative** (leaky-late [L1] / true-late [L2]) phase of the viral replication cycle.

- **Non-Causal Representation**: High predictive weight on a specific amino acid sequence pattern or ESM-2 latent dimension reflects statistical association with viral functional niches (e.g., late structural capsid/envelope components versus early enzymatic/regulatory factors). It does not imply that amino acid composition directly drives transcription timing in the host cell.
- **Transcriptional vs. Translational Kinetics**: Viral temporal classes were historically defined through transcript accumulation under chemical inhibitors (e.g., cycloheximide for IE, phosphonoacetic acid / acyclovir for true late). While steady-state protein abundance generally correlates with viral mRNA cascade kinetics, post-translational modifications, viral proteolysis, and differential protein half-lives can cause divergence between transcript timing and functional protein persistence.

### 1.2 Target Organisms and Taxa
The analyses are calibrated specifically on Human alphaherpesvirus 1 (HSV-1, strain 17) and Human alphaherpesvirus 2 (HSV-2, strain HG52). Findings should not be extrapolated uncritically to betaherpesviruses (e.g., HCMV) or gammaherpesviruses (e.g., EBV, KSHV), which possess divergent genomic architectures, distinct regulatory cascades, and non-homologous accessory proteins.

---

## 2. Benchmark Design & Generalization Limits

### 2.1 Family-Disjoint Evaluation vs. Sequence Leakage
Herpesvirus proteomes contain paralogous gene families and multi-strain variants exhibiting high sequence identity. Evaluating predictive models with naive random splitting produces massive data leakage and inflated performance estimates.
- To counter this, all models were evaluated using **5-fold family-disjoint cross-validation**, where all members of a homology group (gene family) are held out together.
- Consequently, metrics reported at the **family level** (one vote per homology group) represent the true ability of models to generalize across homologous protein clusters.

### 2.2 Sister-Species Ortholog Transfer
The cross-species generalization experiment evaluated an HSV-1-trained model on 7 directly-adjudicated HSV-2 ortholog families:
- **Observed Result**: 7/7 (100%) family concordance.
- **Scientific Limitation**: This demonstrates robust sister-species transfer across highly conserved orthologous pairs. It **does not** constitute evidence that the model can predict replication timing for completely novel, uncharacterized viral protein families that lack any evolutionary kinship to known herpesvirus folds.

---

## 3. Statistical Multiplicity and Significance Nuances

### 3.1 Paired Model Comparisons (ESM-2 vs. Physicochemical LR)
In the 5-fold cross-validation benchmark on 63 gene families:
- **ESM-2 Logistic Regression**: 50/63 concordant families (Macro-F1 = 0.7946, Balanced Accuracy = 0.7946).
- **Physicochemical Logistic Regression**: 44/63 concordant families (Macro-F1 = 0.7107, Balanced Accuracy = 0.7107).
- **Discordant Pairs**: 9 families where ESM-2 was correct and PhysChem was incorrect; 3 families where PhysChem was correct and ESM-2 was incorrect (net gain of $+6$ families for ESM-2).

### 3.2 Unadjusted vs. Multiplicity-Corrected P-Values
- **Unadjusted McNemar Exact Test**: $p = 0.03125$ ($< 0.05$).
- **Holm-Bonferroni Multiplicity Adjustment**: Across the pre-specified family of 5 model comparisons (ESM-2 vs. Fixed Majority, Fold-Specific Majority, PhysChem LR, PhysChem RF, ESM-2 RF), the adjusted p-value becomes:
  $$p_{\text{adj}} = 0.0602$$
- **Rigorous Interpretation**:
  - The hypothesis of statistical superiority at the strict $\alpha = 0.05$ significance level after conservative family-wise error rate control is **not rejected** ($p_{\text{adj}} = 0.0602$).
  - The ESM-2 advantage represents a **strong empirical trend with a moderate effect size** (Cohen's $\kappa = 0.443$, net $+6$ concordant families), but must not be overstated as definitive mathematical proof in the manuscript.

---

## 4. Intrinsic Dimensionality of Representations

Analysis of the ESM-2 650M representation space across the herpesvirus proteome reveals substantial geometric concentration:
- **Ambient Dimension**: $D = 1280$.
- **Raw Unstandardized Covariance Effective Rank**: $R_{\text{eff}} = 5.93$.
- **Standardized Correlation Effective Rank**: $R_{\text{eff}} = 13.63$.
- **Interpretation**: Despite 1,280 available dimensions, viral protein sequences occupy a low-dimensional manifold. Top principal components correlate strongly with basic physicochemical markers (PC1 with molecular weight and length, $r = 0.78$; PC2 with isoelectric point and net charge, $r = -0.65$), explaining why shallow linear models can extract temporal signals without deep non-linear architectures.

---

## 5. Summary of Methodological Decisions

| Potential Vulnerability | Mitigation Implemented | Residual Limitation |
|---|---|---|
| Homology leakage inflating metrics | 5-fold family-disjoint grouped splitting | Family definitions rely on established viral annotations |
| Sequence redundancy skewing class metrics | Macro-averaged family-level evaluation | Benchmark limited to 63 experimentally verified families |
| Spurious p-values from multiple testing | Pre-registered paired tests with Holm correction | Moderate statistical power due to $N=63$ families |
| Black-box model over-interpretation | PCA projection and biophysical correlation analysis | Sequence features capture evolutionary niches, not transcriptional mechanisms |
