# Step 16: Complete Publication-Quality Figure Captions

**Target Manuscript**: [`reports/16_illustrated_manuscript.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/16_illustrated_manuscript.md)  
**Deliverable File**: [`reports/16_figure_captions.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/16_figure_captions.md)  

---

### Figure 1. Comprehensive analytical workflow and multi-tiered cohort selection pipeline.
**(A) Universal Sequence Corpus Accounting**: Initial retrieval of 22,710 protein records from NCBI and UniProtKB for HSV-1 and HSV-2 was subjected to taxonomy curation and quality filtering, purging 1,043 non-target records (retroviral Gag/Pol/Env contaminants and host interactors) and quarantining 284 unverified or cross-species conflicting accessions. Full-length sequences $\ge 70$ amino acids formed a valid primary candidate cohort of 18,172 sequences (80.02% of total corpus), from which CD-HIT clustering at 95% sequence identity yielded 3,111 representative cluster seeds.  
**(B) Strict Benchmark Adjudication**: Sequence representatives were matched against authoritative biochemical literature (Honess & Roizman 1974, 1975; Rutkowski et al. 2015; Kulej et al. 2017). A strict direct-evidence benchmark of 687 sequences spanning 63 canonical gene families was established (excluding 633 homology-inferred HSV-2 sequences and 1,791 unresolved accessions).  
**(C) Feature Representation & Model Training**: 8 bulk physicochemical descriptors and 480-dimensional mean-pooled ESM-2 transformer embeddings (`esm2_t12_35M_UR50D`) were evaluated independently and in combination.  
**(D) Group-Disjoint 5-Fold Cross-Validation**: Sequences were partitioned strictly at the gene-family level (zero train-validation family overlap), with the large `RS1` (ICP4) family (212 sequences) isolated entirely within Fold 0 validation. Feature scalers were fit strictly within training folds.  
**(E) Statistical Inference & Multiplicity Control**: Out-of-fold family predictions were aggregated via mean probability and evaluated using exact McNemar tests and 10,000 paired family bootstrap resamples with Holm-Bonferroni Family-Wise Error Rate control ($lpha = 0.05$).

---

### Figure 2. Universal dataset composition, cohort hierarchy, and benchmark partitioning.
**(A) Universal Corpus Partitioning**: Stacked allocation of all 22,710 HSV protein records into mutually exclusive categories: valid primary candidates (18,172 records; 80.02%, comprising 3,111 unique 95% ID cluster seeds and 15,061 near-identical clinical strain variants), separate fragment cohort (3,211 records; 14.14%), purged non-targets (1,043 records; 4.59%), quarantined unverified constructs (254 records; 1.12%), and quarantined cross-species duplicate conflicts (30 records; 0.13%).  
**(B) Supervised Benchmark Partitioning**: Partition of the 3,111 representative sequence seeds into the strict direct-evidence benchmark (687 records across 63 gene families; 22.08%), homology-inferred candidates (633 records; 20.35%), and unresolved/accessory ORFs (1,791 records; 57.57%). Exact counts and percentages are displayed.

---

### Figure 3. Temporal replication timing benchmark composition and evidence stratification.
**(A) Gene Family Distribution by Kinetic Class**: Class breakdown across the 63 independent gene families: 23 pre-replicative regulatory families (3 immediate-early $lpha$, 20 early $eta$) and 40 post-replicative structural families (17 leaky-late $\gamma_1$, 23 true-late $\gamma_2$).  
**(B) Experimental Evidence Stratification**: Division of the benchmark into 52 directly protein-validated families (82.5%; synthesis timing confirmed by $^{35}	ext{S}$ metabolic radiolabeling or SILAC mass spectrometry under inhibitor arrest) and 11 transcript-only annotated families (17.5%; kinetics inferred exclusively from 4sU nascent RNA-seq).  
**(C) Sequence Representation and Family Size Asymmetry**: Sequence variant count per gene family in the direct benchmark, highlighting extreme asymmetry: the immediate-early transactivator `RS1` contributes 212 sequence variants (30.9% of the benchmark), followed by `UL23` (thymidine kinase, 39 variants), `UL19` (major capsid protein, 25 variants), and `RL2` (ICP0, 21 variants), compared to 20 singleton families.

---

### Figure 4. ESM-2 high-dimensional embedding geometry and low-rank generalization collapse.
**(A) Scree and Cumulative Explained Variance**: Eigenvalue decay across principal components of the 480-dimensional ESM-2 embedding space for 3,111 representative seeds. The first 18 components capture 90.0% variance and 30 components capture 95.1% variance. The Roy-Vetterli effective rank is 13.6 for standardized correlation and 5.9 for raw unstandardized covariance.  
**(B) Two-Dimensional Unsupervised PCA Projection**: Visualization of the 687 benchmark sequences in the first two principal component coordinates (PC1: 34.2% variance, PC2: 12.8% variance). Sequences are colored by true replication timing class (Pre-replicative: blue; Post-replicative: orange) and shaped by experimental evidence tier (Circles: Protein-validated; Triangles: Transcript-only). Unsupervised linear projection shows substantial overlap between classes, demonstrating that temporal kinetic phase cannot be separated linearly in low dimensions.  
**(C) Generalization Collapse of Low-Dimensional PCA Classifiers**: Out-of-fold family Macro-F1 across varying PCA component counts (5, 10, 18, 25, 30) evaluated under group-disjoint 5-fold cross-validation compared to in-sample training performance. While an 18-component PCA classifier achieves 89.1% training accuracy, its held-out family Macro-F1 collapses to 0.4388, demonstrating severe overfitting and proving that discriminative viral kinetic information is non-linear and distributed across high-dimensional coordinates.

---

### Figure 5. Supervised model comparison, paired bootstrap distributions, and multiplicity bounds.
**(A) Family-Level Performance across Evaluated Models**: Point-and-interval plot comparing Family Macro-F1 and Balanced Accuracy across five models evaluated under group-disjoint 5-fold cross-validation on 63 independent gene families: Fixed Majority Baseline (Macro-F1 = 0.4741), Physicochemical Random Forest (Macro-F1 = 0.5102), Physicochemical Logistic Regression (Macro-F1 = 0.5679), Combined Physicochemical + ESM-2 (Macro-F1 = 0.7644), and Ridge-Regularized ESM-2 (Macro-F1 = 0.7946, 51/63 correct). Error bars indicate 95% bootstrap confidence intervals across 10,000 resamples.  
**(B) Paired Bootstrap Difference Distributions ($\Delta 	ext{Macro-F1}$ vs ESM-2)**: Empirical bootstrap distributions of paired performance gains conferred by ESM-2 over comparator models. 95% bootstrap intervals are strictly positive: [+0.0915, +0.4626] vs Random Forest, [+0.0540, +0.3976] vs Logistic Regression, and [+0.1758, +0.4731] vs Majority.  
**(C) Statistical Multiplicity Bounds**: Exact two-tailed McNemar test results. Superiority over the majority baseline is statistically significant ($p = 0.000192$, Holm $p_{	ext{adj}} = 0.000769$). However, advantages over Physicochemical Random Forest ($p = 0.020062$) and Logistic Regression ($p = 0.029449$) yield Holm-adjusted $p_{	ext{adj}} = 0.060185 > 0.05$, missing formal significance under Family-Wise Error Rate control.

---

### Figure 6. Family-level classifier discrimination, confusion matrix, and error decomposition.
**(A) Family-Level Confusion Matrix**: Classification matrix for the Ridge ESM-2 model across all 63 independent gene families evaluated out-of-fold. True Positive (POST): 34; True Negative (PRE): 17; False Positive (PRE predicted as POST): 6; False Negative (POST predicted as PRE): 6. Pre-replicative recall is 73.9% (17/23) and post-replicative recall is 85.0% (34/40).  
**(B) Predicted Probability Distribution and Error Identification**: Out-of-fold predicted probability of pre-replicative timing ($ar{P}_{	ext{PRE}}$) for each of the 63 families, ranked by prediction confidence. True pre-replicative families are colored blue; true post-replicative families are colored orange. The horizontal dashed line indicates the decision threshold ($ar{P}_{	ext{PRE}} = 0.5$). Specific misclassified families are annotated: false positives (`UL2`, `UL3`, `UL4`, `UL24`, `UL35`, `UL45`) and false negatives (`UL14`, `UL21`, `UL50`, `US2`, `US3`, `US11`).

---

### Figure 7. Cross-validation stability, subgroup stratification, and evaluation unit discrepancy.
**(A) Fold-to-Fold Metric Stability**: Family-level Balanced Accuracy and Macro-F1 across the 5 group-disjoint cross-validation folds. Performance ranges from Macro-F1 = 0.636 (Fold 3) to 0.941 (Fold 0), demonstrating stability despite modest fold sample sizes ($N=12$ to 13 families per fold).  
**(B) Biological Evidence-Tier Stratification**: Model performance stratified by experimental annotation endpoint. ESM-2 achieves Macro-F1 = 0.8154 (43/52 correct) on directly protein-validated families ($N=52$), but decreases to Macro-F1 = 0.6857 (8/11 correct) on transcript-only annotated families ($N=11$).  
**(C) Leave-One-Family-Out (LOFO) Sensitivity**: Distribution of family Macro-F1 across 63 LOFO iterations (Mean = 0.7946, Range: [0.778, 0.811]), demonstrating that no single gene family drives the observed predictive performance.  
**(D) Evaluation Unit Inflation**: Comparison between sequence-weighted metrics (Accuracy = 91.12%, Macro-F1 = 0.9112) and gene family-level metrics (Accuracy = 80.95%, Macro-F1 = 0.7946). The sequence-level metric is inflated by +11.7% due to the 212 near-identical `RS1` variants in Fold 0, establishing that family-level aggregation is required for defensible biological evaluation.

---

### Figure 8 (Supplementary). Unsupervised clustering landscape and sister-species ortholog sequence transfer.
**(A) HDBSCAN Unsupervised Clustering**: Cluster size distribution and noise fraction across 3,111 representative sequence seeds. Clusters predominantly reconstruct known homologous gene families, demonstrating that unsupervised representations reflect sequence homology rather than novel functional discovery.  
**(B) Sister-Species Ortholog Transfer Diagnostic**: Out-of-fold classification accuracy of an HSV-1-trained ESM-2 model evaluated on 7 independent HSV-2 gene families (`RL2`, `RS1`, `UL19`, `UL23`, `UL27`, `US4`, `US6`). 7 of 7 families (100%) were correctly classified, establishing high sequence representation transferability between sister species.
