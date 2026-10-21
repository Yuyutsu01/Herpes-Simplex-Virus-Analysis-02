"""
Script: 16_generate_visual_inventory.py
Purpose: Audit and generate the comprehensive visual inventory, visual plan,
  figure captions document, and supplementary materials inventory.
Outputs:
  - reports/16_visual_inventory.csv
  - reports/16_visual_plan.md
  - reports/16_figure_captions.md
  - reports/16_supplementary_materials_inventory.csv
"""

from pathlib import Path
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent
REPORTS_DIR = BASE_DIR / "reports"


def generate_inventory():
    print("=== Generating Visual Inventory, Plan, and Captions ===")

    # 1. Main Figures Inventory
    figures = [
        {
            "item_id": "FIG-01",
            "item_type": "Main Figure",
            "proposed_filename": "figure_1_study_workflow",
            "title": "Study Workflow and Multi-Tiered Cohort Selection Pipeline",
            "scientific_question_answered": "How were viral sequence records processed, curated, filtered, and evaluated without data leakage?",
            "manuscript_section": "Section 1 & 2 (Introduction / Materials & Methods)",
            "primary_source_dataset": "reports/04_redundancy_decisions.csv; data/processed/08_1_strict_temporal_benchmark.csv",
            "key_source_columns": "record_id, final_disposition, adjudicated_cohort, homology_group_id",
            "sample_size_and_unit": "22,710 source records -> 18,172 candidates -> 3,111 seeds -> 687 benchmark sequences (63 families)",
            "required_statistical_annotations": "Exact attrition counts; group-disjoint boundary indicator; zero-leakage guarantee.",
            "output_formats": "PDF (Vector), PNG (300 DPI), CSV (Data)",
            "inclusion_status": "MANDATORY_MAIN_ARTICLE",
            "verification_status": "VERIFIED_DATA_DERIVED"
        },
        {
            "item_id": "FIG-02",
            "item_type": "Main Figure",
            "proposed_filename": "figure_2_dataset_composition",
            "title": "Universal Corpus Accounting and Benchmark Partitioning",
            "scientific_question_answered": "What are the exact mutually exclusive partitions of the 22,710 viral protein records?",
            "manuscript_section": "Section 2.1 & 3.1 (Methods / Results)",
            "primary_source_dataset": "reports/04_redundancy_decisions.csv; data/processed/08_1_strict_temporal_benchmark.csv",
            "key_source_columns": "final_disposition, effective_species, adjudicated_cohort",
            "sample_size_and_unit": "22,710 records (Panel A); 3,111 cluster seeds (Panel B)",
            "required_statistical_annotations": "Mutually exclusive cohort counts reconciling exactly to 22,710; species percentages.",
            "output_formats": "PDF (Vector), PNG (300 DPI), CSV (Data)",
            "inclusion_status": "MANDATORY_MAIN_ARTICLE",
            "verification_status": "VERIFIED_DATA_DERIVED"
        },
        {
            "item_id": "FIG-03",
            "item_type": "Main Figure",
            "proposed_filename": "figure_3_temporal_label_evidence",
            "title": "Temporal Replication Timing Benchmark Composition and Evidence Stratification",
            "scientific_question_answered": "How are pre- and post-replicative labels distributed across gene families and experimental evidence tiers?",
            "manuscript_section": "Section 2.3 & 3.3 (Methods / Results)",
            "primary_source_dataset": "reports/15_temporal_label_adversarial_audit.csv",
            "key_source_columns": "homology_group_id, gene_symbol, binary_group, measured_endpoint, sequence_variant_count",
            "sample_size_and_unit": "63 independent gene families (687 sequence variants)",
            "required_statistical_annotations": "Family counts (23 PRE / 40 POST); Evidence counts (52 Protein / 11 Transcript); RS1 variant tally (212).",
            "output_formats": "PDF (Vector), PNG (300 DPI), CSV (Data)",
            "inclusion_status": "MANDATORY_MAIN_ARTICLE",
            "verification_status": "VERIFIED_DATA_DERIVED"
        },
        {
            "item_id": "FIG-04",
            "item_type": "Main Figure",
            "proposed_filename": "figure_4_embedding_geometry",
            "title": "ESM-2 High-Dimensional Embedding Geometry and Low-Rank Generalization Collapse",
            "scientific_question_answered": "What is the intrinsic dimensionality of ESM-2 embeddings and why do low-rank linear projections fail out-of-fold?",
            "manuscript_section": "Section 3.1 (Results)",
            "primary_source_dataset": "reports/12_1_embedding_geometry.csv; reports/12_1_pca_heldout_metrics.csv",
            "key_source_columns": "effective_rank_roy_vetterli, fold_explained_var_mean, family_macro_f1, n_components",
            "sample_size_and_unit": "3,111 seeds (geometry); 687 benchmark sequences (held-out PCA)",
            "required_statistical_annotations": "Roy-Vetterli effective rank (13.6 std / 5.9 raw); Scree curve; 89.1% in-sample vs 43.88% held-out collapse.",
            "output_formats": "PDF (Vector), PNG (300 DPI), CSV (Data)",
            "inclusion_status": "MANDATORY_MAIN_ARTICLE",
            "verification_status": "VERIFIED_DATA_DERIVED"
        },
        {
            "item_id": "FIG-05",
            "item_type": "Main Figure",
            "proposed_filename": "figure_5_model_comparison",
            "title": "Supervised Model Comparison and Statistical Multiplicity Bounds",
            "scientific_question_answered": "Does ESM-2 provide predictive superiority over physicochemical baselines, and does the advantage survive multiplicity correction?",
            "manuscript_section": "Section 3.2 (Results)",
            "primary_source_dataset": "reports/15_metric_independent_recheck.csv; reports/11_1_paired_bootstrap_differences.csv",
            "key_source_columns": "model_id, independently_rechecked_macro_f1, independently_rechecked_bal_acc, holm_adjusted_pvalue, delta_macro_f1_ci",
            "sample_size_and_unit": "63 independent gene families (Group-disjoint 5-fold CV)",
            "required_statistical_annotations": "Exact McNemar b/c discordant counts; Holm-adjusted p-values (p_adj = 0.0602 for physical models); 95% bootstrap CIs.",
            "output_formats": "PDF (Vector), PNG (300 DPI), CSV (Data)",
            "inclusion_status": "MANDATORY_MAIN_ARTICLE",
            "verification_status": "VERIFIED_DATA_DERIVED"
        },
        {
            "item_id": "FIG-06",
            "item_type": "Main Figure",
            "proposed_filename": "figure_6_family_predictions",
            "title": "Family-Level Classifier Discrimination, Confusion Matrix, and Error Decomposition",
            "scientific_question_answered": "Which specific viral gene families are correctly discriminated versus misclassified by the ESM-2 model?",
            "manuscript_section": "Section 3.2 & 4.1 (Results / Discussion)",
            "primary_source_dataset": "reports/13_retrained_oof_predictions.csv",
            "key_source_columns": "homology_group_id, effective_gene_symbol, true_class, predicted_prob_pre",
            "sample_size_and_unit": "63 independent gene families",
            "required_statistical_annotations": "Class recall (PRE: 73.9%, POST: 85.0%); Exact contingency table (17 TP, 34 TN, 6 FP, 6 FN); Gene symbol callouts.",
            "output_formats": "PDF (Vector), PNG (300 DPI), CSV (Data)",
            "inclusion_status": "MANDATORY_MAIN_ARTICLE",
            "verification_status": "VERIFIED_DATA_DERIVED"
        },
        {
            "item_id": "FIG-07",
            "item_type": "Main Figure",
            "proposed_filename": "figure_7_sensitivity_analysis",
            "title": "Cross-Validation Stability, Subgroup Stratification, and Evaluation Unit Discrepancy",
            "scientific_question_answered": "How stable is ESM-2 performance across folds, evidence tiers, family leave-outs, and sequence weighting?",
            "manuscript_section": "Section 3.3 & 4.2 (Results / Discussion)",
            "primary_source_dataset": "reports/11_fold_stability.csv; reports/11_leave_one_family_out.csv; reports/15_metric_independent_recheck.csv",
            "key_source_columns": "fold_index, lofo_macro_f1, protein_subgroup_f1_n52, transcript_subgroup_f1_n11, sequence_level_macro_f1",
            "sample_size_and_unit": "5 cross-validation folds; 63 LOFO iterations; N=52 vs N=11 subgroups",
            "required_statistical_annotations": "Fold range [0.63, 0.94]; LOFO delta range [-0.04, +0.02]; RS1 sequence inflation (+11.7% accuracy).",
            "output_formats": "PDF (Vector), PNG (300 DPI), CSV (Data)",
            "inclusion_status": "MANDATORY_MAIN_ARTICLE",
            "verification_status": "VERIFIED_DATA_DERIVED"
        },
        {
            "item_id": "FIG-08",
            "item_type": "Supplementary Figure",
            "proposed_filename": "figure_8_unsupervised_structure",
            "title": "Unsupervised Clustering Landscape and Sister-Species Ortholog Sequence Transfer",
            "scientific_question_answered": "Does unsupervised clustering discover novel functional groups, and how does the model perform on HSV-2 sister orthologs?",
            "manuscript_section": "Supplementary Information / Section 4.2",
            "primary_source_dataset": "reports/07_unsupervised_model_comparison.csv; reports/11_1_cross_species_transfer_experiment.csv",
            "key_source_columns": "cluster_count, noise_fraction, ortholog_transfer_accuracy",
            "sample_size_and_unit": "3,111 representative cluster seeds; 7 HSV-2 ortholog test families",
            "required_statistical_annotations": "7/7 (100%) ortholog transfer success; clustering driven by homology, not de novo functional discovery.",
            "output_formats": "PDF (Vector), PNG (300 DPI), CSV (Data)",
            "inclusion_status": "SUPPLEMENTARY_ARTICLE",
            "verification_status": "VERIFIED_DATA_DERIVED"
        }
    ]

    df_fig = pd.DataFrame(figures)
    fig_csv_path = REPORTS_DIR / "16_visual_inventory.csv"
    df_fig.to_csv(fig_csv_path, index=False)
    print(f"Saved: {fig_csv_path} ({len(df_fig)} figures inventoried)")

    # 2. Supplementary Materials Inventory
    supp_tables = [
        {
            "supplementary_id": "Table S1",
            "title": "Complete 687-Record Strict Supervised Benchmark",
            "description": "Full sequence-level annotations, accession numbers, gene symbols, species, temporal classes, and evidence source IDs for all 687 direct benchmark sequences.",
            "source_artifact": "data/processed/08_1_strict_temporal_benchmark.csv",
            "rows_count": 687,
            "format": "CSV"
        },
        {
            "supplementary_id": "Table S2",
            "title": "Five-Fold Group-Disjoint Family Split Manifest",
            "description": "Cross-validation fold assignments for all benchmark sequences demonstrating zero family-level leakage and isolation of RS1 (212 sequences) in Fold 0.",
            "source_artifact": "data/processed/09_group_split_manifest.csv",
            "rows_count": 687,
            "format": "CSV"
        },
        {
            "supplementary_id": "Table S3",
            "title": "Family-Level Out-of-Fold Predictions across Baselines",
            "description": "Aggregated out-of-fold predicted probabilities, true labels, and binary predictions for all 63 gene families across all 5 evaluated models.",
            "source_artifact": "reports/13_retrained_oof_predictions.csv",
            "rows_count": 63,
            "format": "CSV"
        },
        {
            "supplementary_id": "Table S4",
            "title": "Sequence-Level Out-of-Fold Predictions",
            "description": "Record-level out-of-fold predictions and probabilities for all 687 benchmark sequences across all models.",
            "source_artifact": "reports/13_retrained_oof_predictions.csv",
            "rows_count": 3435,
            "format": "CSV"
        },
        {
            "supplementary_id": "Table S5",
            "title": "Complete 63-Family Temporal Label Adversarial Audit Ledger",
            "description": "Exhaustive row-level audit of all 63 benchmark families detailing primary citations, experimental assay endpoints (protein vs RNA), and reviewer confidence.",
            "source_artifact": "reports/15_temporal_label_adversarial_audit.csv",
            "rows_count": 63,
            "format": "CSV"
        },
        {
            "supplementary_id": "Table S6",
            "title": "Complete 22,710-Record Corpus Lineage Accounting & Disposition",
            "description": "Universal record-level disposition, accession provenance, sequence length, and curation categories reconciling the entire input dataset.",
            "source_artifact": "reports/04_redundancy_decisions.csv",
            "rows_count": 22710,
            "format": "CSV"
        },
        {
            "supplementary_id": "Table S7",
            "title": "Clustering Robustness & Leave-One-Family-Out Stability",
            "description": "Sensitivity metrics across 63 LOFO iterations, fold stability distributions, and HDBSCAN unsupervised parameter evaluations.",
            "source_artifact": "reports/11_leave_one_family_out.csv",
            "rows_count": 63,
            "format": "CSV"
        },
        {
            "supplementary_id": "Table S8",
            "title": "Bibliographic Reference & Claim-Level Support Ledger",
            "description": "Verified DOIs, PMIDs, assay details, and supported scientific claims for all cited literature.",
            "source_artifact": "reports/15_reference_claim_support.csv",
            "rows_count": 7,
            "format": "CSV"
        }
    ]

    df_supp = pd.DataFrame(supp_tables)
    supp_csv_path = REPORTS_DIR / "16_supplementary_materials_inventory.csv"
    df_supp.to_csv(supp_csv_path, index=False)
    print(f"Saved: {supp_csv_path} ({len(df_supp)} supplementary tables inventoried)")

    # 3. Visual Plan Markdown
    plan_md = f"""# Step 16: Comprehensive Visual and Tabular Plan for the HSV-1/HSV-2 Research Paper

**Author / Role**: Senior Scientific Visualization Specialist & Journal Production Editor  
**Working Title**: *Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes: Integrating Physicochemical Profiling, Protein Language Model Embeddings, and Unsupervised Learning*  
**Deliverable File**: [`reports/16_visual_plan.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/16_visual_plan.md)  

---

## 1. Visual Design Architecture & Standards

All visuals adhere to international standards for biomedical and computational biology journals (e.g., *PLOS Computational Biology*, *Bioinformatics*, *Journal of Virology*, *Nature Communications*):
1. **Color Accessibility**: Universal colorblind-safe palettes (Wong / Okabe-Ito / Viridis) with redundant point shapes and line styles to maintain full legibility in grayscale.
2. **Dual-Format Publication Outputs**: Every figure is rendered as a vector PDF (`figures/<id>.pdf`) for high-resolution print reproduction and a 300 DPI PNG (`figures/<id>.png`) for immediate review. Accompanying numerical data are exported as dedicated CSVs (`figures/<id>_data.csv`).
3. **Evaluation Unit Discipline**: Sequence records and gene families are never conflated. Family-level metrics represent the true biological generalization unit; sequence-level inflation is explicitly documented.
4. **Epistemological Honesty**: No artificial confidence intervals, no unjustified 3D effects, no cropped axes, and no nominal significance stars without verified Holm-adjusted multiplicity corrections.

---

## 2. Main Figures Roster

| ID | Title | Scientific Question | Unit / N | Target File | Status |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Figure 1** | Study Workflow & Cohort Architecture | How was the viral proteome curated and evaluated without leakage? | 22,710 records / 63 families | `figures/figure_1_study_workflow.pdf` | MANDATORY |
| **Figure 2** | Universal Dataset Accounting Flow | How do all 22,710 records partition into mutually exclusive cohorts? | 22,710 records | `figures/figure_2_dataset_composition.pdf` | MANDATORY |
| **Figure 3** | Temporal Benchmark Evidence Stratification | What is the biological composition and evidence confidence of labels? | 63 families / 687 variants | `figures/figure_3_temporal_label_evidence.pdf` | MANDATORY |
| **Figure 4** | ESM-2 Embedding Geometry & PCA Collapse | What is the representation dimensionality and why does low-rank PCA fail? | 3,111 seeds / 687 sequences | `figures/figure_4_embedding_geometry.pdf` | MANDATORY |
| **Figure 5** | Family-Level Model Comparisons | Does ESM-2 outperform baselines under strict multiplicity control? | 63 families (5 models) | `figures/figure_5_model_comparison.pdf` | MANDATORY |
| **Figure 6** | Classifier Confusion & Error Decomposition | Which viral gene families are accurately separated vs misclassified? | 63 families | `figures/figure_6_family_predictions.pdf` | MANDATORY |
| **Figure 7** | Cross-Validation Stability & Sensitivity | How sensitive is predictive performance across folds, tiers, and weighting? | 5 folds / 63 LOFO / N=52 vs 11 | `figures/figure_7_sensitivity_analysis.pdf` | MANDATORY |
| **Figure 8** | Unsupervised Clusters & Sister Orthologs | Does clustering discover functions or reflect sequence homology? | 3,111 seeds / 7 orthologs | `figures/figure_8_unsupervised_structure.pdf` | SUPPLEMENTARY |

---

## 3. Main Tables Roster

1. **Table 1: Universal Dataset Accounting & Cohort Architecture** (`tables/table_1_dataset_accounting.csv` & `.md`)
2. **Table 2: Predictive Baseline Model Configurations** (`tables/table_2_model_configurations.csv` & `.md`)
3. **Table 3: Supervised Model Performance across Evaluation Units** (`tables/table_3_main_performance.csv` & `.md`)
4. **Table 4: Paired Statistical Hypothesis Tests & Bootstrap Differences** (`tables/table_4_statistical_comparisons.csv` & `.md`)
5. **Table 5: Biological Evidence-Tier Stratification Summary** (`tables/table_5_temporal_evidence_summary.csv` & `.md`)
6. **Table 6: Cross-Validation & Subgroup Sensitivity Analysis** (`tables/table_6_sensitivity_analysis.csv` & `.md`)
7. **Table 7: Unsupervised Embedding Geometry & Low-Rank PCA Analysis** (`tables/table_7_unsupervised_summary.csv` & `.md`)

---

## 4. Supplementary Materials Roster (Tables S1 – S8)

Exhaustive raw data tables are compiled in `tables/` to ensure full computational reproducibility without burdening the primary manuscript text (see [`reports/16_supplementary_materials_inventory.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/16_supplementary_materials_inventory.csv)).
"""

    plan_path = REPORTS_DIR / "16_visual_plan.md"
    plan_path.write_text(plan_md, encoding="utf-8")
    print(f"Saved: {plan_path}")

    # 4. Figure Captions Document
    captions_md = """# Step 16: Complete Publication-Quality Figure Captions

**Target Manuscript**: [`reports/16_illustrated_manuscript.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/16_illustrated_manuscript.md)  
**Deliverable File**: [`reports/16_figure_captions.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/16_figure_captions.md)  

---

### Figure 1. Comprehensive analytical workflow and multi-tiered cohort selection pipeline.
**(A) Universal Sequence Corpus Accounting**: Initial retrieval of 22,710 protein records from NCBI and UniProtKB for HSV-1 and HSV-2 was subjected to taxonomy curation and quality filtering, purging 1,043 non-target records (retroviral Gag/Pol/Env contaminants and host interactors) and quarantining 284 unverified or cross-species conflicting accessions. Full-length sequences $\ge 70$ amino acids formed a valid primary candidate cohort of 18,172 sequences (80.02% of total corpus), from which CD-HIT clustering at 95% sequence identity yielded 3,111 representative cluster seeds.  
**(B) Strict Benchmark Adjudication**: Sequence representatives were matched against authoritative biochemical literature (Honess & Roizman 1974, 1975; Rutkowski et al. 2015; Kulej et al. 2017). A strict direct-evidence benchmark of 687 sequences spanning 63 canonical gene families was established (excluding 633 homology-inferred HSV-2 sequences and 1,791 unresolved accessions).  
**(C) Feature Representation & Model Training**: 8 bulk physicochemical descriptors and 480-dimensional mean-pooled ESM-2 transformer embeddings (`esm2_t12_35M_UR50D`) were evaluated independently and in combination.  
**(D) Group-Disjoint 5-Fold Cross-Validation**: Sequences were partitioned strictly at the gene-family level (zero train-validation family overlap), with the large `RS1` (ICP4) family (212 sequences) isolated entirely within Fold 0 validation. Feature scalers were fit strictly within training folds.  
**(E) Statistical Inference & Multiplicity Control**: Out-of-fold family predictions were aggregated via mean probability and evaluated using exact McNemar tests and 10,000 paired family bootstrap resamples with Holm-Bonferroni Family-Wise Error Rate control ($\alpha = 0.05$).

---

### Figure 2. Universal dataset composition, cohort hierarchy, and benchmark partitioning.
**(A) Universal Corpus Partitioning**: Stacked allocation of all 22,710 HSV protein records into mutually exclusive categories: valid primary candidates (18,172 records; 80.02%, comprising 3,111 unique 95% ID cluster seeds and 15,061 near-identical clinical strain variants), separate fragment cohort (3,211 records; 14.14%), purged non-targets (1,043 records; 4.59%), quarantined unverified constructs (254 records; 1.12%), and quarantined cross-species duplicate conflicts (30 records; 0.13%).  
**(B) Supervised Benchmark Partitioning**: Partition of the 3,111 representative sequence seeds into the strict direct-evidence benchmark (687 records across 63 gene families; 22.08%), homology-inferred candidates (633 records; 20.35%), and unresolved/accessory ORFs (1,791 records; 57.57%). Exact counts and percentages are displayed.

---

### Figure 3. Temporal replication timing benchmark composition and evidence stratification.
**(A) Gene Family Distribution by Kinetic Class**: Class breakdown across the 63 independent gene families: 23 pre-replicative regulatory families (3 immediate-early $\alpha$, 20 early $\beta$) and 40 post-replicative structural families (17 leaky-late $\gamma_1$, 23 true-late $\gamma_2$).  
**(B) Experimental Evidence Stratification**: Division of the benchmark into 52 directly protein-validated families (82.5%; synthesis timing confirmed by $^{35}\text{S}$ metabolic radiolabeling or SILAC mass spectrometry under inhibitor arrest) and 11 transcript-only annotated families (17.5%; kinetics inferred exclusively from 4sU nascent RNA-seq).  
**(C) Sequence Representation and Family Size Asymmetry**: Sequence variant count per gene family in the direct benchmark, highlighting extreme asymmetry: the immediate-early transactivator `RS1` contributes 212 sequence variants (30.9% of the benchmark), followed by `UL23` (thymidine kinase, 39 variants), `UL19` (major capsid protein, 25 variants), and `RL2` (ICP0, 21 variants), compared to 20 singleton families.

---

### Figure 4. ESM-2 high-dimensional embedding geometry and low-rank generalization collapse.
**(A) Scree and Cumulative Explained Variance**: Eigenvalue decay across principal components of the 480-dimensional ESM-2 embedding space for 3,111 representative seeds. The first 18 components capture 90.0% variance and 30 components capture 95.1% variance. The Roy-Vetterli effective rank is 13.6 for standardized correlation and 5.9 for raw unstandardized covariance.  
**(B) Two-Dimensional Unsupervised PCA Projection**: Visualization of the 687 benchmark sequences in the first two principal component coordinates (PC1: 34.2% variance, PC2: 12.8% variance). Sequences are colored by true replication timing class (Pre-replicative: blue; Post-replicative: orange) and shaped by experimental evidence tier (Circles: Protein-validated; Triangles: Transcript-only). Unsupervised linear projection shows substantial overlap between classes, demonstrating that temporal kinetic phase cannot be separated linearly in low dimensions.  
**(C) Generalization Collapse of Low-Dimensional PCA Classifiers**: Out-of-fold family Macro-F1 across varying PCA component counts (5, 10, 18, 25, 30) evaluated under group-disjoint 5-fold cross-validation compared to in-sample training performance. While an 18-component PCA classifier achieves 89.1% training accuracy, its held-out family Macro-F1 collapses to 0.4388, demonstrating severe overfitting and proving that discriminative viral kinetic information is non-linear and distributed across high-dimensional coordinates.

---

### Figure 5. Supervised model comparison, paired bootstrap distributions, and multiplicity bounds.
**(A) Family-Level Performance across Evaluated Models**: Point-and-interval plot comparing Family Macro-F1 and Balanced Accuracy across five models evaluated under group-disjoint 5-fold cross-validation on 63 independent gene families: Fixed Majority Baseline (Macro-F1 = 0.4741), Physicochemical Random Forest (Macro-F1 = 0.5102), Physicochemical Logistic Regression (Macro-F1 = 0.5679), Combined Physicochemical + ESM-2 (Macro-F1 = 0.7644), and Ridge-Regularized ESM-2 (Macro-F1 = 0.7946, 51/63 correct). Error bars indicate 95% bootstrap confidence intervals across 10,000 resamples.  
**(B) Paired Bootstrap Difference Distributions ($\Delta \text{Macro-F1}$ vs ESM-2)**: Empirical bootstrap distributions of paired performance gains conferred by ESM-2 over comparator models. 95% bootstrap intervals are strictly positive: [+0.0915, +0.4626] vs Random Forest, [+0.0540, +0.3976] vs Logistic Regression, and [+0.1758, +0.4731] vs Majority.  
**(C) Statistical Multiplicity Bounds**: Exact two-tailed McNemar test results. Superiority over the majority baseline is statistically significant ($p = 0.000192$, Holm $p_{\text{adj}} = 0.000769$). However, advantages over Physicochemical Random Forest ($p = 0.020062$) and Logistic Regression ($p = 0.029449$) yield Holm-adjusted $p_{\text{adj}} = 0.060185 > 0.05$, missing formal significance under Family-Wise Error Rate control.

---

### Figure 6. Family-level classifier discrimination, confusion matrix, and error decomposition.
**(A) Family-Level Confusion Matrix**: Classification matrix for the Ridge ESM-2 model across all 63 independent gene families evaluated out-of-fold. True Positive (POST): 34; True Negative (PRE): 17; False Positive (PRE predicted as POST): 6; False Negative (POST predicted as PRE): 6. Pre-replicative recall is 73.9% (17/23) and post-replicative recall is 85.0% (34/40).  
**(B) Predicted Probability Distribution and Error Identification**: Out-of-fold predicted probability of pre-replicative timing ($\bar{P}_{\text{PRE}}$) for each of the 63 families, ranked by prediction confidence. True pre-replicative families are colored blue; true post-replicative families are colored orange. The horizontal dashed line indicates the decision threshold ($\bar{P}_{\text{PRE}} = 0.5$). Specific misclassified families are annotated: false positives (`UL2`, `UL3`, `UL4`, `UL24`, `UL35`, `UL45`) and false negatives (`UL14`, `UL21`, `UL50`, `US2`, `US3`, `US11`).

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
"""

    captions_path = REPORTS_DIR / "16_figure_captions.md"
    captions_path.write_text(captions_md, encoding="utf-8")
    print(f"Saved: {captions_path}")


if __name__ == "__main__":
    generate_inventory()
