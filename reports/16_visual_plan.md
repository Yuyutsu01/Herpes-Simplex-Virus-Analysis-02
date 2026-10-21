# Step 16: Comprehensive Visual and Tabular Plan for the HSV-1/HSV-2 Research Paper

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
