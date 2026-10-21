# Journal Submission Package vs. Research Repository Manifest

This document specifies the precise boundary between the **Journal Submission Package** (the curated bundle of manuscript text, vector figures, primary tables, and supplementary materials submitted to an academic journal) and the **Full Research Repository** (the open-science GitHub reproducibility archive).

---

## 1. Submission Package Inventory (For Journal Editorial Office)

The following curated artifacts are designated for journal peer review, editorial production, and publication hosting:

### 1.1 Primary Manuscript Documents
- **Location**: `manuscript/`
  - [main_manuscript.md](file:///manuscript/main_manuscript.md): Full peer-reviewed text containing Title, Abstract, Introduction, Results, Discussion, Methods, Data Availability, and References.

### 1.2 Main Figures (High-Resolution Vector PDF & 300 DPI PNG)
- **Location**: `figures/` and mirrored in `manuscript/figures/`
  - [figure_1_study_workflow.pdf](file:///figures/figure_1_study_workflow.pdf) / `.png`: End-to-end analytical pipeline and multi-cohort validation schematic.
  - [figure_2_dataset_composition.pdf](file:///figures/figure_2_dataset_composition.pdf) / `.png`: Universal corpus accounting (22,710 records) and gene family distribution.
  - [figure_3_temporal_label_evidence.pdf](file:///figures/figure_3_temporal_label_evidence.pdf) / `.png`: Temporal kinetic evidence grading and 63-family benchmark class balance.
  - [figure_4_embedding_geometry.pdf](file:///figures/figure_4_embedding_geometry.pdf) / `.png`: ESM-2 embedding geometry, scree plots, effective rank, and biophysical alignments.
  - [figure_5_model_comparison.pdf](file:///figures/figure_5_model_comparison.pdf) / `.png`: 5-fold family-disjoint cross-validation ROC, PR, and Macro-F1 forest plots with 95% bootstrap CIs.
  - [figure_6_family_predictions.pdf](file:///figures/figure_6_family_predictions.pdf) / `.png`: Family-level confusion matrices and probability distribution profiles.
  - [figure_7_sensitivity_analysis.pdf](file:///figures/figure_7_sensitivity_analysis.pdf) / `.png`: Sensitivity to evidence grading, threshold variation, and dimensional ablations.
  - [figure_8_unsupervised_structure.pdf](file:///figures/figure_8_unsupervised_structure.pdf) / `.png`: Unsupervised clustering, UMAP/t-SNE manifold projections, and silhouette validation.

### 1.3 Main Tables
- **Location**: `tables/` and mirrored in `manuscript/tables/`
  - [table_1_dataset_accounting.csv](file:///manuscript/tables/table_1_dataset_accounting.csv): Universal sequence corpus breakdown and filtering triage.
  - [table_2_model_configurations.csv](file:///manuscript/tables/table_2_model_configurations.csv): Model architectural hyperparameters and regularization settings.
  - [table_3_main_performance.csv](file:///manuscript/tables/table_3_main_performance.csv): Comprehensive comparative model performance across 63 families.
  - [table_4_statistical_comparisons.csv](file:///manuscript/tables/table_4_statistical_comparisons.csv): Paired McNemar hypothesis tests, odds ratios, and Holm-adjusted p-values.
  - [table_5_temporal_evidence_summary.csv](file:///manuscript/tables/table_5_temporal_evidence_summary.csv): Kinetic evidence level distributions and literature concordance.
  - [table_6_sensitivity_analysis.csv](file:///manuscript/tables/table_6_sensitivity_analysis.csv): Multi-endpoint sensitivity, fragment perturbation, and label flip tests.
  - [table_7_unsupervised_summary.csv](file:///manuscript/tables/table_7_unsupervised_summary.csv): Latent dimensionality, effective rank calculations, and variance explained.

### 1.4 Supplementary Materials & Tables
- **Location**: `manuscript/supplementary_materials/` and `tables/`
  - [table_s1_supervised_benchmark_687.csv](file:///manuscript/supplementary_materials/table_s1_supervised_benchmark_687.csv): Complete 687-sequence benchmark ledger with accession, gene, and class.
  - [table_s2_split_manifest_5fold.csv](file:///manuscript/supplementary_materials/table_s2_split_manifest_5fold.csv): Strict family-disjoint 5-fold cross-validation partition assignments.
  - [table_s3_family_oof_predictions.csv](file:///manuscript/supplementary_materials/table_s3_family_oof_predictions.csv): Out-of-fold predicted probabilities for all 63 gene families.
  - [table_s4_sequence_oof_predictions.csv](file:///manuscript/supplementary_materials/table_s4_sequence_oof_predictions.csv): Sequence-level out-of-fold predictions.
  - [table_s5_temporal_evidence_ledger_63.csv](file:///manuscript/supplementary_materials/table_s5_temporal_evidence_ledger_63.csv): PubMed citations, quotes, and experimental assays for each family.
  - [table_s6_corpus_lineage_accounting.csv](file:///manuscript/supplementary_materials/table_s6_corpus_lineage_accounting.csv): Universal 22,710-record provenance ledger.
  - [table_s7_clustering_robustness.csv](file:///manuscript/supplementary_materials/table_s7_clustering_robustness.csv): Silhouette and Davies-Bouldin metrics across clustering algorithms.
  - [table_s8_reference_claim_support.csv](file:///manuscript/supplementary_materials/table_s8_reference_claim_support.csv): Traceability ledger linking manuscript claims to empirical code outputs.
  - [DATA_DICTIONARY.md](file:///docs/DATA_DICTIONARY.md): Complete schema documentation for supplementary data tables.
  - [SCIENTIFIC_SCOPE_AND_LIMITATIONS.md](file:///docs/SCIENTIFIC_SCOPE_AND_LIMITATIONS.md): Formal biological and statistical caveats.

---

## 2. Full Research Repository (For GitHub Open Science Hosting)

The public GitHub repository contains the complete computational environment, lineage history, raw data, intermediate caches, and test suites to guarantee bit-exact reproducibility. These items are **retained in the repository** but **excluded from the journal upload bundle**:

| Component | Path | Role in Repository | Why Excluded from Paper Upload |
|---|---|---|---|
| **Raw FASTA Sequences** | `data/raw/*.fasta` | Immutable input sequence data | Publicly available in UniProt; hashes provided |
| **ESM-2 Embedding Matrices** | `data/processed/*.npy` | High-dimensional tensor cache (~6MB) | Regenerable from sequences using documented script |
| **Comprehensive Test Suite** | `tests/` | 160+ automated scientific and structural tests | Internal verification infrastructure |
| **Phase Documentation** | `docs/phases/*.md` | Detailed 10-phase historical scientific logs | Condensed into manuscript Methods section |
| **Historical Audit Reports** | `reports/` | Intermediate development and peer review trails | Internal audit history; not part of final article |
| **CI / Environment Files** | `requirements.txt`, `.gitignore` | Virtual environment configuration | Standard open source tooling |

---

## 3. Data & Code Availability Statement (For Manuscript)

> **Data Availability**: The complete dataset of curated HSV-1 (strain 17) and HSV-2 (strain HG52) protein sequences, physicochemical features, temporal kinetic annotations, and cross-validation partition manifests are available in the public GitHub repository at [https://github.com/Yuyutsu01/Herpes-Simplex-Virus-Analysis-02](https://github.com/Yuyutsu01/Herpes-Simplex-Virus-Analysis-02) and archived on Zenodo (DOI: pending acceptance). Raw sequences were obtained from UniProt (HSV-1 reference proteome UP000009294; HSV-2 reference proteome UP000001874).
>
> **Code Availability**: All Python source code for data curation, feature extraction, ESM-2 representation inference, cross-validation modeling, statistical hypothesis testing, and publication figure generation is available under the MIT License in the project repository.
