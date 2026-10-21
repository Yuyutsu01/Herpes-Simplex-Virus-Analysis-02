"""
Script to populate manuscript/ directory with main_manuscript.md,
manuscript/figures, manuscript/tables, and manuscript/supplementary_materials.
Also generates manuscript/README.md.
"""
import shutil
from pathlib import Path

def setup_manuscript():
    base_dir = Path(".")
    manuscript_dir = base_dir / "manuscript"
    figures_dest = manuscript_dir / "figures"
    tables_dest = manuscript_dir / "tables"
    supp_dest = manuscript_dir / "supplementary_materials"

    figures_dest.mkdir(parents=True, exist_ok=True)
    tables_dest.mkdir(parents=True, exist_ok=True)
    supp_dest.mkdir(parents=True, exist_ok=True)

    # 1. Copy main manuscript
    src_manuscript = base_dir / "reports" / "16_illustrated_manuscript.md"
    dst_manuscript = manuscript_dir / "main_manuscript.md"
    shutil.copy2(src_manuscript, dst_manuscript)
    print(f"Copied {src_manuscript} -> {dst_manuscript}")

    # 2. Copy figures
    figures_src = base_dir / "figures"
    for f in figures_src.glob("figure_*.*"):
        shutil.copy2(f, figures_dest / f.name)
    print(f"Copied {len(list(figures_dest.glob('*')))} figures to {figures_dest}")

    # 3. Copy main tables (table_1 to table_7)
    tables_src = base_dir / "tables"
    for t in tables_src.glob("table_[0-9]*.*"):
        shutil.copy2(t, tables_dest / t.name)
    print(f"Copied {len(list(tables_dest.glob('*')))} main tables to {tables_dest}")

    # 4. Copy supplementary tables (table_s*)
    for s in tables_src.glob("table_s*.*"):
        shutil.copy2(s, supp_dest / s.name)
    print(f"Copied {len(list(supp_dest.glob('*')))} supplementary tables to {supp_dest}")

    # 5. Generate manuscript/README.md
    readme_content = """# Manuscript Submission Bundle

This directory contains the submission-ready manuscript, high-resolution figures, publication tables, and supplementary materials for:

**"Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes: Integrating Physicochemical Profiling, Protein Language Model Embeddings, and Unsupervised Learning"**

---

## 1. Directory Structure

| Path | Contents | Format |
|---|---|---|
| [main_manuscript.md](file:///manuscript/main_manuscript.md) | Full peer-reviewed manuscript text with integrated figures and tables | Markdown |
| [figures/](file:///manuscript/figures/) | Publication-ready figures (Figures 1-8) | Vector PDF, 300 DPI PNG, Data CSV |
| [tables/](file:///manuscript/tables/) | Main manuscript tables (Tables 1-7) | CSV, Markdown |
| [supplementary_materials/](file:///manuscript/supplementary_materials/) | Supplementary tables (Tables S1-S8) and extended evidence files | CSV |

---

## 2. Main Figures Summary

- [figure_1_study_workflow.pdf](file:///manuscript/figures/figure_1_study_workflow.pdf): End-to-end study workflow and multi-cohort validation architecture.
- [figure_2_dataset_composition.pdf](file:///manuscript/figures/figure_2_dataset_composition.pdf): Universal corpus accounting (22,710 records) and gene family distribution.
- [figure_3_temporal_label_evidence.pdf](file:///manuscript/figures/figure_3_temporal_label_evidence.pdf): Temporal kinetic evidence grading and 63-family benchmark class balance.
- [figure_4_embedding_geometry.pdf](file:///manuscript/figures/figure_4_embedding_geometry.pdf): ESM-2 embedding geometry, scree plots, effective rank, and biophysical alignments.
- [figure_5_model_comparison.pdf](file:///manuscript/figures/figure_5_model_comparison.pdf): Cross-validation ROC, PR, and Macro-F1 forest plots with 95% bootstrap CIs.
- [figure_6_family_predictions.pdf](file:///manuscript/figures/figure_6_family_predictions.pdf): Family-level confusion matrices and probability distribution profiles.
- [figure_7_sensitivity_analysis.pdf](file:///manuscript/figures/figure_7_sensitivity_analysis.pdf): Sensitivity to evidence grading, threshold variation, and dimensional ablations.
- [figure_8_unsupervised_structure.pdf](file:///manuscript/figures/figure_8_unsupervised_structure.pdf): Unsupervised clustering, UMAP/t-SNE manifold projections, and silhouette validation.

---

## 3. Main Tables Summary

- [table_1_dataset_accounting.csv](file:///manuscript/tables/table_1_dataset_accounting.csv): Universal sequence corpus breakdown and filtering triage.
- [table_2_model_configurations.csv](file:///manuscript/tables/table_2_model_configurations.csv): Model architectural hyperparameters and regularization settings.
- [table_3_main_performance.csv](file:///manuscript/tables/table_3_main_performance.csv): Comprehensive comparative model performance across 63 families.
- [table_4_statistical_comparisons.csv](file:///manuscript/tables/table_4_statistical_comparisons.csv): Paired McNemar hypothesis tests, odds ratios, and Holm-adjusted p-values.
- [table_5_temporal_evidence_summary.csv](file:///manuscript/tables/table_5_temporal_evidence_summary.csv): Kinetic evidence level distributions and literature concordance.
- [table_6_sensitivity_analysis.csv](file:///manuscript/tables/table_6_sensitivity_analysis.csv): Multi-endpoint sensitivity, fragment perturbation, and label flip tests.
- [table_7_unsupervised_summary.csv](file:///manuscript/tables/table_7_unsupervised_summary.csv): Latent dimensionality, effective rank calculations, and variance explained.

---

## 4. Supplementary Materials Summary

- `table_s1_supervised_benchmark_687.csv`: Complete 687-sequence benchmark ledger with accession, gene, and class.
- `table_s2_split_manifest_5fold.csv`: Strict family-disjoint 5-fold cross-validation partition assignments.
- `table_s3_family_oof_predictions.csv`: Out-of-fold predicted probabilities for all 63 gene families.
- `table_s4_sequence_oof_predictions.csv`: Out-of-fold sequence-level predictions.
- `table_s5_temporal_evidence_ledger_63.csv`: PubMed citations, quotes, and experimental assays for each family.
- `table_s6_corpus_lineage_accounting.csv`: Universal 22,710-record provenance ledger.
- `table_s7_clustering_robustness.csv`: Silhouette and Davies-Bouldin metrics across clustering algorithms.
- `table_s8_reference_claim_support.csv`: Traceability ledger linking manuscript claims to empirical code outputs.
"""
    with open(manuscript_dir / "README.md", "w", encoding="utf-8") as f:
        f.write(readme_content.strip() + "\n")
    print("Generated manuscript/README.md")

if __name__ == "__main__":
    setup_manuscript()
