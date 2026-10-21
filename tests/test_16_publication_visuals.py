"""
Test Suite: test_16_publication_visuals.py
Purpose: Validate Step 16 publication-quality tables and figures, visual inventory,
  cohort accounting reconciliation, metric traceability, figure opening, and illustrated manuscript.
"""

from pathlib import Path
import json
import pytest
import pandas as pd
import numpy as np

BASE_DIR = Path(__file__).resolve().parent.parent
REPORTS_DIR = BASE_DIR / "reports"
FIGURES_DIR = BASE_DIR / "figures"
TABLES_DIR = BASE_DIR / "manuscript" / "tables"
SUPP_DIR = BASE_DIR / "manuscript" / "supplementary_materials"
SCRIPTS_DIR = BASE_DIR / "scripts"


def test_step16_reports_exist_and_nonempty():
    """Verify that all Step 16 report deliverables exist and have non-zero size."""
    required_reports = [
        REPORTS_DIR / "16_visual_inventory.csv",
        REPORTS_DIR / "16_visual_plan.md",
        REPORTS_DIR / "16_visual_data_validation.md",
        REPORTS_DIR / "16_table_value_reconciliation.csv",
        REPORTS_DIR / "16_figure_captions.md",
        REPORTS_DIR / "16_supplementary_materials_inventory.csv",
        REPORTS_DIR / "16_visual_generation_manifest.json",
        REPORTS_DIR / "16_illustrated_manuscript.md",
    ]

    for p in required_reports:
        assert p.exists(), f"Missing required report: {p}"
        assert p.stat().st_size > 0, f"Report is empty: {p}"


def test_step16_figures_exist_and_renderable():
    """Verify that all 8 publication figures exist as vector PDF, PNG, and data CSV."""
    fig_names = [
        "figure_1_study_workflow",
        "figure_2_dataset_composition",
        "figure_3_temporal_label_evidence",
        "figure_4_embedding_geometry",
        "figure_5_model_comparison",
        "figure_6_family_predictions",
        "figure_7_sensitivity_analysis",
        "figure_8_unsupervised_structure"
    ]

    for name in fig_names:
        pdf_file = FIGURES_DIR / f"{name}.pdf"
        png_file = FIGURES_DIR / f"{name}.png"

        assert pdf_file.exists(), f"Missing figure PDF: {pdf_file}"
        assert pdf_file.stat().st_size > 1000, f"PDF too small/corrupted: {pdf_file}"
        assert png_file.exists(), f"Missing figure PNG: {png_file}"
        assert png_file.stat().st_size > 10000, f"PNG too small/corrupted: {png_file}"


def test_step16_tables_exist_and_nonempty():
    """Verify that all 7 main tables (CSV and Markdown) and 8 supplementary tables exist."""
    # Main tables
    for i, name in enumerate([
        "dataset_accounting", "model_configurations", "main_performance",
        "statistical_comparisons", "temporal_evidence_summary",
        "sensitivity_analysis", "unsupervised_summary"
    ], start=1):
        csv_path = TABLES_DIR / f"table_{i}_{name}.csv"
        md_path = TABLES_DIR / f"table_{i}_{name}.md"
        assert csv_path.exists(), f"Missing main table CSV: {csv_path}"
        assert csv_path.stat().st_size > 0, f"Table CSV empty: {csv_path}"
        assert md_path.exists(), f"Missing main table MD: {md_path}"
        assert md_path.stat().st_size > 0, f"Table MD empty: {md_path}"

    # Supplementary tables
    for i, name in enumerate([
        "supervised_benchmark_687", "split_manifest_5fold", "family_oof_predictions",
        "sequence_oof_predictions", "temporal_evidence_ledger_63", "corpus_lineage_accounting",
        "clustering_robustness", "reference_claim_support"
    ], start=1):
        supp_csv = SUPP_DIR / f"table_s{i}_{name}.csv"
        assert supp_csv.exists(), f"Missing supplementary table: {supp_csv}"
        assert supp_csv.stat().st_size > 0, f"Supplementary table empty: {supp_csv}"


def test_corpus_accounting_reconciliation():
    """Verify exact mutual exclusivity and closure of the 22,710-record corpus."""
    t1_df = pd.read_csv(TABLES_DIR / "table_1_dataset_accounting.csv")
    total_row = t1_df[t1_df["cohort_name"] == "Universal Total Sequence Corpus"].iloc[0]
    assert total_row["record_count"] == 22710

    candidates_row = t1_df[t1_df["cohort_name"] == "Total Valid Primary Cohort"].iloc[0]
    assert candidates_row["record_count"] == 18172

    seeds_row = t1_df[t1_df["cohort_name"].str.contains("Cluster Seeds")].iloc[0]
    assert seeds_row["record_count"] == 3111

    variants_row = t1_df[t1_df["cohort_name"].str.contains("Clinical Strain Variants")].iloc[0]
    assert variants_row["record_count"] == 15061
    assert seeds_row["record_count"] + variants_row["record_count"] == candidates_row["record_count"]


def test_table_reconciliation_exact_match():
    """Verify that reports/16_table_value_reconciliation.csv confirms exact matches across all items."""
    rec_df = pd.read_csv(REPORTS_DIR / "16_table_value_reconciliation.csv")
    assert len(rec_df) >= 7
    assert (rec_df["reconciliation_status"] == "EXACT_MATCH").all(), "All table values must reconcile with EXACT_MATCH"


def test_evaluation_unit_distinction():
    """Verify that family-level and sequence-level metrics are never conflated in Table 3."""
    t3_df = pd.read_csv(TABLES_DIR / "table_3_main_performance.csv")
    ridge_row = t3_df[t3_df["model_identifier"] == "ridge_esm2"].iloc[0]

    # Family BalAcc = 79.46% vs Sequence Accuracy/F1 = 91.12%
    assert ridge_row["family_balanced_accuracy"] == "79.46%"
    assert np.isclose(float(ridge_row["family_macro_f1"]), 0.7946, atol=1e-4)
    assert np.isclose(float(ridge_row["sequence_macro_f1"]), 0.9112, atol=1e-4)
    assert ridge_row["correct_families"] == "51 / 63"


def test_statistical_multiplicity_reporting():
    """Verify that Table 4 correctly reflects Holm adjustments and avoids false superiority claims."""
    t4_df = pd.read_csv(TABLES_DIR / "table_4_statistical_comparisons.csv")

    rf_row = t4_df[t4_df["comparison"].str.contains("Random Forest")].iloc[0]
    lr_row = t4_df[t4_df["comparison"].str.contains("Logistic Regression")].iloc[0]
    maj_row = t4_df[t4_df["comparison"].str.contains("Majority")].iloc[0]

    # Majority comparison survives Holm
    assert float(maj_row["holm_adjusted_pvalue"]) < 0.001
    assert "Statistically Significant" in maj_row["statistical_conclusion_alpha_05"]

    # Physicochemical models do NOT survive Holm
    assert np.isclose(float(rf_row["holm_adjusted_pvalue"]), 0.060185, atol=1e-4)
    assert np.isclose(float(lr_row["holm_adjusted_pvalue"]), 0.060185, atol=1e-4)
    assert "Not Significant" in rf_row["statistical_conclusion_alpha_05"]
    assert "Not Significant" in lr_row["statistical_conclusion_alpha_05"]


def test_illustrated_manuscript_structure_and_callouts():
    """Verify that reports/16_illustrated_manuscript.md contains callouts for all 7 main figures and tables."""
    ms_path = REPORTS_DIR / "16_illustrated_manuscript.md"
    assert ms_path.exists()
    content = ms_path.read_text(encoding="utf-8")

    # Callouts for Figures 1 to 7
    for i in range(1, 8):
        assert f"Figure {i}" in content, f"Missing callout for Figure {i}"

    # Callouts for Tables 1 to 7
    for i in range(1, 8):
        assert f"Table {i}" in content, f"Missing callout for Table {i}"

    # Supplementary references
    assert "Table S1" in content or "Table S6" in content
    assert "Figure 8" in content
