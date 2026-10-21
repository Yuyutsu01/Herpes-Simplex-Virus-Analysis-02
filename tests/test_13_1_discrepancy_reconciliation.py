"""
Unit and integration tests for Step 13.1: Resolution of Benchmark Contradictions
and Independent Result Verification.

Verifies:
1. Deliverable presence, non-emptiness, and format integrity.
2. Record and fold reconciliation (212 RS1 records in Fold 0; Fold 3 Val = 89 seqs; zero leakage).
3. Model configuration audit (C=1.0 confirmed across all linear models).
4. Feature dimension audit (8 bulk descriptors, 480 ESM-2, 488 combined; not 28 / 508).
5. Exact prediction reproduction (0 mismatches, max diff < 1e-12 vs Step 10 and Step 13).
6. Family-level metrics and subgroup evaluations (N=52 protein vs N=11 transcript).
7. Statistical inference and Holm-Bonferroni correction (p_adj = 0.0602 for physicochemical models).
8. Report consistency with underlying CSV deliverables.
"""

from pathlib import Path
import json
import numpy as np
import pandas as pd
import pytest

BASE_DIR = Path(__file__).resolve().parent.parent
REPORTS_DIR = BASE_DIR / "reports"
DATA_DIR = BASE_DIR / "data" / "processed"

RECON_CSV = REPORTS_DIR / "13_1_record_and_fold_reconciliation.csv"
MODEL_CFG_CSV = REPORTS_DIR / "13_1_model_configuration_audit.csv"
FEAT_DIM_CSV = REPORTS_DIR / "13_1_feature_dimension_audit.csv"
HARDCODING_CSV = REPORTS_DIR / "13_1_hardcoding_audit.csv"
METRICS_CSV = REPORTS_DIR / "13_1_independently_recomputed_metrics.csv"
PRED_COMP_CSV = REPORTS_DIR / "13_1_prediction_reproduction_comparison.csv"
STAT_CSV = REPORTS_DIR / "13_1_statistical_reconciliation.csv"
CLAIMS_CSV = REPORTS_DIR / "13_1_scientific_claims_audit.csv"
MANIFEST_JSON = REPORTS_DIR / "13_1_reproducibility_manifest.json"
REPORT_MD = REPORTS_DIR / "13_1_discrepancy_reconciliation.md"


def test_deliverables_exist_and_nonempty():
    """Verify that all 10 Step 13.1 deliverables exist and are non-empty."""
    required = [
        RECON_CSV,
        MODEL_CFG_CSV,
        FEAT_DIM_CSV,
        HARDCODING_CSV,
        METRICS_CSV,
        PRED_COMP_CSV,
        STAT_CSV,
        CLAIMS_CSV,
        MANIFEST_JSON,
        REPORT_MD,
    ]
    for p in required:
        assert p.exists(), f"Missing required deliverable: {p.name}"
        assert p.stat().st_size > 0, f"Deliverable is empty: {p.name}"


def test_rs1_and_fold_reconciliation():
    """Verify that RS1 has exactly 212 records in Fold 0 and Fold 3 Val has 89 sequences."""
    df_recon = pd.read_csv(RECON_CSV)
    
    # Check RS1 audit record
    rs1_row = df_recon[df_recon["item_audited"].str.contains("RS1")].iloc[0]
    assert "212 records" in rs1_row["observed_value"]
    assert "Fold 0" in rs1_row["observed_value"]
    assert rs1_row["disposition"] == "DISCREPANCY_RESOLVED_TEXTUAL_ERROR"

    # Check Fold 0 and Fold 3 allocation
    f0_row = df_recon[df_recon["item_audited"].str.contains("Fold 0")].iloc[0]
    assert "Val=290 seqs" in f0_row["observed_value"]
    assert "Overlap=0" in f0_row["observed_value"]

    f3_row = df_recon[df_recon["item_audited"].str.contains("Fold 3")].iloc[0]
    assert "Val=89 seqs" in f3_row["observed_value"]
    assert "Overlap=0" in f3_row["observed_value"]


def test_model_configuration_audit():
    """Verify that all linear classifiers executed with C=1.0 and typo is resolved."""
    df_model = pd.read_csv(MODEL_CFG_CSV)
    
    # Linear models
    linear_models = ["logreg_physicochemical", "ridge_esm2", "combined_phys_esm2"]
    for mid in linear_models:
        row = df_model[df_model["model_id"] == mid].iloc[0]
        assert row["step10_config_C"] == 1.0
        assert row["step13_executed_C"] == 1.0

    # ESM-2 model typo resolution
    ridge_row = df_model[df_model["model_id"] == "ridge_esm2"].iloc[0]
    assert ridge_row["discrepancy_status"] == "TEXTUAL_TYPO_IN_REPORT"


def test_feature_dimension_audit():
    """Verify that predictive models use 8 bulk features (488 combined) and not 28 / 508."""
    df_feat = pd.read_csv(FEAT_DIM_CSV)
    
    pred_row = df_feat[df_feat["feature_context"].str.contains("Supervised Baseline")].iloc[0]
    assert pred_row["feature_count"] == 8
    assert pred_row["esm2_dimensions"] == "480" or pred_row["esm2_dimensions"] == 480
    assert pred_row["combined_dimension"] == 488

    diag_row = df_feat[df_feat["feature_context"].str.contains("Univariate")].iloc[0]
    assert diag_row["feature_count"] == 28


def test_exact_prediction_reproduction():
    """Verify exact reproduction across all 5 models against Step 10 and Step 13."""
    df_pred_comp = pd.read_csv(PRED_COMP_CSV)
    assert len(df_pred_comp) == 5

    for _, r in df_pred_comp.iterrows():
        assert r["evaluated_records"] == 687
        assert r["class_prediction_mismatches_vs_step10"] == 0
        assert r["class_prediction_mismatches_vs_step13"] == 0
        assert r["max_abs_prob_diff_vs_step10"] < 1e-12
        assert r["max_abs_prob_diff_vs_step13"] < 1e-12
        assert r["reproduction_verdict"] == "EXACT_REPRODUCIBILITY"


def test_recalculated_metrics_and_subgroups():
    """Verify recomputed family-level metrics and evidence-tier subgroup breakdowns."""
    df_metrics = pd.read_csv(METRICS_CSV)
    
    ridge_fam = df_metrics[(df_metrics["model_id"] == "ridge_esm2") & (df_metrics["evaluation_unit"].str.startswith("Homology"))].iloc[0]
    assert np.isclose(ridge_fam["accuracy"], 51 / 63, atol=1e-4)
    assert np.isclose(ridge_fam["balanced_accuracy"], 0.7946, atol=1e-4)
    assert np.isclose(ridge_fam["macro_f1"], 0.7946, atol=1e-4)
    assert np.isclose(ridge_fam["subgroup_protein_macro_f1"], 0.8154, atol=1e-4)
    assert np.isclose(ridge_fam["subgroup_transcript_macro_f1"], 0.6857, atol=1e-4)


def test_statistical_inference_and_holm_adjustment():
    """Verify exact McNemar tests, Holm adjustments, and inconclusive status for phys models."""
    df_stat = pd.read_csv(STAT_CSV)
    assert len(df_stat) == 4

    vs_dummy = df_stat[df_stat["comparison"] == "ridge_esm2 vs dummy_majority"].iloc[0]
    assert vs_dummy["significant_fwer_05"] == True
    assert vs_dummy["holm_adjusted_pvalue"] < 0.05

    vs_rf = df_stat[df_stat["comparison"] == "ridge_esm2 vs rf_physicochemical"].iloc[0]
    assert np.isclose(vs_rf["mcnemar_exact_raw_pvalue"], 0.0201, atol=1e-3)
    assert np.isclose(vs_rf["holm_adjusted_pvalue"], 0.0602, atol=1e-3)
    assert vs_rf["significant_fwer_05"] == False
    assert vs_rf["delta_macro_f1_ci_lower_95"] > 0

    vs_lr = df_stat[df_stat["comparison"] == "ridge_esm2 vs logreg_physicochemical"].iloc[0]
    assert np.isclose(vs_lr["mcnemar_exact_raw_pvalue"], 0.0294, atol=1e-3)
    assert np.isclose(vs_lr["holm_adjusted_pvalue"], 0.0602, atol=1e-3)
    assert vs_lr["significant_fwer_05"] == False
    assert vs_lr["delta_macro_f1_ci_lower_95"] > 0


def test_report_consistency_with_underlying_tables():
    """Verify that report markdown includes all key data values from tables."""
    text = REPORT_MD.read_text(encoding="utf-8")
    
    # RS1 count
    assert "212" in text
    # Feature counts
    assert "488" in text
    # Holm p-value
    assert "0.0602" in text
    # Exact reproduction verdict
    assert "EXACT_REPRODUCIBILITY" in text
    # Family balanced accuracy
    assert "0.7946" in text
