"""
Unit and integration tests for Step 12.1: Validation of Representation Interpretation and Feature Analysis.

Verifies:
1. Deliverable presence and format closure.
2. Data alignment checks across manifest, benchmark, and embedding rows (all PASSED).
3. Correction of PCA generalization claim (in-sample training vs. held-out 5-fold CV).
4. Reassessment of embedding geometry (Roy-Vetterli effective rank and condition number).
5. Sequence composition inference audit (ANCOVA, Cliff's Delta, inconclusive status).
6. Coefficient stability audit (collinear dimension rank swapping).
7. Reproducibility manifest audit verdicts and hashes.
"""

from pathlib import Path
import json
import pandas as pd
import numpy as np
import pytest

BASE_DIR = Path(__file__).resolve().parent.parent
REPORTS_DIR = BASE_DIR / "reports"

ALIGN_CSV = REPORTS_DIR / "12_1_data_alignment_checks.csv"
PCA_CSV = REPORTS_DIR / "12_1_pca_heldout_metrics.csv"
GEOM_CSV = REPORTS_DIR / "12_1_embedding_geometry.csv"
COMP_CSV = REPORTS_DIR / "12_1_composition_inference_audit.csv"
COEF_CSV = REPORTS_DIR / "12_1_coefficient_stability_audit.csv"
REPORT_MD = REPORTS_DIR / "12_1_interpretation_validation_audit.md"
MANIFEST_JSON = REPORTS_DIR / "12_1_reproducibility_manifest.json"


def test_deliverable_files_exist():
    """Verify that all seven required Step 12.1 deliverables exist and are non-empty."""
    for p in [ALIGN_CSV, PCA_CSV, GEOM_CSV, COMP_CSV, COEF_CSV, REPORT_MD, MANIFEST_JSON]:
        assert p.exists(), f"Missing required file: {p.name}"
        assert p.stat().st_size > 0, f"File {p.name} is empty"


def test_data_alignment_checks():
    """Verify that all alignment checks passed with zero row-count or ID drift."""
    df = pd.read_csv(ALIGN_CSV)
    assert len(df) == 6, "Must contain exactly 6 explicit alignment checks"
    assert (df["status"] == "PASSED").all(), "All alignment checks must have status PASSED"

    # Specific check values
    rep_check = df[df["check_name"] == "total_representatives_count"].iloc[0]
    assert "3111" in str(rep_check["observed_value"])

    dir_check = df[df["check_name"] == "direct_cohort_record_count"].iloc[0]
    assert "687" in str(dir_check["observed_value"])

    hom_check = df[df["check_name"] == "homology_groups_count_and_partition"].iloc[0]
    assert "total=63" in str(hom_check["observed_value"])
    assert "pre=23" in str(hom_check["observed_value"])
    assert "post=40" in str(hom_check["observed_value"])


def test_pca_generalization_correction():
    """Verify the correction of PCA generalization: in-sample vs held-out 5-fold CV."""
    df = pd.read_csv(PCA_CSV)
    
    # In-sample model must show high training accuracy (> 85%)
    in_sample = df[df["model_id"] == "pca_18_in_sample"].iloc[0]
    assert in_sample["seq_accuracy"] > 0.85
    assert "IN_SAMPLE_TRAINING" in str(in_sample["evaluation_protocol"])

    # Held-out PCA k=18 model must show collapsed generalization (~50% balanced acc)
    held_out_18 = df[df["model_id"] == "pca_18_held_out"].iloc[0]
    assert held_out_18["family_balanced_acc"] < 0.55
    assert held_out_18["family_macro_f1"] < 0.50
    assert "HELD_OUT" in str(held_out_18["evaluation_protocol"])

    # Full ESM-2 model must outperform PCA held-out models by a large margin
    full_esm = df[df["model_id"] == "ridge_esm2_full_480"].iloc[0]
    assert full_esm["family_balanced_acc"] > 0.75
    assert full_esm["family_macro_f1"] > 0.75


def test_embedding_geometry_and_effective_rank():
    """Verify Roy-Vetterli effective rank and condition number reassessment."""
    df = pd.read_csv(GEOM_CSV)
    assert len(df) == 2, "Must evaluate both raw covariance and correlation matrix"

    corr_row = df[df["matrix_representation"] == "Standardized Correlation Matrix"].iloc[0]
    assert 12.0 < corr_row["effective_rank_roy_vetterli"] < 16.0
    assert corr_row["pcs_for_95pct_variance"] <= 30
    assert corr_row["near_zero_noise_eigenvalues_count"] >= 1

    raw_row = df[df["matrix_representation"] == "Raw Unstandardized Covariance"].iloc[0]
    assert 4.0 < raw_row["effective_rank_roy_vetterli"] < 8.0
    assert raw_row["pcs_for_95pct_variance"] <= 20


def test_sequence_composition_audit():
    """Verify sequence composition inference, FDR adjustments, and inconclusive status."""
    df = pd.read_csv(COMP_CSV)
    assert len(df) == 28, "Must test exactly 28 features"

    # All features must fail FDR correction at 0.05
    assert (df["significant_fdr_05"] == False).all()
    assert (df["bh_fdr_adjusted_pvalue"] > 0.40).all()

    # All features must be classified as inconclusive
    assert df["statistical_inference_conclusion"].str.contains("INCONCLUSIVE_DUE_TO_SAMPLE_SIZE").all()

    # Lysine check: Cliff's delta positive (~0.337), length-adjusted ANCOVA p > 0.05
    k_row = df[df["feature"] == "aa_K"].iloc[0]
    assert np.isclose(k_row["cliffs_delta"], 0.336956, atol=1e-3)
    assert k_row["length_adjusted_ancova_pvalue"] > 0.05


def test_coefficient_stability_and_rank_swapping():
    """Verify detection of collinear dimension rank swapping across folds."""
    df = pd.read_csv(COEF_CSV)
    assert len(df) == 480, "Must audit all 480 dimensions"

    # Verify that rank swapping was detected for multiple dimensions
    swapping_dims = df[df["collinear_rank_swapping_detected"] == True]
    assert len(swapping_dims) >= 5, "Must detect at least 5 rank-swapping dimensions"

    # Verify that Dim 422 has a large rank range
    dim_422 = df[df["embedding_dimension"] == 422].iloc[0]
    assert dim_422["rank_fluctuation_range"] >= 50


def test_manifest_verdicts():
    """Verify reproducibility manifest structure and verdicts."""
    with open(MANIFEST_JSON, "r", encoding="utf-8") as f:
        data = json.load(f)

    assert data["manifest_version"] == "1.0.0_step12_1_validation_audit"
    verdicts = data["audit_verdicts"]
    assert verdicts["data_alignment"] == "PASSED_ZERO_DRIFT"
    assert verdicts["pca_held_out_generalization"] == "FAILED_IN_SAMPLE_OVERFIT_ONLY"
    assert "INCONCLUSIVE" in verdicts["composition_inference"]
