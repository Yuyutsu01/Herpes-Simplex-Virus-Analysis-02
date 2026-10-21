"""
Unit and integration tests for Step 12: Interpretable ESM-2 Representation Analysis.

Verifies:
1. Deliverable file presence and JSON manifest structure.
2. Embedding coefficient stability across cross-validation folds.
3. Multicollinearity, condition number, and low-rank intrinsic geometry.
4. Sequence-level composition, FDR correction, and length confounder adjustment.
5. Stratification across biological evidence tiers (52 protein vs 11 transcript).
6. Compliance with epistemological boundaries and prohibition of token-level fabrication.
"""

from pathlib import Path
import json
import pandas as pd
import numpy as np
import pytest

BASE_DIR = Path(__file__).resolve().parent.parent
REPORTS_DIR = BASE_DIR / "reports"

COEF_CSV = REPORTS_DIR / "12_embedding_coefficient_stability.csv"
FEAT_CSV = REPORTS_DIR / "12_feature_association.csv"
DIAG_CSV = REPORTS_DIR / "12_sequence_level_diagnostics.csv"
REPORT_MD = REPORTS_DIR / "12_representation_interpretation.md"
LIMITATIONS_MD = REPORTS_DIR / "12_validation_and_limitations.md"
MANIFEST_JSON = REPORTS_DIR / "12_reproducibility_manifest.json"


def test_deliverable_files_exist():
    """Verify that all six required Step 12 deliverables exist and are non-empty."""
    for p in [COEF_CSV, FEAT_CSV, DIAG_CSV, REPORT_MD, LIMITATIONS_MD, MANIFEST_JSON]:
        assert p.exists(), f"Missing required file: {p.name}"
        assert p.stat().st_size > 0, f"File {p.name} is empty"


def test_manifest_structure_and_hashes():
    """Verify reproducibility manifest format, hashes, and parameter logging."""
    with open(MANIFEST_JSON, "r", encoding="utf-8") as f:
        data = json.load(f)

    assert data["manifest_version"] == "1.0.0_step12_interpretation"
    assert data["random_seed"] == 42
    assert "input_hashes_sha256" in data

    # Verify all 6 critical inputs have logged hashes
    expected_inputs = [
        "08_1_strict_temporal_benchmark.csv",
        "valid_proteins_metadata.csv",
        "esm2_t12_representatives_metadata.csv",
        "esm2_t12_representatives_embeddings.npy",
        "09_group_split_manifest.csv",
        "10_out_of_fold_predictions.csv"
    ]
    for inp in expected_inputs:
        assert inp in data["input_hashes_sha256"]
        assert len(data["input_hashes_sha256"][inp]) == 64


def test_embedding_stability_table():
    """Verify ESM-2 linear model coefficient stability across 5 folds."""
    df = pd.read_csv(COEF_CSV)
    assert len(df) == 480, "Must contain exactly 480 embedding dimensions"

    expected_cols = [
        "embedding_dimension", "mean_coefficient", "std_coefficient",
        "signal_to_noise_ratio", "direction", "fold_0_coef", "fold_1_coef",
        "fold_2_coef", "fold_3_coef", "fold_4_coef", "mean_magnitude_rank"
    ]
    for col in expected_cols:
        assert col in df.columns

    # Verify fold coefficients are correlated across folds
    fold_matrix = df[[f"fold_{i}_coef" for i in range(5)]].values
    corrs = []
    for i in range(5):
        for j in range(i + 1, 5):
            r = np.corrcoef(fold_matrix[:, i], fold_matrix[:, j])[0, 1]
            corrs.append(r)

    mean_corr = np.mean(corrs)
    assert mean_corr > 0.65, f"Expected mean fold correlation > 0.65, got {mean_corr:.4f}"

    # Verify SNR is non-negative and finite
    assert (df["signal_to_noise_ratio"] >= 0).all()
    assert not df["signal_to_noise_ratio"].isnull().any()


def test_multicollinearity_and_geometry():
    """Verify severe collinearity and low-rank intrinsic dimensionality of ESM-2 embeddings."""
    with open(MANIFEST_JSON, "r", encoding="utf-8") as f:
        data = json.load(f)

    geom = data["multicollinearity_and_geometry"]
    # Condition number must be severe (> 1e12)
    assert geom["covariance_condition_number"] > 1e12
    # 95% of variance must be captured by 25 or fewer PCs out of 480
    assert geom["pcs_for_95_pct_variance"] <= 25
    assert geom["pcs_for_80_pct_variance"] <= 10
    # Mean off-diagonal correlation must be substantial
    assert geom["mean_absolute_pairwise_correlation"] > 0.25


def test_feature_association_fdr_and_confounders():
    """Verify sequence composition tests, BH FDR correction, and length confounder adjustment."""
    df = pd.read_csv(FEAT_CSV)
    # 8 physicochemical + 20 canonical amino acids = 28 features
    assert len(df) == 28

    # Crucial scientific check: ZERO features must survive FDR at alpha=0.05
    assert (df["significant_fdr_05"] == False).all()
    assert (df["significant_fwer_05"] == False).all()
    assert (df["bh_fdr_adjusted_pvalue"] > 0.05).all()

    # Check Lysine (aa_K): unadjusted p may be < 0.05, but length-adjusted Wald p must be > 0.05
    k_row = df[df["feature"] == "aa_K"].iloc[0]
    assert k_row["raw_mann_whitney_pvalue"] < 0.05
    assert k_row["length_adjusted_wald_pvalue"] > 0.05, "aa_K must not be significant after controlling for length"


def test_stratification_by_evidence_tier():
    """Verify performance differences between protein-validated and transcript-only families."""
    df = pd.read_csv(DIAG_CSV)
    assert len(df) == 63, "Diagnostics table must contain all 63 direct benchmark families"

    # 52 protein-validated vs 11 transcript-only
    protein_sub = df[df["measured_molecule"] == "TRANSCRIPT_AND_PROTEIN"]
    transcript_sub = df[df["measured_molecule"] == "TRANSCRIPT"]

    assert len(protein_sub) == 52
    assert len(transcript_sub) == 11

    # Protein-validated families must have higher accuracy (82.7% vs 72.7%)
    acc_protein = protein_sub["is_correct"].mean()
    acc_transcript = transcript_sub["is_correct"].mean()

    assert np.isclose(acc_protein, 43 / 52, atol=1e-4)
    assert np.isclose(acc_transcript, 8 / 11, atol=1e-4)
    assert acc_protein > acc_transcript

    # Overall correctness must be 51/63 (80.95%)
    assert df["is_correct"].sum() == 51


def test_epistemological_boundaries_and_limitations_report():
    """Verify that reports enforce strict boundaries against token fabrication and causal claims."""
    report_text = REPORT_MD.read_text(encoding="utf-8")
    limitations_text = LIMITATIONS_MD.read_text(encoding="utf-8")

    # Verify token-level reality check
    assert "mean-pooled" in report_text.lower()
    assert "token-level hidden states" in report_text or "token-level attributions are not fabricated" in report_text.lower()

    # Verify causal inference warning
    assert "Predictive Association vs. Biological Mechanism" in report_text
    assert "TAATGARAT" in limitations_text or "VP16" in limitations_text
    assert "UNPROVEN" in limitations_text
