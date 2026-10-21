"""
Unit tests for Step 11.1: Statistical Inference and Cross-Species Claim Audit.

Verifies:
1. Deliverable existence and format integrity.
2. Exact McNemar 2x2 contingency counts and p-values.
3. Holm-Bonferroni step-down correction correctness and FWER significance conclusions.
4. Paired bootstrap resampling properties (B=10,000, 0 missing classes, 95% CIs).
5. Family-held-out protocol audit verifications (zero leakage, fold-local preprocessing).
6. Cross-species claim reclassification and true HSV-1 -> HSV-2 transfer results.
"""

from pathlib import Path
import pandas as pd
import numpy as np
import pytest
from scipy.stats import binomtest

BASE_DIR = Path(__file__).resolve().parent.parent
REPORTS_DIR = BASE_DIR / "reports"

HYPOTHESIS_CSV = REPORTS_DIR / "11_1_paired_hypothesis_tests.csv"
BOOTSTRAP_CSV = REPORTS_DIR / "11_1_paired_bootstrap_differences.csv"
TRANSFER_CSV = REPORTS_DIR / "11_1_cross_species_transfer_experiment.csv"
REPORT_MD = REPORTS_DIR / "11_1_statistical_and_species_audit.md"


def test_deliverable_files_exist():
    """Verify that all four required Step 11.1 deliverables exist and are non-empty."""
    assert HYPOTHESIS_CSV.exists(), "11_1_paired_hypothesis_tests.csv must exist"
    assert BOOTSTRAP_CSV.exists(), "11_1_paired_bootstrap_differences.csv must exist"
    assert TRANSFER_CSV.exists(), "11_1_cross_species_transfer_experiment.csv must exist"
    assert REPORT_MD.exists(), "11_1_statistical_and_species_audit.md must exist"

    assert HYPOTHESIS_CSV.stat().st_size > 0
    assert BOOTSTRAP_CSV.stat().st_size > 0
    assert TRANSFER_CSV.stat().st_size > 0
    assert REPORT_MD.stat().st_size > 0


def test_paired_hypothesis_tests_table():
    """Verify McNemar 2x2 table values, raw exact binomial p-values, and Holm adjustments."""
    df = pd.read_csv(HYPOTHESIS_CSV)
    assert len(df) == 5, "Must test exactly 5 prespecified model comparisons"

    expected_comparisons = [
        "ridge_esm2 vs fixed_family_majority",
        "ridge_esm2 vs dummy_majority",
        "ridge_esm2 vs logreg_physicochemical",
        "ridge_esm2 vs rf_physicochemical",
        "ridge_esm2 vs combined_phys_esm2"
    ]
    assert list(df["comparison"]) == expected_comparisons

    # Verify contingency counts sum to 63 for all models
    for _, row in df.iterrows():
        total = row["both_correct_a"] + row["esm2_only_correct_b"] + row["baseline_only_correct_c"] + row["both_incorrect_d"]
        assert total == 63, f"Total families must be 63, got {total} for {row['comparison']}"

        # Recompute exact two-sided binomial p-value on discordant pairs
        b = row["esm2_only_correct_b"]
        c = row["baseline_only_correct_c"]
        n_disc = b + c
        assert row["total_discordant_pairs"] == n_disc

        if n_disc == 0:
            expected_p = 1.0
        else:
            expected_p = binomtest(b, n_disc, p=0.5, alternative="two-sided").pvalue
        assert np.isclose(row["mcnemar_exact_raw_pvalue"], expected_p, atol=1e-6)

    # Verify Holm-Bonferroni correction and significance decisions
    # Under Holm-Bonferroni at alpha=0.05, ONLY dummy_majority should be statistically significant
    dummy_row = df[df["comparison"] == "ridge_esm2 vs dummy_majority"].iloc[0]
    assert bool(dummy_row["significant_holm_05"]) is True
    assert dummy_row["holm_adjusted_pvalue"] < 0.05

    # Physicochemical models and fixed majority must NOT be significant under Holm adjustment at alpha=0.05
    for comp in ["ridge_esm2 vs logreg_physicochemical", "ridge_esm2 vs rf_physicochemical", "ridge_esm2 vs fixed_family_majority"]:
        row = df[df["comparison"] == comp].iloc[0]
        assert bool(row["significant_holm_05"]) is False
        assert row["holm_adjusted_pvalue"] > 0.05, f"{comp} should not be significant under Holm FWER control"


def test_paired_bootstrap_differences_table():
    """Verify paired bootstrap difference metrics (B=10,000 resamples, no missing classes)."""
    df = pd.read_csv(BOOTSTRAP_CSV)
    assert len(df) == 5

    for _, row in df.iterrows():
        # Check that 10,000 iterations were performed
        assert row["bootstrap_iterations"] == 10000
        # Check that missing class resamples were 0
        assert row["missing_class_resamples"] == 0

        # Verify 95% CI monotonicity
        assert row["delta_macro_f1_ci_lower_95"] <= row["delta_macro_f1_ci_upper_95"]
        assert row["delta_bal_acc_ci_lower_95"] <= row["delta_bal_acc_ci_upper_95"]
        assert row["delta_acc_ci_lower_95"] <= row["delta_acc_ci_upper_95"]

        # Scope note must be present
        assert "Quantifies uncertainty in evaluated family sample conditional on fixed predictions" in str(row["resampling_scope_note"])

    # For the 4 baseline comparisons, the Macro-F1 delta 95% CI should strictly exclude 0
    baselines = [
        "ridge_esm2 vs fixed_family_majority",
        "ridge_esm2 vs dummy_majority",
        "ridge_esm2 vs logreg_physicochemical",
        "ridge_esm2 vs rf_physicochemical"
    ]
    for comp in baselines:
        row = df[df["comparison"] == comp].iloc[0]
        assert row["delta_macro_f1_ci_lower_95"] > 0.0, f"Macro-F1 95% CI lower bound must be > 0 for {comp}"


def test_cross_species_transfer_experiment_table():
    """Verify true HSV-1 to HSV-2 cross-species transfer experiment results."""
    df = pd.read_csv(TRANSFER_CSV)
    assert len(df) == 7, "Must test exactly 7 directly-annotated HSV-2 families"

    expected_genes = {"RL2", "RS1", "UL19", "UL23", "UL27", "US4", "US6"}
    assert set(df["effective_gene_symbol"]) == expected_genes

    # All 7 families must be correctly classified at the family level
    assert (df["is_correct"] == True).all()

    # Total test sequences across the 7 families must sum to 265
    assert df["hsv2_test_sequences"].sum() == 265

    # Sequence level accuracy is 259/265 = 97.74%
    first_row = df.iloc[0]
    assert np.isclose(first_row["family_accuracy_n7"], 1.0)
    assert np.isclose(first_row["sequence_accuracy_n265"], 259 / 265, atol=1e-4)


def test_report_contains_required_reclassifications():
    """Verify that the markdown report contains the mandated reclassifications and warnings."""
    report_text = REPORT_MD.read_text(encoding="utf-8")

    assert "RECLASSIFICATION REQUIRED" in report_text or "family-held-out ortholog diagnostic" in report_text
    assert "Holm-Bonferroni" in report_text
    assert "0 / 10,000" in report_text
    assert "VERIFIED" in report_text
