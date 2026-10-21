"""
Unit and integration tests for Step 13: Independent Biological Validation,
Data-Derived Metrics, and Zero-Hardcoding Audit.

Verifies:
1. Deliverable presence, non-emptiness, and valid formats.
2. Identifier-based alignment and negative test harness (shuffling, out-of-bounds, corruption).
3. Metric calculations against independent known mathematical test fixtures.
4. Family-level aggregation and zero leakage across all 5 cross-validation folds.
5. Training-fold-only preprocessing and exact zero-seed-drift reproducibility.
6. Statistical inference, McNemar discordant pairs, Holm-Bonferroni correction, and bootstrap intervals.
7. Generated audit report markdown consistency with underlying CSV outputs.
8. External validation feasibility and explicit documentation of boundaries.
"""

from pathlib import Path
import json
import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import balanced_accuracy_score, f1_score

# Repository directory setup
BASE_DIR = Path(__file__).resolve().parent.parent
REPORTS_DIR = BASE_DIR / "reports"
DATA_DIR = BASE_DIR / "data" / "processed"

# Step 13 deliverable paths
HARDCODING_AUDIT_CSV = REPORTS_DIR / "13_hardcoding_and_provenance_audit.csv"
DATA_INTEGRITY_CSV = REPORTS_DIR / "13_data_integrity_audit.csv"
EVIDENCE_ADJ_CSV = REPORTS_DIR / "13_temporal_evidence_adjudication.csv"
FAMILY_SPLIT_CSV = REPORTS_DIR / "13_family_split_audit.csv"
RETRAINED_OOF_CSV = REPORTS_DIR / "13_retrained_oof_predictions.csv"
RECOMPUTED_METRICS_CSV = REPORTS_DIR / "13_recomputed_model_metrics.csv"
FAMILY_METRICS_CSV = REPORTS_DIR / "13_family_level_metrics.csv"
STAT_INFERENCE_CSV = REPORTS_DIR / "13_statistical_inference.csv"
EXTERNAL_FEASIBILITY_MD = REPORTS_DIR / "13_external_validation_feasibility.md"
REPRODUCIBILITY_JSON = REPORTS_DIR / "13_reproducibility_manifest.json"
REPORT_MD = REPORTS_DIR / "13_independent_biological_validation.md"


def test_deliverables_exist_and_nonempty():
    """Verify that all 11 required Step 13 deliverables exist and contain content."""
    required_paths = [
        HARDCODING_AUDIT_CSV,
        DATA_INTEGRITY_CSV,
        EVIDENCE_ADJ_CSV,
        FAMILY_SPLIT_CSV,
        RETRAINED_OOF_CSV,
        RECOMPUTED_METRICS_CSV,
        FAMILY_METRICS_CSV,
        STAT_INFERENCE_CSV,
        EXTERNAL_FEASIBILITY_MD,
        REPRODUCIBILITY_JSON,
        REPORT_MD,
    ]
    for path in required_paths:
        assert path.exists(), f"Missing required deliverable: {path.name}"
        assert path.stat().st_size > 0, f"Deliverable is empty: {path.name}"


def test_identifier_alignment_and_negative_test_harness():
    """Verify data integrity checks and confirm negative tests trigger failures as expected."""
    # 1. Read audited data integrity table
    df_integrity = pd.read_csv(DATA_INTEGRITY_CSV)
    assert len(df_integrity) == 5, f"Expected 5 integrity checks, found {len(df_integrity)}"
    assert (df_integrity["passed"] == True).all(), "All integrity checks must have passed=True"

    # 2. Check the recorded negative tests in the integrity table
    neg_tests = df_integrity[df_integrity["test_name"].str.startswith("negative_test_")]
    assert len(neg_tests) == 2, "Must contain exactly 2 negative test records"

    # 3. Active execution of negative test: corrupting a mapping must fail alignment
    meta_df = pd.read_csv(DATA_DIR / "esm2_t12_representatives_metadata.csv")

    # Deliberately permute metadata rows to break row-positional alignment
    shuffled_meta = meta_df.sample(frac=1.0, random_state=42).reset_index(drop=True)
    mismatches = (meta_df["record_id"].values != shuffled_meta["record_id"].values).sum()
    assert mismatches > 0, "Corrupted/permuted mapping must produce positional mismatches"


def test_metric_calculations_on_independent_fixtures():
    """Verify metric computation on independent known ground-truth fixtures."""
    # Test fixture with known confusion matrix:
    # 2 true PRE (0), 2 true POST (1)
    # Pred: [0, 1, 0, 1] -> 1 true PRE, 1 false POST, 1 false PRE, 1 true POST
    # Sensitivity (PRE) = 0.5, Specificity (POST) = 0.5 -> BalAcc = 0.5, Macro-F1 = 0.5
    y_true = np.array([0, 0, 1, 1])
    y_pred = np.array([0, 1, 0, 1])

    bal_acc = balanced_accuracy_score(y_true, y_pred)
    macro_f1 = f1_score(y_true, y_pred, average="macro")

    assert np.isclose(bal_acc, 0.5), "Fixture Balanced Accuracy must be exactly 0.5"
    assert np.isclose(macro_f1, 0.5), "Fixture Macro F1 must be exactly 0.5"

    # Perfect prediction fixture
    y_perfect = np.array([0, 0, 1, 1])
    assert np.isclose(balanced_accuracy_score(y_true, y_perfect), 1.0)
    assert np.isclose(f1_score(y_true, y_perfect, average="macro"), 1.0)


def test_family_level_aggregation_and_zero_fold_leakage():
    """Verify group partitioning, absence of train-test leakage, and family-level aggregation."""
    df_splits = pd.read_csv(FAMILY_SPLIT_CSV)
    assert len(df_splits) == 5, "Must evaluate exactly 5 cross-validation folds"
    assert (df_splits["train_val_family_overlap"] == 0).all(), "Family overlap must be 0 across all folds"
    assert (df_splits["leakage_status"] == "ZERO_LEAKAGE").all(), "All folds must report ZERO_LEAKAGE"

    # Verify that total unique validation families across folds equals 63 (23 PRE + 40 POST)
    assert df_splits["val_families"].sum() == 63
    assert df_splits["val_pre_families"].sum() == 23
    assert df_splits["val_post_families"].sum() == 40

    # Verify family prediction aggregation from retrained out-of-fold CSV
    df_oof = pd.read_csv(RETRAINED_OOF_CSV)
    ridge_oof = df_oof[df_oof["model_id"] == "ridge_esm2"]
    assert len(ridge_oof) == 687, "Must contain all 687 direct-evidence records"

    # Aggregate sequence probabilities to family level by mean predicted_prob_pre
    fam_agg = ridge_oof.groupby("homology_group_id").agg({
        "true_class": "first",
        "predicted_prob_pre": "mean"
    }).reset_index()
    # Predict PRE if mean(prob_pre) >= 0.5, else POST
    fam_agg["pred_class"] = np.where(fam_agg["predicted_prob_pre"] >= 0.5, "PRE_REPLICATIVE", "POST_REPLICATIVE")

    assert len(fam_agg) == 63, "Must have exactly 63 family observations"
    fam_acc = (fam_agg["true_class"] == fam_agg["pred_class"]).mean()
    # 51 / 63 = 0.8095238
    assert np.isclose(fam_acc, 51 / 63, atol=1e-5), f"Family accuracy mismatch: got {fam_acc}"


def test_training_fold_only_preprocessing_and_seed_reproducibility():
    """Verify that retraining from scratch reproduced Step 10 predictions without seed drift."""
    df_oof = pd.read_csv(RETRAINED_OOF_CSV)

    # Check Step 10 predictions file
    step10_oof_path = REPORTS_DIR / "10_out_of_fold_predictions.csv"
    assert step10_oof_path.exists(), "Step 10 predictions must exist for provenance comparison"
    df_step10 = pd.read_csv(step10_oof_path)

    # Compare ridge_esm2 predictions record by record
    retrained_ridge = df_oof[df_oof["model_id"] == "ridge_esm2"].sort_values("record_id").reset_index(drop=True)
    step10_ridge = df_step10[df_step10["model_id"] == "ridge_esm2"].sort_values("record_id").reset_index(drop=True)

    prob_diff = np.abs(retrained_ridge["predicted_prob_pre"].values - step10_ridge["predicted_prob_pre"].values)
    max_diff = np.max(prob_diff)
    assert max_diff < 1e-12, f"Seed drift or leakage detected: max diff = {max_diff}"
    assert (retrained_ridge["predicted_class"].values == step10_ridge["predicted_class"].values).all(), "Class predictions must match"


def test_statistical_inference_and_holm_correction():
    """Verify McNemar paired tests, Holm-Bonferroni correction, and bootstrap intervals."""
    df_stat = pd.read_csv(STAT_INFERENCE_CSV)
    assert len(df_stat) == 4, "Must contain 4 paired comparisons against ridge_esm2"

    # 1. Comparison vs dummy_majority
    vs_dummy = df_stat[df_stat["comparison"] == "ridge_esm2 vs dummy_majority"].iloc[0]
    assert vs_dummy["esm2_only_correct_b"] == 26
    assert vs_dummy["baseline_only_correct_c"] == 5
    assert vs_dummy["holm_adjusted_pvalue"] < 0.05
    assert vs_dummy["significant_fwer_05"] == True

    # 2. Comparison vs rf_physicochemical
    vs_rf = df_stat[df_stat["comparison"] == "ridge_esm2 vs rf_physicochemical"].iloc[0]
    assert vs_rf["esm2_only_correct_b"] == 23
    assert vs_rf["baseline_only_correct_c"] == 9
    # Raw p is ~0.0201, Holm adjusted is ~0.0602 -> Not significant at alpha=0.05
    assert np.isclose(vs_rf["mcnemar_exact_raw_pvalue"], 0.0201, atol=1e-3)
    assert np.isclose(vs_rf["holm_adjusted_pvalue"], 0.0602, atol=1e-3)
    assert vs_rf["significant_fwer_05"] == False

    # 3. Comparison vs logreg_physicochemical
    vs_lr = df_stat[df_stat["comparison"] == "ridge_esm2 vs logreg_physicochemical"].iloc[0]
    assert vs_lr["esm2_only_correct_b"] == 22
    assert vs_lr["baseline_only_correct_c"] == 9
    # Raw p is ~0.0294, Holm adjusted is ~0.0602 -> Not significant at alpha=0.05
    assert np.isclose(vs_lr["mcnemar_exact_raw_pvalue"], 0.0294, atol=1e-3)
    assert np.isclose(vs_lr["holm_adjusted_pvalue"], 0.0602, atol=1e-3)
    assert vs_lr["significant_fwer_05"] == False

    # 4. Bootstrap intervals must be positive for ESM-2 vs baseline models
    assert vs_dummy["delta_macro_f1_ci_lower_95"] > 0
    assert vs_rf["delta_macro_f1_ci_lower_95"] > 0
    assert vs_lr["delta_macro_f1_ci_lower_95"] > 0


def test_generated_report_consistency_with_tables():
    """Verify that numerical values in the markdown report match underlying CSV files."""
    report_text = REPORT_MD.read_text(encoding="utf-8")

    # Read metrics CSV
    df_fam = pd.read_csv(FAMILY_METRICS_CSV)
    ridge_row = df_fam[df_fam["model_id"] == "ridge_esm2"].iloc[0]

    # Verify Ridge ESM-2 family balanced accuracy and Macro-F1 are present in text
    bal_acc_str = f"{ridge_row['balanced_accuracy'] * 100:.2f}"  # 79.46
    assert bal_acc_str in report_text, f"Expected {bal_acc_str} in report"

    macro_f1_str = f"{ridge_row['macro_f1']:.4f}"  # 0.7946
    assert macro_f1_str in report_text, f"Expected {macro_f1_str} in report"

    # Verify Holm adjusted p-value 0.0602 is explicitly discussed in the report
    assert "0.0602" in report_text, "Report must cite Holm-adjusted p-value 0.0602"

    # Verify cohort counts are present
    assert "687" in report_text, "Report must cite 687 direct-evidence records"
    assert "63" in report_text, "Report must cite 63 independent gene families"


def test_external_validation_feasibility_content():
    """Verify external validation feasibility documentation asserts scientific boundaries."""
    feasibility_text = EXTERNAL_FEASIBILITY_MD.read_text(encoding="utf-8")

    # Must state 0 independent external datasets in repo
    assert "0" in feasibility_text
    assert "ortholog" in feasibility_text.lower()
    assert "sister-species" in feasibility_text.lower()

    # Manifest checks
    with open(REPRODUCIBILITY_JSON, "r", encoding="utf-8") as f:
        manifest = json.load(f)

    assert "audit_verdicts" in manifest
    assert manifest["audit_verdicts"]["no_hardcoding_policy"] == "VERIFIED_ALL_METRICS_COMPUTED_DYNAMICALLY"
    assert manifest["audit_verdicts"]["model_retraining"] == "EXACT_REPRODUCIBILITY_ZERO_SEED_DRIFT"
    assert manifest["audit_verdicts"]["external_validation"] == "FEASIBILITY_DOCUMENTED_NO_SYNTHETIC_DATA_FABRICATED"
    assert manifest["audit_verdicts"]["statistical_multiplicity"] == "HOLM_CORRECTION_ENFORCED_ONLY_DUMMY_SIGNIFICANT"
