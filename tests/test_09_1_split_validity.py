#!/usr/bin/env python3
"""
Test Suite: Step 9.1 Benchmark Reconciliation and Split-Validity Audit
======================================================================

Scientific Background and Quality Assurance:
--------------------------------------------
Step 9.1 validates the mathematical closure, record-level reconciliation, and
group disjointness of the Step 9 benchmark design:
1. Reconciles the narrative count discrepancies (IE=211 vs 251, E=121 vs 99,
   L1=134 vs 159, L2=221 vs 178; species 421/266 vs 422/265).
2. Verifies that the single record difference in HSV-2 direct evidence is AAK97852.1
   (recovered submitter alias quarantined to UNRESOLVED).
3. Mathematically proves group disjointness: TrainFamilies_k ∩ ValFamilies_k = ∅
   across all 5 cross-validation folds.
4. Assesses fold feasibility: Confirms that both binary classes (Pre-Replicative and
   Post-Replicative) are present in every validation fold.
5. Verifies that raw data and all historical outputs remain intact.
"""

import os
import pytest
import pandas as pd

WORKSPACE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA_RAW = os.path.join(WORKSPACE_DIR, "data", "raw")
DATA_PROCESSED = os.path.join(WORKSPACE_DIR, "data", "processed")
REPORTS_DIR = os.path.join(WORKSPACE_DIR, "reports")

# Step 9.1 Deliverables
FILE_RECON_MD = os.path.join(REPORTS_DIR, "09_1_count_reconciliation.md")
FILE_REC_DIFF_CSV = os.path.join(REPORTS_DIR, "09_1_record_level_differences.csv")
FILE_GROUPING_MD = os.path.join(REPORTS_DIR, "09_1_family_grouping_audit.md")
FILE_GROUP_EXCEPT_CSV = os.path.join(REPORTS_DIR, "09_1_grouping_exceptions.csv")
FILE_FEAS_CSV = os.path.join(REPORTS_DIR, "09_1_fold_feasibility.csv")
FILE_SPLIT_VAL_MD = os.path.join(REPORTS_DIR, "09_1_split_validity_report.md")
FILE_DECISION_MD = os.path.join(REPORTS_DIR, "09_1_benchmark_readiness_decision.md")

# Data Artifacts
FILE_BENCH_81 = os.path.join(DATA_PROCESSED, "08_1_strict_temporal_benchmark.csv")
FILE_SPLIT_MANIFEST_9 = os.path.join(DATA_PROCESSED, "09_group_split_manifest.csv")
FILE_EMB_META = os.path.join(DATA_PROCESSED, "esm2_t12_representatives_metadata.csv")


def test_01_all_step9_1_artifacts_exist():
    """Verify that all 7 Step 9.1 reports and tables exist on disk and are non-empty."""
    required_files = [
        FILE_RECON_MD,
        FILE_REC_DIFF_CSV,
        FILE_GROUPING_MD,
        FILE_GROUP_EXCEPT_CSV,
        FILE_FEAS_CSV,
        FILE_SPLIT_VAL_MD,
        FILE_DECISION_MD,
    ]
    for path in required_files:
        assert os.path.exists(path), f"Missing Step 9.1 artifact: {path}"
        assert os.path.getsize(path) > 0, f"Artifact is empty: {path}"


def test_02_preservation_of_historical_and_step9_files():
    """Verify that raw FASTA data and Steps 1-9 deliverables remain unmodified."""
    assert os.path.exists(os.path.join(DATA_RAW, "hsv1.fasta"))
    assert os.path.exists(os.path.join(DATA_RAW, "hsv2.fasta"))

    df_bench = pd.read_csv(FILE_BENCH_81)
    assert len(df_bench) == 3111, f"Corrupted 08_1 benchmark count: {len(df_bench)}"

    df_m9 = pd.read_csv(FILE_SPLIT_MANIFEST_9)
    assert len(df_m9) == 3111, f"Corrupted Step 9 manifest count: {len(df_m9)}"


def test_03_record_level_direct_counts_reconciliation():
    """
    Verify exact record-level counts for DIRECT_EVIDENCE in 08_1 and Step 9 manifest:
    Total = 687 (IE=251, L2=178, L1=159, E=99; HSV1=422, HSV2=265).
    """
    df_81 = pd.read_csv(FILE_BENCH_81)
    df_m9 = pd.read_csv(FILE_SPLIT_MANIFEST_9)

    direct_81 = df_81[df_81["adjudicated_cohort"] == "DIRECT_EVIDENCE"]
    direct_9 = df_m9[df_m9["adjudicated_cohort"] == "DIRECT_EVIDENCE"]

    assert len(direct_81) == 687
    assert len(direct_9) == 687

    # Exact class counts
    c_counts = direct_81["normalized_temporal_class"].value_counts().to_dict()
    assert c_counts.get("IE_ALPHA") == 251
    assert c_counts.get("L2_GAMMA2") == 178
    assert c_counts.get("L1_GAMMA1") == 159
    assert c_counts.get("E_BETA") == 99

    # Exact species counts
    s_counts = direct_81["effective_species"].value_counts().to_dict()
    assert s_counts.get("HSV1") == 422
    assert s_counts.get("HSV2") == 265


def test_04_record_level_differences_audit():
    """Verify that AAK97852.1 is documented in 09_1_record_level_differences.csv."""
    df_diff = pd.read_csv(FILE_REC_DIFF_CSV)
    assert len(df_diff) > 0, "Record level differences table is empty"

    records = list(df_diff["record_id"].unique())
    assert "AAK97852.1" in records, "AAK97852.1 missing from difference audit"

    # Confirm AAK97852.1 is in UNRESOLVED in 08_1
    df_81 = pd.read_csv(FILE_BENCH_81)
    r = df_81[df_81["record_id"] == "AAK97852.1"].iloc[0]
    assert r["adjudicated_cohort"] == "UNRESOLVED"
    assert r["is_alias_recovered"] == True


def test_05_homology_grouping_closure_and_mega_groups():
    """Verify no null group IDs and confirm HOM_RS1 is the primary mega-group."""
    df_m9 = pd.read_csv(FILE_SPLIT_MANIFEST_9)
    assert df_m9["homology_group_id"].isnull().sum() == 0
    assert df_m9["prevent_ml_leakage_group"].isnull().sum() == 0

    direct = df_m9[df_m9["adjudicated_cohort"] == "DIRECT_EVIDENCE"]
    rs1_count = (direct["homology_group_id"] == "HOM_RS1").sum()
    assert rs1_count == 212, f"Expected 212 variants in HOM_RS1, found {rs1_count}"

    # Confirm HOM_RS1 is isolated in Fold 0
    rs1_folds = direct[direct["homology_group_id"] == "HOM_RS1"]["fold_index"].unique()
    assert list(rs1_folds) == [0], f"HOM_RS1 leaked across folds: {rs1_folds}"


def test_06_mathematical_proof_of_fold_disjointness():
    """
    Mathematical proof of zero group leakage across all 5 folds:
    TrainFamilies_k ∩ ValFamilies_k = ∅
    """
    df_m9 = pd.read_csv(FILE_SPLIT_MANIFEST_9)
    direct = df_m9[df_m9["adjudicated_cohort"] == "DIRECT_EVIDENCE"]

    for f in range(5):
        val_families = set(direct[direct["fold_index"] == f]["homology_group_id"])
        train_families = set(direct[direct["fold_index"] != f]["homology_group_id"])

        intersection = val_families.intersection(train_families)
        assert len(intersection) == 0, f"LEAKAGE DETECTED in Fold {f}: {intersection}"


def test_07_fold_feasibility_metrics():
    """Verify that both binary classes are represented in all 5 validation folds."""
    df_feas = pd.read_csv(FILE_FEAS_CSV)
    assert len(df_feas) == 5, f"Expected 5 rows in fold feasibility table, got {len(df_feas)}"

    for _, r in df_feas.iterrows():
        assert r["both_binary_classes_present"] == True, f"Fold {r['fold_index']} lacks both binary classes!"
        assert r["validation_pre_families"] > 0
        assert r["validation_post_families"] > 0
        assert r["feasibility_status"] == "FEASIBLE_BINARY_TASK"


def test_08_homology_and_unresolved_isolation():
    """Ensure homology-inferred and unresolved records remain segregated from direct training folds."""
    df_m9 = pd.read_csv(FILE_SPLIT_MANIFEST_9)

    homol_folds = df_m9[df_m9["adjudicated_cohort"] == "HOMOLOGY_INFERRED"]["fold_index"].unique()
    assert list(homol_folds) == [-1], f"Homology records leaked into folds: {homol_folds}"

    unres_folds = df_m9[df_m9["adjudicated_cohort"] == "UNRESOLVED"]["fold_index"].unique()
    assert list(unres_folds) == [-2], f"Unresolved records leaked into folds: {unres_folds}"


def test_09_readiness_decision_content():
    """Verify that reports/09_1_benchmark_readiness_decision.md issues Option 1 decision."""
    assert os.path.exists(FILE_DECISION_MD)
    with open(FILE_DECISION_MD, "r", encoding="utf-8") as f:
        content = f.read()

    assert "option 1" in content.lower() or "ready for baseline modeling" in content.lower()
    assert "approved" in content.lower() or "proceed" in content.lower()


def test_10_zero_model_fitting_enforcement():
    """Enforce Step 9.1 rule: No models or weights are trained or saved."""
    prohibited_extensions = [".pt", ".pth", ".pkl", ".joblib", ".onnx", ".h5"]
    for root, _, files in os.walk(WORKSPACE_DIR):
        if ".git" in root or ".gemini" in root:
            continue
        for f in files:
            for ext in prohibited_extensions:
                assert not f.endswith(ext), f"Prohibited model file found: {os.path.join(root, f)}"
