#!/usr/bin/env python3
"""
Test Suite: Step 9 Leakage-Resistant Benchmark Design and Feasibility
====================================================================

Scientific Background and Quality Assurance:
--------------------------------------------
In computational virology, machine-learning benchmarks are frequently compromised
by sequence-level data leakage (pseudoreplication) and severe ascertainment bias:
1. Multiple clinical isolates of a single conserved viral gene (e.g. 212 variants of RS1)
   must not be split across training and validation sets.
2. The independent unit of experimental observation is the canonical gene family
   (N_families = 63 in direct evidence), not the sequence variant count (N_seq = 687).
3. The benchmark must explicitly separate direct wet-lab literature evidence from
   homology-inferred cross-species transfers.
4. Baseline hierarchies (majority class, family lookup, physicochemical, ESM-2)
   and macro-averaged metrics must be strictly prespecified before model fitting.
5. Zero models, classifiers, or hyperparameter tunings are permitted in Step 9.
"""

import os
import pytest
import pandas as pd

WORKSPACE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA_RAW = os.path.join(WORKSPACE_DIR, "data", "raw")
DATA_PROCESSED = os.path.join(WORKSPACE_DIR, "data", "processed")
REPORTS_DIR = os.path.join(WORKSPACE_DIR, "reports")

# Step 9 Deliverables
FILE_RECON_MD = os.path.join(REPORTS_DIR, "09_benchmark_reconstruction.md")
FILE_BENCH_COUNTS_CSV = os.path.join(REPORTS_DIR, "09_benchmark_counts.csv")
FILE_GENE_INDEP_CSV = os.path.join(REPORTS_DIR, "09_gene_family_independence_audit.csv")
FILE_PROV_AUDIT_CSV = os.path.join(REPORTS_DIR, "09_label_provenance_audit.csv")
FILE_FEASIBILITY_MD = os.path.join(REPORTS_DIR, "09_class_feasibility.md")
FILE_FEASIBILITY_CSV = os.path.join(REPORTS_DIR, "09_class_feasibility.csv")
FILE_SPLIT_PROTOCOL_MD = os.path.join(REPORTS_DIR, "09_split_protocol.md")
FILE_SPLIT_MANIFEST_CSV = os.path.join(DATA_PROCESSED, "09_group_split_manifest.csv")
FILE_BASELINE_MD = os.path.join(REPORTS_DIR, "09_baseline_protocol.md")
FILE_METRICS_MD = os.path.join(REPORTS_DIR, "09_evaluation_metrics.md")
FILE_SENSITIVITY_MD = os.path.join(REPORTS_DIR, "09_sensitivity_protocol.md")
FILE_READINESS_MD = os.path.join(REPORTS_DIR, "09_benchmark_readiness_report.md")

# Prior Frozen Inputs
FILE_STRICT_BENCHMARK = os.path.join(DATA_PROCESSED, "08_1_strict_temporal_benchmark.csv")
FILE_EMB_META = os.path.join(DATA_PROCESSED, "esm2_t12_representatives_metadata.csv")


def test_01_all_step9_artifacts_exist():
    """Verify that all 12 Step 9 deliverables and reports exist on disk and are non-empty."""
    required_files = [
        FILE_RECON_MD,
        FILE_BENCH_COUNTS_CSV,
        FILE_GENE_INDEP_CSV,
        FILE_PROV_AUDIT_CSV,
        FILE_FEASIBILITY_MD,
        FILE_FEASIBILITY_CSV,
        FILE_SPLIT_PROTOCOL_MD,
        FILE_SPLIT_MANIFEST_CSV,
        FILE_BASELINE_MD,
        FILE_METRICS_MD,
        FILE_SENSITIVITY_MD,
        FILE_READINESS_MD,
    ]
    for path in required_files:
        assert os.path.exists(path), f"Missing Step 9 artifact: {path}"
        assert os.path.getsize(path) > 0, f"Artifact is empty: {path}"


def test_02_raw_and_historical_inputs_preserved():
    """Verify that raw FASTA data and Steps 1-8.1 deliverables remain uncorrupted."""
    assert os.path.exists(os.path.join(DATA_RAW, "hsv1.fasta"))
    assert os.path.exists(os.path.join(DATA_RAW, "hsv2.fasta"))

    df_bench = pd.read_csv(FILE_STRICT_BENCHMARK)
    assert len(df_bench) == 3111, f"Corrupted strict benchmark count: {len(df_bench)}"

    df_meta = pd.read_csv(FILE_EMB_META)
    assert len(df_meta) == 3111, f"Corrupted representatives metadata count: {len(df_meta)}"


def test_03_strict_benchmark_exact_reconstruction():
    """
    Ensure exact dimensional reconstruction of all 3,111 representative sequences:
    687 DIRECT_EVIDENCE + 633 HOMOLOGY_INFERRED + 1,791 UNRESOLVED = 3,111.
    """
    df_manifest = pd.read_csv(FILE_SPLIT_MANIFEST_CSV)
    assert len(df_manifest) == 3111, f"Expected 3,111 records, got {len(df_manifest)}"

    cohort_counts = df_manifest["adjudicated_cohort"].value_counts().to_dict()
    assert cohort_counts.get("DIRECT_EVIDENCE", 0) == 687
    assert cohort_counts.get("HOMOLOGY_INFERRED", 0) == 633
    assert cohort_counts.get("UNRESOLVED", 0) == 1791


def test_04_independent_unit_accounting():
    """
    Verify that the independent biological units are recognized as gene families,
    with exactly 63 independent gene families in DIRECT_EVIDENCE.
    """
    df_manifest = pd.read_csv(FILE_SPLIT_MANIFEST_CSV)
    direct_sub = df_manifest[df_manifest["adjudicated_cohort"] == "DIRECT_EVIDENCE"]

    n_families = direct_sub["homology_group_id"].nunique()
    assert n_families == 63, f"Expected exactly 63 gene families in DIRECT_EVIDENCE, found {n_families}"

    # Also verify gene family independence audit
    df_indep = pd.read_csv(FILE_GENE_INDEP_CSV)
    assert len(df_indep) == 63, f"Expected 63 rows in gene independence table, found {len(df_indep)}"


def test_05_class_feasibility_power_analysis():
    """
    Verify that IE_ALPHA contains exactly 5 independent gene families and is
    explicitly audited as critically underpowered for generalizable prediction.
    """
    df_feas = pd.read_csv(FILE_FEASIBILITY_CSV)
    ie_row = df_feas[df_feas["temporal_class"] == "IE_ALPHA"].iloc[0]

    assert ie_row["independent_gene_families"] == 5, (
        f"Expected 5 gene families in IE_ALPHA, got {ie_row['independent_gene_families']}"
    )
    assert ie_row["feasibility_assessment"] == "CRITICALLY_UNDERPOWERED"


def test_06_leakage_resistant_split_manifest_integrity():
    """
    Mathematical proof of zero group leakage:
    Verify that NO gene family (homology_group_id) is split across any training and validation fold.
    """
    df_manifest = pd.read_csv(FILE_SPLIT_MANIFEST_CSV)
    direct_sub = df_manifest[df_manifest["adjudicated_cohort"] == "DIRECT_EVIDENCE"]

    # Verify 5 folds exist
    folds = sorted(direct_sub["fold_index"].unique())
    assert folds == [0, 1, 2, 3, 4], f"Expected folds [0, 1, 2, 3, 4], got {folds}"

    # Verify pairwise exclusivity of gene families across folds
    for f in range(5):
        val_families = set(direct_sub[direct_sub["fold_index"] == f]["homology_group_id"])
        train_families = set(direct_sub[direct_sub["fold_index"] != f]["homology_group_id"])

        overlap = val_families.intersection(train_families)
        assert len(overlap) == 0, (
            f"LEAKAGE DETECTED in Fold {f}! Overlapping families: {overlap}"
        )


def test_07_separation_of_direct_and_homology_inferred():
    """
    Ensure homology-inferred records are segregated into secondary sensitivity analysis
    and quarantined records are excluded from training folds.
    """
    df_manifest = pd.read_csv(FILE_SPLIT_MANIFEST_CSV)

    # Homology-inferred records must have fold_index == -1
    homol_sub = df_manifest[df_manifest["adjudicated_cohort"] == "HOMOLOGY_INFERRED"]
    assert (homol_sub["fold_index"] == -1).all(), "Homology records assigned to training folds!"

    # Unresolved records must have fold_index == -2
    unres_sub = df_manifest[df_manifest["adjudicated_cohort"] == "UNRESOLVED"]
    assert (unres_sub["fold_index"] == -2).all(), "Unresolved records assigned to training folds!"


def test_08_transcript_vs_protein_assay_distinction():
    """Verify that transcript-only genes are explicitly audited in the provenance table."""
    df_prov = pd.read_csv(FILE_PROV_AUDIT_CSV)
    assert "transcript_kinetic_class" in df_prov.columns
    assert "protein_kinetic_class" in df_prov.columns

    tx_only = df_prov[df_prov["measured_molecule"] == "TRANSCRIPT"]
    assert len(tx_only) > 0, "No transcript-only records documented"
    assert (tx_only["protein_kinetic_class"] == "NOT_DIRECTLY_MEASURED").all() or (
        tx_only["protein_kinetic_class"] == "INSUFFICIENT_EVIDENCE"
    ).any()


def test_09_prespecified_baselines_and_metrics():
    """Ensure baseline protocols and evaluation metrics are documented and macro-averaged."""
    assert os.path.exists(FILE_BASELINE_MD)
    assert os.path.exists(FILE_METRICS_MD)

    with open(FILE_BASELINE_MD, "r", encoding="utf-8") as f:
        content = f.read()
    assert "Majority-Class" in content
    assert "Gene-Family Label Lookup" in content
    assert "Physicochemical" in content
    assert "ESM-2" in content

    with open(FILE_METRICS_MD, "r", encoding="utf-8") as f:
        metrics_content = f.read()
    assert "Macro-F1" in metrics_content
    assert "Balanced Accuracy" in metrics_content


def test_10_zero_model_training_enforcement():
    """
    Enforce Step 9 mandate: No models, classifiers, weights, or tuning outputs
    may be trained or saved during benchmark design.
    """
    prohibited_extensions = [".pt", ".pth", ".pkl", ".joblib", ".onnx", ".h5"]
    for root, _, files in os.walk(WORKSPACE_DIR):
        # Ignore git and hidden directories
        if ".git" in root or ".gemini" in root:
            continue
        for f in files:
            for ext in prohibited_extensions:
                assert not f.endswith(ext), f"Prohibited model artifact found: {os.path.join(root, f)}"
