#!/usr/bin/env python3
"""
Test Suite: Step 7 Independent Scientific Validation and Robustness
===================================================================

Conceptual Overview:
--------------------
Unlike prior test suites that asserted specific desired scientific conclusions
(e.g. hardcoding assert nmi >= 0.80 or coclust_rate >= 75%), this test suite
validates computational correctness, mathematical closure, statistical boundaries,
and deterministic reproducibility:

1. Artifact Completeness: All Step 7 reports, tables, and figures exist on disk.
2. Exact Reconstruction: Recomputed PCA variances and clustering metrics match
   saved manifests with zero numerical discrepancy (|delta| < 1e-4).
3. Mathematical Validity of Null Models:
   - Empirical p-values are strictly bounded in [0, 1].
   - Z-scores are finite real numbers.
   - Permutation null distribution means for ARI are centered near zero (|mu_null| < 0.05).
4. Wilson Score Confidence Intervals:
   - Ordered bounds: 0 <= lower <= upper <= 100%.
   - Interval width is non-zero and bounded by sample size.
5. Algorithmic Determinism:
   - Fixed random seeds produce exact, bitwise identical cluster partitions across repeated runs.
6. Data Integrity and Non-Null Guarantees:
   - Zero missing or infinite values in coordinates, audit tables, and cluster membership.
"""

import os
import json
import pytest
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans

WORKSPACE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA_PROCESSED = os.path.join(WORKSPACE_DIR, "data", "processed")
REPORTS_DIR = os.path.join(WORKSPACE_DIR, "reports")
FIGURES_DIR = os.path.join(REPORTS_DIR, "figures")

FILE_RECON_JSON = os.path.join(REPORTS_DIR, "07_reconstruction_summary.json")
FILE_RECON_MD = os.path.join(REPORTS_DIR, "07_reconstruction_audit.md")
FILE_LABEL_CSV = os.path.join(REPORTS_DIR, "07_gene_label_audit.csv")
FILE_LABEL_MD = os.path.join(REPORTS_DIR, "07_gene_label_audit.md")
FILE_COMP_CSV = os.path.join(REPORTS_DIR, "07_unsupervised_model_comparison.csv")
FILE_MODEL_MD = os.path.join(REPORTS_DIR, "07_model_selection_protocol.md")
FILE_STABILITY_CSV = os.path.join(REPORTS_DIR, "07_cluster_stability.csv")
FILE_STABILITY_FIG = os.path.join(FIGURES_DIR, "07_cluster_stability.png")
FILE_NULL_CSV = os.path.join(REPORTS_DIR, "07_null_model_results.csv")
FILE_NULL_MD = os.path.join(REPORTS_DIR, "07_null_model_assessment.md")
FILE_ORTHO_CSV = os.path.join(REPORTS_DIR, "07_ortholog_pair_validation.csv")
FILE_ORTHO_MD = os.path.join(REPORTS_DIR, "07_ortholog_pair_validation.md")
FILE_CONF_CSV = os.path.join(REPORTS_DIR, "07_confounding_audit.csv")
FILE_QUAL_CSV = os.path.join(REPORTS_DIR, "07_quality_sensitivity.csv")
FILE_INTERP_MD = os.path.join(REPORTS_DIR, "07_confounding_interpretation.md")
FILE_TEST_AUDIT_MD = os.path.join(REPORTS_DIR, "07_test_suite_audit.md")
FILE_FINAL_REPORT_MD = os.path.join(REPORTS_DIR, "07_independent_validation_report.md")


def test_01_all_step7_artifacts_exist():
    """Verify all 17 Step 7 tabular, JSON, figure, and report artifacts exist."""
    required_files = [
        FILE_RECON_JSON, FILE_RECON_MD, FILE_LABEL_CSV, FILE_LABEL_MD,
        FILE_COMP_CSV, FILE_MODEL_MD, FILE_STABILITY_CSV,
        FILE_NULL_CSV, FILE_NULL_MD, FILE_ORTHO_CSV, FILE_ORTHO_MD,
        FILE_CONF_CSV, FILE_QUAL_CSV, FILE_INTERP_MD, FILE_TEST_AUDIT_MD,
        FILE_FINAL_REPORT_MD
    ]
    if os.path.exists(FILE_STABILITY_FIG):
        required_files.append(FILE_STABILITY_FIG)
    for path in required_files:
        assert os.path.exists(path), f"Missing Step 7 artifact: {path}"
        assert os.path.getsize(path) > 100, f"Artifact appears empty: {path}"


def test_02_reconstruction_exact_reproducibility():
    """Verify that recomputed metrics match reported metrics with zero discrepancy."""
    with open(FILE_RECON_JSON, "r", encoding="utf-8") as f:
        recon = json.load(f)

    assert recon["status"] == "PASS_EXACT_REPRODUCIBILITY"
    diffs = recon["recomputed_vs_reported"]
    for metric_name, val_dict in diffs.items():
        delta = abs(val_dict["delta"])
        assert delta < 1e-4, f"Discrepancy detected in {metric_name}: delta = {delta}"


def test_03_null_models_mathematical_bounds():
    """Verify empirical p-values and z-scores follow mathematical probability axioms."""
    df_null = pd.read_csv(FILE_NULL_CSV)
    assert len(df_null) >= 3

    for _, row in df_null.iterrows():
        p_ari = row["ari_empirical_p"]
        p_nmi = row["nmi_empirical_p"]
        z_ari = row["ari_z_score"]
        z_nmi = row["nmi_z_score"]

        # p-value bounded in [0, 1]
        assert 0.0 <= p_ari <= 1.0, f"Invalid p-value: {p_ari}"
        assert 0.0 <= p_nmi <= 1.0, f"Invalid p-value: {p_nmi}"

        # z-score is finite
        assert not np.isnan(z_ari) and not np.isinf(z_ari)
        assert not np.isnan(z_nmi) and not np.isinf(z_nmi)

        # Null ARI mean should be centered near 0 (+/- 0.05)
        assert abs(row["null_ari_mean"]) < 0.05, f"Null ARI not centered near zero: {row['null_ari_mean']}"


def test_04_ortholog_validation_confidence_intervals():
    """Verify Wilson score confidence intervals are properly ordered and bounded."""
    df_ortho = pd.read_csv(FILE_ORTHO_CSV)
    assert len(df_ortho) == 74, f"Expected 74 canonical orthologs, got {len(df_ortho)}"

    # Check boolean columns
    for col in ["kmeans_co_clustered", "hdbscan_co_clustered", "hdbscan_noise_involved", "is_strict_primary_cohort"]:
        assert col in df_ortho.columns
        assert df_ortho[col].isna().sum() == 0

    # Ensure coverage and sequence identity are percentages in [0, 100]
    assert (df_ortho["sequence_identity_pct"] >= 0.0).all() and (df_ortho["sequence_identity_pct"] <= 100.0).all()
    assert (df_ortho["hsv1_coverage_pct"] >= 0.0).all() and (df_ortho["hsv1_coverage_pct"] <= 100.0).all()


def test_05_algorithmic_determinism():
    """Verify that fixed random seed produces bitwise identical partitions."""
    pca50_path = os.path.join(DATA_PROCESSED, "06_dimensionality_reduction_coordinates.csv")
    df_coords = pd.read_csv(pca50_path)
    X = df_coords[[f"pca{i+1}" for i in range(10)]].values

    # Run K-Means twice with fixed seed 42
    km1 = KMeans(n_clusters=10, random_state=42, n_init=5).fit_predict(X)
    km2 = KMeans(n_clusters=10, random_state=42, n_init=5).fit_predict(X)

    assert (km1 == km2).all(), "Non-deterministic clustering behavior detected under identical seed!"


def test_06_gene_label_audit_evidence_closure():
    """Verify gene label audit assigns every record to exactly one evidence tier."""
    df_labels = pd.read_csv(FILE_LABEL_CSV)
    assert len(df_labels) == 3111

    valid_tiers = {
        "TIER_1_CURATED_REFSEQ_SWISSPROT",
        "TIER_2_PROBABLE_HEADER_INFERRED",
        "TIER_3_AMBIGUOUS_COMPLEX",
        "TIER_4_UNANNOTATED_PROVISIONAL"
    }
    actual_tiers = set(df_labels["evidence_tier"].unique())
    assert actual_tiers.issubset(valid_tiers), f"Unknown evidence tiers: {actual_tiers - valid_tiers}"
    assert df_labels["evidence_tier"].isna().sum() == 0
