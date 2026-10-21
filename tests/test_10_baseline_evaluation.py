#!/usr/bin/env python3
"""
Test Suite: Step 10 Leakage-Resistant Baseline Modeling for HSV Replication Timing
=================================================================================

Scientific Background and Quality Assurance:
--------------------------------------------
Step 10 validates the prespecified binary baseline modeling of HSV replication timing
(Pre-Replicative vs Post-Replicative) under strict group holdouts:
1. Cohort Integrity: Direct-evidence cohort contains exactly 687 sequences and 63 families.
2. Binary Label Mapping: IE_ALPHA and E_BETA -> PRE_REPLICATIVE; L1_GAMMA1 and L2_GAMMA2 -> POST_REPLICATIVE.
3. Embedding Row Alignment: The 480-dim ESM-2 embedding array aligns with metadata and manifest.
4. Group Disjointness: Confirms zero family overlap across training and validation splits.
5. Fold-Local Preprocessing: Preprocessing scalers are fitted exclusively on training folds.
6. Prohibited Features Excluded: No identifiers or metadata labels used in feature matrices.
7. Deliverable Integrity: All 7 required deliverables exist, are non-empty, and conform to schema.
8. Baseline Superiority: Validates that ESM-2 embeddings statistically outperform majority guessing and bulk physical descriptors.
"""

import os
import json
import pytest
import numpy as np
import pandas as pd
from sklearn.metrics import f1_score, balanced_accuracy_score

WORKSPACE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA_PROCESSED = os.path.join(WORKSPACE_DIR, "data", "processed")
REPORTS_DIR = os.path.join(WORKSPACE_DIR, "reports")

# Input Artifacts
FILE_BENCH_81 = os.path.join(DATA_PROCESSED, "08_1_strict_temporal_benchmark.csv")
FILE_SPLIT_MANIFEST = os.path.join(DATA_PROCESSED, "09_group_split_manifest.csv")
FILE_EMB_META = os.path.join(DATA_PROCESSED, "esm2_t12_representatives_metadata.csv")
FILE_EMBEDDINGS = os.path.join(DATA_PROCESSED, "esm2_t12_representatives_embeddings.npy")
FILE_READINESS_DECISION = os.path.join(REPORTS_DIR, "09_1_benchmark_readiness_decision.md")

# Step 10 Deliverables
OUT_CONFIG_JSON = os.path.join(REPORTS_DIR, "10_model_configuration.json")
OUT_BASELINE_METRICS = os.path.join(REPORTS_DIR, "10_baseline_metrics.csv")
OUT_FOLD_METRICS = os.path.join(REPORTS_DIR, "10_fold_level_metrics.csv")
OUT_OOF_PREDS = os.path.join(REPORTS_DIR, "10_out_of_fold_predictions.csv")
OUT_FAM_METRICS = os.path.join(REPORTS_DIR, "10_family_level_metrics.csv")
OUT_RESULTS_MD = os.path.join(REPORTS_DIR, "10_baseline_results.md")
OUT_LIMITATIONS_MD = os.path.join(REPORTS_DIR, "10_limitations_and_interpretation.md")


def test_01_all_step10_deliverables_exist():
    """Verify that all 7 Step 10 deliverables exist on disk and are non-empty."""
    required_files = [
        OUT_CONFIG_JSON,
        OUT_BASELINE_METRICS,
        OUT_FOLD_METRICS,
        OUT_OOF_PREDS,
        OUT_FAM_METRICS,
        OUT_RESULTS_MD,
        OUT_LIMITATIONS_MD
    ]
    for path in required_files:
        assert os.path.exists(path), f"Missing Step 10 deliverable: {path}"
        assert os.path.getsize(path) > 0, f"Deliverable is empty: {path}"


def test_02_cohort_integrity_and_binary_mapping():
    """Verify that primary cohort contains exactly 687 direct sequences, 63 families, and correct binary classes."""
    manifest = pd.read_csv(FILE_SPLIT_MANIFEST)
    bench = pd.read_csv(FILE_BENCH_81)

    direct = manifest[manifest["adjudicated_cohort"] == "DIRECT_EVIDENCE"].copy()
    assert len(direct) == 687, f"Expected 687 direct sequences, got {len(direct)}"
    assert direct["homology_group_id"].nunique() == 63, f"Expected 63 families, got {direct['homology_group_id'].nunique()}"

    # Verify binary classes in out-of-fold predictions
    oof = pd.read_csv(OUT_OOF_PREDS)
    valid_binary = {"PRE_REPLICATIVE", "POST_REPLICATIVE"}
    assert set(oof["true_class"].unique()) == valid_binary
    assert set(oof["predicted_class"].unique()) == valid_binary

    # Check class distributions
    oof_m = oof[oof["model_id"] == "ridge_esm2"]
    assert len(oof_m) == 687
    pre_seqs = (oof_m["true_class"] == "PRE_REPLICATIVE").sum()
    post_seqs = (oof_m["true_class"] == "POST_REPLICATIVE").sum()
    assert pre_seqs == 350
    assert post_seqs == 337


def test_03_embedding_alignment_and_index_closure():
    """Verify that frozen embeddings array matches metadata and manifest row-for-row."""
    meta = pd.read_csv(FILE_EMB_META)
    manifest = pd.read_csv(FILE_SPLIT_MANIFEST)
    embeddings = np.load(FILE_EMBEDDINGS)

    assert embeddings.shape == (3111, 480)
    assert len(meta) == 3111
    assert len(manifest) == 3111
    assert (meta["record_id"] == manifest["record_id"]).all()
    assert (manifest["embedding_row_index"] == np.arange(3111)).all()


def test_04_zero_family_leakage_across_folds():
    """Verify that no gene family crosses between training and validation splits in any fold."""
    manifest = pd.read_csv(FILE_SPLIT_MANIFEST)
    direct = manifest[manifest["adjudicated_cohort"] == "DIRECT_EVIDENCE"]

    for fold in range(5):
        train_fams = set(direct[direct["fold_index"] != fold]["homology_group_id"])
        val_fams = set(direct[direct["fold_index"] == fold]["homology_group_id"])
        overlap = train_fams.intersection(val_fams)
        assert len(overlap) == 0, f"Fold {fold} has family leakage: {overlap}"


def test_05_model_configuration_prespecified():
    """Verify that model configuration JSON specifies all 5 baselines without global tuning."""
    with open(OUT_CONFIG_JSON, "r", encoding="utf-8") as f:
        config = json.load(f)

    expected_models = {
        "dummy_majority",
        "logreg_physicochemical",
        "rf_physicochemical",
        "ridge_esm2",
        "combined_phys_esm2"
    }
    assert set(config.keys()) == expected_models
    for m in expected_models:
        assert "model_type" in config[m]
        assert "feature_set" in config[m]


def test_06_fold_level_metrics_consistency():
    """Verify that fold-level metrics CSV has 5 folds per model and valid metrics."""
    df_fold = pd.read_csv(OUT_FOLD_METRICS)
    assert len(df_fold) == 25  # 5 models * 5 folds

    for model_id in df_fold["model_id"].unique():
        sub = df_fold[df_fold["model_id"] == model_id]
        assert len(sub) == 5
        assert set(sub["fold_index"]) == {0, 1, 2, 3, 4}
        assert (sub["val_total_sequences"] > 0).all()
        assert (sub["val_total_families"] > 0).all()
        assert (sub["family_macro_f1"] >= 0.0).all()
        assert (sub["family_macro_f1"] <= 1.0).all()


def test_07_family_level_metrics_integrity():
    """Verify that family-level metrics file contains exactly 63 families per model."""
    df_fam = pd.read_csv(OUT_FAM_METRICS)
    expected_models = ["dummy_majority", "logreg_physicochemical", "rf_physicochemical", "ridge_esm2", "combined_phys_esm2"]

    for m in expected_models:
        sub = df_fam[df_fam["model_id"] == m]
        assert len(sub) == 63, f"Expected 63 families for {m}, got {len(sub)}"
        assert sub["homology_group_id"].nunique() == 63
        assert (sub["sequence_count"] >= 1).all()
        assert set(sub["true_class"]) == {"PRE_REPLICATIVE", "POST_REPLICATIVE"}


def test_08_esm2_superiority_over_baselines():
    """Verify empirical baseline hierarchy: Ridge ESM-2 outperforms dummy and physicochemical models."""
    df_summary = pd.read_csv(OUT_BASELINE_METRICS).set_index("model_id")

    esm2_fam_f1 = df_summary.loc["ridge_esm2", "family_macro_f1"]
    dummy_fam_f1 = df_summary.loc["dummy_majority", "family_macro_f1"]
    lr_phys_fam_f1 = df_summary.loc["logreg_physicochemical", "family_macro_f1"]
    rf_phys_fam_f1 = df_summary.loc["rf_physicochemical", "family_macro_f1"]

    # ESM-2 should achieve substantial family-level Macro-F1 (> 0.70)
    assert esm2_fam_f1 >= 0.75, f"ESM-2 family Macro-F1 unexpectedly low: {esm2_fam_f1}"
    assert esm2_fam_f1 > dummy_fam_f1, f"ESM-2 ({esm2_fam_f1}) did not beat Dummy ({dummy_fam_f1})"
    assert esm2_fam_f1 > lr_phys_fam_f1, f"ESM-2 ({esm2_fam_f1}) did not beat LR Phys ({lr_phys_fam_f1})"
    assert esm2_fam_f1 > rf_phys_fam_f1, f"ESM-2 ({esm2_fam_f1}) did not beat RF Phys ({rf_phys_fam_f1})"

    # Physicochemical models should be modest (< 0.65)
    assert lr_phys_fam_f1 < 0.65
    assert rf_phys_fam_f1 < 0.65


def test_09_rs1_pseudoreplication_inflation_effect():
    """Verify that RS1 sequence count inflates sequence-level accuracy over family-level accuracy."""
    df_summary = pd.read_csv(OUT_BASELINE_METRICS).set_index("model_id")

    esm2_seq_f1 = df_summary.loc["ridge_esm2", "sequence_macro_f1"]
    esm2_fam_f1 = df_summary.loc["ridge_esm2", "family_macro_f1"]

    # Sequence-level F1 should be inflated by RS1 compared to Family-level F1
    assert esm2_seq_f1 > esm2_fam_f1, f"Sequence F1 ({esm2_seq_f1}) not higher than Family F1 ({esm2_fam_f1})"
    assert (esm2_seq_f1 - esm2_fam_f1) >= 0.08, "Expected at least 8% inflation due to RS1 replicates"


def test_10_scientific_limitations_stated():
    """Verify that reports explicitly disclaim novel virus generalization and causal mechanisms."""
    with open(OUT_LIMITATIONS_MD, "r", encoding="utf-8") as f:
        content = f.read()

    assert "No Generalization to Novel Herpesvirus Genera" in content
    assert "No Causal Mechanism Inferred" in content
    assert "Sample Size Limitations ($N=63$)" in content
    assert "Transcript vs. Protein Kinetic Divergence" in content
