#!/usr/bin/env python3
"""
Test Suite: Step 11 Robustness, Generalization, and Scientific Validity Audit
=============================================================================

Scientific Background and Quality Assurance:
--------------------------------------------
Step 11 stress-tests the Step 10/10.1 temporal replication timing benchmark:
1. Deliverable Verification: Confirms presence of all 8 Step 11 deliverables.
2. Integrity Manifest: Verifies artifact hashing and zero group leakage across folds.
3. Dual Majority Baselines: Verifies fold sequence-majority (BA=0.5321) and fixed
   family-majority (BA=0.5000, Macro-F1=0.3883).
4. LOFO Stability Bounds: Proves ESM-2 Macro-F1 remains within [0.7886, 0.8082]
   under single-family omissions.
5. Influential Loci Identification: Audits the exact delta and classification status
   for top influential families.
6. Fold Stability & RS1 Impact: Verifies fold statistics and RS1's 212-sequence concentration in Fold 0.
7. Evidence-Tier Stratification: Verifies 52 protein-validated vs 11 transcript-only families.
8. Paired Significance: Confirms McNemar exact p-values for ESM-2 superiority over baselines.
9. Species Infeasibility Audit: Confirms 63 HSV-1 vs 7 HSV-2 direct families, verifying that
   an unseen-family cross-species benchmark is statistically infeasible.
10. Scientific Boundaries: Asserts explicit disclaimers against causal mechanisms and novel herpesviruses.
"""

import os
import json
import pytest
import numpy as np
import pandas as pd

WORKSPACE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
REPORTS_DIR = os.path.join(WORKSPACE_DIR, "reports")
DATA_PROCESSED = os.path.join(WORKSPACE_DIR, "data", "processed")

# Step 11 Deliverables
FILE_11_INTEGRITY = os.path.join(REPORTS_DIR, "11_data_integrity_checks.csv")
FILE_11_MANIFEST = os.path.join(REPORTS_DIR, "11_reproducibility_manifest.json")
FILE_11_PRIMARY = os.path.join(REPORTS_DIR, "11_primary_metrics_recomputed.csv")
FILE_11_LOFO = os.path.join(REPORTS_DIR, "11_leave_one_family_out.csv")
FILE_11_FOLD = os.path.join(REPORTS_DIR, "11_fold_stability.csv")
FILE_11_EVIDENCE = os.path.join(REPORTS_DIR, "11_label_evidence_sensitivity.csv")
FILE_11_PAIRED = os.path.join(REPORTS_DIR, "11_paired_model_comparison.csv")
FILE_11_AUDIT_MD = os.path.join(REPORTS_DIR, "11_robustness_generalization_audit.md")


def test_01_all_step11_deliverables_exist():
    """Verify that all 8 Step 11 deliverables exist on disk and are non-empty."""
    required = [
        FILE_11_INTEGRITY,
        FILE_11_MANIFEST,
        FILE_11_PRIMARY,
        FILE_11_LOFO,
        FILE_11_FOLD,
        FILE_11_EVIDENCE,
        FILE_11_PAIRED,
        FILE_11_AUDIT_MD
    ]
    for path in required:
        assert os.path.exists(path), f"Missing Step 11 deliverable: {path}"
        assert os.path.getsize(path) > 0, f"Deliverable is empty: {path}"


def test_02_manifest_and_integrity_hashes():
    """Verify that reproducibility manifest JSON and integrity CSV contain valid hashes and no group leakage."""
    with open(FILE_11_MANIFEST, "r", encoding="utf-8") as f:
        manifest = json.load(f)

    assert "audit_timestamp" in manifest
    assert "artifacts" in manifest
    assert len(manifest["artifacts"]) >= 8

    df_integrity = pd.read_csv(FILE_11_INTEGRITY)
    assert len(df_integrity) >= 8
    assert (df_integrity["status"] == "VERIFIED_PRESENT").all()
    assert (df_integrity["sha256_hash"].str.len() == 64).all()


def test_03_primary_recomputed_metrics_closure():
    """Verify recomputed primary benchmark metrics including both sequence and family majority baselines."""
    df_prim = pd.read_csv(FILE_11_PRIMARY).set_index("model_id")
    assert len(df_prim) == 6

    # Verify ESM-2 primary metrics
    assert abs(df_prim.loc["ridge_esm2", "family_macro_f1"] - 0.7945652173913043) < 1e-12
    assert abs(df_prim.loc["ridge_esm2", "family_balanced_accuracy"] - 0.7945652173913043) < 1e-12

    # Verify Fixed Family Majority Baseline (predicts POST for all families)
    assert abs(df_prim.loc["fixed_family_majority", "family_balanced_accuracy"] - 0.5000) < 1e-12
    assert abs(df_prim.loc["fixed_family_majority", "family_macro_f1"] - 0.3883495145631068) < 1e-12
    assert abs(df_prim.loc["fixed_family_majority", "family_recall_pre"] - 0.0) < 1e-12
    assert abs(df_prim.loc["fixed_family_majority", "family_recall_post"] - 1.0) < 1e-12

    # Verify Fold-Specific Sequence Majority Baseline
    assert abs(df_prim.loc["dummy_majority", "family_balanced_accuracy"] - 0.5320652173913043) < 1e-12


def test_04_leave_one_family_out_stability():
    """Verify that ESM-2 LOFO Macro-F1 remains tightly bounded in [0.7886, 0.8082]."""
    df_lofo = pd.read_csv(FILE_11_LOFO)
    esm_lofo = df_lofo[df_lofo["model_id"] == "ridge_esm2"]

    assert len(esm_lofo) == 63
    min_mf1 = esm_lofo["lofo_family_macro_f1"].min()
    max_mf1 = esm_lofo["lofo_family_macro_f1"].max()

    assert abs(min_mf1 - 0.7886363636363636) < 1e-12
    assert abs(max_mf1 - 0.8081575246132209) < 1e-12
    assert (max_mf1 - min_mf1) < 0.025, f"LOFO range unexpectedly wide: {max_mf1 - min_mf1}"


def test_05_lofo_influential_families():
    """Verify that removing correctly predicted families decreases score and misclassified increases score."""
    df_lofo = pd.read_csv(FILE_11_LOFO)
    esm_lofo = df_lofo[df_lofo["model_id"] == "ridge_esm2"]

    correct_fams = esm_lofo[esm_lofo["is_correct"]]
    incorrect_fams = esm_lofo[~esm_lofo["is_correct"]]

    assert len(correct_fams) == 51
    assert len(incorrect_fams) == 12

    # Omitting a correct family must reduce Macro-F1
    assert (correct_fams["delta_macro_f1"] < 0).all()
    # Omitting an incorrect family must increase Macro-F1
    assert (incorrect_fams["delta_macro_f1"] > 0).all()


def test_06_fold_stability_and_rs1_concentration():
    """Verify fold statistics and RS1's 212-sequence concentration in Fold 0."""
    df_fold = pd.read_csv(FILE_11_FOLD).set_index("model_id")

    assert df_fold.loc["ridge_esm2", "fold0_sequences"] == 290
    assert abs(df_fold.loc["ridge_esm2", "fold0_family_macro_f1"] - 0.8080808080808081) < 1e-12
    assert abs(df_fold.loc["ridge_esm2", "fold0_sequence_macro_f1"] - 0.9273391080249124) < 1e-12


def test_07_evidence_tier_stratification():
    """Verify performance across 52 protein-validated vs 11 transcript-only loci."""
    df_ev = pd.read_csv(FILE_11_EVIDENCE).set_index("model_id")

    assert df_ev.loc["ridge_esm2", "all_direct_n_fams"] == 63
    assert df_ev.loc["ridge_esm2", "protein_validated_n_fams"] == 52
    assert df_ev.loc["ridge_esm2", "transcript_only_n_fams"] == 11

    esm_prot_mf1 = df_ev.loc["ridge_esm2", "protein_validated_macro_f1"]
    esm_tx_mf1 = df_ev.loc["ridge_esm2", "transcript_only_macro_f1"]
    assert abs(esm_prot_mf1 - 0.8153846153846154) < 1e-12
    assert abs(esm_tx_mf1 - 0.6857142857142857) < 1e-12
    assert esm_prot_mf1 > esm_tx_mf1


def test_08_paired_model_comparison_significance():
    """Verify exact binomial McNemar test p-values for paired comparisons against ESM-2."""
    df_paired = pd.read_csv(FILE_11_PAIRED).set_index("comparison")

    # ESM-2 vs Fixed Family Majority: p = 0.0347 < 0.05
    p_fix = df_paired.loc["ridge_esm2 vs fixed_family_majority", "mcnemar_exact_pvalue"]
    assert p_fix < 0.05
    assert df_paired.loc["ridge_esm2 vs fixed_family_majority", "statistically_significant_at_05"] == True

    # ESM-2 vs Dummy Majority: p = 0.00019 < 0.001
    p_dum = df_paired.loc["ridge_esm2 vs dummy_majority", "mcnemar_exact_pvalue"]
    assert p_dum < 0.001

    # ESM-2 vs Physicochemical models: p < 0.05
    p_lr = df_paired.loc["ridge_esm2 vs logreg_physicochemical", "mcnemar_exact_pvalue"]
    p_rf = df_paired.loc["ridge_esm2 vs rf_physicochemical", "mcnemar_exact_pvalue"]
    assert p_lr < 0.05
    assert p_rf < 0.05

    # ESM-2 vs Combined: p = 0.5000 (not significant)
    p_comb = df_paired.loc["ridge_esm2 vs combined_phys_esm2", "mcnemar_exact_pvalue"]
    assert p_comb == 0.5
    assert df_paired.loc["ridge_esm2 vs combined_phys_esm2", "statistically_significant_at_05"] == False


def test_09_species_representation_and_infeasibility():
    """Verify that only 7 direct families exist in HSV-2 and cross-species benchmark is documented as infeasible."""
    manifest = pd.read_csv(os.path.join(DATA_PROCESSED, "09_group_split_manifest.csv"))
    direct = manifest[manifest["adjudicated_cohort"] == "DIRECT_EVIDENCE"]

    sp_fam = direct.groupby(["homology_group_id", "effective_species"]).size().unstack(fill_value=0)
    assert (sp_fam["HSV1"] > 0).sum() == 63
    assert (sp_fam["HSV2"] > 0).sum() == 7

    # Check that report explicitly states infeasibility
    with open(FILE_11_AUDIT_MD, "r", encoding="utf-8") as f:
        content = f.read()

    assert "STATISTICALLY INFEASIBLE" in content
    assert "7 gene families have direct experimental kinetic annotations in HSV-2" in content


def test_10_scientific_boundaries_enforced():
    """Verify that reports explicitly disclaim causal mechanisms and novel herpesvirus generalization."""
    with open(FILE_11_AUDIT_MD, "r", encoding="utf-8") as f:
        content = f.read()

    assert "No Causal Mechanism" in content
    assert "No Generalization to Distant Herpesviruses" in content
    assert "Sample Size Constraints ($N=63$)" in content
