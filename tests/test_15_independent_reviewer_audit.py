"""
Test Suite: test_15_independent_reviewer_audit.py
Purpose: Validate Step 15 independent adversarial peer review deliverables,
  metric rechecks, claim classifications, temporal label audit, and revised manuscript.
"""

from pathlib import Path
import pytest
import pandas as pd
import numpy as np

BASE_DIR = Path(__file__).resolve().parent.parent
REPORTS_DIR = BASE_DIR / "reports"
SCRIPTS_DIR = BASE_DIR / "scripts"


def test_step15_deliverables_exist():
    """Verify that all required Step 15 reports and artifacts exist and are non-empty."""
    required_files = [
        REPORTS_DIR / "15_independent_reviewer_report.md",
        REPORTS_DIR / "15_independent_claim_verification.csv",
        REPORTS_DIR / "15_metric_independent_recheck.csv",
        REPORTS_DIR / "15_temporal_label_adversarial_audit.csv",
        REPORTS_DIR / "15_reference_claim_support.csv",
        REPORTS_DIR / "15_manuscript_completeness_checklist.csv",
        REPORTS_DIR / "15_required_revision_log.csv",
        REPORTS_DIR / "15_submission_readiness_decision.md",
        REPORTS_DIR / "15_revised_manuscript.md",
        SCRIPTS_DIR / "15_independent_reviewer_audit.py",
    ]

    for filepath in required_files:
        assert filepath.exists(), f"Missing required deliverable: {filepath}"
        assert filepath.stat().st_size > 0, f"Deliverable is empty: {filepath}"


def test_metric_independent_recheck():
    """Verify that independent metric recalculations match saved reports with zero discrepancy."""
    recheck_path = REPORTS_DIR / "15_metric_independent_recheck.csv"
    df = pd.read_csv(recheck_path)

    # Check ESM-2 family metrics
    ridge_row = df[df["model_id"] == "ridge_esm2"].iloc[0]
    assert np.isclose(ridge_row["independently_rechecked_bal_acc"], 0.7946, atol=1e-4)
    assert np.isclose(ridge_row["independently_rechecked_macro_f1"], 0.7946, atol=1e-4)
    assert np.isclose(ridge_row["independently_rechecked_accuracy"], 0.8095, atol=1e-4)
    assert ridge_row["discrepancy"] == 0.0
    assert ridge_row["audit_verdict"] == "CONFIRMED_EXACT"

    # Verify sequence level distortion by RS1
    assert np.isclose(ridge_row["sequence_level_macro_f1"], 0.9112, atol=1e-4)

    # Check protein vs transcript subgroups
    assert np.isclose(ridge_row["protein_subgroup_f1_n52"], 0.8154, atol=1e-4)
    assert np.isclose(ridge_row["transcript_subgroup_f1_n11"], 0.6857, atol=1e-4)

    # Verify statistical comparison against majority baseline
    mcnemar_maj = df[df["model_id"] == "ridge_esm2_vs_dummy_majority"].iloc[0]
    assert np.isclose(mcnemar_maj["independently_rechecked_macro_f1"], 0.000192, atol=1e-5)

    # Verify unadjusted p-values against physicochemical models
    mcnemar_rf = df[df["model_id"] == "ridge_esm2_vs_rf_physicochemical"].iloc[0]
    assert np.isclose(mcnemar_rf["independently_rechecked_macro_f1"], 0.02006, atol=1e-4)
    mcnemar_lr = df[df["model_id"] == "ridge_esm2_vs_logreg_physicochemical"].iloc[0]
    assert np.isclose(mcnemar_lr["independently_rechecked_macro_f1"], 0.02945, atol=1e-4)


def test_adversarial_claim_verification():
    """Audit scientific claims to ensure overclaims are appropriately classified and qualified."""
    claims_path = REPORTS_DIR / "15_independent_claim_verification.csv"
    df = pd.read_csv(claims_path)

    assert len(df) == 7

    # ADV-001: Binary framing
    c1 = df[df["claim_id"] == "ADV-001"].iloc[0]
    assert c1["evidentiary_classification"] == "VERIFIED_DIRECTLY_DEMONSTRATED"
    assert "binary" in c1["required_manuscript_qualification"].lower()

    # ADV-003: Physicochemical comparison misses strict Holm correction
    c3 = df[df["claim_id"] == "ADV-003"].iloc[0]
    assert c3["evidentiary_classification"] == "UNPROVEN_OR_BORDERLINE_SIGNIFICANCE"

    # ADV-005: Low-rank PCA held-out collapse
    c5 = df[df["claim_id"] == "ADV-005"].iloc[0]
    assert c5["evidentiary_classification"] == "CONTRADICTED_IN_SAMPLE_OVERFITTING"

    # ADV-006: Sister-species ortholog transfer
    c6 = df[df["claim_id"] == "ADV-006"].iloc[0]
    assert c6["evidentiary_classification"] == "OVERCLAIM_RECLASSIFIED_AS_ORTHOLOG_TRANSFER"

    # ADV-007: External clinical validation
    c7 = df[df["claim_id"] == "ADV-007"].iloc[0]
    assert c7["evidentiary_classification"] == "MISSING_EXTERNAL_DATASET"


def test_adversarial_temporal_label_audit():
    """Verify row-level adversarial audit across all 63 direct benchmark families."""
    audit_path = REPORTS_DIR / "15_temporal_label_adversarial_audit.csv"
    df = pd.read_csv(audit_path)

    # Exactly 63 gene families
    assert len(df) == 63

    # Total sequence variants must equal 687
    assert df["sequence_variant_count"].sum() == 687

    # Protein validated vs transcript only split
    prot_count = (df["measured_endpoint"] == "TRANSCRIPT_AND_PROTEIN").sum()
    tx_count = (df["measured_endpoint"] == "TRANSCRIPT").sum()
    assert prot_count == 52
    assert tx_count == 11

    # Check high vs moderate confidence disposition
    assert (df["label_confidence"] == "HIGH_CONFIDENCE").sum() == 52
    assert (df["label_confidence"] == "MODERATE_CONFIDENCE_QUALIFIED").sum() == 11

    # RS1 family check (212 variants)
    rs1_row = df[df["homology_group_id"] == "HOM_RS1"].iloc[0]
    assert rs1_row["sequence_variant_count"] == 212
    assert rs1_row["gene_symbol"] == "RS1"


def test_reference_claim_support():
    """Audit bibliographic reference verification and claim-level support."""
    ref_path = REPORTS_DIR / "15_reference_claim_support.csv"
    df = pd.read_csv(ref_path)

    assert len(df) == 7
    # All 7 must be verified with primary assay audited
    assert df["metadata_verified"].all()
    assert (df["claim_level_verification_status"] == "VERIFIED_PRIMARY_EVIDENCE_ACCESSIBLE").all()


def test_manuscript_completeness_and_revisions():
    """Verify that manuscript completeness checklist and revision log are comprehensive."""
    check_path = REPORTS_DIR / "15_manuscript_completeness_checklist.csv"
    df_check = pd.read_csv(check_path)
    assert len(df_check) >= 15

    rev_path = REPORTS_DIR / "15_required_revision_log.csv"
    df_rev = pd.read_csv(rev_path)
    assert len(df_rev) >= 6

    # Verify essential revisions are logged
    rev_ids = df_rev["revision_id"].tolist()
    assert "REV-001" in rev_ids  # Binary framing
    assert "REV-003" in rev_ids  # Formatted tables
    assert "REV-006" in rev_ids  # Conclusion


def test_revised_manuscript_structure_and_integrity():
    """Verify that the Step 15 revised manuscript incorporates required corrections and structures."""
    ms_path = REPORTS_DIR / "15_revised_manuscript.md"
    assert ms_path.exists()
    content = ms_path.read_text(encoding="utf-8")

    # Required sections
    assert "## Abstract" in content
    assert "## 1. Introduction" in content
    assert "## 2. Materials and Methods" in content
    assert "## 3. Results" in content
    assert "## 4. Discussion and Limitations" in content
    assert "## 5. Conclusion" in content
    assert "## References" in content

    # Required tables and figure callouts
    assert "### Table 1" in content
    assert "### Table 2" in content
    assert "### Table 3" in content
    assert "Figure 1" in content
    assert "Figure 2" in content
    assert "Figure 3" in content

    # Required epistemological qualifiers
    assert "binary temporal classification task" in content.lower()
    assert "0.0602" in content  # Holm adjusted p-value
    assert "ortholog sequence transfer diagnostic" in content.lower()
    assert "non-linear" in content.lower()


def test_submission_readiness_decision():
    """Verify that the submission readiness decision is unambiguous and separates co-author review."""
    decision_path = REPORTS_DIR / "15_submission_readiness_decision.md"
    content = decision_path.read_text(encoding="utf-8")

    # Must choose exactly option 1
    assert "1. READY FOR CO-AUTHOR REVIEW" in content

    # Must declare separation from journal submission readiness
    assert "NOT equivalent to journal submission readiness" in content
