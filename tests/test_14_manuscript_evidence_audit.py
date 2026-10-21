"""
Unit and integration tests for Step 14: Final Manuscript Evidence Audit,
Scientific Claim Verification, and Publication Readiness.

Verifies:
1. Deliverable presence, non-emptiness, and format integrity.
2. Manuscript inventory and explicit documentation of standalone draft status.
3. Claim-to-evidence ledger coverage, claim types, and audit verdicts.
4. Metric reconciliation (reported vs recomputed exact match for 0.7946, 0.8095).
5. Temporal label evidence audit (52 protein-validated vs 11 transcript-only families).
6. Reference verification (100% bibliographic verification with valid DOIs/PMIDs).
7. Unsupervised and geometry claims (effective rank 13.6, PCA held-out collapse).
8. Statistical and causal language audit (Holm-adjusted p=0.0602 requirement).
9. Audited manuscript draft consistency with underlying data deliverables.
10. Author review memo completeness.
"""

from pathlib import Path
import json
import numpy as np
import pandas as pd
import pytest

BASE_DIR = Path(__file__).resolve().parent.parent
REPORTS_DIR = BASE_DIR / "reports"

INV_CSV = REPORTS_DIR / "14_manuscript_inventory.csv"
LEDGER_CSV = REPORTS_DIR / "14_claim_to_evidence_ledger.csv"
METRIC_RECON_CSV = REPORTS_DIR / "14_metric_reconciliation.csv"
EVIDENCE_AUDIT_CSV = REPORTS_DIR / "14_temporal_label_evidence_audit.csv"
REF_VERIF_CSV = REPORTS_DIR / "14_reference_verification.csv"
UNSUP_CSV = REPORTS_DIR / "14_unsupervised_claim_audit.csv"
LANG_CSV = REPORTS_DIR / "14_statistical_language_audit.csv"
CHANGE_LOG_CSV = REPORTS_DIR / "14_manuscript_change_log.csv"
ASSESSMENT_MD = REPORTS_DIR / "14_publication_readiness_assessment.md"
AUDITED_DRAFT_MD = REPORTS_DIR / "14_audited_manuscript_draft.md"
AUTHOR_MEMO_MD = REPORTS_DIR / "14_unresolved_claims_author_review.md"
MANIFEST_JSON = REPORTS_DIR / "14_reproducibility_manifest.json"


def test_deliverables_exist_and_nonempty():
    """Verify that all 12 Step 14 deliverables exist and are non-empty."""
    required = [
        INV_CSV,
        LEDGER_CSV,
        METRIC_RECON_CSV,
        EVIDENCE_AUDIT_CSV,
        REF_VERIF_CSV,
        UNSUP_CSV,
        LANG_CSV,
        CHANGE_LOG_CSV,
        ASSESSMENT_MD,
        AUDITED_DRAFT_MD,
        AUTHOR_MEMO_MD,
        MANIFEST_JSON,
    ]
    for p in required:
        assert p.exists(), f"Missing required deliverable: {p.name}"
        assert p.stat().st_size > 0, f"Deliverable is empty: {p.name}"


def test_manuscript_inventory_status():
    """Verify inventory records and clear documentation of standalone draft status."""
    df_inv = pd.read_csv(INV_CSV)
    assert len(df_inv) >= 50, f"Expected >= 50 inventoried files, got {len(df_inv)}"
    
    # Verify missing standalone status is documented
    missing_entry = df_inv[df_inv["file_category"] == "MISSING_STANDALONE_DRAFT"]
    assert len(missing_entry) == 1, "Must explicitly document standalone manuscript status"
    assert "reports/14_audited_manuscript_draft.md" in str(missing_entry.iloc[0]["notes"])


def test_claim_ledger_structure_and_verdicts():
    """Verify that claims span all core sections and contain rigorous epistemological verdicts."""
    df_ledger = pd.read_csv(LEDGER_CSV)
    assert len(df_ledger) >= 20, f"Expected >= 20 claims, got {len(df_ledger)}"

    expected_sections = [
        "Data Curation & Sequence Inventory",
        "Representation Analysis & Geometry",
        "Temporal Benchmark & Adjudication",
        "Supervised Predictive Modeling",
        "Statistical Inference & Multiplicity",
        "Biological Stratification",
        "Generalization & Cross-Species Claims"
    ]
    for sec in expected_sections:
        matches = df_ledger[df_ledger["manuscript_section"] == sec]
        assert len(matches) > 0, f"Missing claims for section: {sec}"

    # Verify presence of key verdicts
    verdicts = set(df_ledger["epistemological_verdict"])
    assert "VERIFIED FROM DATA" in verdicts
    assert "INCONCLUSIVE" in verdicts
    assert "CONTRADICTED BY EVIDENCE" in verdicts
    assert "UNSUPPORTED" in verdicts


def test_metric_reconciliation_exactness():
    """Verify reported metrics match primary recomputed outputs down to decimal tolerance."""
    df_recon = pd.read_csv(METRIC_RECON_CSV)
    
    # ESM-2 family balanced accuracy & macro-F1 exact match
    bal_acc_row = df_recon[df_recon["metric_name"] == "ESM-2 Family Balanced Accuracy"].iloc[0]
    assert bal_acc_row["absolute_difference"] == 0.0
    assert bal_acc_row["reconciliation_verdict"] == "CONFIRMED_EXACT_MATCH"

    macro_f1_row = df_recon[df_recon["metric_name"] == "ESM-2 Family Macro-F1"].iloc[0]
    assert macro_f1_row["absolute_difference"] == 0.0
    assert macro_f1_row["reconciliation_verdict"] == "CONFIRMED_EXACT_MATCH"

    # Multiplicity correction flag for physicochemical models
    rf_p_row = df_recon[df_recon["metric_name"] == "ESM-2 vs Physicochemical RF p-value"].iloc[0]
    assert rf_p_row["reconciliation_verdict"] == "REQUIRES_MULTIPLICITY_CORRECTION"
    assert "0.0602" in str(rf_p_row["notes"])


def test_temporal_evidence_audit_subgroups():
    """Verify temporal evidence audit captures 52 protein-validated and 11 transcript-only families."""
    df_ev = pd.read_csv(EVIDENCE_AUDIT_CSV)
    assert len(df_ev) == 63, "Must audit all 63 gene families"

    tier_counts = df_ev["evidence_tier"].value_counts().to_dict()
    assert tier_counts.get("DIRECT_PROTEIN_VALIDATED") == 52
    assert tier_counts.get("TRANSCRIPT_ONLY_ANNOTATED") == 11


def test_reference_verification_completeness():
    """Verify that all cited references have verified DOIs and PMIDs."""
    df_refs = pd.read_csv(REF_VERIF_CSV)
    assert len(df_refs) >= 7

    for _, r in df_refs.iterrows():
        assert r["bibliographic_metadata_status"] == "VERIFIED_EXTERNAL_DOI_PMID"
        assert r["claim_level_evidence_status"] == "VERIFIED_SUPPORTS_SPECIFIC_CLAIM"
        assert pd.notnull(r["doi"])
        assert pd.notnull(r["pmid"]) or "Book Chapter" in str(r["pmid"])


def test_unsupervised_claims_and_effective_rank():
    """Verify effective rank 13.6 and PCA component variance retention."""
    df_unsup = pd.read_csv(UNSUP_CSV)

    rank_row = df_unsup[df_unsup["analysis_component"].str.contains("Effective Rank")].iloc[0]
    assert rank_row["observed_value"] == "13.6" or rank_row["observed_value"] == 13.6
    assert rank_row["claim_status"] == "VERIFIED FROM DATA"


def test_statistical_language_audit_recommendations():
    """Verify that causal/overstated terms have specific required revisions."""
    df_lang = pd.read_csv(LANG_CSV)
    assert len(df_lang) >= 5

    sig_row = df_lang[df_lang["target_term"] == "significant"].iloc[0]
    assert "0.0602" in sig_row["recommended_revision"]

    ext_row = df_lang[df_lang["target_term"] == "independent external validation"].iloc[0]
    assert "ortholog sequence transfer diagnostic" in ext_row["recommended_revision"]


def test_audited_manuscript_draft_consistency():
    """Verify that reports/14_audited_manuscript_draft.md contains reconciled and corrected values."""
    text = AUDITED_DRAFT_MD.read_text(encoding="utf-8")

    # Correct RS1 record count
    assert "212" in text
    # Correct combined feature dimension
    assert "488" in text
    # Correct family balanced accuracy
    assert "79.46%" in text
    # Correct Holm-adjusted p-value
    assert "0.0602" in text
    # Reclassification of HSV-2 transfer
    assert "ortholog sequence transferability" in text


def test_author_review_memo_decisions():
    """Verify that reports/14_unresolved_claims_author_review.md outlines the 4 critical decisions."""
    memo = AUTHOR_MEMO_MD.read_text(encoding="utf-8")
    assert "Decision 1" in memo
    assert "Decision 2" in memo
    assert "Decision 3" in memo
    assert "Decision 4" in memo
    assert "0.0602" in memo
    assert "ortholog sequence transfer diagnostic" in memo.lower()
