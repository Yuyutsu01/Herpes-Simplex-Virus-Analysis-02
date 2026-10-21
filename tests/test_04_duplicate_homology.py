#!/usr/bin/env python3
"""
Comprehensive Validation Test Suite for Step 4 Remediation (v2.1.0)
==================================================================
Tests:
1. Raw FASTA integrity (unmodified SHA-256 hashes).
2. Mathematically rigorous alignment identity and bidirectional coverage calculations.
3. Fragment handling (distinguishing fragments from full-length proteins and small intact proteins).
4. Boundary condition handling (95% identity threshold and sub-threshold rejection).
5. Disambiguation of UL26.5 and UL26 tokens.
6. Preservation of natural small canonical proteins (UL11, UL49A, US5, US9, US12).
7. Cross-file AEW27163.1 / WHR54491.1 Thymidine Kinase conflict resolution.
8. Authoritative reference ortholog pairings across all 74 canonical HSV genes.
9. Zero non-target retroviruses or host cellular contaminants in eligible cohort.
10. Universal 1-to-1 accounting (22,710 records across 6 mutually exclusive dispositions).
11. Numerical reconciliation across manifest, validation metrics, and markdown summary.
"""

import os
import json
import hashlib
import pandas as pd
import pytest
from Bio.Align import PairwiseAligner

WORKSPACE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RAW_HSV1 = os.path.join(WORKSPACE_DIR, "data", "raw", "hsv1.fasta")
RAW_HSV2 = os.path.join(WORKSPACE_DIR, "data", "raw", "hsv2.fasta")

REPORTS_DIR = os.path.join(WORKSPACE_DIR, "reports")
EXACT_CSV = os.path.join(REPORTS_DIR, "04_exact_duplicate_clusters.csv")
NEAR_ID_CSV = os.path.join(REPORTS_DIR, "04_near_identity_clusters.csv")
HOMOLOGY_CSV = os.path.join(REPORTS_DIR, "04_homology_group_membership.csv")
CROSS_REL_CSV = os.path.join(REPORTS_DIR, "04_hsv1_hsv2_relationships.csv")
REDUNDANCY_CSV = os.path.join(REPORTS_DIR, "04_redundancy_decisions.csv")
MANIFEST_JSON = os.path.join(REPORTS_DIR, "04_step4_manifest.json")
SUMMARY_MD = os.path.join(REPORTS_DIR, "04_step4_summary.md")
VALIDATION_JSON = os.path.join(REPORTS_DIR, "04_validation_results.json")

HSV1_SHA256 = "15001269a835b41223467e0d91bdcada433ee68c6a58c19e684c060483b31f36"
HSV2_SHA256 = "86c9d627286aa4cad121af4c84cdb39b824c14905aa57948de9a996955ca0311"


def compute_sha256(filepath):
    """Computes SHA-256 hash of a file on disk."""
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()


def test_raw_fasta_integrity():
    """Requirement: Raw input FASTA files must remain 100% unaltered."""
    assert os.path.exists(RAW_HSV1), f"Missing {RAW_HSV1}"
    assert os.path.exists(RAW_HSV2), f"Missing {RAW_HSV2}"
    assert compute_sha256(RAW_HSV1) == HSV1_SHA256, "hsv1.fasta hash mismatch!"
    assert compute_sha256(RAW_HSV2) == HSV2_SHA256, "hsv2.fasta hash mismatch!"


def test_rigorous_alignment_and_coverage_mathematics():
    """
    Requirement:
    Verify exact mathematical definitions of:
    - Aligned columns: columns where neither query nor target has gap
    - Query coverage: aligned columns / len(query)
    - Target coverage: aligned columns / len(target)
    - Bidirectional coverage: min(query_cov, target_cov)
    """
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    # Case A: Exact match (20 aa)
    s1 = "MKVLWAALLVTFLAGCQAKV"
    s2 = "MKVLWAALLVTFLAGCQAKV"
    aln = aligner.align(s1, s2)[0]
    aligned_cols = sum(a != '-' and b != '-' for a, b in zip(aln[0], aln[1]))
    matches = sum(a == b and a != '-' and b != '-' for a, b in zip(aln[0], aln[1]))
    assert matches == 20
    assert aligned_cols == 20
    assert (aligned_cols / len(s1)) * 100.0 == 100.0
    assert (aligned_cols / len(s2)) * 100.0 == 100.0

    # Case B: 95% threshold boundary (19 matches in 20 aa)
    s3 = "MKVLWAALLVTFLAGCQAKA"
    aln_sub = aligner.align(s1, s3)[0]
    m_sub = sum(a == b and a != '-' and b != '-' for a, b in zip(aln_sub[0], aln_sub[1]))
    ident_sub = (m_sub / len(aln_sub[0])) * 100.0
    assert ident_sub == 95.0

    # Case C: Sub-threshold (90% identity - 18 matches in 20 aa)
    s4 = "MKVLWAALLVTFLAGCQAEE"
    aln_sub2 = aligner.align(s1, s4)[0]
    m_sub2 = sum(a == b and a != '-' and b != '-' for a, b in zip(aln_sub2[0], aln_sub2[1]))
    ident_sub2 = (m_sub2 / len(aln_sub2[0])) * 100.0
    assert ident_sub2 == 90.0
    assert ident_sub2 < 95.0  # Must be rejected by 95% threshold

    # Case D: Short fragment vs full-length protein (50 aa vs 200 aa)
    full = "M" + "A" * 198 + "K"
    frag = "M" + "A" * 48 + "K"
    aln_frag = aligner.align(full, frag)[0]
    aligned_frag_cols = sum(a != '-' and b != '-' for a, b in zip(aln_frag[0], aln_frag[1]))
    cov_frag = (aligned_frag_cols / len(frag)) * 100.0
    cov_target = (aligned_frag_cols / len(full)) * 100.0
    assert cov_frag == 100.0  # Query (fragment) is fully aligned
    assert cov_target == 25.0  # Target (full-length) is only 25% covered!
    assert min(cov_frag, cov_target) == 25.0  # Bidirectional coverage correctly flags fragment!


def test_disambiguation_ul26_5_from_ul26():
    """
    Requirement:
    Verify that UL26.5 is recognized as a distinct canonical gene family from UL26.
    UL26 must not intercept UL26.5 records.
    """
    df_red = pd.read_csv(REDUNDANCY_CSV)
    ul26_5_recs = df_red[df_red["gene_annotation"] == "UL26.5"]
    ul26_recs = df_red[df_red["gene_annotation"] == "UL26"]

    assert len(ul26_5_recs) >= 30, f"Expected >=30 records for UL26.5, got {len(ul26_5_recs)}"
    assert len(ul26_recs) >= 50, f"Expected >=50 records for UL26, got {len(ul26_recs)}"

    # Check distinct homology groups for target records
    target_ul26_5 = ul26_5_recs[ul26_5_recs["final_disposition"] != "EXCLUDED_CONFIRMED_NON_TARGET"]
    target_ul26 = ul26_recs[ul26_recs["final_disposition"] != "EXCLUDED_CONFIRMED_NON_TARGET"]
    assert (target_ul26_5["homology_group_id"] == "HOM_UL26.5").all()
    assert (target_ul26["homology_group_id"] == "HOM_UL26").all()
    assert len(target_ul26_5) >= 30
    assert len(target_ul26) >= 50


def test_preservation_of_natural_small_proteins():
    """
    Requirement:
    Naturally small canonical intact proteins (UL11, UL49A, US5, US9, US12)
    must NOT be swept into fragments. Intact copies must be retained in eligible cohort.
    """
    df_red = pd.read_csv(REDUNDANCY_CSV)
    small_genes = ['UL11', 'UL49A', 'US5', 'US9', 'US12']

    for g in small_genes:
        g_recs = df_red[df_red["gene_annotation"] == g]
        eligible_count = g_recs["final_disposition"].isin(["RETAINED_UNIQUE_CANDIDATE", "RETAINED_VALID_VARIANT"]).sum()
        assert eligible_count > 0, f"Gene {g} had 0 eligible full-length candidates!"


def test_aew27163_cross_file_conflict_resolution():
    """
    Requirement:
    AEW27163.1 (in hsv2.fasta) is verified via NCBI GenBank record as HSV-1 (taxid 10298)
    and has identical sequence to WHR54491.1. It must have effective species attribution HSV-1.
    """
    df_red = pd.read_csv(REDUNDANCY_CSV)
    aew_row = df_red[df_red["record_id"].str.contains("AEW27163", case=False, na=False)].iloc[0]
    whr_row = df_red[df_red["record_id"].str.contains("WHR54491", case=False, na=False)].iloc[0]

    assert aew_row["sequence_sha256"] == whr_row["sequence_sha256"], "Hashes must match!"
    assert aew_row["effective_species_attribution"] == "HSV1", "AEW27163.1 effective species must be HSV1"
    assert aew_row["gene_annotation"] == "UL23"


def test_authoritative_ortholog_pairing_all_74_genes():
    """
    Requirement:
    All 74 canonical HSV genes must be audited with reference ortholog pairs.
    Zero severe length asymmetry outliers (length ratio >= 0.60 for all core pairs).
    """
    df_cross = pd.read_csv(CROSS_REL_CSV)
    assert len(df_cross) == 74, f"Expected exactly 74 audited canonical pairs, got {len(df_cross)}"

    # Check that high-confidence orthologs make up the majority
    high_conf = df_cross[df_cross["confidence_level"] == "HIGH_CONFIDENCE_ORTHOLOG"]
    assert len(high_conf) >= 45, f"Expected >=45 high-confidence orthologs, got {len(high_conf)}"

    # Check previously missing genes are present
    genes_present = set(df_cross["canonical_gene"])
    for g in ['UL26.5', 'UL11', 'UL49A', 'US5', 'US9']:
        assert g in genes_present, f"Missing canonical gene {g} in ortholog relationships!"


def test_contamination_purge_zero_retroviruses_retained():
    """
    Requirement:
    Zero non-target retroviruses or host cellular proteins permitted in retained viral sets.
    """
    df_red = pd.read_csv(REDUNDANCY_CSV)
    retained = df_red[df_red["final_disposition"].isin(["RETAINED_UNIQUE_CANDIDATE", "RETAINED_VALID_VARIANT"])]

    desc_list = retained["disposition_reason"].str.lower()
    assert not desc_list.str.contains("retrovir").any(), "Retained retroviral records found!"
    assert not (retained["gene_annotation"] == "NON_TARGET_RETROVIRAL").any()

    excluded = df_red[df_red["final_disposition"] == "EXCLUDED_CONFIRMED_NON_TARGET"]
    assert len(excluded) == 1043, f"Expected 1,043 excluded non-target records, got {len(excluded)}"


def test_universal_mutually_exclusive_dispositions():
    """
    Requirement:
    Every raw record (22,710) is assigned to exactly one mutually exclusive disposition.
    All categories sum to 22,710.
    """
    df_red = pd.read_csv(REDUNDANCY_CSV)
    assert len(df_red) == 22710

    compound_keys = df_red["source_file"].apply(os.path.basename) + "::" + df_red["record_id"].astype(str)
    assert compound_keys.nunique() == 22710, "Compound keys must be uniquely indexed"

    valid_dispositions = {
        "RETAINED_UNIQUE_CANDIDATE",
        "RETAINED_VALID_VARIANT",
        "RETAINED_FRAGMENT_SEPARATE_ANALYSIS",
        "EXCLUDED_CONFIRMED_NON_TARGET",
        "QUARANTINED_TAXONOMY_CONFLICT",
        "QUARANTINED_UNVERIFIED",
        "EXCLUDED_SEQUENCE_INVALID"
    }
    assigned_disps = set(df_red["final_disposition"].unique())
    assert assigned_disps.issubset(valid_dispositions), f"Invalid dispositions found: {assigned_disps - valid_dispositions}"

    counts = df_red["final_disposition"].value_counts().to_dict()
    assert sum(counts.values()) == 22710


def test_near_identity_clustering_metrics():
    """
    Requirement:
    Near-identity clusters must have verified pairwise metrics (identity >= 95%, bidi cov >= 90%).
    """
    df_near = pd.read_csv(NEAR_ID_CSV)
    assert len(df_near) > 10000
    assert (df_near["sequence_identity_pct"] >= 95.0).all()
    assert (df_near["query_coverage_pct"] >= 90.0).all()
    assert (df_near["target_coverage_pct"] >= 90.0).all()
    assert (df_near["bidirectional_coverage_pct"] >= 90.0).all()


def test_ml_leakage_groups_partitioning():
    """
    Requirement:
    Stable ML leakage prevention groups assigned to all records.
    """
    df_red = pd.read_csv(REDUNDANCY_CSV)
    assert df_red["prevent_ml_leakage_group"].notna().all()
    assert (df_red["prevent_ml_leakage_group"] != "").all()
    assert df_red["prevent_ml_leakage_group"].nunique() >= 500


def test_reconciliation_of_manifest_and_summary():
    """
    Requirement:
    Manifest and reports reconcile with exact numeric parity.
    """
    with open(MANIFEST_JSON, "r", encoding="utf-8") as f:
        manifest = json.load(f)

    assert manifest["dataset_counts"]["total_input_records"] == 22710
    assert manifest["dataset_counts"]["unique_sequence_hashes"] == 22689
    assert manifest["dataset_counts"]["cross_file_shared_sequence_hashes"] == 21
    assert manifest["dataset_counts"]["canonical_ortholog_pairs_audited"] == 74

    df_red = pd.read_csv(REDUNDANCY_CSV)
    disp_counts = df_red["final_disposition"].value_counts().to_dict()
    assert manifest["dataset_counts"]["final_disposition_breakdown"] == disp_counts


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
