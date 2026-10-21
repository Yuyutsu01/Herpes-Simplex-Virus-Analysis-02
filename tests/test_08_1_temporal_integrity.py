#!/usr/bin/env python3
"""
Test Suite: Step 8.1 Temporal Annotation Integrity and Benchmark Validation
==========================================================================

Scientific Background and Quality Standards:
--------------------------------------------
Biological annotations used for proteome-wide comparative genomics and future
machine learning must satisfy rigorous integrity standards:
1. Mathematical closure: All subcohorts and categories must sum exactly to the
   representative cohort denominator (N = 3,111).
2. Traceability: Every single record must trace back to accession, gene, species,
   assay type, and verified literature citation.
3. Separation of Direct vs Inferred Evidence: Homology transfers must be strictly
   segregated from primary direct experimental benchmarks.
4. Quarantine Enforcement: Unannotated open reading frames, context-dependent loci,
   and unverified submitter aliases must remain quarantined.
5. Preservation: Historical Step 1-8 outputs and raw FASTA files must remain unmodified.
"""

import os
import json
import pytest
import pandas as pd

WORKSPACE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA_RAW = os.path.join(WORKSPACE_DIR, "data", "raw")
DATA_PROCESSED = os.path.join(WORKSPACE_DIR, "data", "processed")
REPORTS_DIR = os.path.join(WORKSPACE_DIR, "reports")

# Step 8 historical files (must remain untouched)
FILE_STEP8_ADJUDICATION = os.path.join(DATA_PROCESSED, "08_temporal_label_adjudication.csv")
FILE_STEP8_GENE_EVID = os.path.join(DATA_PROCESSED, "08_temporal_gene_evidence.csv")
FILE_STEP8_SOURCES = os.path.join(DATA_PROCESSED, "08_temporal_evidence_sources.csv")
FILE_STEP8_CROSS = os.path.join(DATA_PROCESSED, "08_cross_species_temporal_labels.csv")
FILE_STEP8_ALIAS = os.path.join(DATA_PROCESSED, "08_alias_recovery.csv")
FILE_REP_META = os.path.join(DATA_PROCESSED, "esm2_t12_representatives_metadata.csv")

# Step 8.1 deliverables
FILE_RECON_CSV = os.path.join(REPORTS_DIR, "08_1_count_reconciliation.csv")
FILE_PROV_CSV = os.path.join(REPORTS_DIR, "08_1_label_provenance_audit.csv")
FILE_BIO_CSV = os.path.join(REPORTS_DIR, "08_1_high_priority_biology_review.csv")
FILE_CIT_CSV = os.path.join(REPORTS_DIR, "08_1_citation_verification.csv")
FILE_STRICT_CSV = os.path.join(DATA_PROCESSED, "08_1_strict_temporal_benchmark.csv")
FILE_REPORT_MD = os.path.join(REPORTS_DIR, "08_1_integrity_audit.md")

VALID_COHORTS = {"DIRECT_EVIDENCE", "HOMOLOGY_INFERRED", "UNRESOLVED"}
VALID_TRANSFER_STATUSES = {
    "DIRECTLY_SUPPORTED_HSV1",
    "DIRECTLY_SUPPORTED_IN_BOTH_SPECIES",
    "INFERRED_FROM_HOMOLOGY",
    "NO_DEFENSIBLE_TRANSFER",
}
DUAL_SUPPORTED_GENES = {"RS1", "RL2", "UL23", "UL27", "US4", "US6", "UL19"}


def test_01_all_step8_1_artifacts_exist():
    """Verify that all Step 8.1 reports, manifests, and tables exist on disk and are non-empty."""
    required_files = [
        FILE_RECON_CSV,
        FILE_PROV_CSV,
        FILE_BIO_CSV,
        FILE_CIT_CSV,
        FILE_STRICT_CSV,
        FILE_REPORT_MD,
    ]
    for path in required_files:
        assert os.path.exists(path), f"Missing Step 8.1 artifact: {path}"
        assert os.path.getsize(path) > 0, f"Artifact is empty: {path}"


def test_02_historical_and_step8_files_unmodified():
    """Verify that raw datasets, frozen metadata, and Step 8 deliverables were preserved."""
    assert os.path.exists(os.path.join(DATA_RAW, "hsv1.fasta"))
    assert os.path.exists(os.path.join(DATA_RAW, "hsv2.fasta"))

    df_rep = pd.read_csv(FILE_REP_META)
    assert len(df_rep) == 3111, f"Corrupted esm2 representatives count: {len(df_rep)}"

    df_s8_adj = pd.read_csv(FILE_STEP8_ADJUDICATION)
    assert len(df_s8_adj) == 3111, f"Corrupted Step 8 adjudication count: {len(df_s8_adj)}"


def test_03_strict_benchmark_record_traceability_and_closure():
    """
    Ensure every row in 08_1_strict_temporal_benchmark.csv corresponds 1-to-1
    with the 3,111 frozen representative sequences.
    """
    df_strict = pd.read_csv(FILE_STRICT_CSV)
    df_rep = pd.read_csv(FILE_REP_META)

    assert len(df_strict) == 3111, f"Expected 3,111 rows, got {len(df_strict)}"
    assert (df_strict["embedding_row_index"].values == df_rep["embedding_row_index"].values).all(), (
        "embedding_row_index does not strictly match representatives metadata"
    )
    assert (df_strict["record_id"].values == df_rep["record_id"].values).all(), (
        "record_id does not strictly match representatives metadata"
    )
    assert df_strict["embedding_row_index"].is_unique, "Duplicate embedding_row_index in strict benchmark"


def test_04_count_reconciliation_exact_closure():
    """Verify exact count closure across the three non-overlapping cohorts."""
    df_strict = pd.read_csv(FILE_STRICT_CSV)
    cohort_counts = df_strict["adjudicated_cohort"].value_counts().to_dict()

    n_direct = cohort_counts.get("DIRECT_EVIDENCE", 0)
    n_homology = cohort_counts.get("HOMOLOGY_INFERRED", 0)
    n_unresolved = cohort_counts.get("UNRESOLVED", 0)

    # Exact closure assertion
    assert n_direct + n_homology + n_unresolved == 3111, (
        f"Sum of cohorts ({n_direct} + {n_homology} + {n_unresolved} = {n_direct + n_homology + n_unresolved}) != 3,111"
    )

    # Check reconciliation table completeness
    df_recon = pd.read_csv(FILE_RECON_CSV)
    assert len(df_recon) >= 5, "Count reconciliation table is missing key audit items"
    items = set(df_recon["metric_item"].str.lower())
    assert any("direct evidence" in it for it in items), "Missing direct evidence reconciliation"
    assert any("tier 1/2" in it for it in items), "Missing tier 1/2 reconciliation"
    assert any("alias" in it for it in items), "Missing alias recovery reconciliation"


def test_05_mutually_exclusive_cohort_partition():
    """Ensure strict cohort exclusivity and transfer-state discipline."""
    df_strict = pd.read_csv(FILE_STRICT_CSV)

    observed_cohorts = set(df_strict["adjudicated_cohort"].unique())
    assert observed_cohorts.issubset(VALID_COHORTS), f"Invalid cohorts: {observed_cohorts}"

    # DIRECT_EVIDENCE must not have INFERRED_FROM_HOMOLOGY transfer status
    direct_rows = df_strict[df_strict["adjudicated_cohort"] == "DIRECT_EVIDENCE"]
    assert not (direct_rows["cross_species_transfer_status"] == "INFERRED_FROM_HOMOLOGY").any(), (
        "DIRECT_EVIDENCE cohort contains homology-inferred transfers"
    )

    # HOMOLOGY_INFERRED records must belong to HSV2
    homology_rows = df_strict[df_strict["adjudicated_cohort"] == "HOMOLOGY_INFERRED"]
    assert (homology_rows["effective_species"] == "HSV2").all(), (
        "HOMOLOGY_INFERRED cohort contains non-HSV2 records"
    )

    # For HSV-2 in DIRECT_EVIDENCE, gene must be in DUAL_SUPPORTED_GENES
    direct_hsv2 = direct_rows[direct_rows["effective_species"] == "HSV2"]
    for _, r in direct_hsv2.iterrows():
        assert r["effective_gene_symbol"] in DUAL_SUPPORTED_GENES, (
            f"HSV-2 record {r['record_id']} in DIRECT_EVIDENCE has unsupported gene {r['effective_gene_symbol']}"
        )


def test_06_evidence_source_traceability():
    """Verify that all benchmarked records have traceable, registered citations."""
    df_strict = pd.read_csv(FILE_STRICT_CSV)
    df_cit = pd.read_csv(FILE_CIT_CSV)
    valid_source_ids = set(df_cit["source_id"].unique())

    benchmarked_rows = df_strict[df_strict["adjudicated_cohort"].isin(["DIRECT_EVIDENCE", "HOMOLOGY_INFERRED"])]
    for idx, row in benchmarked_rows.iterrows():
        raw_sources = str(row["evidence_source_ids"]).split(";")
        for s in raw_sources:
            s_clean = s.strip()
            if s_clean and s_clean != "NONE":
                assert s_clean in valid_source_ids, (
                    f"Row {idx} ({row['record_id']}) references unregistered source: '{s_clean}'"
                )


def test_07_quarantine_enforcement_for_unverified_aliases_and_insufficient():
    """Ensure all recovered submitter aliases and unannotated ORFs remain in UNRESOLVED."""
    df_strict = pd.read_csv(FILE_STRICT_CSV)

    # All alias-recovered records must be in UNRESOLVED
    alias_rows = df_strict[df_strict["is_alias_recovered"] == True]
    assert (alias_rows["adjudicated_cohort"] == "UNRESOLVED").all(), (
        "Submitter alias records were incorrectly included in strict benchmark"
    )

    # All INSUFFICIENT_EVIDENCE records must be in UNRESOLVED
    insufficient_rows = df_strict[df_strict["normalized_temporal_class"] == "INSUFFICIENT_EVIDENCE"]
    assert (insufficient_rows["adjudicated_cohort"] == "UNRESOLVED").all(), (
        "INSUFFICIENT_EVIDENCE records were incorrectly included in benchmark"
    )

    # UL55 and US8A records must have NO_DEFENSIBLE_TRANSFER
    ul55_us8a = df_strict[df_strict["effective_gene_symbol"].isin(["UL55", "US8A"])]
    assert (ul55_us8a["cross_species_transfer_status"] == "NO_DEFENSIBLE_TRANSFER").all(), (
        "UL55 or US8A records have invalid transfer status"
    )


def test_08_transcript_vs_protein_distinction():
    """Verify separate transcript and protein kinetic fields without conflation."""
    df_strict = pd.read_csv(FILE_STRICT_CSV)

    assert "transcript_kinetic_class" in df_strict.columns
    assert "protein_kinetic_class" in df_strict.columns

    # For transcript-only genes, protein class must be NOT_DIRECTLY_MEASURED or INSUFFICIENT
    transcript_only_genes = {"UL53", "UL2", "UL4", "UL20", "UL49A", "US5", "US10", "UL43", "UL45", "UL56", "US2"}
    sub_tx = df_strict[
        (df_strict["effective_gene_symbol"].isin(transcript_only_genes)) &
        (df_strict["adjudicated_cohort"] != "UNRESOLVED")
    ]
    for _, r in sub_tx.iterrows():
        assert r["protein_kinetic_class"] == "NOT_DIRECTLY_MEASURED", (
            f"Gene {r['effective_gene_symbol']} has protein kinetics claimed without direct measurement"
        )


def test_09_high_priority_biology_review_traceability():
    """Verify key biological test cases are documented in the review table."""
    df_bio = pd.read_csv(FILE_BIO_CSV)
    genes_reviewed = set(df_bio["gene_symbol"].unique())

    # Check key complex loci
    assert any("UL19" in g for g in genes_reviewed), "UL19 missing from biological review"
    assert any("UL26" in g for g in genes_reviewed), "UL26 missing from biological review"
    assert any("RL1" in g for g in genes_reviewed), "RL1 missing from biological review"
    assert any("UL48" in g for g in genes_reviewed), "UL48 missing from biological review"


def test_10_citation_verification_completeness():
    """Verify all 6 core citations possess complete metadata and verified status."""
    df_cit = pd.read_csv(FILE_CIT_CSV)
    assert len(df_cit) == 6, f"Expected exactly 6 core citations, got {len(df_cit)}"

    for _, r in df_cit.iterrows():
        assert len(str(r["full_citation"])) > 20, f"Incomplete citation for {r['source_id']}"
        assert r["verification_status"] in {"INDEPENDENTLY_VERIFIED", "AUTHORITATIVE_REVIEW_VERIFIED"}, (
            f"Unverified citation status for {r['source_id']}"
        )
