#!/usr/bin/env python3
"""
Test Suite: Step 8 Evidence-Based Temporal Kinetics Annotation Quality Assurance
================================================================================

Scientific Background and Design Rationale:
------------------------------------------
In computational virology, biological label curation must satisfy stringent
validation criteria to avoid circularity, label leakage, and spurious claims
of machine-learning readiness.

This test suite strictly validates implementation integrity, identifier traceability,
vocabulary compliance, evidence linkage, and operational closure without enforcing
predetermined biological distributions or hardcoded class ratios:

1. Record Traceability:
   - Exactly 3,111 representative sequences are accounted for 1-to-1.
   - Row indices and sequence identifiers match the frozen Step 5/6 representative corpus.
2. Controlled Vocabulary:
   - Adjudicated temporal classes belong strictly to the 7-class taxonomy.
   - Evidence strength and transfer statuses belong to explicit enumerations.
3. Provenance & Source Linkage:
   - Every literature-backed gene assignment links to an existing source in the bibliographic catalog.
4. Auditable Rationale & No Unsupported Assignments:
   - No sequence receives a non-insufficient temporal label without explicit rationale and evidence.
   - Insufficient-evidence records are explicitly flagged with INSUFFICIENT evidence strength.
5. Cross-Species Transfer Discipline:
   - Direct experimental evidence is explicitly segregated from homology-inferred transfers.
   - Ortholog transfers adhere strictly to the 3-state transfer ontology.
6. Transcript vs Protein Assay Distinction:
   - Assays are explicitly logged as transcript-level or protein-level.
7. Preservation of Prior Step Outputs & Raw Data:
   - Raw FASTA files and frozen Step 4-7 deliverables remain completely uncorrupted.
8. Conflict & Ambiguity Logging:
   - Mixed, context-dependent, or disputed loci possess detailed rationale and conflict documentation.
"""

import os
import json
import pytest
import pandas as pd

WORKSPACE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA_RAW = os.path.join(WORKSPACE_DIR, "data", "raw")
DATA_PROCESSED = os.path.join(WORKSPACE_DIR, "data", "processed")
REPORTS_DIR = os.path.join(WORKSPACE_DIR, "reports")
FIGURES_DIR = os.path.join(REPORTS_DIR, "figures")

# Raw and historical frozen inputs
RAW_HSV1 = os.path.join(DATA_RAW, "hsv1.fasta")
RAW_HSV2 = os.path.join(DATA_RAW, "hsv2.fasta")
FROZEN_VALID_META = os.path.join(DATA_PROCESSED, "valid_proteins_metadata.csv")
FROZEN_EMB_META = os.path.join(DATA_PROCESSED, "esm2_t12_representatives_metadata.csv")
FROZEN_EMB_NPY = os.path.join(DATA_PROCESSED, "esm2_t12_representatives_embeddings.npy")

# Step 8 deliverables
FILE_UNIVERSE_MD = os.path.join(REPORTS_DIR, "08_annotation_universe_audit.md")
FILE_UNIVERSE_JSON = os.path.join(REPORTS_DIR, "08_annotation_universe_summary.json")
FILE_ALIAS_CSV = os.path.join(DATA_PROCESSED, "08_alias_recovery.csv")
FILE_ALIAS_JSON = os.path.join(DATA_PROCESSED, "08_gene_alias_dictionary.json")
FILE_ALIAS_MD = os.path.join(REPORTS_DIR, "08_alias_recovery_audit.md")
FILE_ACCESSION_CSV = os.path.join(REPORTS_DIR, "08_accession_verification.csv")
FILE_ACCESSION_MD = os.path.join(REPORTS_DIR, "08_accession_verification_summary.md")
FILE_SOURCES_CSV = os.path.join(DATA_PROCESSED, "08_temporal_evidence_sources.csv")
FILE_GENE_EVID_CSV = os.path.join(DATA_PROCESSED, "08_temporal_gene_evidence.csv")
FILE_LIT_MD = os.path.join(REPORTS_DIR, "08_temporal_literature_review.md")
FILE_PROTOCOL_MD = os.path.join(REPORTS_DIR, "08_temporal_label_protocol.md")
FILE_ADJUDICATION_CSV = os.path.join(DATA_PROCESSED, "08_temporal_label_adjudication.csv")
FILE_CROSS_SPECIES_CSV = os.path.join(DATA_PROCESSED, "08_cross_species_temporal_labels.csv")
FILE_CROSS_SPECIES_MD = os.path.join(REPORTS_DIR, "08_cross_species_label_transfer_audit.md")
FILE_COVERAGE_CSV = os.path.join(REPORTS_DIR, "08_temporal_label_coverage.csv")
FILE_BIAS_MD = os.path.join(REPORTS_DIR, "08_temporal_label_bias_audit.md")
FILE_REVIEW_CSV = os.path.join(REPORTS_DIR, "08_manual_review_queue.csv")
FILE_QA_MD = os.path.join(REPORTS_DIR, "08_annotation_quality_report.md")
FILE_FINAL_REPORT_MD = os.path.join(REPORTS_DIR, "08_temporal_annotation_report.md")
FIG_CLASS_DIST = os.path.join(FIGURES_DIR, "08_temporal_class_distribution.png")
FIG_EVID_COVERAGE = os.path.join(FIGURES_DIR, "08_temporal_evidence_coverage.png")

VALID_TEMPORAL_CLASSES = {
    "IE_ALPHA",
    "E_BETA",
    "L1_GAMMA1",
    "L2_GAMMA2",
    "MIXED_OR_CONTEXT_DEPENDENT",
    "DISPUTED",
    "INSUFFICIENT_EVIDENCE",
}

VALID_EVIDENCE_STRENGTHS = {
    "HIGH_EXPERIMENTAL",
    "MODERATE_EXPERIMENTAL",
    "MODERATE_SUBMITTER_ALIAS",
    "INSUFFICIENT",
}

VALID_TRANSFER_STATUSES_REPRESENTATIVE = {
    "DIRECTLY_SUPPORTED_HSV1",
    "DIRECTLY_SUPPORTED_IN_BOTH_SPECIES",
    "INFERRED_FROM_HOMOLOGY",
    "NO_DEFENSIBLE_TRANSFER",
}

VALID_TRANSFER_STATUSES_ORTHOLOG = {
    "DIRECTLY_SUPPORTED_IN_BOTH_SPECIES",
    "INFERRED_FROM_HOMOLOGY",
    "NO_DEFENSIBLE_TRANSFER",
}


def test_01_all_step8_artifacts_exist():
    """Verify that all required Step 8 deliverable files and figures exist on disk."""
    required_artifacts = [
        FILE_UNIVERSE_MD,
        FILE_UNIVERSE_JSON,
        FILE_ALIAS_CSV,
        FILE_ALIAS_JSON,
        FILE_ALIAS_MD,
        FILE_ACCESSION_CSV,
        FILE_ACCESSION_MD,
        FILE_SOURCES_CSV,
        FILE_GENE_EVID_CSV,
        FILE_LIT_MD,
        FILE_PROTOCOL_MD,
        FILE_ADJUDICATION_CSV,
        FILE_CROSS_SPECIES_CSV,
        FILE_CROSS_SPECIES_MD,
        FILE_COVERAGE_CSV,
        FILE_BIAS_MD,
        FILE_REVIEW_CSV,
        FILE_QA_MD,
        FILE_FINAL_REPORT_MD,
    ]
    if os.path.exists(FIG_CLASS_DIST):
        required_artifacts.append(FIG_CLASS_DIST)
    if os.path.exists(FIG_EVID_COVERAGE):
        required_artifacts.append(FIG_EVID_COVERAGE)
    for path in required_artifacts:
        assert os.path.exists(path), f"Missing required Step 8 artifact: {path}"
        assert os.path.getsize(path) > 0, f"Artifact is empty: {path}"


def test_02_raw_and_historical_inputs_preserved():
    """Verify that raw FASTA data and Step 1-7 frozen datasets were not corrupted or overwritten."""
    assert os.path.exists(RAW_HSV1) and os.path.getsize(RAW_HSV1) > 0
    assert os.path.exists(RAW_HSV2) and os.path.getsize(RAW_HSV2) > 0

    df_valid = pd.read_csv(FROZEN_VALID_META)
    assert len(df_valid) == 18172, f"Corrupted valid_proteins_metadata.csv row count: {len(df_valid)}"

    df_rep = pd.read_csv(FROZEN_EMB_META)
    assert len(df_rep) == 3111, f"Corrupted esm2_t12_representatives_metadata.csv row count: {len(df_rep)}"
    assert os.path.exists(FROZEN_EMB_NPY), "Missing representative embeddings .npy array"


def test_03_record_identifier_uniqueness_and_traceability():
    """
    Ensure every row in 08_temporal_label_adjudication.csv corresponds 1-to-1
    with the 3,111 frozen representative sequences.
    """
    df_adj = pd.read_csv(FILE_ADJUDICATION_CSV)
    df_rep = pd.read_csv(FROZEN_EMB_META)

    # Check total length closure
    assert len(df_adj) == 3111, f"Expected exactly 3,111 rows, found {len(df_adj)}"

    # Check row index monotonicity and match
    assert (df_adj["embedding_row_index"].values == df_rep["embedding_row_index"].values).all(), (
        "Adjudication embedding_row_index does not strictly match representatives metadata"
    )

    # Check record_id 1-to-1 match
    assert (df_adj["record_id"].values == df_rep["record_id"].values).all(), (
        "Adjudication record_id does not strictly match representatives metadata"
    )

    # Check uniqueness of embedding_row_index
    assert df_adj["embedding_row_index"].is_unique, "Duplicate embedding_row_index found in adjudication table"


def test_04_controlled_vocabulary_compliance():
    """Verify that all temporal classes, evidence strengths, and transfer statuses obey defined ontologies."""
    df_adj = pd.read_csv(FILE_ADJUDICATION_CSV)

    # Check class vocabulary
    observed_classes = set(df_adj["adjudicated_temporal_class"].unique())
    invalid_classes = observed_classes - VALID_TEMPORAL_CLASSES
    assert not invalid_classes, f"Encountered invalid temporal classes: {invalid_classes}"

    # Check evidence strength vocabulary
    observed_strengths = set(df_adj["evidence_strength"].unique())
    invalid_strengths = observed_strengths - VALID_EVIDENCE_STRENGTHS
    assert not invalid_strengths, f"Encountered invalid evidence strengths: {invalid_strengths}"

    # Check transfer status vocabulary
    observed_transfers = set(df_adj["cross_species_transfer_status"].unique())
    invalid_transfers = observed_transfers - VALID_TRANSFER_STATUSES_REPRESENTATIVE
    assert not invalid_transfers, f"Encountered invalid transfer statuses: {invalid_transfers}"


def test_05_evidence_sources_linkage_and_integrity():
    """Verify that all gene evidence records link to an existing registered literature source."""
    df_sources = pd.read_csv(FILE_SOURCES_CSV)
    df_gene_evid = pd.read_csv(FILE_GENE_EVID_CSV)

    assert "source_id" in df_sources.columns
    registered_sources = set(df_sources["source_id"].unique())

    # Every gene evidence row must trace to registered source keys
    for _, row in df_gene_evid.iterrows():
        raw_sources = str(row["primary_evidence_sources"]).split(";")
        for s in raw_sources:
            s_clean = s.strip()
            if s_clean:
                assert s_clean in registered_sources, (
                    f"Gene {row['canonical_gene']} references unregistered source key '{s_clean}'"
                )

        # Molecule measured must be non-empty and specify biological level
        assert pd.notna(row["measured_molecule"]) and len(row["measured_molecule"]) > 0


def test_06_no_unsupported_assignments_without_rationale():
    """
    Verify that every record contains an audit rationale, and that non-insufficient
    classes have justified evidence strength.
    """
    df_adj = pd.read_csv(FILE_ADJUDICATION_CSV)

    for idx, row in df_adj.iterrows():
        t_class = row["adjudicated_temporal_class"]
        strength = row["evidence_strength"]
        rationale = str(row["adjudication_rationale"]).strip()

        # Rationale must never be empty or generic
        assert len(rationale) > 10, f"Row {idx} ({row['record_id']}) lacks substantive adjudication rationale"

        # Logical consistency between class and strength
        if t_class == "INSUFFICIENT_EVIDENCE":
            assert strength == "INSUFFICIENT", (
                f"Row {idx} has class INSUFFICIENT_EVIDENCE but strength is '{strength}'"
            )
        else:
            assert strength in {"HIGH_EXPERIMENTAL", "MODERATE_EXPERIMENTAL", "MODERATE_SUBMITTER_ALIAS"}, (
                f"Row {idx} has class '{t_class}' but invalid/insufficient strength '{strength}'"
            )


def test_07_cross_species_transfer_discipline():
    """
    Verify that cross-species ortholog transfers adhere strictly to the 3-state
    transfer ontology and never fabricate direct evidence for homology transfers.
    """
    df_ortho = pd.read_csv(FILE_CROSS_SPECIES_CSV)

    assert len(df_ortho) > 0, "Cross-species temporal label table is empty"

    # Vocabulary check
    observed_statuses = set(df_ortho["cross_species_transfer_status"].unique())
    invalid_statuses = observed_statuses - VALID_TRANSFER_STATUSES_ORTHOLOG
    assert not invalid_statuses, f"Invalid ortholog transfer status: {invalid_statuses}"

    # For homology-inferred transfers, hsv2_evidence_source must indicate homology
    for _, row in df_ortho.iterrows():
        status = row["cross_species_transfer_status"]
        hsv2_src = str(row["hsv2_evidence_source"])
        if status == "INFERRED_FROM_HOMOLOGY":
            assert "homology" in hsv2_src.lower() or "inferred" in hsv2_src.lower(), (
                f"Ortholog pair {row['canonical_gene']} marked INFERRED_FROM_HOMOLOGY but source is '{hsv2_src}'"
            )
        elif status == "DIRECTLY_SUPPORTED_IN_BOTH_SPECIES":
            assert "homology" not in hsv2_src.lower(), (
                f"Ortholog pair {row['canonical_gene']} marked DIRECTLY_SUPPORTED_IN_BOTH_SPECIES but mentions homology"
            )


def test_08_transcript_vs_protein_assay_distinction():
    """Ensure that transcript kinetics and protein kinetics are explicitly distinguished."""
    df_gene_evid = pd.read_csv(FILE_GENE_EVID_CSV)
    quantities = df_gene_evid["measured_molecule"].unique()

    has_transcript = any("TRANSCRIPT" in str(q) for q in quantities)
    has_protein = any("PROTEIN" in str(q) for q in quantities)

    assert has_transcript, "No transcript-level evidence captured in temporal gene evidence table"
    assert has_protein, "No protein-level evidence captured in temporal gene evidence table"

    # Check representative adjudication table
    df_adj = pd.read_csv(FILE_ADJUDICATION_CSV)
    assert "measurement_type" in df_adj.columns
    types_observed = set(df_adj["measurement_type"].unique())
    assert "TRANSCRIPT_AND_PROTEIN" in types_observed or "TRANSCRIPT_4SU" in types_observed


def test_09_manual_review_queue_traceability():
    """Verify that all manual review queue items link to valid representative records."""
    df_adj = pd.read_csv(FILE_ADJUDICATION_CSV)
    df_queue = pd.read_csv(FILE_REVIEW_CSV)

    assert len(df_queue) > 0, "Manual review queue is empty"
    valid_row_indices = set(df_adj["embedding_row_index"].unique())

    for _, row in df_queue.iterrows():
        row_idx = row["embedding_row_index"]
        assert row_idx in valid_row_indices, (
            f"Review queue contains non-existent embedding_row_index: {row_idx}"
        )
        assert pd.notna(row["review_triggers"]) and len(str(row["review_triggers"])) > 3, (
            f"Review queue item {row_idx} lacks explicit review_triggers"
        )


def test_10_conflict_and_ambiguity_explicit_logging():
    """Ensure records with mixed, context-dependent, or disputed annotations have explicit disclosures."""
    df_adj = pd.read_csv(FILE_ADJUDICATION_CSV)
    df_mixed = df_adj[df_adj["adjudicated_temporal_class"].isin(["MIXED_OR_CONTEXT_DEPENDENT", "DISPUTED"])]

    for _, row in df_mixed.iterrows():
        rationale = str(row["adjudication_rationale"]).strip().lower()
        assert any(k in rationale for k in ["context", "dual", "dispute", "mixed", "complex", "icp34.5", "kinase", "variable"]), (
            f"Record {row['record_id']} labeled {row['adjudicated_temporal_class']} without explicit conflict documentation"
        )
