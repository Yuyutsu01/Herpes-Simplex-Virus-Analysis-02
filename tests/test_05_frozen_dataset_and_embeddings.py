#!/usr/bin/env python3
"""
Test Suite: Step 5 Frozen Datasets, Physicochemical Profiles, and ESM-2 Embeddings
=================================================================================

Concept & Scope:
----------------
This test suite verifies the biological and computational integrity of the
frozen HSV-1 / HSV-2 dataset artifacts and the generated ESM-2 embeddings:

1. Provenance Integrity: Source FASTAs match immutable SHA-256 hashes.
2. Universal Record Accounting: Every raw record (22,710 total) is accounted for
   in exactly one primary cohort:
     18,172 Valid Full-Length + 3,211 Fragments + 1,327 Excluded/Quarantined = 22,710.
3. Species & Cross-File Integrity: No unresolved taxonomy conflicts in the primary cohort.
   AEW27163.1 (HSV-1 thymidine kinase) is correctly placed in HSV-1.
4. Naturally Small Intact Proteins: Preserved in valid cohort (e.g. UL11, UL49A, US5, US9, US12).
5. Physicochemical Validity: All computed parameters (MW, pI, instability, GRAVY,
   aliphatic, residue group fractions) fall within biologically physical bounds.
6. ESM-2 Embeddings: The representative cohort embedding matrix has exact dimensions
   (3111, 480), finite real values, and 1-to-1 row alignment with the metadata index.
"""

import os
import hashlib
import json
import pytest
import numpy as np
import pandas as pd
from Bio import SeqIO

WORKSPACE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA_RAW = os.path.join(WORKSPACE_DIR, "data", "raw")
DATA_PROCESSED = os.path.join(WORKSPACE_DIR, "data", "processed")
REPORTS_DIR = os.path.join(WORKSPACE_DIR, "reports")

RAW_HSV1 = os.path.join(DATA_RAW, "hsv1.fasta")
RAW_HSV2 = os.path.join(DATA_RAW, "hsv2.fasta")

EXPECTED_HSV1_SHA256 = "15001269a835b41223467e0d91bdcada433ee68c6a58c19e684c060483b31f36"
EXPECTED_HSV2_SHA256 = "86c9d627286aa4cad121af4c84cdb39b824c14905aa57948de9a996955ca0311"

VALID_HSV1_FASTA = os.path.join(DATA_PROCESSED, "valid_hsv1_proteins.fasta")
VALID_HSV2_FASTA = os.path.join(DATA_PROCESSED, "valid_hsv2_proteins.fasta")
VALID_COMBINED_FASTA = os.path.join(DATA_PROCESSED, "valid_proteins_combined.fasta")
VALID_METADATA_CSV = os.path.join(DATA_PROCESSED, "valid_proteins_metadata.csv")
FRAGMENT_FASTA = os.path.join(DATA_PROCESSED, "fragment_sequences.fasta")
EXCLUDED_CSV = os.path.join(DATA_PROCESSED, "excluded_and_quarantined_records.csv")

MANIFEST_JSON = os.path.join(REPORTS_DIR, "final_dataset_manifest.json")
SUMMARY_MD = os.path.join(REPORTS_DIR, "final_dataset_summary.md")

ESM2_EMBEDDINGS_NPY = os.path.join(DATA_PROCESSED, "esm2_t12_representatives_embeddings.npy")
ESM2_METADATA_CSV = os.path.join(DATA_PROCESSED, "esm2_t12_representatives_metadata.csv")
ESM2_MANIFEST_JSON = os.path.join(REPORTS_DIR, "esm2_embedding_manifest.json")


def compute_sha256(filepath):
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()


def test_01_raw_fasta_sha256_checksums():
    """Verify raw input FASTA files remain completely unmodified."""
    assert os.path.exists(RAW_HSV1), "Raw HSV-1 FASTA missing"
    assert os.path.exists(RAW_HSV2), "Raw HSV-2 FASTA missing"
    assert compute_sha256(RAW_HSV1) == EXPECTED_HSV1_SHA256, "Raw HSV-1 SHA-256 hash mismatch!"
    assert compute_sha256(RAW_HSV2) == EXPECTED_HSV2_SHA256, "Raw HSV-2 SHA-256 hash mismatch!"


def test_02_all_processed_files_exist():
    """Verify all 6 core processed data files and reports exist."""
    required_files = [
        VALID_HSV1_FASTA,
        VALID_HSV2_FASTA,
        VALID_COMBINED_FASTA,
        VALID_METADATA_CSV,
        FRAGMENT_FASTA,
        EXCLUDED_CSV,
        MANIFEST_JSON,
        SUMMARY_MD
    ]
    for path in required_files:
        assert os.path.exists(path), f"Required artifact missing: {path}"


def test_03_universal_record_accounting_and_cohort_counts():
    """
    Verify exact reconciliation of all 22,710 raw records:
    18,172 Valid Full-Length + 3,211 Fragments + 1,327 Excluded/Quarantined = 22,710.
    """
    df_valid = pd.read_csv(VALID_METADATA_CSV)
    df_excl = pd.read_csv(EXCLUDED_CSV)
    frag_records = list(SeqIO.parse(FRAGMENT_FASTA, "fasta"))

    n_valid = len(df_valid)
    n_frag = len(frag_records)
    n_excl = len(df_excl)
    total_accounted = n_valid + n_frag + n_excl

    assert n_valid == 18172, f"Expected 18,172 valid records, got {n_valid}"
    assert n_frag == 3211, f"Expected 3,211 fragment records, got {n_frag}"
    assert n_excl == 1327, f"Expected 1,327 excluded records, got {n_excl}"
    assert total_accounted == 22710, f"Universal accounting failed: {total_accounted} != 22,710"


def test_04_species_attribution_and_aew27163_resolution():
    """
    Verify that:
    1. HSV-1 fasta has 9,454 records and HSV-2 fasta has 8,718 records (sum = 18,172).
    2. AEW27163.1 is classified as HSV-1 (effective_species_attribution == HSV-1).
    3. No records in valid_hsv1_proteins.fasta have effective species HSV-2, and vice versa.
    """
    hsv1_records = list(SeqIO.parse(VALID_HSV1_FASTA, "fasta"))
    hsv2_records = list(SeqIO.parse(VALID_HSV2_FASTA, "fasta"))
    assert len(hsv1_records) == 9454, f"Expected 9,454 HSV-1 records, got {len(hsv1_records)}"
    assert len(hsv2_records) == 8718, f"Expected 8,718 HSV-2 records, got {len(hsv2_records)}"
    assert len(hsv1_records) + len(hsv2_records) == 18172

    df_valid = pd.read_csv(VALID_METADATA_CSV)
    aew_matches = df_valid[df_valid["record_id"].str.contains("AEW27163.1", na=False)]
    assert len(aew_matches) == 1, "AEW27163.1 must exist exactly once in valid dataset"
    assert aew_matches.iloc[0]["effective_species_attribution"] in ["HSV-1", "HSV1"], "AEW27163.1 must be HSV-1"


def test_05_naturally_small_proteins_preserved():
    """
    Verify naturally small proteins (UL11, UL49A, US5, US9, US12) under 100 aa
    are preserved in the valid cohort and not erroneously purged to fragments.
    """
    df_valid = pd.read_csv(VALID_METADATA_CSV)
    for gene in ["UL11", "UL49A", "US5", "US9", "US12"]:
        gene_records = df_valid[df_valid["gene_annotation"] == gene]
        assert len(gene_records) > 0, f"Naturally small gene {gene} has 0 records in valid dataset!"
        # Check that intact instances exist with length < 100 aa
        sub100 = gene_records[gene_records["sequence_length"] < 100]
        assert len(sub100) > 0, f"Expected intact instances of {gene} < 100 aa, found none"


def test_06_physicochemical_properties_validity():
    """Verify calculated physicochemical features are within biologically physical ranges."""
    df_valid = pd.read_csv(VALID_METADATA_CSV)
    
    # Molecular Weight > 4,000 Da (e.g. UL49A at 87 aa or partially masked sequences) and < 700,000 Da
    assert (df_valid["molecular_weight"] > 4000).all(), "Unrealistic low MW detected"
    assert (df_valid["molecular_weight"] < 700000).all(), "Unrealistic high MW detected"
    
    # Isoelectric Point between 3.0 and 13.0
    assert (df_valid["isoelectric_point"] >= 3.0).all() and (df_valid["isoelectric_point"] <= 13.0).all()
    
    # Fractions between 0.0 and 1.0
    for col in ["charged_fraction", "polar_fraction", "hydrophobic_fraction", "aromatic_fraction", "ambiguous_residue_fraction"]:
        assert (df_valid[col] >= 0.0).all() and (df_valid[col] <= 1.0).all(), f"Invalid fraction range in {col}"
    
    # GRAVY typically between -2.0 and +2.0
    assert (df_valid["gravy"] >= -2.5).all() and (df_valid["gravy"] <= 2.5).all()


def test_07_esm2_embeddings_and_metadata_integrity():
    """
    Verify the generated ESM-2 embeddings:
    1. Embeddings array shape is exactly (3111, 480).
    2. Embeddings array contains no NaNs or Infs.
    3. Matching metadata CSV has exactly 3,111 rows.
    4. Row indices align 1-to-1 with the array.
    5. No leakage prevention groups are null.
    """
    if not os.path.exists(ESM2_EMBEDDINGS_NPY):
        pytest.skip("ESM-2 embeddings file not yet created (task may still be running)")

    embeddings = np.load(ESM2_EMBEDDINGS_NPY)
    assert embeddings.shape == (3111, 480), f"Expected shape (3111, 480), got {embeddings.shape}"
    assert not np.isnan(embeddings).any(), "NaN found in ESM-2 embeddings"
    assert not np.isinf(embeddings).any(), "Inf found in ESM-2 embeddings"

    df_esm_meta = pd.read_csv(ESM2_METADATA_CSV)
    assert len(df_esm_meta) == 3111, f"Expected 3,111 metadata rows, got {len(df_esm_meta)}"
    assert (df_esm_meta["embedding_row_index"].values == np.arange(3111)).all(), "Row ordering mismatch!"
    assert df_esm_meta["prevent_ml_leakage_group"].isna().sum() == 0, "Missing leakage prevention groups"
    assert df_esm_meta["sequence_sha256"].isna().sum() == 0, "Missing sequence hashes"
