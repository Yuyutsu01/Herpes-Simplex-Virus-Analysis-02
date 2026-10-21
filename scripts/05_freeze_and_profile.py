#!/usr/bin/env python3
"""
Step 5 Part 1: Dataset Freezing & Physicochemical Profiling
==========================================================
This script:
1. Ingests raw FASTAs and Step 4 validated redundancy decisions.
2. Generates the final frozen dataset files:
   - data/processed/valid_hsv1_proteins.fasta
   - data/processed/valid_hsv2_proteins.fasta
   - data/processed/valid_proteins_combined.fasta
   - data/processed/valid_proteins_metadata.csv
   - data/processed/fragment_sequences.fasta
   - data/processed/excluded_and_quarantined_records.csv
3. Computes comprehensive physicochemical features:
   - Length, molecular weight, isoelectric point (pI), instability index,
     aliphatic index, GRAVY, amino acid frequencies, residue group fractions,
     ambiguous residue fractions.
4. Generates reports/final_dataset_manifest.json and reports/final_dataset_summary.md.
"""

import os
import sys
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(line_buffering=True)
import re
import json
import hashlib
import datetime
from collections import Counter
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

WORKSPACE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA_RAW = os.path.join(WORKSPACE_DIR, "data", "raw")
FILE_HSV1 = os.path.join(DATA_RAW, "hsv1.fasta")
FILE_HSV2 = os.path.join(DATA_RAW, "hsv2.fasta")

DATA_PROCESSED = os.path.join(WORKSPACE_DIR, "data", "processed")
os.makedirs(DATA_PROCESSED, exist_ok=True)

REPORTS_DIR = os.path.join(WORKSPACE_DIR, "reports")
REDUNDANCY_CSV = os.path.join(REPORTS_DIR, "04_redundancy_decisions.csv")
HOMOLOGY_CSV = os.path.join(REPORTS_DIR, "04_homology_group_membership.csv")
CROSS_REL_CSV = os.path.join(REPORTS_DIR, "04_hsv1_hsv2_relationships.csv")

# Output dataset files
OUT_HSV1_FASTA = os.path.join(DATA_PROCESSED, "valid_hsv1_proteins.fasta")
OUT_HSV2_FASTA = os.path.join(DATA_PROCESSED, "valid_hsv2_proteins.fasta")
OUT_COMBINED_FASTA = os.path.join(DATA_PROCESSED, "valid_proteins_combined.fasta")
OUT_METADATA_CSV = os.path.join(DATA_PROCESSED, "valid_proteins_metadata.csv")
OUT_FRAGMENT_FASTA = os.path.join(DATA_PROCESSED, "fragment_sequences.fasta")
OUT_EXCLUDED_CSV = os.path.join(DATA_PROCESSED, "excluded_and_quarantined_records.csv")

OUT_MANIFEST_JSON = os.path.join(REPORTS_DIR, "final_dataset_manifest.json")
OUT_SUMMARY_MD = os.path.join(REPORTS_DIR, "final_dataset_summary.md")

HSV1_SHA256 = "15001269a835b41223467e0d91bdcada433ee68c6a58c19e684c060483b31f36"
HSV2_SHA256 = "86c9d627286aa4cad121af4c84cdb39b824c14905aa57948de9a996955ca0311"


def compute_sha256(filepath):
    """Computes SHA-256 hash of a file on disk."""
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()


def compute_physicochemical_properties(seq):
    """
    Computes rigorous biophysical features for a normalized amino acid sequence.
    Handles non-standard / ambiguous characters cleanly without throwing errors.
    """
    clean_seq = "".join([c for c in seq.upper() if c in "ACDEFGHIKLMNPQRSTVWY"])
    ambig_count = len(seq) - len(clean_seq)
    ambig_fraction = round(ambig_count / len(seq), 5) if len(seq) > 0 else 0.0

    if len(clean_seq) == 0:
        return {
            "molecular_weight": np.nan, "isoelectric_point": np.nan,
            "instability_index": np.nan, "aliphatic_index": np.nan,
            "gravy": np.nan, "charged_fraction": np.nan,
            "polar_fraction": np.nan, "hydrophobic_fraction": np.nan,
            "aromatic_fraction": np.nan, "ambiguous_residue_fraction": ambig_fraction
        }

    analysis = ProteinAnalysis(clean_seq)
    
    # 1. Molecular Weight (Daltons)
    try:
        mw = round(analysis.molecular_weight(), 2)
    except Exception:
        mw = np.nan

    # 2. Theoretical Isoelectric Point (pI)
    try:
        pi = round(analysis.isoelectric_point(), 3)
    except Exception:
        pi = np.nan

    # 3. Instability Index (Guruprasad et al. 1990: <40 stable, >40 unstable)
    try:
        instability = round(analysis.instability_index(), 2)
    except Exception:
        instability = np.nan

    # 4. GRAVY (Grand average of hydropathicity, Kyte & Doolittle 1982)
    try:
        gravy = round(analysis.gravy(), 3)
    except Exception:
        gravy = np.nan

    # 5. Amino Acid Frequencies
    aa_counts = analysis.count_amino_acids()
    l_clean = len(clean_seq)

    # 6. Aliphatic Index (Ikai 1980: X_Ala + 2.9*X_Val + 3.9*(X_Ile + X_Leu))
    x_ala = aa_counts.get('A', 0) / l_clean * 100.0
    x_val = aa_counts.get('V', 0) / l_clean * 100.0
    x_ile = aa_counts.get('I', 0) / l_clean * 100.0
    x_leu = aa_counts.get('L', 0) / l_clean * 100.0
    aliphatic_index = round(x_ala + 2.9 * x_val + 3.9 * (x_ile + x_leu), 2)

    # 7. Functional Residue Group Fractions
    charged = sum(aa_counts.get(a, 0) for a in ['D', 'E', 'R', 'K', 'H']) / l_clean
    polar = sum(aa_counts.get(a, 0) for a in ['S', 'T', 'N', 'Q', 'C', 'Y']) / l_clean
    hydrophobic = sum(aa_counts.get(a, 0) for a in ['A', 'V', 'I', 'L', 'M', 'F', 'W', 'P']) / l_clean
    aromatic = sum(aa_counts.get(a, 0) for a in ['F', 'W', 'Y']) / l_clean

    return {
        "molecular_weight": mw,
        "isoelectric_point": pi,
        "instability_index": instability,
        "aliphatic_index": aliphatic_index,
        "gravy": gravy,
        "charged_fraction": round(charged, 4),
        "polar_fraction": round(polar, 4),
        "hydrophobic_fraction": round(hydrophobic, 4),
        "aromatic_fraction": round(aromatic, 4),
        "ambiguous_residue_fraction": ambig_fraction
    }


def main():
    print("================================================================================")
    print("STEP 5: FREEZING DATASET & GENERATING PHYSICOCHEMICAL PROFILES")
    print("================================================================================\n")

    # 1. Verify Raw Input Integrity
    print("--- 1. VERIFYING RAW FASTA INTEGRITY ---")
    sha1 = compute_sha256(FILE_HSV1)
    sha2 = compute_sha256(FILE_HSV2)
    assert sha1 == HSV1_SHA256, f"Corrupted {FILE_HSV1}"
    assert sha2 == HSV2_SHA256, f"Corrupted {FILE_HSV2}"
    print(f"  [OK] hsv1.fasta SHA-256: {sha1}")
    print(f"  [OK] hsv2.fasta SHA-256: {sha2}")

    # 2. Ingest Validated Redundancy Decisions & Sequences
    print("\n--- 2. INGESTING STEP 4 REDUNDANCY DECISIONS & FASTA RECORDS ---")
    df_red = pd.read_csv(REDUNDANCY_CSV)
    print(f"  Loaded {len(df_red):,} total records from {REDUNDANCY_CSV}.")
    assert len(df_red) == 22710, "Record accounting mismatch!"

    # Read raw sequences directly from FASTAs
    raw_seq_dict = {}
    for fpath in [FILE_HSV1, FILE_HSV2]:
        fname = os.path.basename(fpath)
        for rec in SeqIO.parse(fpath, "fasta"):
            s_norm = str(rec.seq).strip().replace(" ", "").replace("\r", "").replace("\n", "").upper()
            raw_seq_dict[(fname, rec.id)] = (s_norm, rec.description)

    print(f"  Indexed {len(raw_seq_dict):,} sequence records from raw FASTAs.")

    # 3. Partition Cohorts
    print("\n--- 3. PARTITIONING COHORTS ---")
    # Dispositions:
    # 1. Primary Full-Length Cohort: RETAINED_UNIQUE_CANDIDATE (3,111) + RETAINED_VALID_VARIANT (15,061) = 18,172
    # 2. Fragment Cohort: RETAINED_FRAGMENT_SEPARATE_ANALYSIS = 3,211
    # 3. Excluded & Quarantined: EXCLUDED_CONFIRMED_NON_TARGET (1,043) + QUARANTINED_UNVERIFIED (254) + QUARANTINED_TAXONOMY_CONFLICT (30) = 1,327

    valid_mask = df_red["final_disposition"].isin(["RETAINED_UNIQUE_CANDIDATE", "RETAINED_VALID_VARIANT"])
    frag_mask = df_red["final_disposition"] == "RETAINED_FRAGMENT_SEPARATE_ANALYSIS"
    excl_mask = ~valid_mask & ~frag_mask

    df_valid = df_red[valid_mask].copy()
    df_frag = df_red[frag_mask].copy()
    df_excl = df_red[excl_mask].copy()

    print(f"  Primary Full-Length Target Cohort: {len(df_valid):,} records")
    print(f"    - Unique Representatives: {(df_valid['final_disposition'] == 'RETAINED_UNIQUE_CANDIDATE').sum():,}")
    print(f"    - Valid Strains / Variants: {(df_valid['final_disposition'] == 'RETAINED_VALID_VARIANT').sum():,}")
    print(f"    - Annotated Canonical Genes: {(df_valid['gene_annotation'] != 'UNANNOTATED').sum():,}")
    print(f"    - Provisional Unannotated: {(df_valid['gene_annotation'] == 'UNANNOTATED').sum():,}")
    print(f"  Fragment Cohort: {len(df_frag):,} records")
    print(f"  Excluded & Quarantined Cohort: {len(df_excl):,} records")
    assert len(df_valid) + len(df_frag) + len(df_excl) == 22710

    # 4. Compute Physicochemical Features for Valid Cohort
    print("\n--- 4. COMPUTING PHYSICOCHEMICAL PROPERTIES (18,172 RECORDS) ---")
    valid_records_data = []
    hsv1_fasta_records = []
    hsv2_fasta_records = []
    combined_fasta_records = []

    for idx, row in df_valid.iterrows():
        key = (row["source_file"], row["record_id"])
        seq_str, full_header = raw_seq_dict[key]
        
        # Calculate biophysical features
        props = compute_physicochemical_properties(seq_str)

        # Build clean header for FASTA outputs
        eff_sp = row["effective_species_attribution"]
        disp = row["final_disposition"]
        gene = row["gene_annotation"]
        cid = row["near_identity_cluster_id_95"]
        leak = row["prevent_ml_leakage_group"]
        
        fasta_header = f"{row['record_id']} species={eff_sp} gene={gene} cluster={cid} leakage_group={leak} disposition={disp} len={len(seq_str)}"
        fasta_entry = f">{fasta_header}\n{seq_str}\n"

        combined_fasta_records.append(fasta_entry)
        if eff_sp == "HSV1":
            hsv1_fasta_records.append(fasta_entry)
        else:
            hsv2_fasta_records.append(fasta_entry)

        entry_data = dict(row)
        entry_data["sequence"] = seq_str
        entry_data["original_header"] = full_header
        entry_data.update(props)
        valid_records_data.append(entry_data)

        if len(valid_records_data) % 3000 == 0:
            print(f"    Processed {len(valid_records_data):,} / {len(df_valid):,} sequences...", flush=True)

    df_valid_meta = pd.DataFrame(valid_records_data)

    # 5. Write Processed FASTAs and Metadata
    print("\n--- 5. WRITING FROZEN FASTA FILES & METADATA TABLE ---")
    with open(OUT_HSV1_FASTA, "w", encoding="utf-8") as f:
        f.writelines(hsv1_fasta_records)
    print(f"  Saved: {OUT_HSV1_FASTA} ({len(hsv1_fasta_records):,} records)")

    with open(OUT_HSV2_FASTA, "w", encoding="utf-8") as f:
        f.writelines(hsv2_fasta_records)
    print(f"  Saved: {OUT_HSV2_FASTA} ({len(hsv2_fasta_records):,} records)")

    with open(OUT_COMBINED_FASTA, "w", encoding="utf-8") as f:
        f.writelines(combined_fasta_records)
    print(f"  Saved: {OUT_COMBINED_FASTA} ({len(combined_fasta_records):,} records)")

    df_valid_meta.to_csv(OUT_METADATA_CSV, index=False)
    print(f"  Saved: {OUT_METADATA_CSV} ({len(df_valid_meta):,} rows, {len(df_valid_meta.columns)} features)")

    # 6. Write Fragment FASTA and Excluded Table
    print("\n--- 6. WRITING FRAGMENT FASTA & EXCLUDED/QUARANTINED TABLE ---")
    frag_fasta_records = []
    for idx, row in df_frag.iterrows():
        key = (row["source_file"], row["record_id"])
        seq_str, full_header = raw_seq_dict[key]
        header = f"{row['record_id']} species={row['effective_species_attribution']} gene={row['gene_annotation']} status={row['fragment_status']} len={len(seq_str)}"
        frag_fasta_records.append(f">{header}\n{seq_str}\n")

    with open(OUT_FRAGMENT_FASTA, "w", encoding="utf-8") as f:
        f.writelines(frag_fasta_records)
    print(f"  Saved: {OUT_FRAGMENT_FASTA} ({len(frag_fasta_records):,} records)")

    # Excluded & Quarantined records with original headers and reasons
    excl_rows = []
    for idx, row in df_excl.iterrows():
        key = (row["source_file"], row["record_id"])
        seq_str, full_header = raw_seq_dict[key]
        item = dict(row)
        item["original_header"] = full_header
        excl_rows.append(item)

    df_excl_meta = pd.DataFrame(excl_rows)
    df_excl_meta.to_csv(OUT_EXCLUDED_CSV, index=False)
    print(f"  Saved: {OUT_EXCLUDED_CSV} ({len(df_excl_meta):,} records)")

    # 7. Generate Comprehensive Manifest JSON
    print("\n--- 7. GENERATING FINAL DATASET MANIFEST ---")
    # Compute summary statistics
    length_dist = {
        "min": int(df_valid_meta["sequence_length"].min()),
        "max": int(df_valid_meta["sequence_length"].max()),
        "mean": round(float(df_valid_meta["sequence_length"].mean()), 2),
        "median": int(df_valid_meta["sequence_length"].median()),
        "std": round(float(df_valid_meta["sequence_length"].std()), 2)
    }

    manifest = {
        "manifest_version": "3.0.0_frozen",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "study_title": "Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes",
        "provenance_statement": (
            "Primary FASTA files (data/raw/hsv1.fasta, data/raw/hsv2.fasta) were downloaded from authoritative public "
            "repositories. Working analysis cohorts are established under the working assumption of database veracity, "
            "with full accession-by-accession manual curation deferred. Plausible unannotated viral proteins are preserved "
            "with explicit uncertainty flags."
        ),
        "raw_inputs": {
            "hsv1_fasta": {"path": FILE_HSV1, "records": 12170, "sha256": sha1},
            "hsv2_fasta": {"path": FILE_HSV2, "records": 10540, "sha256": sha2},
            "total_raw_records": 22710,
            "unique_raw_hashes": 22689
        },
        "frozen_cohort_counts": {
            "primary_valid_full_length_records": len(df_valid),
            "hsv1_primary_records": len(hsv1_fasta_records),
            "hsv2_primary_records": len(hsv2_fasta_records),
            "representative_cluster_seeds": int((df_valid["final_disposition"] == "RETAINED_UNIQUE_CANDIDATE").sum()),
            "strain_variant_records": int((df_valid["final_disposition"] == "RETAINED_VALID_VARIANT").sum()),
            "annotated_canonical_records": int((df_valid["gene_annotation"] != "UNANNOTATED").sum()),
            "provisional_unannotated_records": int((df_valid["gene_annotation"] == "UNANNOTATED").sum()),
            "fragment_cohort_records": len(df_frag),
            "excluded_non_target_records": int((df_red["final_disposition"] == "EXCLUDED_CONFIRMED_NON_TARGET").sum()),
            "quarantined_unverified_records": int((df_red["final_disposition"] == "QUARANTINED_UNVERIFIED").sum()),
            "quarantined_taxonomy_conflict_records": int((df_red["final_disposition"] == "QUARANTINED_TAXONOMY_CONFLICT").sum()),
            "total_accounted_records": 22710
        },
        "sequence_length_distribution": length_dist,
        "physicochemical_summary": {
            "mean_molecular_weight_da": round(float(df_valid_meta["molecular_weight"].mean()), 2),
            "mean_isoelectric_point": round(float(df_valid_meta["isoelectric_point"].mean()), 3),
            "mean_gravy": round(float(df_valid_meta["gravy"].mean()), 3),
            "mean_instability_index": round(float(df_valid_meta["instability_index"].mean()), 2)
        },
        "output_files": {
            "valid_hsv1_proteins_fasta": OUT_HSV1_FASTA,
            "valid_hsv2_proteins_fasta": OUT_HSV2_FASTA,
            "valid_proteins_combined_fasta": OUT_COMBINED_FASTA,
            "valid_proteins_metadata_csv": OUT_METADATA_CSV,
            "fragment_sequences_fasta": OUT_FRAGMENT_FASTA,
            "excluded_and_quarantined_records_csv": OUT_EXCLUDED_CSV,
            "manifest_json": OUT_MANIFEST_JSON,
            "summary_md": OUT_SUMMARY_MD
        }
    }

    with open(OUT_MANIFEST_JSON, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
    print(f"  Saved: {OUT_MANIFEST_JSON}")

    # 8. Generate Summary Markdown
    print("\n--- 8. GENERATING FINAL DATASET SUMMARY REPORT ---")
    summary_md = f"""# Final Frozen Analysis Dataset & Physicochemical Profiling Summary
**Study Title**: *Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes: Integrating Physicochemical Profiling, Protein Language Model Embeddings, and Unsupervised Learning*  
**Date**: {datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")}  
**Author**: Computational Virology Research Team  
**Manifest**: `reports/final_dataset_manifest.json`  
**Dataset Version**: 3.0.0 (Frozen)  

---

## 1. Executive Summary & Provenance Statement

The complete analysis dataset for Herpes Simplex Virus Type 1 (HSV-1) and Type 2 (HSV-2) has been formally frozen. All **22,710 original records** from the source FASTA archives have been assigned to mutually exclusive, biologically defensible cohorts.

### Provenance Assumption & Working Baseline:
Both raw FASTA archives were retrieved from trusted international databases (NCBI GenBank / RefSeq / UniProtKB). In accordance with the study protocol:
- Individual accession-level live database retrieval has been deferred.
- The working assumption of FASTA provenance is established, while acknowledging that 97.4% of eligible records depend on NCBI GenBank submitter descriptors rather than experimental Swiss-Prot records.
- Confirmed non-target retroviruses, host cellular interactors, and contradictory species conflicts have been purged or quarantined.
- Plausible unannotated viral open reading frames (5,997 records) have been preserved as provisional candidates rather than prematurely discarded.

---

## 2. Reconciled Cohort Architecture

| Analysis Cohort | Record Count | % of Corpus | File Artifact | Analytical Purpose |
| :--- | :---: | :---: | :--- | :--- |
| **Primary HSV-1 Full-Length** | **{len(hsv1_fasta_records):,}** | 44.59% | `data/processed/valid_hsv1_proteins.fasta` | Full-length HSV-1 candidates (includes AEW27163.1) |
| **Primary HSV-2 Full-Length** | **{len(hsv2_fasta_records):,}** | 35.43% | `data/processed/valid_hsv2_proteins.fasta` | Full-length HSV-2 candidates |
| **Combined Primary Cohort** | **{len(df_valid):,}** | **80.02%** | `data/processed/valid_proteins_combined.fasta` | Complete study-ready viral dataset for pLM & ML |
| **- Unique Representatives** | **{int((df_valid['final_disposition'] == 'RETAINED_UNIQUE_CANDIDATE').sum()):,}** | 13.70% | `metadata: RETAINED_UNIQUE_CANDIDATE` | 95% near-identity cluster seeds (non-redundant core) |
| **- Multi-Strain Variants** | **{int((df_valid['final_disposition'] == 'RETAINED_VALID_VARIANT').sum()):,}** | 66.32% | `metadata: RETAINED_VALID_VARIANT` | Clinical isolate strain variants ($\ge 95\%$ ID) |
| **- Annotated Canonical Genes** | **{int((df_valid['gene_annotation'] != 'UNANNOTATED').sum()):,}** | 53.64% | `metadata: gene != UNANNOTATED` | Core 74 canonical HSV gene families |
| **- Provisional Unannotated** | **{int((df_valid['gene_annotation'] == 'UNANNOTATED').sum()):,}** | 26.38% | `metadata: gene == UNANNOTATED` | Plausible viral candidates lacking canonical symbols |
| **Fragment Sensitivity Cohort** | **{len(df_frag):,}** | **14.14%** | `data/processed/fragment_sequences.fasta` | Truncated sequences and short unannotated peptides |
| **Excluded Non-Targets** | **{int((df_red['final_disposition'] == 'EXCLUDED_CONFIRMED_NON_TARGET').sum()):,}** | 4.59% | `data/processed/excluded_and_quarantined_records.csv` | Retroviruses (Gag/Pol/Env), host receptors, non-HSV |
| **Quarantined Unverified** | **{int((df_red['final_disposition'] == 'QUARANTINED_UNVERIFIED').sum()):,}** | 1.12% | `data/processed/excluded_and_quarantined_records.csv` | Synthetic patent claims, PDB crystallization chains |
| **Quarantined Species Conflict**| **{int((df_red['final_disposition'] == 'QUARANTINED_TAXONOMY_CONFLICT').sum()):,}** | 0.13% | `data/processed/excluded_and_quarantined_records.csv` | Contradictory cross-species exact duplicate pairs |
| **Total Accounted Records** | **22,710** | **100.00%** | **Universal Record Accounting** | **100.00% 1-to-1 traceability** |

---

## 3. Sequence Length & Physicochemical Profiles

Physicochemical properties were computed across all **18,172 valid full-length sequences**:

| Feature | Mean | Median | Min | Max | Standard Deviation |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **Sequence Length (aa)** | {length_dist['mean']:,} | {length_dist['median']:,} | {length_dist['min']:,} | {length_dist['max']:,} | {length_dist['std']:,} |
| **Molecular Weight (kDa)** | {round(df_valid_meta['molecular_weight'].mean() / 1000.0, 2)} | {round(df_valid_meta['molecular_weight'].median() / 1000.0, 2)} | {round(df_valid_meta['molecular_weight'].min() / 1000.0, 2)} | {round(df_valid_meta['molecular_weight'].max() / 1000.0, 2)} | {round(df_valid_meta['molecular_weight'].std() / 1000.0, 2)} |
| **Isoelectric Point (pI)** | {round(df_valid_meta['isoelectric_point'].mean(), 2)} | {round(df_valid_meta['isoelectric_point'].median(), 2)} | {round(df_valid_meta['isoelectric_point'].min(), 2)} | {round(df_valid_meta['isoelectric_point'].max(), 2)} | {round(df_valid_meta['isoelectric_point'].std(), 2)} |
| **Instability Index** | {round(df_valid_meta['instability_index'].mean(), 2)} | {round(df_valid_meta['instability_index'].median(), 2)} | {round(df_valid_meta['instability_index'].min(), 2)} | {round(df_valid_meta['instability_index'].max(), 2)} | {round(df_valid_meta['instability_index'].std(), 2)} |
| **Aliphatic Index** | {round(df_valid_meta['aliphatic_index'].mean(), 2)} | {round(df_valid_meta['aliphatic_index'].median(), 2)} | {round(df_valid_meta['aliphatic_index'].min(), 2)} | {round(df_valid_meta['aliphatic_index'].max(), 2)} | {round(df_valid_meta['aliphatic_index'].std(), 2)} |
| **GRAVY (Hydropathicity)** | {round(df_valid_meta['gravy'].mean(), 3)} | {round(df_valid_meta['gravy'].median(), 3)} | {round(df_valid_meta['gravy'].min(), 3)} | {round(df_valid_meta['gravy'].max(), 3)} | {round(df_valid_meta['gravy'].std(), 3)} |
| **Charged Fraction** | {round(df_valid_meta['charged_fraction'].mean() * 100.0, 2)}% | {round(df_valid_meta['charged_fraction'].median() * 100.0, 2)}% | {round(df_valid_meta['charged_fraction'].min() * 100.0, 2)}% | {round(df_valid_meta['charged_fraction'].max() * 100.0, 2)}% | — |
| **Polar Fraction** | {round(df_valid_meta['polar_fraction'].mean() * 100.0, 2)}% | {round(df_valid_meta['polar_fraction'].median() * 100.0, 2)}% | {round(df_valid_meta['polar_fraction'].min() * 100.0, 2)}% | {round(df_valid_meta['polar_fraction'].max() * 100.0, 2)}% | — |
| **Hydrophobic Fraction** | {round(df_valid_meta['hydrophobic_fraction'].mean() * 100.0, 2)}% | {round(df_valid_meta['hydrophobic_fraction'].median() * 100.0, 2)}% | {round(df_valid_meta['hydrophobic_fraction'].min() * 100.0, 2)}% | {round(df_valid_meta['hydrophobic_fraction'].max() * 100.0, 2)}% | — |

---

## 4. Machine Learning Leakage Prevention Guidelines

To maintain scientific validity in downstream machine learning tasks:
1. **Redundancy Reduction**: Primary models should be trained on the **3,111 unique cluster representatives** (`RETAINED_UNIQUE_CANDIDATE`). Evaluating models on all 18,172 records without cluster partitioning results in severe optimism bias due to near-identical strain variants.
2. **Partition Scheme**:
   - For variant-level generalization: use `near_identity_cluster_id_95`.
   - For gene-family holdout and kinetic class evaluation: use `prevent_ml_leakage_group` (holds out entire homologous families across HSV-1 and HSV-2).
"""

    with open(OUT_SUMMARY_MD, "w", encoding="utf-8") as f:
        f.write(summary_md)
    print(f"  Saved: {OUT_SUMMARY_MD}")

    print("\n================================================================================")
    print("STEP 5 PART 1 COMPLETE: DATASET FROZEN & PROFILED SUCCESSFULLY.")
    print("================================================================================")


if __name__ == "__main__":
    main()
