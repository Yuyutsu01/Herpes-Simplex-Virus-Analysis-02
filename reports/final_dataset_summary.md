# Final Frozen Analysis Dataset & Physicochemical Profiling Summary
**Study Title**: *Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes: Integrating Physicochemical Profiling, Protein Language Model Embeddings, and Unsupervised Learning*  
**Date**: 2026-10-21 12:02:58 UTC  
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
| **Primary HSV-1 Full-Length** | **9,454** | 44.59% | `data/processed/valid_hsv1_proteins.fasta` | Full-length HSV-1 candidates (includes AEW27163.1) |
| **Primary HSV-2 Full-Length** | **8,718** | 35.43% | `data/processed/valid_hsv2_proteins.fasta` | Full-length HSV-2 candidates |
| **Combined Primary Cohort** | **18,172** | **80.02%** | `data/processed/valid_proteins_combined.fasta` | Complete study-ready viral dataset for pLM & ML |
| **- Unique Representatives** | **3,111** | 13.70% | `metadata: RETAINED_UNIQUE_CANDIDATE` | 95% near-identity cluster seeds (non-redundant core) |
| **- Multi-Strain Variants** | **15,061** | 66.32% | `metadata: RETAINED_VALID_VARIANT` | Clinical isolate strain variants ($\ge 95\%$ ID) |
| **- Annotated Canonical Genes** | **12,156** | 53.64% | `metadata: gene != UNANNOTATED` | Core 74 canonical HSV gene families |
| **- Provisional Unannotated** | **6,016** | 26.38% | `metadata: gene == UNANNOTATED` | Plausible viral candidates lacking canonical symbols |
| **Fragment Sensitivity Cohort** | **3,211** | **14.14%** | `data/processed/fragment_sequences.fasta` | Truncated sequences and short unannotated peptides |
| **Excluded Non-Targets** | **1,043** | 4.59% | `data/processed/excluded_and_quarantined_records.csv` | Retroviruses (Gag/Pol/Env), host receptors, non-HSV |
| **Quarantined Unverified** | **254** | 1.12% | `data/processed/excluded_and_quarantined_records.csv` | Synthetic patent claims, PDB crystallization chains |
| **Quarantined Species Conflict**| **30** | 0.13% | `data/processed/excluded_and_quarantined_records.csv` | Contradictory cross-species exact duplicate pairs |
| **Total Accounted Records** | **22,710** | **100.00%** | **Universal Record Accounting** | **100.00% 1-to-1 traceability** |

---

## 3. Sequence Length & Physicochemical Profiles

Physicochemical properties were computed across all **18,172 valid full-length sequences**:

| Feature | Mean | Median | Min | Max | Standard Deviation |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **Sequence Length (aa)** | 677.39 | 512 | 74 | 6,165 | 578.2 |
| **Molecular Weight (kDa)** | 71.34 | 55.04 | 4.49 | 645.73 | 60.65 |
| **Isoelectric Point (pI)** | 7.92 | 7.64 | 4.23 | 12.0 | 1.98 |
| **Instability Index** | 52.88 | 52.05 | 17.3 | 135.65 | 12.61 |
| **Aliphatic Index** | 78.45 | 81.39 | 20.66 | 125.41 | 14.74 |
| **GRAVY (Hydropathicity)** | -0.241 | -0.214 | -2.054 | 0.96 | 0.315 |
| **Charged Fraction** | 22.94% | 22.96% | 10.0% | 52.5% | — |
| **Polar Fraction** | 21.46% | 21.6% | 3.85% | 41.38% | — |
| **Hydrophobic Fraction** | 47.19% | 47.12% | 20.67% | 62.26% | — |

---

## 4. Machine Learning Leakage Prevention Guidelines

To maintain scientific validity in downstream machine learning tasks:
1. **Redundancy Reduction**: Primary models should be trained on the **3,111 unique cluster representatives** (`RETAINED_UNIQUE_CANDIDATE`). Evaluating models on all 18,172 records without cluster partitioning results in severe optimism bias due to near-identical strain variants.
2. **Partition Scheme**:
   - For variant-level generalization: use `near_identity_cluster_id_95`.
   - For gene-family holdout and kinetic class evaluation: use `prevent_ml_leakage_group` (holds out entire homologous families across HSV-1 and HSV-2).
