# Phase 1: Data Inventory, Taxonomic Curation & Redundancy Reduction

**Consolidated Steps**: Step 1 (FASTA Inventory), Step 2 (Taxonomy Audit), Step 3 (Target Proteome Curation), Step 4 (Homology & Redundancy Reduction)  
**Primary Deliverables**: `reports/04_redundancy_decisions.csv`, `reports/04_final_validation.md`, `data/processed/curated_protein_candidates.fasta`  

---

## 1. Scientific Objective
Establish a clean, non-redundant, authoritative viral protein sequence corpus from raw public databases (NCBI GenBank, RefSeq, UniProtKB) for Herpes Simplex Virus Type 1 (HSV-1) and Type 2 (HSV-2). The primary challenge in viral computational proteomics is that public sequence repositories contain massive clinical strain redundancy, truncated sequencing fragments, host interactors, and recombinant vector contaminants (such as HIV-1 Gag/Pol).

## 2. Inputs & Source Data
- `data/raw/hsv1.fasta` (SHA-256: `15001269a835b41223467e0d91bdcada433ee68c6a58c19e684c060483b31f36`, Size: 7,954,295 bytes, 12,028 records)
- `data/raw/hsv2.fasta` (SHA-256: `86c9d627286aa4cad121af4c84cdb39b824c14905aa57948de9a996955ca0311`, Size: 7,095,643 bytes, 10,682 records)
- Universal total: **22,710 records**

## 3. Methodology & Transformations
1. **Initial Sequence Parsing & Exact Hashing**: Each FASTA entry was digested with SHA-256 to identify identical amino acid sequences across accessions.
2. **Taxonomic & Host Contaminant Audit**:
   - Identified 1,043 non-target records including HIV-1 Gag/Pol/Env structural proteins from lentiviral gene therapy vectors and host ribosomal interactors.
   - Quarantined 254 unverified entries (synthetic patent claims, PDB constructs).
   - Quarantined 30 records (15 duplicate pairs) exhibiting identical sequences with contradictory species labels (`[Human alphaherpesvirus 1]` vs `[Human alphaherpesvirus 2]`).
3. **Length Thresholding**: Sequestered 3,211 truncated fragments (< 70 amino acids) into a separate fragment cohort, leaving 18,172 valid full-length candidates (HSV-1: 9,454; HSV-2: 8,718).
4. **95% Identity Redundancy Reduction (CD-HIT)**:
   - Partitioned the 18,172 candidates into **3,111 representative cluster seeds** (13.70% of corpus) and **15,061 near-identical strain variants** (66.32% of corpus).
   - Mapped all 3,111 representative seeds to canonical gene families (74 ortholog pairs).

## 4. Key Outputs & Provenance
- `reports/04_redundancy_decisions.csv`: Exhaustive 22,710-row ledger documenting the disposition and cluster assignment of every sequence.
- `data/processed/valid_proteins_combined.fasta`: 18,172 valid full-length viral candidates.
- `data/processed/esm2_t12_representatives_metadata.csv`: Metadata for the 3,111 representative cluster seeds.

## 5. Validation & Test Outcomes
- Validated under `tests/test_04_duplicate_homology.py` (10 passing tests).
- Mathematical closure confirmed: 18,172 (valid) + 3,211 (fragments) + 1,043 (purged) + 254 (unverified) + 30 (conflicts) = exactly 22,710 records.
