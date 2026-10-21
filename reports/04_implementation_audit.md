# Step 4 Implementation Audit & Remediation Plan

**Project:** `Herpes-Simplex-Virus-Analysis-02`  
**Step:** 4 — Exact Duplicate Detection, Homology Clustering, and Protein Redundancy Analysis  
**Audit Date:** 2026-10-21  
**Status:** In Remediation (Audit Complete, Refactoring Active)

---

## 1. Executive Summary & Audit Scope

This audit provides a forensic evaluation of the initial implementation of Step 4 (`scripts/04_duplicate_homology_analysis.py`, its outputs in `reports/`, and the test suite in `tests/test_04_duplicate_homology.py`).

The initial implementation correctly established the raw sequence inventory (22,710 total records, 22,689 cryptographically distinct amino acid SHA-256 hashes, 0 intra-file exact duplicates, and 21 cross-file shared hashes). However, it suffered from severe mathematical, algorithmic, and biological failures in near-identity clustering, alignment calculations, orthology verification, non-target contamination filtering, and leakage-group assignment.

This document itemizes every identified failure, classifies its severity, details the root cause, and specifies the required remediation.

---

## 2. Itemized Defect Matrix

| ID | Issue Description | Severity | Affected Files & Lines | Root Cause | Required Correction |
| :--- | :--- | :---: | :--- | :--- | :--- |
| **DEF-01** | **Fabricated Near-Identity Clusters via String Concatenation** | **CRITICAL** | `scripts/04_duplicate_homology_analysis.py`: L375–377 | Cluster IDs were generated as `f"NEAR95_{hom_group_id}_{s_len // 50 * 50}"` rather than inferred from sequence similarity. No pairwise alignments or sequence comparisons were performed to form these clusters. | Implement genuine sequence-based greedy incremental clustering (CD-HIT algorithm) across eligible viral sequences using explicit alignment identity and coverage thresholds (95%, 90%, 70%). |
| **DEF-02** | **Invalid Near-Identity Sampling (702 pairs only)** | **CRITICAL** | `scripts/04_duplicate_homology_analysis.py`: L415–435 | Only the first 10 records of each annotation group were compared against record 0. 702 pairs cannot support corpus-wide clustering claims for 22,710 records. | Perform comprehensive clustering across all eligible viral sequences. Generate candidate pairs via k-mer indexing and length constraints, followed by exact alignment validation. |
| **DEF-03** | **Flawed Alignment Mathematics** | **CRITICAL** | `scripts/04_duplicate_homology_analysis.py`: L165–173 | Identity was calculated as `score / min_len * 100`, where `score` is an alignment score (including mismatch and gap penalties), not identity. Coverage was calculated as `min_len / len`, not aligned-residue coverage. Produced 0% identity for standard homologs. | Replace with exact aligned-residue identity: $N_{\text{identical}} / L_{\text{alignment}} \times 100\%$ (global) and $N_{\text{identical}} / \min(L_1, L_2) \times 100\%$ (local). Calculate true query and target coverage based on aligned columns. |
| **DEF-04** | **Unvetted Orthology Pair Generation & Erroneous 0% Identities** | **CRITICAL** | `scripts/04_duplicate_homology_analysis.py`: L470–515, `reports/04_hsv1_hsv2_relationships.csv` | Ortholog pairs paired the arbitrary first record matching a gene regex from HSV-1 and HSV-2. In many cases, fragments were aligned against full-length proteins under global alignment, yielding 0% identity (e.g. UL27/gB, US4/gG, US10, US12). | Base cross-species orthology on curated, authoritative reference genomes (HSV-1 strain 17 RefSeq NC_001806.2 and HSV-2 strain HG52 RefSeq NC_001798.2 / Swiss-Prot). Align full-length canonical proteins and validate bidirectional coverage. |
| **DEF-05** | **Contamination Leakage: 626 Retroviral Records Retained as Primary** | **CRITICAL** | `reports/04_redundancy_decisions.csv`, `scripts/04_duplicate_homology_analysis.py`: L540–580 | 626 records belonging to `HOM_NON_TARGET_RETROVIRAL` (Gag, Pol, Env) were assigned `RETAIN_PRIMARY_REPRESENTATIVE` because their NCBI headers erroneously stated `[Human alphaherpesvirus 1]`. | Re-audit all records for non-target retroviruses, host proteins, synthetic patents, and vectors. Strictly assign them to `EXCLUDED_CONFIRMED_NON_TARGET`. Zero non-target records permitted in primary viral cohorts. |
| **DEF-06** | **Unsupported Claims of "Biopython Alignment & K-mer Clustering"** | **HIGH** | `reports/04_step4_manifest.json`, `reports/04_step4_summary.md` | The manifest claimed k-mer clustering was applied, but no k-mer indexing, inverted index, or candidate search was implemented in the code. | Implement a genuine k-mer indexing candidate filter and document the exact parameters, data structures, and algorithmic complexity. |
| **DEF-07** | **Discrepancy Between Narrative Report and Generated Tables** | **HIGH** | `reports/04_step4_summary.md` vs `reports/04_hsv1_hsv2_relationships.csv` | Narrative text cited gB 87.1%, UL30 85.8%, UL5 88.4%, ICP0 61.3%, and ICP47 69.4%, while the actual CSV table recorded 0.00%, 82.19%, 79.06%, 17.78%, and 0.00%. | All narrative summaries must be programmatically derived from the generated output tables with 100% numerical consistency. |
| **DEF-08** | **Flawed ML Leakage Prevention Architecture** | **HIGH** | `scripts/04_duplicate_homology_analysis.py`: L545–560 | `prevent_ml_leakage_group` was set directly to the broad header regex. All 3,121 unannotated records were lumped into a single mega-group `HOM_UNANNOTATED`, making family-disjoint train/test splits useless. | Construct leakage groups hierarchically: gene-family level for verified homologs; sequence-cluster level (at 70% or 90% identity) for unannotated records. |
| **DEF-09** | **Fragment Accounting Discrepancies** | **MEDIUM** | `reports/04_step4_summary.md`, `reports/04_redundancy_decisions.csv` | Summary stated 1,745 fragments; decision table reported 1,636; raw FASTAs contain 1,790 sequences $<100$ aa. | Reconcile fragment definitions: distinguish "all short sequences $<100$ aa" (1,790) from "retained flagged viral fragments" (after excluding non-target short peptides and quarantined fragments). |
| **DEF-10** | **Superficial Test Suite** | **HIGH** | `tests/test_04_duplicate_homology.py` | Tests checked file existence, row counts, and weak non-empty conditions, failing to catch broken alignment math, fabricated cluster labels, and retained retroviral contaminants. | Build a comprehensive test suite including synthetic sequence benchmarks with hand-calculated expected values, non-target exclusion assertions, and algorithmic verification. |

---

## 3. Systematic Answers to the 10 Audit Questions (Section 1)

1. **Accountability of 22,710 raw records:**  
   Every raw record is present in the record table, but their classification was biologically flawed (e.g., retaining 626 retroviruses as primary representatives).
2. **Consistency of sequence hashes:**  
   Sequence hashes are computed consistently from normalized amino acid strings (uppercase, whitespace stripped). 22,689 unique hashes and 21 shared hashes are verified.
3. **Basis of near-identity clusters:**  
   Previously based on length bins and annotation text, NOT sequence similarity. This is being replaced with greedy incremental clustering based on pairwise alignments.
4. **Mathematical correctness of alignment identity and coverage:**  
   Previously incorrect (`score / min_len`). Being replaced with exact aligned-residue counting: $N_{\text{matches}} / L_{\text{alignment}}$.
5. **Transitive over-grouping / greedy order dependence:**  
   Greedy incremental clustering will sort sequences by length descending so full-length sequences act as cluster seeds, preventing fragments from acting as cluster centroids.
6. **Derivation of 78 homology groups:**  
   Previously assigned via regex matching on FASTA header strings, erroneously assigning retroviral Gag proteins to homology groups. Must be grounded in alignment to canonical reference proteins.
7. **Evidence for proposed ortholog pairs:**  
   Previously unvalidated, leading to fragment vs full-length pairings and 0% identity. Must be grounded in RefSeq / Swiss-Prot reference strain annotations.
8. **Handling of ambiguous/contradictory taxonomy:**  
   Previously, 8,429 records with provisional flags and 2,113 unannotated records were marked as primary representatives without qualifying their provisional nature.
9. **Classification of fragments, patents, PDB, and non-targets:**  
   PDB and patents were partially quarantined, but retroviruses and host proteins leaked into primary viral sets.
10. **Reconciliation of counts:**  
    Narrative text contradicted generated CSV tables. All counts must be strictly bound to computed output dataframes.

---

## 4. Remediation Implementation Plan

1. **Core Alignment Engine (`Bio.Align.PairwiseAligner`):**
   - Exact match counting, true query coverage, true target coverage.
   - Verified on synthetic test suite with known hand-calculated expected values.
2. **Eligible Target Cohort & Non-Target Purge:**
   - Immediately exclude 685 retroviral contaminants and 20 host cellular receptors.
   - Isolate 254 unverified records (PDB/patents) and 40 taxonomy conflicts.
3. **CD-HIT Style Greedy Incremental Clustering:**
   - Sort eligible candidates by length descending.
   - Candidate filtering via length difference ($\le 10\%$) and k-mer matching ($k=5$).
   - Exact pairwise alignment to assign sequences to clusters at 95%, 90%, and 70% identity.
4. **Authoritative Ortholog Alignment:**
   - Curate canonical HSV-1 (strain 17) and HSV-2 (strain HG52) reference proteins.
   - Perform global alignments and report true identity, coverage, and evidence tiers.
5. **Universal Mutually Exclusive Redundancy Decisions:**
   - Ensure sum across 7 mutually exclusive categories equals exactly 22,710.
6. **Comprehensive Automated Test Suite (`tests/test_04_duplicate_homology.py`):**
   - Test synthetic sequences, alignment math, non-target exclusion, and count consistency.
