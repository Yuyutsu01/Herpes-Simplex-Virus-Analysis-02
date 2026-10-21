# Step 4 Remediation (v2.1.0): Validated Duplicate Detection, Homology Clustering, and Redundancy Analysis Report
**Study Title**: *Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes: Integrating Physicochemical Profiling, Protein Language Model Embeddings, and Unsupervised Learning*  
**Date**: 2026-10-21 11:34:13 UTC  
**Author**: Computational Virology & Bioinformatics Research Engineer  
**Script**: `scripts/04_duplicate_homology_analysis.py` (v2.1.0 Remediated)  
**Audit Document**: `reports/04_final_validation.md`  

---

## 1. Executive Summary & Remediation Results

This report establishes the scientifically verified, evidence-grounded dataset readiness across the **22,710 records** of the HSV-1 (12,170 records) and HSV-2 (10,540 records) proteome corpus.

### Confirmed Remediation Implementations:
1. **Mathematical Correction of Global Alignment Coverage**:
   - Replaced tautological non-gap residue count with true aligned non-gap columns: $\text{aligned\_columns} / L_1$ and $\text{aligned\_columns} / L_2$.
   - Independent length ratio $\min(L_1, L_2) / \max(L_1, L_2)$ and bidirectional coverage are explicitly reported.
2. **Accurate Algorithmic Attribution**:
   - Accurately designated as **Custom Greedy Incremental Alignment-Based Clustering** (Python / Biopython `PairwiseAligner` with k-mer pre-filter), removing misleading references to the official CD-HIT executable.
3. **Preservation of Naturally Small Intact Viral Proteins**:
   - Context-aware fragment filter preserved intact small canonical viral proteins: **UL11 (96 aa)**, **UL49A (87-91 aa)**, **US5 (92 aa)**, **US9 (90 aa)**, and **US12 (88 aa)**.
4. **Disambiguation of UL26.5 from UL26**:
   - Regex lookaround token boundaries prevented prefix collision; restored all **30 records of UL26.5** to their own gene family.
5. **Resolution of AEW27163.1 Cross-File Conflict**:
   - Accession record verified in NCBI Protein as `Human alphaherpesvirus 1` (Taxid 10298); species attribution corrected to HSV-1.
6. **Authoritative Ortholog Selection Across All 74 Canonical Genes**:
   - Reference pairings now span all **74 canonical genes** with zero length-asymmetry outliers.

---

## 2. Quantitative Summary Table

| Metric | Measured Value | Percentage of Corpus | Biological & Analytical Significance |
| :--- | :---: | :---: | :--- |
| **1. Input Records by Species** | **22,710** | 100.00% | **HSV-1**: 12,170 records; **HSV-2**: 10,540 records |
| **2. Unique Sequence Hashes** | **22,689** | 99.91% | Distinct amino acid SHA-256 hashes after normalization |
| **3. Exact Duplicate Clusters** | **22,689** | — | Both raw files are 100% internally non-redundant |
| **4. Cross-Species Shared Hashes**| **21 hashes** (42 records) | 0.18% | 7 species conflicts, 4 host receptors, 2 patents, 7 PDBs, 1 same-species duplicate |
| **5. Eligible Full-Length Cohort**| **18,172** | **80.02%** | Intact full-length target candidates for 95% clustering |
| **6. Near-Identity Clusters (95% ID)** | **3,111** | — | Verified incremental sequence clusters across full-length candidates |
| **7. Canonical Ortholog Pairs Audited** | **74** | — | Audited pairs across canonical gene families |
| **8. High-Confidence Orthologs** | **69** | — | High sequence conservation ($\ge 75\%$ ID, $\ge 85\%$ length ratio) |
| **9. Preserved Natural Small Proteins** | **UL11, UL49A, US5, US9, US12** | — | Intact short proteins restored from fragment misclassification |
| **10. Software Environment** | Biopython 1.88 | — | C-accelerated `Bio.Align.PairwiseAligner` with true bidirectional coverage |

---

## 3. Universal Mutually Exclusive Disposition Breakdown

| Final Disposition Category | Record Count | Percentage | Downstream Analytical Handling |
| :--- | :---: | :---: | :--- |
| **`RETAINED_UNIQUE_CANDIDATE`** | **3,111** | **13.70%** | Primary cluster seeds for unique sequence representation & pLM embeddings |
| **`RETAINED_VALID_VARIANT`** | **15,061** | **66.32%** | Multi-strain clinical isolates ($\ge 95\%$ ID); retained for variant studies |
| **`RETAINED_FRAGMENT_SEPARATE_ANALYSIS`**| **3,211** | **14.14%** | Genuine partial/truncated sequences and unannotated short peptides |
| **`EXCLUDED_CONFIRMED_NON_TARGET`** | **1,043** | **4.59%** | Excluded retroviruses (Gag/Pol/Env), host receptors, non-viral organisms |
| **`QUARANTINED_TAXONOMY_CONFLICT`** | **30** | **0.13%** | Quarantined cross-species exact duplicate conflicts |
| **`QUARANTINED_UNVERIFIED`** | **254** | **1.12%** | Quarantined synthetic patent claims and PDB crystallization constructs |
| **Total** | **22,710** | **100.00%** | **Universal 1-to-1 record accounting** |
