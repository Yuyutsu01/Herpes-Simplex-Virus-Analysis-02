# Phase B Gene Label Audit: Annotation Provenance and Evidence Tiers

**Audit Date**: 2026-10-21 13:08:49 UTC  
**Target Corpus**: Frozen Representative Cohort (N = 3,111 sequences)  
**Deliverable**: `reports/07_gene_label_audit.csv`  

---

## 1. Executive Summary & Provenance Disclosure

In accordance with Step 7 requirements, we audited the provenance of the gene labels used in all prior analyses. Because live individual accession-by-accession web lookups were deferred:
1. **Submitter Header Dependency**: 96.8% of annotated representative records rely on GenBank submitter descriptors rather than curated Swiss-Prot experimental entries.
2. **Four Distinct Evidence Tiers**:
   - **Tier 1 (Curated RefSeq / Swiss-Prot)**: **5 records (0.2%)**. Authoritative reference standards.
   - **Tier 2 (Probable Submitter Header)**: **1,211 records (38.9%)**. Unambiguous standard gene tokens in GenBank headers.
   - **Tier 3 (Ambiguous / Complex Loci)**: **608 records (19.5%)**. Overlapping open reading frames (`UL26` vs `UL26.5`), diploid terminal repeats (`RL1`, `RL2`, `RS1`), and unannotated records containing descriptive product text (e.g. "alpha trans-inducing factor" for VP16/UL48).
   - **Tier 4 (Provisional Unannotated)**: **1,287 records (41.4%)**. Truly uncharacterized open reading frames or hypothetical proteins.

---

## 2. Evidence Tier Distribution Table

| Evidence Tier | Record Count | % Cohort | Confidence Level | Description |
| :--- | :---: | :---: | :--- | :--- |
| **Tier 1: Curated Reference** | **5** | 0.2% | High | NCBI RefSeq (NP_/YP_) or Swiss-Prot entries |
| **Tier 2: Probable Header** | **1,211** | 38.9% | Moderate | Explicit submitter gene symbols in standard GenBank records |
| **Tier 3: Ambiguous / Complex** | **608** | 19.5% | Moderate-Low | Nested genes (UL26/UL26.5), repeats (RL1/RS1), or missed product text |
| **Tier 4: Provisional Unannotated**| **1,287** | 41.4% | Low | Hypothetical proteins lacking canonical gene symbols |
| **Total Evaluated** | **3,111** | **100.0%** | Universal | Complete audit closure |

---

## 3. Critical Findings on Gene Labeling

1. **Nested Reading Frames**: `UL26` (capsid protease) and `UL26.5` (capsid scaffold) share an identical C-terminal sequence. In clinical sequencing without transcript mapping, submitters frequently label full-length scaffold precursors as `UL26` or vice versa.
2. **Hidden Product Aliases**: 199 records currently classified as `UNANNOTATED` contain explicit functional product descriptions (e.g. `alpha trans-inducing factor` for `UL48`, `major capsid protein` for `UL19`). While these do not invalidate Step 6 clustering, they indicate that true unannotated novel proteins are fewer than 1,489.
3. **Absence of Independent Verification**: We explicitly reiterate that individual accession-level database lookups remain deferred; gene labels should be treated as provisional evidence tiers rather than absolute ground truth.
