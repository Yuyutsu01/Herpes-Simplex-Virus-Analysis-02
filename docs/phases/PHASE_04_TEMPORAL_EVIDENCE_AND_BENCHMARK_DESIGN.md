# Phase 4: Temporal Evidence Curation & Benchmark Split Design

**Consolidated Steps**: Step 8 & 8.1 (Temporal Evidence Curation & Integrity Audit), Step 9 & 9.1 (Supervised Benchmark Design & Feasibility)  
**Primary Deliverables**: `data/processed/08_1_strict_temporal_benchmark.csv`, `data/processed/09_group_split_manifest.csv`  

---

## 1. Scientific Objective
Establish an authoritative, leakage-resistant benchmark for viral replication timing by curating primary biochemical literature, recovering submitter gene aliases, adjudicating experimental evidence tiers, and designing a group-disjoint cross-validation split that prevents homology leakage.

## 2. Primary Literature Sources & Evidence Tiers
- **Honess & Roizman (1974, 1975)**: Fundamental metabolic arrest assays using cycloheximide (CHX; protein synthesis inhibitor) and phosphonoacetic acid (PAA; viral DNA polymerase inhibitor).
- **Rutkowski et al. (2015)**: Time-resolved 4sU nascent RNA sequencing.
- **Kulej et al. (2017)**: SILAC quantitative mass spectrometry.
- **Evidence Tiers**:
  1. *Direct Protein-Validated* ($N=52$ families): Synthesis timing verified at the protein level under metabolic arrest.
  2. *Transcript-Only Annotated* ($N=11$ families): Kinetics inferred exclusively from 4sU RNA-seq.

## 3. Alias Recovery & Benchmark Partitioning
- A 10-gene regex dictionary resolved 270 submitter aliases (e.g., "alpha-trans-inducing factor" $	o$ `UL48`, "major capsid protein" $	o$ `UL19`).
- Partition of 3,111 representative cluster seeds:
  - **Strict Direct Evidence**: 687 sequences across 63 independent gene families.
  - **Homology Inferred**: 633 sequences (HSV-2 sequences assigned labels solely via HSV-1 orthology; excluded from direct benchmark).
  - **Unresolved / Accessory ORFs**: 1,791 sequences.

## 4. Group-Disjoint 5-Fold Cross-Validation Design
- All sequences belonging to the same gene family were allocated to exactly one validation fold ($	ext{Overlap} = 0$).
- **RS1 (ICP4) Family Isolation**: `RS1` contains 212 sequence variants (30.9% of the benchmark). It was allocated strictly to Fold 0 validation. Zero `RS1` sequences exist in Fold 0 training.

## 5. Validation & Test Outcomes
- Validated under `tests/test_08_temporal_annotation.py`, `test_08_1_temporal_integrity.py`, `test_09_benchmark_design.py`, `test_09_1_split_validity.py` (40 passing tests).
