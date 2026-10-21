# Step 9.1 Audit Report: Gene Family Grouping and Independence Audit

**Project**: *Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes*  
**Date**: 2026-10-21 14:36:02 UTC  
**Deliverable**: `reports/09_1_family_grouping_audit.md`  

---

## 1. Grouping Column Architecture

The project maintains three nested levels of sequence grouping:
1. `near_identity_cluster_id_95`: 95% near-identity clustering (CD-HIT/MMseqs2 equivalent).
2. `prevent_ml_leakage_group`: Clustered sequence groups ensuring near-identical proteins cannot cross splits.
3. `homology_group_id`: Canonical ortholog gene-family identifier (`HOM_RS1`, `HOM_UL19`, etc.).

---

## 2. Mega-Group Analysis: The `HOM_RS1` Asymmetry

- **`HOM_RS1` (ICP4)**: Contains **212 representative sequence records** in `DIRECT_EVIDENCE` (30.86% of the primary benchmark).
- **Impact on Splitting**:
  * Whichever cross-validation fold receives `HOM_RS1` will contain ~290 sequences, whereas folds without `HOM_RS1` contain ~90–107 sequences.
  * In [09_group_split_manifest.csv](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/data/processed/09_group_split_manifest.csv), all 212 records of `HOM_RS1` are strictly assigned to **Fold 0**.
  * **Zero Leakage**: No variant of `HOM_RS1` appears in Folds 1, 2, 3, or 4.

---

## 3. Mathematical Proof of Group Disjointness

Across all 5 folds in [09_group_split_manifest.csv](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/data/processed/09_group_split_manifest.csv):
`TrainFamilies_k ∩ ValFamilies_k = ∅` for all folds k in (0, 1, 2, 3, 4).
- Zero sequence records sharing canonical gene identity or >= 95% similarity cross any fold boundary.
