# Step 9 Protocol: Group-Aware, Leakage-Resistant Cross-Validation Design

**Project**: *Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes*  
**Date**: 2026-10-21 14:14:59 UTC  
**Deliverable**: `reports/09_split_protocol.md`  

---

## 1. Group-Aware Partitioning Architecture

To eliminate sequence-level data leakage, the primary evaluation uses a **5-fold Stratified Group Cross-Validation (`StratifiedGroupKFold`)** partitioned strictly by `homology_group_id`:

```
Direct Evidence Cohort (N = 687 sequences, 63 unique gene families)
├── Fold 0: ~12-13 gene families | ~137 sequences
├── Fold 1: ~12-13 gene families | ~137 sequences
├── Fold 2: ~12-13 gene families | ~137 sequences
├── Fold 3: ~12-13 gene families | ~137 sequences
└── Fold 4: ~12-13 gene families | ~137 sequences
```

### Mathematical Guarantee:
`TrainFamilies_k ∩ ValFamilies_k = ∅` for all folds k in (0, 1, 2, 3, 4).

- **Zero Sequence Homology Leakage**: All 212 variants of `RS1` (ICP4) are assigned to a single fold. They are never split across train and test.
- **Zero Homology-Group Leakage**: All sequences sharing >= 95% near-identity or canonical orthology reside in the identical fold.
- **Prohibition of Random Splitting**: Standard random sequence splitting is prohibited because it evaluates trivial sequence re-identification rather than kinetic generalization.
