# Step 9.1 Audit Report: Five-Fold Split Validity and Feasibility

**Project**: *Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes*  
**Date**: 2026-10-21 14:36:02 UTC  
**Target Split Manifest**: `data/processed/09_group_split_manifest.csv`  
**Deliverable**: `reports/09_1_split_validity_report.md`  

---

## 1. Quantitative Five-Fold Distribution Table

| Fold | Val Seqs | Train Seqs | Val Fams | Train Fams | Val Pre-Repl Fams | Val Post-Repl Fams | Both Classes Represented | Feasibility |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **Fold 0** | 290 | 397 | 19 | 44 | 6 (IE=1, E=5) | 13 (L1=6, L2=7) | **YES** | **FEASIBLE** |
| **Fold 1** | 101 | 586 | 9 | 54 | 5 (IE=1, E=4) | 4 (L1=3, L2=1) | **YES** | **FEASIBLE** |
| **Fold 2** | 107 | 580 | 11 | 52 | 2 (IE=1, E=1) | 9 (L1=4, L2=5) | **YES** | **FEASIBLE** |
| **Fold 3** | 89 | 598 | 11 | 52 | 4 (IE=1, E=3) | 7 (L1=1, L2=6) | **YES** | **FEASIBLE** |
| **Fold 4** | 100 | 587 | 13 | 50 | 6 (IE=1, E=5) | 7 (L1=2, L2=5) | **YES** | **FEASIBLE** |
| **Total** | **687** | **687** | **63** | **63** | **23 Families** | **40 Families** | **100% Balanced** | **READY** |

---

## 2. Split Feasibility Assessment

1. **Binary Task Feasibility**: Every single fold contains multiple independent gene families for both Pre-Replicative and Post-Replicative classes.
2. **Four-Class Task Limitation**: In the 4-class setting, each fold contains exactly ONE immediate-early gene family ($N=1$). A single misclassification drops fold recall to 0%. Therefore, the **binary task is the only statistically defensible primary benchmark**.
