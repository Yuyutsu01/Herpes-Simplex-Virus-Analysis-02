# Step 9 Report: Class Feasibility and Statistical Generalizability

**Project**: *Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes*  
**Date**: 2026-10-21 14:14:58 UTC  
**Deliverable**: `reports/09_class_feasibility.md`  

---

## 1. Sequence-Level vs Gene-Family Discrepancy

Evaluating class feasibility at both sequence and gene-family levels reveals massive ascertainment distortion:

| Temporal Class | Sequences ($N = 687$) | % Seqs | Gene Families ($N = 63$) | % Families | Ratio Seqs/Family | Feasibility Status |
| :--- | :---: | :---: | :---: | :---: | :---: | :--- |
| **`IE_ALPHA`** | 251 | 36.54% | **5** | **7.94%** | **50.2 : 1** | **CRITICALLY UNDERPOWERED** |
| **`E_BETA`** | 99 | 14.41% | **18** | **28.57%** | 5.5 : 1 | **FEASIBLE WITH GROUPING** |
| **`L1_GAMMA1`** | 159 | 23.14% | **16** | **25.40%** | 9.9 : 1 | **FEASIBLE WITH GROUPING** |
| **`L2_GAMMA2`** | 178 | 25.91% | **24** | **38.10%** | 7.4 : 1 | **FEASIBLE WITH GROUPING** |

---

## 2. Statistical Consequences for Supervised Machine Learning

1. **The Immediate-Early ($N=5$) Generalization Paradox**:
   - `IE_ALPHA` appears to be the largest class at the sequence level (251 sequences, 36.5%), but it consists of **only 5 gene families in the entire virus**: `RS1` (ICP4), `RL2` (ICP0), `US1` (ICP22), `UL54` (ICP27), `US12` (ICP47).
   - In a 5-fold cross-validation held out at the gene family level, **each validation fold will contain exactly ONE immediate-early gene family**.
   - If the model misses a single gene family (e.g. ICP4), the validation recall for `IE_ALPHA` drops instantaneously from 100% to 0%.
   - **Conclusion**: A model cannot learn generalizable sequence features for "unseen immediate-early viral proteins" from a training set of only 4 genes.

2. **Feasibility of Binary / Multi-Class Tasks**:
   - **Full 4-Class Generalization Task (`IE` vs `E` vs `L1` vs `L2`)**: Statistically fragile due to N_families = 5 in `IE_ALPHA`.
   - **Binary Replication Timing Task (`Early / Pre-Replicative` vs `Late / Post-Replicative`)**:
     * Pre-Replicative (`IE_ALPHA` + `E_BETA`): 5 + 18 = **23 gene families** (350 sequences).
     * Post-Replicative (`L1_GAMMA1` + `L2_GAMMA2`): 16 + 24 = **40 gene families** (337 sequences).
     * This binary task is **substantially more balanced, statistically robust, and biologically meaningful** for evaluating whether ESM-2 embeddings encode promoter kinetics.
