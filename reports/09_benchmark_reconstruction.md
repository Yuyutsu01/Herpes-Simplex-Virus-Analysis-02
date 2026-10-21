# Step 9 Report: Benchmark Reconstruction and Dimensional Accounting

**Project**: *Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes*  
**Date**: 2026-10-21 14:14:58 UTC  
**Deliverable**: `reports/09_benchmark_reconstruction.md`  

---

## 1. Primary Benchmark Cohort Accounting

Reconstruction of `data/processed/08_1_strict_temporal_benchmark.csv` confirms exact 1-to-1 alignment with the frozen representative embedding array (3,111 rows, 480 dimensions):

```
Total Representatives: N = 3,111 sequences
├── 1. DIRECT_EVIDENCE:   687 records (22.08%) | 63 unique gene families | 63 homology groups
├── 2. HOMOLOGY_INFERRED: 633 records (20.35%) | 56 unique gene families | 56 homology groups
└── 3. UNRESOLVED:      1,791 records (57.57%) | Quarantined from supervised training
```

### Direct Evidence Primary Benchmark Breakdown:
Across the **687 sequence records** with direct literature validation:
- **`IE_ALPHA`**: 251 sequence records, but **ONLY 5 unique gene families** (`RS1`, `RL2`, `US1`, `UL54`, `US12`).
  * `RS1` alone accounts for **212 sequences (84.5% of all IE sequences)**!
- **`E_BETA`**: 99 sequence records across **18 unique gene families**.
- **`L1_GAMMA1`**: 159 sequence records across **16 unique gene families**.
- **`L2_GAMMA2`**: 178 sequence records across **24 unique gene families**.
- **Total Independent Biological Units**: Exactly **63 unique gene families** (homology groups).

> [!CRITICAL]
> **Pseudoreplication Discovery**: While the benchmark contains 687 representative sequences, it contains only **63 independent viral gene families**. Treating each sequence as an independent observation would artificially inflate the effective sample size by over 10-fold and mask severe class imbalance.

---

## 2. Species-Level Provenance & Balance

- **HSV-1 Representation**: **422 records** across 63 gene families.
- **HSV-2 Representation**: **265 records** strictly restricted to the 7 dual-supported gene families (`RS1`, `RL2`, `UL23`, `UL27`, `US4`, `US6`, `UL19`).
- **Homology-Inferred Cohort (HSV-2 only)**: **633 records** across 56 gene families where kinetics are transferred from HSV-1 homologs.
