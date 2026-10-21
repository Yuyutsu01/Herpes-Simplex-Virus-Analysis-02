# Phase F Audit Report: Independent Ortholog-Pair Validation

**Audit Date**: 2026-10-21 13:16:18 UTC  
**Target Corpus**: 74 Canonical Core Alphaherpesvirus Ortholog Gene Pairs  
**Deliverable**: `reports/07_ortholog_pair_validation.csv`  

---

## 1. Executive Summary & Scientific Findings

Purported HSV-1 / HSV-2 ortholog pairs were reconstructed independently from the Step 4 canonical relationship table, without using cluster assignments to infer orthology.

### Reconstructed Ortholog Co-Clustering Rates:
1. **Universal Canonical Cohort (All 74 Genes)**:
   - Co-clustered under K-Means ($k=74$): **61 of 74 (82.43%)** [95% Wilson CI: **72.23% – 89.44%**].
2. **Strict Primary Cohort (Tiers 1 & 2 Curated/Probable, excluding complex nested loci)**:
   - Evaluated Genes: **69**
   - Co-clustered: **58 of 69 (84.06%)** [95% Wilson CI: **73.67% – 90.86%**].
3. **HDBSCAN Density Partitioning (Noise-Excluded Pairs)**:
   - Co-clustered: **46 of 47 (97.87%)** [95% Wilson CI: **88.89% – 99.62%**].

---

## 2. Analysis of Non-Co-Clustering Gene Pairs

Thirteen ortholog gene families split into different primary clusters between HSV-1 and HSV-2:
1. **Extreme Sequence Divergence**:
   - `US4` (Glycoprotein G): Only 34.6% sequence identity and 59.8% bidirectional alignment coverage. HSV-1 gG and HSV-2 gG possess completely different C-terminal secreted domains, which the ESM-2 model correctly captures as topologically distinct.
   - `US11` (RNA-binding tegument): 52.8% sequence identity with high C-terminal poly-proline repeat variation.
2. **Nested / Overlapping Loci**:
   - `UL26` (Protease) vs `UL26.5` (Capsid scaffold): Cross-contamination between protease precursors and cleaved scaffold forms in submitter records splits cluster centroids.
3. **Diploid Terminal Repeats**:
   - `RL1` (ICP34.5) and `RS1` (ICP4): Extensive GC-rich repeat expansions and sequencing gap masking (`X` runs) cause variant representations to disperse.
