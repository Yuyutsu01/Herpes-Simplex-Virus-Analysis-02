# Step 9 Final Report: Benchmark Readiness Decision and Task Feasibility Assessment

**Project**: *Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes*  
**Date**: 2026-10-21 14:14:59 UTC  
**Author**: Computational Virology Machine Learning Governance Committee  
**Target Benchmark**: Direct Evidence Primary Cohort (N = 687 sequences, 63 independent gene families)  
**Deliverables Generated**:
- `reports/09_benchmark_reconstruction.md`: SHA-256 `458e3c48f3965be89e521c7bdb82a8ccb898f55378c4a1f692b2880f1353934f`
- `reports/09_benchmark_counts.csv`: SHA-256 `369d01e044b865ca5011ca6a91dfc1cdff765ec3f43cbaba473eb0d3b3df20b4`
- `reports/09_gene_family_independence_audit.csv`: SHA-256 `aa6bf4e1000059cbf7760c166fcf5ef8b0c21af1fa1ffba3d9a659c069031377`
- `reports/09_label_provenance_audit.csv`: SHA-256 `c8b8308cb70f881b6a687c960aea09723088526f6d3572569b3bc8f6c66bf725`
- `reports/09_class_feasibility.md`: SHA-256 `0cfcbf6ab74b3fdbb96e389c7d9ce1e25d8e27c8e0bf97389c17fe650b346a07`
- `reports/09_class_feasibility.csv`: SHA-256 `362171c3983e2df3714e5d3059b182672d8f732a2e67407f0da439e9cb5d8ffd`
- `reports/09_split_protocol.md`: SHA-256 `bb71fc209fffa8a4ff5331cfcb7c862d14df0abca1f3dcdc4f3688b93c4efdaa`
- `data/processed/09_group_split_manifest.csv`: SHA-256 `44c46b4bdd13022b04ea6e6352119b0cbb7923eb259c01e35dd29f0a0aa38065`
- `reports/09_baseline_protocol.md`: SHA-256 `b5f756cde93e548e82a79ee4275ef03ec2e57696ef99ffd9ca1d0947fd4b7a0a`
- `reports/09_evaluation_metrics.md`: SHA-256 `9cd12f1908e8fb0fa435e486bbefe8b01aebb23c2e5e6be6785c232ddb9c18bb`
- `reports/09_sensitivity_protocol.md`: SHA-256 `de0ed4a465a20416ad79990168cf60817a534ba5f5e120db3eeae9f5e89339db`

---

## 1. Executive Summary & Mandatory Scientific Rules Maintained

Step 9 evaluated whether the audited HSV-1/HSV-2 temporal kinetics annotations support a scientifically defensible supervised prediction task:
- **No Classifiers Trained**: In strict adherence to protocol, **zero classifiers were trained or tuned**, and no predictive performance claims were manufactured.
- **Exact Reconstruction**: Verified exact closure of the 3,111 representative sequences (687 direct + 633 homology + 1,791 unresolved).
- **Leakage Elimination**: Formulated an audited, group-aware cross-validation manifest ([09_group_split_manifest.csv](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/data/processed/09_group_split_manifest.csv)) guaranteeing zero gene-family overlap across splits.

---

## 2. Answers to the Seven Governance Questions

### Question 1: How many independent gene families are available per class?
- **`IE_ALPHA`**: Exactly **5 gene families** (`RS1`, `RL2`, `US1`, `UL54`, `US12`).
- **`E_BETA`**: Exactly **18 gene families**.
- **`L1_GAMMA1`**: Exactly **16 gene families**.
- **`L2_GAMMA2`**: Exactly **24 gene families**.
- **Total Independent Units**: Exactly **63 gene families**.

### Question 2: Can a held-out test set contain genuinely unseen gene families?
**Yes**, but with severe statistical power constraints for `IE_ALPHA`. In a 5-fold cross-validation partitioned by gene family, each fold contains exactly ONE immediate-early gene family. A single classification error produces a 0% recall for that class in that fold.

### Question 3: Can related sequences be kept in the same split?
**Yes**. The Stratified Group K-Fold design implemented in [09_group_split_manifest.csv](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/data/processed/09_group_split_manifest.csv) groups sequences strictly by `homology_group_id`. All 212 variants of `RS1` are isolated into a single fold, mathematically precluding sequence-level leakage.

### Question 4: Are labels sufficiently independent and traceable?
**Yes, but at the gene-family level**. Every sequence traces back to primary literature (Honess & Roizman 1974/1975, Rutkowski 2015, Kulej 2017). However, multiple sequences within a gene family share the identical experimental literature source.

### Question 5: Is the task prediction of a known gene's temporal class, or generalization to unseen gene families?
This is the **critical scientific distinction**:
- **Known Gene / Variant Annotation**: A trivial family lookup table achieves near-100% accuracy because clinical variants of known genes retain their temporal class. This does not require complex deep learning.
- **De Novo Kinetic Generalization**: Predicting temporal class for novel viral gene families purely from sequence/ESM-2 embeddings is severely constrained by the small sample size (N_families = 63).

### Question 6: What performance estimates can the available sample size realistically support?
With only 63 independent families:
- Macro-F1 estimates will exhibit wide 95% confidence intervals (+/- 12-18%).
- Generalization to unseen immediate-early proteins cannot be validated with statistical confidence from N=5 training instances.

### Question 7: Should supervised learning proceed, be narrowed, or be postponed?

#### **Official Governance Decision: PROCEED WITH A NARROWED, SCIENTIFICALLY DEFENSIBLE TASK**

1. **Reject Naive 4-Class Generalization as Primary**: We explicitly reject claiming that an ESM-2 classifier can reliably learn a generalizable 4-class predictor for unseen viral immediate-early proteins from 5 examples.
2. **Approve Binary Replication-Timing Benchmark**:
   - Supervised classification is approved for a **Binary Replication Timing Task**:
     * **Pre-Replicative (IE + E)**: 23 independent gene families.
     * **Post-Replicative (L1 + L2)**: 40 independent gene families.
   - This task tests the genuine biological hypothesis of whether protein language model embeddings encode structural/functional constraints associated with pre-replicative replication enzymes vs post-replicative structural virion components.
3. **Mandatory Group-Aware Evaluation**: Any future model training must execute strictly against [09_group_split_manifest.csv](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/data/processed/09_group_split_manifest.csv). Sequence-random splitting is permanently prohibited.
