# Step 9.1 Governance Decision: Benchmark Readiness and Claim Verification

**Project**: *Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes*  
**Date**: 2026-10-21 14:36:02 UTC  
**Author**: Computational Virology Machine Learning Governance Committee  
**Deliverable Checksums**:
- `reports/09_1_count_reconciliation.md`: SHA-256 `50b337b4efa68a1d042efaa90eb0616f630bc1691cd3e230e412a9e55897cf27`
- `reports/09_1_record_level_differences.csv`: SHA-256 `9e45df4145afa31e768d3e91e2449b6ba547d1db89c79bd5076065c131cc835b`
- `reports/09_1_family_grouping_audit.md`: SHA-256 `5084a3ae138f7a84a0234c7854ee23480180bba3babbde595d72659c2b024f7a`
- `reports/09_1_grouping_exceptions.csv`: SHA-256 `bae72323dd7089598901ff2fcb51c1189fbe420c02b870e63b4f94af9b3ad466`
- `reports/09_1_fold_feasibility.csv`: SHA-256 `30728b18b4898e4fe58084ca1411471a4104b5b4f71b97faa4a16717dc92c2f9`
- `reports/09_1_split_validity_report.md`: SHA-256 `38472707fd3d52f972dbefd4a540e5e9b5a2d4b63499b714346227b3b3cb8691`

---

## 1. Scientific Claim Verification and Boundaries

In compliance with Phase D guidelines, the project establishes the following strict scientific boundaries:
1. **Unseen HSV Gene-Family vs. New Virus Generalization**:
   - The HSV genome encodes a fixed evolutionary repertoire of 63 canonical characterized genes. A group-aware fold tests whether a model can predict timing for a held-out HSV gene family; it **does NOT prove generalization across distantly related beta- or gammaherpesviruses**.
2. **Sequence Sample Size ($N=687$) vs. Family Sample Size ($N=63$)**:
   - The effective sample size for generalization is **63**, not 687. Sequence variants represent clinical diversity, not independent evolutionary experiments.
3. **Transcript vs. Directly Measured Protein Kinetics**:
   - 11 gene families have only nascent 4sU-seq RNA data; their protein accumulation kinetics are not directly measured and are segregated in sensitivity analyses.
4. **Direct Evidence vs. Homology-Inferred Labels**:
   - All primary benchmarks are restricted to the 687 direct-evidence records. Homology-inferred labels ($N=633$) remain strictly secondary.
5. **Prohibition of Premature Predictive Claims**:
   - We **explicitly refrain from claiming that ESM-2 encodes viral temporal kinetics** until empirical baseline and embedding models are evaluated on the leakage-resistant benchmark.

---

## 2. Official Readiness Decision

### **OFFICIAL DECISION: OPTION 1 — READY FOR BASELINE MODELING**

#### Justification:
1. **Mathematical Closure Verified**: All 3,111 representative sequences are accounted for without contradiction.
2. **Discrepancy Fully Resolved**: The discrepancy between Step 8.1 narrative text and record-level CSV data was investigated and documented down to individual accessions (e.g. `AAK97852.1`).
3. **Zero Group Leakage Proven**: The Stratified Group K-Fold design implemented in [09_group_split_manifest.csv](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/data/processed/09_group_split_manifest.csv) mathematically guarantees zero gene-family overlap across splits.
4. **Feasibility Confirmed**: All 5 folds contain adequate independent gene families for the binary replication timing task (`Pre-Replicative` vs. `Post-Replicative`).
5. **No Model Fitting Executed**: Zero classifiers were trained, and all historical artifacts remain intact.

The project is formally approved to proceed to Step 10: Baseline Modeling and Linear Evaluation.
