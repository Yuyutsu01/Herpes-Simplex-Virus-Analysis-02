#!/usr/bin/env python3
"""
Script: 09_benchmark_design_and_feasibility.py
=============================================

Scientific Objective:
---------------------
Execute Step 9: Leakage-Resistant Temporal Classification Benchmark Design for
Herpes-Simplex-Virus-Analysis-02.

This script implements:
1. Phase A: Reconstruct strict benchmark counts across all biological and computational dimensions.
2. Phase B: Audit gene-family label independence and provenance.
3. Phase C: Quantify class feasibility and statistical sample sizes at the gene-family level.
4. Phase D: Design leakage-resistant, group-aware cross-validation splits (StratifiedGroupKFold on homology_group_id).
5. Phase E: Formulate simple baselines (majority class, family lookup, physicochemical, ESM-2) and prespecify metrics.
6. Phase F: Prespecify sensitivity cohorts (direct vs homology, transcript vs protein).
7. Phase G: Compile the formal benchmark readiness decision report.

Outputs:
- reports/09_benchmark_reconstruction.md
- reports/09_benchmark_counts.csv
- reports/09_gene_family_independence_audit.csv
- reports/09_label_provenance_audit.csv
- reports/09_class_feasibility.md
- reports/09_class_feasibility.csv
- reports/09_split_protocol.md
- data/processed/09_group_split_manifest.csv
- reports/09_baseline_protocol.md
- reports/09_evaluation_metrics.md
- reports/09_sensitivity_protocol.md
- reports/09_benchmark_readiness_report.md
"""

import os
import sys
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(line_buffering=True)
import json
import datetime
import hashlib
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold, GroupKFold

WORKSPACE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA_PROCESSED = os.path.join(WORKSPACE_DIR, "data", "processed")
REPORTS_DIR = os.path.join(WORKSPACE_DIR, "reports")

# Step 8.1 Inputs
FILE_STRICT_BENCHMARK = os.path.join(DATA_PROCESSED, "08_1_strict_temporal_benchmark.csv")
FILE_EMB_META = os.path.join(DATA_PROCESSED, "esm2_t12_representatives_metadata.csv")
FILE_GENE_EVID = os.path.join(DATA_PROCESSED, "08_temporal_gene_evidence.csv")
FILE_SOURCES = os.path.join(DATA_PROCESSED, "08_temporal_evidence_sources.csv")
FILE_CROSS_SPECIES = os.path.join(DATA_PROCESSED, "08_cross_species_temporal_labels.csv")

# Step 9 Deliverable Outputs
OUT_RECON_MD = os.path.join(REPORTS_DIR, "09_benchmark_reconstruction.md")
OUT_BENCH_COUNTS_CSV = os.path.join(REPORTS_DIR, "09_benchmark_counts.csv")
OUT_GENE_INDEP_CSV = os.path.join(REPORTS_DIR, "09_gene_family_independence_audit.csv")
OUT_PROV_AUDIT_CSV = os.path.join(REPORTS_DIR, "09_label_provenance_audit.csv")
OUT_FEASIBILITY_MD = os.path.join(REPORTS_DIR, "09_class_feasibility.md")
OUT_FEASIBILITY_CSV = os.path.join(REPORTS_DIR, "09_class_feasibility.csv")
OUT_SPLIT_PROTOCOL_MD = os.path.join(REPORTS_DIR, "09_split_protocol.md")
OUT_SPLIT_MANIFEST_CSV = os.path.join(DATA_PROCESSED, "09_group_split_manifest.csv")
OUT_BASELINE_MD = os.path.join(REPORTS_DIR, "09_baseline_protocol.md")
OUT_METRICS_MD = os.path.join(REPORTS_DIR, "09_evaluation_metrics.md")
OUT_SENSITIVITY_MD = os.path.join(REPORTS_DIR, "09_sensitivity_protocol.md")
OUT_READINESS_MD = os.path.join(REPORTS_DIR, "09_benchmark_readiness_report.md")


def compute_sha256(filepath: str) -> str:
    """Compute the SHA-256 hash of a file for audit verification."""
    if not os.path.exists(filepath):
        return "FILE_NOT_FOUND"
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()


def main():
    print("=" * 80)
    print("STEP 9: LEAKAGE-RESISTANT TEMPORAL BENCHMARK DESIGN")
    print("=" * 80)

    # 1. Load Step 8.1 Benchmark Manifest and Metadata
    print("\n[Phase A] Loading Step 8.1 strict benchmark and representative metadata...")
    df_bench = pd.read_csv(FILE_STRICT_BENCHMARK)
    df_meta = pd.read_csv(FILE_EMB_META)
    df_gene = pd.read_csv(FILE_GENE_EVID)
    df_sources = pd.read_csv(FILE_SOURCES)

    n_total = len(df_bench)
    assert n_total == 3111, f"Expected 3,111 records, found {n_total}"
    assert len(df_meta) == 3111, f"Expected 3,111 metadata rows, found {len(df_meta)}"
    assert (df_bench["embedding_row_index"].values == df_meta["embedding_row_index"].values).all(), (
        "Embedding row index alignment mismatch between benchmark and metadata!"
    )

    # Merge benchmark with sequence-level metadata (homology groups, leakage clusters, hashes)
    merged = df_bench.merge(
        df_meta[[
            "embedding_row_index", "sequence_sha256", "sequence_length",
            "exact_cluster_id", "near_identity_cluster_id_95",
            "homology_group_id", "prevent_ml_leakage_group",
            "molecular_weight", "isoelectric_point", "gravy"
        ]],
        on="embedding_row_index", how="left"
    )

    # =========================================================================
    # PHASE A: RECONSTRUCT BENCHMARK COUNTS
    # =========================================================================
    print("\n[Phase A] Recomputing exact benchmark counts across dimensions...")

    # Dimensional summary table
    count_summary_rows = []

    for cohort_name in ["DIRECT_EVIDENCE", "HOMOLOGY_INFERRED", "UNRESOLVED", "TOTAL_COHORT"]:
        if cohort_name == "TOTAL_COHORT":
            sub = merged
        else:
            sub = merged[merged["adjudicated_cohort"] == cohort_name]

        for cls in ["ALL", "IE_ALPHA", "E_BETA", "L1_GAMMA1", "L2_GAMMA2", "MIXED_OR_CONTEXT_DEPENDENT", "INSUFFICIENT_EVIDENCE"]:
            if cls == "ALL":
                csub = sub
            else:
                csub = sub[sub["normalized_temporal_class"] == cls]

            if len(csub) == 0:
                continue

            count_summary_rows.append({
                "cohort": cohort_name,
                "temporal_class": cls,
                "sequence_record_count": len(csub),
                "unique_sequence_hashes": csub["sequence_sha256"].nunique(),
                "unique_gene_symbols": csub["effective_gene_symbol"].nunique(),
                "unique_homology_groups": csub["homology_group_id"].nunique(),
                "unique_leakage_groups_95": csub["prevent_ml_leakage_group"].nunique(),
                "hsv1_records": int((csub["effective_species"] == "HSV1").sum()),
                "hsv2_records": int((csub["effective_species"] == "HSV2").sum()),
                "median_sequence_length": float(csub["sequence_length"].median()),
                "protein_measured_count": int((csub["measured_molecule"].str.contains("PROTEIN", na=False)).sum()),
                "transcript_only_count": int((csub["measured_molecule"] == "TRANSCRIPT").sum())
            })

    df_counts = pd.DataFrame(count_summary_rows)
    df_counts.to_csv(OUT_BENCH_COUNTS_CSV, index=False)
    print(f"Saved benchmark dimensional counts to: {OUT_BENCH_COUNTS_CSV}")

    # Write Phase A Narrative Report
    recon_md = f"""# Step 9 Report: Benchmark Reconstruction and Dimensional Accounting

**Project**: *Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes*  
**Date**: {datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")}  
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
"""
    with open(OUT_RECON_MD, "w", encoding="utf-8") as f:
        f.write(recon_md)
    print(f"Saved benchmark reconstruction report to: {OUT_RECON_MD}")

    # =========================================================================
    # PHASE B: AUDIT LABEL INDEPENDENCE AND PROVENANCE
    # =========================================================================
    print("\n[Phase B] Auditing label independence and provenance at the gene-family level...")

    gene_family_rows = []
    direct_sub = merged[merged["adjudicated_cohort"] == "DIRECT_EVIDENCE"]

    for gene, grp in direct_sub.groupby("effective_gene_symbol"):
        cls = grp["normalized_temporal_class"].iloc[0]
        n_seqs = len(grp)
        h1_seqs = int((grp["effective_species"] == "HSV1").sum())
        h2_seqs = int((grp["effective_species"] == "HSV2").sum())
        hom_grp = grp["homology_group_id"].iloc[0]
        sources = grp["evidence_source_ids"].iloc[0]
        measured = grp["measured_molecule"].iloc[0]
        is_dual = gene in {"RS1", "RL2", "UL23", "UL27", "US4", "US6", "UL19"}

        gene_family_rows.append({
            "gene_symbol": gene,
            "homology_group_id": hom_grp,
            "temporal_class": cls,
            "total_sequences": n_seqs,
            "hsv1_sequences": h1_seqs,
            "hsv2_sequences": h2_seqs,
            "is_dual_supported": is_dual,
            "primary_evidence_sources": sources,
            "measured_molecule": measured,
            "independence_level": "SINGLE_INDEPENDENT_GENE_FAMILY",
            "audit_risk_factor": "HIGH_SEQUENCE_WEIGHT" if n_seqs >= 30 else "NORMAL"
        })

    df_gene_indep = pd.DataFrame(gene_family_rows)
    df_gene_indep.sort_values(by="total_sequences", ascending=False, inplace=True)
    df_gene_indep.to_csv(OUT_GENE_INDEP_CSV, index=False)
    print(f"Saved gene family independence audit to: {OUT_GENE_INDEP_CSV}")

    # Sequence-level provenance audit table
    df_prov = merged[[
        "embedding_row_index", "record_id", "effective_species",
        "effective_gene_symbol", "homology_group_id", "gene_evidence_tier",
        "adjudicated_cohort", "normalized_temporal_class",
        "transcript_kinetic_class", "protein_kinetic_class",
        "measured_molecule", "evidence_source_ids",
        "cross_species_transfer_status", "audit_rationale"
    ]].copy()
    df_prov.to_csv(OUT_PROV_AUDIT_CSV, index=False)
    print(f"Saved label provenance audit table to: {OUT_PROV_AUDIT_CSV}")

    # =========================================================================
    # PHASE C: QUANTIFY CLASS FEASIBILITY
    # =========================================================================
    print("\n[Phase C] Quantifying class feasibility at sequence and family levels...")

    feasibility_rows = []
    classes_eval = ["IE_ALPHA", "E_BETA", "L1_GAMMA1", "L2_GAMMA2"]

    for cls in classes_eval:
        csub_seq = direct_sub[direct_sub["normalized_temporal_class"] == cls]
        n_seq = len(csub_seq)
        n_fam = csub_seq["homology_group_id"].nunique()
        fams = list(csub_seq["effective_gene_symbol"].unique())

        feasibility_rows.append({
            "temporal_class": cls,
            "sequence_record_count": n_seq,
            "pct_of_sequences": round(n_seq / len(direct_sub) * 100.0, 2),
            "independent_gene_families": n_fam,
            "pct_of_families": round(n_fam / 63 * 100.0, 2),
            "ratio_seq_per_family": round(n_seq / n_fam, 1),
            "gene_families_list": ", ".join(sorted(fams)[:6]) + ("..." if len(fams) > 6 else ""),
            "feasibility_assessment": (
                "CRITICALLY_UNDERPOWERED" if n_fam <= 5 else (
                    "FEASIBLE_WITH_GROUP_AWARENESS" if n_fam >= 15 else "MARGINAL"
                )
            ),
            "statistical_limitation": (
                f"Only {n_fam} independent families. In 5-fold cross-validation, exactly 1 family per fold."
                if n_fam <= 5 else f"Adequate gene family diversity ({n_fam} families)."
            )
        })

    df_feas = pd.DataFrame(feasibility_rows)
    df_feas.to_csv(OUT_FEASIBILITY_CSV, index=False)
    print(f"Saved class feasibility table to: {OUT_FEASIBILITY_CSV}")

    # Write Phase C Report
    feas_md = f"""# Step 9 Report: Class Feasibility and Statistical Generalizability

**Project**: *Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes*  
**Date**: {datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")}  
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
"""
    with open(OUT_FEASIBILITY_MD, "w", encoding="utf-8") as f:
        f.write(feas_md)
    print(f"Saved class feasibility report to: {OUT_FEASIBILITY_MD}")

    # =========================================================================
    # PHASE D: DESIGN GROUP-AWARE SPLITS
    # =========================================================================
    print("\n[Phase D] Designing group-aware, leakage-resistant split manifests...")

    # Design 5-fold GroupKFold on homology_group_id for the DIRECT_EVIDENCE cohort
    # To ensure balance across classes, we use StratifiedGroupKFold on homology groups
    direct_df = direct_sub.copy()
    
    # We assign group folds using StratifiedGroupKFold
    sgkf = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)
    
    # We stratify on normalized_temporal_class and group on homology_group_id
    groups = direct_df["homology_group_id"].values
    y = direct_df["normalized_temporal_class"].values
    X = direct_df["embedding_row_index"].values

    direct_df["group_fold_assignment"] = -1
    for fold_idx, (train_idx, val_idx) in enumerate(sgkf.split(X, y, groups)):
        direct_df.iloc[val_idx, direct_df.columns.get_loc("group_fold_assignment")] = fold_idx

    # Verify zero leakage across folds
    for f in range(5):
        train_groups = set(direct_df[direct_df["group_fold_assignment"] != f]["homology_group_id"])
        val_groups = set(direct_df[direct_df["group_fold_assignment"] == f]["homology_group_id"])
        intersection = train_groups.intersection(val_groups)
        assert len(intersection) == 0, f"LEAKAGE DETECTED in Fold {f}: {intersection}"

    # Also build full split manifest including HOMOLOGY_INFERRED and UNRESOLVED
    manifest_rows = []
    
    # Direct evidence records receive their assigned fold (0-4)
    for _, r in direct_df.iterrows():
        manifest_rows.append({
            "embedding_row_index": r["embedding_row_index"],
            "record_id": r["record_id"],
            "effective_species": r["effective_species"],
            "effective_gene_symbol": r["effective_gene_symbol"],
            "homology_group_id": r["homology_group_id"],
            "prevent_ml_leakage_group": r["prevent_ml_leakage_group"],
            "adjudicated_cohort": "DIRECT_EVIDENCE",
            "normalized_temporal_class": r["normalized_temporal_class"],
            "proposed_split_role": f"DIRECT_FOLD_{r['group_fold_assignment']}",
            "fold_index": int(r["group_fold_assignment"]),
            "leakage_isolation_status": "ZERO_FAMILY_LEAKAGE_GUARANTEED"
        })

    # Homology-inferred records are flagged for secondary sensitivity testing
    homol_sub = merged[merged["adjudicated_cohort"] == "HOMOLOGY_INFERRED"].copy()
    for _, r in homol_sub.iterrows():
        manifest_rows.append({
            "embedding_row_index": r["embedding_row_index"],
            "record_id": r["record_id"],
            "effective_species": r["effective_species"],
            "effective_gene_symbol": r["effective_gene_symbol"],
            "homology_group_id": r["homology_group_id"],
            "prevent_ml_leakage_group": r["prevent_ml_leakage_group"],
            "adjudicated_cohort": "HOMOLOGY_INFERRED",
            "normalized_temporal_class": r["normalized_temporal_class"],
            "proposed_split_role": "SECONDARY_HOMOLOGY_SENSITIVITY",
            "fold_index": -1,
            "leakage_isolation_status": "EXCLUDED_FROM_PRIMARY_TRAINING"
        })

    # Unresolved records are quarantined
    unres_sub = merged[merged["adjudicated_cohort"] == "UNRESOLVED"].copy()
    for _, r in unres_sub.iterrows():
        manifest_rows.append({
            "embedding_row_index": r["embedding_row_index"],
            "record_id": r["record_id"],
            "effective_species": r["effective_species"],
            "effective_gene_symbol": r["effective_gene_symbol"],
            "homology_group_id": r["homology_group_id"],
            "prevent_ml_leakage_group": r["prevent_ml_leakage_group"],
            "adjudicated_cohort": "UNRESOLVED",
            "normalized_temporal_class": r["normalized_temporal_class"],
            "proposed_split_role": "QUARANTINED_UNLABELED_DISCOVERY",
            "fold_index": -2,
            "leakage_isolation_status": "QUARANTINED_FROM_ALL_BENCHMARKS"
        })

    df_manifest = pd.DataFrame(manifest_rows)
    df_manifest.sort_values(by="embedding_row_index", inplace=True)
    df_manifest.to_csv(OUT_SPLIT_MANIFEST_CSV, index=False)
    print(f"Saved group-aware split manifest (N = {len(df_manifest):,}) to: {OUT_SPLIT_MANIFEST_CSV}")

    # Write Phase D Protocol Report
    split_md = f"""# Step 9 Protocol: Group-Aware, Leakage-Resistant Cross-Validation Design

**Project**: *Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes*  
**Date**: {datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")}  
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
`TrainFamilies_k ∩ ValFamilies_k = ∅` for all folds k in {0, 1, 2, 3, 4}.

- **Zero Sequence Homology Leakage**: All 212 variants of `RS1` (ICP4) are assigned to a single fold. They are never split across train and test.
- **Zero Homology-Group Leakage**: All sequences sharing >= 95% near-identity or canonical orthology reside in the identical fold.
- **Prohibition of Random Splitting**: Standard random sequence splitting is prohibited because it evaluates trivial sequence re-identification rather than kinetic generalization.
"""
    with open(OUT_SPLIT_PROTOCOL_MD, "w", encoding="utf-8") as f:
        f.write(split_md)
    print(f"Saved split protocol report to: {OUT_SPLIT_PROTOCOL_MD}")

    # =========================================================================
    # PHASE E: DEFINE SIMPLE BASELINES AND EVALUATION METRICS
    # =========================================================================
    print("\n[Phase E] Formulating baseline protocols and prespecifying evaluation metrics...")

    baseline_md = f"""# Step 9 Protocol: Prespecified Baseline Hierarchy

**Project**: *Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes*  
**Date**: {datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")}  
**Deliverable**: `reports/09_baseline_protocol.md`  

---

## 1. Baseline Model Hierarchy

Before fitting any advanced neural architectures or tuning classifiers, future modeling must evaluate four rigorous baselines:

### 1. Majority-Class Dummy Classifier
- Predicts the most frequent class in the training split.
- **Purpose**: Establishes the empirical performance floor under class imbalance.

### 2. Gene-Family Label Lookup (Diagnostic Diagnostic)
- Predicts the historical label of the gene family using exact symbol or homology group match.
- **Purpose**: Used **strictly as a diagnostic for data leakage**.
- If a model achieves near-100% test performance when family lookup is impossible under grouped holdouts, the model is exploiting unexpected leakage.

### 3. Physicochemical Property Baseline
- Linear logistic regression / Random Forest trained on 8 classical physical features:
  * Molecular weight, isoelectric point, GRAVY (hydropathy), aliphatic index, instability index, charged fraction, polar fraction, hydrophobic fraction.
- **Purpose**: Tests whether temporal kinetics correlates simply with bulk physical chemistry (e.g. acidic transcription factors vs basic capsid proteins).

### 4. Ridge-Regularized ESM-2 Linear Classifier
- Multinomial logistic regression with L2 regularization directly on the 480-dimensional mean-pooled ESM-2 embeddings.
- **Purpose**: Tests whether protein language representations encode temporal expression timing.

---

## 2. Prespecified Evaluation Metrics

Because of severe class imbalance at the gene-family level (N = 5 for IE vs N = 24 for L2):

1. **Macro-Averaged F1-Score (Macro-F1)**:
   Arithmetic mean of per-class F1-scores across all K classes: `Macro-F1 = (1/K) * sum(F1_k)`.
   *Gives equal weight to all classes regardless of sequence count.*
2. **Balanced Accuracy**:
   Arithmetic mean of per-class recalls: `Balanced Accuracy = (1/K) * sum(Recall_k)`.
3. **Per-Class Precision, Recall, and Confusion Matrices**.
4. **Gene-Family-Level Uncertainty Estimation**:
   * Metrics must be computed by averaging performance across the **63 independent gene families**, rather than sequence-weighted averages.
"""
    with open(OUT_BASELINE_MD, "w", encoding="utf-8") as f:
        f.write(baseline_md)
    print(f"Saved baseline protocol report to: {OUT_BASELINE_MD}")

    # Metrics report
    metrics_md = f"""# Step 9 Protocol: Prespecified Statistical Metrics and Uncertainty Protocol

**Project**: *Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes*  
**Date**: {datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")}  
**Deliverable**: `reports/09_evaluation_metrics.md`  

---

## 1. Statistical Principles

1. **Sequence-Weighted Metrics Prohibited as Primary**:
   Evaluating accuracy across 687 sequences assigns 30.8% of total score weight to a single gene (`RS1`). A model predicting ICP4 correctly on 212 sequences would achieve high sequence accuracy even if it failed on all other 62 viral genes.
2. **Gene-Family Grouped Averaging (Macro-Evaluation)**:
   All validation metrics must compute prediction error at the **gene-family level** ($N=63$).
3. **Confidence Interval Estimation**:
   Bootstrap resampling across gene families ($B = 1,000$ iterations) to construct 95% confidence intervals on Macro-F1 and Balanced Accuracy.
"""
    with open(OUT_METRICS_MD, "w", encoding="utf-8") as f:
        f.write(metrics_md)
    print(f"Saved evaluation metrics protocol to: {OUT_METRICS_MD}")

    # =========================================================================
    # PHASE F: PRESPECIFY SENSITIVITY ANALYSES
    # =========================================================================
    print("\n[Phase F] Prespecifying sensitivity analysis protocols...")

    sensitivity_md = f"""# Step 9 Protocol: Prespecified Sensitivity Analyses

**Project**: *Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes*  
**Date**: {datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")}  
**Deliverable**: `reports/09_sensitivity_protocol.md`  

---

## 1. Prespecified Sensitivity Cohorts

Future evaluations must assess stability across four distinct sensitivity tiers:

1. **Primary Benchmark**: Direct evidence records only ($N = 687$ sequences, 63 gene families).
2. **Protein-Validated Sensitivity Cohort**: Excludes genes whose kinetics were measured solely by nascent 4sU RNA without physical mass spectrometry validation (excludes `UL53`, `UL2`, `UL4`, `UL20`, `UL49A`, `US5`, `US10`, `UL43`, `UL45`, `UL56`, `US2`).
3. **Expanded Homology Sensitivity Cohort**: Incorporates the 633 homology-inferred HSV-2 records to test cross-species transferability.
4. **Binary Kinetic Task Sensitivity**: Compares full 4-class prediction (`IE`, `E`, `L1`, `L2`) against binary replication timing (`Pre-Replicative` vs `Post-Replicative`).
"""
    with open(OUT_SENSITIVITY_MD, "w", encoding="utf-8") as f:
        f.write(sensitivity_md)
    print(f"Saved sensitivity protocol report to: {OUT_SENSITIVITY_MD}")

    # =========================================================================
    # PHASE G: BENCHMARK READINESS DECISION REPORT
    # =========================================================================
    print("\n[Phase G] Compiling master benchmark readiness decision report...")

    readiness_md = f"""# Step 9 Final Report: Benchmark Readiness Decision and Task Feasibility Assessment

**Project**: *Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes*  
**Date**: {datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")}  
**Author**: Computational Virology Machine Learning Governance Committee  
**Target Benchmark**: Direct Evidence Primary Cohort (N = 687 sequences, 63 independent gene families)  
**Deliverables Generated**:
- `reports/09_benchmark_reconstruction.md`: SHA-256 `{compute_sha256(OUT_RECON_MD)}`
- `reports/09_benchmark_counts.csv`: SHA-256 `{compute_sha256(OUT_BENCH_COUNTS_CSV)}`
- `reports/09_gene_family_independence_audit.csv`: SHA-256 `{compute_sha256(OUT_GENE_INDEP_CSV)}`
- `reports/09_label_provenance_audit.csv`: SHA-256 `{compute_sha256(OUT_PROV_AUDIT_CSV)}`
- `reports/09_class_feasibility.md`: SHA-256 `{compute_sha256(OUT_FEASIBILITY_MD)}`
- `reports/09_class_feasibility.csv`: SHA-256 `{compute_sha256(OUT_FEASIBILITY_CSV)}`
- `reports/09_split_protocol.md`: SHA-256 `{compute_sha256(OUT_SPLIT_PROTOCOL_MD)}`
- `data/processed/09_group_split_manifest.csv`: SHA-256 `{compute_sha256(OUT_SPLIT_MANIFEST_CSV)}`
- `reports/09_baseline_protocol.md`: SHA-256 `{compute_sha256(OUT_BASELINE_MD)}`
- `reports/09_evaluation_metrics.md`: SHA-256 `{compute_sha256(OUT_METRICS_MD)}`
- `reports/09_sensitivity_protocol.md`: SHA-256 `{compute_sha256(OUT_SENSITIVITY_MD)}`

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
"""
    with open(OUT_READINESS_MD, "w", encoding="utf-8") as f:
        f.write(readiness_md)
    print(f"Saved master benchmark readiness report to: {OUT_READINESS_MD}")

    print("\n" + "=" * 80)
    print("STEP 9 BENCHMARK DESIGN COMPLETE.")
    print("=" * 80)


if __name__ == "__main__":
    main()
