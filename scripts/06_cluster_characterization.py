#!/usr/bin/env python3
"""
Step 6 Phase 5: Physicochemical Characterization and Statistical Correlation
=============================================================================

Conceptual Overview:
--------------------
While ESM-2 representations are derived from unsupervised masked language modeling
over amino acid sequences, protein structures and evolutionary constraints are
fundamentally constrained by physical chemistry (mass, charge, hydropathicity,
and structural stability).

To characterize the clusters without circular reasoning or causal over-claims:
1. Feature Profiling: For every discovered cluster, we compute the median and
   interquartile range (IQR) across length, molecular weight, isoelectric point (pI),
   GRAVY (hydropathicity), instability index, aliphatic index, and residue group fractions.
2. Dominant Annotation & Purity: We assess the primary gene family associated with each
   cluster and determine cluster purity (% belonging to dominant gene family).
3. Statistical Correlation with Embedding Geometry: We compute Spearman rank
   correlations between the top 5 principal components and the 10 physicochemical
   variables, strictly adjusting for multiple hypothesis testing using the
   conservative Bonferroni correction (m = 50 hypothesis tests, alpha = 0.05 / 50 = 0.001).
4. Strict Scientific Guardrails: Physicochemical correlations describe the biophysical
   properties associated with latent dimensions; THEY DO NOT PROVE CAUSATION or provide
   evidence for temporal gene expression classes.
"""

import os
import sys
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(line_buffering=True)
import json
import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

WORKSPACE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA_PROCESSED = os.path.join(WORKSPACE_DIR, "data", "processed")
REPORTS_DIR = os.path.join(WORKSPACE_DIR, "reports")
FIGURES_DIR = os.path.join(REPORTS_DIR, "figures")

FILE_MEMBERSHIP_CSV = os.path.join(REPORTS_DIR, "06_cluster_membership.csv")
OUT_CHARACTERIZATION_CSV = os.path.join(REPORTS_DIR, "06_cluster_characterization.csv")

PHYSICO_COLS = [
    "sequence_length", "molecular_weight", "isoelectric_point",
    "instability_index", "aliphatic_index", "gravy",
    "charged_fraction", "polar_fraction", "hydrophobic_fraction",
    "aromatic_fraction"
]


def run_characterization():
    print("=" * 80)
    print("STEP 6 PHASE 5: CLUSTER PHYSICOCHEMICAL CHARACTERIZATION")
    print("=" * 80)

    # 1. Load Cluster Membership Data
    df = pd.read_csv(FILE_MEMBERSHIP_CSV)
    print(f"Loaded cluster membership table: {len(df):,} records")

    # 2. Cluster-Level Summaries for Primary Clustering Model (K-Means k=74)
    cluster_summaries = []
    unique_clusters = sorted(df["primary_cluster_id"].unique())
    print(f"\nComputing physicochemical profiles across {len(unique_clusters)} clusters...")

    for cid in unique_clusters:
        sub = df[df["primary_cluster_id"] == cid]
        size = len(sub)
        
        # Species composition
        n_hsv1 = int((sub["effective_species_attribution"] == "HSV1").sum())
        n_hsv2 = int((sub["effective_species_attribution"] == "HSV2").sum())
        pct_hsv1 = round(n_hsv1 / size * 100.0, 1)

        # Gene family composition
        gene_counts = sub["gene_annotation"].value_counts()
        top_gene = gene_counts.index[0]
        top_gene_count = int(gene_counts.iloc[0])
        purity_overall = round(top_gene_count / size * 100.0, 1)

        # Annotated subset purity
        ann_sub = sub[sub["gene_annotation"] != "UNANNOTATED"]
        if len(ann_sub) > 0:
            ann_counts = ann_sub["gene_annotation"].value_counts()
            top_ann_gene = ann_counts.index[0]
            ann_purity = round(ann_counts.iloc[0] / len(ann_sub) * 100.0, 1)
        else:
            top_ann_gene = "NONE"
            ann_purity = 0.0

        # Physicochemical medians and IQRs
        row_dict = {
            "cluster_id": cid,
            "cluster_size": size,
            "hsv1_count": n_hsv1,
            "hsv2_count": n_hsv2,
            "pct_hsv1": pct_hsv1,
            "top_gene": top_gene,
            "top_gene_count": top_gene_count,
            "cluster_purity_pct": purity_overall,
            "top_annotated_gene": top_ann_gene,
            "annotated_purity_pct": ann_purity,
            "median_length": round(float(sub["sequence_length"].median()), 1),
            "median_mw_kda": round(float(sub["molecular_weight"].median()) / 1000.0, 2),
            "median_pi": round(float(sub["isoelectric_point"].median()), 2),
            "median_gravy": round(float(sub["gravy"].median()), 3),
            "median_instability": round(float(sub["instability_index"].median()), 2),
            "median_aliphatic": round(float(sub["aliphatic_index"].median()), 2),
            "median_charged_pct": round(float(sub["charged_fraction"].median()) * 100.0, 1),
            "median_hydrophobic_pct": round(float(sub["hydrophobic_fraction"].median()) * 100.0, 1),
            "median_aromatic_pct": round(float(sub["aromatic_fraction"].median()) * 100.0, 1)
        }
        cluster_summaries.append(row_dict)

    df_char = pd.DataFrame(cluster_summaries)
    df_char.to_csv(OUT_CHARACTERIZATION_CSV, index=False)
    print(f"Saved cluster characterization table to: {OUT_CHARACTERIZATION_CSV}")
    print(f"  Processed {len(df_char)} clusters. Mean cluster purity: {df_char['cluster_purity_pct'].mean():.1f}%")

    # 3. Correlation between Principal Components and Physicochemical Features
    print("\n--- Computing Spearman Rank Correlations (PCs vs Physicochemical Features) ---")
    pc_cols = [f"pca{i+1}" for i in range(5)]
    m_tests = len(pc_cols) * len(PHYSICO_COLS)
    alpha_bonferroni = 0.05 / m_tests
    print(f"  Hypothesis tests: {m_tests} (Bonferroni significance threshold: alpha = {alpha_bonferroni:.6f})")

    corr_matrix = pd.DataFrame(index=pc_cols, columns=PHYSICO_COLS, dtype=float)
    pval_matrix = pd.DataFrame(index=pc_cols, columns=PHYSICO_COLS, dtype=float)

    for pc in pc_cols:
        for feat in PHYSICO_COLS:
            # Drop NaNs if any (handled gracefully)
            valid_mask = ~(df[pc].isna() | df[feat].isna())
            rho, pval = stats.spearmanr(df.loc[valid_mask, pc], df.loc[valid_mask, feat])
            corr_matrix.loc[pc, feat] = round(rho, 3)
            pval_matrix.loc[pc, feat] = pval

    print("\nSpearman Correlation Matrix (PC1-PC5 vs Features):")
    print(corr_matrix)

    # Report statistically significant correlations
    print("\nStatistically Significant Associations after Bonferroni Correction:")
    sig_count = 0
    for pc in pc_cols:
        for feat in PHYSICO_COLS:
            pval = pval_matrix.loc[pc, feat]
            rho = corr_matrix.loc[pc, feat]
            if pval < alpha_bonferroni:
                sig_count += 1
                print(f"  - {pc} vs {feat:20s}: rho = {rho:+.3f} (p = {pval:.2e} < {alpha_bonferroni:.6f})")
    print(f"  Total significant associations: {sig_count}/{m_tests}")

    # 4. Generate Figure 7: Cluster-Level Physicochemical Comparison Plots
    print("\n--- Generating Figure 7: Physicochemical Comparison Plots ---")
    plt.figure(figsize=(14, 10))

    # Panel A: Isoelectric Point (pI) by Broad Gene Category
    plt.subplot(2, 2, 1)
    sns.boxplot(
        data=df, x="broad_category", y="isoelectric_point",
        palette="Set2", showfliers=False
    )
    plt.xticks(rotation=30, ha="right")
    plt.title("A. Theoretical Isoelectric Point (pI) by Functional Category")
    plt.xlabel("")
    plt.ylabel("Isoelectric Point (pH)")
    plt.grid(True, linestyle="--", alpha=0.4)

    # Panel B: GRAVY (Hydropathicity) by Broad Gene Category
    plt.subplot(2, 2, 2)
    sns.boxplot(
        data=df, x="broad_category", y="gravy",
        palette="Set2", showfliers=False
    )
    plt.axhline(0, color="gray", linestyle="--", alpha=0.7)
    plt.xticks(rotation=30, ha="right")
    plt.title("B. Hydropathicity Index (GRAVY) by Functional Category\n(Positive = Hydrophobic / Membrane-associated)")
    plt.xlabel("")
    plt.ylabel("GRAVY Score")
    plt.grid(True, linestyle="--", alpha=0.4)

    # Panel C: Molecular Weight (kDa) vs Sequence Length
    plt.subplot(2, 2, 3)
    sns.scatterplot(
        data=df, x="sequence_length", y=df["molecular_weight"] / 1000.0,
        hue="broad_category", palette="tab10", alpha=0.6, s=25, edgecolor="none"
    )
    plt.title("C. Protein Length vs Molecular Weight (kDa)")
    plt.xlabel("Sequence Length (amino acids)")
    plt.ylabel("Molecular Weight (kDa)")
    plt.legend([],[], frameon=False)  # Suppress redundant legend
    plt.grid(True, linestyle="--", alpha=0.4)

    # Panel D: Instability Index vs Aliphatic Index
    plt.subplot(2, 2, 4)
    sns.scatterplot(
        data=df, x="aliphatic_index", y="instability_index",
        hue="broad_category", palette="tab10", alpha=0.6, s=25, edgecolor="none"
    )
    plt.axhline(40, color="red", linestyle="--", label="Instability Threshold (40)")
    plt.title("D. Aliphatic Index vs Instability Index\n(>40 indicates predicted in vitro instability)")
    plt.xlabel("Aliphatic Index")
    plt.ylabel("Instability Index")
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", title="Category")
    plt.grid(True, linestyle="--", alpha=0.4)

    plt.tight_layout()
    fig7_path = os.path.join(FIGURES_DIR, "07_cluster_physicochemical_distributions.png")
    plt.savefig(fig7_path, dpi=300)
    plt.close()
    print(f"  Saved: {fig7_path}")

    print("\n================================================================================")
    print("STEP 6 PHASE 5 COMPLETE: PHYSICOCHEMICAL INTEGRATION SUCCEEDED.")
    print("================================================================================")


if __name__ == "__main__":
    run_characterization()
