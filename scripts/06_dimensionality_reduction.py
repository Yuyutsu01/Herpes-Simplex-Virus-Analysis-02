#!/usr/bin/env python3
"""
Step 6 Phase 3: Dimensionality Reduction Analysis (PCA, UMAP, and t-SNE)
========================================================================

Conceptual Overview:
--------------------
ESM-2 embeddings map protein sequences into a continuous 480-dimensional latent
space learned from masked language modeling over UniRef50. To inspect how the
HSV-1 and HSV-2 proteomes are organized within this space:

1. Principal Component Analysis (PCA):
   Identifies the dominant linear orthogonal axes of variance. PCA reveals the
   global variance spectrum (scree plot), determines the effective dimensionality
   of the proteome embedding space, and tests whether linear components capture
   broad biological categories (e.g. envelope glycoproteins vs capsids vs tegument)
   or technical artifacts (e.g. sequence length or ambiguity).

2. Uniform Manifold Approximation and Projection (UMAP):
   Captures non-linear local neighborhood topology and manifold structure. We evaluate
   hyperparameter stability across n_neighbors in [15, 30, 50] and min_dist in [0.1, 0.3].
   We explicitly evaluate whether UMAP islands persist in the high-quality sensitivity
   cohort (Cohort B, <=5% ambiguous residues) to rule out sequencing-gap artifacts.

3. t-Distributed Stochastic Neighbor Embedding (t-SNE):
   Provides a non-linear local density representation to complement UMAP and verify
   that observed clusters are not artifacts of a specific graph layout algorithm.
"""

import os
import sys
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(line_buffering=True)
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # Non-interactive backend for headless execution
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap

WORKSPACE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA_PROCESSED = os.path.join(WORKSPACE_DIR, "data", "processed")
REPORTS_DIR = os.path.join(WORKSPACE_DIR, "reports")
FIGURES_DIR = os.path.join(REPORTS_DIR, "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)

FILE_EMB_NPY = os.path.join(DATA_PROCESSED, "esm2_t12_representatives_embeddings.npy")
FILE_EMB_META = os.path.join(DATA_PROCESSED, "esm2_t12_representatives_metadata.csv")
OUT_REDUCTION_CSV = os.path.join(DATA_PROCESSED, "06_dimensionality_reduction_coordinates.csv")

RANDOM_SEED = 42


def run_dimensionality_reduction():
    print("=" * 80)
    print("STEP 6 PHASE 3: DIMENSIONALITY REDUCTION (PCA, UMAP, t-SNE)")
    print("=" * 80)

    # 1. Load Embeddings and Metadata
    print(f"Loading embeddings: {FILE_EMB_NPY}")
    embeddings = np.load(FILE_EMB_NPY)
    df_meta = pd.read_csv(FILE_EMB_META)
    num_records, embedding_dim = embeddings.shape
    print(f"  Loaded {num_records:,} sequences with embedding dimension {embedding_dim}")

    # 2. Define Cohorts
    # Cohort A: Full Representative Cohort (all 3,111 sequences)
    # Cohort B: High-Quality Sensitivity Cohort (<= 5% ambiguous residues, removes severe gap artifacts)
    df_meta["is_high_quality"] = df_meta["ambiguous_residue_fraction"] <= 0.05
    n_cohort_a = len(df_meta)
    n_cohort_b = int(df_meta["is_high_quality"].sum())
    print(f"\nCohort A (Full Representatives): {n_cohort_a:,} records")
    print(f"Cohort B (High-Quality Sensitivity, ambig <= 0.05): {n_cohort_b:,} records ({n_cohort_b/n_cohort_a*100:.1f}%)")
    print(f"  Removed {n_cohort_a - n_cohort_b:,} sequences with > 5% ambiguous residues (e.g. X-runs in RL1/ICP34.5)")

    # 3. Standardization (without outcome labels)
    print("\nStandardizing 480-dimensional embeddings (zero mean, unit variance)...")
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(embeddings)

    # 4. PCA Analysis (Top 50 Components)
    print("\n--- 4A. Computing PCA (n_components=50) ---")
    pca = PCA(n_components=50, random_state=RANDOM_SEED)
    X_pca = pca.fit_transform(X_scaled)

    evr = pca.explained_variance_ratio_
    cum_evr = np.cumsum(evr)
    print(f"  PC1 Explained Variance: {evr[0]*100:.2f}%")
    print(f"  PC2 Explained Variance: {evr[1]*100:.2f}%")
    print(f"  Top 10 PCs Cumulative Explained Variance: {cum_evr[9]*100:.2f}%")
    print(f"  Top 50 PCs Cumulative Explained Variance: {cum_evr[49]*100:.2f}%")

    # Figure 1: PCA Explained Variance (Scree Plot)
    plt.figure(figsize=(10, 5))
    plt.subplot(1, 2, 1)
    plt.plot(np.arange(1, 51), evr * 100, "o-", color="#1f77b4", markersize=4)
    plt.title("PCA Scree Plot (Top 50 PCs)")
    plt.xlabel("Principal Component")
    plt.ylabel("Individual Explained Variance (%)")
    plt.grid(True, linestyle="--", alpha=0.5)

    plt.subplot(1, 2, 2)
    plt.plot(np.arange(1, 51), cum_evr * 100, "s-", color="#2ca02c", markersize=4)
    plt.axhline(y=50, color="r", linestyle="--", label="50% Threshold")
    plt.title("Cumulative Explained Variance")
    plt.xlabel("Principal Component")
    plt.ylabel("Cumulative Explained Variance (%)")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    fig1_path = os.path.join(FIGURES_DIR, "01_pca_explained_variance.png")
    plt.savefig(fig1_path, dpi=300)
    plt.close()
    print(f"  Saved: {fig1_path}")

    # Figure 2: PCA Projection colored by Species (HSV1 vs HSV2)
    plt.figure(figsize=(8, 7))
    palette_species = {"HSV1": "#1f77b4", "HSV2": "#ff7f0e"}
    sns.scatterplot(
        x=X_pca[:, 0], y=X_pca[:, 1],
        hue=df_meta["effective_species_attribution"],
        palette=palette_species, alpha=0.6, s=25, edgecolor="none"
    )
    plt.title(f"PCA Projection: PC1 vs PC2 by Species (N={n_cohort_a:,})\nPC1: {evr[0]*100:.2f}% var, PC2: {evr[1]*100:.2f}% var")
    plt.xlabel("Principal Component 1")
    plt.ylabel("Principal Component 2")
    plt.legend(title="Species")
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    fig2_path = os.path.join(FIGURES_DIR, "02_pca_species_projection.png")
    plt.savefig(fig2_path, dpi=300)
    plt.close()
    print(f"  Saved: {fig2_path}")

    # Figure 3: PCA Projection colored by Major Gene Categories
    # Define broad structural/functional categories for HSV genes
    def categorize_gene(gene):
        if gene == "UNANNOTATED":
            return "Unannotated"
        # Glycoproteins (Envelope)
        if gene in ["UL27", "UL44", "UL53", "US4", "US6", "US7", "US8", "UL49A", "US5"]:
            return "Glycoprotein (Env)"
        # Capsid
        if gene in ["UL18", "UL19", "UL26", "UL26.5", "UL35", "UL38"]:
            return "Capsid Structural"
        # Tegument
        if gene in ["UL11", "UL16", "UL21", "UL36", "UL37", "UL46", "UL47", "UL48", "UL49"]:
            return "Tegument"
        # Replication / Enzymes
        if gene in ["UL5", "UL8", "UL9", "UL23", "UL29", "UL30", "UL42", "UL50", "UL52"]:
            return "Replication / Enzyme"
        # Regulatory (Immediate-Early / Latency)
        if gene in ["RS1", "RL2", "US1", "US11", "US12", "RL1"]:
            return "Regulatory / Latency"
        return "Other Canonical"

    df_meta["broad_category"] = df_meta["gene_annotation"].apply(categorize_gene)
    category_counts = df_meta["broad_category"].value_counts()
    print("\nBroad Gene Category Counts in Representative Cohort:")
    for cat, count in category_counts.items():
        print(f"  - {cat}: {count:,} ({count/n_cohort_a*100:.1f}%)")

    plt.figure(figsize=(9, 7))
    sns.scatterplot(
        x=X_pca[:, 0], y=X_pca[:, 1],
        hue=df_meta["broad_category"],
        palette="tab10", alpha=0.7, s=28, edgecolor="none"
    )
    plt.title(f"PCA Projection: PC1 vs PC2 by Gene Category (N={n_cohort_a:,})")
    plt.xlabel("Principal Component 1")
    plt.ylabel("Principal Component 2")
    plt.legend(title="Gene Category", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    fig3_path = os.path.join(FIGURES_DIR, "03_pca_gene_family_projection.png")
    plt.savefig(fig3_path, dpi=300)
    plt.close()
    print(f"  Saved: {fig3_path}")

    # 5. UMAP Analysis Across Multiple Hyperparameters
    print("\n--- 4B. Computing UMAP Across Parameter Grid ---")
    umap_grid = [
        {"n_neighbors": 15, "min_dist": 0.1, "label": "nn15_md01"},
        {"n_neighbors": 30, "min_dist": 0.1, "label": "nn30_md01_primary"},
        {"n_neighbors": 50, "min_dist": 0.3, "label": "nn50_md03"}
    ]

    umap_coords = {}
    for p in umap_grid:
        print(f"  Running UMAP with n_neighbors={p['n_neighbors']}, min_dist={p['min_dist']}...")
        reducer = umap.UMAP(
            n_neighbors=p["n_neighbors"],
            min_dist=p["min_dist"],
            n_components=2,
            metric="cosine",  # Cosine metric aligns with transformer latent spaces
            random_state=RANDOM_SEED
        )
        embedding_umap = reducer.fit_transform(X_scaled)
        umap_coords[p["label"]] = embedding_umap
        df_meta[f"umap1_{p['label']}"] = embedding_umap[:, 0]
        df_meta[f"umap2_{p['label']}"] = embedding_umap[:, 1]

    # Primary UMAP coordinates (nn=30, md=0.1, cosine metric)
    primary_umap = umap_coords["nn30_md01_primary"]

    # Figure 4: UMAP Projection colored by Species (HSV1 vs HSV2)
    plt.figure(figsize=(8, 7))
    sns.scatterplot(
        x=primary_umap[:, 0], y=primary_umap[:, 1],
        hue=df_meta["effective_species_attribution"],
        palette=palette_species, alpha=0.6, s=25, edgecolor="none"
    )
    plt.title(f"UMAP Projection (n_neighbors=30, min_dist=0.1, metric=cosine)\nColored by Species (N={n_cohort_a:,})")
    plt.xlabel("UMAP Dimension 1")
    plt.ylabel("UMAP Dimension 2")
    plt.legend(title="Species")
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    fig4_path = os.path.join(FIGURES_DIR, "04_umap_species_projection.png")
    plt.savefig(fig4_path, dpi=300)
    plt.close()
    print(f"  Saved: {fig4_path}")

    # Figure 6: UMAP Projection colored by Sequence Quality (Ambiguous Residue Fraction)
    plt.figure(figsize=(8.5, 7))
    scatter = plt.scatter(
        primary_umap[:, 0], primary_umap[:, 1],
        c=df_meta["ambiguous_residue_fraction"],
        cmap="viridis", alpha=0.7, s=25, edgecolors="none"
    )
    cbar = plt.colorbar(scatter)
    cbar.set_label("Ambiguous Residue Fraction (X / total)")
    plt.title(f"UMAP Projection Colored by Sequence Ambiguity Fraction (N={n_cohort_a:,})\nEvaluating Sequencing-Gap Artifact Clustering")
    plt.xlabel("UMAP Dimension 1")
    plt.ylabel("UMAP Dimension 2")
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    fig6_path = os.path.join(FIGURES_DIR, "06_umap_quality_projection.png")
    plt.savefig(fig6_path, dpi=300)
    plt.close()
    print(f"  Saved: {fig6_path}")

    # Figure 8: Sensitivity Comparison: Cohort A (All) vs Cohort B (High-Quality <= 5% ambig)
    mask_hq = df_meta["is_high_quality"].values
    plt.figure(figsize=(14, 6))
    plt.subplot(1, 2, 1)
    sns.scatterplot(
        x=primary_umap[:, 0], y=primary_umap[:, 1],
        hue=df_meta["effective_species_attribution"],
        palette=palette_species, alpha=0.5, s=20, edgecolor="none"
    )
    plt.title(f"Cohort A: All Representatives (N={n_cohort_a:,})\nIncludes Sequences with Gap Masking")
    plt.xlabel("UMAP Dimension 1")
    plt.ylabel("UMAP Dimension 2")
    plt.legend(title="Species")
    plt.grid(True, linestyle="--", alpha=0.4)

    plt.subplot(1, 2, 2)
    sns.scatterplot(
        x=primary_umap[mask_hq, 0], y=primary_umap[mask_hq, 1],
        hue=df_meta.loc[mask_hq, "effective_species_attribution"],
        palette=palette_species, alpha=0.6, s=25, edgecolor="none"
    )
    plt.title(f"Cohort B: High-Quality Subset (N={n_cohort_b:,})\nExcludes Sequences with >5% Ambiguous Residues")
    plt.xlabel("UMAP Dimension 1")
    plt.ylabel("UMAP Dimension 2")
    plt.legend(title="Species")
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    fig8_path = os.path.join(FIGURES_DIR, "08_sensitivity_comparison_cohort_a_vs_b.png")
    plt.savefig(fig8_path, dpi=300)
    plt.close()
    print(f"  Saved: {fig8_path}")

    # 6. Optional t-SNE Analysis
    print("\n--- 4C. Computing t-SNE (perplexity=30 on 50 PCs) ---")
    tsne = TSNE(n_components=2, perplexity=30, random_state=RANDOM_SEED, max_iter=1000)
    X_tsne = tsne.fit_transform(X_pca)
    df_meta["tsne1"] = X_tsne[:, 0]
    df_meta["tsne2"] = X_tsne[:, 1]
    print("  t-SNE computation completed.")

    # 7. Store Principal Components in Output Table
    for i in range(10):
        df_meta[f"pca{i+1}"] = X_pca[:, i]

    # Save comprehensive dimensionality reduction coordinate table
    print(f"\nSaving dimensionality reduction coordinates to: {OUT_REDUCTION_CSV}")
    df_meta.to_csv(OUT_REDUCTION_CSV, index=False)
    print(f"  Saved {len(df_meta):,} rows and {len(df_meta.columns)} columns.")

    print("\n================================================================================")
    print("STEP 6 PHASE 3 COMPLETE: DIMENSIONALITY REDUCTION PIPELINE SUCCEEDED.")
    print("================================================================================")


if __name__ == "__main__":
    run_dimensionality_reduction()
