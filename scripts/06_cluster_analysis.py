#!/usr/bin/env python3
"""
Step 6 Phase 4: Rigorous Protein Cluster Analysis
=================================================

Conceptual Overview:
--------------------
Clustering protein representations allows us to discover whether proteins naturally
partition into distinct structural, functional, or evolutionary modules in the
ESM-2 embedding space, without relying on subjective visual inspection of 2D plots.

Methodological Principles:
1. Primary Representation Space: Clustering is performed directly in the validated
   standardized embedding space (480 dimensions, or top 50 PCs preserving 96.11% variance),
   NOT on the 2D UMAP projection. Non-linear 2D projections warp distances and densities;
   clustering in 2D would introduce severe artificial boundary artifacts.
2. Multi-Algorithm Evaluation:
   - K-Means: Evaluated across k in [10, 20, 30, 40, 50, 74, 100]. (k=74 is biologically
     grounded, matching the 74 canonical HSV core gene families).
   - Agglomerative Hierarchical Clustering: Evaluated with Ward linkage across k in [20, 50, 74].
   - HDBSCAN: Density-based clustering (min_cluster_size in [15, 25, 50]), explicitly identifying
     dense protein cores while labeling divergent or fragmented provisional candidates as noise (-1).
3. Comprehensive Metrics:
   - Internal geometry: Silhouette Score, Davies-Bouldin Index (DBI), Calinski-Harabasz Index (CHI).
   - Annotation concordance: Adjusted Rand Index (ARI) and Normalized Mutual Information (NMI)
     measured against known canonical gene families (Cohort C).
   - Sensitivity: Metric recalculation on Cohort B (high quality, <=5% ambiguity).
"""

import os
import sys
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(line_buffering=True)
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, AgglomerativeClustering, HDBSCAN
from sklearn.metrics import silhouette_score, davies_bouldin_score, calinski_harabasz_score
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

WORKSPACE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA_PROCESSED = os.path.join(WORKSPACE_DIR, "data", "processed")
REPORTS_DIR = os.path.join(WORKSPACE_DIR, "reports")
FIGURES_DIR = os.path.join(REPORTS_DIR, "figures")

FILE_EMB_NPY = os.path.join(DATA_PROCESSED, "esm2_t12_representatives_embeddings.npy")
FILE_COORD_CSV = os.path.join(DATA_PROCESSED, "06_dimensionality_reduction_coordinates.csv")

OUT_METRICS_CSV = os.path.join(REPORTS_DIR, "06_cluster_quality_metrics.csv")
OUT_MEMBERSHIP_CSV = os.path.join(REPORTS_DIR, "06_cluster_membership.csv")
OUT_SENSITIVITY_CSV = os.path.join(REPORTS_DIR, "06_sensitivity_analysis.csv")

RANDOM_SEED = 42


def evaluate_clustering():
    print("=" * 80)
    print("STEP 6 PHASE 4: CLUSTERING EVALUATION AND SENSITIVITY AUDIT")
    print("=" * 80)

    # 1. Load Embeddings and Coordinates
    print(f"Loading embeddings: {FILE_EMB_NPY}")
    embeddings = np.load(FILE_EMB_NPY)
    df_coords = pd.read_csv(FILE_COORD_CSV)
    num_records, embedding_dim = embeddings.shape
    print(f"  Loaded {num_records:,} records, {embedding_dim} dimensions")

    # Standardize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(embeddings)

    # Top 50 PCA components (96.11% variance, de-noised space)
    pca = PCA(n_components=50, random_state=RANDOM_SEED)
    X_pca50 = pca.fit_transform(X_scaled)

    # Define masks
    mask_hq = df_coords["is_high_quality"].values  # Cohort B
    mask_annotated = (df_coords["gene_annotation"] != "UNANNOTATED").values  # Cohort C
    ground_truth_genes = df_coords["gene_annotation"].values

    quality_results = []
    cluster_labels_dict = {}

    # -------------------------------------------------------------
    # 2. Evaluate K-Means Across Range of k
    # -------------------------------------------------------------
    k_range = [10, 20, 30, 40, 50, 74, 100]
    print("\n--- Evaluating K-Means Clustering on 50-PC Space ---")
    for k in k_range:
        km = KMeans(n_clusters=k, random_state=RANDOM_SEED, n_init=10)
        labels = km.fit_predict(X_pca50)
        cluster_labels_dict[f"kmeans_k{k}"] = labels

        # Compute internal metrics on full representative cohort
        sil = silhouette_score(X_pca50, labels)
        dbi = davies_bouldin_score(X_pca50, labels)
        chi = calinski_harabasz_score(X_pca50, labels)

        # Concordance with ground truth on annotated subset
        ari = adjusted_rand_score(ground_truth_genes[mask_annotated], labels[mask_annotated])
        nmi = normalized_mutual_info_score(ground_truth_genes[mask_annotated], labels[mask_annotated])

        # Sensitivity on High-Quality Cohort B
        sil_hq = silhouette_score(X_pca50[mask_hq], labels[mask_hq])
        dbi_hq = davies_bouldin_score(X_pca50[mask_hq], labels[mask_hq])

        # Cluster size stats
        counts = pd.Series(labels).value_counts()
        singletons = (counts == 1).sum()

        print(f"  K-Means (k={k:3d}): Silhouette={sil:.4f} (HQ={sil_hq:.4f}), DBI={dbi:.3f}, CHI={chi:.1f}, ARI={ari:.3f}, NMI={nmi:.3f}, Singletons={singletons}")

        quality_results.append({
            "algorithm": "K-Means",
            "parameters": f"k={k}",
            "n_clusters": k,
            "noise_points": 0,
            "singleton_clusters": int(singletons),
            "silhouette_full": round(sil, 4),
            "silhouette_hq": round(sil_hq, 4),
            "davies_bouldin_full": round(dbi, 4),
            "davies_bouldin_hq": round(dbi_hq, 4),
            "calinski_harabasz": round(chi, 2),
            "ari_annotated": round(ari, 4),
            "nmi_annotated": round(nmi, 4),
            "min_cluster_size": int(counts.min()),
            "median_cluster_size": float(counts.median()),
            "max_cluster_size": int(counts.max())
        })

    # -------------------------------------------------------------
    # 3. Evaluate Agglomerative Hierarchical Clustering
    # -------------------------------------------------------------
    print("\n--- Evaluating Agglomerative Hierarchical Clustering (Ward Linkage) ---")
    for k in [20, 50, 74, 100]:
        agg = AgglomerativeClustering(n_clusters=k, linkage="ward")
        labels = agg.fit_predict(X_pca50)
        cluster_labels_dict[f"agglom_ward_k{k}"] = labels

        sil = silhouette_score(X_pca50, labels)
        dbi = davies_bouldin_score(X_pca50, labels)
        chi = calinski_harabasz_score(X_pca50, labels)
        ari = adjusted_rand_score(ground_truth_genes[mask_annotated], labels[mask_annotated])
        nmi = normalized_mutual_info_score(ground_truth_genes[mask_annotated], labels[mask_annotated])

        sil_hq = silhouette_score(X_pca50[mask_hq], labels[mask_hq])
        dbi_hq = davies_bouldin_score(X_pca50[mask_hq], labels[mask_hq])
        counts = pd.Series(labels).value_counts()
        singletons = (counts == 1).sum()

        print(f"  Agglom Ward (k={k:3d}): Silhouette={sil:.4f} (HQ={sil_hq:.4f}), DBI={dbi:.3f}, CHI={chi:.1f}, ARI={ari:.3f}, NMI={nmi:.3f}, Singletons={singletons}")

        quality_results.append({
            "algorithm": "Agglomerative (Ward)",
            "parameters": f"k={k}",
            "n_clusters": k,
            "noise_points": 0,
            "singleton_clusters": int(singletons),
            "silhouette_full": round(sil, 4),
            "silhouette_hq": round(sil_hq, 4),
            "davies_bouldin_full": round(dbi, 4),
            "davies_bouldin_hq": round(dbi_hq, 4),
            "calinski_harabasz": round(chi, 2),
            "ari_annotated": round(ari, 4),
            "nmi_annotated": round(nmi, 4),
            "min_cluster_size": int(counts.min()),
            "median_cluster_size": float(counts.median()),
            "max_cluster_size": int(counts.max())
        })

    # -------------------------------------------------------------
    # 4. Evaluate HDBSCAN Density-Based Clustering
    # -------------------------------------------------------------
    print("\n--- Evaluating HDBSCAN Density-Based Clustering ---")
    hdbscan_params = [
        {"min_cluster_size": 15, "min_samples": 5},
        {"min_cluster_size": 25, "min_samples": 5},
        {"min_cluster_size": 50, "min_samples": 5}
    ]
    for p in hdbscan_params:
        hdb = HDBSCAN(min_cluster_size=p["min_cluster_size"], min_samples=p["min_samples"], metric="euclidean")
        labels = hdb.fit_predict(X_pca50)
        label_key = f"hdbscan_mcs{p['min_cluster_size']}"
        cluster_labels_dict[label_key] = labels

        non_noise = labels != -1
        n_clusters = len(set(labels[non_noise]))
        n_noise = int((labels == -1).sum())

        if n_clusters > 1 and non_noise.sum() > n_clusters:
            sil = silhouette_score(X_pca50[non_noise], labels[non_noise])
            dbi = davies_bouldin_score(X_pca50[non_noise], labels[non_noise])
            chi = calinski_harabasz_score(X_pca50[non_noise], labels[non_noise])
            
            mask_ann_non_noise = mask_annotated & non_noise
            ari = adjusted_rand_score(ground_truth_genes[mask_ann_non_noise], labels[mask_ann_non_noise])
            nmi = normalized_mutual_info_score(ground_truth_genes[mask_ann_non_noise], labels[mask_ann_non_noise])
        else:
            sil, dbi, chi, ari, nmi = np.nan, np.nan, np.nan, np.nan, np.nan

        counts = pd.Series(labels[non_noise]).value_counts() if non_noise.any() else pd.Series([0])
        singletons = (counts == 1).sum()

        print(f"  HDBSCAN (min_size={p['min_cluster_size']}): Clusters={n_clusters}, Noise={n_noise} ({n_noise/len(labels)*100:.1f}%), Silhouette={sil:.4f}, ARI={ari:.3f}, NMI={nmi:.3f}")

        quality_results.append({
            "algorithm": "HDBSCAN",
            "parameters": f"min_cluster_size={p['min_cluster_size']},min_samples={p['min_samples']}",
            "n_clusters": n_clusters,
            "noise_points": n_noise,
            "singleton_clusters": int(singletons),
            "silhouette_full": round(sil, 4) if not np.isnan(sil) else np.nan,
            "silhouette_hq": np.nan,
            "davies_bouldin_full": round(dbi, 4) if not np.isnan(dbi) else np.nan,
            "davies_bouldin_hq": np.nan,
            "calinski_harabasz": round(chi, 2) if not np.isnan(chi) else np.nan,
            "ari_annotated": round(ari, 4) if not np.isnan(ari) else np.nan,
            "nmi_annotated": round(nmi, 4) if not np.isnan(nmi) else np.nan,
            "min_cluster_size": int(counts.min()) if non_noise.any() else 0,
            "median_cluster_size": float(counts.median()) if non_noise.any() else 0,
            "max_cluster_size": int(counts.max()) if non_noise.any() else 0
        })

    # Save Quality Metrics Table
    df_metrics = pd.DataFrame(quality_results)
    df_metrics.to_csv(OUT_METRICS_CSV, index=False)
    print(f"\nSaved cluster quality metrics to: {OUT_METRICS_CSV}")

    # -------------------------------------------------------------
    # 5. Selection of Primary Model for Downstream Characterization
    # -------------------------------------------------------------
    # K-Means at k=74 represents the biologically grounded core model,
    # matching the 74 canonical HSV gene families with high NMI (0.763) and 0 noise.
    # HDBSCAN (mcs=25) provides the complementary density/noise model.
    primary_model_key = "kmeans_k74"
    primary_labels = cluster_labels_dict[primary_model_key]
    hdb_labels = cluster_labels_dict["hdbscan_mcs25"]

    df_membership = df_coords.copy()
    df_membership["primary_cluster_id"] = primary_labels
    df_membership["hdbscan_cluster_id"] = hdb_labels
    df_membership.to_csv(OUT_MEMBERSHIP_CSV, index=False)
    print(f"Saved cluster membership table to: {OUT_MEMBERSHIP_CSV}")

    # -------------------------------------------------------------
    # 6. Generate Figure 5: UMAP Colored by Primary Cluster Assignments
    # -------------------------------------------------------------
    plt.figure(figsize=(10, 8))
    # Pick categorical colormap with enough distinct colors
    cmap = plt.get_cmap("tab20", 74)
    scatter = plt.scatter(
        df_membership["umap1_nn30_md01_primary"],
        df_membership["umap2_nn30_md01_primary"],
        c=primary_labels,
        cmap="turbo",
        alpha=0.7,
        s=25,
        edgecolors="none"
    )
    plt.colorbar(scatter, label="Primary Cluster ID (K-Means, k=74)")
    plt.title("UMAP Projection Colored by Primary Cluster Assignments\n(K-Means k=74 on 50-PC ESM-2 space, N=3,111)")
    plt.xlabel("UMAP Dimension 1")
    plt.ylabel("UMAP Dimension 2")
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    fig5_path = os.path.join(FIGURES_DIR, "05_umap_cluster_assignment.png")
    plt.savefig(fig5_path, dpi=300)
    plt.close()
    print(f"  Saved: {fig5_path}")

    # -------------------------------------------------------------
    # 7. Sensitivity Comparison Table (Cohort A vs Cohort B)
    # -------------------------------------------------------------
    sensitivity_records = []
    for k in [20, 50, 74]:
        labels_full = cluster_labels_dict[f"kmeans_k{k}"]
        # Retrain K-Means strictly on Cohort B
        km_b = KMeans(n_clusters=k, random_state=RANDOM_SEED, n_init=10)
        labels_b = km_b.fit_predict(X_pca50[mask_hq])

        sil_a = silhouette_score(X_pca50, labels_full)
        sil_b = silhouette_score(X_pca50[mask_hq], labels_b)
        dbi_a = davies_bouldin_score(X_pca50, labels_full)
        dbi_b = davies_bouldin_score(X_pca50[mask_hq], labels_b)
        
        # Ground truth agreement on annotated subsets
        mask_ann_a = mask_annotated
        mask_ann_b = mask_annotated[mask_hq]
        nmi_a = normalized_mutual_info_score(ground_truth_genes[mask_ann_a], labels_full[mask_ann_a])
        nmi_b = normalized_mutual_info_score(ground_truth_genes[mask_hq][mask_ann_b], labels_b[mask_ann_b])

        sensitivity_records.append({
            "model": f"K-Means (k={k})",
            "cohort_a_records": len(X_pca50),
            "cohort_a_silhouette": round(sil_a, 4),
            "cohort_a_dbi": round(dbi_a, 4),
            "cohort_a_nmi_annotated": round(nmi_a, 4),
            "cohort_b_records": int(mask_hq.sum()),
            "cohort_b_silhouette": round(sil_b, 4),
            "cohort_b_dbi": round(dbi_b, 4),
            "cohort_b_nmi_annotated": round(nmi_b, 4),
            "silhouette_delta": round(sil_b - sil_a, 4),
            "nmi_delta": round(nmi_b - nmi_a, 4)
        })

    df_sens = pd.DataFrame(sensitivity_records)
    df_sens.to_csv(OUT_SENSITIVITY_CSV, index=False)
    print(f"Saved sensitivity analysis table to: {OUT_SENSITIVITY_CSV}")

    print("\n================================================================================")
    print("STEP 6 PHASE 4 COMPLETE: CLUSTERING ANALYSIS SUCCEEDED.")
    print("================================================================================")


if __name__ == "__main__":
    evaluate_clustering()
