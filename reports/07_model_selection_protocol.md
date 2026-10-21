# Phase C & D Audit Report: Label-Free Model Selection and Stability Audit

**Audit Date**: 2026-10-21 13:14:20 UTC  
**Auditor**: Independent Computational Biologist  
**Target Manifest**: `reports/07_unsupervised_model_comparison.csv`  

---

## 1. Executive Summary & Anti-Circularity Protocol

To prevent confirmation bias, Step 7 evaluated a pre-specified grid of **67 distinct clustering experiments** across 5 input representations and 3 algorithm families **WITHOUT USING ANY GENE LABELS OR SPECIES LABELS**.

### Core Findings on Model Selection:
1. **Internal Geometric Trajectory**:
   - Silhouette scores increase monotonically from $k=10$ ($0.1837$) to $k=150$ ($0.4124$) as clusters become smaller and more compact.
   - Davies-Bouldin index shows optimal separation (lowest values) in the $k=74 \dots 100$ range ($1.357 \dots 1.223$).
   - K-Means with $k=74$ is scientifically defensible based purely on internal stability (seed ARI = $0.941 \pm 0.015$, 80% subsample ARI = $0.812 \pm 0.024$).
2. **HDBSCAN Density-Based Natural Partitioning**:
   - Across `min_cluster_size` in $[10, 15, 20, 25, 35, 50]$, HDBSCAN with `min_cluster_size=15` natively detects **74 dense clusters** with a silhouette score of **0.4747** and 697 noise points (22.4%).
   - At `min_cluster_size >= 50`, the dataset collapses into 2 macro-clusters, confirming that viral proteomes are modular and multi-family rather than binary.
3. **Representation Hierarchy**:
   - PCA-reduced embeddings (top 50 PCs, 96.1% variance) outperform raw 480d embeddings in cluster compactness and computational efficiency by eliminating high-frequency noise.
   - Clustering purely on 10-dimensional physicochemical features yields inferior separation (Silhouette $\approx 0.16$), demonstrating that ESM-2 language model representations capture rich structural homology beyond gross physical descriptors.

---

## 2. Stability Audit Summary Table

| Perturbation Type | Model Evaluated | Metric | Mean Stability | Verdict |
| :--- | :--- | :--- | :---: | :--- |
| **5 Random Seeds** | K-Means ($k=74$) | Pairwise ARI | **0.941** | High deterministic stability |
| **80% Subsampling (10 iters)** | K-Means ($k=74$) | Overlap ARI | **0.812** | Robust to sequence dropout |
| **Scaling (Standard vs Robust)** | K-Means ($k=74$) | Cross-Method ARI | **0.985** | Insensitive to outlier scaling |
| **Scaling (Standard vs L2)** | K-Means ($k=74$) | Cross-Method ARI | **0.892** | Cosine and Euclidean alignment |
| **PCA Cutoff (20 vs 50 PCs)** | K-Means ($k=74$) | Cross-Dimension ARI| **0.865** | Dominant geometry in top 20 PCs |
| **PCA Cutoff (50 vs 100 PCs)** | K-Means ($k=74$) | Cross-Dimension ARI| **0.948** | Trailing PCs contribute negligible signal |
| **Cohort B Sensitivity (<=5% ambig)**| K-Means ($k=74$) | Shared Record ARI | **0.884** | Gap artifacts do not drive clustering |
