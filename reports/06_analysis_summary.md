# Step 6 Analysis Summary: Unsupervised Learning, Dimensionality Reduction, and Protein Cluster Analysis

**Project**: *Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes*  
**Date**: 2026-10-21 12:57:26 UTC  
**Dataset**: Frozen Representative Cohort (N = 3,111 sequences)  
**Embeddings**: ESM-2 `facebook/esm2_t12_35M_UR50D` (480-dimensional latent representations)  
**Manifest**: `reports/06_analysis_manifest.json`  

---

## 1. Executive Summary & Core Scientific Findings

Step 6 investigated the topological and structural organization of the HSV-1 and HSV-2 proteomes by integrating ESM-2 protein language model embeddings with sequence-level physicochemical features, canonical gene annotations, and viral species labels.

### Key Conclusions:

1. **Homology Dominates Over Viral Species**:
   - The ESM-2 embedding space is overwhelmingly structured by **gene family and structural function** (Normalized Mutual Information $\text{NMI} = 0.867$, Adjusted Rand Index $\text{ARI} = 0.658$), rather than viral species ($\text{NMI} = 0.084$, $\text{ARI} = 0.020$).
   - **82.4% of canonical ortholog gene families (61 of 74)** co-cluster into the **exact same** primary cluster between HSV-1 and HSV-2.
   - Homologous proteins across species are vastly closer in embedding space than non-homologous proteins within the same viral species.

2. **PCA Latent Axes Align with Physical Chemistry**:
   - Although ESM-2 was trained without explicit biophysical objectives, its principal linear components correlate strongly and significantly with physical chemistry:
     - **PC1 (22.14% variance)**: Correlates positively with isoelectric point ($\rho = +0.761$, $p < 10^{-300}$) and instability index ($\rho = +0.598$), and negatively with aliphatic index ($\rho = -0.734$) and GRAVY ($\rho = -0.585$).
     - **PC2 (16.14% variance)**: Correlates negatively with polar residue fraction ($\rho = -0.726$, $p < 10^{-300}$) and positively with sequence length/molecular weight ($\rho = +0.474$).

3. **HDBSCAN Native Cluster Discovery Matches Canonical Biology**:
   - Without supervision, HDBSCAN (`min_cluster_size=15`) natively discovered **exactly 74 dense clusters**, matching the 74 canonical core gene families with an NMI of **0.967** and ARI of **0.821**.
   - HDBSCAN categorized 697 divergent/provisional proteins (22.4%) as noise (`cluster = -1`), effectively isolating unannotated or degraded candidates from the canonical viral core.

4. **Robustness to Sequence Quality**:
   - Comparing the full representative cohort (Cohort A: 3,111 records) against the high-quality sensitivity cohort (Cohort B: 1,615 records with $\le 5\%$ ambiguous residues) demonstrated that observed cluster topologies are stable and not driven by sequencing-gap artifacts ($\Delta \text{NMI} < 0.03$).

---

## 2. Analysis Cohorts Definition

| Cohort | Description | Filter Criteria | Record Count | % Cohort | Analytical Role |
| :--- | :--- | :--- | :---: | :---: | :--- |
| **Cohort A** | Main Representative Cohort | `final_disposition == 'RETAINED_UNIQUE_CANDIDATE'` | **3,111** | 100.0% | Primary unsupervised manifold learning and clustering |
| **Cohort B** | High-Quality Sensitivity Cohort | `ambiguous_residue_fraction <= 0.05` | **1,615** | 51.9% | Sensitivity audit excluding clinical sequencing gap masking |
| **Cohort C** | Annotated-Protein Core | `gene_annotation != 'UNANNOTATED'` | **1,622** | 52.1% | External validation against canonical HSV gene families |
| **Cohort D** | Provisional Unannotated Set | `gene_annotation == 'UNANNOTATED'` | **1,489** | 47.9% | Uncharacterized ORF cluster assignment and discovery |

---

## 3. Dimensionality Reduction Metrics

### A. Principal Component Analysis (PCA)
- Standardized 480-dimensional ESM-2 embedding space (zero mean, unit variance).
- **PC1 Explained Variance**: 22.14%
- **PC2 Explained Variance**: 16.14%
- **Top 10 PCs Cumulative Variance**: 72.58%
- **Top 50 PCs Cumulative Variance**: 96.11%
- *Figure Artifacts*:
  - `reports/figures/01_pca_explained_variance.png`: Scree and cumulative variance plots.
  - `reports/figures/02_pca_species_projection.png`: PC1 vs PC2 colored by species (intermixed).
  - `reports/figures/03_pca_gene_family_projection.png`: PC1 vs PC2 colored by functional gene category.

### B. Uniform Manifold Approximation and Projection (UMAP)
- Computed on standardized 50-PC space using cosine distance metric.
- Evaluated across grid: `(nn=15, md=0.1)`, `(nn=30, md=0.1)`, `(nn=50, md=0.3)`.
- *Figure Artifacts*:
  - `reports/figures/04_umap_species_projection.png`: Broad intermixing of HSV-1 and HSV-2.
  - `reports/figures/05_umap_cluster_assignment.png`: 74 primary clusters projected on UMAP.
  - `reports/figures/06_umap_quality_projection.png`: Ambiguity fraction gradient.
  - `reports/figures/08_sensitivity_comparison_cohort_a_vs_b.png`: Manifold stability across Cohorts A and B.

---

## 4. Multi-Algorithm Clustering Comparison

Evaluated on 50-PC standardized embedding space (96.11% variance retained):

| Algorithm | Parameters | Clusters | Noise Points | Silhouette (Full) | Silhouette (HQ) | Davies-Bouldin | Calinski-Harabasz | ARI (Genes) | NMI (Genes) |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **K-Means** | $k=20$ | 20 | 0 | 0.1951 | 0.1591 | 1.613 | 314.8 | 0.332 | 0.658 |
| **K-Means** | $k=50$ | 50 | 0 | 0.2812 | 0.2421 | 1.458 | 218.7 | 0.653 | 0.850 |
| **K-Means (Primary)** | **$k=74$** | **74** | **0** | **0.3220** | **0.2824** | **1.357** | **195.6** | **0.658** | **0.867** |
| **K-Means** | $k=100$ | 100 | 0 | 0.3638 | 0.3305 | 1.223 | 184.8 | 0.716 | 0.910 |
| **Agglomerative** | Ward ($k=74$) | 74 | 0 | 0.3296 | 0.2947 | 1.398 | 192.7 | 0.715 | 0.889 |
| **HDBSCAN** | mcs=15, ms=5 | **74** | 697 (22.4%) | **0.4747** | — | — | — | **0.821** | **0.967** |
| **HDBSCAN** | mcs=25, ms=5 | 50 | 1,073 (34.5%)| 0.4716 | — | — | — | 0.811 | 0.960 |

---

## 5. Physicochemical Associations (Bonferroni-Corrected)

Spearman rank correlation ($m = 50$ tests, significance threshold $\alpha = 0.001$):
- **PC1 (22.14% var)**:
  - Strongly positive with Isoelectric Point ($\rho = +0.761$) and Instability Index ($\rho = +0.598$).
  - Strongly negative with Aliphatic Index ($\rho = -0.734$), GRAVY ($\rho = -0.585$), and Aromatic Fraction ($\rho = -0.570$).
- **PC2 (16.14% var)**:
  - Strongly negative with Polar Fraction ($\rho = -0.726$).
  - Strongly positive with Sequence Length ($\rho = +0.474$) and Molecular Weight ($\rho = +0.439$).
- *Figure Artifact*:
  - `reports/figures/07_cluster_physicochemical_distributions.png`: Boxplots of pI, GRAVY, and length across functional classes.

---

## 6. HSV-1 versus HSV-2 Species Comparison

- **Species Distribution in Representative Set**: HSV-1 = 1,124 (36.13%), HSV-2 = 1,987 (63.87%).
- **Cluster Species Skew**: 46 of 74 clusters (62.2%) show no statistically significant deviation from baseline species proportions ($p_{\text{FDR}} \ge 0.05$).
- **Ortholog Co-Clustering**: 61 of 74 dual-species canonical gene families (82.4%) share the exact same primary cluster across HSV-1 and HSV-2.
- **Species vs Gene Association**:
  - Species NMI: **0.0843** (negligible species-specific partitioning).
  - Gene Family NMI: **0.8670** (dominant functional partitioning).

---

## 7. Leakage Prevention and Data Traceability

- All 3,111 records analyzed belong strictly to `RETAINED_UNIQUE_CANDIDATE` ($<95\%$ identity across clusters).
- The `near_identity_cluster_id_95` and `prevent_ml_leakage_group` fields are mapped for all records.
- In future supervised tasks, cross-validation folds MUST split on `prevent_ml_leakage_group` to avoid homologous data leakage.

---

## 8. Limitations & Guardrails

1. **No Temporal Class Over-Interpretation**: Cluster assignments must **not** be equated with Immediate-Early ($\alpha$), Early ($\beta$), or Late ($\gamma$) kinetic classes without independent experimental labeling.
2. **Provisional Unannotated Proteins**: 1,489 records lack canonical gene symbols. In HDBSCAN, 697 records were flagged as noise (`-1`), suggesting that many unannotated ORFs are truncated, atypical, or divergent.
3. **No Causality**: Physicochemical correlations with PCs demonstrate alignment between language model geometry and physical chemistry, but do not imply a causal mechanism.
