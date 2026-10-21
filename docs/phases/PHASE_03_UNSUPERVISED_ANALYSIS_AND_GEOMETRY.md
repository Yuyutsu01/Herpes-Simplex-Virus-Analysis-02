# Phase 3: Unsupervised Representation Analysis, Geometry & Cluster Validation

**Consolidated Steps**: Step 6 (Dimensionality Reduction & Clustering), Step 7 (Unsupervised Validation & Confounding Audit)  
**Primary Deliverables**: `data/processed/06_dimensionality_reduction_coordinates.csv`, `reports/07_unsupervised_model_comparison.csv`  

---

## 1. Scientific Objective
Evaluate whether the latent geometry of ESM-2 embeddings reveals intrinsic biological structure (e.g. gene families, functional classes) without supervision, and determine whether unsupervised clustering discovers novel functional categories or simply reconstructs known sequence homology.

## 2. Inputs
- `data/processed/esm2_t12_representatives_embeddings.npy` (3,111 sequences $	imes$ 480 dimensions)
- Canonical gene family annotations and homology groups

## 3. Methodology
1. **Dimensionality Reduction**:
   - PCA (Principal Component Analysis): Computed top 50 components.
   - UMAP (Uniform Manifold Approximation and Projection): Parameters $n\_neighbors \in \{15, 30, 50\}$, $min\_dist \in \{0.1, 0.3\}$.
   - t-SNE: Perplexity = 30, iterations = 1000.
2. **Clustering**:
   - HDBSCAN density clustering across varying minimum cluster sizes ($min\_cluster\_size \in \{5, 10, 15\}$).
   - K-Means clustering ($K \in [5, 74]$).
3. **Ortholog Co-Clustering & Confounding Audits**:
   - Evaluated co-clustering of 74 canonical HSV-1/HSV-2 ortholog pairs.
   - Evaluated sequence length and GC content as potential confounders.

## 4. Key Findings & Epistemological Boundaries
- **Homology Re-Identification**: Clusters predominantly re-identified known homologous sequence families (Adjusted Rand Index $> 0.85$ vs gene symbol). Unsupervised clustering does **not** discover de novo viral functional classes.
- **Ortholog Alignment**: Conserved core replication enzymes (e.g., DNA polymerase UL30, helicase-primase UL5/UL8/UL52) cluster tightly by sister-species orthology.

## 5. Validation & Test Outcomes
- Validated under `tests/test_06_unsupervised_analysis.py` and `tests/test_07_independent_validation.py` (20 passing tests).
