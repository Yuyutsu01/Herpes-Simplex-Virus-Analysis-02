# Phase G Audit Report: Confounding Audit, Selection Bias, and Noise Anatomy

**Audit Date**: 2026-10-21 13:18:34 UTC  
**Target Corpus**: Frozen Representative Cohort (N = 3,111 sequences)  
**Deliverables**: `reports/07_confounding_audit.csv`, `reports/07_quality_sensitivity.csv`  

---

## 1. Executive Summary & Confounder Assessment

We audited whether the unsupervised manifold and clustering structures observed in Step 6 are driven by biological homology or technical confounders:

1. **Sequence Length**:
   - Kruskal-Wallis test across 74 clusters: $H = 2649.8$, $p < 10^{-100}$, $\eta^2 = 0.848$.
   - **Interpretation**: While protein length differs significantly across clusters, this reflects genuine biological diversity (e.g. UL11 at 96 aa vs UL36 at 3,164 aa). Length alone does not explain sub-family partitioning, as proteins of identical length (e.g. UL27 and UL19) separate cleanly based on sequence semantics.
2. **Ambiguity Fraction (Sequencing Gap Confounding)**:
   - Kruskal-Wallis test: $H = 691.4$, $p < 10^{-100}$, $\eta^2 = 0.204$.
   - Ambiguity fraction is heavily localized to specific GC-rich repeat loci (`RL1` and `RS1`), where sequencing gaps (`X` runs) cause local centroid fragmentation.
3. **Correlation vs Causation Guardrail**:
   - Principal component correlations with molecular weight ($\rho = +0.439$), isoelectric point ($\rho = +0.761$), and GRAVY ($\rho = -0.585$) demonstrate that language model embeddings capture physical chemistry, but **do not prove that these properties caused cluster formation**.

---

## 2. Selection Bias in Quality Filtering (Cohort B)

Excluding sequences with $>5\%$ ambiguous residues (Cohort B) removed **1,496 records (48.1%)**. An audit of the removed records reveals severe biological selection bias:

| Gene Locus | Total Representatives | Retained in Cohort B | Removed in Cohort B | % Depletion | Biological Mechanism |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **RL1 (ICP34.5)** | **135** | **1** | **134** | **99.3%** | Extreme GC-rich repeat masking in clinical isolates |
| **RS1 (ICP4)** | **212** | **7** | **205** | **96.7%** | Diploid terminal repeat sequencing gap runs |
| **UNANNOTATED** | **1,489** | **992** | **497** | **33.4%** | Degraded / short provisional open reading frames |
| **UL27 (gB)** | **32** | **31** | **1** | **3.1%** | Core conserved essential glycoprotein (clean sequencing) |
| **UL19 (MCP)** | **23** | **22** | **1** | **4.3%** | Major capsid protein (clean sequencing) |

> [!WARNING]
> **Methodological Selection Bias Finding**: Filtering out sequences with $>5\%$ ambiguous residues disproportionately eliminates **over 95% of all RL1 and RS1 records** from the analysis! While Cohort B is chemically pure, it blindspots the diploid regulatory and neurovirulence factors of the herpesvirus genome.

---

## 3. Anatomy of HDBSCAN Noise Points (N = 697)

HDBSCAN (`min_cluster_size=15`) assigned **697 sequences (22.4%)** to noise (`cluster = -1`):
1. **Unannotated Over-Representation**: **536 of 697 noise points (76.9%)** are `UNANNOTATED` proteins. Their median length is 269.0 aa (vs 500.0 aa for clustered proteins).
2. **Biological Reality of Noise**: These noise points do not represent software failures; they capture atypical, truncated, or divergent provisional open reading frames that lack sufficient sequence density to form cohesive homologous families.
