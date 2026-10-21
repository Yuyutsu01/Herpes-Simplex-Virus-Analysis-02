"""
Script: reorganize_docs_phases.py
Purpose: Consolidate project documentation from Steps 1-16 into 10 comprehensive,
  authoritative Markdown documents in docs/phases/.
"""

from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
PHASES_DIR = BASE_DIR / "docs" / "phases"
PHASES_DIR.mkdir(parents=True, exist_ok=True)


def create_phase_documents():
    print("=== Consolidating Documentation into docs/phases/ ===")

    # -------------------------------------------------------------
    # PHASE 1: Data Inventory, Taxonomic Curation & Redundancy Reduction
    # -------------------------------------------------------------
    p1_content = """# Phase 1: Data Inventory, Taxonomic Curation & Redundancy Reduction

**Consolidated Steps**: Step 1 (FASTA Inventory), Step 2 (Taxonomy Audit), Step 3 (Target Proteome Curation), Step 4 (Homology & Redundancy Reduction)  
**Primary Deliverables**: `reports/04_redundancy_decisions.csv`, `reports/04_final_validation.md`, `data/processed/curated_protein_candidates.fasta`  

---

## 1. Scientific Objective
Establish a clean, non-redundant, authoritative viral protein sequence corpus from raw public databases (NCBI GenBank, RefSeq, UniProtKB) for Herpes Simplex Virus Type 1 (HSV-1) and Type 2 (HSV-2). The primary challenge in viral computational proteomics is that public sequence repositories contain massive clinical strain redundancy, truncated sequencing fragments, host interactors, and recombinant vector contaminants (such as HIV-1 Gag/Pol).

## 2. Inputs & Source Data
- `data/raw/hsv1.fasta` (SHA-256: `15001269a835b41223467e0d91bdcada433ee68c6a58c19e684c060483b31f36`, Size: 7,954,295 bytes, 12,028 records)
- `data/raw/hsv2.fasta` (SHA-256: `86c9d627286aa4cad121af4c84cdb39b824c14905aa57948de9a996955ca0311`, Size: 7,095,643 bytes, 10,682 records)
- Universal total: **22,710 records**

## 3. Methodology & Transformations
1. **Initial Sequence Parsing & Exact Hashing**: Each FASTA entry was digested with SHA-256 to identify identical amino acid sequences across accessions.
2. **Taxonomic & Host Contaminant Audit**:
   - Identified 1,043 non-target records including HIV-1 Gag/Pol/Env structural proteins from lentiviral gene therapy vectors and host ribosomal interactors.
   - Quarantined 254 unverified entries (synthetic patent claims, PDB constructs).
   - Quarantined 30 records (15 duplicate pairs) exhibiting identical sequences with contradictory species labels (`[Human alphaherpesvirus 1]` vs `[Human alphaherpesvirus 2]`).
3. **Length Thresholding**: Sequestered 3,211 truncated fragments (< 70 amino acids) into a separate fragment cohort, leaving 18,172 valid full-length candidates (HSV-1: 9,454; HSV-2: 8,718).
4. **95% Identity Redundancy Reduction (CD-HIT)**:
   - Partitioned the 18,172 candidates into **3,111 representative cluster seeds** (13.70% of corpus) and **15,061 near-identical strain variants** (66.32% of corpus).
   - Mapped all 3,111 representative seeds to canonical gene families (74 ortholog pairs).

## 4. Key Outputs & Provenance
- `reports/04_redundancy_decisions.csv`: Exhaustive 22,710-row ledger documenting the disposition and cluster assignment of every sequence.
- `data/processed/valid_proteins_combined.fasta`: 18,172 valid full-length viral candidates.
- `data/processed/esm2_t12_representatives_metadata.csv`: Metadata for the 3,111 representative cluster seeds.

## 5. Validation & Test Outcomes
- Validated under `tests/test_04_duplicate_homology.py` (10 passing tests).
- Mathematical closure confirmed: 18,172 (valid) + 3,211 (fragments) + 1,043 (purged) + 254 (unverified) + 30 (conflicts) = exactly 22,710 records.
"""
    (PHASES_DIR / "PHASE_01_DATA_INVENTORY_AND_CURATION.md").write_text(p1_content, encoding="utf-8")

    # -------------------------------------------------------------
    # PHASE 2: Feature Engineering & ESM-2 Embeddings
    # -------------------------------------------------------------
    p2_content = """# Phase 2: Feature Engineering & Protein Language Model Embeddings

**Consolidated Step**: Step 5 (Frozen Dataset & ESM-2 Embedding Generation)  
**Primary Deliverables**: `data/processed/esm2_t12_representatives_embeddings.npy`, `data/processed/esm2_t12_representatives_metadata.csv`  

---

## 1. Scientific Objective
Extract multi-scale representations for all 3,111 representative viral proteins: classical bulk physicochemical descriptors and deep contextual transformer embeddings from the Evolutionary Scale Modeling (ESM-2) protein language model.

## 2. Inputs
- `data/processed/valid_proteins_combined.fasta` (Filtered representative seeds)
- Pre-trained model weights: `esm2_t12_35M_UR50D` (12 layers, 480 embedding dimensions, 35M parameters)

## 3. Methodology & Feature Extraction
1. **Bulk Physicochemical Descriptors ($D=8$)**:
   - Sequence length (aa), molecular weight (Da), isoelectric point (pI), aliphatic index, instability index, GRAVY hydropathicity, charged residue fraction, aromatic residue fraction.
2. **ESM-2 Embeddings ($D=480$)**:
   - Computed on GPU/CPU using fairseq / ESM PyTorch library.
   - Extracted residue representations from layer 12 and applied mean-pooling across the length of each protein sequence.
   - Output stored as a dense floating-point matrix of shape `(3111, 480)`.
3. **Combined Representation ($D=488$)**:
   - Concatenation of the 8 physical descriptors with the 480 ESM-2 dimensions.

## 4. Key Outputs & Provenance
- `data/processed/esm2_t12_representatives_embeddings.npy`: Binary NumPy array `(3111, 480)` (5,973,248 bytes).
- `data/processed/esm2_t12_representatives_metadata.csv`: Complete metadata joining FASTA headers, accessions, lengths, and embedding row indices.

## 5. Validation & Test Outcomes
- Validated under `tests/test_05_frozen_dataset_and_embeddings.py` (10 passing tests).
- Row-index closure confirmed: `row_index` in metadata strictly maps 1-to-1 to row index of `.npy` array without off-by-one errors or NaN entries.
"""
    (PHASES_DIR / "PHASE_02_FEATURE_ENGINEERING_AND_EMBEDDINGS.md").write_text(p2_content, encoding="utf-8")

    # -------------------------------------------------------------
    # PHASE 3: Unsupervised Analysis, Geometry & Cluster Validation
    # -------------------------------------------------------------
    p3_content = """# Phase 3: Unsupervised Representation Analysis, Geometry & Cluster Validation

**Consolidated Steps**: Step 6 (Dimensionality Reduction & Clustering), Step 7 (Unsupervised Validation & Confounding Audit)  
**Primary Deliverables**: `data/processed/06_dimensionality_reduction_coordinates.csv`, `reports/07_unsupervised_model_comparison.csv`  

---

## 1. Scientific Objective
Evaluate whether the latent geometry of ESM-2 embeddings reveals intrinsic biological structure (e.g. gene families, functional classes) without supervision, and determine whether unsupervised clustering discovers novel functional categories or simply reconstructs known sequence homology.

## 2. Inputs
- `data/processed/esm2_t12_representatives_embeddings.npy` (3,111 sequences $\times$ 480 dimensions)
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
"""
    (PHASES_DIR / "PHASE_03_UNSUPERVISED_ANALYSIS_AND_GEOMETRY.md").write_text(p3_content, encoding="utf-8")

    # -------------------------------------------------------------
    # PHASE 4: Temporal Evidence Curation & Benchmark Split Design
    # -------------------------------------------------------------
    p4_content = """# Phase 4: Temporal Evidence Curation & Benchmark Split Design

**Consolidated Steps**: Step 8 & 8.1 (Temporal Evidence Curation & Integrity Audit), Step 9 & 9.1 (Supervised Benchmark Design & Feasibility)  
**Primary Deliverables**: `data/processed/08_1_strict_temporal_benchmark.csv`, `data/processed/09_group_split_manifest.csv`  

---

## 1. Scientific Objective
Establish an authoritative, leakage-resistant benchmark for viral replication timing by curating primary biochemical literature, recovering submitter gene aliases, adjudicating experimental evidence tiers, and designing a group-disjoint cross-validation split that prevents homology leakage.

## 2. Primary Literature Sources & Evidence Tiers
- **Honess & Roizman (1974, 1975)**: Fundamental metabolic arrest assays using cycloheximide (CHX; protein synthesis inhibitor) and phosphonoacetic acid (PAA; viral DNA polymerase inhibitor).
- **Rutkowski et al. (2015)**: Time-resolved 4sU nascent RNA sequencing.
- **Kulej et al. (2017)**: SILAC quantitative mass spectrometry.
- **Evidence Tiers**:
  1. *Direct Protein-Validated* ($N=52$ families): Synthesis timing verified at the protein level under metabolic arrest.
  2. *Transcript-Only Annotated* ($N=11$ families): Kinetics inferred exclusively from 4sU RNA-seq.

## 3. Alias Recovery & Benchmark Partitioning
- A 10-gene regex dictionary resolved 270 submitter aliases (e.g., "alpha-trans-inducing factor" $\to$ `UL48`, "major capsid protein" $\to$ `UL19`).
- Partition of 3,111 representative cluster seeds:
  - **Strict Direct Evidence**: 687 sequences across 63 independent gene families.
  - **Homology Inferred**: 633 sequences (HSV-2 sequences assigned labels solely via HSV-1 orthology; excluded from direct benchmark).
  - **Unresolved / Accessory ORFs**: 1,791 sequences.

## 4. Group-Disjoint 5-Fold Cross-Validation Design
- All sequences belonging to the same gene family were allocated to exactly one validation fold ($\text{Overlap} = 0$).
- **RS1 (ICP4) Family Isolation**: `RS1` contains 212 sequence variants (30.9% of the benchmark). It was allocated strictly to Fold 0 validation. Zero `RS1` sequences exist in Fold 0 training.

## 5. Validation & Test Outcomes
- Validated under `tests/test_08_temporal_annotation.py`, `test_08_1_temporal_integrity.py`, `test_09_benchmark_design.py`, `test_09_1_split_validity.py` (40 passing tests).
"""
    (PHASES_DIR / "PHASE_04_TEMPORAL_EVIDENCE_AND_BENCHMARK_DESIGN.md").write_text(p4_content, encoding="utf-8")

    # -------------------------------------------------------------
    # PHASE 5: Predictive Modeling & Baseline Evaluation
    # -------------------------------------------------------------
    p5_content = """# Phase 5: Predictive Modeling & Baseline Evaluation

**Consolidated Steps**: Step 10 (Baseline Modeling), Step 10.1 (Independent Metrics Audit)  
**Primary Deliverables**: `reports/13_retrained_oof_predictions.csv`, `reports/10_baseline_metrics.csv`, `reports/10_family_level_metrics.csv`  

---

## 1. Scientific Objective
Evaluate whether protein representations (physicochemical properties vs. ESM-2 transformer embeddings) contain discriminative signal distinguishing pre-replicative regulatory enzymes ($\\text{PRE} = \\text{IE} + \\text{E}$) from post-replicative structural virion proteins ($\\text{POST} = \\text{L1} + \\text{L2}$) across unseen viral gene families.

## 2. Models Evaluated
1. **Majority Baseline**: `DummyClassifier(strategy='most_frequent')`
2. **Physicochemical Logistic Regression**: `LogisticRegression(C=1.0, penalty='l2', solver='lbfgs', max_iter=1000)` on 8 physical features.
3. **Physicochemical Random Forest**: `RandomForestClassifier(n_estimators=100, max_depth=5, min_samples_leaf=3)` on 8 physical features.
4. **Ridge-Regularized ESM-2**: `LogisticRegression(C=1.0, penalty='l2', solver='lbfgs', max_iter=1000)` on 480 ESM-2 dimensions.
5. **Combined Model**: `LogisticRegression(C=1.0, penalty='l2', solver='lbfgs', max_iter=1000)` on 488 concatenated features.

## 3. Evaluation Protocol
- **Generalization Unit**: Gene Family ($N=63$). Family out-of-fold probability is the mean sequence probability: $\\bar{P}_{\\text{PRE}} \\ge 0.5 \\implies \\text{PRE}$.
- Preprocessing: `StandardScaler` fit strictly inside training folds.

## 4. Re-Verified Benchmark Results
| Model | Family Accuracy | Family BalAcc | Family Macro-F1 | Sequence Macro-F1 | Correct Families |
| :--- | :---: | :---: | :---: | :---: | :---: |
| Majority Baseline | 47.62% | 53.21% | 0.4741 | 0.2651 | 30 / 63 |
| Physicochemical RF | 58.73% | 51.79% | 0.5102 | 0.5399 | 37 / 63 |
| Physicochemical LogReg | 60.32% | 56.74% | 0.5679 | 0.5747 | 38 / 63 |
| Combined (Phys + ESM-2) | 77.78% | 76.96% | 0.7644 | 0.8966 | 49 / 63 |
| **Ridge ESM-2** | **80.95%** | **79.46%** | **0.7946** | **0.9112** | **51 / 63** |

## 5. Validation & Test Outcomes
- Validated under `tests/test_10_baseline_evaluation.py` and `test_10_1_independent_verification.py` (20 passing tests).
- Mathematical equality verified: For $N=63$ with 17 true PRE and 34 true POST correct, Balanced Accuracy $\\equiv$ Macro-F1 $\\equiv 0.7946$.
"""
    (PHASES_DIR / "PHASE_05_PREDICTIVE_MODELING_AND_BASELINES.md").write_text(p5_content, encoding="utf-8")

    # -------------------------------------------------------------
    # PHASE 6: Robustness, Generalization & Statistical Inference
    # -------------------------------------------------------------
    p6_content = """# Phase 6: Robustness, Generalization & Statistical Inference

**Consolidated Steps**: Step 11 (Robustness Audit), Step 11.1 (Statistical Multiplicity Audit)  
**Primary Deliverables**: `reports/11_1_paired_hypothesis_tests.csv`, `reports/11_1_paired_bootstrap_differences.csv`, `reports/11_leave_one_family_out.csv`  

---

## 1. Scientific Objective
Evaluate the statistical significance and sensitivity of ESM-2 predictive performance under leave-one-family-out cross-validation, sister-species transfer, paired exact McNemar tests, 10,000 paired bootstrap resamples, and Holm-Bonferroni Family-Wise Error Rate (FWER) multiplicity control.

## 2. Statistical Multiplicity Audit
Paired comparisons against Ridge ESM-2 on family predictions ($N=63$):
- **vs. Majority Baseline**: Discordant $b=26, c=5$, raw $p = 0.000192$, Holm **$p_{\\text{adj}} = 0.000769$** (Statistically Significant).
- **vs. Physicochemical Random Forest**: Discordant $b=23, c=9$, raw $p = 0.020062$, Holm **$p_{\\text{adj}} = 0.060185$** (Borderline / Not Significant at $\\alpha=0.05$).
- **vs. Physicochemical Logistic Regression**: Discordant $b=22, c=9$, raw $p = 0.029449$, Holm **$p_{\\text{adj}} = 0.060185$** (Borderline / Not Significant at $\\alpha=0.05$).
- **vs. Combined Model**: Discordant $b=2, c=0$, raw $p = 0.500000$, Holm **$p_{\\text{adj}} = 0.500000$** (No difference).

## 3. Paired Bootstrap Confidence Intervals (10,000 Resamples)
- $\\Delta \\text{Macro-F1}$ vs Random Forest: Mean $= +0.2844$, 95% CI $[+0.0915, +0.4626]$.
- $\\Delta \\text{Macro-F1}$ vs Logistic Regression: Mean $= +0.2267$, 95% CI $[+0.0540, +0.3976]$.
- Strictly positive bootstrap intervals confirm a robust empirical trend despite missing strict FWER significance due to $N=63$.

## 4. Sensitivity Analyses
- **Leave-One-Family-Out (LOFO)**: 63 LOFO runs yielded Mean Macro-F1 = 0.7946 (Range: [0.778, 0.811]), proving no individual family drives performance.
- **Sister-Species Ortholog Transfer**: Evaluating an HSV-1-trained model on 7 HSV-2 ortholog families yielded 7/7 (100%) correct classification.

## 5. Validation & Test Outcomes
- Validated under `tests/test_11_robustness_audit.py` and `test_11_1_statistical_audit.py` (15 passing tests).
"""
    (PHASES_DIR / "PHASE_06_ROBUSTNESS_AND_STATISTICAL_INFERENCE.md").write_text(p6_content, encoding="utf-8")

    # -------------------------------------------------------------
    # PHASE 7: Representation Interpretation & Embedding Geometry
    # -------------------------------------------------------------
    p7_content = """# Phase 7: Representation Interpretation & Embedding Geometry

**Consolidated Steps**: Step 12 (Representation Interpretation), Step 12.1 (Geometry & PCA Audit)  
**Primary Deliverables**: `reports/12_1_embedding_geometry.csv`, `reports/12_1_pca_heldout_metrics.csv`  

---

## 1. Scientific Objective
Investigate the intrinsic dimensionality of the 480-dimensional ESM-2 embedding space, test whether low-dimensional linear projections generalize out-of-fold, and determine whether embeddings reflect causal kinetic mechanisms or distributed structural signatures.

## 2. Embedding Geometry & Effective Rank
- **Roy-Vetterli Effective Rank**:
  - *Standardized Correlation Matrix*: **13.63** (captures distributed degrees of freedom).
  - *Raw Unstandardized Covariance*: **5.93** (dominated by top 1-5 large variance directions).
- **Variance Distribution**: The first 18 principal components capture 90.0% of variance; 30 components capture 95.1%.

## 3. Low-Rank PCA Generalization Collapse
- In Step 12, an 18-component PCA classifier achieved 89.1% training accuracy.
- In Step 12.1, held-out 5-fold cross-validation revealed that out-of-fold Macro-F1 collapsed to **0.4388**.
- **Scientific Conclusion**: Unsupervised linear dimensionality reduction discards subtle discriminative coordinates necessary for out-of-fold viral phase discrimination. Viral replication timing is non-linear and distributed across high-dimensional latent space.

## 4. Non-Causal Boundary
ESM-2 embeddings correlate with globular catalytic domains in regulatory enzymes and transmembrane/repeat regions in structural proteins. The model leverages structural topology correlates, **not** causal transcriptional cascade mechanisms.

## 5. Validation & Test Outcomes
- Validated under `tests/test_12_representation_interpretation.py` and `test_12_1_interpretation_audit.py` (14 passing tests).
"""
    (PHASES_DIR / "PHASE_07_REPRESENTATION_INTERPRETATION_AND_PCA.md").write_text(p7_content, encoding="utf-8")

    # -------------------------------------------------------------
    # PHASE 8: Independent Biological Validation & Reconciliation
    # -------------------------------------------------------------
    p8_content = """# Phase 8: Independent Biological Validation & Reconciliation

**Consolidated Steps**: Step 13 (Biological Validation & Zero-Hardcoding), Step 13.1 (Discrepancy Reconciliation)  
**Primary Deliverables**: `reports/13_1_discrepancy_reconciliation.md`, `reports/13_1_independently_recomputed_metrics.csv`  

---

## 1. Scientific Objective
Audit all code, inputs, and outputs across the research pipeline to enforce zero-hardcoding, confirm deterministic reproducibility, and resolve historical narrative discrepancies that arose across earlier project steps.

## 2. Reconciled Discrepancies
1. **RS1 Sequence Count & Fold 0 Isolation**: Resolved ambiguity between 212 direct benchmark sequences vs. 129 cluster metadata entries. Proven that exactly 212 valid full-length sequences belong to `RS1` and all are strictly isolated in Fold 0.
2. **Model Hyperparameters**: Confirmed $C=1.0$, `penalty='l2'`, `solver='lbfgs'`, `max_iter=1000`, `random_state=42` across all Logistic Regression baselines.
3. **Feature Dimensionality**: Confirmed $D=8$ for physical, $D=480$ for ESM-2, and $D=488$ for combined models.
4. **Zero-Hardcoding**: All metrics derived by dynamically executing model training and prediction evaluation.

## 3. Validation & Test Outcomes
- Validated under `tests/test_13_independent_validation.py` and `test_13_1_discrepancy_reconciliation.py` (16 passing tests).
- All 687 out-of-fold prediction rows reproduced with absolute floating-point agreement ($< 1.11 \\times 10^{-16}$).
"""
    (PHASES_DIR / "PHASE_08_INDEPENDENT_VALIDATION_AND_RECONCILIATION.md").write_text(p8_content, encoding="utf-8")

    # -------------------------------------------------------------
    # PHASE 9: Manuscript Evidence Audit & Adversarial Peer Review
    # -------------------------------------------------------------
    p9_content = """# Phase 9: Manuscript Evidence Audit & Adversarial Peer Review

**Consolidated Steps**: Step 14 (Manuscript Evidence Audit), Step 15 (Adversarial Peer Review Challenge)  
**Primary Deliverables**: `reports/15_independent_reviewer_report.md`, `reports/15_revised_manuscript.md`, `reports/15_submission_readiness_decision.md`  

---

## 1. Scientific Objective
Execute an independent, adversarial peer review challenging central claims, verifying reference support, testing manuscript completeness, logging required revisions, and revising the research draft.

## 2. Adversarial Claim Verification
Seven claims evaluated in `reports/15_independent_claim_verification.csv`:
- Binary classification task (PRE vs POST) confirmed; fine-grained 4-class discovery refuted.
- Majority baseline superiority verified ($p_{\\text{adj}} = 0.000769$).
- Physicochemical superiority classified as borderline / unproven at $\\alpha=0.05$ ($p_{\\text{adj}} = 0.0602$).
- Low-rank PCA generalization collapse documented as an essential finding.
- Sister-species transfer reclassified strictly as an ortholog sequence transfer diagnostic.
- External clinical isolate validation disclosed as absent.

## 3. Bibliographic Verification
All 7 cited primary references (Honess & Roizman 1974, 1975; Rutkowski et al. 2015; Kulej et al. 2017; Dremel & DeLuca 2019; Roizman et al. 2013; Lin et al. 2023) verified with valid DOIs, PMIDs, and accessible primary assays.

## 4. Readiness Decision
- Outcome: **1. READY FOR CO-AUTHOR REVIEW**.
- Explicitly separated from journal submission readiness pending publication figure compilation and supplementary packaging.

## 5. Validation & Test Outcomes
- Validated under `tests/test_14_manuscript_evidence_audit.py` and `test_15_independent_reviewer_audit.py` (18 passing tests).
"""
    (PHASES_DIR / "PHASE_09_EVIDENCE_AUDIT_AND_PEER_REVIEW.md").write_text(p9_content, encoding="utf-8")

    # -------------------------------------------------------------
    # PHASE 10: Publication Visuals & Illustrated Manuscript
    # -------------------------------------------------------------
    p10_content = """# Phase 10: Publication Visuals & Illustrated Manuscript

**Consolidated Step**: Step 16 (Publication-Quality Tables & Figures)  
**Primary Deliverables**: `figures/` (Figures 1-8 in PDF, PNG, CSV), `tables/` (Tables 1-7, S1-S8 in CSV & MD), `reports/16_illustrated_manuscript.md`  

---

## 1. Scientific Objective
Produce a complete, publication-quality suite of vector graphics, formatted data tables, and an illustrated research manuscript integrating all verified callouts and supplementary materials.

## 2. Visual Figures Suite (`figures/`)
- **Figure 1**: Study Workflow & Cohort Pipeline (22,710 $\to$ 18,172 $\to$ 3,111 $\to$ 687 benchmark).
- **Figure 2**: Dataset Composition & Universal Corpus Accounting Flow.
- **Figure 3**: Temporal Benchmark Evidence & Sequence Asymmetry (RS1 = 212 variants).
- **Figure 4**: ESM-2 Embedding Geometry & Low-Rank PCA Generalization Collapse.
- **Figure 5**: Family-Level Model Comparisons & Multiplicity Bounds ($p_{\\text{adj}} = 0.0602$).
- **Figure 6**: Classifier Confusion Matrix (PRE: 73.9%, POST: 85.0%) & Predicted Probability Spectrum.
- **Figure 7**: Cross-Validation Stability, Subgroup Stratification, and Evaluation Unit Inflation (+11.7%).
- **Figure 8 (Supplementary)**: Unsupervised Clustering Landscape & 7/7 Sister-Species Ortholog Transfer.

## 3. Main & Supplementary Tables (`tables/`)
- Tables 1 to 7: Universal dataset accounting, model configurations, main performance, statistical comparisons, temporal evidence summary, sensitivity analysis, unsupervised geometry.
- Tables S1 to S8: 687 benchmark records, 5-fold split manifest, family out-of-fold predictions, sequence out-of-fold predictions, 63-family temporal evidence ledger, 22,710 corpus lineage, LOFO stability, reference claim support.

## 4. Illustrated Manuscript
- Authored `reports/16_illustrated_manuscript.md` integrating all figure callouts, embedded tables, solver configurations, and conclusion sections.

## 5. Validation & Test Outcomes
- Validated under `tests/test_16_publication_visuals.py` (8 passing tests).
"""
    (PHASES_DIR / "PHASE_10_PUBLICATION_VISUALS_AND_ILLUSTRATED_MANUSCRIPT.md").write_text(p10_content, encoding="utf-8")

    print(f"Generated 10 phase documents in {PHASES_DIR}")


if __name__ == "__main__":
    create_phase_documents()
