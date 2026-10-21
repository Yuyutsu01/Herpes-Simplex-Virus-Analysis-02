# Reproducibility Guide: HSV-1 & HSV-2 Temporal Proteomics Pipeline

This guide provides end-to-end instructions for reproducing all computational analyses, machine learning models, statistical tests, publication tables, and figures from a clean checkout of the repository.

---

## 1. Environment & Computational Prerequisites

### 1.1 Operating Systems & Python Architecture
- **Validated Environments**:
  - Windows 11 (64-bit), Python 3.11.8
  - Linux (Ubuntu 22.04 LTS), Python 3.11.x
  - macOS (Ventura / Sonoma, Apple Silicon via Rosetta or Native), Python 3.11.x
- **Python Version**: Python 3.10+ (Python 3.11 recommended).

### 1.2 Hardware Requirements
- **Feature Extraction & Modeling Pipeline**:
  - CPU: Standard modern multi-core processor (x86_64 or ARM64).
  - RAM: Minimum 8 GB (16 GB recommended for full matrix operations).
  - Disk Space: ~500 MB for repository, scripts, and precomputed embeddings.
  - Execution Time: Under 2 minutes for complete model training, statistical tests, and figure generation.
- **ESM-2 650M Embedding Recomputation (Optional)**:
  - If regenerating embeddings from raw sequence files from scratch rather than using saved artifacts:
    - GPU: NVIDIA GPU with $\ge 12$ GB VRAM (e.g., RTX 3060, RTX 4080, A100) or CPU execution (requires ~16 GB RAM and ~20 minutes).
    - CUDA: 11.8 or 12.x.

### 1.3 Dependencies Installation
To configure a clean virtual environment:

```bash
# Create a virtual environment
python -m venv venv

# Activate virtual environment
# Windows (PowerShell):
.\venv\Scripts\Activate.ps1
# Linux / macOS:
source venv/bin/activate

# Upgrade pip
pip install --upgrade pip

# Install project dependencies
pip install -r requirements.txt
```

---

## 2. Determinism & Random Seeds

All stochastic procedures (cross-validation fold splitting, logistic regression solver initialization, random forest bootstrap sampling, and empirical percentile confidence interval bootstrap) use a fixed seed:
- **Global Seed**: `42`
- **Scikit-Learn Solvers**: `random_state=42`
- **NumPy Random Generator**: `np.random.seed(42)`
- **Bootstrap Replicates**: $B = 10,000$ iterations with fixed generator states.

### Numerical Tolerance Bounds
- Point estimates (e.g., ESM-2 family Macro-F1 = 0.7946, Balanced Accuracy = 0.7946, Concordant Families = 50 / 63) are fully deterministic and bit-exact across platforms.
- 95% bootstrap confidence intervals are stable to within $\pm 0.005$ across standard scientific Python distributions.

---

## 3. Verifying Raw Source Data Integrity

Before executing any analytical pipeline, verify the integrity of the immutable raw FASTA datasets:

```bash
# Verify cryptographic hashes
python scripts/verify_raw_checksums.py
```
Expected output:
- `data/raw/hsv1.fasta`: `15001269a835b41223467e0d91bdcada433ee68c6a58c19e684c060483b31f36` [OK]
- `data/raw/hsv2.fasta`: `86c9d627286aa4cad121af4c84cdb39b824c14905aa57948de9a996955ca0311` [OK]

---

## 4. End-to-End Execution Workflow

To reproduce the study from sequence ingestion through to final publication artifacts, execute scripts in the following canonical sequence:

```bash
# -------------------------------------------------------------
# Phase 1: Data Inventory, Curation, and Quality Control
# -------------------------------------------------------------
# Ingest raw FASTA sequences, filter non-canonical residues, curate proteomes
python scripts/01_curate_proteomes.py

# -------------------------------------------------------------
# Phase 2: Feature Extraction & Embedding Generation
# -------------------------------------------------------------
# Compute physicochemical indices (pI, GRAVY, instability, charge, secondary structure)
python scripts/02_extract_physicochemical_properties.py

# (Optional: Only if re-extracting ESM-2 650M representations)
python scripts/03_extract_esm2_embeddings.py

# -------------------------------------------------------------
# Phase 3 & 4: Benchmark Curation & Split Design
# -------------------------------------------------------------
# Adjudicate literature evidence and generate family-disjoint 5-fold partitions
python scripts/04_build_temporal_benchmark.py

# -------------------------------------------------------------
# Phase 5: Supervised Model Training & Evaluation
# -------------------------------------------------------------
# Train baseline majority, physicochemical, and ESM-2 models under 5-fold CV
python scripts/05_train_models.py

# -------------------------------------------------------------
# Phase 6: Statistical Inference & Paired Hypothesis Testing
# -------------------------------------------------------------
# Compute exact McNemar tests, bootstrap CIs, and Holm-adjusted p-values
python scripts/06_evaluate_models.py

# -------------------------------------------------------------
# Phase 7: Representation Interpretation & PCA
# -------------------------------------------------------------
# Decompose embedding covariance/correlation matrices and assess alignment
python scripts/07_interpretability_analysis.py

# -------------------------------------------------------------
# Phase 8: Publication Visuals & Manuscript Tables
# -------------------------------------------------------------
# Generate all main and supplementary publication figures (PDF/PNG) and tables (CSV)
python scripts/16_generate_publication_visuals.py
```

---

## 5. Automated Test Suite Validation

The repository includes a comprehensive automated test suite enforcing scientific invariants, data alignment, zero hardcoding, statistical validity, and repository structural integrity.

To execute all tests:

```bash
# Run the complete test suite with verbose reporting
pytest tests/ -v
```

Expected result:
- **100% passing tests** (no failures, no skipped tests).
- Automated checks verify that:
  - No sequence leaks across folds within the same homology group.
  - Metrics in manuscript tables strictly match out-of-fold prediction ledgers.
  - Zero hardcoded results exist in test verification scripts.
  - File links, manifests, and documentation references resolve cleanly.

---

## 6. Clean-Clone Quickstart Summary

For a collaborator or reviewer inspecting the repository for the first time:

```bash
git clone https://github.com/Yuyutsu01/Herpes-Simplex-Virus-Analysis-02.git
cd Herpes-Simplex-Virus-Analysis-02
python -m venv venv
.\venv\Scripts\Activate.ps1  # or source venv/bin/activate
pip install -r requirements.txt
pytest tests/ -v
```
