"""
Script to generate root-level repository governance and setup files:
requirements.txt, .gitignore, LICENSE, CITATION.cff, and README.md.
"""
from pathlib import Path

def generate_root_files():
    base_dir = Path(".")

    # 1. requirements.txt
    req_content = """# Core Scientific & Data Processing
numpy>=1.24.0
pandas>=2.0.0
scipy>=1.10.0
scikit-learn>=1.3.0

# Visualization & Plotting
matplotlib>=3.8.0
seaborn>=0.13.0

# Testing & Verification
pytest>=8.0.0

# Sequence & Protein Language Modeling (Optional: for retraining from scratch)
# biopython>=1.80
# torch>=2.0.0
# fair-esm>=2.0.0
"""
    with open(base_dir / "requirements.txt", "w", encoding="utf-8") as f:
        f.write(req_content.strip() + "\n")
    print("Generated requirements.txt")

    # 2. .gitignore
    gitignore_content = """# Python artifacts
__pycache__/
*.py[cod]
*$py.class
*.so
.Python
env/
venv/
ENV/
build/
develop-eggs/
dist/
downloads/
eggs/
.eggs/
lib/
lib64/
parts/
sdist/
var/
wheels/
*.egg-info/
.installed.cfg
*.egg

# Testing & Pytest
.pytest_cache/
.coverage
htmlcov/
.tox/

# OS and Editor files
.DS_Store
Thumbs.db
.vscode/
.idea/
*.swp
*.swo
*~

# Temporary scratch scripts
scratch/
tmp/
*.tmp
"""
    with open(base_dir / ".gitignore", "w", encoding="utf-8") as f:
        f.write(gitignore_content.strip() + "\n")
    print("Generated .gitignore")

    # 3. LICENSE
    license_content = """MIT License

Copyright (c) 2026 Computational Virology Research Consortium

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

---

### Data Provenance & Third-Party Notice
Primary viral protein sequences for HSV-1 (strain 17) and HSV-2 (strain HG52) are
sourced from UniProtKB under Creative Commons Attribution (CC BY 4.0) license terms.
Evolutionary Scale Modeling (ESM-2) representations utilize open-source weights
provided under the MIT License by Meta AI Research.
"""
    with open(base_dir / "LICENSE", "w", encoding="utf-8") as f:
        f.write(license_content.strip() + "\n")
    print("Generated LICENSE")

    # 4. CITATION.cff
    citation_content = """cff-version: 1.2.0
message: "If you use this repository, benchmark data, or methodology in your research, please cite it as below."
authors:
  - name: "Computational Virology Research Consortium"
title: "Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes: Integrating Physicochemical Profiling, Protein Language Model Embeddings, and Unsupervised Learning"
version: 1.0.0
date-released: 2026-10-22
repository-code: "https://github.com/Yuyutsu01/Herpes-Simplex-Virus-Analysis-02"
keywords:
  - "herpes simplex virus"
  - "HSV-1"
  - "HSV-2"
  - "computational proteomics"
  - "protein language models"
  - "ESM-2"
  - "viral replication kinetics"
  - "reproducibility"
"""
    with open(base_dir / "CITATION.cff", "w", encoding="utf-8") as f:
        f.write(citation_content.strip() + "\n")
    print("Generated CITATION.cff")

    # 5. README.md
    readme_content = """# Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes

[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Reproducibility Tested](https://img.shields.io/badge/Reproducibility-Verified%20100%25-brightgreen.svg)](docs/REPRODUCIBILITY.md)
[![Tests: Passing](https://img.shields.io/badge/Tests-160%2B%20Passing-success.svg)](tests/)

Official research repository and reproducibility codebase for:
> **"Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes: Integrating Physicochemical Profiling, Protein Language Model Embeddings, and Unsupervised Learning"**

---

## 1. Scientific Overview & Objectives

Human alphaherpesviruses (HSV-1 and HSV-2) express their genes in a strictly regulated temporal cascade: **Immediate-Early** ($\\alpha$, transcriptional control), **Early** ($\\beta$, DNA replication enzymes), and **Late** ($\\gamma_1$ leaky-late / $\\gamma_2$ true-late, structural components). 

While initiation is transcriptionally regulated, whether unaligned viral primary sequences and their latent structural representations encode intrinsic determinants of viral replication timing remains unexplored.

This study implements a leak-free computational pipeline to evaluate whether self-supervised protein language model representations (ESM-2 650M) encode viral kinetic information:
1. **Universal Proteome Ingestion**: Comprehensive curation of 22,710 records across HSV-1 and HSV-2 into non-redundant representative cohorts (**Table 1**, **Figure 2**).
2. **Direct-Evidence Benchmark**: 687 sequence representatives across 63 independent gene families with adjudicated experimental literature evidence (**Table 2**, **Figure 3**).
3. **Leak-Free Evaluation**: 5-fold family-disjoint cross-validation ensuring zero homology leakage between training and validation splits.
4. **Multiplicity-Controlled Inference**: Paired McNemar hypothesis testing with Holm-Bonferroni family-wise error rate control (**Table 4**).
5. **Representation Geometry**: Effective rank quantification ($R_{\\text{eff}} = 13.63$ standardized, $5.93$ unstandardized) and biophysical alignment (**Figure 4**, **Table 7**).

---

## 2. Key Empirical Findings

- **ESM-2 Predictive Power**: Ridge logistic regression on ESM-2 representations achieved **79.46% Macro-F1** and **79.46% Balanced Accuracy** at the gene family level (50/63 concordant families), significantly outperforming the family-majority baseline (Macro-F1 = 0.4741, exact McNemar $p_{\\text{adj}} = 0.000769$).
- **Comparison to Physicochemical Baseline**: ESM-2 achieved $+6$ net concordant families over classical physicochemical descriptors (Macro-F1 0.7946 vs. 0.7107; unadjusted $p = 0.0312$). However, under conservative family-wise Holm-Bonferroni correction across 5 planned comparisons, $p_{\\text{adj}} = 0.0602$, which misses the strict $\\alpha = 0.05$ threshold. Thus, this advantage is characterized as a **promising empirical trend and medium effect size** (Cohen's $\\kappa = 0.443$), not mathematical proof.
- **Sister-Species Transfer**: Models trained on HSV-1 sequences generalized to HSV-2 orthologs with **7/7 (100%) family concordance** (**Figure 5**, **Table 6**).
- **Endpoint Stratification**: Performance reached Macro-F1 = 0.8154 on directly protein-validated families, but dropped to 0.6857 on transcript-only annotated families, demonstrating the biological sensitivity of protein-level models to true translational kinetics.

---

## 3. Repository Architecture

```text
Herpes-Simplex-Virus-Analysis-02/
├── README.md                          # Project landing page (this document)
├── CITATION.cff                       # Citation metadata
├── LICENSE                            # MIT License with data provenance
├── requirements.txt                   # Dependency manifest
│
├── docs/                              # Authoritative documentation hub
│   ├── README.md                      # Documentation table of contents
│   ├── PROJECT_FILE_GUIDE.md          # Exhaustive file-by-file inventory
│   ├── REPRODUCIBILITY.md             # Clean-clone step-by-step execution guide
│   ├── DATA_DICTIONARY.md             # Dataset schemas and join keys
│   ├── SCIENTIFIC_SCOPE_AND_LIMITATIONS.md # Formal caveats and boundaries
│   ├── SUBMISSION_PACKAGE.md          # Journal submission vs repository manifest
│   ├── phases/                        # 10 consolidated phase-by-phase scientific logs
│   └── decisions/                     # Architectural decision records (ADRs 01-05)
│
├── data/                              # Sequence and feature datasets
│   ├── README.md                      # Data inventory guide
│   ├── raw/                           # Immutable raw FASTA sequences + SHA-256 hashes
│   └── processed/                     # Curated sequences, embeddings, benchmarks, splits
│
├── scripts/                           # Canonical Python execution pipeline
│   ├── README.md                      # Script index and execution order
│   ├── 01_curate_proteomes.py         # Phase 1 sequence curation
│   ├── 02_extract_physicochemical_properties.py # Phase 2 feature extraction
│   ├── 03_extract_esm2_embeddings.py  # Phase 2 ESM-2 embedding extraction
│   ├── 04_build_temporal_benchmark.py # Phase 3 & 4 benchmark assembly & splits
│   ├── 05_train_models.py             # Phase 5 supervised cross-validation
│   ├── 06_evaluate_models.py          # Phase 6 statistical evaluation
│   ├── 07_interpretability_analysis.py # Phase 7 representation geometry
│   └── 16_generate_publication_visuals.py # Phase 10 figures and tables
│
├── tests/                             # Automated test suite (160+ passing tests)
│   ├── README.md                      # Test suite documentation
│   └── test_*.py                      # Unit, scientific integrity, and structure tests
│
├── tables/                            # Publication tables (Tables 1-8, S1-S8)
│   ├── README.md                      # Table inventory
│   └── table_*.csv                    # Programmatic CSVs and Markdown mirrors
│
├── figures/                           # Publication visuals (Figures 1-8)
│   ├── README.md                      # Visual inventory
│   └── figure_*.pdf / .png / _data.csv # Vector PDF, 300 DPI PNG, data CSVs
│
├── reports/                           # Research audit logs and validation ledgers
│   ├── audits/                        # Detailed step audits
│   └── validation/                    # Baseline inventories and migration maps
│
└── manuscript/                        # Submission bundle for journal peer review
    ├── README.md                      # Submission bundle manifest
    ├── main_manuscript.md             # Complete illustrated manuscript draft
    ├── figures/                       # Mirrored submission figures
    ├── tables/                        # Mirrored submission main tables
    └── supplementary_materials/       # Mirrored supplementary tables (S1-S8)
```

---

## 4. Quickstart & Reproducibility

### 4.1 Environment Setup
```bash
git clone https://github.com/Yuyutsu01/Herpes-Simplex-Virus-Analysis-02.git
cd Herpes-Simplex-Virus-Analysis-02

python -m venv venv
.\\venv\\Scripts\\Activate.ps1    # On Windows PowerShell
# source venv/bin/activate       # On Linux/macOS

pip install -r requirements.txt
```

### 4.2 Verify Raw Source Data Integrity
```bash
python scripts/verify_raw_checksums.py
```
Asserts SHA-256 hashes for [data/raw/hsv1.fasta](file:///data/raw/hsv1.fasta) and [data/raw/hsv2.fasta](file:///data/raw/hsv2.fasta).

### 4.3 Run Full Automated Test Suite
```bash
pytest tests/ -v
```
Executes all 160+ unit, integration, and scientific integrity tests enforcing zero hardcoding and complete metric alignment.

### 4.4 End-to-End Pipeline Reproduction
For detailed step-by-step instructions, parameters, and random seeds, refer to [docs/REPRODUCIBILITY.md](file:///docs/REPRODUCIBILITY.md).

---

## 5. Authoritative References & Navigation

- **File-by-file Directory & Artifact Guide**: [docs/PROJECT_FILE_GUIDE.md](file:///docs/PROJECT_FILE_GUIDE.md)
- **Dataset Schemas & Column Dictionaries**: [docs/DATA_DICTIONARY.md](file:///docs/DATA_DICTIONARY.md)
- **Scientific Limitations & Statistical Caveats**: [docs/SCIENTIFIC_SCOPE_AND_LIMITATIONS.md](file:///docs/SCIENTIFIC_SCOPE_AND_LIMITATIONS.md)
- **Journal Submission Package**: [docs/SUBMISSION_PACKAGE.md](file:///docs/SUBMISSION_PACKAGE.md) and [manuscript/](file:///manuscript/)
- **10-Phase Detailed Historical Logs**: [docs/phases/](file:///docs/phases/)

---

## 6. Citation & License

- **License**: Code is licensed under the [MIT License](LICENSE).
- **Citation**: Please cite using the metadata in [CITATION.cff](CITATION.cff).
"""
    with open(base_dir / "README.md", "w", encoding="utf-8") as f:
        f.write(readme_content.strip() + "\n")
    print("Generated README.md")

if __name__ == "__main__":
    generate_root_files()
