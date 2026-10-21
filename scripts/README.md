# Scripts Directory Guide

This directory houses the executable Python scripts supporting the entire scientific pipeline. Scripts follow a consistent phase-based naming convention and preserve end-to-end execution order.

---

## Canonical Pipeline Execution Order

| Phase | Canonical Script | Primary Function | Key Output |
|---|---|---|---|
| **01** | `01_curate_proteomes.py` | Ingest raw FASTA sequences and filter non-canonical residues | `data/processed/all_curated_proteins.csv` |
| **02** | `02_extract_physicochemical_properties.py` | Extract 29 biophysical indices (pI, GRAVY, charge, etc.) | `data/processed/physicochemical_properties.csv` |
| **02** | `03_extract_esm2_embeddings.py` | Extract 1280-dim embeddings from ESM-2 650M | `data/processed/esm2_t33_650M_embeddings.npz` |
| **03** | `04_build_temporal_benchmark.py` | Adjudicate literature evidence and build 5-fold CV splits | `data/processed/benchmark_dataset_687.csv` |
| **05** | `05_train_models.py` | Train baseline majority, PhysChem LR/RF, ESM-2 LR/RF | `results/predictions/oof_predictions_*.csv` |
| **06** | `06_evaluate_models.py` | Compute family Macro-F1, bootstrap CIs, exact McNemar | `tables/table_3_main_performance.csv` |
| **07** | `07_interpretability_analysis.py` | PCA decomposition and correlation with biophysical features | `tables/table_7_pca_interpretation.csv` |
| **08** | `08_run_independent_validation.py` | Zero-hardcoding verification and pipeline reconciliation | `reports/validation/` |
| **09** | `09_run_adversarial_challenge.py` | Evaluate sensitivity to effective rank and test choices | `tables/table_8_adversarial_reconciliation.csv` |
| **10** | `16_generate_publication_visuals.py` | Generate publication figures (PDF/PNG) and tables | `figures/` and `tables/` |

## Auxiliary Utility Scripts
- `verify_raw_checksums.py`: Verifies cryptographic SHA-256 hashes of immutable raw FASTA data.
- `setup_manuscript.py`: Bundles submission files into `manuscript/`.

All scripts use deterministic random seeds (`seed=42`) and relative path resolution from the project root.
