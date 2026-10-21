# Automated Test Suite

This directory contains the automated test suite for the `Herpes-Simplex-Virus-Analysis-02` project. The test suite enforces scientific invariants, data integrity, statistical validity, zero hardcoding, and repository structure.

---

## Test Organization

- `test_01_raw_fasta_integrity.py`: Validates raw FASTA checksums, accession formats, and record counts.
- `test_02_sequence_curation_qc.py`: Asserts quality filtering, non-canonical residue handling, and length thresholds.
- `test_03_physicochemical_properties.py`: Asserts non-null values and valid biophysical domains for all 29 indices.
- `test_04_esm2_embeddings.py`: Asserts embedding matrix dimensions, finite values, and manifest alignment.
- `test_05_unsupervised_analysis.py`: Validates PCA coordinates, scree calculations, and manifold invariants.
- `test_06_temporal_evidence_audit.py`: Asserts PubMed citations, evidence grading, and kinetic classification logic.
- `test_07_benchmark_dataset_integrity.py`: Asserts 687 rows, 63 homology groups, and class balance.
- `test_08_family_disjoint_splits.py`: Asserts zero homology group overlap across the 5 cross-validation folds.
- `test_09_model_training_reproducibility.py`: Validates deterministic model convergence and seed stability.
- `test_10_statistical_inference_audit.py`: Recomputes McNemar tables, discordant pairs, and Holm-Bonferroni math.
- `test_11_robustness_analysis.py`: Asserts model sensitivity under feature ablation and label flips.
- `test_12_representation_interpretation.py`: Validates PCA effective rank math and biophysical correlation matrices.
- `test_13_independent_validation.py`: Asserts zero hardcoding in pipeline results and ledgers.
- `test_14_evidence_audit_consistency.py`: Asserts metric alignment between manuscript tables and prediction files.
- `test_15_independent_reviewer_audit.py`: Validates adversarial challenges, error reporting, and claim qualification.
- `test_16_publication_visuals.py`: Asserts existence, non-zero size, and integrity of all publication figures and tables.
- `test_repository_structure.py`: Validates repository file layout, internal markdown links, and raw data immutability.

## Running Tests

Execute the full suite using pytest:

```bash
pytest tests/ -v
```
