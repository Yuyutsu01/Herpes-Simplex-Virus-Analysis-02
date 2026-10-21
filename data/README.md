# Data Directory Guide

This directory houses the sequence data, physicochemical features, protein language model embeddings, benchmark datasets, and split manifests for the HSV-1 and HSV-2 temporal proteomics study.

---

## Directory Organization

```text
data/
├── raw/                      # Immutable primary source FASTA files and SHA-256 manifest
│   ├── hsv1.fasta           # Reference proteome for HSV-1 (strain 17)
│   ├── hsv2.fasta           # Reference proteome for HSV-2 (strain HG52)
│   └── RAW_DATA_CHECKSUMS.sha256 # Cryptographic hashes for raw FASTAs
├── processed/                # Curated sequences, embeddings, features, and benchmarks
└── cache/                    # Local processing caches (if retained)
```

## Key Files in `data/processed/`

- `all_curated_proteins.csv`: Consolidated inventory of quality-filtered protein sequences.
- `physicochemical_properties.csv`: 29 sequence-derived biophysical indices (pI, GRAVY, instability, charge, MW, AA frequencies).
- `esm2_t33_650M_embeddings.npz` / `.csv`: Mean-pooled 1280-dim ESM-2 650M sequence embeddings.
- `esm2_embeddings_manifest.json`: Array index to UniProt accession mapping.
- `temporal_evidence_adjudication.json`: Adjudicated kinetic classifications with PubMed citations and experimental quotes.
- `benchmark_dataset_687.csv`: Authoritative 687-sequence, 63-family direct-evidence benchmark.
- `family_disjoint_5fold_splits.csv`: 5-fold cross-validation partition assignments grouped by homology group.

For full schemas and column definitions, see [docs/DATA_DICTIONARY.md](file:///docs/DATA_DICTIONARY.md).
