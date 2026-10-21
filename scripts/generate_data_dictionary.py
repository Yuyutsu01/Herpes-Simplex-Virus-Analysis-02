"""
Script to generate docs/DATA_DICTIONARY.md for the Herpes-Simplex-Virus-Analysis-02 project.
Documents schemas, types, descriptions, valid values, and join keys for all critical datasets.
"""
from pathlib import Path

def generate_data_dictionary():
    doc_path = Path("docs/DATA_DICTIONARY.md")
    
    content = """# Data Dictionary: HSV-1 & HSV-2 Proteome Analysis

This document provides the authoritative data dictionary and schema specification for all raw, processed, benchmark, and tabular artifacts in the `Herpes-Simplex-Virus-Analysis-02` repository.

---

## 1. Primary Identifiers & Linkage Architecture

The analytical pipeline links sequence data, physicochemical features, protein language model embeddings, temporal adjudication evidence, and cross-validation folds using the following primary and secondary keys:

| Identifier | Format / Example | Semantic Definition | Authoritative Source |
|---|---|---|---|
| `record_id` / `uniprot_id` | `sp|P04294|KITH_HHV11` / `P04294` | Primary FASTA header / UniProt accession number for individual protein sequence | [data/raw/hsv1.fasta](file:///data/raw/hsv1.fasta), [data/raw/hsv2.fasta](file:///data/raw/hsv2.fasta) |
| `effective_gene_symbol` | `UL23`, `US3`, `RL2` | Standardized canonical viral gene symbol | [data/processed/valid_proteins_metadata.csv](file:///data/processed/valid_proteins_metadata.csv) |
| `effective_species` | `HSV-1` or `HSV-2` | Viral species designation | [data/processed/valid_proteins_metadata.csv](file:///data/processed/valid_proteins_metadata.csv) |
| `homology_group_id` | `UL23`, `US3`, `ICP0` | Homology gene-family cluster ID used for family-disjoint splitting | [tables/table_s2_split_manifest_5fold.csv](file:///tables/table_s2_split_manifest_5fold.csv) |
| `fold_index` | Integer `0` to `4` | 5-fold cross-validation partition assignment (grouped by `homology_group_id`) | [tables/table_s2_split_manifest_5fold.csv](file:///tables/table_s2_split_manifest_5fold.csv) |
| `normalized_temporal_class`| `IE_ALPHA`, `E_BETA`, `L1_GAMMA1`, `L2_GAMMA2` | 4-class detailed kinetic transcription/accumulation stage | [tables/table_s1_supervised_benchmark_687.csv](file:///tables/table_s1_supervised_benchmark_687.csv) |
| `true_class` | `PRE` or `POST` | Binary operationalized replication phase (`PRE` = IE + E; `POST` = L1 + L2) | [tables/table_s3_family_oof_predictions.csv](file:///tables/table_s3_family_oof_predictions.csv) |
| `embedding_row_index` | Integer ($0$ to $3110$) | Row index aligning metadata with ESM-2 embedding tensor | [data/processed/esm2_t12_representatives_metadata.csv](file:///data/processed/esm2_t12_representatives_metadata.csv) |

---

## 2. Raw Sequence Files

### 2.1 [hsv1.fasta](file:///data/raw/hsv1.fasta) & [hsv2.fasta](file:///data/raw/hsv2.fasta)
- **Directory**: `data/raw/`
- **Format**: Standard FASTA (`.fasta`)
- **Immutability**: Immutable source data locked with SHA-256 checksums recorded in [RAW_DATA_CHECKSUMS.sha256](file:///data/raw/RAW_DATA_CHECKSUMS.sha256).
- **Header Structure**: `>sp|Accession|EntryName ProteinDescription OS=OrganismName OX=TaxID GN=GeneName PE=EvidenceLevel SV=Version`
  - Example: `>sp|P04294|KITH_HHV11 Thymidine kinase OS=Human herpesvirus 1 (strain 17) OX=10298 GN=UL23 PE=1 SV=1`
- **Sequence Validation Rules**:
  - Valid IUPAC amino acid characters: standard 20 standard amino acids (`A`, `C`, `D`, `E`, `F`, `G`, `H`, `I`, `K`, `L`, `M`, `N`, `P`, `Q`, `R`, `S`, `T`, `V`, `W`, `Y`).
  - Terminal stop asterisks (`*`) stripped during ingestion.
  - Zero non-canonical or ambiguous characters (`X`, `B`, `Z`, `J`) permitted in the final benchmark.

---

## 3. Curated Proteome & Feature Datasets

### 3.1 [data/processed/curated_protein_candidates.csv](file:///data/processed/curated_protein_candidates.csv)
- **Purpose**: Consolidated inventory of all curated HSV-1 and HSV-2 protein sequences passing length ($\ge 70$ aa) and quality thresholds.
- **Row Grain**: One record per candidate protein sequence.

| Column Name | Type | Definition | Valid Values / Range | Nullable |
|---|---|---|---|---|
| `record_id` | String | Full FASTA sequence header identifier | e.g. `sp\|P04294\|KITH_HHV11` | No |
| `species` | String | Virus species | `HSV-1`, `HSV-2` | No |
| `gene_symbol` | String | Annotated gene symbol | e.g., `UL23` | Yes |
| `sequence_length` | Integer | Length in amino acid residues | $70 - 3,164$ | No |
| `sequence` | String | Validated amino acid sequence | Standard 20 amino acids | No |

### 3.2 [data/processed/valid_proteins_metadata.csv](file:///data/processed/valid_proteins_metadata.csv)
- **Purpose**: Metadata ledger for 18,172 valid primary candidate viral sequences.
- **Row Grain**: One record per valid viral sequence.

| Column Name | Type | Definition | Valid Values / Units | Nullable |
|---|---|---|---|---|
| `record_id` | String | Primary join key to sequence registry | Matches accession header | No |
| `effective_species` | String | Viral species designation | `HSV-1`, `HSV-2` | No |
| `effective_gene_symbol` | String | Canonical gene symbol resolved via alias recovery | e.g. `UL23`, `US3` | No |
| `sequence_length` | Integer | Sequence length in residues | $\ge 70$ | No |
| `is_fragment` | Boolean | True if truncated fragment | `False` | No |

### 3.3 [data/processed/esm2_t12_representatives_embeddings.npy](file:///data/processed/esm2_t12_representatives_embeddings.npy)
- **Purpose**: Mean-pooled sequence representations extracted from pre-trained ESM-2 (`esm2_t12_35M_UR50D`).
- **Row Grain**: One row per non-redundant cluster representative ($N = 3,111$) in matching order to [data/processed/esm2_t12_representatives_metadata.csv](file:///data/processed/esm2_t12_representatives_metadata.csv).
- **Dimensionality**: Matrix of shape $(3111, 480)$, where 480 is the embedding hidden dimension.

---

## 4. Benchmark & Validation Datasets

### 4.1 [tables/table_s1_supervised_benchmark_687.csv](file:///tables/table_s1_supervised_benchmark_687.csv)
- **Purpose**: Authoritative direct-evidence benchmark dataset containing 687 curated sequence representatives across 63 homology groups.
- **Row Grain**: One sequence representative.

| Column Name | Type | Definition | Valid Values / Constraints | Nullable |
|---|---|---|---|---|
| `embedding_row_index` | Integer | Index into `esm2_t12_representatives_embeddings.npy` | $0 - 3110$ | No |
| `record_id` | String | Unique sequence accession header | Unique identifier | No |
| `effective_species` | String | Organism of origin | `HSV-1`, `HSV-2` | No |
| `effective_gene_symbol`| String | Canonical viral gene symbol / homology group | 63 unique families | No |
| `adjudicated_cohort` | String | Benchmark inclusion category | `STRICT_DIRECT_TEMPORAL` | No |
| `normalized_temporal_class` | String | 4-class temporal kinetic classification | `IE_ALPHA`, `E_BETA`, `L1_GAMMA1`, `L2_GAMMA2` | No |
| `transcript_kinetic_class` | String | Transcript-level timing annotation | `IE`, `E`, `L1`, `L2`, `UNASSIGNED` | No |
| `protein_kinetic_class` | String | Direct protein accumulation timing annotation | `IE`, `E`, `L1`, `L2`, `UNASSIGNED` | No |
| `measured_molecule` | String | Experimental assay endpoint | `PROTEIN`, `TRANSCRIPT` | No |
| `evidence_source_ids` | String | Key literature reference / PubMed accession IDs | e.g. `REF_HONESS_1974` | No |

### 4.2 [tables/table_s2_split_manifest_5fold.csv](file:///tables/table_s2_split_manifest_5fold.csv)
- **Purpose**: Strict family-disjoint partition definitions ensuring zero homology group leakage between training and validation splits.
- **Row Grain**: One entry per sequence in the 687-sequence benchmark.

| Column Name | Type | Definition | Constraints |
|---|---|---|---|
| `record_id` | String | Sequence identifier | Matches `table_s1` |
| `effective_gene_symbol`| String | Gene family symbol | 63 families |
| `homology_group_id` | String | Homology cluster key | Homology group isolation |
| `fold_index` | Integer | Fold assignment | $0, 1, 2, 3, 4$ balanced across classes |
| `leakage_isolation_status` | String | Leakage audit result | `CONFIRMED_DISJOINT` |

---

## 5. Model Predictions & Output Ledgers

### 5.1 [tables/table_s3_family_oof_predictions.csv](file:///tables/table_s3_family_oof_predictions.csv)
- **Purpose**: Out-of-fold prediction ledger for all 63 gene families across the 5 evaluated models ($63 \\times 5 = 315$ rows).

| Column Name | Type | Definition | Range / Values |
|---|---|---|---|
| `homology_group_id` | String | Gene family identifier | 63 unique families |
| `effective_gene_symbol`| String | Canonical gene symbol | e.g. `UL23`, `US3` |
| `true_class` | String | Ground-truth binary class | `PRE`, `POST` |
| `predicted_prob_pre` | Float | Predicted probability for `PRE` class | $[0.0, 1.0]$ |
| `measured_molecule` | String | Evidence basis for family | `PROTEIN` ($N=52$), `TRANSCRIPT` ($N=11$) |
| `fold_index` | Integer | Cross-validation fold during which family was held out | $0 - 4$ |
| `model_id` | String | Evaluated machine learning model architecture | `dummy_majority`, `lr_physicochemical`, `rf_physicochemical`, `ridge_esm2`, `rf_esm2` |
| `predicted_class` | String | Thresholded operational prediction ($\tau = 0.5$) | `PRE`, `POST` |
| `is_correct` | Boolean / Int | Family-level concordance indicator | `True` (1) or `False` (0) |

---

## 6. Publication Tables

Located in [tables/](file:///tables/):
- [table_1_dataset_accounting.csv](file:///tables/table_1_dataset_accounting.csv): Universal corpus accounting across 22,710 records.
- [table_2_model_configurations.csv](file:///tables/table_2_model_configurations.csv): Model architectural hyperparameters and regularization settings.
- [table_3_main_performance.csv](file:///tables/table_3_main_performance.csv): Comprehensive comparative model performance across 63 families with bootstrap 95% CIs.
- [table_4_statistical_comparisons.csv](file:///tables/table_4_statistical_comparisons.csv): Paired McNemar hypothesis tests, odds ratios, and Holm-adjusted p-values.
- [table_5_temporal_evidence_summary.csv](file:///tables/table_5_temporal_evidence_summary.csv): Temporal evidence levels, experimental assays, and literature concordance.
- [table_6_sensitivity_analysis.csv](file:///tables/table_6_sensitivity_analysis.csv): Multi-endpoint sensitivity, fragment perturbation, and label flip tests.
- [table_7_unsupervised_summary.csv](file:///tables/table_7_unsupervised_summary.csv): Latent dimensionality, effective rank calculations, and variance explained.
"""
    
    with open(doc_path, "w", encoding="utf-8") as f:
        f.write(content.strip() + "\n")
    print(f"Generated {doc_path}")

if __name__ == "__main__":
    generate_data_dictionary()
