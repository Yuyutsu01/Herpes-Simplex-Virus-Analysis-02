"""
Script to generate docs/README.md for the Herpes-Simplex-Virus-Analysis-02 project.
Serves as the documentation landing page and navigation hub.
"""
from pathlib import Path

def generate_docs_readme():
    doc_path = Path("docs/README.md")
    
    content = """# Documentation Hub: HSV-1 & HSV-2 Proteome Research

Welcome to the documentation repository for the computational proteomics study:
**"Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes: Integrating Physicochemical Profiling, Protein Language Model Embeddings, and Unsupervised Learning."**

This directory houses the complete scientific, technical, operational, and architectural documentation for the project.

---

## 1. Core Reference Documentation

| Document | Purpose |
|---|---|
| [PROJECT_FILE_GUIDE.md](file:///docs/PROJECT_FILE_GUIDE.md) | **Authoritative master inventory** detailing every file, script, dataset, table, and figure in the repository, including inputs, outputs, and reproduction status. |
| [REPRODUCIBILITY.md](file:///docs/REPRODUCIBILITY.md) | Step-by-step instructions for reproducing the entire pipeline from scratch, covering dependencies, hardware, random seeds, and execution sequence. |
| [DATA_DICTIONARY.md](file:///docs/DATA_DICTIONARY.md) | Detailed schema specifications, data types, value constraints, and primary join keys for all CSV, JSON, FASTA, and NPZ datasets. |
| [SCIENTIFIC_SCOPE_AND_LIMITATIONS.md](file:///docs/SCIENTIFIC_SCOPE_AND_LIMITATIONS.md) | Explicit definition of biological scope, observational vs. causal boundaries, statistical multiplicity considerations ($p_{\\text{adj}} = 0.0602$), and effective rank nuances. |
| [SUBMISSION_PACKAGE.md](file:///docs/SUBMISSION_PACKAGE.md) | Clear specification distinguishing the curated journal submission package from the complete open-science reproducibility repository. |

---

## 2. Consolidated Phase Documentation (`docs/phases/`)

The scientific workflow was executed and audited across 10 structured phases. Each document below consolidates all technical steps, rationale, scripts, inputs, outputs, validation checks, and findings for that phase:

1. [PHASE_01_DATA_INVENTORY_AND_CURATION.md](file:///docs/phases/PHASE_01_DATA_INVENTORY_AND_CURATION.md)
   - Primary FASTA ingestion, UniProt accession validation, quality filtering, and non-canonical residue handling.
2. [PHASE_02_FEATURE_ENGINEERING_AND_EMBEDDINGS.md](file:///docs/phases/PHASE_02_FEATURE_ENGINEERING_AND_EMBEDDINGS.md)
   - Bulk physicochemical descriptors and mean-pooled ESM-2 representation extraction.
3. [PHASE_03_UNSUPERVISED_ANALYSIS_AND_GEOMETRY.md](file:///docs/phases/PHASE_03_UNSUPERVISED_ANALYSIS_AND_GEOMETRY.md)
   - Manifold geometry, intrinsic dimensionality estimation (effective rank), and unsupervised clustering.
4. [PHASE_04_TEMPORAL_EVIDENCE_AND_BENCHMARK_DESIGN.md](file:///docs/phases/PHASE_04_TEMPORAL_EVIDENCE_AND_BENCHMARK_DESIGN.md)
   - Literature evidence curation, kinetic adjudication, and strict 5-fold family-disjoint partition construction.
5. [PHASE_05_PREDICTIVE_MODELING_AND_BASELINES.md](file:///docs/phases/PHASE_05_PREDICTIVE_MODELING_AND_BASELINES.md)
   - Training baseline majority, physicochemical LR/RF, and ESM-2 LR/RF models with family-level pooling.
6. [PHASE_06_ROBUSTNESS_AND_STATISTICAL_INFERENCE.md](file:///docs/phases/PHASE_06_ROBUSTNESS_AND_STATISTICAL_INFERENCE.md)
   - Paired McNemar hypothesis testing, 10,000-fold bootstrap confidence intervals, and Holm-Bonferroni correction.
7. [PHASE_07_REPRESENTATION_INTERPRETATION_AND_PCA.md](file:///docs/phases/PHASE_07_REPRESENTATION_INTERPRETATION_AND_PCA.md)
   - Principal component decomposition, linear biophysical correlation analysis, and feature alignment.
8. [PHASE_08_INDEPENDENT_VALIDATION_AND_RECONCILIATION.md](file:///docs/phases/PHASE_08_INDEPENDENT_VALIDATION_AND_RECONCILIATION.md)
   - Full pipeline audit, zero-hardcoding verification, and reconciliation of metric discrepancies.
9. [PHASE_09_EVIDENCE_AUDIT_AND_PEER_REVIEW.md](file:///docs/phases/PHASE_09_EVIDENCE_AUDIT_AND_PEER_REVIEW.md)
   - Independent adversarial peer review challenge, stress testing, and claim qualification.
10. [PHASE_10_PUBLICATION_VISUALS_AND_ILLUSTRATED_MANUSCRIPT.md](file:///docs/phases/PHASE_10_PUBLICATION_VISUALS_AND_ILLUSTRATED_MANUSCRIPT.md)
    - Production of high-resolution publication figures (1-8), publication tables (1-7, S1-S8), and the illustrated manuscript draft.

---

## 3. Architectural Decision Records (`docs/decisions/`)

Key scientific and methodological choices are formally documented in our Architectural Decision Records:

- [DECISION_01_GROUP_DISJOINT_CV_AND_RS1_ISOLATION.md](file:///docs/decisions/DECISION_01_GROUP_DISJOINT_CV_AND_RS1_ISOLATION.md): Strict grouping by homology group to prevent sequence leakage.
- [DECISION_02_BINARY_REPLICATION_TIMING_FORMULATION.md](file:///docs/decisions/DECISION_02_BINARY_REPLICATION_TIMING_FORMULATION.md): Formulation of binary replication timing (PRE vs POST).
- [DECISION_03_HOLM_BONFERRONI_MULTIPLICITY_CONTROL.md](file:///docs/decisions/DECISION_03_HOLM_BONFERRONI_MULTIPLICITY_CONTROL.md): Exact McNemar testing and Holm-Bonferroni correction ($p_{\\text{adj}} = 0.0602$).
- [DECISION_04_DISCLOSURE_OF_PCA_GENERALIZATION_COLLAPSE.md](file:///docs/decisions/DECISION_04_DISCLOSURE_OF_PCA_GENERALIZATION_COLLAPSE.md): Transparent disclosure of PCA performance collapse under held-out cross-validation.
- [DECISION_05_SISTER_SPECIES_ORTHOLOG_TRANSFER_FRAMING.md](file:///docs/decisions/DECISION_05_SISTER_SPECIES_ORTHOLOG_TRANSFER_FRAMING.md): Framing 7/7 sister-species transfer as ortholog conservation rather than novel family generalization.

---

## 4. Documentation Maintenance and Contributing Guidelines

All documentation files are maintained in standard GitHub-flavored Markdown. Any code modifications that alter data schemas, model metrics, or file paths must update the corresponding phase documents, data dictionaries, and [PROJECT_FILE_GUIDE.md](file:///docs/PROJECT_FILE_GUIDE.md) concurrently.
"""

    with open(doc_path, "w", encoding="utf-8") as f:
        f.write(content.strip() + "\n")
    print(f"Generated {doc_path}")

if __name__ == "__main__":
    generate_docs_readme()
