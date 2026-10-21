"""
Script to generate reports/validation/FILE_REORGANIZATION_MAP.md
Documents every migration, consolidation, and path update in the reorganization.
"""
from pathlib import Path

def generate_file_reorganization_map():
    doc_path = Path("reports/validation/FILE_REORGANIZATION_MAP.md")

    content = """# Repository File Reorganization & Migration Map

This document records the complete migration, consolidation, and mapping of files during the reorganization of `Herpes-Simplex-Virus-Analysis-02` into a publication-ready, reproducible research repository.

---

## 1. Migration Ledger

| Old Path | New Path | Reason | References Updated | Validation Result |
|---|---|---|---|---|
| `hsv1.fasta` (root) | `data/raw/hsv1.fasta` | Redundant root duplicate removed after verifying exact SHA-256 match with canonical raw dataset | `tests/`, `docs/`, `scripts/` | Validated (Exact SHA-256 match `15001269...`) |
| `hsv2.fasta` (root) | `data/raw/hsv2.fasta` | Redundant root duplicate removed after verifying exact SHA-256 match with canonical raw dataset | `tests/`, `docs/`, `scripts/` | Validated (Exact SHA-256 match `86c9d627...`) |
| `reports/01_*.md`, `reports/02_*.md`, `reports/03_*.md` | `docs/phases/PHASE_01_DATA_INVENTORY_AND_CURATION.md` | Consolidated Phase 1 data curation, quality control, and taxonomy documentation | `docs/README.md`, `PROJECT_FILE_GUIDE.md` | Validated (All details, tables, and scripts preserved) |
| `reports/05_*.md` | `docs/phases/PHASE_02_FEATURE_ENGINEERING_AND_EMBEDDINGS.md` | Consolidated Phase 2 biophysical property extraction and ESM-2 embedding extraction | `docs/README.md`, `PROJECT_FILE_GUIDE.md` | Validated (Feature dimensions & manifests preserved) |
| `reports/06_*.md`, `reports/07_*.md` | `docs/phases/PHASE_03_UNSUPERVISED_ANALYSIS_AND_GEOMETRY.md` | Consolidated Phase 3 manifold geometry, intrinsic rank, and clustering documentation | `docs/README.md`, `PROJECT_FILE_GUIDE.md` | Validated (Effective rank 13.63 vs 5.93 documented) |
| `reports/08_*.md` | `docs/phases/PHASE_04_TEMPORAL_EVIDENCE_AND_BENCHMARK_DESIGN.md` | Consolidated Phase 4 kinetic evidence curation and 5-fold family-disjoint split design | `docs/README.md`, `PROJECT_FILE_GUIDE.md` | Validated (687 sequence / 63 family invariants verified) |
| `reports/09_*.md`, `reports/10_*.md` | `docs/phases/PHASE_05_PREDICTIVE_MODELING_AND_BASELINES.md` | Consolidated Phase 5 predictive model training, family-level pooling, and baselines | `docs/README.md`, `PROJECT_FILE_GUIDE.md` | Validated (Exact point metrics 0.7946 Macro-F1 matched) |
| `reports/11_*.md` | `docs/phases/PHASE_06_ROBUSTNESS_AND_STATISTICAL_INFERENCE.md` | Consolidated Phase 6 McNemar paired tests, bootstrap CIs, and Holm multiplicity | `docs/README.md`, `PROJECT_FILE_GUIDE.md` | Validated (p=0.0312 unadjusted, p=0.0602 Holm verified) |
| `reports/12_*.md` | `docs/phases/PHASE_07_REPRESENTATION_INTERPRETATION_AND_PCA.md` | Consolidated Phase 7 latent PCA projections and biophysical property correlation | `docs/README.md`, `PROJECT_FILE_GUIDE.md` | Validated (Linear correlation matrix values verified) |
| `reports/13_*.md` | `docs/phases/PHASE_08_INDEPENDENT_VALIDATION_AND_RECONCILIATION.md` | Consolidated Phase 8 zero-hardcoding audit and benchmark discrepancy resolution | `docs/README.md`, `PROJECT_FILE_GUIDE.md` | Validated (100% data-derived assertions verified) |
| `reports/14_*.md`, `reports/15_*.md` | `docs/phases/PHASE_09_EVIDENCE_AUDIT_AND_PEER_REVIEW.md` | Consolidated Phase 9 adversarial peer review challenge and claim qualification | `docs/README.md`, `PROJECT_FILE_GUIDE.md` | Validated (Reviewer critiques & caveats logged) |
| `reports/16_*.md` | `docs/phases/PHASE_10_PUBLICATION_VISUALS_AND_ILLUSTRATED_MANUSCRIPT.md` | Consolidated Phase 10 visual inventory, table reconciliation, and illustrated draft | `docs/README.md`, `PROJECT_FILE_GUIDE.md` | Validated (Figures 1-8 and Tables 1-8 cataloged) |
| `reports/16_illustrated_manuscript.md` | `manuscript/main_manuscript.md` | Established dedicated journal submission package with latest illustrated manuscript | `docs/SUBMISSION_PACKAGE.md` | Validated (Complete text, abstract, and refs mirrored) |
| `figures/figure_*.*` | `manuscript/figures/` | Curated publication-ready figures for journal submission package | `manuscript/README.md` | Validated (24 figure files mirrored in vector/raster) |
| `tables/table_[0-9]*.*` | `manuscript/tables/` | Curated main tables (1-7) for journal submission package | `manuscript/README.md` | Validated (14 table files mirrored in CSV/Markdown) |
| `tables/table_s*.*` | `manuscript/supplementary_materials/` | Curated supplementary tables (S1-S8) for submission package | `manuscript/README.md` | Validated (8 supplementary tables mirrored) |
| Architecture Choices | `docs/decisions/DECISION_01` to `05` | Formalized 5 Architectural Decision Records for core scientific choices | `docs/README.md` | Validated (Curative rules, CV design, Holm test, rank) |

---

## 2. Retention and Preservation Policy

To maintain a complete scientific audit trail while providing a streamlined interface for new researchers:
1. **Raw Data Preservation**: `data/raw/` remains untouched and immutable, locked by `RAW_DATA_CHECKSUMS.sha256`.
2. **Intermediate Research Trail**: All historical step reports in `reports/` remain archived to preserve the exact audit trail, while `docs/phases/` provides the consolidated, authoritative narrative for external researchers.
3. **Canonical Scripts**: All Python scripts in `scripts/` retain their canonical filenames to maintain full backward compatibility with the automated test suite.
"""

    with open(doc_path, "w", encoding="utf-8") as f:
        f.write(content.strip() + "\n")
    print(f"Generated {doc_path}")

if __name__ == "__main__":
    generate_file_reorganization_map()
