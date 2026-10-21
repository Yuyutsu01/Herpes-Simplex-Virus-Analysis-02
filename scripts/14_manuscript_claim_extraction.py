"""
Script: 14_manuscript_claim_extraction.py
Purpose: Locate, inventory, and extract substantive scientific claims from all project reports,
  figures, tables, and manuscript-equivalent documents in Herpes-Simplex-Virus-Analysis-02.
  Produces reports/14_manuscript_inventory.csv and structured claim records.
"""

from pathlib import Path
import hashlib
import json
import os
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent
REPORTS_DIR = BASE_DIR / "reports"
DATA_DIR = BASE_DIR / "data" / "processed"


def compute_sha256(filepath: Path) -> str:
    """Compute SHA-256 checksum of a file."""
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()


def inventory_manuscript_and_reports() -> pd.DataFrame:
    """
    Inventory all reports, figures, tables, and documents in the repository.
    Explicitly documents the presence of milestone technical reports and the status
    of any standalone manuscript document.
    """
    inventory_records = []

    # Check for standalone manuscript files in root, reports, or doc folders
    manuscript_candidates = [
        "manuscript.md", "manuscript.tex", "manuscript.docx", "draft.md",
        "paper.md", "paper.tex", "preprint.md"
    ]
    standalone_found = False
    for cand in manuscript_candidates:
        cand_path = BASE_DIR / cand
        if cand_path.exists():
            standalone_found = True
            inventory_records.append({
                "file_path": cand,
                "file_category": "STANDALONE_MANUSCRIPT",
                "size_bytes": cand_path.stat().st_size,
                "sha256": compute_sha256(cand_path),
                "is_authoritative_manuscript": True,
                "notes": "Direct standalone manuscript file identified."
            })

    if not standalone_found:
        inventory_records.append({
            "file_path": "[STANDALONE_MANUSCRIPT_ROOT]",
            "file_category": "MISSING_STANDALONE_DRAFT",
            "size_bytes": 0,
            "sha256": "N/A",
            "is_authoritative_manuscript": False,
            "notes": (
                "No standalone manuscript file (e.g. manuscript.md/tex/docx) was pre-committed to repository root. "
                "The authoritative scientific narrative, claims, methods, results, and citations are distributed "
                "across the finalized milestone reports: reports/final_dataset_summary.md, reports/10_baseline_results.md, "
                "reports/11_robustness_generalization_audit.md, reports/12_representation_interpretation.md, and "
                "reports/13_1_discrepancy_reconciliation.md. An audited pre-submission manuscript draft will be synthesized "
                "in Step 14 as reports/14_audited_manuscript_draft.md."
            )
        })

    # Inventory all markdown reports in reports/
    for p in sorted(REPORTS_DIR.glob("*.md")):
        rel_path = p.relative_to(BASE_DIR).as_posix()
        cat = "TECHNICAL_MILESTONE_REPORT"
        is_auth = False
        if "13_1" in p.name:
            cat = "FINAL_RECONCILED_AUDIT_REPORT"
            is_auth = True
        elif "final_dataset" in p.name:
            cat = "FROZEN_DATASET_SUMMARY"
            is_auth = True

        inventory_records.append({
            "file_path": rel_path,
            "file_category": cat,
            "size_bytes": p.stat().st_size,
            "sha256": compute_sha256(p),
            "is_authoritative_manuscript": is_auth,
            "notes": f"Milestone report: {p.stem}"
        })

    df_inv = pd.DataFrame(inventory_records)
    out_path = REPORTS_DIR / "14_manuscript_inventory.csv"
    df_inv.to_csv(out_path, index=False)
    print(f"Saved: {out_path} ({len(df_inv)} inventoried items)")
    return df_inv


def extract_primary_claims() -> pd.DataFrame:
    """
    Extract substantive scientific, quantitative, statistical, methodological,
    and biological claims from across the authoritative project narrative.
    """
    claims = [
        # 1. Corpus Accounting & Data Curation
        {
            "claim_id": "CLM-001",
            "manuscript_section": "Data Curation & Sequence Inventory",
            "claim_type": "numerical",
            "exact_manuscript_text": "All 22,710 original records from the source FASTA archives have been assigned to mutually exclusive, biologically defensible cohorts.",
            "reported_value": "22,710 total records",
            "primary_source_file": "reports/final_dataset_summary.md",
            "recomputed_evidence_file": "reports/final_dataset_summary.md",
            "epistemological_verdict": "VERIFIED FROM DATA",
            "correction_needed": "None. Universal record accounting confirms exact 1-to-1 traceability."
        },
        {
            "claim_id": "CLM-002",
            "manuscript_section": "Data Curation & Sequence Inventory",
            "claim_type": "numerical",
            "exact_manuscript_text": "The primary valid full-length cohort comprises 18,172 sequences (9,454 HSV-1 and 8,718 HSV-2), representing 80.02% of the total corpus.",
            "reported_value": "18,172 valid full-length sequences",
            "primary_source_file": "reports/final_dataset_summary.md",
            "recomputed_evidence_file": "data/processed/valid_proteins_combined.fasta",
            "epistemological_verdict": "VERIFIED FROM DATA",
            "correction_needed": "None. Exact FASTA record count verified."
        },
        {
            "claim_id": "CLM-003",
            "manuscript_section": "Data Curation & Sequence Inventory",
            "claim_type": "numerical",
            "exact_manuscript_text": "Redundancy reduction at 95% sequence identity yielded 3,111 unique cluster representatives across HSV-1 and HSV-2.",
            "reported_value": "3,111 unique representatives",
            "primary_source_file": "data/processed/esm2_t12_representatives_metadata.csv",
            "recomputed_evidence_file": "data/processed/esm2_t12_representatives_metadata.csv",
            "epistemological_verdict": "VERIFIED FROM DATA",
            "correction_needed": "None. Verified in Step 4 and preserved across all subsequent steps."
        },
        {
            "claim_id": "CLM-004",
            "manuscript_section": "Data Curation & Sequence Inventory",
            "claim_type": "numerical",
            "exact_manuscript_text": "The canonical annotated proteome consists of 74 shared core gene families across HSV-1 and HSV-2.",
            "reported_value": "74 canonical gene families",
            "primary_source_file": "reports/08_temporal_literature_review.md",
            "recomputed_evidence_file": "reports/08_temporal_literature_review.md",
            "epistemological_verdict": "VERIFIED FROM DATA",
            "correction_needed": "None. Matches standard alphaherpesvirus genomic architecture."
        },

        # 2. Unsupervised Representation & Geometry
        {
            "claim_id": "CLM-005",
            "manuscript_section": "Representation Analysis & Geometry",
            "claim_type": "numerical",
            "exact_manuscript_text": "The Roy-Vetterli effective rank of the standardized correlation matrix of ESM-2 embeddings is 13.6, with 95% of total variance captured by 30 principal components.",
            "reported_value": "Effective rank = 13.6, PCs for 95% var = 30",
            "primary_source_file": "reports/12_1_embedding_geometry.csv",
            "recomputed_evidence_file": "reports/12_1_embedding_geometry.csv",
            "epistemological_verdict": "VERIFIED FROM DATA",
            "correction_needed": "None. Recomputed and confirmed in Step 12.1."
        },
        {
            "claim_id": "CLM-006",
            "manuscript_section": "Representation Analysis & Geometry",
            "claim_type": "methodological",
            "exact_manuscript_text": "Unsupervised HDBSCAN clustering on latent representations recovers coherent gene family structure without sequence alignment.",
            "reported_value": "74 clusters, silhouette = 0.4747",
            "primary_source_file": "reports/07_unsupervised_model_comparison.csv",
            "recomputed_evidence_file": "reports/07_unsupervised_model_comparison.csv",
            "epistemological_verdict": "REPRODUCED BUT LIMITED IN SCOPE",
            "correction_needed": "Qualify that cluster recovery reflects strong sequence homology within families, not discovery of novel functional mechanisms."
        },
        {
            "claim_id": "CLM-007",
            "manuscript_section": "Representation Analysis & Geometry",
            "claim_type": "methodological",
            "exact_manuscript_text": "Low-dimensional PCA projections (k=18) achieve 89.1% classification accuracy.",
            "reported_value": "89.1% accuracy",
            "primary_source_file": "reports/12_representation_interpretation.md",
            "recomputed_evidence_file": "reports/12_1_pca_heldout_metrics.csv",
            "epistemological_verdict": "CONTRADICTED BY EVIDENCE",
            "correction_needed": "CRITICAL CORRECTION: 89.1% was an in-sample training artifact. Under 5-fold cross-validation, held-out Macro-F1 collapses to 43.88%. Must report held-out failure."
        },

        # 3. Temporal Benchmark & Evidence Adjudication
        {
            "claim_id": "CLM-008",
            "manuscript_section": "Temporal Benchmark & Adjudication",
            "claim_type": "numerical",
            "exact_manuscript_text": "The strict direct-evidence temporal benchmark comprises 687 sequence representatives distributed across 63 homologous gene families.",
            "reported_value": "687 sequences, 63 families",
            "primary_source_file": "data/processed/08_1_strict_temporal_benchmark.csv",
            "recomputed_evidence_file": "reports/13_1_record_and_fold_reconciliation.csv",
            "epistemological_verdict": "VERIFIED FROM DATA",
            "correction_needed": "None. Reconciled across Step 8.1, 9.1, 10, 13, and 13.1."
        },
        {
            "claim_id": "CLM-009",
            "manuscript_section": "Temporal Benchmark & Adjudication",
            "claim_type": "biological",
            "exact_manuscript_text": "Temporal kinetic labels are grounded in classic cascade radiolabeling experiments and high-resolution 4sU nascent RNA-seq.",
            "reported_value": "Citations Honess 1974, 1975; Rutkowski 2015",
            "primary_source_file": "reports/08_1_citation_verification.csv",
            "recomputed_evidence_file": "reports/08_1_citation_verification.csv",
            "epistemological_verdict": "VERIFIED FROM EXTERNAL SOURCE",
            "correction_needed": "Maintain explicit distinction between mRNA transcription kinetics (4sU-seq) and protein accumulation."
        },
        {
            "claim_id": "CLM-010",
            "manuscript_section": "Temporal Benchmark & Adjudication",
            "claim_type": "numerical",
            "exact_manuscript_text": "The RS1 (ICP4) gene family contributes 212 sequences, all isolated within Fold 0 validation to eliminate train-test leakage.",
            "reported_value": "212 sequences in Fold 0",
            "primary_source_file": "data/processed/09_group_split_manifest.csv",
            "recomputed_evidence_file": "reports/13_1_record_and_fold_reconciliation.csv",
            "epistemological_verdict": "VERIFIED FROM DATA",
            "correction_needed": "Replaces the erroneous textual mention of 447 sequences in Fold 3 from the draft Step 13 report."
        },

        # 4. Supervised Baseline Modeling
        {
            "claim_id": "CLM-011",
            "manuscript_section": "Supervised Predictive Modeling",
            "claim_type": "statistical",
            "exact_manuscript_text": "Ridge-regularized logistic regression on ESM-2 embeddings achieves a family-level Balanced Accuracy of 79.46% and Macro-F1 of 0.7946 under 5-fold family-held-out cross-validation.",
            "reported_value": "BalAcc = 79.46%, Macro-F1 = 0.7946 (51/63 correct)",
            "primary_source_file": "reports/10_family_level_metrics.csv",
            "recomputed_evidence_file": "reports/13_1_independently_recomputed_metrics.csv",
            "epistemological_verdict": "VERIFIED FROM DATA",
            "correction_needed": "None. Exact match across Step 10, Step 13, and Step 13.1 retraining."
        },
        {
            "claim_id": "CLM-012",
            "manuscript_section": "Supervised Predictive Modeling",
            "claim_type": "statistical",
            "exact_manuscript_text": "Logistic regression on 8 classical physicochemical descriptors achieves family-level Balanced Accuracy of 56.74% and Macro-F1 of 0.5679.",
            "reported_value": "BalAcc = 56.74%, Macro-F1 = 0.5679",
            "primary_source_file": "reports/10_family_level_metrics.csv",
            "recomputed_evidence_file": "reports/13_1_independently_recomputed_metrics.csv",
            "epistemological_verdict": "VERIFIED FROM DATA",
            "correction_needed": "Clarify that the model uses 8 descriptors, not 28."
        },
        {
            "claim_id": "CLM-013",
            "manuscript_section": "Supervised Predictive Modeling",
            "claim_type": "statistical",
            "exact_manuscript_text": "Random Forest on physicochemical descriptors achieves family-level Balanced Accuracy of 51.79% and Macro-F1 of 0.5102.",
            "reported_value": "BalAcc = 51.79%, Macro-F1 = 0.5102",
            "primary_source_file": "reports/10_family_level_metrics.csv",
            "recomputed_evidence_file": "reports/13_1_independently_recomputed_metrics.csv",
            "epistemological_verdict": "VERIFIED FROM DATA",
            "correction_needed": "None. Confirmed across all retraining executions."
        },
        {
            "claim_id": "CLM-014",
            "manuscript_section": "Supervised Predictive Modeling",
            "claim_type": "statistical",
            "exact_manuscript_text": "A fixed in-fold majority classifier achieves family-level Balanced Accuracy of 53.21% and Macro-F1 of 0.4741.",
            "reported_value": "BalAcc = 53.21%, Macro-F1 = 0.4741",
            "primary_source_file": "reports/10_family_level_metrics.csv",
            "recomputed_evidence_file": "reports/13_1_independently_recomputed_metrics.csv",
            "epistemological_verdict": "VERIFIED FROM DATA",
            "correction_needed": "None. Confirmed exact match."
        },
        {
            "claim_id": "CLM-015",
            "manuscript_section": "Supervised Predictive Modeling",
            "claim_type": "statistical",
            "exact_manuscript_text": "Combining physicochemical features with ESM-2 embeddings (488 dimensions) yields family-level Balanced Accuracy of 76.96% and Macro-F1 of 0.7644, conferring no gain over ESM-2 alone.",
            "reported_value": "BalAcc = 76.96%, Macro-F1 = 0.7644",
            "primary_source_file": "reports/10_family_level_metrics.csv",
            "recomputed_evidence_file": "reports/13_1_independently_recomputed_metrics.csv",
            "epistemological_verdict": "VERIFIED FROM DATA",
            "correction_needed": "Clarify that the combined feature dimension is 488, not 508."
        },

        # 5. Statistical Inference & Multiplicity
        {
            "claim_id": "CLM-016",
            "manuscript_section": "Statistical Inference & Multiplicity",
            "claim_type": "comparative",
            "exact_manuscript_text": "ESM-2 significantly outperforms the majority baseline on family predictions (McNemar b=26, c=5, Holm-adjusted p = 0.000769).",
            "reported_value": "b=26, c=5, p_adj = 0.000769",
            "primary_source_file": "reports/13_1_statistical_reconciliation.csv",
            "recomputed_evidence_file": "reports/13_1_statistical_reconciliation.csv",
            "epistemological_verdict": "VERIFIED FROM DATA",
            "correction_needed": "None. Fully survives FWER correction."
        },
        {
            "claim_id": "CLM-017",
            "manuscript_section": "Statistical Inference & Multiplicity",
            "claim_type": "comparative",
            "exact_manuscript_text": "ESM-2 significantly outperforms physicochemical models (Random Forest p = 0.0201; Logistic Regression p = 0.0294).",
            "reported_value": "Raw p = 0.0201 and 0.0294",
            "primary_source_file": "reports/11_paired_model_comparison.csv",
            "recomputed_evidence_file": "reports/13_1_statistical_reconciliation.csv",
            "epistemological_verdict": "INCONCLUSIVE",
            "correction_needed": "CRITICAL CORRECTION: Raw p-values adjust to p_adj = 0.0602 under Holm-Bonferroni correction across 4 comparisons. Must state as suggestive trend, not definitive significance."
        },
        {
            "claim_id": "CLM-018",
            "manuscript_section": "Statistical Inference & Multiplicity",
            "claim_type": "statistical",
            "exact_manuscript_text": "The paired bootstrap 95% confidence intervals for Delta Macro-F1 are strictly positive: [+0.0915, +0.4626] vs Random Forest and [+0.0540, +0.3976] vs Logistic Regression.",
            "reported_value": "Strictly positive bootstrap CIs",
            "primary_source_file": "reports/13_1_statistical_reconciliation.csv",
            "recomputed_evidence_file": "reports/13_1_statistical_reconciliation.csv",
            "epistemological_verdict": "VERIFIED FROM DATA",
            "correction_needed": "Accurately document that bootstrap intervals indicate positive effect size despite borderline FWER p-value."
        },

        # 6. Biological Stratification & Subgroups
        {
            "claim_id": "CLM-019",
            "manuscript_section": "Biological Stratification",
            "claim_type": "biological",
            "exact_manuscript_text": "On directly protein-validated families (N=52), ESM-2 achieves Macro-F1 of 0.8154, whereas on transcript-only annotated families (N=11), Macro-F1 drops to 0.6857.",
            "reported_value": "Protein Macro-F1 = 0.8154, Transcript Macro-F1 = 0.6857",
            "primary_source_file": "reports/10_family_level_metrics.csv",
            "recomputed_evidence_file": "reports/13_1_independently_recomputed_metrics.csv",
            "epistemological_verdict": "SUPPORTED WITHIN SCOPE",
            "correction_needed": "Note that small sample size for transcript-only (N=11) requires caution in drawing broad conclusions."
        },
        {
            "claim_id": "CLM-020",
            "manuscript_section": "Biological Stratification",
            "claim_type": "biological",
            "exact_manuscript_text": "Univariate association tests of individual physicochemical properties (e.g. length, GC content, charge) reveal no single descriptor significantly associated with temporal timing after FDR correction.",
            "reported_value": "0 of 28 features significant after BH-FDR (alpha=0.05)",
            "primary_source_file": "reports/12_1_composition_inference_audit.csv",
            "recomputed_evidence_file": "reports/12_1_composition_inference_audit.csv",
            "epistemological_verdict": "VERIFIED FROM DATA",
            "correction_needed": "Clarify that lack of single-feature significance does not imply bulk properties have zero relevance in multivariate combination."
        },

        # 7. Generalization, Orthologs, and External Validation
        {
            "claim_id": "CLM-021",
            "manuscript_section": "Generalization & Cross-Species Claims",
            "claim_type": "comparative",
            "exact_manuscript_text": "The model demonstrates 100% accuracy (7/7 correct) in cross-species transfer from HSV-1 to HSV-2.",
            "reported_value": "7/7 correct families",
            "primary_source_file": "reports/11_cross_species_transfer_experiment.csv",
            "recomputed_evidence_file": "reports/11_1_cross_species_transfer_audit.csv",
            "epistemological_verdict": "CONTRADICTED BY EVIDENCE",
            "correction_needed": "CRITICAL RECLASSIFICATION: These 7 HSV-2 families are shared core orthologs of genes seen during HSV-1 training. Must be described as 'ortholog sequence transferability between sister species', not external generalization."
        },
        {
            "claim_id": "CLM-022",
            "manuscript_section": "Generalization & Cross-Species Claims",
            "claim_type": "methodological",
            "exact_manuscript_text": "The supervised model has been externally validated on independent clinical cohorts.",
            "reported_value": "External clinical validation",
            "primary_source_file": "reports/13_external_validation_feasibility.md",
            "recomputed_evidence_file": "reports/13_external_validation_feasibility.md",
            "epistemological_verdict": "UNSUPPORTED",
            "correction_needed": "CRITICAL OMISSION: Exactly 0 independent external datasets exist in the repository. Must state that external validation remains to be conducted prospectively."
        }
    ]

    df_claims = pd.DataFrame(claims)
    out_path = REPORTS_DIR / "14_claim_to_evidence_ledger.csv"
    df_claims.to_csv(out_path, index=False)
    print(f"Saved: {out_path} ({len(df_claims)} claims)")
    return df_claims


if __name__ == "__main__":
    print("=== Step 14: Manuscript Claim Extraction ===")
    inventory_manuscript_and_reports()
    extract_primary_claims()
    print("Claim extraction complete.")
