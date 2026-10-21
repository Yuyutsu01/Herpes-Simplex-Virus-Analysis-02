"""
Script: generate_decision_records.py
Purpose: Create authoritative decision records in docs/decisions/ documenting
  key scientific, architectural, and methodological decisions.
"""

from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
DECISIONS_DIR = BASE_DIR / "docs" / "decisions"
DECISIONS_DIR.mkdir(parents=True, exist_ok=True)


def create_decisions():
    print("=== Generating Decision Records in docs/decisions/ ===")

    # Decision 1
    d1 = """# Decision 01: Group-Disjoint Cross-Validation & RS1 Isolation

**Status**: ADOPTED AND ENFORCED  
**Date**: Phase 4 / Step 9  
**Decision**: Enforce strict group-disjoint cross-validation at the homologous gene family level and isolate the entire RS1 (ICP4) gene family within Fold 0 validation.

## Context
Viral genomes encode gene duplications and numerous clinical strain variants in public repositories. In our dataset, the immediate-early transactivator `RS1` contains 212 near-identical sequence variants, representing 30.9% of the direct benchmark (212/687). Naive sequence-level random splitting places near-identical variants of `RS1` in both training and test folds, creating catastrophic pseudoreplication leakage where the model merely memorizes strain variants.

## Decision Rationale
1. **Zero Family Overlap**: Partition all sequences strictly by homologous gene family (`homology_group_id`) across 5 folds ($\text{Overlap} = 0$).
2. **RS1 Isolation**: Confine all 212 `RS1` sequences to Fold 0 validation. Zero `RS1` variants exist in Fold 0 training.
3. **Primary Evaluation Unit**: Report performance aggregated at the gene-family level ($N=63$) as the primary metric, preventing sequence-count weighting from distorting generalization estimates.
"""
    (DECISIONS_DIR / "DECISION_01_GROUP_DISJOINT_CV_AND_RS1_ISOLATION.md").write_text(d1, encoding="utf-8")

    # Decision 2
    d2 = """# Decision 02: Binary Replication Timing Formulation

**Status**: ADOPTED AND ENFORCED  
**Date**: Phase 4 / Step 8 & 9  
**Decision**: Formulate the machine learning task strictly as a binary classification distinguishing pre-replicative regulatory enzymes (IE + E) from post-replicative structural virion components (L1 + L2).

## Context
The canonical herpesvirus cascade comprises four kinetic phases: immediate-early ($\alpha$), early ($\beta$), leaky-late ($\gamma_1$), and true-late ($\gamma_2$). However, the number of independent gene families in the human alphaherpesvirus proteome is small ($N=63$). Immediate-early proteins comprise only 3 families (`RL2`, `RS1`, `US1`), while leaky-late comprises 17 and true-late 23. A 4-class multi-class benchmark with $N=3$ in one class suffers from extreme class scarcity and unreliable multi-class generalization.

## Decision Rationale
1. **Biochemical Partition**: Delineated by viral DNA replication. Pre-replicative proteins ($\alpha + \beta$) are synthesized prior to genome replication and function as regulatory transactivators, helicases, primases, and polymerases. Post-replicative proteins ($\gamma_1 + \gamma_2$) require or are accelerated by DNA synthesis and form structural capsids, tegument matrix, and envelope spikes.
2. **Balanced Sample Distribution**: Yields 23 pre-replicative vs 40 post-replicative families, providing statistically viable class support for 5-fold cross-validation.
3. **Refusal of Overclaims**: Prohibit claiming direct 4-class fine-grained kinetic discovery from sequence alone.
"""
    (DECISIONS_DIR / "DECISION_02_BINARY_REPLICATION_TIMING_FORMULATION.md").write_text(d2, encoding="utf-8")

    # Decision 3
    d3 = """# Decision 03: Holm-Bonferroni Multiplicity Control

**Status**: ADOPTED AND ENFORCED  
**Date**: Phase 6 / Step 11  
**Decision**: Enforce Holm-Bonferroni step-down Family-Wise Error Rate (FWER) correction across all paired model comparisons at $\alpha = 0.05$.

## Context
When comparing ESM-2 against classical baselines (Majority, Physicochemical Random Forest, Physicochemical Logistic Regression, Combined), nominal unadjusted p-values showed $p = 0.0201$ vs RF and $p = 0.0294$ vs LogReg. In informal machine learning literature, these are frequently reported with asterisks ($*$) claiming "statistical superiority".

## Decision Rationale
1. **Multiplicity Discipline**: Testing 4 paired hypotheses on the same 63 gene families inflates false positive risk.
2. **Holm Adjustment**: Ordering the p-values and adjusting by rank yields $p_{\text{adj}} = 0.060185$ for both physicochemical models.
3. **Scientific Boundary**: Because $0.0602 > 0.05$, formal superiority cannot be claimed at $\alpha=0.05$. It must be described as a strong empirical trend with positive bootstrap intervals ($[+0.0915, +0.4626]$ vs RF; $[+0.0540, +0.3976]$ vs LogReg).
"""
    (DECISIONS_DIR / "DECISION_03_HOLM_BONFERRONI_MULTIPLICITY_CONTROL.md").write_text(d3, encoding="utf-8")

    # Decision 4
    d4 = """# Decision 04: Disclosure of PCA Generalization Collapse

**Status**: ADOPTED AND ENFORCED  
**Date**: Phase 7 / Step 12  
**Decision**: Transparently report the out-of-fold generalization collapse of low-rank PCA classifiers ($0.4388$ Macro-F1) rather than suppressing the negative finding or reporting in-sample training accuracy ($89.10\%$).

## Context
In Step 12, an exploratory logistic regression model fit on an 18-component PCA projection achieved 89.10% training accuracy. However, evaluating this same PCA pipeline under 5-fold family-held-out cross-validation caused family Macro-F1 to collapse to $0.4388$ (worse than majority baseline).

## Decision Rationale
1. **Epistemological Honesty**: Reporting $89.10\%$ training accuracy as a "linear low-rank compression" would be a fraudulent misrepresentation of model generalization.
2. **Scientific Insight**: The held-out collapse provides valuable biological evidence: viral replication timing coordinates cannot be compressed into linear top-variance components. The discriminative information is distributed across non-linear high-dimensional embedding space.
"""
    (DECISIONS_DIR / "DECISION_04_DISCLOSURE_OF_PCA_GENERALIZATION_COLLAPSE.md").write_text(d4, encoding="utf-8")

    # Decision 5
    d5 = """# Decision 05: Sister-Species Ortholog Transfer Framing

**Status**: ADOPTED AND ENFORCED  
**Date**: Phase 6 / Step 11 & Phase 9 / Step 14  
**Decision**: Reclassify the 100% (7/7) cross-species prediction on HSV-2 strictly as a "sister-species ortholog sequence transfer diagnostic", explicitly refusing to frame it as universal cross-species generalization to unobserved viral clades.

## Context
A model trained on HSV-1 was evaluated on 7 HSV-2 gene families (`RL2`, `RS1`, `UL19`, `UL23`, `UL27`, `US4`, `US6`) and achieved 7/7 (100%) accuracy. Early drafts proposed describing this as "proof of cross-species generalization".

## Decision Rationale
1. **Ortholog Identity**: The 7 HSV-2 genes are direct sister-species orthologs of core HSV-1 proteins present during training, sharing $60\text{--}90\%$ sequence identity.
2. **Scope Limitation**: Transfer across sister species within the same *Simplexvirus* genus does not prove generalization to distant betaherpesviruses (e.g. CMV), gammaherpesviruses (e.g. EBV), or novel unobserved viral families.
3. **Defensible Phrasing**: Frame as an ortholog sequence transfer diagnostic confirming high representation fidelity between closely related sister proteomes.
"""
    (DECISIONS_DIR / "DECISION_05_SISTER_SPECIES_ORTHOLOG_TRANSFER_FRAMING.md").write_text(d5, encoding="utf-8")

    print(f"Generated 5 decision records in {DECISIONS_DIR}")


if __name__ == "__main__":
    create_decisions()
