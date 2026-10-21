# Decision 01: Group-Disjoint Cross-Validation & RS1 Isolation

**Status**: ADOPTED AND ENFORCED  
**Date**: Phase 4 / Step 9  
**Decision**: Enforce strict group-disjoint cross-validation at the homologous gene family level and isolate the entire RS1 (ICP4) gene family within Fold 0 validation.

## Context
Viral genomes encode gene duplications and numerous clinical strain variants in public repositories. In our dataset, the immediate-early transactivator `RS1` contains 212 near-identical sequence variants, representing 30.9% of the direct benchmark (212/687). Naive sequence-level random splitting places near-identical variants of `RS1` in both training and test folds, creating catastrophic pseudoreplication leakage where the model merely memorizes strain variants.

## Decision Rationale
1. **Zero Family Overlap**: Partition all sequences strictly by homologous gene family (`homology_group_id`) across 5 folds ($	ext{Overlap} = 0$).
2. **RS1 Isolation**: Confine all 212 `RS1` sequences to Fold 0 validation. Zero `RS1` variants exist in Fold 0 training.
3. **Primary Evaluation Unit**: Report performance aggregated at the gene-family level ($N=63$) as the primary metric, preventing sequence-count weighting from distorting generalization estimates.
