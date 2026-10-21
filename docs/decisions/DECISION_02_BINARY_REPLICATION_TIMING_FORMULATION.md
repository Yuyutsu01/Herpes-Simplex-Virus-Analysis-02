# Decision 02: Binary Replication Timing Formulation

**Status**: ADOPTED AND ENFORCED  
**Date**: Phase 4 / Step 8 & 9  
**Decision**: Formulate the machine learning task strictly as a binary classification distinguishing pre-replicative regulatory enzymes (IE + E) from post-replicative structural virion components (L1 + L2).

## Context
The canonical herpesvirus cascade comprises four kinetic phases: immediate-early ($lpha$), early ($eta$), leaky-late ($\gamma_1$), and true-late ($\gamma_2$). However, the number of independent gene families in the human alphaherpesvirus proteome is small ($N=63$). Immediate-early proteins comprise only 3 families (`RL2`, `RS1`, `US1`), while leaky-late comprises 17 and true-late 23. A 4-class multi-class benchmark with $N=3$ in one class suffers from extreme class scarcity and unreliable multi-class generalization.

## Decision Rationale
1. **Biochemical Partition**: Delineated by viral DNA replication. Pre-replicative proteins ($lpha + eta$) are synthesized prior to genome replication and function as regulatory transactivators, helicases, primases, and polymerases. Post-replicative proteins ($\gamma_1 + \gamma_2$) require or are accelerated by DNA synthesis and form structural capsids, tegument matrix, and envelope spikes.
2. **Balanced Sample Distribution**: Yields 23 pre-replicative vs 40 post-replicative families, providing statistically viable class support for 5-fold cross-validation.
3. **Refusal of Overclaims**: Prohibit claiming direct 4-class fine-grained kinetic discovery from sequence alone.
