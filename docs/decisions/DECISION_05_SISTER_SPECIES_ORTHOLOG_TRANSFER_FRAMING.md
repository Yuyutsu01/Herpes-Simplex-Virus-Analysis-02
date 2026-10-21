# Decision 05: Sister-Species Ortholog Transfer Framing

**Status**: ADOPTED AND ENFORCED  
**Date**: Phase 6 / Step 11 & Phase 9 / Step 14  
**Decision**: Reclassify the 100% (7/7) cross-species prediction on HSV-2 strictly as a "sister-species ortholog sequence transfer diagnostic", explicitly refusing to frame it as universal cross-species generalization to unobserved viral clades.

## Context
A model trained on HSV-1 was evaluated on 7 HSV-2 gene families (`RL2`, `RS1`, `UL19`, `UL23`, `UL27`, `US4`, `US6`) and achieved 7/7 (100%) accuracy. Early drafts proposed describing this as "proof of cross-species generalization".

## Decision Rationale
1. **Ortholog Identity**: The 7 HSV-2 genes are direct sister-species orthologs of core HSV-1 proteins present during training, sharing $60	ext{--}90\%$ sequence identity.
2. **Scope Limitation**: Transfer across sister species within the same *Simplexvirus* genus does not prove generalization to distant betaherpesviruses (e.g. CMV), gammaherpesviruses (e.g. EBV), or novel unobserved viral families.
3. **Defensible Phrasing**: Frame as an ortholog sequence transfer diagnostic confirming high representation fidelity between closely related sister proteomes.
