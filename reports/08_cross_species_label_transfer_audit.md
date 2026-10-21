# Step 8 Audit Report: HSV-1 vs HSV-2 Cross-Species Label Transfer

**Audit Date**: 2026-10-21 13:32:33 UTC  
**Deliverable**: `data/processed/08_cross_species_temporal_labels.csv`  

---

## 1. Executive Summary & Transfer Policy

Historically, over 85% of molecular kinetic studies in herpesvirology were performed using **HSV-1** (primarily strains 17, F, or KOS). While HSV-1 and HSV-2 share 74 orthologous genes and conserve the general $\alpha \to \beta \to \gamma$ cascade architecture, assuming 100% identical kinetics across all genes without documentation is scientifically unsound.

### Transfer Categories:
1. **`DIRECTLY_SUPPORTED_IN_BOTH_SPECIES`**:
   - 7 key canonical genes have independently documented kinetics in both HSV-1 and HSV-2:
     - `RS1` (ICP4), `RL2` (ICP0), `UL23` (TK), `UL27` (gB), `UL19` (VP5), `US4` (gG), `US6` (gD).
2. **`INFERRED_FROM_HOMOLOGY`**:
   - 61 core genes have robust experimental data in HSV-1 and strong sequence orthology ($>70\%$ identity), but direct HSV-2 metabolic arrest data is sparse.
   - Recommended for secondary sensitivity benchmarks, not strict primary benchmarks.
3. **`NO_DEFENSIBLE_TRANSFER`**:
   - Genes lacking reliable kinetic data in either species (e.g. `UL55`, `US8A`).
