# Step 8 Protocol: Temporal Label Adjudication and Evidence Hierarchy

**Protocol Version**: 1.0.0 (Step 8 Standard)  
**Deliverable**: `data/processed/08_temporal_label_adjudication.csv`  

---

## 1. Evidence Hierarchy Rules

To ensure scientific integrity:
1. **Rule 1 (No Unsupervised Inference)**: No temporal label may be assigned based on ESM-2 embeddings, PCA coordinates, cluster membership, protein length, or molecular weight.
2. **Rule 2 (Distinguish Transcript from Protein)**: mRNA kinetics (e.g. 4sU-seq) are classified as `TRANSCRIPT`, while mass spectrometry / Western blots are classified as `PROTEIN`.
3. **Rule 3 (Honest Insufficient Evidence)**: Records lacking canonical gene mapping are designated `INSUFFICIENT_EVIDENCE`. We reject fabricating predictions to inflate benchmark size.
4. **Rule 4 (Cross-Species Transparency)**: HSV-1 experimental evidence is not falsely claimed as HSV-2 direct evidence; homology-transferred labels are flagged explicitly.

---

## 2. Label Adjudication Decision Tree

```
               [Representative Sequence]
                           |
             Does it have a canonical gene?
             (Direct or Phase B Alias)
                     /            \
                   YES             NO  --> [INSUFFICIENT_EVIDENCE]
                    |
      Is it in Canonical Kinetics Catalog?
              /              \
            YES               NO  --> [INSUFFICIENT_EVIDENCE]
             |
   Check Inhibitor Response:
   - CHX Resistant? -----------> [IE_ALPHA]
   - CHX Sensitive, PAA Res? --> [E_BETA]
   - PAA Partially Inhibited? -> [L1_GAMMA1]
   - PAA Abolished? -----------> [L2_GAMMA2]
   - Variable / Complex Locus? -> [MIXED_OR_CONTEXT_DEPENDENT]
```
