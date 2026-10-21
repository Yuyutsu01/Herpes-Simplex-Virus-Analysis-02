# Step 13: External Validation Feasibility and Generalization Audit

**Project**: *Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes*  
**Date**: 2026-10-21  
**Deliverable**: [`reports/13_external_validation_feasibility.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/13_external_validation_feasibility.md)  

---

## 1. Executive Summary & Epistemological Stance

A central tenet of scientific machine learning is that performance on development data—even under cross-validation—does not guarantee generalization to independent biological populations. Step 13 evaluates whether genuine external validation data exist within the repository or publicly, and outlines the experimental protocols necessary to achieve true independent validation.

### Principal Findings:
1. **Zero Truly Independent External Datasets in Repository**: All 3,111 sequences in this repository derive from historical HSV-1 (strain 17, strain F, KOS) and HSV-2 (strain HG52, 333) FASTA records and associated kinetic literature. No separate, prospectively collected temporal proteomic cohort was frozen.
2. **Reclassification of the HSV-2 Cohort**: The 7 HSV-2 families with direct literature evidence (`RL2`, `RS1`, `UL19`, `UL23`, `UL27`, `US4`, `US6`) are sister-species orthologs of genes seen during HSV-1 training. Evaluating performance on these 7 families tests **ortholog sequence transferability between closely related sister species for known viral genes**, NOT generalization to novel or unseen viral gene families across species.
3. **Prohibition of Synthetic or Inferred 'Validation'**: Using homology-inferred temporal labels as a synthetic test set is rejected as circular. Computational predictions cannot validate computational predictions.

---

## 2. Requirements for a Defensible External Validation Study

To claim genuine external validation in computational virology, a validation dataset must satisfy four strict criteria:

### Criteria for Independence:
- **Independent Experimental Provenance**: Viral infections and temporal profiling must be conducted in an independent laboratory on distinct biological replicates (e.g. fresh primary clinical isolates rather than heavily passaged laboratory strains).
- **Direct Molecule Measurement**: Protein accumulation must be measured directly via quantitative proteomics (e.g., SILAC, TMT, or targeted PRM/SRM mass spectrometry) or de novo translation via Ribo-seq, avoiding reliance on transcript kinetics alone.
- **Independent Label Adjudication**: Kinetic classes must be adjudicated without knowledge of or influence by the model's predictions.
- **Explicit Homology / Closeness Boundaries**: Sequence identity and Pfam domain overlap between external viral proteins and the development set must be quantified.

---

## 3. Potential External Validation Candidate Systems

| Candidate Viral System | Biological Relevance | Technical Feasibility | Current Limitation |
| :--- | :--- | :---: | :--- |
| **Varicella-Zoster Virus (VZV / HHV-3)** | Human neurotropic alphaherpesvirus; shares 60+ core orthologs with HSV. | Medium | Lower direct kinetic literature density; requires independent manual adjudication. |
| **Pseudorabies Virus (PRV)** | Veterinary alphaherpesvirus model; extensively sequenced with time-course data. | High | Significant sequence divergence; requires domain-level alignment. |
| **Prospective HSV Clinical Isolates** | Directly assesses within-species evolutionary robustness and variant drift. | High (Lab dependent) | Requires new mass spectrometry profiling under CHX and PAA metabolic arrest. |

---

## 4. Scientific Verdict on External Generalization

At present, claims of model generalization are **strictly supported within scope** to unseen gene families within HSV-1 and HSV-2. Broad claims that the model generalizes across all herpesviruses or discovers universal viral replication timing mechanisms are scientifically unsupported and premature.
