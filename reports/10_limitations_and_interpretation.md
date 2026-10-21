# Step 10: Limitations, Methodological Safeguards, and Scientific Interpretation

**Project**: *Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes*  
**Date**: 2026-10-21  
**Deliverable**: [`reports/10_limitations_and_interpretation.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/10_limitations_and_interpretation.md)  

---

## 1. What Step 10 Establishes (Scientific Scope)

1. **Language Models Capture Viral Replication Timing Information**:
   - Linear evaluation of pre-trained ESM-2 representations achieves **~80% balanced accuracy** across held-out HSV gene families without fine-tuning.
   - This significantly surpasses majority guessing (~50%) and bulk physicochemical descriptors (~51-57%), proving that evolutionary embeddings encode constraints correlated with pre- vs. post-replicative lifecycle roles.
2. **Leakage-Resistant Generalization**:
   - Because cross-validation strictly held out entire homology groups, the performance represents genuine generalization to *unseen viral gene families* within the HSV proteome.

## 2. What Step 10 Does NOT Claim (Explicit Boundary Constraints)

In accordance with scientific integrity guidelines, the following boundaries must be maintained:

1. **No Generalization to Novel Herpesvirus Genera**:
   - The evaluation was conducted exclusively on Alphaherpesvirinae (*Simplexvirus*: HSV-1 and HSV-2). It does **NOT** prove that this linear classifier can predict replication timing in Betaherpesvirinae (e.g., CMV) or Gammaherpesvirinae (e.g., EBV, KSHV), which possess diverged temporal cascades and lineage-specific genes.
2. **No Causal Mechanism Inferred**:
   - Predictive association does not imply that ESM-2 has discovered the mechanistic cause of transcription timing. Viral transcription timing is primarily regulated by genomic promoter architectures, viral transactivator binding motifs (e.g., VP16/Oct-1 response elements for IE genes), host RNA Polymerase II recruitment, and chromatin accessibility, rather than intrinsic protein coding sequence alone.
3. **Sample Size Limitations ($N=63$)**:
   - Although 687 sequence records were analyzed, they correspond to only **63 independent evolutionary loci**. With $N=63$, the 95% bootstrap confidence interval for ESM-2 Macro-F1 spans **[0.6865, 0.8931]**. The true population generalization error remains subject to this fundamental sample size constraint.
4. **Transcript vs. Protein Kinetic Divergence**:
   - While 52 families have concordance between protein and nascent transcript kinetics, 11 families rely solely on 4sU-seq RNA data. In viruses, translation efficiency and post-translational stability can decouple protein accumulation from mRNA abundance.

## 3. Preserved Artifacts and Next Steps

All raw FASTA files, frozen embeddings (`esm2_t12_representatives_embeddings.npy`), and historical Step 1-9.1 reports remain unmodified.
The complete results are serialized in:
- [`reports/10_baseline_metrics.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/10_baseline_metrics.csv)
- [`reports/10_fold_level_metrics.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/10_fold_level_metrics.csv)
- [`reports/10_out_of_fold_predictions.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/10_out_of_fold_predictions.csv)
- [`reports/10_family_level_metrics.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/10_family_level_metrics.csv)
- [`reports/10_model_configuration.json`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/10_model_configuration.json)
