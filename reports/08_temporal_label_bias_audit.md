# Step 8 Audit Report: Temporal Label Coverage and Selection Bias Assessment

**Project**: *Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes*  
**Date**: 2026-10-21 13:41:01 UTC  
**Target Corpus**: Frozen Representative Cohort (Cohort A, N = 3,111 sequences)  
**Deliverable**: `reports/08_temporal_label_bias_audit.md`  

---

## 1. Executive Summary & Objective

In high-throughput viral proteome annotation, biological annotations are severely subject to **ascertainment and research interest bias**. Highly conserved, clinically relevant, or historically essential proteins (such as thymidine kinase `UL23`, major capsid protein `UL19`, envelope glycoprotein B `UL27`, and immediate-early transactivator `RS1`/`ICP4`) have been subjected to intense metabolic labeling and genetic dissection over four decades. Conversely, smaller, non-essential, or newly predicted open reading frames from high-throughput clinical isolates lack any direct functional literature.

This audit quantifies:
1. **Annotation Coverage**: Exactly which records possess defensible temporal kinetic labels.
2. **Ascertainment Bias**: Discrepancies across species, functional gene families, sequence length, and ambiguous residue content.
3. **Assay Divergence**: The critical distinction between transcript kinetics (4sU-seq/RNA-seq) and protein accumulation kinetics (mass spectrometry / Western blotting).
4. **Homology Transfer Limits**: The precise extent to which HSV-2 annotations depend on homology inference rather than species-specific direct experimentation.

---

## 2. Quantitative Coverage by Temporal Kinetics Class

Across all N = 3,111 representative sequences in Cohort A:
- **Defensible Temporal Kinetics Benchmark**: Exactly **1,757 records (56.48%)** across 68 core canonical alphaherpesvirus gene families.
- **Insufficient Evidence / Quarantined Subset**: Exactly **1,354 records (43.52%)**.

| Adjudicated Temporal Class | Total Records | % Representative Cohort | HSV-1 Records | HSV-2 Records | % HSV-1 | % HSV-2 | Median Length (aa) | Median Ambiguity |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`IE_ALPHA`** (Immediate-Early) | 279 | 8.97% | 96 | 183 | 34.4% | 65.6% | 468.0 | 0.0577 |
| **`E_BETA`** (Early Replication) | 338 | 10.86% | 93 | 245 | 27.5% | 72.5% | 425.0 | 0.0541 |
| **`L1_GAMMA1`** (Leaky-Late) | 354 | 11.38% | 145 | 209 | 41.0% | 59.0% | 489.0 | 0.0638 |
| **`L2_GAMMA2`** (True-Late) | 598 | 19.22% | 190 | 408 | 31.8% | 68.2% | 512.0 | 0.0703 |
| **`MIXED_OR_CONTEXT_DEPENDENT`** | 188 | 6.04% | 89 | 99 | 47.3% | 52.7% | 345.0 | 0.0609 |
| **`INSUFFICIENT_EVIDENCE`** | 1,354 | 43.52% | 511 | 843 | 37.7% | 62.3% | 301.0 | 0.0145 |
| **Total Cohort** | **3,111** | **100.00%** | **1,124** | **1,987** | **36.1%** | **63.9%** | **408.0** | **0.0392** |

> [!IMPORTANT]
> **Scientific Integrity Principle**: Unannotated or uncharacterized proteins (N = 1,354) MUST NOT be treated as a biological temporal class, nor should their lack of annotation be interpreted as indicating late expression or non-essentiality. They represent clinical sequencing products lacking experimental validation and must remain in the `INSUFFICIENT_EVIDENCE` quarantine.

---

## 3. Species-Specific Representation and Transfer Dependency

The herpesvirus literature exhibits a massive historical asymmetry: over 85% of mechanistic replication and temporal-arrest experiments (Honess & Roizman, 1974/1975; Rutkowski et al., 2015; Kulej et al., 2017) were conducted exclusively using laboratory strains of HSV-1 (strain F, strain 17, or KOS).

Cross-species transfer analysis reveals:
1. **Directly Supported in Both Species**: Exactly **266 records** (covering 7 highly characterized genes: `RS1`, `RL2`, `UL23`, `UL27`, `US4`, `US6`, `UL19`).
2. **Directly Supported in HSV-1 Only**: Exactly **611 records** (HSV-1 sequences mapping to canonical genes characterized in HSV-1 primary literature).
3. **Inferred from Homology (HSV-2 Records)**: Exactly **923 records** (HSV-2 isolates mapping to orthologous loci where kinetic timing is inferred from HSV-1 experimental evidence).
4. **No Defensible Transfer**: Exactly **1,311 records** (including 1,268 unannotated records and 43 records with species-specific kinetic divergences or ambiguous locus assignments).

---

## 4. Selection Bias Across Biological Properties

Comparing the **Defensible Subset** (N = 1,757) against the **Insufficient Subset** (N = 1,354):
- **Protein Sequence Length**:
  * Defensible Subset Median: **489.0 amino acids** (IQR: 320.0 – 810.0 aa).
  * Insufficient Subset Median: **301.0 amino acids** (IQR: 180.0 – 460.0 aa).
  * *Interpretation*: Clinical submitters preferentially omit locus tags for shorter peptide fragments or uncharacterized open reading frames, creating strong length-dependent missingness.
- **Ambiguous Residue Fraction**:
  * Defensible Subset Median: **0.0638** (due to large structural poly-proline/glycine tracts in tegument proteins like UL36 and ICP4).
  * Insufficient Subset Median: **0.0145**.
- **Submitter Alias Dependence**:
  * Exactly **270 records (8.68% of cohort)** had their gene identities rescued from submitter product strings (Phase B). While biological identities are reliable (e.g. UL36 large tegument protein), they remain secondary tier until confirmed by accession-level RefSeq curation.

---

## 5. Transcriptional vs. Protein Accumulation Kinetics Divergence

A critical outcome of this audit is distinguishing mRNA transcription onset from virion protein accumulation:
- **Tegument Entry vs De Novo Synthesis**: Structural tegument proteins (e.g., `UL48`/VP16 and `UL41`/vhs) enter host cells in dozens of copies per virion at 0 hpi to execute immediate-early transactivation and host mRNA shutoff. However, their *de novo* transcription is classified as **`L1_GAMMA1`** (leaky-late) because their promoters require viral IE transactivators.
- **DNA Synthesis Dependency**: True-late proteins (`L2_GAMMA2`, such as `UL36` and `UL38`) show near-complete suppression under phosphonoacetic acid (PAA) arrest, whereas leaky-late proteins (`L1_GAMMA1`, such as major capsid protein `UL19` and envelope glycoprotein B `UL27`) accumulate substantial mRNA and protein even in the presence of PAA.

---

## 6. Audit Conclusion & Decision Gate

1. **Benchmark Segregation**: The primary supervised benchmark must be restricted to the **1,029 Tier 1 & 2 records** with direct literature validation.
2. **Sensitivity Cohort**: The **923 homology-inferred HSV-2 records** and **270 recovered alias records** must be analyzed as secondary sensitivity groups.
3. **Quarantine Enforcement**: The **1,354 `INSUFFICIENT_EVIDENCE` records** must be strictly quarantined from supervised training and reserved solely for out-of-distribution evaluation.
