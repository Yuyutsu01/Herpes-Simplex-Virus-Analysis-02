# Step 8 Literature Review: HSV-1 & HSV-2 Temporal Expression Kinetics

**Audit Date**: 2026-10-21 13:32:33 UTC  
**Auditor**: Independent Computational Biologist  
**Deliverables**: `data/processed/08_temporal_evidence_sources.csv`, `data/processed/08_temporal_gene_evidence.csv`  

---

## 1. Executive Summary & Experimental Definitions

Herpes simplex virus gene expression operates via an orderly, tightly regulated cascade originally elucidated by Honess & Roizman (1974, 1975) and validated by modern high-throughput transcriptomics (Rutkowski et al. 2015) and quantitative mass spectrometry (Kulej et al. 2017).

### Operational Definitions of Temporal Kinetic Classes:

1. **Immediate-Early (IE / $\alpha$)**:
   - **Operational Hallmark**: Transcribed in the absence of de novo viral protein synthesis.
   - **Experimental Definition**: Resistant to translation inhibition by cycloheximide (CHX, 50–100 $\mu$g/ml). Upon removal of CHX in the presence of actinomycin D, $\alpha$ mRNAs are selectively translated.
   - **Regulatory Mechanism**: Promoters contain TAATGARAT consensus motifs recognized by host Oct-1 and the virion tegument transactivator VP16 (UL48).
   - **Canonical Core Genes**: **`RL2` (ICP0), `RS1` (ICP4), `US1` (ICP22), `UL54` (ICP27), `US12` (ICP47)**.

2. **Early (E / $\beta$)**:
   - **Operational Hallmark**: Transcribed prior to viral DNA replication; strictly dependent on functional $\alpha$ proteins (primarily ICP4 and ICP27).
   - **Experimental Definition**: Inhibited by CHX pre-treatment (sensitive to translation block), but **resistant to viral DNA synthesis inhibitors** (phosphonoacetic acid [PAA, 100–300 $\mu$g/ml], acyclovir [ACV], or phosphonoformic acid [PFA]).
   - **Biological Roles**: Enzymes and factors required for nucleotide metabolism and viral DNA synthesis.
   - **Canonical Core Genes**: **`UL23` (TK), `UL29` (ICP8), `UL30` (DNA Pol), `UL42` (processivity), `UL5`, `UL8`, `UL52` (helicase-primase), `UL9` (OBP), `UL39`, `UL40` (RNR)**.

3. **Leaky-Late (L1 / $\gamma_1$)**:
   - **Operational Hallmark**: Initiated early before DNA replication, but expression increases substantially (5- to 20-fold) following viral DNA synthesis.
   - **Experimental Definition**: Partially attenuated but **detectable** in the presence of PAA or ACV.
   - **Canonical Core Genes**: **`UL19` (major capsid VP5), `UL27` (gB), `UL48` (VP16), `US6` (gD), `UL26` / `UL26.5` (capsid scaffold)**.

4. **True-Late (L2 / $\gamma_2$)**:
   - **Operational Hallmark**: Strictly dependent on viral DNA replication; promoters require ongoing replication fork progression and TATT-like non-canonical TATA boxes.
   - **Experimental Definition**: **Completely blocked** by PAA, ACV, or Ara-C.
   - **Canonical Core Genes**: **`UL36` (large tegument), `UL38` (VP19C), `UL44` (gC), `US4` (gG), `UL11`, `US11`**.

5. **Mixed / Context-Dependent / Insufficient Evidence**:
   - Genes exhibiting cell-type-dependent kinetics (e.g. `RL1`/ICP34.5), overlapping reading frames with discordant promoter usage, or genes lacking clear metabolic arrest documentation.

---

## 2. Primary Evidence Sources Table

| Source ID | Study Citation | Assays / Treatments | Target System | Key Finding |
| :--- | :--- | :--- | :--- | :--- |
| **SRC_HONESS_1974** | Honess & Roizman, *J Virol* 14:8–19 | $^{35}\text{S}$-Met, CHX + ActD | HSV-1 (F) / HEp-2 | Foundational definition of $\alpha$ cascade |
| **SRC_HONESS_1975** | Honess & Roizman, *J Virol* 16:1308–1326 | PAA, Ara-C, SDS-PAGE | HSV-1 (F) / HEp-2 | Established $\beta$ (PAA-res) vs $\gamma$ (PAA-sens) |
| **SRC_RUTKOWSKI_2015** | Rutkowski et al., *PLoS Pathog* 11:e1005245 | 4sU-seq time-course, CHX, PAA | HSV-1 (17) / primary HFF | Genome-wide quantitative transcription kinetics |
| **SRC_KULEJ_2017** | Kulej et al., *J Proteome Res* 16:922–939 | SILAC mass spectrometry | HSV-1 (17) / HeLa | Time-resolved protein accumulation kinetics |
| **SRC_DREMEL_2019** | Dremel & DeLuca, *J Virol* 93:e00452-19 | ChIP-seq RNA Pol II, mutants | HSV-1 wt & d109 | Confirmed ICP4/27 recruitment to $\beta$ and $\gamma$ |
| **SRC_ROIZMAN_KNIPE_2013**| Fields Virology 6th ed., Ch 39 | Multi-decade consensus | HSV-1 & HSV-2 | Standard reference catalog across 74 genes |
