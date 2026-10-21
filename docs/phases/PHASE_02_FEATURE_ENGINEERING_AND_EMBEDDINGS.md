# Phase 2: Feature Engineering & Protein Language Model Embeddings

**Consolidated Step**: Step 5 (Frozen Dataset & ESM-2 Embedding Generation)  
**Primary Deliverables**: `data/processed/esm2_t12_representatives_embeddings.npy`, `data/processed/esm2_t12_representatives_metadata.csv`  

---

## 1. Scientific Objective
Extract multi-scale representations for all 3,111 representative viral proteins: classical bulk physicochemical descriptors and deep contextual transformer embeddings from the Evolutionary Scale Modeling (ESM-2) protein language model.

## 2. Inputs
- `data/processed/valid_proteins_combined.fasta` (Filtered representative seeds)
- Pre-trained model weights: `esm2_t12_35M_UR50D` (12 layers, 480 embedding dimensions, 35M parameters)

## 3. Methodology & Feature Extraction
1. **Bulk Physicochemical Descriptors ($D=8$)**:
   - Sequence length (aa), molecular weight (Da), isoelectric point (pI), aliphatic index, instability index, GRAVY hydropathicity, charged residue fraction, aromatic residue fraction.
2. **ESM-2 Embeddings ($D=480$)**:
   - Computed on GPU/CPU using fairseq / ESM PyTorch library.
   - Extracted residue representations from layer 12 and applied mean-pooling across the length of each protein sequence.
   - Output stored as a dense floating-point matrix of shape `(3111, 480)`.
3. **Combined Representation ($D=488$)**:
   - Concatenation of the 8 physical descriptors with the 480 ESM-2 dimensions.

## 4. Key Outputs & Provenance
- `data/processed/esm2_t12_representatives_embeddings.npy`: Binary NumPy array `(3111, 480)` (5,973,248 bytes).
- `data/processed/esm2_t12_representatives_metadata.csv`: Complete metadata joining FASTA headers, accessions, lengths, and embedding row indices.

## 5. Validation & Test Outcomes
- Validated under `tests/test_05_frozen_dataset_and_embeddings.py` (10 passing tests).
- Row-index closure confirmed: `row_index` in metadata strictly maps 1-to-1 to row index of `.npy` array without off-by-one errors or NaN entries.
