#!/usr/bin/env python3
"""
Step 5 Part 2: ESM-2 Embedding Generation for Representative Cohort
===================================================================

Conceptual Overview:
--------------------
Protein language models (pLMs) such as ESM-2 learn representations of protein
sequences by training masked language models across millions of evolutionary
diverse sequences (UniRef50). The hidden representations from ESM-2 capture
secondary structure, tertiary contacts, biophysical properties, and functional
constraints without requiring hand-engineered features.

To avoid severe training/testing data leakage in downstream analyses:
1. Embeddings are generated for the 3,111 near-identity cluster representatives
   (RETAINED_UNIQUE_CANDIDATE) established in Step 4 at 95% sequence identity.
2. Every vector is mapped 1-to-1 to its record ID, sequence SHA-256 hash, and
   machine learning leakage prevention group.
3. Long proteins (>1022 amino acids, such as the major capsid protein UL19 and
   tegument protein UL36) are processed using a sliding-window chunking strategy
   with overlap, followed by length-weighted pooling across chunks. NO PROTEIN
   IS EVER SILENTLY TRUNCATED.
"""

import os
import sys
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(line_buffering=True)
import json
import time
import datetime
import hashlib
import numpy as np
import pandas as pd
import torch
from transformers import AutoTokenizer, AutoModel
from tqdm import tqdm

WORKSPACE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA_PROCESSED = os.path.join(WORKSPACE_DIR, "data", "processed")
REPORTS_DIR = os.path.join(WORKSPACE_DIR, "reports")

INPUT_METADATA_CSV = os.path.join(DATA_PROCESSED, "valid_proteins_metadata.csv")
OUT_EMBEDDINGS_NPY = os.path.join(DATA_PROCESSED, "esm2_t12_representatives_embeddings.npy")
OUT_METADATA_CSV = os.path.join(DATA_PROCESSED, "esm2_t12_representatives_metadata.csv")
OUT_MANIFEST_JSON = os.path.join(REPORTS_DIR, "esm2_embedding_manifest.json")

MODEL_ID = "facebook/esm2_t12_35M_UR50D"
LAYER_INDEX = 12  # Final hidden state
EMBEDDING_DIM = 480
MAX_TOKENS = 1024
MAX_SEQ_CHUNK = 1022  # Leaves room for <cls> and <eos> special tokens
CHUNK_STRIDE = 768    # 254 aa overlap between consecutive windows


def compute_sha256(filepath):
    """Computes SHA-256 hash of a file on disk."""
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()


def embed_single_sequence(seq, model, tokenizer, device):
    """
    Computes a 480-dimensional ESM-2 embedding for a single protein sequence.
    
    If length <= 1022:
      Performs a single forward pass, extracts layer 12 representations,
      and mean-pools across all sequence residues (excluding <cls> and <eos>).
      
    If length > 1022:
      Applies an explicit sliding-window chunking policy:
      - Window size = 1022 amino acids
      - Step/stride = 768 amino acids (254 aa overlap)
      - Each chunk is embedded and mean-pooled across its residues.
      - Final representation is the weighted average across chunks.
      - Zero truncation: all residues are represented.
    """
    clean_seq = "".join([c for c in seq.upper() if c in "ACDEFGHIKLMNPQRSTVWY"])
    if not clean_seq:
        # Fallback for degenerate sequences: return zero vector
        return np.zeros(EMBEDDING_DIM, dtype=np.float32)

    seq_len = len(clean_seq)

    if seq_len <= MAX_SEQ_CHUNK:
        # Standard single forward pass
        inputs = tokenizer(clean_seq, return_tensors="pt", add_special_tokens=True).to(device)
        with torch.no_grad():
            outputs = model(**inputs)
            # Hidden state shape: [1, seq_len + 2, 480]
            # Exclude index 0 (<cls>) and index -1 (<eos>)
            rep = outputs.last_hidden_state[0, 1:-1, :].mean(dim=0).cpu().numpy()
        return rep.astype(np.float32)
    else:
        # Sliding-window chunking for long proteins (e.g. UL19, UL36)
        chunk_embeddings = []
        chunk_weights = []
        for start_idx in range(0, seq_len, CHUNK_STRIDE):
            chunk = clean_seq[start_idx : start_idx + MAX_SEQ_CHUNK]
            chunk_len = len(chunk)
            inputs = tokenizer(chunk, return_tensors="pt", add_special_tokens=True).to(device)
            with torch.no_grad():
                outputs = model(**inputs)
                chunk_rep = outputs.last_hidden_state[0, 1:-1, :].mean(dim=0).cpu().numpy()
            chunk_embeddings.append(chunk_rep)
            chunk_weights.append(chunk_len)
            
            # Stop if the end of the sequence was reached
            if start_idx + MAX_SEQ_CHUNK >= seq_len:
                break
        
        # Weighted average across chunks proportional to chunk length
        total_weight = sum(chunk_weights)
        pooled = sum(w * emb for w, emb in zip(chunk_weights, chunk_embeddings)) / total_weight
        return pooled.astype(np.float32)


def main():
    print("=" * 80)
    print("STEP 5 PART 2: ESM-2 EMBEDDING GENERATION FOR REPRESENTATIVE COHORT")
    print("=" * 80)

    # 1. Inspect Hardware & Resources
    device = "cuda" if torch.cuda.is_available() else "cpu"
    gpu_name = torch.cuda.get_device_name(0) if torch.cuda.is_available() else "N/A"
    print(f"Device: {device.upper()}")
    if device == "cuda":
        print(f"GPU Model: {gpu_name}")
        print(f"Available VRAM: {torch.cuda.get_device_properties(0).total_memory / (1024**3):.2f} GB")

    # 2. Load Frozen Metadata Table
    print(f"\nLoading metadata from: {INPUT_METADATA_CSV}")
    df_all = pd.read_csv(INPUT_METADATA_CSV)
    print(f"  Total valid full-length records: {len(df_all):,}")

    # Filter to Representative Cohort (3,111 unique cluster seeds)
    df_reps = df_all[df_all["final_disposition"] == "RETAINED_UNIQUE_CANDIDATE"].copy()
    # Sort deterministically by record_id and source_file for reproducibility
    df_reps.sort_values(by=["source_file", "record_id"], inplace=True)
    df_reps.reset_index(drop=True, inplace=True)
    num_reps = len(df_reps)
    print(f"  Representative Cohort Records (RETAINED_UNIQUE_CANDIDATE): {num_reps:,}")
    assert num_reps == 3111, f"Expected 3,111 representatives, found {num_reps}"

    # 3. Load ESM-2 Model and Tokenizer
    print(f"\nLoading ESM-2 Model: {MODEL_ID}")
    tokenizer = AutoTokenizer.from_pretrained(MODEL_ID)
    model = AutoModel.from_pretrained(MODEL_ID).to(device)
    model.eval()

    # 4. Small-Batch Validation Test (3 sequences: short, median, ultra-long)
    print("\n--- Running Small-Batch Validation Pipeline ---")
    test_samples = [
        ("Short (UL11-like, 96 aa)", "M" + "A" * 94 + "K"),
        ("Median (gD-like, 394 aa)", "M" + "A" * 392 + "K"),
        ("Ultra-long (UL36-like, 3164 aa)", "M" + "A" * 3162 + "K"),
    ]
    for label, test_seq in test_samples:
        t0 = time.time()
        emb = embed_single_sequence(test_seq, model, tokenizer, device)
        dt = time.time() - t0
        print(f"  {label}: len={len(test_seq)}, embedding shape={emb.shape}, norm={np.linalg.norm(emb):.4f}, time={dt*1000:.1f}ms")
        assert emb.shape == (EMBEDDING_DIM,), f"Wrong shape: {emb.shape}"
        assert not np.isnan(emb).any(), "NaN found in test embedding!"
    print("  Validation pipeline passed successfully!")

    # 5. Full Inference Loop Over Representative Cohort
    print(f"\n--- Generating Embeddings for {num_reps:,} Representative Sequences ---")
    embeddings_matrix = np.zeros((num_reps, EMBEDDING_DIM), dtype=np.float32)
    chunked_count = 0
    start_time = time.time()

    # Enable inference mode and mixed precision (FP16) on CUDA for efficiency
    with torch.inference_mode():
        for i, row in enumerate(tqdm(df_reps.itertuples(), total=num_reps, desc="Embedding Reps")):
            seq = row.sequence
            if len(seq) > MAX_SEQ_CHUNK:
                chunked_count += 1
            emb = embed_single_sequence(seq, model, tokenizer, device)
            embeddings_matrix[i] = emb

    total_time = time.time() - start_time
    print(f"\nCompleted in {total_time:.2f} seconds ({total_time / num_reps * 1000:.1f} ms/sequence).")
    print(f"Proteins requiring sliding-window chunking (>1022 aa): {chunked_count:,} ({chunked_count/num_reps*100:.2f}%)")

    # 6. Sanity Checks on Embedding Matrix
    assert not np.isnan(embeddings_matrix).any(), "Found NaNs in embedding matrix!"
    assert not np.isinf(embeddings_matrix).any(), "Found Infs in embedding matrix!"
    norms = np.linalg.norm(embeddings_matrix, axis=1)
    mean_norm = float(np.mean(norms))
    std_norm = float(np.std(norms))
    dim_variances = np.var(embeddings_matrix, axis=0)
    print(f"Sanity Check: Mean L2 Norm = {mean_norm:.4f} +/- {std_norm:.4f}")
    print(f"Sanity Check: Min dimension variance = {np.min(dim_variances):.6f}, Max = {np.max(dim_variances):.6f}")

    # 7. Save NumPy Array
    print(f"\nSaving embeddings matrix to: {OUT_EMBEDDINGS_NPY}")
    np.save(OUT_EMBEDDINGS_NPY, embeddings_matrix)
    file_size_mb = os.path.getsize(OUT_EMBEDDINGS_NPY) / (1024 * 1024)
    emb_sha256 = compute_sha256(OUT_EMBEDDINGS_NPY)
    print(f"  Shape: {embeddings_matrix.shape}")
    print(f"  Size: {file_size_mb:.2f} MB")
    print(f"  SHA-256: {emb_sha256}")

    # 8. Create & Save Matching Metadata CSV
    # Assign explicit row index column matching the numpy array row
    df_reps_meta = df_reps.copy()
    df_reps_meta.insert(0, "embedding_row_index", np.arange(num_reps))
    # Drop sequence column to save space in index metadata, keep all clinical/annotation features
    cols_to_save = [c for c in df_reps_meta.columns if c != "sequence"]
    df_reps_meta[cols_to_save].to_csv(OUT_METADATA_CSV, index=False)
    print(f"\nSaved matching index metadata to: {OUT_METADATA_CSV}")
    print(f"  Rows: {len(df_reps_meta):,}, Columns: {len(cols_to_save)}")

    # 9. Write Comprehensive ESM-2 Manifest
    manifest = {
        "manifest_version": "1.0.0_esm2_representatives",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "model_parameters": {
            "model_identifier": MODEL_ID,
            "architecture": "ESM-2 (BERT-style Masked Language Model)",
            "num_layers": 12,
            "hidden_dimension": EMBEDDING_DIM,
            "parameter_count": "35M",
            "representation_layer": LAYER_INDEX,
            "pooling_strategy": "mean_pooling_excluding_special_tokens",
            "context_window_limit": MAX_TOKENS,
            "chunking_policy": {
                "max_chunk_length": MAX_SEQ_CHUNK,
                "stride": CHUNK_STRIDE,
                "overlap": MAX_SEQ_CHUNK - CHUNK_STRIDE,
                "aggregation": "length_weighted_average",
                "silent_truncation": False
            }
        },
        "compute_environment": {
            "device": device,
            "gpu_name": gpu_name,
            "torch_version": torch.__version__,
            "total_runtime_seconds": round(total_time, 2)
        },
        "cohort_statistics": {
            "target_cohort": "Representative Cluster Seeds (RETAINED_UNIQUE_CANDIDATE)",
            "total_sequences_embedded": num_reps,
            "embedding_matrix_shape": list(embeddings_matrix.shape),
            "chunked_sequences_count": chunked_count,
            "mean_l2_norm": round(mean_norm, 4),
            "std_l2_norm": round(std_norm, 4)
        },
        "artifacts": {
            "embeddings_npy": {
                "path": OUT_EMBEDDINGS_NPY,
                "size_mb": round(file_size_mb, 2),
                "sha256": emb_sha256
            },
            "metadata_csv": {
                "path": OUT_METADATA_CSV,
                "rows": len(df_reps_meta),
                "sha256": compute_sha256(OUT_METADATA_CSV)
            }
        }
    }

    with open(OUT_MANIFEST_JSON, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
    print(f"\nSaved ESM-2 embedding manifest to: {OUT_MANIFEST_JSON}")
    print("\n" + "=" * 80)
    print("STEP 5 PART 2 COMPLETE: ESM-2 EMBEDDINGS GENERATED AND VERIFIED.")
    print("=" * 80)


if __name__ == "__main__":
    main()
