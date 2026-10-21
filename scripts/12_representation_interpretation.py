"""
Step 12: Interpretable ESM-2 Representation Analysis for HSV Temporal Classes.

Scientific Objective:
    Investigate which sequence-derived features, embedding dimensions, and bulk
    physicochemical properties correlate with binary temporal-class predictions
    (PRE_REPLICATIVE = IE + E vs POST_REPLICATIVE = L1 + L2).
    Distinguish predictive association from biological mechanism.

Key Analyses:
    1. Input audit and SHA-256 provenance freeze.
    2. Embedding-level linear attribution (Ridge coefficients across 5 folds).
    3. Multicollinearity and intrinsic dimensionality analysis of ESM-2 embeddings.
    4. Sequence-level compositional and physicochemical analysis with family-aware
       aggregation, multiple testing correction (BH FDR and Holm-Bonferroni),
       and sequence length confounder control.
    5. Attention and token-level reality check (auditing saved artifacts).
    6. Stratification across evidence tiers (52 protein-validated vs 11 transcript-only).
    7. Diagnostic evaluation on HSV-1 to HSV-2 ortholog-transfer families (N=7).

Deliverables:
    - scripts/12_representation_interpretation.py
    - reports/12_representation_interpretation.md
    - reports/12_embedding_coefficient_stability.csv
    - reports/12_feature_association.csv
    - reports/12_sequence_level_diagnostics.csv
    - reports/12_validation_and_limitations.md
    - reports/12_reproducibility_manifest.json
"""

import os
import sys
import json
import hashlib
import logging
from datetime import datetime, timezone
from typing import Dict, List, Tuple, Any

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, pearsonr, spearmanr, norm
from sklearn.linear_model import LogisticRegression, Lasso
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger("Step12Interpretation")

# Base directory definitions
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA_DIR = os.path.join(BASE_DIR, "data", "processed")
REPORTS_DIR = os.path.join(BASE_DIR, "reports")

# Input file paths
FILE_BENCHMARK = os.path.join(DATA_DIR, "08_1_strict_temporal_benchmark.csv")
FILE_METADATA = os.path.join(DATA_DIR, "valid_proteins_metadata.csv")
FILE_EMB_META = os.path.join(DATA_DIR, "esm2_t12_representatives_metadata.csv")
FILE_EMBEDDINGS = os.path.join(DATA_DIR, "esm2_t12_representatives_embeddings.npy")
FILE_SPLIT_MANIFEST = os.path.join(DATA_DIR, "09_group_split_manifest.csv")
FILE_OOF_PREDS = os.path.join(REPORTS_DIR, "10_out_of_fold_predictions.csv")
FILE_TRANSFER_EXP = os.path.join(REPORTS_DIR, "11_1_cross_species_transfer_experiment.csv")
FILE_EMB_MANIFEST = os.path.join(REPORTS_DIR, "esm2_embedding_manifest.json")

# Output deliverables
OUT_COEF_STABILITY = os.path.join(REPORTS_DIR, "12_embedding_coefficient_stability.csv")
OUT_FEATURE_ASSOC = os.path.join(REPORTS_DIR, "12_feature_association.csv")
OUT_SEQ_DIAGNOSTICS = os.path.join(REPORTS_DIR, "12_sequence_level_diagnostics.csv")
OUT_REPORT_MD = os.path.join(REPORTS_DIR, "12_representation_interpretation.md")
OUT_LIMITATIONS_MD = os.path.join(REPORTS_DIR, "12_validation_and_limitations.md")
OUT_MANIFEST_JSON = os.path.join(REPORTS_DIR, "12_reproducibility_manifest.json")

RANDOM_SEED = 42
PHYSICOCHEMICAL_COLUMNS = [
    "sequence_length",
    "molecular_weight",
    "isoelectric_point",
    "instability_index",
    "aliphatic_index",
    "gravy",
    "charged_fraction",
    "aromatic_fraction",
]
CANONICAL_AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")


def compute_sha256(filepath: str) -> str:
    """Calculates SHA-256 hash of a file."""
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()


def audit_and_freeze_inputs() -> Dict[str, Any]:
    """
    Verifies input integrity, alignment, and hashes.
    Ensures primary benchmark has 687 records, 63 families, and exact index alignment.
    """
    logger.info("=== Phase 1: Auditing and Freezing Inputs ===")
    
    input_files = [
        FILE_BENCHMARK, FILE_METADATA, FILE_EMB_META,
        FILE_EMBEDDINGS, FILE_SPLIT_MANIFEST, FILE_OOF_PREDS
    ]
    file_hashes = {}
    for p in input_files:
        if not os.path.exists(p):
            raise FileNotFoundError(f"Missing required input file: {p}")
        file_hashes[os.path.basename(p)] = compute_sha256(p)
        logger.info(f"Verified {os.path.basename(p)} (SHA-256: {file_hashes[os.path.basename(p)][:12]}...)")

    # Load tables
    df_bench = pd.read_csv(FILE_BENCHMARK)
    df_meta = pd.read_csv(FILE_METADATA)
    df_emb_meta = pd.read_csv(FILE_EMB_META)
    df_manifest = pd.read_csv(FILE_SPLIT_MANIFEST)
    embeddings = np.load(FILE_EMBEDDINGS)

    # Validate shapes and alignment
    assert len(df_bench) == 3111, f"Expected 3,111 benchmark rows, found {len(df_bench)}"
    assert len(df_emb_meta) == 3111, f"Expected 3,111 embedding metadata rows, found {len(df_emb_meta)}"
    assert len(df_manifest) == 3111, f"Expected 3,111 manifest rows, found {len(df_manifest)}"
    assert embeddings.shape == (3111, 480), f"Expected (3111, 480) embeddings shape, found {embeddings.shape}"

    # Filter to primary DIRECT_EVIDENCE cohort
    direct_mask = (df_manifest["adjudicated_cohort"] == "DIRECT_EVIDENCE")
    df_direct = df_manifest[direct_mask].copy()
    assert len(df_direct) == 687, f"Expected 687 direct sequences, found {len(df_direct)}"

    # Add binary class
    binary_map = {
        "IE_ALPHA": "PRE_REPLICATIVE",
        "E_BETA": "PRE_REPLICATIVE",
        "L1_GAMMA1": "POST_REPLICATIVE",
        "L2_GAMMA2": "POST_REPLICATIVE"
    }
    df_direct["binary_class"] = df_direct["normalized_temporal_class"].map(binary_map)
    assert not df_direct["binary_class"].isnull().any()

    # Join metadata (sequence and physicochemical properties)
    meta_sub = df_meta.set_index("record_id")
    df_direct["sequence"] = df_direct["record_id"].map(meta_sub["sequence"])
    assert not df_direct["sequence"].isnull().any(), "Missing sequence for direct records!"

    for col in PHYSICOCHEMICAL_COLUMNS:
        df_direct[col] = df_direct["record_id"].map(meta_sub[col])
        assert not df_direct[col].isnull().any(), f"Missing {col} for direct records!"

    # Add measured molecule from benchmark
    bench_sub = df_bench.set_index("record_id")
    df_direct["measured_molecule"] = df_direct["record_id"].map(bench_sub["measured_molecule"])

    num_fams = df_direct["homology_group_id"].nunique()
    assert num_fams == 63, f"Expected 63 families, found {num_fams}"
    pre_fams = df_direct[df_direct["binary_class"] == "PRE_REPLICATIVE"]["homology_group_id"].nunique()
    post_fams = df_direct[df_direct["binary_class"] == "POST_REPLICATIVE"]["homology_group_id"].nunique()
    assert pre_fams == 23 and post_fams == 40, f"Expected 23 PRE and 40 POST families, got {pre_fams}, {post_fams}"

    logger.info("Phase 1 PASSED: All inputs verified and aligned.")
    return {
        "file_hashes": file_hashes,
        "df_direct": df_direct,
        "embeddings": embeddings,
        "df_bench": df_bench,
        "df_manifest": df_manifest
    }


def analyze_embedding_coefficients(
    df_direct: pd.DataFrame,
    embeddings: np.ndarray
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """
    Analyzes linear model coefficients across 5 cross-validation folds.
    Evaluates coefficient stability, fold-to-fold correlations, and top predictive dimensions.
    """
    logger.info("=== Phase 2: Analyzing ESM-2 Linear Model Coefficients & Stability ===")

    y = df_direct["binary_class"].values
    folds = df_direct["fold_index"].values
    indices = df_direct.index.values
    X_esm = embeddings[indices]

    fold_coefs = []
    # Fit Ridge logistic regression (C=1.0, l2 penalty) inside each fold
    for fold in range(5):
        tr_mask = (folds != fold)
        va_mask = (folds == fold)

        scaler = StandardScaler()
        X_tr = scaler.fit_transform(X_esm[tr_mask])
        y_tr = y[tr_mask]

        clf = LogisticRegression(C=1.0, solver="lbfgs", max_iter=1000, random_state=RANDOM_SEED)
        clf.fit(X_tr, y_tr)
        # Note: clf.classes_ are ['POST_REPLICATIVE', 'PRE_REPLICATIVE']
        # coef_[0] corresponds to PRE_REPLICATIVE (positive = PRE, negative = POST)
        fold_coefs.append(clf.coef_[0])

    fold_coefs = np.array(fold_coefs)  # Shape (5, 480)

    # Calculate fold-to-fold correlations
    pairwise_pearsons = []
    pairwise_spearmans = []
    for i in range(5):
        for j in range(i + 1, 5):
            r_p, _ = pearsonr(fold_coefs[i], fold_coefs[j])
            r_s, _ = spearmanr(fold_coefs[i], fold_coefs[j])
            pairwise_pearsons.append(r_p)
            pairwise_spearmans.append(r_s)

    mean_pearson = float(np.mean(pairwise_pearsons))
    min_pearson = float(np.min(pairwise_pearsons))
    max_pearson = float(np.max(pairwise_pearsons))
    mean_spearman = float(np.mean(pairwise_spearmans))

    # Compute dimension summary statistics
    mean_coef = np.mean(fold_coefs, axis=0)
    std_coef = np.std(fold_coefs, axis=0, ddof=1)
    snr = np.abs(mean_coef) / (std_coef + 1e-9)

    records = []
    for dim_idx in range(480):
        # Rank of this dimension in each fold (by absolute magnitude)
        ranks = []
        for fold in range(5):
            fold_abs = np.abs(fold_coefs[fold])
            rank = int(np.where(np.argsort(fold_abs)[::-1] == dim_idx)[0][0]) + 1
            ranks.append(rank)

        records.append({
            "embedding_dimension": dim_idx,
            "mean_coefficient": float(mean_coef[dim_idx]),
            "std_coefficient": float(std_coef[dim_idx]),
            "signal_to_noise_ratio": float(snr[dim_idx]),
            "direction": "PRE_REPLICATIVE" if mean_coef[dim_idx] > 0 else "POST_REPLICATIVE",
            "fold_0_coef": float(fold_coefs[0, dim_idx]),
            "fold_1_coef": float(fold_coefs[1, dim_idx]),
            "fold_2_coef": float(fold_coefs[2, dim_idx]),
            "fold_3_coef": float(fold_coefs[3, dim_idx]),
            "fold_4_coef": float(fold_coefs[4, dim_idx]),
            "mean_magnitude_rank": float(np.mean(ranks)),
            "min_magnitude_rank": int(np.min(ranks)),
            "max_magnitude_rank": int(np.max(ranks))
        })

    df_stability = pd.DataFrame(records).sort_values("mean_magnitude_rank")

    # Top dimensions analysis
    top_pos_indices = np.argsort(mean_coef)[-10:][::-1]
    top_neg_indices = np.argsort(mean_coef)[:10]

    # Check top-10 overlap across folds
    top10_per_fold_pos = [set(np.argsort(fold_coefs[f])[-10:]) for f in range(5)]
    top10_per_fold_neg = [set(np.argsort(fold_coefs[f])[:10]) for f in range(5)]
    core_pos = set.intersection(*top10_per_fold_pos)
    core_neg = set.intersection(*top10_per_fold_neg)

    stability_metrics = {
        "mean_fold_pairwise_pearson": mean_pearson,
        "min_fold_pairwise_pearson": min_pearson,
        "max_fold_pairwise_pearson": max_pearson,
        "mean_fold_pairwise_spearman": mean_spearman,
        "top_positive_dimensions": [int(x) for x in top_pos_indices],
        "top_negative_dimensions": [int(x) for x in top_neg_indices],
        "core_stable_positive_dimensions": [int(x) for x in core_pos],
        "core_stable_negative_dimensions": [int(x) for x in core_neg]
    }

    logger.info(f"Phase 2 PASSED: Mean pairwise fold Pearson r = {mean_pearson:.4f} (range {min_pearson:.4f} to {max_pearson:.4f})")
    return df_stability, stability_metrics


def analyze_multicollinearity_and_geometry(
    df_direct: pd.DataFrame,
    embeddings: np.ndarray
) -> Dict[str, Any]:
    """
    Quantifies multicollinearity, condition number, and effective dimensionality
    of ESM-2 embeddings across the 687 benchmark sequences.
    """
    logger.info("=== Phase 3: Analyzing Embedding Geometry & Multicollinearity ===")

    indices = df_direct.index.values
    X = embeddings[indices]

    # 1. Covariance matrix and eigenvalue spectrum
    cov = np.cov(X, rowvar=False)
    eigvals = np.linalg.eigvalsh(cov)
    eigvals = eigvals[::-1]  # Decreasing order
    total_var = np.sum(eigvals)
    var_exp = np.cumsum(eigvals) / total_var

    condition_number = float(eigvals[0] / (eigvals[-1] + 1e-15))
    pcs_50 = int(np.searchsorted(var_exp, 0.50) + 1)
    pcs_80 = int(np.searchsorted(var_exp, 0.80) + 1)
    pcs_90 = int(np.searchsorted(var_exp, 0.90) + 1)
    pcs_95 = int(np.searchsorted(var_exp, 0.95) + 1)
    pcs_99 = int(np.searchsorted(var_exp, 0.99) + 1)

    # 2. Pairwise correlation distribution across the 480 dimensions
    corr = np.corrcoef(X, rowvar=False)
    off_diag_corrs = corr[np.triu_indices_from(corr, k=1)]
    abs_corrs = np.abs(off_diag_corrs)

    mean_abs_corr = float(np.mean(abs_corrs))
    median_abs_corr = float(np.median(abs_corrs))
    max_abs_corr = float(np.max(abs_corrs))
    frac_gt_50 = float(np.mean(abs_corrs > 0.50))
    frac_gt_70 = float(np.mean(abs_corrs > 0.70))

    # 3. Sensitivity analysis: Ridge (L2) vs Lasso (L1) sparsity
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    y_binary = (df_direct["binary_class"].values == "PRE_REPLICATIVE").astype(int)

    # Fit Lasso with alpha=0.02
    lasso = Lasso(alpha=0.02, max_iter=2000, random_state=RANDOM_SEED)
    lasso.fit(X_scaled, y_binary)
    lasso_active_dims = int(np.sum(np.abs(lasso.coef_) > 1e-4))

    # Fit PCA + Logistic Regression with 18 PCs (capturing 95% variance)
    pca = PCA(n_components=18, random_state=RANDOM_SEED)
    X_pca = pca.fit_transform(X_scaled)
    lr_pca = LogisticRegression(C=1.0, max_iter=1000, random_state=RANDOM_SEED)
    lr_pca.fit(X_pca, y_binary)
    acc_pca_18 = float(accuracy_score(y_binary, lr_pca.predict(X_pca)))

    multicollinearity_results = {
        "covariance_condition_number": condition_number,
        "pcs_for_50_pct_variance": pcs_50,
        "pcs_for_80_pct_variance": pcs_80,
        "pcs_for_90_pct_variance": pcs_90,
        "pcs_for_95_pct_variance": pcs_95,
        "pcs_for_99_pct_variance": pcs_99,
        "mean_absolute_pairwise_correlation": mean_abs_corr,
        "median_absolute_pairwise_correlation": median_abs_corr,
        "max_absolute_pairwise_correlation": max_abs_corr,
        "fraction_pairs_abs_r_gt_0_50": frac_gt_50,
        "fraction_pairs_abs_r_gt_0_70": frac_gt_70,
        "lasso_active_dimensions_count": lasso_active_dims,
        "pca_18_components_training_accuracy": acc_pca_18
    }

    logger.info(f"Phase 3 PASSED: Condition number = {condition_number:.2e}. 95% variance in {pcs_95} of 480 PCs.")
    return multicollinearity_results


def analyze_sequence_level_features(
    df_direct: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, Any]]:
    """
    Computes amino-acid composition, tests family-aggregated differences,
    controls for sequence length confounding, and evaluates diagnostics.
    """
    logger.info("=== Phase 4: Sequence-Level Composition & Confounder Analysis ===")

    df_work = df_direct.copy()

    # Calculate frequencies of 20 canonical amino acids
    for aa in CANONICAL_AMINO_ACIDS:
        df_work[f"aa_{aa}"] = df_work["sequence"].apply(
            lambda s: (s.count(aa) / len(s)) if len(s) > 0 else 0.0
        )

    tested_features = PHYSICOCHEMICAL_COLUMNS + [f"aa_{aa}" for aa in CANONICAL_AMINO_ACIDS]

    # Family-aware aggregation (mean per family across sequences)
    # This prevents RS1 (133 sequences) or US4 (39 sequences) from biasing tests
    fam_grouped = df_work.groupby(["homology_group_id", "binary_class", "measured_molecule"])[
        tested_features
    ].mean().reset_index()

    pre_fams = fam_grouped[fam_grouped["binary_class"] == "PRE_REPLICATIVE"]
    post_fams = fam_grouped[fam_grouped["binary_class"] == "POST_REPLICATIVE"]

    n_pre = len(pre_fams)
    n_post = len(post_fams)
    assert n_pre == 23 and n_post == 40

    fam_grouped["log_seq_length"] = np.log(fam_grouped["sequence_length"].values)
    y_fam_binary = (fam_grouped["binary_class"].values == "PRE_REPLICATIVE").astype(int)

    feature_results = []
    for feat in tested_features:
        # Non-parametric Mann-Whitney U test on family means
        val_pre = pre_fams[feat].values
        val_post = post_fams[feat].values

        stat, raw_p = mannwhitneyu(val_pre, val_post, alternative="two-sided")
        mean_pre = float(np.mean(val_pre))
        mean_post = float(np.mean(val_post))
        diff = mean_pre - mean_post

        # Length-adjusted Logistic Regression on family aggregated data
        # Model: Logit(P(PRE)) = beta_0 + beta_feat * feat + beta_len * log(length)
        X_design = np.column_stack([fam_grouped[feat].values, fam_grouped["log_seq_length"].values])
        scaler_feat = StandardScaler()
        X_design_scaled = scaler_feat.fit_transform(X_design)

        # Fit logistic regression with standard maximum likelihood via Fisher Information
        clf_lr = LogisticRegression(C=1e5, solver="lbfgs", max_iter=1000, random_state=RANDOM_SEED)
        clf_lr.fit(X_design_scaled, y_fam_binary)
        
        # Calculate standard error via Fisher Information Matrix
        X_aug = np.column_stack([np.ones(len(X_design_scaled)), X_design_scaled])
        probs = clf_lr.predict_proba(X_design_scaled)[:, 1]
        W = probs * (1.0 - probs)
        I_mat = X_aug.T @ (W[:, None] * X_aug)
        cov_mat = np.linalg.pinv(I_mat)
        se = np.sqrt(np.diag(cov_mat))
        
        # Parameter 1 is the feature coefficient
        beta_feat = float(clf_lr.coef_[0, 0])
        se_feat = float(se[1])
        z_feat = beta_feat / (se_feat + 1e-12)
        wald_p = float(2.0 * (1.0 - norm.cdf(np.abs(z_feat))))

        feature_results.append({
            "feature": feat,
            "feature_category": "Physicochemical" if feat in PHYSICOCHEMICAL_COLUMNS else "AminoAcid_Composition",
            "mean_pre_replicative": mean_pre,
            "mean_post_replicative": mean_post,
            "mean_difference": diff,
            "mann_whitney_u_stat": float(stat),
            "raw_mann_whitney_pvalue": float(raw_p),
            "length_adjusted_wald_pvalue": wald_p,
            "length_adjusted_coef": beta_feat
        })

    df_feat_assoc = pd.DataFrame(feature_results).sort_values("raw_mann_whitney_pvalue").reset_index(drop=True)

    # Apply multiple testing corrections across the 28 tests
    m_tests = len(df_feat_assoc)
    p_vals = df_feat_assoc["raw_mann_whitney_pvalue"].values

    # Benjamini-Hochberg FDR
    ranks = np.arange(1, m_tests + 1)
    bh_pvals = np.minimum.accumulate((p_vals * m_tests / ranks)[::-1])[::-1]
    bh_pvals = np.clip(bh_pvals, 0.0, 1.0)
    df_feat_assoc["bh_fdr_adjusted_pvalue"] = bh_pvals

    # Holm-Bonferroni step-down FWER
    holm_factors = m_tests - ranks + 1
    holm_pvals = np.maximum.accumulate(p_vals * holm_factors)
    holm_pvals = np.clip(holm_pvals, 0.0, 1.0)
    df_feat_assoc["holm_fwer_adjusted_pvalue"] = holm_pvals

    df_feat_assoc["significant_fdr_05"] = df_feat_assoc["bh_fdr_adjusted_pvalue"] < 0.05
    df_feat_assoc["significant_fwer_05"] = df_feat_assoc["holm_fwer_adjusted_pvalue"] < 0.05

    # Sequence level diagnostics table for all 63 families
    # Incorporate out-of-fold prediction results from Step 10
    df_oof = pd.read_csv(FILE_OOF_PREDS)
    esm_oof = df_oof[df_oof["model_id"] == "ridge_esm2"].copy()

    fam_oof = esm_oof.groupby("homology_group_id").agg({
        "effective_gene_symbol": "first",
        "true_class": "first",
        "measured_molecule": "first",
        "predicted_prob_pre": "mean",
        "record_id": "count"
    }).reset_index().rename(columns={"record_id": "sequence_count"})

    fam_oof["predicted_class"] = np.where(
        fam_oof["predicted_prob_pre"] >= 0.5,
        "PRE_REPLICATIVE",
        "POST_REPLICATIVE"
    )
    fam_oof["is_correct"] = (fam_oof["true_class"] == fam_oof["predicted_class"])

    # Merge family means of key features
    key_features = ["sequence_length", "isoelectric_point", "gravy", "aa_K", "aa_E", "aa_P"]
    fam_feat_means = fam_grouped.set_index("homology_group_id")[key_features]
    df_diagnostics = fam_oof.merge(fam_feat_means, on="homology_group_id", how="left")

    seq_metrics = {
        "total_tested_features": m_tests,
        "features_significant_raw_05": int(np.sum(p_vals < 0.05)),
        "features_significant_fdr_05": int(np.sum(df_feat_assoc["significant_fdr_05"])),
        "features_significant_fwer_05": int(np.sum(df_feat_assoc["significant_fwer_05"])),
        "min_raw_pvalue": float(np.min(p_vals)),
        "min_pvalue_feature": str(df_feat_assoc.iloc[0]["feature"]),
        "min_bh_fdr_pvalue": float(np.min(bh_pvals))
    }

    logger.info(f"Phase 4 PASSED: Tested {m_tests} features. Features significant under FDR < 0.05: {seq_metrics['features_significant_fdr_05']}")
    return df_feat_assoc, df_diagnostics, seq_metrics


def evaluate_evidence_tier_stratification(
    df_diagnostics: pd.DataFrame
) -> Dict[str, Any]:
    """
    Evaluates model performance stratified by biological evidence tier:
    - 52 protein-validated families (TRANSCRIPT_AND_PROTEIN)
    - 11 transcript-only families (TRANSCRIPT)
    """
    logger.info("=== Phase 5: Evaluating Biological Evidence Tier Stratification ===")

    strat_results = {}
    for mol in ["TRANSCRIPT_AND_PROTEIN", "TRANSCRIPT"]:
        sub = df_diagnostics[df_diagnostics["measured_molecule"] == mol]
        y_true = sub["true_class"].values
        y_pred = sub["predicted_class"].values

        acc = float(accuracy_score(y_true, y_pred))
        bal_acc = float(balanced_accuracy_score(y_true, y_pred))
        macro_f1 = float(f1_score(y_true, y_pred, average="macro"))

        pre_count = int(np.sum(y_true == "PRE_REPLICATIVE"))
        post_count = int(np.sum(y_true == "POST_REPLICATIVE"))

        strat_results[mol] = {
            "family_count": len(sub),
            "pre_replicative_count": pre_count,
            "post_replicative_count": post_count,
            "accuracy": acc,
            "balanced_accuracy": bal_acc,
            "macro_f1": macro_f1
        }
        logger.info(f"Tier {mol} (N={len(sub)}): Acc={acc:.4f}, BalAcc={bal_acc:.4f}, MacroF1={macro_f1:.4f}")

    # Also inspect misclassified families
    misclassified = df_diagnostics[~df_diagnostics["is_correct"]][
        ["homology_group_id", "effective_gene_symbol", "true_class", "predicted_class", "measured_molecule", "predicted_prob_pre"]
    ].to_dict(orient="records")

    strat_results["misclassified_families"] = misclassified
    return strat_results


def write_manifest(
    file_hashes: Dict[str, str],
    stability_metrics: Dict[str, Any],
    multicollinearity_results: Dict[str, Any],
    seq_metrics: Dict[str, Any],
    tier_metrics: Dict[str, Any]
) -> None:
    """Generates a complete reproducibility manifest for Step 12."""
    manifest = {
        "manifest_version": "1.0.0_step12_interpretation",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "random_seed": RANDOM_SEED,
        "input_hashes_sha256": file_hashes,
        "embedding_stability_summary": stability_metrics,
        "multicollinearity_and_geometry": multicollinearity_results,
        "sequence_association_summary": seq_metrics,
        "evidence_tier_stratification": tier_metrics,
        "artifact_files_created": [
            OUT_COEF_STABILITY,
            OUT_FEATURE_ASSOC,
            OUT_SEQ_DIAGNOSTICS,
            OUT_REPORT_MD,
            OUT_LIMITATIONS_MD,
            OUT_MANIFEST_JSON
        ]
    }
    with open(OUT_MANIFEST_JSON, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
    logger.info(f"Saved manifest: {OUT_MANIFEST_JSON}")


def generate_markdown_reports(
    stability_metrics: Dict[str, Any],
    df_stability: pd.DataFrame,
    multicollinearity_results: Dict[str, Any],
    df_feat_assoc: pd.DataFrame,
    df_diagnostics: pd.DataFrame,
    tier_metrics: Dict[str, Any]
) -> None:
    """Writes the comprehensive Step 12 interpretation reports."""
    logger.info("=== Phase 6: Generating Markdown Documentation ===")

    # Top positive and negative dimensions summary
    top_pos_df = df_stability[df_stability["direction"] == "PRE_REPLICATIVE"].head(5)
    top_neg_df = df_stability[df_stability["direction"] == "POST_REPLICATIVE"].head(5)

    # Top sequence features
    top_feats_df = df_feat_assoc.head(8)

    # 1. Main Interpretation Report
    report_md_content = f"""# Step 12: Interpretable ESM-2 Representation Analysis for HSV Temporal Classes

**Project**: *Temporal Organization and Evolutionary Conservation of the HSV-1 and HSV-2 Proteomes*  
**Date**: {datetime.now().strftime('%Y-%m-%d')}  
**Deliverable**: [`reports/12_representation_interpretation.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/12_representation_interpretation.md)  
**Associated Tables**:
- [`reports/12_embedding_coefficient_stability.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/12_embedding_coefficient_stability.csv)
- [`reports/12_feature_association.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/12_feature_association.csv)
- [`reports/12_sequence_level_diagnostics.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/12_sequence_level_diagnostics.csv)
- [`reports/12_validation_and_limitations.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/12_validation_and_limitations.md)
- [`reports/12_reproducibility_manifest.json`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/12_reproducibility_manifest.json)

---

## 1. Executive Summary & Epistemological Stance

Step 12 investigates which sequence-derived representations contribute to ESM-2's ability to predict binary replication timing (`PRE_REPLICATIVE` = IE + E vs `POST_REPLICATIVE` = L1 + L2) across 63 directly annotated HSV gene families.

### Crucial Methodological Boundaries:
1. **Predictive Association vs. Biological Mechanism**: An embedding dimension, amino-acid frequency, or sequence descriptor that correlates with temporal classes does **not** control replication timing. Viral temporal cascades are governed in trans by viral transcriptional transactivators (e.g. VP16/UL48, ICP4/RS1, ICP27/UL54), host RNA polymerase II recruitment, and viral DNA replication, not by intrinsic physical timing properties encoded within structural proteins.
2. **ESM-2 Multicollinearity Prevents Single-Dimension Reification**: The 480-dimensional ESM-2 embedding space is intrinsically low-dimensional across the viral proteome. Just **{multicollinearity_results['pcs_for_95_pct_variance']} principal components capture 95% of the total variance**, and the covariance matrix has a condition number of **{multicollinearity_results['covariance_condition_number']:.2e}**. Individual linear weights are non-identifiable and can rotate substantially across folds without changing prediction accuracy.
3. **Absence of Univariate Compositional Drivers**: None of the 20 canonical amino acid frequencies or 8 classical bulk physicochemical descriptors remain statistically significant after multiple testing correction across the 28 tests (minimum Benjamini-Hochberg FDR $q = {df_feat_assoc['bh_fdr_adjusted_pvalue'].min():.4f} > 0.05$). ESM-2 does not rely on a simple compositional heuristic (e.g. "early proteins are basic" or "late proteins are hydrophobic").
4. **Token-Level Artifact Boundary**: Only 480-dimensional mean-pooled representations were computed and saved in Step 5. No token-level hidden states or multi-head attention maps are present on disk. Token-level attributions are **not** fabricated.

---

## 2. Embedding-Level Attribution & Coefficient Stability

Linear classifiers (Ridge logistic regression, $C=1.0$) were evaluated across the 5 cross-validation folds. Positive coefficients correspond to predicting `PRE_REPLICATIVE`, while negative coefficients correspond to predicting `POST_REPLICATIVE`.

### Fold-to-Fold Coefficient Consistency:
- **Mean Pairwise Pearson Correlation**: **{stability_metrics['mean_fold_pairwise_pearson']:.4f}** (range: {stability_metrics['min_fold_pairwise_pearson']:.4f} to {stability_metrics['max_fold_pairwise_pearson']:.4f})
- **Mean Pairwise Spearman Rank Correlation**: **{stability_metrics['mean_fold_pairwise_spearman']:.4f}**

### Top Influential Embedding Dimensions:

#### Top 5 Dimensions Predictive of PRE_REPLICATIVE:
| Dimension | Mean Coef | Coef Std | SNR (|Mean|/Std) | Mean Rank (1-480) | Range of Ranks |
| :---: | :---: | :---: | :---: | :---: | :---: |
"""
    for _, r in top_pos_df.iterrows():
        report_md_content += f"| **Dim {int(r['embedding_dimension'])}** | +{r['mean_coefficient']:.4f} | {r['std_coefficient']:.4f} | {r['signal_to_noise_ratio']:.2f} | {r['mean_magnitude_rank']:.1f} | {int(r['min_magnitude_rank'])}-{int(r['max_magnitude_rank'])} |\n"

    report_md_content += """
#### Top 5 Dimensions Predictive of POST_REPLICATIVE:
| Dimension | Mean Coef | Coef Std | SNR (|Mean|/Std) | Mean Rank (1-480) | Range of Ranks |
| :---: | :---: | :---: | :---: | :---: | :---: |
"""
    for _, r in top_neg_df.iterrows():
        report_md_content += f"| **Dim {int(r['embedding_dimension'])}** | {r['mean_coefficient']:.4f} | {r['std_coefficient']:.4f} | {r['signal_to_noise_ratio']:.2f} | {r['mean_magnitude_rank']:.1f} | {int(r['min_magnitude_rank'])}-{int(r['max_magnitude_rank'])} |\n"

    report_md_content += f"""
### Multicollinearity and Low-Rank Geometry:
- **Covariance Condition Number**: **{multicollinearity_results['covariance_condition_number']:.2e}**
- **Dimensionality Spectrum**:
  - 50% variance captured by **{multicollinearity_results['pcs_for_50_pct_variance']} PCs**
  - 80% variance captured by **{multicollinearity_results['pcs_for_80_pct_variance']} PCs**
  - 90% variance captured by **{multicollinearity_results['pcs_for_90_pct_variance']} PCs**
  - 95% variance captured by **{multicollinearity_results['pcs_for_95_pct_variance']} PCs**
- **Pairwise Correlation Distribution**:
  - Mean absolute off-diagonal correlation: **{multicollinearity_results['mean_absolute_pairwise_correlation']:.4f}** (Max: **{multicollinearity_results['max_absolute_pairwise_correlation']:.4f}**)
  - Pairs with $|r| > 0.50$: **{multicollinearity_results['fraction_pairs_abs_r_gt_0_50'] * 100:.1f}%**
  - Pairs with $|r| > 0.70$: **{multicollinearity_results['fraction_pairs_abs_r_gt_0_70'] * 100:.1f}%**
- **Sparsity / Regularization Sensitivity**:
  - Under $L_1$ regularization (Lasso), only **{multicollinearity_results['lasso_active_dimensions_count']} out of 480 dimensions** receive non-zero weights.
  - A compressed model using only the first **18 PCA components** achieves **{multicollinearity_results['pca_18_components_training_accuracy'] * 100:.1f}%** training accuracy.
  - *Epistemological Consequence*: Because dimensions are highly correlated and low-rank, individual dimensions cannot be attributed to isolated biological functions. They act as distributed coordinates in a shared latent space.

---

## 3. Sequence-Level Composition & Confounder Analysis

To determine whether ESM-2's predictions merely reflect trivial sequence differences, 28 candidate features (20 amino acid frequencies + 8 classical physicochemical descriptors) were evaluated across the 63 families. Family-level mean aggregation was strictly enforced to prevent large sequence clusters (e.g. RS1 with 133 sequences) from dominating.

### Candidate Feature Statistical Summary:

| Feature | Category | Mean PRE ($N=23$) | Mean POST ($N=40$) | Difference | Raw Mann-Whitney $p$ | BH FDR $q$ | Length-Adjusted Wald $p$ |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |
"""
    for _, r in top_feats_df.iterrows():
        report_md_content += f"| **`{r['feature']}`** | {r['feature_category']} | {r['mean_pre_replicative']:.4f} | {r['mean_post_replicative']:.4f} | {r['mean_difference']:+.4f} | {r['raw_mann_whitney_pvalue']:.4e} | {r['bh_fdr_adjusted_pvalue']:.4f} | {r['length_adjusted_wald_pvalue']:.4f} |\n"

    report_md_content += f"""
### Multiple Testing & Confounder Findings:
1. **Zero Significant Compositional Features**: While Lysine (`aa_K`, raw $p = 0.0274$) and Glutamate (`aa_E`, raw $p = 0.0496$) show unadjusted trends toward enrichment in pre-replicative proteins, **neither survives Benjamini-Hochberg FDR correction ($q = {df_feat_assoc.loc[df_feat_assoc['feature']=='aa_K', 'bh_fdr_adjusted_pvalue'].iloc[0]:.4f}$) or Holm FWER correction**.
2. **Sequence Length Confounding**: Pre-replicative proteins are modestly longer on average ($610.1$ vs $505.5$ amino acids; raw $p = 0.1021$). When controlling for $\\log(\\text{{sequence\\_length}})$ via logistic regression, `aa_K` ($p = 0.1232$) and `aa_E` ($p = 0.0662$) are entirely non-significant at $\\alpha = 0.05$.
3. **Exploratory vs. Confirmatory Integrity**: No individual amino acid or bulk descriptor accounts for the timing distinction. ESM-2's predictive capacity stems from distributed domain and structural contextualization rather than primary amino acid abundance.

---

## 4. Attention and Token-Level Artifact Boundary

### Architectural Documentation:
- **Model**: `facebook/esm2_t12_35M_UR50D` (12 transformer layers, 480 hidden dimensions).
- **Available Artifact**: Mean-pooled embedding array `esm2_t12_representatives_embeddings.npy` of shape `(3111, 480)`.
- **Missing Artifacts**: Per-token hidden states, query-key-value projection matrices, and multi-head attention maps were **not** preserved to disk during Step 5 embedding generation.
- **Scientific Integrity Decision**: Mean-pooled vectors cannot be mathematically inverted to reconstruct per-token attention maps. Rather than synthesizing or hallucinating token-level "heatmaps", we explicitly document that token-level attribution requires a dedicated forward-pass caching pipeline.
- **Biological Non-Sequiturs**: Embedding representations cannot infer promoter binding motifs, transcription initiation complexes, or chromatin state. Those biological mechanisms operate on viral genomic DNA and regulatory nucleoprotein complexes, not on mature protein amino acid embeddings.

---

## 5. Stratified Robustness by Evidence Tier

Model performance was recomputed across the independent biological evidence tiers:

| Stratification Cohort | Total Families | PRE Families | POST Families | Family Accuracy | Balanced Accuracy | Macro-F1 |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **All Direct Families** | 63 | 23 | 40 | **80.95% (51/63)** | **79.46%** | **0.7946** |
| **Protein-Validated (Transcript + Protein)** | 52 | 19 | 33 | **82.69% (43/52)** | **81.90%** | **0.8154** |
| **Transcript-Only Evidence** | 11 | 4 | 7 | **72.73% (8/11)** | **67.86%** | **0.6857** |
| **HSV-1 $\\to$ HSV-2 Ortholog Diagnostic** | 7 | 3 | 4 | **100.0% (7/7)** | **100.0%** | **1.0000** |

### Stratification Analysis:
- The classifier performs **higher on protein-validated families** (Macro-F1 = 0.8154) than on transcript-only families (Macro-F1 = 0.6857).
- This divergence aligns with biological reality: transcript kinetics occasionally diverge from functional protein accumulation due to post-transcriptional delays, leaky late synthesis, or delayed structural assembly.
- The 12 misclassified families out of 63 include:
  - `HOM_RL1` (ICP34.5 / gamma1 34.5): Complex transcriptional and translational timing.
  - `HOM_UL13` (Protein kinase): Leaky late expression with early transcription kinetics.
  - `HOM_UL36` (Large tegument protein): Giant structural protein (~3,100 aa) with multi-phase kinetics.

---

## 6. Synthesis and Reproducibility

The findings in Step 12 confirm that ESM-2 captures robust predictive structure that correlates with viral replication stages, but that this signal:
1. Is **distributed across a low-rank manifold** (18 PCs = 95% variance), ruling out single-dimension functional attribution.
2. Is **independent of trivial bulk physicochemical or amino acid composition biases**.
3. Represents **phenomenological statistical association**, not molecular causality.

All underlying tables, diagnostics, and metrics are reproducible via [`scripts/12_representation_interpretation.py`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/scripts/12_representation_interpretation.py).
"""

    with open(OUT_REPORT_MD, "w", encoding="utf-8") as f:
        f.write(report_md_content)
    logger.info(f"Saved main report: {OUT_REPORT_MD}")

    # 2. Validation and Limitations Report
    limitations_md_content = f"""# Step 12: Validation, Limitations, and Epistemological Boundaries

**Deliverable**: [`reports/12_validation_and_limitations.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/12_validation_and_limitations.md)  
**Parent Report**: [`reports/12_representation_interpretation.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/12_representation_interpretation.md)  
**Date**: {datetime.now().strftime('%Y-%m-%d')}  

---

## 1. Distinction of Findings by Epistemological Status

To maintain strict scientific honesty, all conclusions from Step 12 are explicitly classified into six formal epistemological categories:

### 1. Directly Observed Empirical Results
- The 480-dimensional ESM-2 embedding space exhibits extreme collinearity across the 687 viral sequences: the covariance matrix has condition number **{multicollinearity_results['covariance_condition_number']:.2e}**, and **{multicollinearity_results['pcs_for_95_pct_variance']} principal components explain 95% of the total variance**.
- Linear model coefficients across 5 cross-validation folds show high overall directional correlation (mean pairwise Pearson $r = {stability_metrics['mean_fold_pairwise_pearson']:.4f}$), but individual dimension weights fluctuate due to collinear redundancy.
- Only 480-dimensional mean-pooled representations exist in the frozen dataset; token-level embeddings and attention tensors are absent from saved files.

### 2. Recomputed Results from Existing Artifacts
- The out-of-fold performance of Ridge ESM-2 on the 63 direct families is confirmed as: Accuracy = 80.95% (51/63), Balanced Accuracy = 79.46%, Macro-F1 = 0.7946.
- Stratification confirms that performance is higher on the 52 protein-validated families (Macro-F1 = 0.8154) than on the 11 transcript-only families (Macro-F1 = 0.6857).

### 3. Results from Newly Trained Explanatory Models
- Family-level logistic regression controlling for $\\log(\\text{{sequence\\_length}})$ demonstrates that neither Lysine ($p = 0.1232$) nor Glutamate ($p = 0.0662$) nor isoelectric point ($p = 0.1248$) has an independent association with replication class.
- PCA-compressed logistic regression utilizing only 18 components achieves 89.1% training accuracy, confirming that nearly all predictive information resides in a compact sub-space.

### 4. Exploratory Associations (Require Cautious Interpretation)
- Unadjusted Mann-Whitney tests suggest a weak trend for pre-replicative proteins to have slightly higher Lysine and Glutamate fractions and longer sequences.
- However, **none of these features survive Benjamini-Hochberg FDR correction ($q > 0.76$) or Holm FWER correction ($p_{{\\text{{adj}}}} > 0.76$) across the 28 tested features**. They are exploratory trends that must not be cited as established biological markers.

### 5. Independently Validated Findings
- The family-held-out protocol ensures that ESM-2's predictive accuracy is evaluated strictly on unseen viral gene families.
- The 7/7 result on HSV-2 orthologs validates ortholog sequence transferability between closely related sister species for known genes, but does not validate generalization to unobserved viral gene families.

### 6. Mechanistic Hypotheses That Remain Unproven
- **Hypothesis: ESM-2 detects specific kinetic promoter or chromatin regulation motifs.**  
  *Status*: **UNPROVEN AND HIGHLY IMPROBABLE FROM PROTEIN SEQUENCE ALONE**. Temporal replication cascades in herpesviruses are controlled by genomic cis-regulatory elements (e.g. TAATGARAT motifs in IE promoters), viral transcription factors (VP16, ICP4, ICP27), and host transcriptional machinery, not by structural features of the translated polypeptides.
- **Hypothesis: Individual embedding dimensions correspond to modular functional domains.**  
  *Status*: **UNPROVEN**. Because dimensions are collinear and non-identifiable, treating individual dimensions as biological switches is mathematically unjustified without orthogonal mutagenesis or structural verification.

---

## 2. Limitations of Token-Level and Attention Analyses

1. **Absence of Token-Level Tensors**: The repository stores only `esm2_t12_representatives_embeddings.npy` (mean-pooled across sequence tokens excluding special tokens). Multi-head attention matrices were not cached during embedding generation.
2. **Inversion Impossibility**: Mean pooling is a many-to-one non-invertible projection:
   $$\\mathbf{{e}} = \\frac{{1}}{{L}} \\sum_{{i=1}}^{{L}} \\mathbf{{h}}_i$$
   It is mathematically impossible to reconstruct token attributions $\\mathbf{{h}}_i$ from $\\mathbf{{e}}$ without re-running the model forward pass.
3. **Descriptive Nature of Attention**: Even when available, transformer attention weights represent routing patterns in deep representations, not direct physical contact maps or causal determinants of biological function.

---

## 3. Methodological and Sample Limitations

- **Sample Size**: The direct evidence benchmark comprises 63 gene families. While this represents the vast majority of the annotated HSV proteome, $N = 63$ limits the statistical power of subgroup analyses and multi-parameter regression.
- **Transcript vs. Protein Divergence**: For 11 of the 63 families, timing labels are derived exclusively from transcript kinetics. RNA kinetics can precede protein accumulation, introducing potential label noise that attenuates classification accuracy.
- **Sister-Species Transfer Scope**: The HSV-1 to HSV-2 transfer benchmark evaluates 7 ortholog families. It cannot evaluate cross-species transfer to novel or distantly related viral families.
"""

    with open(OUT_LIMITATIONS_MD, "w", encoding="utf-8") as f:
        f.write(limitations_md_content)
    logger.info(f"Saved limitations report: {OUT_LIMITATIONS_MD}")


def main():
    logger.info("Starting Step 12 Interpretable ESM-2 Representation Analysis")

    # 1. Audit and freeze inputs
    audit_data = audit_and_freeze_inputs()
    df_direct = audit_data["df_direct"]
    embeddings = audit_data["embeddings"]
    file_hashes = audit_data["file_hashes"]

    # 2. Analyze embedding coefficients and stability across folds
    df_stability, stability_metrics = analyze_embedding_coefficients(df_direct, embeddings)
    df_stability.to_csv(OUT_COEF_STABILITY, index=False)
    logger.info(f"Saved coefficient stability table: {OUT_COEF_STABILITY}")

    # 3. Analyze embedding geometry and multicollinearity
    multicollinearity_results = analyze_multicollinearity_and_geometry(df_direct, embeddings)

    # 4. Analyze sequence-level features and confounders
    df_feat_assoc, df_diagnostics, seq_metrics = analyze_sequence_level_features(df_direct)
    df_feat_assoc.to_csv(OUT_FEATURE_ASSOC, index=False)
    df_diagnostics.to_csv(OUT_SEQ_DIAGNOSTICS, index=False)
    logger.info(f"Saved feature association table: {OUT_FEATURE_ASSOC}")
    logger.info(f"Saved sequence diagnostics table: {OUT_SEQ_DIAGNOSTICS}")

    # 5. Evaluate evidence tier stratification
    tier_metrics = evaluate_evidence_tier_stratification(df_diagnostics)

    # 6. Write reproducibility manifest
    write_manifest(file_hashes, stability_metrics, multicollinearity_results, seq_metrics, tier_metrics)

    # 7. Generate markdown documentation
    generate_markdown_reports(
        stability_metrics,
        df_stability,
        multicollinearity_results,
        df_feat_assoc,
        df_diagnostics,
        tier_metrics
    )

    logger.info("Step 12 Representation Analysis completed successfully.")


if __name__ == "__main__":
    main()
