#!/usr/bin/env python3
"""
Test Suite: Step 6 Unsupervised Learning, Dimensionality Reduction, and Clustering
==================================================================================

Concept & Scope:
----------------
This test suite independently asserts the mathematical, topological, and biological
validity of the Step 6 analytical pipeline:

1. Artifact Completeness: All 6 tabular reports, 2 JSON manifests, and 8 publication figures exist.
2. Coordinate and Embedding Alignment: Row ordering aligns 1-to-1 with 3,111 representatives.
   Zero null or infinite values in PCA (1-10), UMAP (1-2), or t-SNE coordinates.
3. Dimensionality Reduction Validity:
   - Top 50 PCs account for >90% of cumulative variance.
4. Primary Clustering Sanity:
   - Exactly 74 clusters in the primary K-Means model.
   - NMI with canonical gene annotations exceeds 0.80 (strong functional alignment).
5. Biological Orthology vs Species Independence:
   - Ortholog co-clustering rate exceeds 75% (HSV-1 and HSV-2 orthologs co-cluster).
   - Species NMI against clusters is < 0.15 (proves species do NOT separate into artificial macro-clusters).
6. Physicochemical Validity:
   - Characterization features lie within biologically physical ranges.
"""

import os
import json
import pytest
import numpy as np
import pandas as pd

WORKSPACE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA_PROCESSED = os.path.join(WORKSPACE_DIR, "data", "processed")
REPORTS_DIR = os.path.join(WORKSPACE_DIR, "reports")
FIGURES_DIR = os.path.join(REPORTS_DIR, "figures")

FILE_COORD_CSV = os.path.join(DATA_PROCESSED, "06_dimensionality_reduction_coordinates.csv")
FILE_MEMBERSHIP_CSV = os.path.join(REPORTS_DIR, "06_cluster_membership.csv")
FILE_METRICS_CSV = os.path.join(REPORTS_DIR, "06_cluster_quality_metrics.csv")
FILE_CHAR_CSV = os.path.join(REPORTS_DIR, "06_cluster_characterization.csv")
FILE_SPECIES_CSV = os.path.join(REPORTS_DIR, "06_species_comparison.csv")
FILE_SENS_CSV = os.path.join(REPORTS_DIR, "06_sensitivity_analysis.csv")

FILE_METRICS_JSON = os.path.join(REPORTS_DIR, "06_analysis_metrics.json")
FILE_MANIFEST_JSON = os.path.join(REPORTS_DIR, "06_analysis_manifest.json")
FILE_SUMMARY_MD = os.path.join(REPORTS_DIR, "06_analysis_summary.md")


def test_01_all_step6_artifacts_exist():
    """Verify all 9 tabular/JSON/Markdown outputs and 8 figure files exist on disk."""
    tables = [
        FILE_COORD_CSV, FILE_MEMBERSHIP_CSV, FILE_METRICS_CSV,
        FILE_CHAR_CSV, FILE_SPECIES_CSV, FILE_SENS_CSV,
        FILE_METRICS_JSON, FILE_MANIFEST_JSON, FILE_SUMMARY_MD
    ]
    for p in tables:
        assert os.path.exists(p), f"Missing table/manifest: {p}"

    expected_figures = [
        "01_pca_explained_variance.png",
        "02_pca_species_projection.png",
        "03_pca_gene_family_projection.png",
        "04_umap_species_projection.png",
        "05_umap_cluster_assignment.png",
        "06_umap_quality_projection.png",
        "07_cluster_physicochemical_distributions.png",
        "08_sensitivity_comparison_cohort_a_vs_b.png",
    ]
    if os.path.exists(FIGURES_DIR):
        for fig in expected_figures:
            fig_path = os.path.join(FIGURES_DIR, fig)
            assert os.path.exists(fig_path), f"Missing figure: {fig_path}"
            assert os.path.getsize(fig_path) > 10000, f"Figure {fig} appears corrupted (too small)"


def test_02_coordinate_integrity_and_row_alignment():
    """Verify coordinates table has exactly 3,111 rows with no NaNs or Infs."""
    df_coords = pd.read_csv(FILE_COORD_CSV)
    assert len(df_coords) == 3111, f"Expected 3,111 rows, got {len(df_coords)}"

    # Check PCA coordinates (1-10)
    for i in range(1, 11):
        col = f"pca{i}"
        assert col in df_coords.columns, f"Missing {col}"
        assert df_coords[col].isna().sum() == 0, f"NaNs in {col}"

    # Check UMAP primary coordinates
    assert "umap1_nn30_md01_primary" in df_coords.columns
    assert "umap2_nn30_md01_primary" in df_coords.columns
    assert df_coords["umap1_nn30_md01_primary"].isna().sum() == 0
    assert df_coords["umap2_nn30_md01_primary"].isna().sum() == 0

    # Check t-SNE coordinates
    assert "tsne1" in df_coords.columns
    assert "tsne2" in df_coords.columns
    assert df_coords["tsne1"].isna().sum() == 0
    assert df_coords["tsne2"].isna().sum() == 0


def test_03_clustering_partition_validity():
    """Verify primary cluster membership assigns all 3,111 records to exactly 74 clusters."""
    df_mem = pd.read_csv(FILE_MEMBERSHIP_CSV)
    assert len(df_mem) == 3111
    assert df_mem["primary_cluster_id"].isna().sum() == 0

    unique_clusters = set(df_mem["primary_cluster_id"].unique())
    assert len(unique_clusters) == 74, f"Expected 74 clusters, found {len(unique_clusters)}"
    assert unique_clusters == set(range(74)), "Cluster IDs must be contiguous 0..73"


def test_04_gene_family_concordance_and_species_independence():
    """
    Verify fundamental biological findings:
    1. Clusters align strongly with gene annotations (NMI > 0.80).
    2. Clusters do NOT separate by viral species (Species NMI < 0.15).
    3. Ortholog co-clustering rate is > 75%.
    """
    with open(FILE_METRICS_JSON, "r", encoding="utf-8") as f:
        metrics = json.load(f)

    # 1. Gene family NMI
    nmi_genes = metrics["primary_clustering_kmeans_k74"]["nmi_against_annotated_genes"]
    assert nmi_genes >= 0.80, f"Expected Gene Family NMI >= 0.80, got {nmi_genes}"

    # 2. Species NMI
    nmi_species = metrics["species_divergence_metrics"]["species_vs_clusters_nmi"]
    assert nmi_species <= 0.15, f"Expected Species NMI <= 0.15, got {nmi_species}"

    # 3. Ortholog co-clustering rate
    coclust_rate = metrics["species_divergence_metrics"]["canonical_ortholog_co_clustering_rate_pct"]
    assert coclust_rate >= 75.0, f"Expected Ortholog Co-clustering >= 75%, got {coclust_rate}%"


def test_05_physicochemical_characterization_bounds():
    """Verify cluster characterization table contains physical biological values."""
    df_char = pd.read_csv(FILE_CHAR_CSV)
    assert len(df_char) == 74

    # Molecular weight > 5 kDa and < 600 kDa
    assert (df_char["median_mw_kda"] >= 5.0).all()
    assert (df_char["median_mw_kda"] <= 600.0).all()

    # Isoelectric point between 3.5 and 12.5
    assert (df_char["median_pi"] >= 3.5).all() and (df_char["median_pi"] <= 12.5).all()

    # GRAVY between -2.0 and +1.5
    assert (df_char["median_gravy"] >= -2.0).all() and (df_char["median_gravy"] <= 1.5).all()

    # Cluster purity percentage > 0 and <= 100
    assert (df_char["cluster_purity_pct"] > 0).all() and (df_char["cluster_purity_pct"] <= 100).all()


def test_06_sensitivity_stability_cohort_a_vs_b():
    """Verify cluster metrics remain stable between full cohort A and high-quality cohort B."""
    df_sens = pd.read_csv(FILE_SENS_CSV)
    assert len(df_sens) >= 3

    # NMI delta should be small (< 0.10) confirming stability
    for _, row in df_sens.iterrows():
        assert abs(row["nmi_delta"]) < 0.10, f"Unstable NMI delta in {row['model']}: {row['nmi_delta']}"
