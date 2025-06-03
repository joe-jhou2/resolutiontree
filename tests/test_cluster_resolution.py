# tests/test_cluster_resolution.py
from __future__ import annotations

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import re
import pandas as pd
import pytest
import scanpy as sc
from src import cluster_resolution_finder
from testing.scanpy._helpers.data import pbmc68k_reduced


@pytest.fixture
def adata_for_test():
    """Fixture to provide a preprocessed AnnData object for testing."""
    import scanpy as sc

    adata = pbmc68k_reduced()
    sc.pp.neighbors(adata)
    return adata


# Test 1: Basic functionality
def test_cluster_resolution_finder_basic(adata_for_test):
    """Test that cluster_resolution_finder runs without errors and modifies adata."""
    adata = adata_for_test.copy()  # Create a copy to avoid modifying the fixture
    resolutions = [0.1, 0.5]
    result = cluster_resolution_finder(
        adata,
        resolutions,
        prefix="leiden_res_",
        method="wilcoxon",
        n_top_genes=2,
        min_cells=2,
        deg_mode="within_parent",
        flavor="igraph",
        n_iterations=2,
    )

    # Check that clustering columns were added to adata.obs
    for res in resolutions:
        assert f"leiden_res_{res}" in result.obs

    # Check that top_genes_dict was added to adata.uns
    assert "cluster_resolution_top_genes" in result.uns
    top_genes_dict = result.uns["cluster_resolution_top_genes"]
    assert isinstance(top_genes_dict, dict)
    assert len(top_genes_dict) > 0

    for k, genes in top_genes_dict.items():
        parent, child = k.split("_", 1)  # Split the combined key back into tuple
        assert isinstance(parent, str)
        assert isinstance(child, str)
        assert isinstance(genes, list)
        assert all(isinstance(g, str) for g in genes)

    # Check that cluster_data was added to adata.uns
    assert "cluster_resolution_cluster_data" in result.uns
    cluster_data = result.uns["cluster_resolution_cluster_data"]
    assert isinstance(cluster_data, pd.DataFrame)
    for res in resolutions:
        assert f"leiden_res_{res}" in cluster_data.columns


# Test 2: Conflicting arguments (invalid deg_mode)
def test_cluster_resolution_finder_invalid_deg_mode(adata_for_test):
    """Test that an invalid deg_mode raises a ValueError."""
    adata = adata_for_test.copy()
    with pytest.raises(
        ValueError, match=r"deg_mode must be 'within_parent' or 'per_resolution'"
    ):
        cluster_resolution_finder(
            adata,
            resolutions=[0.1],
            deg_mode="invalid_mode",  # type: ignore[arg-type]
        )


# Test 3: Input values that should cause an error (empty resolutions)
def test_cluster_resolution_finder_empty_resolutions(adata_for_test):
    """Test that an empty resolutions list raises a ValueError."""
    adata = adata_for_test.copy()
    with pytest.raises(ValueError, match=r"resolutions list cannot be empty"):
        cluster_resolution_finder(
            adata,
            resolutions=[],
        )


# Test 4: Input values that should cause an error (negative resolutions)
def test_cluster_resolution_finder_negative_resolutions(adata_for_test):
    """Test that negative resolutions raise a ValueError."""
    adata = adata_for_test.copy()
    with pytest.raises(
        ValueError, match="All resolutions must be non-negative numbers"
    ):
        cluster_resolution_finder(
            adata,
            resolutions=[0.1, -0.5],
        )


# Test 5: Input values that should cause an error (missing neighbors)
def test_cluster_resolution_finder_missing_neighbors():
    """Test that an adata object without neighbors raises a ValueError."""
    adata = sc.datasets.pbmc68k_reduced()  # Create a fresh adata
    # Remove neighbors if they exist
    if "neighbors" in adata.uns:
        del adata.uns["neighbors"]
    # Also remove connectivities and distances to ensure leiden doesn't recompute
    if "connectivities" in adata.obsp:
        del adata.obsp["connectivities"]
    if "distances" in adata.obsp:
        del adata.obsp["distances"]
    with pytest.raises(
        ValueError,
        match=re.escape(
            "adata must have precomputed neighbors (run sc.pp.neighbors first)."
        ),
    ):
        cluster_resolution_finder(
            adata,
            resolutions=[0.1],
        )


# Test 6: Helpful error message (unsupported method)
def test_cluster_resolution_finder_unsupported_method(adata_for_test):
    """Test that an unsupported method raises a ValueError with a helpful message."""
    adata = adata_for_test.copy()
    with pytest.raises(ValueError, match="Only method='wilcoxon' is supported"):
        cluster_resolution_finder(
            adata,
            resolutions=[0.1],
            method="t-test",  # type: ignore[arg-type]
        )


# Test 7: Bounds on returned values (n_top_genes)
@pytest.mark.parametrize("n_top_genes", [1, 3])
def test_cluster_resolution_finder_n_top_genes(adata_for_test, n_top_genes):
    """Test that n_top_genes bounds the number of genes stored in adata.uns."""
    adata = adata_for_test.copy()
    resolutions = [0.1, 0.5]
    result = cluster_resolution_finder(
        adata,
        resolutions,
        n_top_genes=n_top_genes,
    )

    # Check the number of genes in adata.uns["cluster_resolution_top_genes"]
    top_genes_dict = result.uns["cluster_resolution_top_genes"]
    for genes in top_genes_dict.values():
        assert len(genes) <= n_top_genes
