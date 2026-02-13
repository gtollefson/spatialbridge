"""
Foundational loading and efficient access for Xenium spatial data.
Single entry point for downstream pipeline: load, subset, query by gene or region.
"""
from pathlib import Path
from typing import Optional, Sequence, Union

import numpy as np
from anndata import AnnData

from src.config import V1_BUNDLE, PRIME_BUNDLE, V1_EXTRACTED, PRIME_EXTRACTED
from src.load_xenium import ensure_extracted, load_xenium_minimal, load_xenium_genes_subset

DATASET_NAMES = ("v1", "prime")
_BUNDLES = {"v1": V1_BUNDLE, "prime": PRIME_BUNDLE}
_EXTRACTED = {"v1": V1_EXTRACTED, "prime": PRIME_EXTRACTED}


def get_dataset_path(name: str) -> Path:
    """Resolve path to extracted Xenium output dir; extract from zip if needed."""
    if name not in DATASET_NAMES:
        raise ValueError(f"name must be one of {DATASET_NAMES}")
    bundle = Path(V1_BUNDLE if name == "v1" else PRIME_BUNDLE)
    out_dir = Path(V1_EXTRACTED if name == "v1" else PRIME_EXTRACTED)
    return ensure_extracted(bundle, out_dir)


def load_dataset(
    name: str,
    genes: Optional[Sequence[str]] = None,
    n_cells: Optional[int] = None,
    use_gene_subset_io: bool = True,
) -> AnnData:
    """
    Load v1 or Prime spatial data. Efficient when only a subset of genes is needed.

    Parameters
    ----------
    name : 'v1' | 'prime'
    genes : optional list of gene names; if set and use_gene_subset_io True,
        only these genes are read from disk (saves memory for large panels).
    n_cells : optional; if set, subsample to this many cells after load (for quick tests).
    use_gene_subset_io : if True and genes is provided, read only those genes from H5.

    Returns
    -------
    AnnData with X (cells x genes), obs, obsm['spatial'].
    """
    path = get_dataset_path(name)
    if genes is not None and use_gene_subset_io:
        try:
            adata = load_xenium_genes_subset(path, list(genes))
        except Exception:
            adata = load_xenium_minimal(path)
            adata = adata[:, [g for g in genes if g in adata.var_names]].copy()
    else:
        adata = load_xenium_minimal(path)
        if genes is not None:
            adata = adata[:, [g for g in genes if g in adata.var_names]].copy()
    if n_cells is not None and adata.n_obs > n_cells:
        rng = np.random.default_rng(0)
        idx = rng.choice(adata.n_obs, size=n_cells, replace=False)
        adata = adata[idx].copy()
    return adata


def get_expression(
    adata: AnnData,
    genes: Union[str, Sequence[str]],
    layer: Optional[str] = None,
) -> np.ndarray:
    """
    Return expression matrix for given gene(s). (n_cells,) for one gene, (n_cells, n_genes) for multiple.
    """
    if isinstance(genes, str):
        genes = [genes]
    genes = [g for g in genes if g in adata.var_names]
    if not genes:
        return np.full((adata.n_obs, 0), np.nan)
    if layer and layer in adata.layers:
        X = adata[:, genes].layers[layer]
    else:
        X = adata[:, genes].X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X, dtype=float)
    return X.ravel() if X.shape[1] == 1 else X


def get_spatial(adata: AnnData) -> np.ndarray:
    """Return (n_cells, 2) array of x, y coordinates in microns."""
    return np.asarray(adata.obsm["spatial"], dtype=float)


def cells_in_bbox(
    adata: AnnData,
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
) -> np.ndarray:
    """Boolean mask of cells whose centroid lies inside the given bounding box."""
    xy = get_spatial(adata)
    return (
        (xy[:, 0] >= x_min) & (xy[:, 0] <= x_max) &
        (xy[:, 1] >= y_min) & (xy[:, 1] <= y_max)
    )


def subset_by_bbox(
    adata: AnnData,
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
) -> AnnData:
    """Return subset of adata for cells inside the bounding box."""
    mask = cells_in_bbox(adata, x_min, x_max, y_min, y_max)
    return adata[mask].copy()


def shared_genes(adata_v1: AnnData, adata_prime: AnnData) -> list:
    """Sorted list of gene names present in both datasets."""
    out = sorted(set(adata_v1.var_names) & set(adata_prime.var_names))
    return out


def prime_only_genes(adata_v1: AnnData, adata_prime: AnnData) -> list:
    """Sorted list of gene names in Prime but not in v1."""
    out = sorted(set(adata_prime.var_names) - set(adata_v1.var_names))
    return out


def expression_weighted_centroid(
    adata: AnnData,
    gene: str,
    layer: Optional[str] = None,
) -> np.ndarray:
    """(x, y) centroid of cell positions weighted by expression of the given gene."""
    v = get_expression(adata, gene, layer=layer)
    v = np.maximum(np.asarray(v, float).ravel(), 0)
    xy = get_spatial(adata)
    if v.sum() <= 0:
        return np.nanmean(xy, axis=0)
    return np.average(xy, axis=0, weights=v)
