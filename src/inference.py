"""KNN-based imputation: train on Prime (shared -> Prime-only), predict on v1.
Memory-efficient: one NearestNeighbors fit, then stream over genes (no per-gene model storage).
"""
import numpy as np
from sklearn.neighbors import NearestNeighbors
from anndata import AnnData

_EPS = 1e-10


def _to_dense(X):
    if hasattr(X, "toarray"):
        return X.toarray()
    return np.asarray(X, dtype=float)


def _knn_weighted_predict(neighbor_indices: np.ndarray, neighbor_distances: np.ndarray, y_prime: np.ndarray) -> np.ndarray:
    """Predict for each row as distance-weighted mean of y_prime at neighbor indices. Same as KNeighborsRegressor(weights='distance')."""
    # neighbor_indices (n_v1, k), neighbor_distances (n_v1, k), y_prime (n_prime,)
    n_v1, k = neighbor_indices.shape
    weights = 1.0 / (neighbor_distances + _EPS)
    y_at_neighbors = y_prime[neighbor_indices]
    return (weights * y_at_neighbors).sum(axis=1) / weights.sum(axis=1)


def add_imputed_layer(
    adata_v1: AnnData,
    adata_prime: AnnData,
    shared_genes: list[str],
    prime_only_genes: list[str],
    layer: str = "normalized",
    k: int = 15,
    weights: str = "distance",
    imputed_layer_name: str = "imputed_5k",
) -> AnnData:
    """
    Impute Prime-only genes for v1 using one shared KNN in Prime shared-gene space.
    Memory-efficient: single NearestNeighbors fit, then stream over genes.
    """
    X_prime = _to_dense(adata_prime[:, shared_genes].layers[layer] if layer in adata_prime.layers else adata_prime[:, shared_genes].X)
    X_v1 = _to_dense(adata_v1[:, shared_genes].layers[layer] if layer in adata_v1.layers else adata_v1[:, shared_genes].X)

    nn = NearestNeighbors(n_neighbors=k, metric="euclidean")
    nn.fit(X_prime)
    neighbor_distances, neighbor_indices = nn.kneighbors(X_v1)

    n_prime_only = len(prime_only_genes)
    imputed = np.zeros((adata_v1.n_obs, n_prime_only), dtype=float)
    for j, g in enumerate(prime_only_genes):
        if g not in adata_prime.var_names:
            continue
        idx = adata_prime.var_names.get_loc(g)
        y_prime = _to_dense(adata_prime.layers[layer][:, idx] if layer in adata_prime.layers else adata_prime.X[:, idx]).ravel()
        imputed[:, j] = _knn_weighted_predict(neighbor_indices, neighbor_distances, y_prime)

    shared_arr = _to_dense(adata_v1.layers[layer] if layer in adata_v1.layers else adata_v1.X)
    full = np.hstack([shared_arr, imputed])
    all_genes = shared_genes + prime_only_genes
    import pandas as pd
    expanded = AnnData(
        full,
        obs=adata_v1.obs,
        var=pd.DataFrame(index=all_genes),
        obsm=dict(adata_v1.obsm),
    )
    expanded.layers[imputed_layer_name] = full
    if layer in adata_v1.layers:
        shared_layer = _to_dense(adata_v1.layers[layer])
        expanded.layers[layer] = np.hstack([shared_layer, np.zeros((adata_v1.n_obs, n_prime_only), dtype=float)])
    return expanded
