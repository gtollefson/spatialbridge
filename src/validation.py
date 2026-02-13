"""Holdout validation (Pearson, MSE) and virtual stain 2x2 plot."""
from pathlib import Path
from typing import List, Optional
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from sklearn.neighbors import NearestNeighbors
from anndata import AnnData

_EPS = 1e-10


def holdout_validation(
    adata_v1: AnnData,
    adata_prime: AnnData,
    shared_genes: list[str],
    n_holdout: int = 50,
    layer: str = "normalized",
    k: int = 15,
    weights: str = "distance",
    seed: int = 42,
) -> tuple[float, float, np.ndarray, np.ndarray]:
    """
    Hold out n_holdout shared genes; one NN in Prime (train_genes) space, predict holdout via weighted mean.
    Returns (pearson_r, mse, real_flat, imputed_flat). Memory-efficient (no per-gene models).
    """
    rng = np.random.default_rng(seed)
    n_shared = len(shared_genes)
    if n_holdout >= n_shared:
        n_holdout = max(1, n_shared // 2)
    holdout_idx = rng.choice(n_shared, size=n_holdout, replace=False)
    holdout_genes = [shared_genes[i] for i in holdout_idx]
    train_genes = [g for i, g in enumerate(shared_genes) if i not in set(holdout_idx)]

    def _to_dense(X):
        if hasattr(X, "toarray"):
            return X.toarray()
        return np.asarray(X, dtype=float)

    X_prime = _to_dense(adata_prime[:, train_genes].layers[layer] if layer in adata_prime.layers else adata_prime[:, train_genes].X)
    X_v1 = _to_dense(adata_v1[:, train_genes].layers[layer] if layer in adata_v1.layers else adata_v1[:, train_genes].X)

    nn = NearestNeighbors(n_neighbors=k, metric="euclidean")
    nn.fit(X_prime)
    neighbor_distances, neighbor_indices = nn.kneighbors(X_v1)
    w = 1.0 / (neighbor_distances + _EPS)

    imputed_list = []
    for g in holdout_genes:
        idx = adata_prime.var_names.get_loc(g)
        y = _to_dense(adata_prime.layers[layer][:, idx] if layer in adata_prime.layers else adata_prime.X[:, idx]).ravel()
        y_at_nn = y[neighbor_indices]
        imputed_list.append((w * y_at_nn).sum(axis=1) / w.sum(axis=1))
    imputed_vals = np.column_stack(imputed_list).ravel()

    real_list = []
    for g in holdout_genes:
        idx = adata_v1.var_names.get_loc(g)
        r = _to_dense(adata_v1.layers[layer][:, idx] if layer in adata_v1.layers else adata_v1.X[:, idx]).ravel()
        real_list.append(r)
    real_vals = np.column_stack(real_list).ravel()

    valid = np.isfinite(real_vals) & np.isfinite(imputed_vals)
    real_vals = real_vals[valid]
    imputed_vals = imputed_vals[valid]
    if real_vals.size < 2:
        return 0.0, np.nan, real_vals, imputed_vals
    r, _ = pearsonr(real_vals, imputed_vals)
    mse = np.mean((real_vals - imputed_vals) ** 2)
    return float(r), float(mse), real_vals, imputed_vals


def plot_validation_scatter(
    real: np.ndarray,
    imputed: np.ndarray,
    pearson_r: float,
    out_path: Optional[Path] = None,
) -> None:
    """Scatter real vs imputed with regression line and r in title."""
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(real, imputed, s=1, alpha=0.5)
    if real.size > 1:
        z = np.polyfit(real, imputed, 1)
        xl = np.array([real.min(), real.max()])
        ax.plot(xl, np.polyval(z, xl), "r-", lw=2, label=f"r = {pearson_r:.4f}")
    ax.set_xlabel("Real (held-out)")
    ax.set_ylabel("Imputed")
    ax.set_title(f"Validation: Pearson r = {pearson_r:.4f}")
    ax.legend()
    ax.set_aspect("equal")
    if out_path:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()


def _get_expression(adata: AnnData, gene: str, layer: str, imputed_layer: str, imputed_gene_idx: Optional[int]) -> np.ndarray:
    """Get per-cell expression for a gene from raw, normalized, or imputed layer."""
    if gene in adata.var_names:
        idx = adata.var_names.get_loc(gene)
        arr = adata.layers[layer] if layer in adata.layers else adata.X
        v = np.asarray(arr[:, idx]).ravel()
        if hasattr(v, "toarray"):
            v = v.toarray().ravel()
        return v
    if imputed_gene_idx is not None and imputed_layer in adata.layers:
        return np.asarray(adata.layers[imputed_layer][:, imputed_gene_idx]).ravel()
    return np.zeros(adata.n_obs)


def plot_virtual_stain(
    adata_v1: AnnData,
    genes: tuple[str, str, str, str] = ("PAX8", "PTPRC", "IDO1", "LAG3"),
    layer: str = "normalized",
    imputed_layer: str = "imputed_5k",
    shared_genes: Optional[List[str]] = None,
    prime_only_genes: Optional[List[str]] = None,
    out_path: Optional[Path] = None,
) -> None:
    """
    2x2 spatial scatter: (a) Real PAX8, (b) Real PTPRC, (c) Imputed IDO1, (d) Imputed LAG3.
    For imputed genes, index into imputed_layer is n_shared + index_in_prime_only.
    """
    if shared_genes is None:
        shared_genes = []
    if prime_only_genes is None:
        prime_only_genes = []
    n_shared = len(shared_genes)
    fig, axes = plt.subplots(2, 2, figsize=(10, 10))
    axes = axes.ravel()
    for ax, gene in zip(axes, genes):
        imputed_idx = None
        if gene not in adata_v1.var_names and gene in prime_only_genes:
            imputed_idx = n_shared + prime_only_genes.index(gene)
        vals = _get_expression(adata_v1, gene, layer, imputed_layer, imputed_idx)
        xy = adata_v1.obsm["spatial"]
        label = "Imputed " + gene if imputed_idx is not None else "Real " + gene
        ax.scatter(xy[:, 0], xy[:, 1], c=vals, s=0.5, cmap="viridis")
        ax.set_title(label)
        ax.set_aspect("equal")
    plt.tight_layout()
    if out_path:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()


def write_results_summary(
    out_path: Path,
    n_shared: int,
    n_prime_only: int,
    pearson_r: float,
    mse: float,
    translation: np.ndarray,
    concordance_target: Optional[float] = 0.8,
) -> None:
    """Write results_summary.txt with key metrics. Concordance = Pearson r on held-out shared genes."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        f.write("SpatialBridge results summary\n")
        f.write("============================\n")
        f.write(f"Shared genes: {n_shared}\n")
        f.write(f"Prime-only genes imputed: {n_prime_only}\n")
        f.write(f"Shared-gene concordance (holdout): Pearson r = {pearson_r:.4f}\n")
        f.write(f"Validation MSE: {mse:.4f}\n")
        if concordance_target is not None:
            met = "met" if pearson_r >= concordance_target else "not met"
            f.write(f"Target r >= {concordance_target:.2f}: {met}\n")
        f.write(f"Alignment translation (dx, dy): ({translation[0]:.2f}, {translation[1]:.2f})\n")
