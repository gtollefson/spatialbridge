"""Spatial registration: PAX8/PTPRC density centroids, translation of Prime to align with v1."""
from pathlib import Path
from typing import Optional
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from anndata import AnnData


def _density_centroid(adata: AnnData, gene: str, layer: str = "normalized") -> np.ndarray:
    """(x, y) centroid of expression-weighted cell positions."""
    if gene not in adata.var_names:
        return np.array([np.nan, np.nan])
    idx = adata.var_names.get_loc(gene)
    if layer in adata.layers:
        w = np.asarray(adata.layers[layer][:, idx]).ravel()
    else:
        w = np.asarray(adata.X[:, idx]).ravel()
    w = np.maximum(w, 0)
    xy = adata.obsm["spatial"]
    if w.sum() <= 0:
        return np.nanmean(xy, axis=0)
    return np.average(xy, axis=0, weights=w)


def align_sections(
    adata_v1: AnnData,
    adata_prime: AnnData,
    anchor_genes: tuple[str, str] = ("PAX8", "PTPRC"),
    layer: str = "normalized",
) -> tuple[AnnData, np.ndarray]:
    """
    Compute translation (dx, dy) from Prime density centroids to v1 for anchor genes;
    apply to prime.obsm['spatial']. Return (adata_prime_aligned, translation_vector).
    """
    prime = adata_prime.copy()
    trans = np.zeros(2)
    for gene in anchor_genes:
        c_v1 = _density_centroid(adata_v1, gene, layer)
        c_prime = _density_centroid(prime, gene, layer)
        if np.any(np.isnan(c_v1)) or np.any(np.isnan(c_prime)):
            continue
        trans += (c_v1 - c_prime)
    trans /= len(anchor_genes)
    prime.obsm["spatial"] = np.asarray(prime.obsm["spatial"], dtype=float) + trans
    return prime, trans


def plot_alignment_check(
    adata_v1: AnnData,
    adata_prime_aligned: AnnData,
    anchor_genes: tuple[str, str] = ("PAX8", "PTPRC"),
    layer: str = "normalized",
    out_path: Optional[Path] = None,
) -> None:
    """2-panel or overlay of density contours / scatter for anchor genes (v1 vs aligned Prime)."""
    out_path = Path(out_path) if out_path is not None else None
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, gene in zip(axes, anchor_genes):
        if gene not in adata_v1.var_names:
            ax.set_title(f"{gene} (v1 only)")
            continue
        idx_v1 = adata_v1.var_names.get_loc(gene)
        if layer in adata_v1.layers:
            v1_vals = np.asarray(adata_v1.layers[layer][:, idx_v1]).ravel()
        else:
            v1_vals = np.asarray(adata_v1.X[:, idx_v1]).ravel()
        xy_v1 = adata_v1.obsm["spatial"]
        ax.scatter(xy_v1[:, 0], xy_v1[:, 1], c=v1_vals, s=0.5, cmap="viridis", alpha=0.7, label="v1")
        if gene in adata_prime_aligned.var_names:
            idx_p = adata_prime_aligned.var_names.get_loc(gene)
            if layer in adata_prime_aligned.layers:
                p_vals = np.asarray(adata_prime_aligned.layers[layer][:, idx_p]).ravel()
            else:
                p_vals = np.asarray(adata_prime_aligned.X[:, idx_p]).ravel()
            xy_p = adata_prime_aligned.obsm["spatial"]
            ax.scatter(xy_p[:, 0], xy_p[:, 1], c=p_vals, s=0.5, cmap="magma", alpha=0.5, label="Prime")
        ax.set_title(gene)
        ax.set_aspect("equal")
        ax.legend(loc="upper right", fontsize=8)
    plt.tight_layout()
    if out_path:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
