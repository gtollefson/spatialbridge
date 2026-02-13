"""
Discovery utilities: compute a niche enrichment score for "virtual stains".

Definition (MVP):
- Define a tumor mask using an anchor gene (default: PAX8) in v1, using a high-expression quantile.
- Define tumor boundary cells as tumor cells that neighbor at least one non-tumor cell (kNN in XY space).
- Compute each cell's distance to the nearest boundary cell.
- Define a boundary band (e.g. 50 µm) and compute enrichment of target genes inside vs outside the band.

Enrichment score:
  log2( (mean_in_band + eps) / (mean_outside + eps) )

Optionally compute a permutation p-value by shuffling expression across cells.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.spatial import cKDTree


@dataclass(frozen=True)
class NicheEnrichmentParams:
    boundary_gene: str = "PAX8"
    tumor_quantile: float = 0.90
    band_um: float = 50.0
    k_neighbors: int = 10
    n_permutations: int = 200
    eps: float = 1e-9
    boundary_layer: str = "normalized"
    expression_layer: str = "imputed_5k"


def _get_vector(adata: AnnData, gene: str, layer: str) -> np.ndarray:
    if gene not in adata.var_names:
        raise KeyError(f"Gene '{gene}' not in adata.var_names")
    idx = adata.var_names.get_loc(gene)
    arr = adata.layers[layer] if layer in adata.layers else adata.X
    v = np.asarray(arr[:, idx]).ravel()
    if hasattr(v, "toarray"):
        v = v.toarray().ravel()
    return np.asarray(v, dtype=float)


def _tumor_mask_from_quantile(expr: np.ndarray, q: float) -> Tuple[np.ndarray, float]:
    if expr.size == 0:
        return np.zeros(0, dtype=bool), np.nan
    thr = float(np.quantile(expr, q))
    return expr >= thr, thr


def _boundary_mask_knn(xy: np.ndarray, tumor_mask: np.ndarray, k: int) -> np.ndarray:
    if xy.shape[0] == 0:
        return np.zeros(0, dtype=bool)
    if tumor_mask.sum() == 0 or tumor_mask.sum() == xy.shape[0]:
        return np.zeros(xy.shape[0], dtype=bool)
    k_eff = int(min(max(2, k + 1), xy.shape[0]))  # +1 to include self
    tree = cKDTree(xy)
    _, nn = tree.query(xy, k=k_eff)
    # Ensure 2D for k=1 cases
    if nn.ndim == 1:
        nn = nn[:, None]
    is_tumor = tumor_mask[nn]
    # boundary: tumor cell with at least one neighbor that is not tumor
    boundary = tumor_mask & (~np.all(is_tumor, axis=1))
    return boundary


def _distance_to_boundary(xy: np.ndarray, boundary_mask: np.ndarray) -> np.ndarray:
    if xy.shape[0] == 0:
        return np.zeros(0, dtype=float)
    if boundary_mask.sum() == 0:
        return np.full(xy.shape[0], np.nan, dtype=float)
    bxy = xy[boundary_mask]
    tree = cKDTree(bxy)
    dist, _ = tree.query(xy, k=1)
    return np.asarray(dist, dtype=float)


def compute_niche_enrichment(
    adata_v1: AnnData,
    genes: Iterable[str],
    params: Optional[NicheEnrichmentParams] = None,
    out_tsv: Optional[Path] = None,
    seed: int = 0,
) -> pd.DataFrame:
    """
    Compute niche enrichment scores for given genes in `adata_v1`.

    Requires:
    - `adata_v1.obsm['spatial']` (micron coordinates)
    - `params.boundary_layer` available for `params.boundary_gene` (typically 'normalized')
    - `params.expression_layer` available for scoring genes (typically 'imputed_5k')
    """
    params = params or NicheEnrichmentParams()
    genes = [g for g in genes if g in adata_v1.var_names]
    if not genes:
        return pd.DataFrame()

    xy = np.asarray(adata_v1.obsm["spatial"], dtype=float)
    boundary_expr = _get_vector(adata_v1, params.boundary_gene, params.boundary_layer)
    tumor_mask, tumor_thr = _tumor_mask_from_quantile(boundary_expr, params.tumor_quantile)
    boundary_mask = _boundary_mask_knn(xy, tumor_mask, params.k_neighbors)
    dist_to_boundary = _distance_to_boundary(xy, boundary_mask)
    niche_mask = dist_to_boundary <= float(params.band_um)

    rng = np.random.default_rng(seed)

    rows: List[Dict] = []
    for g in genes:
        v = _get_vector(adata_v1, g, params.expression_layer)
        in_band = v[niche_mask]
        out_band = v[~niche_mask]
        mean_in = float(np.nanmean(in_band)) if in_band.size else np.nan
        mean_out = float(np.nanmean(out_band)) if out_band.size else np.nan
        score = float(np.log2((mean_in + params.eps) / (mean_out + params.eps))) if np.isfinite(mean_in) and np.isfinite(mean_out) else np.nan

        p_perm = np.nan
        if params.n_permutations and np.isfinite(score):
            perm_scores = np.zeros(int(params.n_permutations), dtype=float)
            for i in range(int(params.n_permutations)):
                vp = rng.permutation(v)
                mi = float(np.nanmean(vp[niche_mask])) if niche_mask.any() else np.nan
                mo = float(np.nanmean(vp[~niche_mask])) if (~niche_mask).any() else np.nan
                perm_scores[i] = float(np.log2((mi + params.eps) / (mo + params.eps)))
            # two-sided permutation p-value
            p_perm = float((np.sum(np.abs(perm_scores) >= abs(score)) + 1.0) / (perm_scores.size + 1.0))

        rows.append(
            {
                "gene": g,
                "n_cells": int(adata_v1.n_obs),
                "boundary_gene": params.boundary_gene,
                "tumor_quantile": float(params.tumor_quantile),
                "tumor_threshold": float(tumor_thr),
                "band_um": float(params.band_um),
                "k_neighbors": int(params.k_neighbors),
                "n_boundary_cells": int(boundary_mask.sum()),
                "n_band_cells": int(np.sum(niche_mask)),
                "mean_in_band": mean_in,
                "mean_outside_band": mean_out,
                "log2_enrichment": score,
                "p_perm_two_sided": p_perm,
                "n_permutations": int(params.n_permutations),
            }
        )

    df = pd.DataFrame(rows).sort_values("log2_enrichment", ascending=False)
    if out_tsv is not None:
        out_tsv = Path(out_tsv)
        out_tsv.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_tsv, sep="\t", index=False)
    return df


def format_niche_summary(df: pd.DataFrame, top_n: int = 5) -> List[str]:
    """Format a short multi-line summary suitable for appending to results_summary.txt."""
    if df is None or df.empty:
        return ["Niche enrichment: not computed (no genes scored)."]
    lines = ["Niche enrichment (boundary band):"]
    show = df.head(int(top_n))
    for _, r in show.iterrows():
        gene = r["gene"]
        score = r["log2_enrichment"]
        p = r.get("p_perm_two_sided", np.nan)
        lines.append(f"- {gene}: log2_enrichment={score:.4f}, p_perm={p:.4g}")
    return lines

