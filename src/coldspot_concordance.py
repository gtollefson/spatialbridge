"""
Cold-spot concordance between Prime (real 5K) and v1 (imputed 5K).
Zones: tumor boundary band, distance from vasculature, necrotic cores.
Outputs: side-by-side cold-spot maps, side-by-side by context, zone-level concordance scatter.
"""
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from anndata import AnnData
from scipy.spatial import cKDTree
from scipy import stats

from src.discovery import (
    _get_vector,
    _tumor_mask_from_quantile,
    _boundary_mask_knn,
    _distance_to_boundary,
)


def _to_dense(arr):
    if hasattr(arr, "toarray"):
        return arr.toarray().ravel()
    return np.asarray(arr, dtype=float).ravel()


def cold_spot_score(
    adata: AnnData,
    immune_gene: str,
    layer: str = "normalized",
    imputed_layer: str = "imputed_5k",
) -> np.ndarray:
    """Per-cell cold score: immune-only. High = low immune = cold. score = -z(immune)."""
    if immune_gene not in adata.var_names:
        return np.full(adata.n_obs, np.nan)
    arr_i = adata.layers[imputed_layer] if imputed_layer in adata.layers else (adata.layers[layer] if layer in adata.layers else adata.X)
    idx_i = adata.var_names.get_loc(immune_gene)
    i = _to_dense(arr_i[:, idx_i])
    i = np.asarray(i, dtype=float)
    z_i = (i - np.nanmean(i)) / (np.nanstd(i) + 1e-10)
    return -z_i  # high score = low immune = cold


def cold_spot_score_aggregate(
    adata: AnnData,
    genes: List[str],
    layer: str = "normalized",
    imputed_layer: str = "imputed_5k",
) -> Tuple[np.ndarray, int]:
    """Per-cell cold score = -mean(z) over genes (high = low immune = cold). Returns (score, n_genes_used)."""
    genes = [g for g in genes if g in adata.var_names]
    if not genes:
        return np.full(adata.n_obs, np.nan), 0
    arr = adata.layers[imputed_layer] if imputed_layer in adata.layers else (adata.layers[layer] if layer in adata.layers else adata.X)
    neg_z_list = []
    for g in genes:
        idx = adata.var_names.get_loc(g)
        x = _to_dense(arr[:, idx])
        x = np.asarray(x, dtype=float)
        z = (x - np.nanmean(x)) / (np.nanstd(x) + 1e-10)
        neg_z_list.append(-z)
    score = np.nanmean(np.stack(neg_z_list, axis=0), axis=0)
    return score, len(genes)


def boundary_band_mask(
    xy: np.ndarray,
    tumor_expr: np.ndarray,
    q: float = 0.9,
    k: int = 10,
    band_um: float = 50.0,
) -> np.ndarray:
    """Boolean mask: cells within band_um of tumor boundary."""
    tumor_mask, _ = _tumor_mask_from_quantile(tumor_expr, q)
    boundary = _boundary_mask_knn(xy, tumor_mask, k)
    dist = _distance_to_boundary(xy, boundary)
    return (dist <= band_um) & np.isfinite(dist)


def vessel_distance_bins(
    xy: np.ndarray,
    vessel_expr: np.ndarray,
    bins_um: Tuple[float, ...],
) -> np.ndarray:
    """Per-cell distance to nearest vessel (vessel = top quantile of vessel_expr). Return bin index 0..len(bins)-1, or -1 if beyond last bin."""
    q = 0.90
    vessel_mask = vessel_expr >= np.nanquantile(vessel_expr, q)
    if vessel_mask.sum() == 0:
        return np.full(xy.shape[0], -1, dtype=int)
    vxy = xy[vessel_mask]
    tree = cKDTree(vxy)
    dist, _ = tree.query(xy, k=1)
    dist = np.asarray(dist, dtype=float).ravel()
    out = np.full(xy.shape[0], -1, dtype=int)
    for i in range(len(bins_um) - 1):
        out[(dist >= bins_um[i]) & (dist <= bins_um[i + 1])] = i
    return out


def necrosis_mask_from_low_rna(adata: AnnData, layer: str = "normalized", bottom_quantile: float = 0.02) -> np.ndarray:
    """Cells in bottom total-expression quantile as proxy for necrotic/low viability."""
    if layer in adata.layers:
        X = adata.layers[layer]
    else:
        X = adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    total = np.asarray(X.sum(axis=1)).ravel()
    thr = np.nanquantile(total, bottom_quantile)
    return total <= thr


def compute_zone_means(
    score: np.ndarray,
    zone_labels: List[Tuple[str, np.ndarray]],
) -> List[Tuple[str, float]]:
    """For each (name, boolean mask), return (name, mean score)."""
    out = []
    for name, mask in zone_labels:
        if mask.sum() == 0:
            out.append((name, np.nan))
        else:
            out.append((name, float(np.nanmean(score[mask]))))
    return out


def concordance_ccc(x: np.ndarray, y: np.ndarray) -> float:
    """Lin's concordance correlation coefficient."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    if valid.sum() < 2:
        return np.nan
    x, y = x[valid], y[valid]
    r = np.corrcoef(x, y)[0, 1]
    sx, sy = np.std(x), np.std(y)
    mx, my = np.mean(x), np.mean(y)
    return float(2 * r * sx * sy / (sx**2 + sy**2 + (mx - my)**2))


def run_concordance(
    adata_prime: AnnData,
    adata_v1: AnnData,
    tumor_gene: str = "PAX8",
    immune_gene: str = "CD8A",
    native_immune_gene: Optional[str] = "PTPRC",
    imputed_only_immune_gene: Optional[str] = "LAG3",
    shared_immune_genes: Optional[List[str]] = None,
    prime_only_immune_genes: Optional[List[str]] = None,
    vessel_gene: Optional[str] = "PECAM1",
    necrosis_gene: Optional[str] = None,
    boundary_band_um: float = 50.0,
    vessel_bins_um: Tuple[float, ...] = (0, 50, 100, 200, 500),
    layer: str = "normalized",
    imputed_layer: str = "imputed_5k",
    out_dir: Optional[Path] = None,
) -> pd.DataFrame:
    """
    Compute cold-spot scores in Prime and v1; define zones (boundary, vessel bins, necrosis);
    produce side-by-side maps, context panels, and concordance scatter. Return zone-level table.
    """
    out_dir = Path(out_dir) if out_dir else Path("results")
    out_dir.mkdir(parents=True, exist_ok=True)

    xy_prime = np.asarray(adata_prime.obsm["spatial"], dtype=float)
    xy_v1 = np.asarray(adata_v1.obsm["spatial"], dtype=float)
    t_prime = _get_vector(adata_prime, tumor_gene, layer)
    t_v1 = _get_vector(adata_v1, tumor_gene, layer)

    # Cold scores: aggregate (two lists) or single-gene fallback
    use_aggregate = bool(shared_immune_genes) and bool(prime_only_immune_genes)
    n_shared_used = 0
    n_prime_only_used = 0
    gene_imputed = imputed_only_immune_gene if (imputed_only_immune_gene and imputed_only_immune_gene in adata_prime.var_names and imputed_only_immune_gene in adata_v1.var_names) else immune_gene
    immune_native = native_immune_gene if (native_immune_gene and native_immune_gene in adata_v1.var_names) else immune_gene

    if use_aggregate:
        genes_native = [g for g in shared_immune_genes if g in adata_v1.var_names]
        genes_prime = [g for g in prime_only_immune_genes if g in adata_prime.var_names]
        genes_v1_imputed = [g for g in prime_only_immune_genes if g in adata_v1.var_names]
        if genes_native:
            score_v1_native, n_shared_used = cold_spot_score_aggregate(adata_v1, genes_native, layer=layer, imputed_layer=layer)
        else:
            score_v1_native = cold_spot_score(adata_v1, immune_native, layer=layer, imputed_layer=layer)
        if genes_prime and genes_v1_imputed:
            score_prime, n_prime_only_used = cold_spot_score_aggregate(adata_prime, genes_prime, layer=layer, imputed_layer=layer)
            score_v1, _ = cold_spot_score_aggregate(adata_v1, genes_v1_imputed, layer=layer, imputed_layer=imputed_layer)
        else:
            score_prime = cold_spot_score(adata_prime, gene_imputed, layer=layer, imputed_layer=layer)
            score_v1 = cold_spot_score(adata_v1, gene_imputed, layer=layer, imputed_layer=imputed_layer)
        if not genes_native and (not genes_prime or not genes_v1_imputed):
            use_aggregate = False
    if not use_aggregate:
        score_prime = cold_spot_score(adata_prime, gene_imputed, layer=layer, imputed_layer=layer)
        score_v1 = cold_spot_score(adata_v1, gene_imputed, layer=layer, imputed_layer=imputed_layer)
        score_v1_native = cold_spot_score(adata_v1, immune_native, layer=layer, imputed_layer=layer)

    # Tumor masks (PAX8 high) for visualization: color only tumor cells by cold score; background = light grey
    tumor_p, _ = _tumor_mask_from_quantile(t_prime, 0.9)
    tumor_v1, _ = _tumor_mask_from_quantile(t_v1, 0.9)

    v_prime = v_v1 = None
    if vessel_gene and vessel_gene in adata_prime.var_names and vessel_gene in adata_v1.var_names:
        v_prime = _get_vector(adata_prime, vessel_gene, layer)
        arr_v1 = adata_v1.layers[imputed_layer] if imputed_layer in adata_v1.layers else adata_v1.X
        v_v1 = _to_dense(arr_v1[:, adata_v1.var_names.get_loc(vessel_gene)])

    zone_names: List[str] = []
    mean_prime: List[float] = []
    mean_v1: List[float] = []
    n_cells_prime: List[int] = []
    n_cells_v1: List[int] = []

    # Boundary band
    band_p = boundary_band_mask(xy_prime, t_prime, q=0.9, k=10, band_um=boundary_band_um)
    band_v1 = boundary_band_mask(xy_v1, t_v1, q=0.9, k=10, band_um=boundary_band_um)
    zone_names.append("Tumor boundary band")
    mean_prime.append(float(np.nanmean(score_prime[band_p])) if band_p.sum() else np.nan)
    mean_v1.append(float(np.nanmean(score_v1[band_v1])) if band_v1.sum() else np.nan)
    n_cells_prime.append(int(band_p.sum()))
    n_cells_v1.append(int(band_v1.sum()))

    # Vessel distance bins
    if v_prime is not None and v_v1 is not None:
        bins_p = vessel_distance_bins(xy_prime, v_prime, vessel_bins_um)
        bins_v1 = vessel_distance_bins(xy_v1, v_v1, vessel_bins_um)
        for i in range(len(vessel_bins_um) - 1):
            label = f"Vessel {vessel_bins_um[i]}-{vessel_bins_um[i+1]} µm"
            zone_names.append(label)
            mp = score_prime[bins_p == i]
            mv = score_v1[bins_v1 == i]
            mean_prime.append(float(np.nanmean(mp)) if len(mp) else np.nan)
            mean_v1.append(float(np.nanmean(mv)) if len(mv) else np.nan)
            n_cells_prime.append(int((bins_p == i).sum()))
            n_cells_v1.append(int((bins_v1 == i).sum()))

    # Necrosis: low-RNA proxy or gene
    if necrosis_gene and necrosis_gene in adata_prime.var_names and necrosis_gene in adata_v1.var_names:
        nec_p = _get_vector(adata_prime, necrosis_gene, layer) >= np.nanquantile(
            _get_vector(adata_prime, necrosis_gene, layer), 0.98
        )
        arr_n = adata_v1.layers[imputed_layer] if imputed_layer in adata_v1.layers else adata_v1.X
        nv = _to_dense(arr_n[:, adata_v1.var_names.get_loc(necrosis_gene)])
        nec_v1 = nv >= np.nanquantile(nv, 0.98)
    else:
        nec_p = necrosis_mask_from_low_rna(adata_prime, layer)
        nec_v1 = necrosis_mask_from_low_rna(adata_v1, layer)
    zone_names.append("Necrotic / low-RNA")
    mean_prime.append(float(np.nanmean(score_prime[nec_p])) if nec_p.sum() else np.nan)
    mean_v1.append(float(np.nanmean(score_v1[nec_v1])) if nec_v1.sum() else np.nan)
    n_cells_prime.append(int(nec_p.sum()))
    n_cells_v1.append(int(nec_v1.sum()))

    df = pd.DataFrame({
        "zone": zone_names,
        "n_cells_Prime": n_cells_prime,
        "n_cells_v1_imputed": n_cells_v1,
        "mean_cold_score_Prime": mean_prime,
        "mean_cold_score_v1_imputed": mean_v1,
    })
    df.to_csv(out_dir / "coldspot_zone_means.csv", index=False)

    # Concordance stats (only over zones where both Prime and v1 have values)
    p_vals = np.array(mean_prime)
    v_vals = np.array(mean_v1)
    valid = np.isfinite(p_vals) & np.isfinite(v_vals)
    if valid.sum() >= 2:
        r, r_pvalue = stats.pearsonr(p_vals[valid], v_vals[valid])
        ccc = concordance_ccc(p_vals, v_vals)
    else:
        r, ccc, r_pvalue = np.nan, np.nan, np.nan
    concordance_summary = pd.DataFrame([{"pearson_r": r, "pearson_p": r_pvalue, "ccc": ccc, "n_zones_with_both": int(valid.sum())}])
    concordance_summary.to_csv(out_dir / "coldspot_concordance_stats.csv", index=False)

    # Short interpretation for "is imputed data useful?"
    with open(out_dir / "coldspot_interpretation.txt", "w") as f:
        f.write("Cold-spot concordance: Is imputed v1 data useful vs real 5K (Prime)?\n")
        f.write("=" * 60 + "\n")
        if use_aggregate and (n_shared_used or n_prime_only_used):
            f.write(f"Cold score = aggregate immune-only (mean of -z per gene). Native: {n_shared_used} shared genes; Prime/imputed: {n_prime_only_used} 5K-only genes.\n")
        else:
            f.write("Cold score = immune-only (high = low immune = cold). Maps show cold score in tumor (PAX8+) cells only; background = light grey.\n")
        f.write(f"Zone-level: Pearson r = {r:.3f}, p = {r_pvalue:.3f}, CCC = {ccc:.3f} (n = {valid.sum()} zones with both scores).\n")
        if valid.sum() >= 2 and np.isfinite(r_pvalue):
            f.write("Significant concordance (p < 0.05) supports that imputed data recapitulates zone-level cold-spot geography.\n")
        else:
            f.write("Limited or non-significant concordance; interpret with caution (few zones or missing values).\n")
        f.write("Blank cells in zone_means = no cells in that zone for that dataset (see n_cells_Prime, n_cells_v1_imputed).\n")
        f.write("Gene-level usefulness: see results_summary.txt (holdout Pearson r) and validation_scatter.png.\n")
        f.write("v1 improvement: coldspot_v1_native_vs_imputed_sidebyside.png (v1 without imputation vs v1 with 5K imputation).\n")

    # Boundary edge (tumor–stroma interface) for crisp outline
    edge_p = _boundary_mask_knn(xy_prime, tumor_p, 10)
    edge_v1 = _boundary_mask_knn(xy_v1, tumor_v1, 10)
    # Vessel masks (PECAM1+)
    if v_prime is not None and v_v1 is not None:
        vessel_p = v_prime >= np.nanquantile(v_prime, 0.90)
        vessel_v1 = v_v1 >= np.nanquantile(v_v1, 0.90)
    else:
        vessel_p = np.zeros(adata_prime.n_obs, dtype=bool)
        vessel_v1 = np.zeros(adata_v1.n_obs, dtype=bool)

    # Figure 1: Side-by-side cold-spot maps (immune-only; only tumor cells colored, background light grey)
    all_scores = np.concatenate([np.ravel(score_prime), np.ravel(score_v1)])
    all_scores = all_scores[np.isfinite(all_scores)]
    lim = float(np.percentile(np.abs(all_scores), 99)) if all_scores.size else 1.0
    lim = max(lim, 0.5)
    vmin, vmax = -lim, lim
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    for ax, xy, score, tumor_mask, title in [
        (ax1, xy_prime, score_prime, tumor_p, "Prime (real 5K)"),
        (ax2, xy_v1, score_v1, tumor_v1, "v1 (imputed 5K)"),
    ]:
        ax.scatter(xy[~tumor_mask, 0], xy[~tumor_mask, 1], c="#e8e8e8", s=0.25, rasterized=True)
        sc = ax.scatter(xy[tumor_mask, 0], xy[tumor_mask, 1], c=score[tumor_mask], s=0.4, cmap="coolwarm", vmin=vmin, vmax=vmax, rasterized=True)
        ax.set_title(title)
        ax.set_aspect("equal")
    plt.colorbar(sc, ax=[ax1, ax2], label="Cold spot score (low immune = high)")
    if use_aggregate and n_prime_only_used:
        fig.suptitle(f"Cold spots in tumor (PAX8+); immune-only (aggregate, {n_prime_only_used} genes)", fontsize=11)
    else:
        fig.suptitle(f"Cold spots in tumor (PAX8+); immune-only ({gene_imputed})", fontsize=11)
    plt.tight_layout()
    plt.savefig(out_dir / "coldspot_maps_sidebyside.png", dpi=150, bbox_inches="tight")
    plt.close()

    # Figure: v1 native (no 5K imputation) vs v1 imputed — show improvement from imputation (tumor only colored)
    all_v1_scores = np.concatenate([np.ravel(score_v1_native), np.ravel(score_v1)])
    all_v1_scores = all_v1_scores[np.isfinite(all_v1_scores)]
    lim_v1 = float(np.percentile(np.abs(all_v1_scores), 99)) if all_v1_scores.size else 1.0
    lim_v1 = max(lim_v1, 0.5)
    vmin_v1, vmax_v1 = -lim_v1, lim_v1
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    ax1.scatter(xy_v1[~tumor_v1, 0], xy_v1[~tumor_v1, 1], c="#e8e8e8", s=0.25, rasterized=True)
    ax1.scatter(xy_v1[tumor_v1, 0], xy_v1[tumor_v1, 1], c=score_v1_native[tumor_v1], s=0.4, cmap="coolwarm", vmin=vmin_v1, vmax=vmax_v1, rasterized=True)
    if use_aggregate and n_shared_used:
        ax1.set_title(f"v1 native (no imputation)\naggregate ({n_shared_used} shared)")
    else:
        ax1.set_title(f"v1 native (no imputation)\n{immune_native}")
    ax1.set_aspect("equal")
    ax2.scatter(xy_v1[~tumor_v1, 0], xy_v1[~tumor_v1, 1], c="#e8e8e8", s=0.25, rasterized=True)
    sc = ax2.scatter(xy_v1[tumor_v1, 0], xy_v1[tumor_v1, 1], c=score_v1[tumor_v1], s=0.4, cmap="coolwarm", vmin=vmin_v1, vmax=vmax_v1, rasterized=True)
    if use_aggregate and n_prime_only_used:
        ax2.set_title(f"v1 imputed (5K)\naggregate ({n_prime_only_used} 5K-only)")
    else:
        ax2.set_title(f"v1 imputed (5K)\n{gene_imputed} (Prime-only)")
    ax2.set_aspect("equal")
    plt.colorbar(sc, ax=[ax1, ax2], label="Cold spot score (low immune = high)")
    if use_aggregate and (n_shared_used or n_prime_only_used):
        fig.suptitle("v1 cold-spot: aggregate shared (native) vs aggregate 5K-only (imputed); similar patterns = imputation accuracy", fontsize=10)
    else:
        fig.suptitle(f"v1 cold-spot: {immune_native} (native) vs {gene_imputed} (imputed); similar patterns = imputation accuracy", fontsize=10)
    plt.tight_layout()
    plt.savefig(out_dir / "coldspot_v1_native_vs_imputed_sidebyside.png", dpi=150, bbox_inches="tight")
    plt.close()

    # Figure: Tumor boundary (v1 vs 5K) — band + interface line
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    for ax, xy, band, edge, title in [
        (ax1, xy_prime, band_p, edge_p, "Prime (real 5K)"),
        (ax2, xy_v1, band_v1, edge_v1, "v1 (imputed 5K)"),
    ]:
        ax.scatter(xy[:, 0], xy[:, 1], c="#e8e8e8", s=0.25, rasterized=True, label="All cells")
        ax.scatter(xy[band, 0], xy[band, 1], c="#0d9488", s=0.8, alpha=0.5, rasterized=True, label="Boundary band (≤50 µm)")
        ax.scatter(xy[edge, 0], xy[edge, 1], c="#115e59", s=1.2, rasterized=True, label="Tumor–stroma interface")
        ax.set_title(title)
        ax.set_aspect("equal")
        ax.legend(loc="upper right", fontsize=7)
    fig.suptitle("Tumor boundary", fontsize=12)
    plt.tight_layout()
    plt.savefig(out_dir / "coldspot_feature_tumor_boundary_sidebyside.png", dpi=150, bbox_inches="tight")
    plt.close()

    # Figure: Necrotic / low-RNA (v1 vs 5K)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    for ax, xy, nec, title in [
        (ax1, xy_prime, nec_p, "Prime (real 5K)"),
        (ax2, xy_v1, nec_v1, "v1 (imputed 5K)"),
    ]:
        ax.scatter(xy[:, 0], xy[:, 1], c="#e8e8e8", s=0.25, rasterized=True, label="All cells")
        ax.scatter(xy[nec, 0], xy[nec, 1], c="#3d2817", s=1.0, alpha=0.85, rasterized=True, label="Necrotic / low-RNA")
        ax.set_title(title)
        ax.set_aspect("equal")
        ax.legend(loc="upper right", fontsize=7)
    fig.suptitle("Necrotic regions (low-RNA proxy)", fontsize=12)
    plt.tight_layout()
    plt.savefig(out_dir / "coldspot_feature_necrosis_sidebyside.png", dpi=150, bbox_inches="tight")
    plt.close()

    # Figure: Vasculature PECAM1+ (v1 vs 5K)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    for ax, xy, vessel, title in [
        (ax1, xy_prime, vessel_p, "Prime (real 5K)"),
        (ax2, xy_v1, vessel_v1, "v1 (imputed 5K)"),
    ]:
        ax.scatter(xy[:, 0], xy[:, 1], c="#e8e8e8", s=0.25, rasterized=True, label="All cells")
        ax.scatter(xy[vessel, 0], xy[vessel, 1], c="#b91c1c", s=1.0, alpha=0.9, rasterized=True, label="Vasculature (PECAM1+)")
        ax.set_title(title)
        ax.set_aspect("equal")
        ax.legend(loc="upper right", fontsize=7)
    fig.suptitle("Vasculature (PECAM1+)", fontsize=12)
    plt.tight_layout()
    plt.savefig(out_dir / "coldspot_feature_vasculature_sidebyside.png", dpi=150, bbox_inches="tight")
    plt.close()

    # Figure 2: Side-by-side by context (bar chart: Prime vs v1 per zone)
    x = np.arange(len(zone_names))
    w = 0.35
    fig, ax = plt.subplots(figsize=(max(8, len(zone_names) * 0.8), 5))
    ax.bar(x - w/2, mean_prime, w, label="Prime (real 5K)", color="C0")
    ax.bar(x + w/2, mean_v1, w, label="v1 (imputed 5K)", color="C1")
    ax.set_xticks(x)
    ax.set_xticklabels(zone_names, rotation=45, ha="right")
    ax.set_ylabel("Mean cold spot score")
    ax.legend()
    ax.set_title("Cold spot score by zone: Prime vs imputed v1")
    plt.tight_layout()
    plt.savefig(out_dir / "coldspot_by_context_sidebyside.png", dpi=150, bbox_inches="tight")
    plt.close()

    # Figure 3: Concordance scatter (zone-level) — one color/marker per zone, legend, shared axis limits
    fig, ax = plt.subplots(figsize=(7, 6))
    colors = plt.cm.tab10(np.linspace(0, 1, max(len(zone_names), 10)))
    markers = ["o", "s", "^", "D", "v", "P", "*", "X", "h", "p"]
    for i, name in enumerate(zone_names):
        ax.scatter(
            mean_prime[i], mean_v1[i],
            s=100, c=[colors[i % len(colors)]], marker=markers[i % len(markers)],
            edgecolors="k", linewidths=0.5, zorder=3, label=name,
        )
    all_vals = [x for x in mean_prime + mean_v1 if np.isfinite(x)]
    if all_vals:
        lim_lo = min(all_vals) - 0.15
        lim_hi = max(all_vals) + 0.15
    else:
        lim_lo, lim_hi = -1.6, 0.5
    ax.plot([lim_lo, lim_hi], [lim_lo, lim_hi], "k--", alpha=0.5, label="y=x", zorder=1)
    ax.set_xlim(lim_lo, lim_hi)
    ax.set_ylim(lim_lo, lim_hi)
    ax.set_xlabel("Prime (real 5K) mean cold score")
    ax.set_ylabel("v1 (imputed 5K) mean cold score")
    sig_str = f"; p = {r_pvalue:.3f}" if np.isfinite(r_pvalue) else ""
    ax.set_title(f"Zone-level concordance (n={len(zone_names)} zones)\nPearson r = {r:.3f}{sig_str}; CCC = {ccc:.3f}")
    ax.set_aspect("equal")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=7)
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.savefig(out_dir / "coldspot_concordance_scatter.png", dpi=150, bbox_inches="tight")
    plt.close()

    return df
