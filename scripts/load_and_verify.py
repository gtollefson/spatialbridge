#!/usr/bin/env python3
"""
Load Xenium v1 and Prime via foundational API and confirm spatial data is queryable.
Demonstrates efficient access helpers. Run from project root.
"""
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from src.spatial_data import (
    get_dataset_path,
    load_dataset,
    get_expression,
    get_spatial,
    cells_in_bbox,
    shared_genes,
    prime_only_genes,
    expression_weighted_centroid,
)


def main():
    print("Resolving paths (extract if needed)...")
    v1_path = get_dataset_path("v1")
    prime_path = get_dataset_path("prime")
    print(f"  v1:   {v1_path}")
    print(f"  Prime: {prime_path}")

    print("Loading v1 (full)...")
    adata_v1 = load_dataset("v1")
    print(f"  shape: {adata_v1.n_obs} x {adata_v1.n_vars}")

    print("Loading Prime (full)...")
    adata_prime = load_dataset("prime")
    print(f"  shape: {adata_prime.n_obs} x {adata_prime.n_vars}")

    print("\n--- Query checks ---")
    xy = get_spatial(adata_v1)
    print(f"  get_spatial(v1): shape {xy.shape}, x [{xy[:, 0].min():.1f}, {xy[:, 0].max():.1f}], y [{xy[:, 1].min():.1f}, {xy[:, 1].max():.1f}]")

    gene = adata_v1.var_names[0]
    v = get_expression(adata_v1, gene)
    print(f"  get_expression(v1, '{gene}'): shape {v.shape}, sum={v.sum():.0f}, n_nonzero={(v > 0).sum()}")

    shared = shared_genes(adata_v1, adata_prime)
    prime_only = prime_only_genes(adata_v1, adata_prime)
    print(f"  shared_genes: {len(shared)}")
    print(f"  prime_only_genes: {len(prime_only)}")

    if "PAX8" in adata_v1.var_names:
        c = expression_weighted_centroid(adata_v1, "PAX8")
        print(f"  expression_weighted_centroid(v1, 'PAX8'): ({c[0]:.1f}, {c[1]:.1f})")

    x_min, x_max = float(xy[:, 0].mean() - 200), float(xy[:, 0].mean() + 200)
    y_min, y_max = float(xy[:, 1].mean() - 200), float(xy[:, 1].mean() + 200)
    mask = cells_in_bbox(adata_v1, x_min, x_max, y_min, y_max)
    print(f"  cells_in_bbox(v1, 400x400 region): {mask.sum()} cells")

    print("\n--- Load subset (efficient) ---")
    adata_v1_sub = load_dataset("v1", genes=shared[:10], n_cells=1000)
    print(f"  load_dataset('v1', genes=shared[:10], n_cells=1000): {adata_v1_sub.n_obs} x {adata_v1_sub.n_vars}")

    print("\nDone.")


if __name__ == "__main__":
    main()
