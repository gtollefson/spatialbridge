"""Synthetic-data test for pipeline (no Xenium load). Requires: anndata, numpy, scipy, sklearn, matplotlib."""
import sys
from pathlib import Path

import numpy as np
from anndata import AnnData

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

# Minimal norm/log so we don't need scanpy
def _normalize_total(adata, layer, target_sum=1e4):
    X = adata.layers[layer]
    if hasattr(X, "toarray"):
        X = X.toarray()
    s = X.sum(axis=1, keepdims=True)
    s[s == 0] = 1
    adata.layers[layer] = X / s * target_sum

def _log1p(adata, layer):
    X = adata.layers[layer]
    if hasattr(X, "toarray"):
        X = X.toarray()
    adata.layers[layer] = np.log1p(X)

def main():
    np.random.seed(42)
    n = 300
    n_shared, n_prime_only = 40, 25
    shared = [f"g{i}" for i in range(n_shared)]
    prime_only = [f"p{i}" for i in range(n_prime_only)]

    v1 = AnnData(np.random.poisson(3, (n, n_shared)).astype(np.float64))
    v1.var_names = shared
    v1.obsm["spatial"] = np.random.rand(n, 2) * 200
    v1.layers["normalized"] = v1.X.copy()
    _normalize_total(v1, "normalized")
    _log1p(v1, "normalized")

    prime = AnnData(np.random.poisson(3, (n, n_shared + n_prime_only)).astype(np.float64))
    prime.var_names = shared + prime_only
    prime.obsm["spatial"] = np.random.rand(n, 2) * 200
    prime.layers["normalized"] = prime.X.copy()
    _normalize_total(prime, "normalized")
    _log1p(prime, "normalized")

    from src.alignment import align_sections, plot_alignment_check
    from src.inference import add_imputed_layer
    from src.validation import holdout_validation, plot_validation_scatter, plot_virtual_stain, write_results_summary

    prime, trans = align_sections(v1, prime, anchor_genes=("g0", "g1"), layer="normalized")
    assert prime.obsm["spatial"].shape == (n, 2)

    add_imputed_layer(v1, prime, shared, prime_only, layer="normalized", k=10, imputed_layer_name="imputed_5k")
    assert v1.layers["imputed_5k"].shape == (n, n_shared + n_prime_only)

    r, mse, real_vals, imputed_vals = holdout_validation(
        v1, prime, shared, n_holdout=min(15, len(shared) - 1), k=10
    )
    assert np.isfinite(r) or real_vals.size < 2
    assert np.isfinite(mse) or real_vals.size < 2

    out_dir = PROJECT_ROOT / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_validation_scatter(real_vals, imputed_vals, r, out_path=out_dir / "validation_scatter.png")
    plot_virtual_stain(
        v1, genes=("g0", "g1", "p0", "p1"),
        layer="normalized", imputed_layer="imputed_5k",
        shared_genes=shared, prime_only_genes=prime_only,
        out_path=out_dir / "virtual_stain_2x2.png",
    )
    write_results_summary(out_dir / "results_summary.txt", len(shared), len(prime_only), r, mse, trans)
    print("Synthetic pipeline test OK: r={:.4f}, mse={:.4f}".format(r, mse))


if __name__ == "__main__":
    main()
