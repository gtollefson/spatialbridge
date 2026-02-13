"""Load Xenium v1 and Prime data, find shared genes, normalize."""
from pathlib import Path
from typing import Optional
import numpy as np
from anndata import AnnData

from src.load_xenium import ensure_extracted, load_xenium_minimal

try:
    import scanpy as sc
except ImportError:
    sc = None

try:
    from spatialdata_io import xenium
    import spatialdata as sd
except ImportError:
    xenium = None
    sd = None


def _normalize_total_layer(adata: AnnData, layer: str, target_sum: float = 1e4) -> None:
    X = adata.layers[layer]
    if hasattr(X, "toarray"):
        X = X.toarray()
    s = np.array(X.sum(axis=1)).ravel()
    s[s == 0] = 1
    adata.layers[layer] = X / s[:, np.newaxis] * target_sum


def _log1p_layer(adata: AnnData, layer: str) -> None:
    X = adata.layers[layer]
    if hasattr(X, "toarray"):
        X = X.toarray()
    adata.layers[layer] = np.log1p(np.asarray(X, dtype=float))


def _ensure_extracted(zip_path: Path, out_dir: Path) -> Path:
    """Extract zip if needed; return path to directory with cell_feature_matrix.h5 and cells.parquet."""
    if zip_path.suffix.lower() != ".zip":
        return zip_path
    return ensure_extracted(zip_path, out_dir)


def load_xenium(path: Path) -> AnnData:
    """Load Xenium output (extracted directory) into AnnData with obsm['spatial']."""
    path = Path(path)
    if path.suffix.lower() == ".zip":
        raise ValueError("Pass an extracted directory or use ingest_datasets() which extracts zips.")
    if xenium is not None and sd is not None:
        sdata = xenium(str(path))
        return sdata.tables["table"].copy()
    return load_xenium_minimal(path)


def ingest_datasets(
    v1_bundle: Path,
    prime_bundle: Path,
    v1_extracted: Path,
    prime_extracted: Path,
) -> tuple[AnnData, AnnData]:
    """Extract bundles if zips, load v1 and Prime, return (adata_v1, adata_prime)."""
    v1_path = _ensure_extracted(v1_bundle, v1_extracted)
    prime_path = _ensure_extracted(prime_bundle, prime_extracted)
    adata_v1 = load_xenium(v1_path)
    adata_prime = load_xenium(prime_path)
    return adata_v1, adata_prime


def harmonize(
    adata_v1: AnnData,
    adata_prime: AnnData,
    normalize_total: float = 1e4,
) -> tuple[AnnData, AnnData, list[str], list[str]]:
    """Subset v1 to shared genes; prime to shared+prime_only. Normalize (total 1e4, log1p) -> layer 'normalized'."""
    shared = list(set(adata_v1.var_names) & set(adata_prime.var_names))
    shared.sort()
    prime_only = [g for g in adata_prime.var_names if g not in set(adata_v1.var_names)]
    prime_only.sort()
    v1 = adata_v1[:, shared].copy()
    prime = adata_prime[:, shared + prime_only].copy()

    for adata in (v1, prime):
        adata.layers["normalized"] = adata.X.copy()
        if sc is not None:
            sc.pp.normalize_total(adata, target_sum=normalize_total, layer="normalized")
            sc.pp.log1p(adata, layer="normalized")
        else:
            _normalize_total_layer(adata, "normalized", normalize_total)
            _log1p_layer(adata, "normalized")

    return v1, prime, shared, prime_only


def subset_small(
    adata_v1: AnnData,
    adata_prime: AnnData,
    shared_genes: list[str],
    prime_only_genes: list[str],
    n_cells: int = 2000,
    n_shared: Optional[int] = 100,
    n_prime_only: Optional[int] = 200,
) -> tuple[AnnData, AnnData, list[str], list[str]]:
    """Subset for fast testing: fewer cells and optionally fewer genes. Returns (v1, prime, shared, prime_only)."""
    v1 = adata_v1.copy()
    prime = adata_prime.copy()
    shared = list(shared_genes)
    prime_only = list(prime_only_genes)
    if n_cells < v1.n_obs:
        idx = np.random.default_rng(0).choice(v1.n_obs, size=n_cells, replace=False)
        v1 = v1[idx].copy()
    if n_cells < prime.n_obs:
        idx = np.random.default_rng(0).choice(prime.n_obs, size=n_cells, replace=False)
        prime = prime[idx].copy()
    if n_shared is not None and len(shared) > n_shared:
        # Keep cold-spot genes (PAX8, PTPRC, PECAM1, etc.) for demo compatibility
        from src.config import ANCHOR_GENES, COLDSPOT_TUMOR_GENE, COLDSPOT_VESSEL_GENE
        required = set(ANCHOR_GENES) | {COLDSPOT_TUMOR_GENE, COLDSPOT_VESSEL_GENE}
        keep = [g for g in shared if g in required]
        rest = [g for g in shared if g not in keep]
        shared = keep + rest[: max(0, n_shared - len(keep))]
    if n_prime_only is not None and len(prime_only) > n_prime_only:
        from src.config import (
            COLDSPOT_IMPUTED_ONLY_IMMUNE_GENE,
            COLDSPOT_TUMOR_GENE,
            COLDSPOT_VESSEL_GENE,
            VIRTUAL_STAIN_GENES,
        )
        keep = [
            g for g in prime_only
            if g in {COLDSPOT_IMPUTED_ONLY_IMMUNE_GENE, COLDSPOT_TUMOR_GENE, COLDSPOT_VESSEL_GENE}
            or g in VIRTUAL_STAIN_GENES
        ]
        rest = [g for g in prime_only if g not in keep]
        prime_only = keep + rest[: max(0, n_prime_only - len(keep))]
    v1 = v1[:, shared].copy()
    prime = prime[:, shared + prime_only].copy()
    return v1, prime, shared, prime_only
