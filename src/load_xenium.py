"""
Minimal Xenium loader: read cell_feature_matrix.h5 and cells.parquet into AnnData
with obsm['spatial'] and queryable var/obs. No spatialdata dependency.
"""
from pathlib import Path
import zipfile

import numpy as np
import pandas as pd
from anndata import AnnData

try:
    import h5py
except ImportError:
    h5py = None

try:
    import scanpy as sc
except ImportError:
    sc = None


def _read_matrix_h5(h5_path: Path):
    """Read 10x-style HDF5 matrix; return (counts sparse, barcodes, gene_names)."""
    if sc is not None:
        adata = sc.read_10x_h5(str(h5_path), gex_only=True)
        return adata.X, np.array(adata.obs_names), np.array(adata.var_names)
    if h5py is None:
        raise ImportError("Need scanpy or h5py to read cell_feature_matrix.h5")
    with h5py.File(h5_path, "r") as f:
        g = f["matrix"]
        shape = tuple(g["shape"][:])
        barcodes = g["barcodes"][:].astype(str)
        if barcodes.dtype.kind == "U" and barcodes.size and hasattr(barcodes[0], "decode"):
            barcodes = np.array([x.decode() if isinstance(x, bytes) else x for x in barcodes])
        features = g["features"]
        gene_names = features["name"][:].astype(str)
        if gene_names.size and hasattr(gene_names[0], "decode"):
            gene_names = np.array([x.decode() if isinstance(x, bytes) else x for x in gene_names])
        ft = features["feature_type"][:].astype(str)
        if ft.size and hasattr(ft[0], "decode"):
            ft = np.array([x.decode() if isinstance(x, bytes) else x for x in ft])
        keep = ft == "Gene Expression"
        if not np.all(keep):
            gene_names = gene_names[keep]
        data = g["data"][:]
        indices = g["indices"][:]
        indptr = g["indptr"][:]
        n_rows = shape[0]
        n_cols = shape[1]
        import scipy.sparse
        m = scipy.sparse.csc_matrix((data, indices, indptr), shape=(n_rows, n_cols))
        if not np.all(keep):
            m = m[np.where(keep)[0], :]
        m = m.T
        return m, barcodes, gene_names


def load_xenium_minimal(outs_dir: Path) -> AnnData:
    """
    Load Xenium output from a directory containing cell_feature_matrix.h5 and cells.parquet.
    Returns AnnData with X (cells x genes), var_names, obs, and obsm['spatial'] (x_centroid, y_centroid).
    """
    outs_dir = Path(outs_dir)
    h5_path = outs_dir / "cell_feature_matrix.h5"
    cells_path = outs_dir / "cells.parquet"
    if not h5_path.exists():
        raise FileNotFoundError(f"Missing {h5_path}")
    if not cells_path.exists():
        raise FileNotFoundError(f"Missing {cells_path}")

    X, barcodes, gene_names = _read_matrix_h5(h5_path)
    cells = pd.read_parquet(cells_path)

    if "cell_id" not in cells.columns:
        raise ValueError("cells.parquet must contain 'cell_id'")
    x_col = "x_centroid" if "x_centroid" in cells.columns else "x"
    y_col = "y_centroid" if "y_centroid" in cells.columns else "y"
    if x_col not in cells.columns or y_col not in cells.columns:
        raise ValueError(f"cells.parquet must contain {x_col} and {y_col} (or x_centroid, y_centroid)")

    cell_ids = cells["cell_id"].astype(str).values
    spatial = np.column_stack([cells[x_col].values.astype(float), cells[y_col].values.astype(float)])

    if len(barcodes) != len(cell_ids):
        idx = pd.Index(cell_ids).get_indexer(barcodes)
        if np.any(idx < 0):
            raise ValueError("Barcode/cell_id mismatch: some matrix barcodes not in cells.parquet")
        spatial = spatial[idx]
        cells = cells.iloc[idx].reset_index(drop=True)
    else:
        order = pd.Index(cell_ids).get_indexer(barcodes)
        if np.any(order < 0):
            raise ValueError("Barcode/cell_id mismatch")
        if not np.all(order == np.arange(len(order))):
            spatial = spatial[order]
            cells = cells.iloc[order].reset_index(drop=True)

    adata = AnnData(X, obs=cells.copy(), var=pd.DataFrame(index=gene_names))
    adata.obs_names = barcodes
    adata.obsm["spatial"] = spatial
    return adata


def load_xenium_genes_subset(outs_dir: Path, genes: list) -> AnnData:
    """
    Load Xenium output keeping only requested genes. Use when panel is large and
    only a subset of genes is needed (e.g. shared genes only) to reduce memory.
    """
    adata = load_xenium_minimal(outs_dir)
    found = [g for g in genes if g in adata.var_names]
    if not found:
        return adata[:, []].copy()
    return adata[:, found].copy()


def ensure_extracted(zip_path: Path, out_dir: Path) -> Path:
    """Extract zip to out_dir if needed; return out_dir. Only extracts h5 and parquet if possible."""
    zip_path = Path(zip_path)
    out_dir = Path(out_dir)
    if zip_path.suffix.lower() != ".zip":
        return zip_path
    out_dir.mkdir(parents=True, exist_ok=True)
    if (out_dir / "cell_feature_matrix.h5").exists() and (out_dir / "cells.parquet").exists():
        return out_dir
    with zipfile.ZipFile(zip_path, "r") as z:
        names = z.namelist()
        for name in ("cell_feature_matrix.h5", "cells.parquet"):
            for n in names:
                if n == name or n.endswith("/" + name):
                    z.extract(n, out_dir)
                    break
    h5 = next(out_dir.rglob("cell_feature_matrix.h5"), None)
    parquet = next(out_dir.rglob("cells.parquet"), None)
    if h5 and parquet:
        return h5.parent
    if (out_dir / "cell_feature_matrix.h5").exists() and (out_dir / "cells.parquet").exists():
        return out_dir
    return out_dir
