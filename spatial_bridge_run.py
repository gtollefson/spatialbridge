#!/usr/bin/env python3
"""
SpatialBridge: full pipeline from Xenium v1 + Prime to aligned, imputed v1 and validation.
Run from project root. On HPCC start an interactive session first, e.g.:
  interact -n 8 -t 24:00:00 -m 200g
"""
import argparse
import sys
from pathlib import Path

# Project root and src on path
PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.config import (
    PROJECT_ROOT,
    RESULTS_DIR,
    PROCESSED_DATA_DIR,
    V1_BUNDLE,
    PRIME_BUNDLE,
    V1_EXTRACTED,
    PRIME_EXTRACTED,
    NORMALIZE_TOTAL,
    K_NEIGHBORS,
    KNN_WEIGHTS,
    N_HOLDOUT_GENES,
    ANCHOR_GENES,
    VIRTUAL_STAIN_GENES,
    CONCORDANCE_TARGET_R,
)
from src.ingest import ingest_datasets, harmonize, subset_small
from src.alignment import align_sections, plot_alignment_check
from src.inference import add_imputed_layer
from src.validation import (
    holdout_validation,
    plot_validation_scatter,
    plot_virtual_stain,
    write_results_summary,
)


def main():
    parser = argparse.ArgumentParser(description="SpatialBridge pipeline")
    parser.add_argument("--small", action="store_true", help="Use small subset for testing")
    parser.add_argument("--n-cells", type=int, default=2000, help="If --small, max cells per dataset")
    parser.add_argument("--n-shared", type=int, default=100, help="If --small, max shared genes")
    parser.add_argument("--n-prime-only", type=int, default=200, help="If --small, max Prime-only genes")
    parser.add_argument("--out-dir", type=Path, default=RESULTS_DIR, help="Results directory")
    parser.add_argument("--log", type=Path, default=PROJECT_ROOT / "AGENT_LOG.md", help="Log file")
    args = parser.parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    def log(msg: str) -> None:
        with open(args.log, "a") as f:
            f.write(msg + "\n")

    log("SpatialBridge run started.")

    # Phase 1: ingest and harmonize
    try:
        adata_v1, adata_prime = ingest_datasets(V1_BUNDLE, PRIME_BUNDLE, V1_EXTRACTED, PRIME_EXTRACTED)
    except Exception as e:
        log(f"Ingest failed: {e}")
        raise
    log("Ingest OK.")

    adata_v1, adata_prime, shared_genes, prime_only_genes = harmonize(
        adata_v1, adata_prime, normalize_total=NORMALIZE_TOTAL
    )
    log(f"Harmonize OK: {len(shared_genes)} shared, {len(prime_only_genes)} Prime-only.")

    if args.small:
        adata_v1, adata_prime, shared_genes, prime_only_genes = subset_small(
            adata_v1, adata_prime, shared_genes, prime_only_genes,
            n_cells=args.n_cells, n_shared=args.n_shared, n_prime_only=args.n_prime_only,
        )
        log(f"Subset: {adata_v1.n_obs} cells, {len(shared_genes)} shared, {len(prime_only_genes)} Prime-only.")

    # Phase 2: alignment
    adata_prime, translation = align_sections(
        adata_v1, adata_prime, anchor_genes=ANCHOR_GENES, layer="normalized"
    )
    log(f"Alignment translation: {translation}")

    plot_alignment_check(
        adata_v1, adata_prime, anchor_genes=ANCHOR_GENES, layer="normalized",
        out_path=out_dir / "alignment_check.png",
    )
    log("alignment_check.png written.")

    adata_v1 = add_imputed_layer(
        adata_v1, adata_prime, shared_genes, prime_only_genes,
        layer="normalized", k=K_NEIGHBORS, weights=KNN_WEIGHTS, imputed_layer_name="imputed_5k",
    )
    log("Imputation OK.")

    # Phase 4: validation and reporting
    pearson_r, mse, real_vals, imputed_vals = holdout_validation(
        adata_v1, adata_prime, shared_genes,
        n_holdout=min(N_HOLDOUT_GENES, max(1, len(shared_genes) - 1)),
        layer="normalized", k=K_NEIGHBORS, weights=KNN_WEIGHTS,
    )
    log(f"Validation: r={pearson_r:.4f}, MSE={mse:.4f}")

    plot_validation_scatter(real_vals, imputed_vals, pearson_r, out_path=out_dir / "validation_scatter.png")
    plot_virtual_stain(
        adata_v1, genes=("PAX8", "PTPRC") + VIRTUAL_STAIN_GENES,
        layer="normalized", imputed_layer="imputed_5k",
        shared_genes=shared_genes, prime_only_genes=prime_only_genes,
        out_path=out_dir / "virtual_stain_2x2.png",
    )
    write_results_summary(
        out_dir / "results_summary.txt",
        n_shared=len(shared_genes), n_prime_only=len(prime_only_genes),
        pearson_r=pearson_r, mse=mse, translation=translation,
        concordance_target=CONCORDANCE_TARGET_R,
    )
    processed_dir = Path(PROCESSED_DATA_DIR)
    processed_dir.mkdir(parents=True, exist_ok=True)
    adata_v1.write_h5ad(processed_dir / "adata_v1_imputed.h5ad")
    adata_prime.write_h5ad(processed_dir / "adata_prime.h5ad")
    log("results_summary.txt and plots written.")
    log("Adatas saved to processed_data/ (adata_v1_imputed.h5ad, adata_prime.h5ad).")
    log("SpatialBridge run finished.")


if __name__ == "__main__":
    main()
