"""
Run cold-spot concordance: load saved adatas and produce side-by-side maps,
context bar chart, and zone-level concordance scatter (Pearson r, CCC).
Run after spatial_bridge_run.py (expects processed_data/adata_v1_imputed.h5ad, processed_data/adata_prime.h5ad).
"""
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

import anndata
from src.config import (
    PROCESSED_DATA_DIR,
    RESULTS_DIR,
    COLDSPOT_TUMOR_GENE,
    COLDSPOT_IMMUNE_GENE,
    COLDSPOT_NATIVE_IMMUNE_GENE,
    COLDSPOT_IMPUTED_ONLY_IMMUNE_GENE,
    COLDSPOT_SHARED_IMMUNE_GENES,
    COLDSPOT_PRIME_ONLY_IMMUNE_GENES,
    COLDSPOT_VESSEL_GENE,
    COLDSPOT_NECROSIS_GENE,
    COLDSPOT_BOUNDARY_BAND_UM,
    COLDSPOT_VESSEL_BINS_UM,
    COLDSPOT_IMPUTED_LAYER,
)
from src.coldspot_concordance import run_concordance


def main():
    adata_v1_path = PROCESSED_DATA_DIR / "adata_v1_imputed.h5ad"
    adata_prime_path = PROCESSED_DATA_DIR / "adata_prime.h5ad"
    if not adata_v1_path.exists() or not adata_prime_path.exists():
        print("Run spatial_bridge_run.py first to generate processed_data/adata_v1_imputed.h5ad and processed_data/adata_prime.h5ad")
        sys.exit(1)
    adata_v1 = anndata.read_h5ad(adata_v1_path)
    adata_prime = anndata.read_h5ad(adata_prime_path)
    run_concordance(
        adata_prime=adata_prime,
        adata_v1=adata_v1,
        tumor_gene=COLDSPOT_TUMOR_GENE,
        immune_gene=COLDSPOT_IMMUNE_GENE,
        native_immune_gene=COLDSPOT_NATIVE_IMMUNE_GENE,
        imputed_only_immune_gene=COLDSPOT_IMPUTED_ONLY_IMMUNE_GENE,
        shared_immune_genes=list(COLDSPOT_SHARED_IMMUNE_GENES) if COLDSPOT_SHARED_IMMUNE_GENES else None,
        prime_only_immune_genes=list(COLDSPOT_PRIME_ONLY_IMMUNE_GENES) if COLDSPOT_PRIME_ONLY_IMMUNE_GENES else None,
        vessel_gene=COLDSPOT_VESSEL_GENE,
        necrosis_gene=COLDSPOT_NECROSIS_GENE,
        boundary_band_um=COLDSPOT_BOUNDARY_BAND_UM,
        vessel_bins_um=COLDSPOT_VESSEL_BINS_UM,
        layer="normalized",
        imputed_layer=COLDSPOT_IMPUTED_LAYER,
        out_dir=RESULTS_DIR,
    )
    print("Cold-spot concordance written to results/: coldspot_maps_sidebyside.png, coldspot_by_context_sidebyside.png, coldspot_concordance_scatter.png, coldspot_zone_means.csv")
    print("Feature maps: coldspot_feature_tumor_boundary_sidebyside.png, coldspot_feature_necrosis_sidebyside.png, coldspot_feature_vasculature_sidebyside.png")
    print("v1 improvement: coldspot_v1_native_vs_imputed_sidebyside.png")


if __name__ == "__main__":
    main()
