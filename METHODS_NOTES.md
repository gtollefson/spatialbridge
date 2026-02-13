# Methods notes for manuscript

Use these notes when writing the SpatialBridge methods section.

## Data

- Adjacent FFPE sections: Xenium v1 (Human Multi-Tissue and Cancer + 100 custom genes, ~477 genes) and Xenium Prime 5K (Pan Tissue and Pathways, ~5000 genes). Ovarian cancer (Ovarian Papillary Serous Carcinoma, II-A).
- Xenium Onboard Analysis / Xenium Ranger (e.g. v4.0) output: cell feature matrix, cells (coordinates and metadata), transcripts. Loaded via spatialdata-io Xenium reader.

## Preprocessing and harmonization

- Shared genes between v1 and Prime were identified. Total-count normalization (target sum 10,000 per cell) followed by log1p transformation was applied to both datasets; normalized values stored in a dedicated layer for downstream steps.

## Spatial registration

- To align the Prime section to the v1 section, a shared-anchor strategy was used. For each of two anchor genes—PAX8 (tumor) and PTPRC (immune)—a 2D expression-weighted centroid was computed from cell-level normalized expression and spatial coordinates. A single translation vector (dx, dy) was derived as the average of the per-gene centroid differences (v1 minus Prime). This translation was applied to all Prime cell coordinates. No rotation or scaling was applied.

## Gene expression imputation (inference)

- For each Prime-only gene, a k-nearest-neighbors regressor (k = 15, distance-weighted) was trained on the Prime dataset using the shared-gene normalized expression as predictors and the Prime-only gene’s normalized expression as the response. The fitted model was used to predict that gene’s expression for every v1 cell from its shared-gene expression. The resulting v1 dataset thus had a layer with dimensions (cells × [shared + imputed Prime-only genes]).

## Validation

- Holdout validation: 50 shared genes were held out at random. For each held-out gene, a KNN regressor was trained on Prime (remaining shared genes → held-out gene) and used to predict that gene in v1. Pearson correlation and mean squared error (MSE) between observed and imputed expression across all v1 cells and held-out genes were reported.
- Virtual stains: spatial plots of real PAX8 and PTPRC and imputed IDO1 and LAG3 were generated to illustrate imputed Prime-only expression in v1 space.

## Software and parameters

- Python 3.10+; scanpy, squidpy, spatialdata, spatialdata-io, scikit-learn. Key parameters: normalization target 10,000; k = 15; weights = distance; 50 holdout genes; anchor genes PAX8, PTPRC.
