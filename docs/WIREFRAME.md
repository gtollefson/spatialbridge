# SpatialBridge Wireframe

## Objective

Upscale 10x Xenium v1 (477 genes) to Prime 5K (5000 genes) using adjacent ovarian cancer sections via spatial registration and gene-expression imputation.

## Component Overview

| Component | Technical logic | Success metric |
|-----------|-----------------|----------------|
| **I. Ingestion** | Load transcripts/cells from v1 and Prime bundles; shared gene set | Common gene count ≈ 477 |
| **II. Alignment** | Global translation from PAX8 (tumor) and PTPRC (immune) density centroids | Overlap Dice > 0.7 |
| **III. Inference** | KNeighborsRegressor: train on Prime (477 → Prime-only), predict on v1 | Pearson r (held-out) > 0.75 |
| **IV. Validation** | Leave-one-out (or leave-N-out) on shared genes | MSE < 0.1 |
| **V. Discovery** | Spatial “virtual stains” (e.g. LAG3, IDO1) vs v1 tumor boundaries | Niche enrichment |

## Data Flow

```
[v1 outs] + [Prime outs]
       → Ingest (shared genes, normalize total 1e4, log1p)
       → Align (PAX8/PTPRC centroid translation on Prime)
       → Infer (KNN: Prime 477→Prime-only; impute v1 → imputed_5k)
       → Validate (holdout correlation, virtual stain 2×2)
       → results_summary.txt, alignment_check.png, validation_scatter.png
```

## Key Parameters

- Normalization: total count 1e4, then log1p → `adata.layers['normalized']`
- Alignment: 2D density centroids for PAX8 and PTPRC; translation (dx, dy) applied to Prime `obsm['spatial']`
- Inference: `KNeighborsRegressor(k=15, weights='distance')`, train on Prime, predict for v1 cells
- Validation: hold out 50 shared genes, impute, report Pearson r and MSE

## Deliverables

- `spatial_bridge_run.py`: single entry point for the pipeline
- `results/alignment_check.png`: PAX8/PTPRC overlap
- `results/validation_scatter.png`: real vs imputed correlation
- `results_summary.txt`: executive summary
- Optional: `AGENT_LOG.md` for run notes and errors
