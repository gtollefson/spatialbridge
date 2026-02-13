# SpatialBridge: Complete Project Workflow and Wireframe

This document is the **resurrection guide** for the SpatialBridge project. Use it to understand the full workflow, where everything lives, and how to make edits or extend the pipeline.

---

## 1. Purpose and scope

**SpatialBridge** does two things:

1. **Upsample Xenium v1 (477 genes) to a Prime 5K–like panel** using an adjacent Xenium Prime 5K section: spatial registration (PAX8/PTPRC anchors) + KNN-based imputation of Prime-only genes into v1 space. Outputs: aligned Prime, v1 with `imputed_5k` layer, validation metrics, virtual stains.
2. **Cold-spot concordance**: Compare cold-spot geography between Prime (real 5K) and v1 (imputed 5K), and between v1 native (no imputation) and v1 imputed. Uses **aggregate cold scores** from two gene lists (shared vs 5K-only), spatial zones (tumor boundary, vessel distance, necrosis), and zone-level concordance (Pearson r, CCC, p-value).

**Data**: Adjacent FFPE sections—v1 and Prime 5K—from the same block (e.g. ovarian cancer). Inputs are Xenium output bundles (zips or extracted dirs).

---

## 2. End-to-end workflow (two phases)

```
Phase 1: SpatialBridge pipeline (spatial_bridge_run.py)
  data/ (v1 + Prime bundles)
    → Ingest (load, extract if needed)
    → Harmonize (shared genes, normalize 1e4 + log1p → layer 'normalized')
    → Align (PAX8/PTPRC centroids, translate Prime to v1)
    → Impute (KNN: shared → Prime-only on Prime; predict for v1 → layer 'imputed_5k')
    → Validate (holdout shared genes, Pearson r, MSE, virtual stain)
    → Write: results/*.png, results/results_summary.txt
    → Write: processed_data/adata_v1_imputed.h5ad, processed_data/adata_prime.h5ad

Phase 2: Cold-spot concordance (scripts/coldspot_concordance_run.py)
  processed_data/ (adata_v1_imputed.h5ad, adata_prime.h5ad)
    → Load adatas
    → Compute aggregate cold scores (List A = shared for native v1; List B = 5K-only for Prime and imputed v1)
    → Define zones (tumor boundary band, vessel distance bins, necrotic/low-RNA)
    → Zone-level means, concordance (r, CCC, p-value)
    → Write: results/coldspot_*.png, coldspot_zone_means.csv, coldspot_concordance_stats.csv, coldspot_interpretation.txt
```

**Order**: Run Phase 1 first (creates processed_data). Then run Phase 2 (reads processed_data, writes only to results/).

---

## 3. Data flow (detailed)

| Step | Input | Action | Output |
|------|--------|--------|--------|
| Ingest | `data/*.zip` or extracted dirs | Load v1 and Prime via spatialdata-io; get shared and Prime-only gene lists | `adata_v1`, `adata_prime` (raw then normalized in next step) |
| Harmonize | adatas | Normalize total 1e4, log1p → `layers['normalized']`; subset to shared + prime_only | `shared_genes`, `prime_only_genes`; adatas with normalized layer |
| Align | adatas, ANCHOR_GENES | PAX8/PTPRC expression-weighted centroids; translation (dx, dy); apply to Prime `obsm['spatial']` | `adata_prime` with updated spatial coords |
| Impute | adatas, shared, prime_only | One KNN per Prime-only gene (shared → gene on Prime); predict for v1; stack shared + imputed | `adata_v1` with `layers['imputed_5k']`, vars = shared + prime_only |
| Validate | adatas, holdout genes | Train KNN for each held-out shared gene; predict in v1; Pearson r, MSE; plot virtual stain | results_summary.txt, validation_scatter.png, virtual_stain_2x2.png |
| Save adatas | adatas | Write to processed_data | adata_v1_imputed.h5ad, adata_prime.h5ad |
| Cold-spot | processed_data adatas | Aggregate scores (List A native, List B Prime/imputed); zones; figures and CSVs | results/coldspot_* |

---

## 4. File map (where to edit what)

| Path | Role | Edit when |
|------|------|-----------|
| **Entry points** | | |
| `spatial_bridge_run.py` | Phase 1: ingest → align → impute → validate → save adatas | Change pipeline order, add steps, change CLI |
| `scripts/coldspot_concordance_run.py` | Phase 2: load adatas, call `run_concordance`, print paths | Change cold-spot inputs or output paths |
| **Configuration** | | |
| `src/config.py` | All paths and parameters (pipelines, cold-spot, discovery) | Change data paths, K, gene lists, zone params |
| **Pipeline (Phase 1)** | | |
| `src/ingest.py` | Load v1/Prime, extract zips; `harmonize()`, `subset_small()` | Change normalization, subset logic |
| `src/alignment.py` | PAX8/PTPRC centroids, translation, `plot_alignment_check` | Change anchor genes or alignment model |
| `src/inference.py` | `add_imputed_layer()`: one KNN space, stream Prime-only genes | Change KNN (k, weights), memory strategy |
| `src/validation.py` | Holdout validation, `plot_validation_scatter`, `plot_virtual_stain`, `write_results_summary` | Change metrics, plot layout, summary format |
| **Cold-spot (Phase 2)** | | |
| `src/coldspot_concordance.py` | `cold_spot_score()`, `cold_spot_score_aggregate()`, zones, `run_concordance()`, all cold-spot figures and CSVs | Change score formula, zone definitions, new figures |
| **Data and discovery** | | |
| `src/spatial_data.py` | High-level load, get_expression, get_spatial, shared_genes, prime_only_genes, centroids | Change data API or paths |
| `src/load_xenium.py` | Low-level H5/parquet load from Xenium bundles | Change raw read logic |
| `src/discovery.py` | Niche enrichment (boundary band, log2 enrichment, permutation p); helpers: _get_vector, _tumor_mask_from_quantile, _boundary_mask_knn, _distance_to_boundary | Change niche definition or reuse helpers |
| **Tests and docs** | | |
| `tests/test_pipeline_synthetic.py` | Synthetic adatas, pipeline smoke test | Add tests for new steps |
| `scripts/load_and_verify.py` | Quick check: load data, shared/prime_only counts | Verify data without full run |
| `README.md` | User-facing install, run, outputs | Update when commands or outputs change |
| `docs/WIREFRAME.md` | Original high-level pipeline design | Legacy; PROJECT_WORKFLOW is the single wireframe |
| `METHODS_NOTES.md` | Methods text for manuscript | Update when methods change |
| `presentation/slide_outline.txt` | Slide headings, contents for talks | Update when results or approach change |

---

## 5. Configuration reference (`src/config.py`)

All tunable parameters live here. Key groups:

**Paths**
- `DATA_DIR`, `RESULTS_DIR`, `PROCESSED_DATA_DIR` — where data lives, where results/plots go, where large h5ad files go.
- `V1_BUNDLE`, `PRIME_BUNDLE`, `V1_EXTRACTED`, `PRIME_EXTRACTED` — Xenium bundles and extracted dirs.

**Pipeline (Phase 1)**
- `NORMALIZE_TOTAL=1e4`, `K_NEIGHBORS=15`, `KNN_WEIGHTS='distance'`, `N_HOLDOUT_GENES=50`, `ANCHOR_GENES=('PAX8','PTPRC')`, `CONCORDANCE_TARGET_R=0.8`, `VIRTUAL_STAIN_GENES=('IDO1','LAG3')`.

**Cold-spot (Phase 2)**
- `COLDSPOT_TUMOR_GENE='PAX8'` — tumor mask for “tumor-only” coloring and boundary.
- `COLDSPOT_SHARED_IMMUNE_GENES` — tuple of immune genes in both v1 and 5K; used for **native v1** aggregate score only.
- `COLDSPOT_PRIME_ONLY_IMMUNE_GENES` — tuple of 5K-only immune genes (10–30); used for **Prime** and **imputed v1** aggregate scores.
- `COLDSPOT_NATIVE_IMMUNE_GENE`, `COLDSPOT_IMPUTED_ONLY_IMMUNE_GENE`, `COLDSPOT_IMMUNE_GENE` — single-gene fallbacks when aggregate lists are empty or no genes pass filter.
- `COLDSPOT_VESSEL_GENE='PECAM1'`, `COLDSPOT_BOUNDARY_BAND_UM=50`, `COLDSPOT_VESSEL_BINS_UM=(0,50,100,200,500)`, `COLDSPOT_IMPUTED_LAYER='imputed_5k'`.

**Discovery (optional)**
- `NICHE_*` — boundary band, tumor quantile, permutations for niche enrichment (used in discovery, not in main run script by default).

---

## 6. Run commands (resurrection checklist)

From project root (`Xenium_profiler/`):

```bash
# Optional: verify data load only (no pipeline)
python scripts/load_and_verify.py

# 1. Start interactive session on HPCC (recommended for full run)
interact -n 8 -t 24:00:00 -m 200g

# 2. Phase 1: full pipeline (writes results/ and processed_data/)
python spatial_bridge_run.py

# Or quick test (subset cells and genes)
python spatial_bridge_run.py --small

# 3. Phase 2: cold-spot concordance (requires processed_data/ from step 2)
python scripts/coldspot_concordance_run.py
```

**Requirements**: Python 3.10+; install deps with `pip install -r requirements.txt`. Data: place v1 and Prime Xenium bundles in `data/` (see README for names).

---

## 7. Outputs (where everything lands)

| Output | Location | Produced by |
|--------|----------|-------------|
| alignment_check.png | results/ | Phase 1 |
| validation_scatter.png | results/ | Phase 1 |
| virtual_stain_2x2.png | results/ | Phase 1 |
| results_summary.txt | results/ | Phase 1 |
| adata_v1_imputed.h5ad, adata_prime.h5ad | processed_data/ | Phase 1 |
| coldspot_maps_sidebyside.png | results/ | Phase 2 |
| coldspot_v1_native_vs_imputed_sidebyside.png | results/ | Phase 2 |
| coldspot_by_context_sidebyside.png | results/ | Phase 2 |
| coldspot_concordance_scatter.png | results/ | Phase 2 |
| coldspot_feature_tumor_boundary_sidebyside.png | results/ | Phase 2 |
| coldspot_feature_necrosis_sidebyside.png | results/ | Phase 2 |
| coldspot_feature_vasculature_sidebyside.png | results/ | Phase 2 |
| coldspot_zone_means.csv | results/ | Phase 2 |
| coldspot_concordance_stats.csv | results/ | Phase 2 |
| coldspot_interpretation.txt | results/ | Phase 2 |
| AGENT_LOG.md | project root | Phase 1 (append log) |

---

## 8. Cold-spot logic (current approach)

- **Score**: Immune-only. High score = low immune = cold. Formula: per cell, `score = -mean(z)` over a list of immune genes (z = z-score of expression).
- **Two lists**:
  - **List A (shared)**: Genes in both v1 and 5K. Used only for **native v1** score (normalized layer; no imputation). Fewer genes → noisier, more extreme scores.
  - **List B (5K-only)**: Genes only in 5K panel. Used for **Prime** map and **imputed v1** map. More genes → smoother, more neutral-looking map.
- **Tumor mask**: PAX8-high cells (e.g. top 90th percentile). Only tumor cells are colored by cold score; non-tumor shown as light grey.
- **Zones** (for concordance): (1) Tumor boundary band (cells within 50 µm of PAX8-boundary); (2) Vessel distance bins (distance to PECAM1-high cells, bins 0–50, 50–100, 100–200, 200–500 µm); (3) Necrotic / low-RNA (bottom 2% total expression). Per zone: mean cold score in Prime, mean in v1 imputed → scatter and bar chart; Pearson r, CCC, p-value.
- **Fallback**: If aggregate lists are empty or no genes pass the presence filter, single-gene cold score is used (config: native and imputed-only gene names).

---

## 9. Extension points (where to add or change behavior)

- **New pipeline step**: Add a step in `spatial_bridge_run.py` (after impute or after save); implement in a new or existing module under `src/`.
- **New cold-spot figure or metric**: In `src/coldspot_concordance.py`, inside `run_concordance()`, add computation and a new `plt.savefig(out_dir / "coldspot_*.png")` (and optionally extend `coldspot_interpretation.txt`).
- **New zone type**: In `coldspot_concordance.py`, add a new mask (e.g. from a new gene or distance), append to `zone_names`, `mean_prime`, `mean_v1`, `n_cells_*`; figures and concordance will include it automatically.
- **Different aggregate formula**: Change `cold_spot_score_aggregate()` (e.g. median instead of mean, or weighted mean); keep the same signature so `run_concordance()` does not need changes.
- **New config parameter**: Add in `src/config.py` and pass through from `scripts/coldspot_concordance_run.py` or `spatial_bridge_run.py` as needed.
- **Data from a different assay or path**: Adjust `config.py` paths and, if needed, `src/load_xenium.py` or ingest (e.g. different reader) while keeping the same AnnData layout (obsm['spatial'], layers['normalized'], etc.).

---

## 10. Dependencies and environment

- **Python**: 3.10+.
- **Key packages**: anndata, scanpy, pandas, numpy, scipy, scikit-learn, matplotlib; spatialdata, spatialdata-io for Xenium load. See `requirements.txt`.
- **Optional**: `scripts/load_and_verify.py` can run with minimal deps (anndata, pandas, h5py, scipy) for a quick data check.
- **Large runs**: Full pipeline and cold-spot script can be memory-heavy; use an interactive session with sufficient RAM (e.g. 200 GB) and run from project root so `src` and `data` resolve correctly.

---

## 11. References to other docs

- **README.md** — Install, data layout, run commands, project layout, limitations.
- **docs/WIREFRAME.md** — Original pipeline wireframe (ingest → align → infer → validate); superseded by this workflow for full scope.
- **METHODS_NOTES.md** — Text for manuscript methods (data, preprocessing, registration, imputation, validation).
- **presentation/slide_outline.txt** — Slide headings, subtitles, and content outline for talks (biology, math, pipeline, cold-spot, results).

---

*Last updated to reflect: two-phase workflow (pipeline + cold-spot), processed_data for h5ad, aggregate cold scores (List A / List B), tumor-only coloring, zone-level concordance, and all current output files.*
