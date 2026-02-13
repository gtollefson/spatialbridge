"""Paths and parameters for SpatialBridge pipeline."""
from pathlib import Path

# Paths (relative to project root)
PROJECT_ROOT = Path(__file__).resolve().parents[1]
# Demo: use original data from parent Xenium_profiler/data/ (too large for GitHub); else own data/
_ORIGINAL_DATA = PROJECT_ROOT.parent / "data"
DATA_DIR = _ORIGINAL_DATA if _ORIGINAL_DATA.exists() else PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"
PROCESSED_DATA_DIR = PROJECT_ROOT / "processed_data"  # large outputs (h5ad, etc.); keep results/ for plots/summaries

# Xenium bundles: zips or extracted dirs
V1_BUNDLE = DATA_DIR / "Xenium_V1_Human_Ovary_Cancer_FF_outs.zip"
PRIME_BUNDLE = DATA_DIR / "Xenium_Prime_Human_Ovary_Cancer_FF_outs.zip"
V1_EXTRACTED = DATA_DIR / "v1_outs"
PRIME_EXTRACTED = DATA_DIR / "prime_outs"

# Pipeline parameters
NORMALIZE_TOTAL = 1e4
K_NEIGHBORS = 15
KNN_WEIGHTS = "distance"
N_HOLDOUT_GENES = 50
ANCHOR_GENES = ("PAX8", "PTPRC")
# Optional: report "Target r > X: met/not met" in results summary (None = do not check)
CONCORDANCE_TARGET_R = 0.8
VIRTUAL_STAIN_GENES = ("IDO1", "LAG3")

# Niche enrichment (discovery)
# Define a "tumor" region using PAX8 (top quantile) and score enrichment of selected genes
# within a boundary band around the tumor boundary.
NICHE_BOUNDARY_GENE = "PAX8"
NICHE_TUMOR_QUANTILE = 0.90
NICHE_BAND_UM = 50.0
NICHE_K_NEIGHBORS = 10
NICHE_N_PERMUTATIONS = 200
NICHE_EPS = 1e-9

# Cold-spot concordance (Prime vs imputed v1)
COLDSPOT_TUMOR_GENE = "PAX8"
COLDSPOT_IMMUNE_GENE = "CD8A"  # fallback if imputed-only missing
COLDSPOT_NATIVE_IMMUNE_GENE = "CD8A"  # left panel: native v1 marker (in v1 panel)
COLDSPOT_IMPUTED_ONLY_IMMUNE_GENE = "LAG3"  # right panel: Prime-only, imputed; similar pattern to native shows imputation works
# Aggregate cold score: only genes present in adatas are used at runtime
COLDSPOT_SHARED_IMMUNE_GENES = ("PTPRC", "CD3D", "CD3E", "CD4", "CD8A", "CD8B")  # immune genes in both v1 and 5K; used for native v1
COLDSPOT_PRIME_ONLY_IMMUNE_GENES = (
    "LAG3", "IDO1", "PDCD1", "HAVCR2", "CTLA4", "TIGIT", "LAYN", "ENTPD1", "CD274", "CD80",
    "CD86", "ICOS", "TNFRSF18", "TNFRSF9", "CXCL13", "IL10", "IFNG", "GZMB", "PRF1", "NKG7",
)  # 5K-only immune/cold-spot genes (10-30); used for Prime map and imputed v1
COLDSPOT_VESSEL_GENE = "PECAM1"
COLDSPOT_NECROSIS_GENE = None  # set to gene name if in panel, else use low-RNA proxy
COLDSPOT_BOUNDARY_BAND_UM = 50.0
COLDSPOT_VESSEL_BINS_UM = (0, 50, 100, 200, 500)  # distance bins: (0,50], (50,100], ...
COLDSPOT_IMPUTED_LAYER = "imputed_5k"
