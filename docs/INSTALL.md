# Installation and Build

## Requirements

- Python 3.10+
- See `requirements.txt` for dependencies (anndata, scanpy, squidpy, spatialdata, spatialdata-io, scikit-learn, matplotlib, etc.)

## Quick install

```bash
git clone <repository-url>
cd clean_repo_xenium
pip install -r requirements.txt
```

## Virtual environment (recommended)

```bash
python -m venv .venv
source .venv/bin/activate   # Linux/macOS
# or: .venv\Scripts\activate   # Windows
pip install -r requirements.txt
```

## HPCC / cluster

If running on a high-performance cluster (e.g. HPCC), start an interactive session with sufficient resources before running the pipeline:

```bash
interact -n 8 -t 24:00:00 -m 200g
```

Then install and run from the project root.

## Optional: minimal load-and-verify

The `scripts/load_and_verify.py` script can run with fewer dependencies (anndata, pandas, scanpy or h5py+scipy) for a quick data check without the full pipeline.
