# Xenium Spatial Transcriptomics Analysis Pipeline

A Jupyter-notebook-based analysis pipeline for 10x Xenium spatial
transcriptomics data with optional Phenocycler (CODEX) integration.
Optimized for Linux workstations and the University of Florida HiPerGator
HPC cluster, with scanpy / squidpy / scvi-tools on the Python side and
Seurat on the R side.

## Table of Contents
- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Workflow](#workflow)
- [Exporting Clusters / Phenotypes to Xenium Explorer](#exporting-clusters--phenotypes-to-xenium-explorer)
- [HiPerGator Usage](#hipergator-usage)
- [Directory Structure](#directory-structure)
- [Notebooks](#notebooks)
- [Utility Package (`utils`)](#utility-package-utils)
- [Continuous Integration and Testing](#continuous-integration-and-testing)
- [Configuration](#configuration)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)

## Features

### Analysis
- **Preprocessing**: QC, filtering, normalization from Xenium h5/h5ad
- **Phenotyping**: marker-gene scoring + scVI-based clustering
- **Spatial analysis**: neighborhood enrichment, co-occurrence, Moran's I,
  Ripley's L, ligand-receptor
- **Group / tissue comparisons**: composition tests, DEG with fixed
  log-fold-change on Xenium counts, volcano plots
- **Phenocycler integration**: coordinate alignment, nearest-neighbor
  cross-modal cell mapping
- **Xenium Explorer export**: push any `.obs` categorical (clusters,
  phenotypes, scVI labels) back to Xenium Explorer 3.0+ as Cell Groups
  CSV or `analysis.zarr.zip`

### Tooling
- **Python 3.12** (pinned in `environment.yml`), with a declared minimum of
  3.10 in `pyproject.toml` so contributor machines can still run the lint +
  test suite
- **R 4.3 / Seurat** via `rpy2` for statistical follow-up
- **SLURM** scripts for every pipeline stage
- **CI/CD**: lightweight GitHub Actions (ruff + pytest + shellcheck) plus a
  HiPerGator-side SLURM runner for the full heavy suite
- **Pre-commit hooks** that match CI exactly

## Installation

### 1. Clone

```bash
git clone https://github.com/smith6jt-cop/Xenium_Analysis.git
cd Xenium_Analysis
```

### 2. Create the conda environment

```bash
conda env create -f environment.yml
conda activate xenium_analysis
```

This creates the full runtime environment:
- Python 3.12, numpy, pandas, scipy, scikit-learn
- Single-cell: scanpy, squidpy, scvi-tools, anndata
- Spatial: shapely, spatialdata (+ `spatialdata-io`)
- R: `r-base`, Seurat, ggplot2, tidyverse + `rpy2`
- Visualization: matplotlib, seaborn, plotly, napari
- I/O: h5py, zarr, tifffile, openpyxl

### 3. Verify

```bash
python -c "import scanpy as sc, squidpy as sq, scvi; print('ok')"
```

### 4. (Optional) install developer hooks

```bash
pip install pre-commit && pre-commit install
```

This wires up ruff, shellcheck, and line-ending hooks so commits match what
CI enforces.

## Quick Start

### 1. Drop data in `data/raw/`

See [DATA_README.md](DATA_README.md) for supported formats (`.h5ad`, 10x
`.h5`, or a Xenium output directory).

### 2. Run locally

```bash
conda activate xenium_analysis
jupyter lab
# Open notebooks in order: 01_preprocessing.ipynb -> ... -> 06_...
```

### 3. Run on HiPerGator

```bash
# One stage
sbatch scripts/slurm/01_run_preprocessing.sh

# Full pipeline
sbatch scripts/slurm/run_full_pipeline.sh
```

## Workflow

```
Raw h5 / Xenium output
          │
          ▼
01_preprocessing.ipynb            ──▶ *_preprocessed.h5ad
          │
          ▼
02_phenotyping.ipynb              ──▶ *_annotated.h5ad
          │                                │
          │                                └──▶ export_for_xenium_explorer()
          │                                      ▶ analysis.zarr.zip + CSV
          ▼
03_spatial_analysis.ipynb         ──▶ *_spatial_analysis.h5ad
          │
          ▼
04_group_comparisons.ipynb        ──▶ DEG tables, volcanos
          │
          ▼
05_tissue_comparisons.ipynb       ──▶ integrated_tissues.h5ad
          │
          ▼
06_xenium_phenocycler_integration ──▶ *_integrated.h5ad
```

## Exporting Clusters / Phenotypes to Xenium Explorer

After phenotyping, any `.obs` categorical (Leiden clusters, scVI labels,
curated phenotypes, manual annotations, ...) can be pushed back to
**Xenium Explorer 3.0+** so cells are colored by that grouping in the
viewer.

The `utils/xenium_explorer_export` module supports both Xenium Explorer
import paths and only depends on `numpy`, `pandas`, `anndata`, and `zarr`
(no scanpy), so you can call it from any environment — including plain
CI.

### When to use which format

| Format | Loaded via | Matched by | Best for |
|---|---|---|---|
| **Cell Groups CSV** | Explorer → Cell → *Add cell categorization* | `cell_id` | After any filtering / reordering — always safe |
| **`analysis.zarr.zip`** | Drop next to `experiment.xenium` → auto-picked up in the *Clusters* dropdown | Row position | When the AnnData row order matches the original `cells.parquet` |

Most notebook callers should just use the top-level
`export_for_xenium_explorer()` — it writes both and also dumps a JSON
palette so matplotlib figures match the Explorer view.

### From a notebook

```python
from utils import export_for_xenium_explorer

artifacts = export_for_xenium_explorer(
    adata,
    group_keys=["leiden_0.5", "celltype", "phenotype"],
    output_dir="data/xenium_explorer",
    sample_name="THYHDL065",
    cell_id_key="cell_id",          # None -> auto-detects from adata.obs or uses the index
    palette="tab20",                # any matplotlib categorical colormap
    rename={"leiden_0.5": "Clusters"},  # optional: prettier names in Explorer
)
# artifacts = {"csv": ..., "zarr": ..., "palette": ...}
```

This writes:

```
data/xenium_explorer/
├── THYHDL065_cell_groups.csv      # Explorer → Cell → Add cell categorization
├── THYHDL065/
│   └── analysis.zarr.zip          # drop next to experiment.xenium
└── THYHDL065_color_palette.json   # {grouping: {category: '#RRGGBB'}}
```

### From the CLI

```bash
python scripts/export_xenium_explorer_groups.py \
    --input   data/processed/THYHDL065_annotated.h5ad \
    --output  data/xenium_explorer \
    --sample  THYHDL065 \
    --groups  leiden_0.5 celltype phenotype \
    --cell-id-key cell_id
```

Flags:

| Flag | Purpose |
|---|---|
| `--input` / `-i` | Path to the `.h5ad` |
| `--output` / `-o` | Output directory |
| `--sample` / `-s` | Prefix for all emitted files |
| `--groups` / `-g` | One or more `.obs` column names |
| `--cell-id-key` | Column with the Xenium `cell_id` (defaults to common candidates) |
| `--palette` | Matplotlib colormap name (default `tab20`) |
| `--no-csv` / `--no-zarr` | Skip one of the two artifacts |

### Loading in Xenium Explorer

**CSV** — open the Xenium sample, go to the **Cell** panel →
*Add cell categorization*, pick the CSV. Every non-`cell_id` column becomes
a selectable layer in the Cell ▸ Groups dropdown.

**analysis.zarr.zip** — copy or symlink it to the same directory as the
sample's `experiment.xenium` before opening. Xenium Explorer picks it up
automatically; your groupings then appear under the built-in **Clusters**
dropdown with the colors stored in the zarr attrs.

### Public API

Importable from `utils.xenium_explorer_export` or the top-level `utils`
package:

| Function | Purpose |
|---|---|
| `export_for_xenium_explorer(adata, group_keys, output_dir, sample_name, ...)` | High-level: writes CSV + `analysis.zarr.zip` + palette JSON |
| `export_groups_to_csv(adata, group_keys, output_path, cell_id_key=None, rename=None)` | Cell Groups CSV only |
| `export_groups_to_zarr(adata, group_keys, output_path, cell_id_key=None, rename=None, colors=None, palette='tab20')` | `analysis.zarr.zip` only |
| `generate_color_palette(categories, palette='tab20')` | `{category: '#RRGGBB'}` dict, with a matplotlib-free fallback |

### Troubleshooting

- **Xenium Explorer does not see my categorization** — the CSV's first
  column must be exactly `cell_id` and values must match the sample's
  `cells.parquet`. Run `adata.obs['cell_id'].head()` and compare.
- **Colors don't match my Python figures** — use the emitted
  `*_color_palette.json` as the `palette` argument in your plotting code
  (e.g. `sc.pl.umap(adata, palette=palette_from_json)`).
- **`analysis.zarr.zip` shows the wrong cells** — this format is
  position-indexed and breaks if the AnnData row order differs from the
  original `cells.parquet`. Use the CSV instead, which matches by
  `cell_id`.

## HiPerGator Usage

### Initial setup

```bash
module load conda
conda env create -f environment.yml
```

Edit the SLURM headers in `scripts/slurm/*.sh` with your values:

```bash
#SBATCH --mail-user=your_email@ufl.edu
#SBATCH --qos=your_qos
#SBATCH --account=your_account
```

### Submitting jobs

```bash
mkdir -p logs

# One stage
sbatch scripts/slurm/01_run_preprocessing.sh
sbatch scripts/slurm/02_run_phenotyping.sh
sbatch scripts/slurm/03_run_spatial_analysis.sh

# Full pipeline
sbatch scripts/slurm/run_full_pipeline.sh

# Full test suite (see "CI" below)
sbatch scripts/slurm/run_tests.sh
```

### Monitoring

```bash
squeue -u $USER
tail -f logs/preprocess_<JOBID>.out
scancel <JOBID>
```

### Default resource allocations

| Step | CPUs | Memory | Time |
|------|-----:|-------:|-----:|
| Preprocessing | 8 | 64 GB | 24 h |
| Phenotyping | 8 | 64 GB | 24 h |
| Spatial analysis | 16 | 128 GB | 48 h |
| Full pipeline | 16 | 128 GB | 96 h |
| Test suite | 4 | 32 GB | 2 h |

## Directory Structure

```
Xenium_Analysis/
├── README.md                          # this file
├── DATA_README.md                     # data organization
├── IMPLEMENTATION_SUMMARY.md          # high-level overview
├── environment.yml                    # conda runtime (scanpy, scvi, R, ...)
├── pyproject.toml                     # ruff + pytest config, project metadata
├── requirements-ci.txt                # minimal deps for GitHub Actions
├── config.ini                         # pipeline parameters template
├── .pre-commit-config.yaml            # local developer hooks
├── .github/
│   └── workflows/
│       └── ci.yml                     # GitHub Actions CI
├── notebooks/                         # analysis notebooks
│   ├── 01_preprocessing.ipynb
│   ├── 01_preprocessing_v2.ipynb
│   ├── 01_preprocessing_panc.ipynb
│   ├── 02_phenotyping.ipynb
│   ├── 03_spatial_analysis.ipynb
│   ├── 04_group_comparisons.ipynb
│   ├── 05_tissue_comparisons.ipynb
│   └── 06_xenium_phenocycler_integration.ipynb
├── scripts/
│   ├── export_xenium_explorer_groups.py  # CLI: .h5ad -> Xenium Explorer
│   ├── R/
│   │   └── xenium_analysis.R
│   ├── slurm/
│   │   ├── 01_run_preprocessing.sh
│   │   ├── 02_run_phenotyping.sh
│   │   ├── 03_run_spatial_analysis.sh
│   │   ├── run_full_pipeline.sh
│   │   └── run_tests.sh                  # HiPerGator test runner
│   └── (one-off refactor helpers)
├── utils/
│   ├── __init__.py                    # lazy-imports both submodules
│   ├── analysis_utils.py              # scanpy-based helpers
│   └── xenium_explorer_export.py      # Xenium Explorer export (scanpy-free)
├── tests/
│   ├── conftest.py                    # --run-slow / --run-hipergator flags
│   ├── test_repo_health.py            # env.yml, SLURM, notebooks, config.ini
│   ├── test_xenium_explorer_export.py
│   ├── test_deg_fixes.py              # structural + @slow functional
│   └── test_preprocessing_refactor.py # structural
├── data/         (git-ignored)
│   ├── raw/                           # input Xenium / Phenocycler
│   ├── processed/                     # pipeline outputs
│   └── phenocycler/
└── figures/      (git-ignored)
```

## Notebooks

| Notebook | Input | Output |
|---|---|---|
| `01_preprocessing.ipynb` | raw h5 / h5ad | `*_preprocessed.h5ad` |
| `02_phenotyping.ipynb` | `*_preprocessed.h5ad` | `*_annotated.h5ad` (pairs naturally with the [Xenium Explorer exporter](#exporting-clusters--phenotypes-to-xenium-explorer)) |
| `03_spatial_analysis.ipynb` | `*_annotated.h5ad` | `*_spatial_analysis.h5ad` |
| `04_group_comparisons.ipynb` | `*_spatial_analysis.h5ad` | DE tables, composition, volcanos |
| `05_tissue_comparisons.ipynb` | multiple `*_annotated.h5ad` | `integrated_tissues.h5ad` |
| `06_xenium_phenocycler_integration.ipynb` | Xenium + Phenocycler | `*_integrated.h5ad` |

See [IMPLEMENTATION_SUMMARY.md](IMPLEMENTATION_SUMMARY.md) for per-notebook
processing details.

## Utility Package (`utils`)

`utils/__init__.py` lazy-imports its submodules so the heavy scientific
stack is only loaded when you actually need it. That means
`utils.xenium_explorer_export` can be used in the lightweight CI
environment (no scanpy), while `utils.analysis_utils` pulls scanpy only
when referenced.

### `utils.analysis_utils` (scanpy / squidpy)

```python
from utils import (
    load_xenium_data,
    add_spatial_coordinates,
    calculate_qc_metrics,
    filter_cells_and_genes,
    normalize_and_hvg,
    run_standard_workflow,
    export_to_csv,
    create_summary_report,
)

adata = load_xenium_data(data_path, "sample_01")
adata = calculate_qc_metrics(adata)
adata = filter_cells_and_genes(adata)
adata = normalize_and_hvg(adata)
adata = run_standard_workflow(adata)
```

### `utils.xenium_explorer_export` (no scanpy)

See [Exporting Clusters / Phenotypes to Xenium Explorer](#exporting-clusters--phenotypes-to-xenium-explorer).

## Continuous Integration and Testing

Two complementary layers guard every commit.

### GitHub Actions — fast, every push and PR

Defined in `.github/workflows/ci.yml`. Runs on `ubuntu-latest` against
Python 3.10 and 3.12 with only `requirements-ci.txt` installed. Total
runtime is about one minute.

| Job | What it does |
|---|---|
| **lint** | `ruff check .` |
| **test** | `pytest --cov=utils` across the 3.10 / 3.12 matrix |
| **shellcheck** | Static analysis of `scripts/slurm/*.sh` |
| **environment-yml** | Parses `environment.yml` and checks the Python pin |

### HiPerGator — full suite, on demand

```bash
sbatch scripts/slurm/run_tests.sh
```

Activates the `xenium_analysis` conda env, verifies the core scientific
imports, then runs:

```bash
ruff check .
pytest tests/ -v
pytest tests/ -v --run-slow --run-hipergator
```

The `--run-slow` / `--run-hipergator` flags unlock tests that need real
Xenium data and the full conda env (see `tests/conftest.py`).

### Test layout

| File | Dependencies | Purpose |
|---|---|---|
| `tests/test_repo_health.py` | stdlib + pytest | Validates `environment.yml`, every SLURM script (shebang, LF line endings, resource directives, conda env), every notebook (valid JSON, `nbformat=4`, kernelspec), `utils` package imports cleanly without scanpy, `config.ini` sections |
| `tests/test_xenium_explorer_export.py` | `anndata`, `zarr` | CSV + `analysis.zarr.zip` roundtrip, color palette, rename, numeric clusters, error paths |
| `tests/test_preprocessing_refactor.py` | stdlib + pytest | Structural checks on preprocessing-v2 notebooks (skipped if notebooks not in checkout) |
| `tests/test_deg_fixes.py` | numpy; scanpy for `@pytest.mark.slow` tests | Structural + functional DEG fix verification |

### Running tests locally

```bash
pytest                                   # fast tests, ~1 s
pytest --run-slow --run-hipergator       # everything (needs real data)
pytest -k xenium_explorer -v             # a subset
pytest --cov=utils --cov-report=html     # coverage report -> htmlcov/
```

### Local pre-commit

```bash
pip install pre-commit
pre-commit install
```

On every `git commit` runs `ruff check`, `shellcheck` on SLURM scripts,
and trailing-whitespace / LF-line-ending hooks — a strict superset of
what CI does, so clean commits imply green CI.

## Configuration

`config.ini` is a template for per-project settings (sample list, QC
thresholds, marker genes, HiPerGator account, ...). Copy it and adjust:

```bash
cp config.ini my_project.ini
$EDITOR my_project.ini
```

The test suite validates the template parses and that numeric thresholds
are actually numeric — if you edit `config.ini` and CI starts failing,
`tests/test_repo_health.py::TestConfigIni` is where to look.

### Common knobs

**Cell-type markers** (`config.ini` → `[phenotyping]`, or `02_phenotyping.ipynb`):
```python
marker_genes = {
    "T cells": ["CD3D", "CD3E", "CD8A", "CD4"],
    "B cells": ["CD19", "CD79A", "MS4A1"],
    # add tissue-specific markers here
}
```

**QC thresholds** (`01_preprocessing.ipynb`):
```python
MIN_GENES = 10
MIN_CELLS = 3
MAX_MT_PCT = 20
MIN_COUNTS = 50
```

**Spatial parameters** (`03_spatial_analysis.ipynb`):
```python
sq.gr.spatial_neighbors(adata, n_neighs=10)
n_domains = 5
```

## Troubleshooting

**Memory errors on HiPerGator**
```bash
#SBATCH --mem=256gb  # or higher
```

**SLURM job fails immediately**
```bash
sacctmgr show assoc user=$USER  # confirm account / QOS
```

**Conda environment issues**
```bash
conda env remove -n xenium_analysis
conda env create -f environment.yml
```

**"No module named scanpy" in a script/notebook**
```bash
conda activate xenium_analysis
which python   # should point into the env
```

**Xenium Explorer does not show my Cell Groups CSV**  
The first column must be exactly `cell_id` and values must match the
original `cells.parquet`. See the [Exporting](#exporting-clusters--phenotypes-to-xenium-explorer)
section for details.

**CI failing on a PR but passing locally**  
Run the exact CI commands locally:
```bash
pip install -r requirements-ci.txt
ruff check .
pytest
```

## Citation

- **Scanpy** — Wolf et al., Genome Biology (2018)
- **Squidpy** — Palla et al., Nature Methods (2022)
- **scvi-tools** — Lopez et al., Nature Methods (2018)
- **Seurat** — Hao et al., Cell (2021)
- **10x Xenium Explorer** — 10x Genomics (https://www.10xgenomics.com/support/software/xenium-explorer)

## Support

- Issues: https://github.com/smith6jt-cop/Xenium_Analysis/issues
- HiPerGator docs: https://help.rc.ufl.edu/
- Scanpy tutorials: https://scanpy-tutorials.readthedocs.io/

## License

Provided as-is for research purposes.

## Acknowledgments

- University of Florida Research Computing and the HiPerGator cluster
- The scanpy, squidpy, and scvi-tools teams
- The `sopa` / `spatialdata-io` projects, whose `analysis.zarr.zip`
  layout this repo's export module mirrors
