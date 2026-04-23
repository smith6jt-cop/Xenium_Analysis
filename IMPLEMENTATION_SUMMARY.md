# Xenium Analysis Pipeline — Implementation Summary

## Project Overview

This repository provides a production-ready pipeline for analyzing 10x
Xenium spatial transcriptomics data, with optional Phenocycler (CODEX)
integration, tooling for Xenium Explorer, R/Seurat follow-up analysis,
and SLURM infrastructure for the University of Florida HiPerGator cluster.

## What Has Been Implemented

### Core Analysis Pipeline (6 Jupyter Notebooks)

#### 1. Preprocessing (`01_preprocessing.ipynb`, plus `_v2` and `_panc` variants)
- Input: raw h5 / h5ad files from Xenium
- Processing: QC metrics (genes/cell, counts, mt%), cell/gene filtering,
  normalization (log1p), highly-variable genes, PCA, UMAP, Leiden
- Output: `*_preprocessed.h5ad`

#### 2. Phenotyping (`02_phenotyping.ipynb`)
- Marker-gene identification and scoring
- Automated annotation with tunable markers
- scVI-based advanced clustering
- Output: `*_annotated.h5ad`
- Recommended follow-up: call `utils.export_for_xenium_explorer()` to push
  clusters / phenotypes straight into Xenium Explorer (see the Xenium
  Explorer Export Module section below)

#### 3. Spatial Analysis (`03_spatial_analysis.ipynb`)
- Spatial neighborhood graph
- Neighborhood enrichment, co-occurrence, Moran's I
- Spatial domain identification
- Ligand-receptor analysis, Ripley's L
- Output: `*_spatial_analysis.h5ad`

#### 4. Group Comparisons (`04_group_comparisons.ipynb`)
- Cell-type composition, statistical testing
- Differential expression between groups (Xenium-appropriate LFC from
  log-normalized counts, no artificial p-value cap, `adjustText` labels)
- Volcano plots, gene-set export
- Output: DE tables, composition tables

#### 5. Tissue Comparisons (`05_tissue_comparisons.ipynb`)
- Multi-sample integration with scVI batch correction
- Cross-tissue DE
- Output: `integrated_tissues.h5ad`

#### 6. Phenocycler Integration (`06_xenium_phenocycler_integration.ipynb`)
- Coordinate alignment, nearest-neighbor cross-modal cell mapping
- Output: `*_integrated.h5ad`

### Xenium Explorer Export Module

`utils/xenium_explorer_export.py` (new) — 400 LOC, scanpy-free.

Exports any `.obs` categorical (clusters, phenotypes, scVI labels, ...)
to the two formats Xenium Explorer 3.0+ accepts:

| Function | Emits | Loaded in Explorer via |
|---|---|---|
| `export_groups_to_csv` | Cell Groups CSV matched by `cell_id` | Cell → *Add cell categorization* |
| `export_groups_to_zarr` | `analysis.zarr.zip` (zarr v2, position-indexed) | Drop next to `experiment.xenium` — auto-loaded |
| `export_for_xenium_explorer` | Both of the above + palette JSON | — |
| `generate_color_palette` | `{category: '#RRGGBB'}` map with matplotlib-free fallback | — |

A CLI wrapper (`scripts/export_xenium_explorer_groups.py`) exposes the
same functionality for non-notebook use.

### CI/CD and Testing Infrastructure

Two complementary layers:

#### GitHub Actions (`.github/workflows/ci.yml`, ~60s runtime)
Runs on every push / PR using only `requirements-ci.txt` — no scanpy:
- **lint** — `ruff check .`
- **test** — `pytest --cov=utils` on Python 3.10 and 3.12
- **shellcheck** — static analysis of `scripts/slurm/*.sh`
- **environment-yml** — validates `environment.yml`

#### HiPerGator test runner (`scripts/slurm/run_tests.sh`)
SLURM job that activates the full conda env and runs both fast tests
and `@pytest.mark.slow` / `@pytest.mark.hipergator` tests against real
Xenium data.

#### Test modules
- `tests/test_repo_health.py` — validates `environment.yml`, SLURM
  scripts (shebang, LF endings, resource directives, conda env), every
  notebook (valid JSON, `nbformat=4`, kernelspec), `utils` imports, and
  `config.ini` sections
- `tests/test_xenium_explorer_export.py` — CSV + `.zarr.zip` roundtrip,
  palette, rename, numeric clusters, error paths
- `tests/test_deg_fixes.py` — structural (always) + functional
  `@pytest.mark.slow` (HiPerGator only)
- `tests/test_preprocessing_refactor.py` — structural notebook checks
  with graceful skips when optional notebooks aren't present
- `tests/conftest.py` — `--run-slow` / `--run-hipergator` opt-in flags

#### Developer ergonomics
- `pyproject.toml` centralizes ruff + pytest config
- `.pre-commit-config.yaml` runs the same checks CI does, on every commit

### HiPerGator Integration

SLURM scripts under `scripts/slurm/`:

| Script | CPUs | Memory | Time |
|---|---:|---:|---:|
| `01_run_preprocessing.sh` | 8 | 64 GB | 24 h |
| `02_run_phenotyping.sh` | 8 | 64 GB | 24 h |
| `03_run_spatial_analysis.sh` | 16 | 128 GB | 48 h |
| `run_full_pipeline.sh` | 16 | 128 GB | 96 h |
| `run_tests.sh` (new) | 4 | 32 GB | 2 h |

All scripts include automatic conda activation, resource allocation,
email notifications, and error logging.

### R Integration (`scripts/R/xenium_analysis.R`)

- `load_xenium_h5ad()` — h5ad → Seurat
- `run_de_analysis()` — Seurat-based DE
- `plot_volcano()`, `plot_spatial_features()`, `plot_composition()`

### Utility Package (`utils/`)

`utils/__init__.py` uses **lazy imports** so heavy deps are only loaded
when actually needed. This means `utils.xenium_explorer_export` can run
in the lightweight CI environment without scanpy.

- `utils/analysis_utils.py` — scanpy / squidpy helpers (load, QC, filter,
  normalize+HVG, standard workflow, CSV export, summary report)
- `utils/xenium_explorer_export.py` — Xenium Explorer export (scanpy-free)

### Environment

`environment.yml` pins **Python 3.12** and all runtime dependencies
(scanpy, squidpy, scvi-tools, R/Seurat, rpy2, napari, spatialdata, ...).
`requirements-ci.txt` lists the minimal lightweight set for CI.

`pyproject.toml` declares Python 3.10 as the project minimum so
contributors with older systems can still run the lint and test suite.

### Documentation

| File | Purpose |
|---|---|
| `README.md` | Main guide: install, workflow, Xenium Explorer export, HiPerGator, CI, config, troubleshooting |
| `DATA_README.md` | Data organization, expected formats, conversion, Xenium Explorer outputs |
| `IMPLEMENTATION_SUMMARY.md` | This file — high-level inventory |
| `config.ini` | Parameter template (60+ knobs) |

## Technology Stack

### Python ecosystem
- **scanpy**, **squidpy**, **scvi-tools**, **anndata** — analysis
- **numpy / pandas / scipy / scikit-learn** — base
- **matplotlib / seaborn / plotly / napari** — visualization
- **zarr** — Xenium Explorer `analysis.zarr.zip` format

### R ecosystem
- **Seurat**, **ggplot2**, **dplyr/tidyr** via **rpy2**

### HPC + Dev
- **SLURM** — scheduling (HiPerGator)
- **conda** — environment management
- **GitHub Actions** — CI
- **ruff**, **pytest**, **shellcheck**, **pre-commit** — quality gates

## Key Features

### End-to-end analysis
Raw Xenium data → publication-ready figures, with every intermediate
step reproducible from saved h5ad files.

### Xenium Explorer round-trip
Analysis results in Python → back into the viewer that biologists
actually use, with two format options covering every filtering scenario.

### HPC-ready
SLURM scripts with sensible defaults, optimized resource allocations,
and a dedicated test-runner job.

### Tested and linted
Repo health checked on every commit: notebook JSON validity, SLURM
script correctness, environment pins, config integrity, utility
imports. No silent drift.

### Reproducible
Pinned conda environment, version-controlled pipeline, documented
parameters. CI enforces format and structure on every PR.

## Quality Assurance

### Code
- Every notebook validated as JSON with `nbformat=4` and a Python kernelspec
- `ruff` lint clean across `utils/`, `scripts/`, `tests/`
- SLURM scripts pass `shellcheck`
- `utils/` package imports without scanpy (lazy-imported)

### Tests
- Repo-health tests catch environment / notebook / SLURM regressions
  before they reach HiPerGator
- Xenium Explorer export covered by 8 dedicated unit tests
- DEG functional tests gated behind `@pytest.mark.slow` for
  HiPerGator-only execution

### Documentation
- Feature-level README covering each subsystem
- First-class section for the Xenium Explorer export feature with API
  reference, format comparison, CLI + notebook usage, and
  troubleshooting
- Data organization guide
- Inline docstrings on every public function

## Usage Patterns

### Local analysis
```bash
conda activate xenium_analysis
jupyter lab
```

### HiPerGator
```bash
sbatch scripts/slurm/run_full_pipeline.sh
squeue -u $USER
```

### Xenium Explorer export (post-phenotyping)
```python
from utils import export_for_xenium_explorer
export_for_xenium_explorer(
    adata,
    group_keys=["leiden_0.5", "celltype"],
    output_dir="data/xenium_explorer",
    sample_name="THYHDL065",
)
```

### CI locally
```bash
pip install -r requirements-ci.txt
ruff check .
pytest
```

## File Inventory

- Notebooks: 8 (`01_preprocessing` + two variants, `02`–`06`)
- SLURM scripts: 5 (4 pipeline + `run_tests.sh`)
- Python packages: `utils/` (2 submodules), `scripts/` (CLI + refactor
  helpers), `tests/` (4 test modules + conftest)
- CI: GitHub Actions workflow, pre-commit config, pyproject.toml,
  requirements-ci.txt
- R scripts: `scripts/R/xenium_analysis.R`
- Documentation: `README.md`, `DATA_README.md`,
  `IMPLEMENTATION_SUMMARY.md`, `config.ini`

## Future Extensions

- Additional spatial analysis methods (spatially-aware DE, CellCharter)
- Interactive dashboards (holoviews / panel)
- Docker / Apptainer image with the full conda env pre-built
- Example Xenium dataset for the test suite

## Support

- GitHub Issues — bug reports, feature requests
- HiPerGator documentation — https://help.rc.ufl.edu/
