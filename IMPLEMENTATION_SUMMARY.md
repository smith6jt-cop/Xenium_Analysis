# Xenium Analysis Pipeline — Implementation Summary

## Project Overview

This repository provides a production-ready pipeline for analyzing 10x
Xenium spatial transcriptomics data, with optional Phenocycler (CODEX)
integration, tooling for Xenium Explorer, R/Seurat follow-up analysis,
and SLURM infrastructure for the University of Florida HiPerGator cluster.

## What Has Been Implemented

### Core Analysis Pipeline

> **For pickup after disconnect, see [`HANDOFF.md`](HANDOFF.md) first.**
> For architectural decisions and structure, see [`CLAUDE.md`](CLAUDE.md).

#### 1. Preprocessing (`01_preprocessing_v2.ipynb`) — **OPTIONAL**
- Input: spatialdata zarr (`data/raw/{sample}.zarr`)
- Processing (10x-aligned QC): drop control probes, `min_counts=50`,
  `min_genes=20`, adaptive `max_counts=98%-ile`, `min_cells=100`;
  `normalize_total` (median target), `log1p`; HVG (seurat_v3, metadata
  only — PCA uses all panel genes); **GPU `rsc.pp.scale(zero_center=False)`**;
  **GPU `rsc.pp.pca(n_comps=30)`**; **GPU `rsc.pp.neighbors(k=15, cosine)`**;
  CPU `sc.tl.umap(min_dist=0.5)` (rapids UMAP outliers — see CLAUDE.md);
  CPU `sc.tl.leiden(resolution=1.5, flavor='igraph', n_iterations=-1)`
  (cugraph leiden has dask-cuda 24.12 vs dask 2026.1 version drift).
  Squidpy decorative stats (centrality / nhood_enrichment permutations /
  Moran's I) are **skipped** at this scale.
- Cell 1 enables B200 TF32 (`torch.set_float32_matmul_precision('high')`).
- Output: `data/processed/{sample}/{sample}_preprocessed.h5ad` (~100 GB)
  with `leiden_1.5` in obs, `counts`/`lognorm` in layers, `X_pca`/`X_umap`
  in obsm — verified output contract for stage 02.
- **Not required**: stage 02 cell 3 ingests directly from zarr if the
  preprocessed h5ad is missing.
- Legacy variants (`01_preprocessing.ipynb`,
  `01_preprocessing_panc.ipynb`) live under `notebooks/legacy/`.

#### 2. Phenotyping (`02_phenotyping.ipynb`) — **PRIMARY ENTRY POINT**
- Input: zarr OR `_preprocessed.h5ad` (auto-detected in cell 3)
- Panel-aware pancreas marker dict (cell 10) — explicitly omits
  `INS / GCG / SST / PPY / TTR / IRX2` and most acinar zymogens that are
  absent from the T1D add-on panel; uses surrogates for Beta/Alpha/
  Delta/Epsilon. `LOW_CONFIDENCE_TYPES = {Alpha, Delta, Epsilon, Acinar}`.
- Score-based assignment with **anti-acinar override** (cell 15) —
  the Xenium panel has only `AMY1A` for Acinar, so an extra pass
  reassigns Acinar-called cells whenever curated specific markers of
  another lineage (Beta/Endocrine/Ductal/Stellate/Endothelial/
  T_cells/B_cells/Myeloid/Schwann) exceed their per-type baseline.
- Confidence-based Indeterminate flagging (median - 1·MAD on
  `celltype_confidence`) for low-coverage types — runs AFTER override.
- scVI training (cell 17): `layer='counts'`, `n_latent=30`,
  `gene_likelihood='nb'`, `continuous_covariate_keys=['cell_area']`
  (NO `total_counts` — NB likelihood already models library size),
  **GPU on B200 with TF32 / Tensor Cores enabled** (cell 1 sets
  `torch.set_float32_matmul_precision('high')`, ~30-50% scVI speedup),
  `max_epochs=200` with early stopping.
- **GPU `rsc.pp.neighbors(use_rep='X_scvi', cosine)`** (cuML kNN —
  ~15s on 1.2M cells vs ~15 min CPU pynndescent), then CPU
  `sc.tl.umap` + CPU `sc.tl.leiden(resolution=1.0, flavor='igraph')`.
- Quantitative validation: Spearman ρ between each UMAP axis and
  `total_counts` rendered to `umap_qc_overlay.png`. Both samples
  currently < 0.3 (passing).
- **GPU `rsc.pp.neighbors(use_rep='X_pca')`** for the PCA-fallback
  embedding stored as `X_umap_pca` / `leiden_pca`.
- Output: `data/processed/{sample}/{sample}_phenotyped.h5ad` (~60 GB)
- Recommended follow-up: call `utils.export_for_xenium_explorer()` to push
  clusters / phenotypes straight into Xenium Explorer (see the Xenium
  Explorer Export Module section below)

#### 3. Spatial Analysis (`03_spatial_analysis.ipynb`)
- Cell 1 imports rapids/cuML (cuKMeans, cuNearestNeighbors, rsc).
- Cell 2 slim adata (drop heavy obsm/obsp/uns we don't need at this stage)
- Cell 4 `sq.gr.spatial_neighbors(n_neighs=10, delaunay=True)` on full
  data — CPU, already fast.
- Cell 6 `sq.gr.nhood_enrichment(n_perms=100, n_jobs=1, seed=0)` —
  `n_jobs>1` hits a numba/joblib readonly-array bug; sequential with
  reduced perms gives same statistical resolution
- Cell 8 co-occurrence on 150 k subsample (CPU squidpy)
- Cell 10 **GPU Moran's I (`rsc.gr.spatial_autocorr`)** — attempts full
  1.2M cells; 200k subsample retained as graceful fallback if OOM.
- Cell 12 spatial niches via **GPU `cuml.cluster.KMeans` (n_clusters=12,
  scalable-k-means++)** on composition vectors built from **GPU
  `cuml.neighbors.NearestNeighbors`** k=50 spatial graph + 2-pass
  numpy majority-vote smoothing (Banksy-style).
- Cell 14 `sq.gr.ligrec` on 100 k subsample, n_perms=50
- Cell 16 `sq.gr.ripley` on 100 k subsample
- Cell 18 defensive try/except CSV exports + h5ad save
- Output: `data/processed/{sample}/{sample}_spatial_analysis.h5ad` (~10 GB)
  + per-sample CSVs (ligrec means/pvalues/metadata, neighborhood_enrichment,
  moranI, spatial_niches)

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

#### GPU stack (HiPerGator B200 node)

- `rapids_singlecell == 0.14.1` (required for scanpy 1.12 API compat)
- `cudf-cu12 == 24.12`, `cuml-cu12 == 24.12`, `rmm-cu12 == 24.12.1`
- `cupy-cuda12x == 14.0.1`, `numba-cuda == 0.0.17.1`,
  `pynvjitlink-cu12 == 0.7.0`
- `torch == 2.10.0` with CUDA 12.9, cudnn 91002 (B200 native)

**Critical env vars set by `scripts/run_local_pipeline.sh`** —
required for the rapids stack to import:

```bash
export CUDA_HOME=/apps/compilers/cuda/12.8.1
export LD_LIBRARY_PATH="$CUDA_HOME/lib64:$CUDA_HOME/nvvm/lib64:$LD_LIBRARY_PATH"
export LD_PRELOAD="$CONDA_PREFIX/lib/libstdc++.so.6:$CUDA_HOME/lib64/libcudart.so.12"
```

The `LD_PRELOAD` line pins runtime libcudart to 12.8 (matching the
driver), which makes `cudf 24.12._setup_numba()` skip the
MVC patch path that would otherwise call `pynvjitlink.patch_numba_linker()`
and raise `RuntimeError: numba_cuda includes patches from pynvjitlink`.
Conda's `libstdc++.so.6` is preloaded first to avoid the older system
`/lib64/libstdc++.so.6` (no CXXABI_1.3.15) being pulled in alongside
the cuda toolkit lib dir, which would break matplotlib.

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

## File Inventory (post-2026-04-25 cleanup)

- Notebooks (active): 7 (`00_ingest`, `01_preprocessing_v2`, `02`–`06`)
- Notebooks (legacy, deprecated): 2 under `notebooks/legacy/`
  (`01_preprocessing.ipynb`, `01_preprocessing_panc.ipynb`)
- SLURM scripts: 5 (4 pipeline + `run_tests.sh`)
- Python scripts (active): 3 — `run_local_pipeline.sh`,
  `slim_phenotyped.py`, `export_xenium_explorer_groups.py`
- Python scripts (legacy, deprecated): 13 under `scripts/legacy/`
  (one-shot fix-ups from iterative debugging — see
  `scripts/legacy/README.md` for what each replaced)
- Utils: `utils/` (2 submodules)
- Tests: `tests/` (4 test modules + conftest)
- CI: GitHub Actions workflow, pre-commit config, pyproject.toml,
  requirements-ci.txt
- R scripts: `scripts/R/xenium_analysis.R`
- Documentation:
  - `HANDOFF.md` — per-session pickup notes (read first after disconnect)
  - `CLAUDE.md` — repo guide for Claude / new contributors
  - `README.md` — user-facing overview
  - `DATA_README.md` — expected data formats + Xenium Explorer outputs
  - `IMPLEMENTATION_SUMMARY.md` — this file
  - `config.ini` — parameter template (panel-aware pancreas markers)

## Future Extensions

- Additional spatial analysis methods (spatially-aware DE, CellCharter)
- Interactive dashboards (holoviews / panel)
- Docker / Apptainer image with the full conda env pre-built
- Example Xenium dataset for the test suite

## Support

- GitHub Issues — bug reports, feature requests
- HiPerGator documentation — https://help.rc.ufl.edu/
