# Xenium Spatial Transcriptomics Analysis Pipeline

A comprehensive, jupyter-notebook-based analysis pipeline for Xenium spatial transcriptomics data with Phenocycler integration. Optimized for HiPerGator HPC cluster and includes scvi-tools, scanpy, squidpy, and R integration.

## Table of Contents
- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Workflow](#workflow)
- [HiPerGator Usage](#hipergator-usage)
- [Directory Structure](#directory-structure)
- [Notebooks](#notebooks)
- [Scripts](#scripts)
- [Citation](#citation)

## Features

### Analysis Capabilities
- ✅ **Preprocessing**: QC, filtering, normalization from h5 files
- ✅ **Phenotyping**: Cell type annotation with marker genes and scVI
- ✅ **Spatial Analysis**: Neighborhood enrichment, cell-cell interactions, spatial statistics
- ✅ **Group Comparisons**: Differential expression and composition analysis
- ✅ **Tissue Comparisons**: Multi-sample integration and comparison
- ✅ **Phenocycler Integration**: Multi-modal data integration

### Technologies
- **Python**: scanpy, squidpy, scvi-tools, spatialdata
- **R**: Seurat, ggplot2, rpy2 integration
- **HPC**: SLURM scripts for HiPerGator
- **Visualization**: matplotlib, seaborn, plotly, napari

## Installation

### 1. Clone Repository

```bash
git clone https://github.com/smith6jt-cop/Xenium_Analysis.git
cd Xenium_Analysis
```

### 2. Create Conda Environment

```bash
conda env create -f environment.yml
conda activate xenium_analysis
```

This creates an environment with all required packages including:
- Python 3.10
- scanpy, squidpy, scvi-tools
- R and Seurat
- Jupyter Lab
- All dependencies

### 3. Verify Installation

```bash
python -c "import scanpy as sc; import squidpy as sq; import scvi; print('All packages loaded successfully!')"
```

## Quick Start

### 1. Prepare Your Data

Place your Xenium h5 files in the `data/raw/` directory:

```bash
data/raw/
├── xenium_sample_01.h5ad
├── xenium_sample_02.h5ad
└── ...
```

### 2. Run Analysis Locally

```bash
# Activate environment
conda activate xenium_analysis

# Launch Jupyter Lab
jupyter lab

# Open and run notebooks in order:
# 01_preprocessing.ipynb
# 02_phenotyping.ipynb
# 03_spatial_analysis.ipynb
# ...
```

### 3. Run on HiPerGator

```bash
# Submit individual jobs
sbatch scripts/slurm/01_run_preprocessing.sh

# Or run complete pipeline
sbatch scripts/slurm/run_full_pipeline.sh
```

## Workflow

The analysis pipeline consists of 6 sequential notebooks:

```
Raw h5 files
     ↓
[01_preprocessing.ipynb]
     ↓
Filtered & normalized data
     ↓
[02_phenotyping.ipynb]
     ↓
Annotated cell types
     ↓
[03_spatial_analysis.ipynb]
     ↓
Spatial metrics & interactions
     ↓
[04_group_comparisons.ipynb]
     ↓
Differential expression
     ↓
[05_tissue_comparisons.ipynb]
     ↓
Integrated multi-tissue data
     ↓
[06_xenium_phenocycler_integration.ipynb]
     ↓
Final integrated dataset
```

## HiPerGator Usage

### Setup on HiPerGator

1. **Load conda module**:
```bash
module load conda
```

2. **Create environment**:
```bash
conda env create -f environment.yml
```

3. **Update SLURM scripts**:

Edit the following in each SLURM script:
- `--mail-user`: Your UFL email
- `--qos`: Your QOS allocation
- `--account`: Your account/group

### Submit Jobs

**Individual notebooks**:
```bash
# Create logs directory
mkdir -p logs

# Submit jobs
sbatch scripts/slurm/01_run_preprocessing.sh
sbatch scripts/slurm/02_run_phenotyping.sh
sbatch scripts/slurm/03_run_spatial_analysis.sh
```

**Full pipeline**:
```bash
sbatch scripts/slurm/run_full_pipeline.sh
```

### Monitor Jobs

```bash
# Check job status
squeue -u $USER

# Check job output
tail -f logs/preprocess_JOBID.out

# Cancel job
scancel JOBID
```

### Resource Requirements

| Step | CPUs | Memory | Time | Notes |
|------|------|--------|------|-------|
| Preprocessing | 8 | 64GB | 24h | Basic QC and filtering |
| Phenotyping | 8 | 64GB | 24h | Includes scVI training |
| Spatial Analysis | 16 | 128GB | 48h | Computationally intensive |
| Full Pipeline | 16 | 128GB | 96h | All steps sequentially |

## Directory Structure

```
Xenium_Analysis/
├── README.md                          # This file
├── environment.yml                    # Conda environment
├── notebooks/                         # Analysis notebooks
│   ├── 01_preprocessing.ipynb
│   ├── 02_phenotyping.ipynb
│   ├── 03_spatial_analysis.ipynb
│   ├── 04_group_comparisons.ipynb
│   ├── 05_tissue_comparisons.ipynb
│   └── 06_xenium_phenocycler_integration.ipynb
├── scripts/
│   ├── R/                            # R analysis scripts
│   │   └── xenium_analysis.R
│   └── slurm/                        # HiPerGator SLURM scripts
│       ├── 01_run_preprocessing.sh
│       ├── 02_run_phenotyping.sh
│       ├── 03_run_spatial_analysis.sh
│       └── run_full_pipeline.sh
├── utils/                            # Utility functions
│   └── analysis_utils.py
├── data/
│   ├── raw/                          # Raw h5 files (not tracked)
│   ├── processed/                    # Processed data (not tracked)
│   └── phenocycler/                  # Phenocycler data (not tracked)
└── figures/                          # Generated figures (not tracked)
    ├── 01_preprocessing/
    ├── 02_phenotyping/
    ├── 03_spatial_analysis/
    ├── 04_group_comparisons/
    ├── 05_tissue_comparisons/
    └── 06_integration/
```

## Notebooks

### 01_preprocessing.ipynb
- Load Xenium h5 files
- Quality control metrics
- Cell and gene filtering
- Normalization and log-transformation
- Highly variable gene selection
- Dimensionality reduction (PCA, UMAP)
- Initial clustering

**Input**: Raw h5/h5ad files  
**Output**: `*_preprocessed.h5ad`

### 02_phenotyping.ipynb
- Marker gene analysis
- Cell type scoring
- Automated annotation
- scVI-based clustering
- Manual curation support
- Final cell type assignment

**Input**: `*_preprocessed.h5ad`  
**Output**: `*_annotated.h5ad`

### 03_spatial_analysis.ipynb
- Spatial neighborhood graph
- Neighborhood enrichment
- Co-occurrence analysis
- Spatial autocorrelation (Moran's I)
- Spatial domains
- Ligand-receptor interactions
- Ripley's statistics

**Input**: `*_annotated.h5ad`  
**Output**: `*_spatial_analysis.h5ad`

### 04_group_comparisons.ipynb
- Cell type composition analysis
- Differential expression between groups
- Cell type-specific DE
- Volcano plots
- Statistical testing

**Input**: `*_spatial_analysis.h5ad`  
**Output**: DE results, composition tables

### 05_tissue_comparisons.ipynb
- Multi-tissue integration
- Batch correction with scVI
- Cross-tissue differential expression
- Tissue-specific spatial patterns
- Comparative analysis

**Input**: Multiple `*_annotated.h5ad` files  
**Output**: `integrated_tissues.h5ad`

### 06_xenium_phenocycler_integration.ipynb
- Phenocycler data preprocessing
- Spatial alignment
- Multi-modal cell mapping
- Cross-modal analysis
- Integrated visualization

**Input**: Xenium + Phenocycler data  
**Output**: `*_integrated.h5ad`

## Scripts

### R Integration (`scripts/R/xenium_analysis.R`)

Advanced statistical analysis and visualization:

```r
source("scripts/R/xenium_analysis.R")

# Load data
seurat_obj <- load_xenium_h5ad("data/processed/sample_annotated.h5ad")

# Differential expression
de_results <- run_de_analysis(seurat_obj, "condition", "Treatment", "Control")

# Visualizations
plot_volcano(de_results)
plot_spatial_features(seurat_obj, c("CD3D", "CD8A"))
plot_composition(seurat_obj, "celltype", "condition")
```

### SLURM Scripts

All SLURM scripts support:
- Email notifications
- Resource allocation
- Error logging
- Automatic output organization

Modify these variables in each script:
```bash
#SBATCH --mail-user=your_email@ufl.edu
#SBATCH --qos=your_qos
#SBATCH --account=your_account
```

## Exporting Clusters / Phenotypes to Xenium Explorer

After phenotyping, you can push any `.obs` categorical (Leiden clusters,
scVI labels, curated phenotypes, ...) back to Xenium Explorer 3.0+ so cells
are colored by that grouping in the viewer.

### From a notebook

```python
from utils import export_for_xenium_explorer

export_for_xenium_explorer(
    adata,
    group_keys=["leiden_0.5", "celltype", "phenotype"],
    output_dir="data/xenium_explorer",
    sample_name="THYHDL065",
    cell_id_key="cell_id",   # defaults to common candidates or the index
)
```

Artifacts written:

- `THYHDL065_cell_groups.csv` – load via **Cell → Add cell categorization**.
- `THYHDL065/analysis.zarr.zip` – drop next to `experiment.xenium` so the
  groupings show up automatically in the **Clusters** dropdown.
- `THYHDL065_color_palette.json` – the category colors used, to keep
  matplotlib figures consistent with the Explorer view.

### From the CLI

```bash
python scripts/export_xenium_explorer_groups.py \
    --input data/processed/THYHDL065_annotated.h5ad \
    --output data/xenium_explorer \
    --sample THYHDL065 \
    --groups leiden_0.5 celltype phenotype \
    --cell-id-key cell_id
```

> The CSV format matches cells by `cell_id` (safe after filtering).  The
> `analysis.zarr.zip` is position-indexed, so use it only when the AnnData
> row order matches the original `cells.parquet` order.

## Continuous Integration

Two complementary layers of automated testing guard the repo:

### GitHub Actions (lightweight, every push / PR)

Defined in `.github/workflows/ci.yml`.  Runs on `ubuntu-latest` against
Python 3.10 and 3.12 (matches the `environment.yml` pin).  Four jobs:

| Job | What it does |
|-----|--------------|
| **lint** | `ruff check` on `utils/`, `scripts/`, `tests/` |
| **test** | `pytest --cov=utils` with only `requirements-ci.txt` installed (no scanpy/scvi-tools) |
| **shellcheck** | Static analysis of the SLURM scripts under `scripts/slurm/` |
| **environment-yml** | Verifies `environment.yml` parses and pins Python 3.10+ |

Total runtime: ~60s.  Heavy deps are deliberately not installed so PRs
don't wait ten minutes on every push.

### HiPerGator (full integration, on demand)

```bash
sbatch scripts/slurm/run_tests.sh
```

Activates the `xenium_analysis` conda env, runs lint + the fast suite,
then the `@pytest.mark.slow` + `@pytest.mark.hipergator` tests that need
real Xenium data.

### Local pre-commit

```bash
pip install pre-commit && pre-commit install
```

Runs the same `ruff check`, shellcheck, trailing-whitespace, and
LF-line-ending hooks on every commit — matches what CI enforces.

### Running tests locally

```bash
# Fast tests only (what GitHub Actions runs)
pytest

# Include functional tests that need real data (HiPerGator)
pytest --run-slow --run-hipergator
```

## Utility Functions

The `utils/analysis_utils.py` module provides reusable functions:

```python
from utils.analysis_utils import *

# Load data
adata = load_xenium_data(data_path, "sample_01")

# Add spatial coordinates
adata = add_spatial_coordinates(adata, coords_file)

# Calculate QC metrics
adata = calculate_qc_metrics(adata)

# Filter
adata = filter_cells_and_genes(adata)

# Normalize and identify HVGs
adata = normalize_and_hvg(adata)

# Run standard workflow
adata = run_standard_workflow(adata)

# Export results
export_to_csv(adata, output_dir, "sample_01")

# Create summary
summary = create_summary_report(adata, "sample_01", "summary.csv")
```

## Customization

### Cell Type Markers

Edit marker genes in `02_phenotyping.ipynb`:

```python
marker_genes = {
    'T cells': ['CD3D', 'CD3E', 'CD8A', 'CD4'],
    'B cells': ['CD19', 'CD79A', 'MS4A1'],
    # Add your tissue-specific markers
}
```

### QC Thresholds

Adjust in `01_preprocessing.ipynb`:

```python
MIN_GENES = 10
MIN_CELLS = 3
MAX_MT_PCT = 20
MIN_COUNTS = 50
```

### Spatial Parameters

Modify in `03_spatial_analysis.ipynb`:

```python
# Spatial neighbors
sq.gr.spatial_neighbors(adata, n_neighs=10)

# Spatial domains
n_domains = 5
```

## Troubleshooting

### Common Issues

**1. Memory errors on HiPerGator**
```bash
# Increase memory allocation
#SBATCH --mem=256gb
```

**2. SLURM job fails immediately**
```bash
# Check account and QOS
sacctmgr show assoc user=$USER
```

**3. Conda environment issues**
```bash
# Remove and recreate
conda env remove -n xenium_analysis
conda env create -f environment.yml
```

**4. Module not found errors**
```bash
# Verify environment is activated
conda activate xenium_analysis
which python
```

## Best Practices

1. **Start with small test dataset** to verify pipeline
2. **Save intermediate results** after each major step
3. **Use version control** for custom modifications
4. **Document parameters** used for reproducibility
5. **Keep raw data separate** from processed outputs
6. **Regular backups** of analysis results

## Data Organization

Recommended file naming:
```
<project>_<tissue>_<condition>_<replicate>_<step>.h5ad

Examples:
study1_liver_control_rep1_preprocessed.h5ad
study1_liver_treatment_rep1_annotated.h5ad
```

## Citation

If you use this pipeline, please cite:

- **Scanpy**: Wolf et al., Genome Biology (2018)
- **Squidpy**: Palla et al., Nature Methods (2022)
- **scvi-tools**: Lopez et al., Nature Methods (2018)
- **Seurat**: Hao et al., Cell (2021)

## Support

For questions and issues:
- Open an issue on GitHub
- Check HiPerGator documentation: https://help.rc.ufl.edu/
- Scanpy tutorials: https://scanpy-tutorials.readthedocs.io/

## License

This pipeline is provided as-is for research purposes.

## Acknowledgments

- University of Florida Research Computing
- HiPerGator HPC cluster
- Scanpy, Squidpy, and scvi-tools development teams
