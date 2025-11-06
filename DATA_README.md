# Data Directory Structure

This document describes the expected data structure for the Xenium analysis pipeline.

## Directory Overview

```
data/
├── raw/              # Raw Xenium and Phenocycler files
├── processed/        # Intermediate and final processed files
└── phenocycler/      # Phenocycler-specific data
```

## Raw Data (`data/raw/`)

### Xenium Data

Place your Xenium h5 files here. The pipeline supports multiple formats:

**Option 1: AnnData h5ad format (Recommended)**
```
data/raw/
└── xenium_sample_01.h5ad
```

**Option 2: 10X h5 format**
```
data/raw/
└── xenium_sample_01.h5
```

**Option 3: Xenium output directory**
```
data/raw/
└── xenium_sample_01/
    ├── cells.csv.gz
    ├── cell_feature_matrix.h5
    ├── cell_feature_matrix/
    └── ...
```

### File Naming Convention

Use descriptive, consistent naming:
```
<project>_<tissue>_<condition>_<replicate>.h5ad

Examples:
study1_liver_control_rep1.h5ad
study1_liver_treatment_rep1.h5ad
study2_kidney_healthy_rep1.h5ad
```

## Processed Data (`data/processed/`)

The pipeline automatically generates processed files at each step:

```
data/processed/
├── sample_01_preprocessed.h5ad          # After 01_preprocessing.ipynb
├── sample_01_annotated.h5ad             # After 02_phenotyping.ipynb
├── sample_01_spatial_analysis.h5ad      # After 03_spatial_analysis.ipynb
├── sample_01_celltype_assignments.csv   # Cell type annotations
├── sample_01_marker_genes.csv           # Marker genes per cluster
├── sample_01_de_genes_by_group.csv      # Differential expression results
└── integrated_tissues.h5ad              # After multi-tissue integration
```

## Phenocycler Data (`data/phenocycler/`)

### Expected Format

Phenocycler/CODEX data can be in various formats:

**Option 1: AnnData format (Recommended)**
```
data/phenocycler/
└── sample_01_phenocycler.h5ad
```

**Option 2: CSV/TSV format**
```
data/phenocycler/
├── sample_01_expression.csv       # Protein expression matrix
├── sample_01_metadata.csv         # Cell metadata
└── sample_01_coordinates.csv      # Spatial coordinates
```

### CSV Format Specifications

If providing CSV files, use the following structure:

**expression.csv**
```csv
cell_id,CD3,CD4,CD8,CD19,CD68,...
cell_0001,15.2,8.3,1.2,0.5,2.1,...
cell_0002,2.1,1.5,18.7,0.3,1.8,...
...
```

**metadata.csv**
```csv
cell_id,cluster,region,sample_id
cell_0001,T_cells,tumor,sample_01
cell_0002,B_cells,stroma,sample_01
...
```

**coordinates.csv**
```csv
cell_id,x,y
cell_0001,1523.4,2847.6
cell_0002,1528.9,2851.2
...
```

## Data Requirements

### Minimum Requirements

For Xenium data:
- Gene expression matrix (cells × genes)
- Spatial coordinates (x, y) for each cell
- Cell IDs

For Phenocycler data:
- Protein expression matrix (cells × proteins)
- Spatial coordinates (x, y) for each cell
- Cell IDs

### Recommended Additional Data

- Cell type annotations (if available)
- Sample metadata (condition, replicate, etc.)
- QC metrics
- Segmentation boundaries (for visualization)

## Example Data

For testing the pipeline, you can use publicly available datasets:

### Xenium Public Datasets

1. **10x Genomics Demo Data**
   - https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast
   - Download the h5 file to `data/raw/`

2. **Chan Zuckerberg CELLxGENE**
   - https://cellxgene.cziscience.com/
   - Search for Xenium datasets

### Converting Data Formats

If you have data in other formats, use the conversion utilities:

```python
import scanpy as sc
import pandas as pd

# From CSV
expr = pd.read_csv('expression.csv', index_col=0)
metadata = pd.read_csv('metadata.csv', index_col=0)
coords = pd.read_csv('coordinates.csv', index_col=0)

# Create AnnData
adata = sc.AnnData(X=expr.values)
adata.obs = metadata
adata.var_names = expr.columns
adata.obs_names = expr.index
adata.obsm['spatial'] = coords[['x', 'y']].values

# Save
adata.write_h5ad('data/raw/sample_01.h5ad')
```

## Storage Considerations

### Disk Space Requirements

Typical file sizes:
- Raw Xenium data: 1-10 GB per sample
- Preprocessed data: 2-20 GB per sample
- Figures: 10-100 MB per sample
- Total per sample: ~5-30 GB

### Backup Strategy

Important data to backup:
1. ✅ Raw data files (always keep originals)
2. ✅ Final processed files (annotated, integrated)
3. ✅ Analysis reports and figures
4. ⚠️ Intermediate files (can be regenerated)

## Data Privacy and Compliance

### For Human Data

- Ensure proper IRB approval
- De-identify patient data
- Follow HIPAA/GDPR guidelines
- Secure data storage

### For Collaborative Projects

- Use consistent naming conventions
- Document metadata thoroughly
- Version control analysis parameters
- Share preprocessing scripts

## Troubleshooting

### Common Issues

**1. File not found errors**
```bash
# Check file exists
ls -lh data/raw/

# Check file permissions
chmod 644 data/raw/your_file.h5ad
```

**2. Format errors**
```python
# Verify AnnData structure
import scanpy as sc
adata = sc.read_h5ad('data/raw/sample.h5ad')
print(adata)  # Should show cells × genes
print(adata.obs.columns)  # Check metadata
print('spatial' in adata.obsm)  # Check coordinates
```

**3. Large file handling**
```python
# For very large files, use backed mode
adata = sc.read_h5ad('data/raw/large_file.h5ad', backed='r')
```

## Contact

For data format questions or issues:
- Open an issue on GitHub
- Check the README.md for general documentation
- See example notebooks for data loading patterns
