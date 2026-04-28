# Data layout

Project-specific reference. Generic Xenium data-loading docs live in the [10x
Xenium documentation](https://www.10xgenomics.com/support/in-situ-gene-expression).

## What this project actually has

Two human pancreas FFPE samples on the **hAtlas v1.1 + 100-gene custom T1D
add-on** panel (Xenium 5K base + custom probes). Both samples are FFPE,
probe-based, chemistry v2.

| Sample | Cells (post-QC) | Donor sex | Notes |
|---|---:|---|---|
| 0041323 | ~1,145,295 | female | (KDM5D ~0%) |
| 0041326 | ~1,168,473 | male | (KDM5D ~10%, more advanced T1D progression by IAPP loss) |

`donor_sex` is auto-derived from KDM5D + DDX3Y detection — see
`scripts/finalize_panel.py`.

## Directory structure

```
data/
├── raw/
│   ├── 0041323.zarr/        # spatialdata zarr ingested via notebook 00
│   └── 0041326.zarr/
├── processed/
│   ├── 0041323/
│   │   ├── 0041323_phenotyped.h5ad         # PRIMARY output of stage 02
│   │   ├── 0041323_spatial_analysis.h5ad   # PRIMARY output of stage 03
│   │   ├── 0041323_immune_phenotyped.h5ad  # gated immune subset (notebook 07)
│   │   └── *.csv                           # per-stage tables
│   ├── 0041326/                            # same layout
│   ├── panel_audit.csv                     # candidate panel exclusions
│   ├── islet_*.csv                         # insulitis headline tables
│   ├── immune_proximity_summary.csv
│   └── legacy/                             # archived pre-rebuild artifacts
└── xenium_explorer/
    ├── {sample}/analysis.zarr.zip          # drop next to experiment.xenium
    ├── {sample}_cell_groups.csv            # Cell → Add cell categorization
    └── {sample}_color_palette.json
```

## h5ad layer / obs / obsm contracts

`*_phenotyped.h5ad` (stage 02 output) MUST contain:

- `layers['counts']` — raw integer counts
- `layers['lognorm']` — `log1p(normalize_total)` (sample-comparable scale)
- `obsm['spatial']` — (n_cells, 2) xy in microns
- `obsm['X_scvi']` — scVI latent (30 dims)
- `obsm['X_umap']` — scVI-derived UMAP
- `obs['celltype']`, `obs['celltype_lineage']`, `obs['leiden_scvi']`
- `obs['donor_sex']` — `'male'` / `'female'` (added by `scripts/finalize_panel.py`)
- `var['panel_for_scoring']`, `var['panel_for_embedding']` — panel-audit boolean masks

`*_spatial_analysis.h5ad` (stage 03 output) adds:
- `obs['spatial_niche']` — KMeans cluster ID
- `uns['nhood_enrichment']`, `uns['co_occurrence']`, `uns['ligrec']` — squidpy outputs

## Headline CSV outputs

| File | Owner script | Schema |
|---|---|---|
| `panel_audit.csv` | `scripts/audit_panel_exclusions.py` | gene-by-gene flagging w/ exclude category |
| `islet_infiltration_per100endo.csv` | `scripts/insulitis_analysis.py` | per-islet rows; per-phenotype zone counts, rates, grades; composition |
| `islet_insulitis_grades.csv` | `scripts/insulitis_analysis.py` | long format: per-(sample × subtype) grade prevalence |
| `islet_insulitis_thresholds.csv` | `scripts/insulitis_analysis.py` | per-(sample × size_class × subtype) rotation-null thresholds |
| `islet_insulitis_regression.csv` | `scripts/insulitis_analysis.py` | per-(sample × subtype) odds ratios, no p-values (n=2) |
| `immune_proximity_summary.csv` | `scripts/insulitis_analysis.py` | density enrichment per immune subtype × sample |
| `rare_celltype_verification.csv` | `scripts/verify_rare_celltypes.py` | per-(sample × type) dotplot/coexpression checks |

`scripts/insulitis_analysis.py` is the **single writer** of the islet/immune
headline CSVs. Do not regenerate them from `scripts/rerun_immune_pipeline.py` or
`scripts/build_notebook_07.py` — those have early-exit redirects pointing here.

## Xenium Explorer round-trip

`utils/xenium_explorer_export.py::export_for_xenium_explorer(adata,
group_keys=[...], output_dir=..., sample_name=...)` writes both formats and a
color palette JSON. The CLI wrapper is
`scripts/export_xenium_explorer_groups.py`.

## Storage footprint per sample

- raw zarr: ~5–10 GB
- `_phenotyped.h5ad`: ~60 GB
- `_spatial_analysis.h5ad`: ~10 GB
- `_immune_phenotyped.h5ad`: ~1 GB

## Data privacy

Human pancreas tissue from de-identified donors (consent on file). Standard
HIPAA/IRB practices apply; no PHI in any committed file.
