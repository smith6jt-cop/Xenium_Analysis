# Legacy / deprecated scripts

These scripts were one-shot fix-ups created during iterative debugging of
the original (pre-2026-04-25) pipeline. They are kept for git history but
are NOT part of the current pipeline. **Do not run them.**

## What replaced them

| Legacy script | Replaced by | Reason |
|---|---|---|
| `fix_umap_pacmap.py` | scVI in `notebooks/02_phenotyping.ipynb` cell 17 | PaCMAP was a workaround for a library-size leak that scVI now corrects via the NB likelihood + `cell_area` covariate. |
| `fix_umap_qc.py` | Inline QC in `notebooks/01_preprocessing_v2.ipynb` cell 14 and `02_phenotyping.ipynb` cell 3 | QC moved into the standard pipeline (`min_counts=50`, `min_genes=20`, 98%-percentile `max_counts`, `min_cells=100`). |
| `redo_umap_from_scvi.py` | `02_phenotyping.ipynb` cell 17 (scanpy CPU UMAP on `X_scvi`) | Rapids/cuML UMAP produced extreme outliers; scanpy CPU UMAP with cosine metric is now the primary embedding. |
| `fix_preprocessing_params.py`, `fix_clustering_params.py`, `refactor_preprocessing_pipeline.py` | Direct edits to `01_preprocessing_v2.ipynb` | Parameter tuning is in-notebook now. |
| `fix_volcano_*.py`, `fix_deg_*.py` | Stage 04 (`04_group_comparisons.ipynb`) | DE plotting now lives in the analysis notebooks. |
| `finalize_03_niche.py`, `finalize_05_figures.py` | Inline plot rendering in `03_spatial_analysis.ipynb` and `05_tissue_comparisons.ipynb` | One-shot finalize scripts merged into the notebooks. |
| `edit_notebook_04_subsample.py` | Inline subsampling in `04_group_comparisons.ipynb` | One-shot notebook editor no longer needed. |
