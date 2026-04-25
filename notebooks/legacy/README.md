# Legacy / deprecated notebooks

These early preprocessing notebooks have been superseded by
`notebooks/01_preprocessing_v2.ipynb`. Kept for git history; do not run.

| Legacy notebook | Why deprecated |
|---|---|
| `01_preprocessing.ipynb` | Original v1: PCA on log-normalized values without HVG mask, no `random_state` on UMAP, leiden via `leidenalg` (slow). |
| `01_preprocessing_panc.ipynb` | Pancreas-specific variant: missed `sc.pp.log1p` between `normalize_total` and PCA, producing arpack failures and NaN propagation. |

The current pipeline does NOT need stage 01 at all when running stage 02 from
zarr (see `02_phenotyping.ipynb` cell 3, the standalone-from-zarr branch).
