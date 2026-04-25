# Repository guide for Claude

## Read first
Before doing ANYTHING in this repo, read [`HANDOFF.md`](HANDOFF.md) — it
documents the current pipeline state, what was just changed, and how to
resume work. The information below is the durable architectural guide;
HANDOFF.md is the per-session pickup note.

## What this repo does

End-to-end analysis of 10x Xenium 5k-panel spatial transcriptomics on
two human pancreas FFPE samples:

| Sample | Cells | Tissue | Slide |
|---|---:|---|---|
| 0041323 | 1,218,859 | pancreas (T1D add-on panel) | hAtlas_v1.1+100 custom |
| 0041326 | 1,254,240 | pancreas (T1D add-on panel) | hAtlas_v1.1+100 custom |

Both panels are FFPE / probe-based / chemistry v2. The custom add-on is
T1D-focused — **canonical hormones (INS, GCG, SST, PPY, TTR, IRX2) and
most acinar zymogens (PRSS1/2, CPA1/2, CTRB1/2, CELA1/2A/3A, REG1A/3A)
are NOT on the panel**. Phenotyping uses surrogate markers — see the
panel-aware marker dict in `notebooks/02_phenotyping.ipynb` cell 10.

## Pipeline (current, post-2026-04-25)

```
zarr (data/raw/{sample}.zarr)
  │
  ├─► [optional] notebooks/01_preprocessing_v2.ipynb
  │       → data/processed/{sample}/{sample}_preprocessed.h5ad
  │       (NOT required — stage 02 ingests from zarr directly)
  │
  └─► notebooks/02_phenotyping.ipynb  ← PRIMARY ENTRY POINT
        Cell 3:    ingest zarr (or load preprocessed if present),
                   drop control probes, QC filter (min_counts=50,
                   min_genes=20, max_counts=98%-ile, min_cells=100),
                   counts layer preserved.
        Cell 5:    normalize_total + log1p + lognorm layer.
                   Marker discovery only if leiden_1.5 already exists
                   (Xenium-builtin clustering from zarr).
        Cell 10:   panel-aware pancreas marker dict (NO INS/GCG/SST/PPY).
                   LOW_CONFIDENCE_TYPES = {Alpha, Delta, Epsilon, Acinar}.
        Cell 13:   sc.tl.score_genes for each cell type.
                   Plot deferred to cell 17 (so it lands on scVI UMAP).
        Cell 15:   argmax assignment → ANTI-ACINAR OVERRIDE
                   (per-cell mean of curated specific anti-acinar markers
                   in lognorm; reassign Acinar→alt if margin > 0)
                   → Indeterminate flagging for low-coverage types
                   (median - 1·MAD on celltype_confidence).
        Cell 17:   scVI training (n_latent=30, gene_likelihood='nb',
                   continuous_covariate_keys=['cell_area'], NO total_counts)
                   → sc.pp.neighbors(use_rep='X_scvi', metric='cosine')
                   → sc.tl.umap(min_dist=0.5, random_state=0)
                   → sc.tl.leiden(resolution=1.0, flavor='igraph',
                                  n_iterations=-1, random_state=0)
                   → umap_qc_overlay.png (validation: Spearman ρ < 0.3
                     between each UMAP axis and total_counts)
                   → PCA-fallback embedding (X_umap_pca, leiden_pca)
                   → celltype_scores_umap.png + predicted_celltypes.png
                     rendered HERE on scVI UMAP.
        Cell 20:   per-cluster consensus annotation
                   (manual_celltype = leiden_scvi → majority predicted_celltype).
        Cell 22:   final celltypes overview, save to {sample}_phenotyped.h5ad.
        │
        └─► notebooks/03_spatial_analysis.ipynb
              Cell 2:    load + slim adata (drop X_pca, X_scvi, X_umap_pca,
                         scvi_connectivities, etc. — keeps spatial workflow
                         RAM bounded).
              Cell 4:    sq.gr.spatial_neighbors(n_neighs=10, delaunay) on full.
              Cell 6:    sq.gr.nhood_enrichment(n_perms=100, n_jobs=1)
                         (n_jobs>1 hits a numba/joblib readonly-array bug)
              Cell 8:    sq.gr.co_occurrence on 150k subsample.
              Cell 10:   Moran's I on 200k subsample, n_jobs=16, top 100
                         genes by mean expression.
              Cell 12:   spatial niches via composition kmeans + 2-pass
                         spatial smoothing (Banksy-style).
              Cell 14:   sq.gr.ligrec on 100k subsample, n_perms=50.
              Cell 16:   sq.gr.ripley on 100k subsample.
              Cell 18:   defensive try/except CSV exports + h5ad save.
```

Stages 04, 05, 06 (group comparisons / cross-tissue / Phenocycler) are
present in `notebooks/` but not part of the per-sample default flow.

## Critical architectural decisions — do not regress

1. **scVI covariate = `cell_area` ONLY.** Do not pass `total_counts` as
   a continuous covariate — NB likelihood already models library size
   via a size factor.
2. **No HVG masking on PCA / scVI input.** The 5k panel is already
   curated. HVG is metadata-only.
3. **Squidpy graph stats are subsampled.** Full 1.2M-cell ops with
   permutation tests or all-pairs Dijkstra are infeasible.
   `centrality_scores` skipped entirely.
4. **scanpy CPU UMAP, not rapids/cuML.** Rapids UMAP produced extreme
   outliers on this data.
5. **`leiden` flavor=`igraph`, `n_iterations=-1`** — orders of
   magnitude faster than leidenalg on 1 M cells.
6. **Pin random_state=0 everywhere** — PCA, neighbors, UMAP, leiden,
   scVI training. Reproducibility is non-negotiable.

## Repo layout (current)

```
notebooks/
├── 00_ingest_xenium.ipynb       # Xenium output dir → spatialdata zarr
├── 01_preprocessing_v2.ipynb    # OPTIONAL — produces _preprocessed.h5ad
├── 02_phenotyping.ipynb         # PRIMARY — ingests zarr, scVI, score, assign
├── 03_spatial_analysis.ipynb    # squidpy spatial stats
├── 04_group_comparisons.ipynb   # group-level DE
├── 05_tissue_comparisons.ipynb  # cross-tissue integration (NOT used here)
├── 06_xenium_phenocycler_integration.ipynb
├── *.bak                        # pre-pivot backups (kept for git safety)
└── legacy/                      # deprecated v1 + panc preprocessing variants

scripts/
├── run_local_pipeline.sh         # primary driver — STAGES env var selects steps
├── slim_phenotyped.py            # post-process to shrink h5ad sizes
├── export_xenium_explorer_groups.py  # CLI for utils/xenium_explorer_export
├── slurm/                        # SLURM job scripts
└── legacy/                       # deprecated one-shot fix-ups (see README.md inside)

utils/
├── analysis_utils.py             # scanpy/squidpy helpers
└── xenium_explorer_export.py     # CSV + zarr export to Xenium Explorer

data/
├── raw/{sample}.zarr             # Xenium spatialdata zarr stores
└── processed/{sample}/           # h5ads + per-sample CSVs

figures/{sample}/{stage}/         # all rendered figures
logs/                             # per-stage .err, .out, plus pipeline_driver_*.log
notebooks/executed/{sample}/      # nbconvert-executed notebook copies (debug)
```

## Common operations

### Re-run stage 02 only (most common after a notebook edit)
```bash
cd /blue/maigan/smith6jt/Xenium_Analysis
STAGES="02" nohup bash scripts/run_local_pipeline.sh \
    > logs/pipeline_driver_$(date +%Y%m%d_%H%M%S).log 2>&1 &
```

### Run a single stage on a single sample
```bash
bash scripts/run_local_pipeline.sh 02 0041323
```

### Re-render figures from saved h5ad without re-running anything
```bash
source /blue/maigan/smith6jt/miniforge3/etc/profile.d/conda.sh
conda activate xenium_analysis
python3 /tmp/regen_phenotyping_figs.py 0041323  # if the script still exists
```

### Validate a phenotyped h5ad
```python
import scanpy as sc
from scipy.stats import spearmanr
a = sc.read_h5ad("data/processed/0041326/0041326_phenotyped.h5ad")
um = a.obsm.get("X_umap_scvi", a.obsm["X_umap"])
print(spearmanr(um[:, 0], a.obs["total_counts"]))  # |ρ| should be < 0.3
print(a.obs["leiden_scvi"].nunique())              # ≈ 15-20 at resolution 1.0
print(a.obs["celltype"].value_counts())            # cell-type breakdown
```

## Background on the original failure mode

The pre-pivot pipeline (April 22-24) had every UMAP dominated by a
library-size gradient — `total_counts` varied monotonically along the
axes. Root causes that were fixed:

1. `MIN_COUNTS=5 / MAX_COUNTS=1500` (cell 14 of old preprocessing) —
   too loose at the bottom (single-transcript spike cells survived) and
   truncating the upper biological tail.
2. HVG selection AFTER PCA — dead code; PCA used all genes.
3. scVI without `cell_area` covariate — segmentation-area effects bled
   into the latent.
4. `svd_solver='arpack'` on 1.2 M cells (initial fix attempt) — Lanczos
   inner loop is single-threaded; one run wedged for 50+ min before
   producing anything. Switched to `svd_solver='randomized'` with
   pinned `random_state=0`.
5. Rapids/cuML UMAP produced extreme outliers — moved to scanpy CPU.
6. Default `n_jobs=1` on Moran's I and centrality — single-threaded
   permutation tests on 1.2 M cells were infeasible. Centrality
   skipped, Moran's I subsampled to 200 k.

Quantitative validation: Spearman correlation of each UMAP axis vs
`total_counts` is now {0.067, 0.185} for 0041323 and {0.125, 0.140}
for 0041326 — both well under the 0.3 threshold for "library size is
NOT driving the embedding."

## When in doubt

- HANDOFF.md → "what just happened, how to resume"
- This file (CLAUDE.md) → "how the repo is structured"
- README.md → user-facing docs
- IMPLEMENTATION_SUMMARY.md → high-level inventory (CI, Xenium Explorer
  export, etc.)
- DATA_README.md → expected data formats and Xenium Explorer outputs
