# Pickup handoff — 2026-04-25 09:55 EDT

This document is the single source of truth for resuming work on the
Xenium analysis pipeline after disconnect. **Read this top-to-bottom
before doing anything else.**

## Current state

### What is on disk and ready to use

| Sample | Phenotyped h5ad (stage 02) | Spatial h5ad (stage 03) | Validated UMAP |
|---|---|---|---|
| 0041323 | `data/processed/0041323/0041323_phenotyped.h5ad` (59 GB, 2026-04-24 20:15) | `data/processed/0041323/0041323_spatial_analysis.h5ad` (10.4 GB, 2026-04-25 09:32) | ✓ scVI, ρ_leak={0.067, 0.185} |
| 0041326 | `data/processed/0041326/0041326_phenotyped.h5ad` (60 GB, 2026-04-24 22:15) | `data/processed/0041326/0041326_spatial_analysis.h5ad` (10.5 GB, 2026-04-25 09:50) | ✓ scVI, ρ_leak={0.125, 0.140} |

Both samples have full sets of figures, ligrec/Moran/niche/nhood CSVs.
Library-size correction quantitatively validated (Spearman correlation of
each UMAP axis vs `total_counts` < 0.3 for both samples).

### What changed in the notebooks (not yet executed)

`notebooks/02_phenotyping.ipynb` was patched after the runs above
finished. The patches are STAGED IN THE NOTEBOOK but the saved h5ads
above were produced from the PREVIOUS version. To get the improved
clustering the notebook now describes, **rerun stage 02**.

The three patches:

1. **Cell 17 — leiden resolution 0.5 → 1.0** (both `leiden_scvi` and
   `leiden_pca`). Was producing ~10 clusters that merged biologically
   distinct populations (T_cells+Myeloid landed in one cluster,
   Ductal+Stellate+Endothelial in another). At 1.0 we expect ~15-20
   clusters that better separate the 14 marker categories.

2. **Cell 15 — anti-acinar override** added after the initial
   `idxmax`-based assignment. The Xenium panel has only `AMY1A` for
   Acinar, so any cell with non-zero `AMY1A` can argmax-out as Acinar
   even if it's strongly expressing specific markers of another lineage.
   The override pass:
   - Computes per-cell mean expression of curated specific anti-acinar
     marker sets in `layers["lognorm"]`:
     - Beta: `IAPP, MAFA, NKX6-1, PDX1, SLC30A8`
     - Endocrine: `CHGA, INSM1, ISL1, NEUROD1, FEV`
     - Ductal: `KRT19, CFTR, MUC1`
     - Stellate: `ACTA2, PDGFRB, RGS5`
     - Endothelial: `PECAM1, CDH5, CLDN5, PLVAP`
     - T_cells: `CD3D, CD3E, CD8A, PTPRC`
     - B_cells: `MS4A1, CD79A`
     - Myeloid: `CD68, CD163, MARCO`
     - Schwann: `MPZ, SOX10`
   - For each Acinar-called cell, finds the alt type with the largest
     positive `(alt_marker_mean - alt_marker_threshold)` margin and
     reassigns. Threshold = median of that mean across cells already
     called the alt type (or 75th-percentile fallback if too few).
   - Per-type breakdown of reassignments printed to the cell output.

3. **Cell 17 — `predicted_celltypes.png` rendering moved to the end**
   so it lands on the scVI UMAP. Previously it was in cell 15 which
   ran BEFORE scVI created `X_umap`, so the figure was on the
   stale Xenium-builtin UMAP. Same problem as `celltype_scores_umap.png`
   that we already fixed.

`notebooks/03_spatial_analysis.ipynb` cells 6 and 18 were also patched
during yesterday's runs and the saved spatial-analysis h5ads above
already reflect those patches:
- Cell 6: `nhood_enrichment(n_jobs=1, n_perms=100, seed=0)` (full
  defaults hung; `n_jobs=-1` hit a numba TypeError on readonly arrays).
- Cell 18: defensive `try/except` around each CSV export so a single
  failure doesn't kill the rest.

## How to resume — exactly

### 1. Submit the GPU node

```bash
sbatch /home/smith6jt/vscode_gpu.sh
```

`vscode_gpu.sh` is already configured with: B200 GPU, 64 CPUs, 512 GB
RAM, 72 h walltime. That's plenty for stage 02+03 × 2 samples (we used
~80 GB peak and the previous run finished in <6 h).

### 2. Connect VS Code and open a terminal in the tunnel

`cd /blue/maigan/smith6jt/Xenium_Analysis`

### 3. Re-run stage 02 only (most important)

The patched 02 notebook still ingests directly from zarr — stage 01
is not required.

```bash
STAGES="02" nohup bash scripts/run_local_pipeline.sh \
    > logs/pipeline_driver_$(date +%Y%m%d_%H%M%S).log 2>&1 &
```

Per-sample timing from the previous run:
- 0041323: 2 h 10 min (18:05 → 20:15)
- 0041326: 2 h (20:15 → 22:15)

Total stage 02 × 2 samples ≈ **4 h 10 min**.

### 4. (Optional) Re-run stage 03 if the new clusters change spatial niches

Stage 03's niche analysis (`spatial_niche` in cell 12) is keyed on
`celltype` (the manual_celltype mapping), which is downstream of leiden
resolution and the anti-acinar override. After stage 02 completes, the
existing stage 03 outputs are stale.

```bash
STAGES="03" nohup bash scripts/run_local_pipeline.sh \
    > logs/pipeline_driver_$(date +%Y%m%d_%H%M%S).log 2>&1 &
```

Stage 03 × 2 samples ≈ **30 min** (it's mostly subsampled).

Or run both stages back-to-back:
```bash
STAGES="02 03" nohup bash scripts/run_local_pipeline.sh \
    > logs/pipeline_driver_$(date +%Y%m%d_%H%M%S).log 2>&1 &
```

### 5. Monitor

```bash
# tail the driver log for stage transitions
tail -F $(ls -t logs/pipeline_driver_*.log | head -1)

# OR launch the py-spy frame monitor (requires conda env active)
source /blue/maigan/smith6jt/miniforge3/etc/profile.d/conda.sh
conda activate xenium_analysis
# (optional one-time install) pip install py-spy

# Find the live nbconvert kernel:
KPID=$(ps -eo pid,etime,cmd --sort=-etime \
    | awk '/ipykernel_launcher.*scratch.*tmp.*\.json/ && !/awk/ {print $1}' \
    | head -1)
py-spy dump --pid $KPID
```

`nbconvert` buffers all cell stdout/stderr into the executed-notebook
output, so the `.err` log only shows the nbconvert header until the
run ends. Use `py-spy` (or check file mtimes on `figures/` and
`data/processed/`) to see live progress.

## Verification once stage 02 finishes

```bash
source /blue/maigan/smith6jt/miniforge3/etc/profile.d/conda.sh
conda activate xenium_analysis

python3 << 'EOF'
import scanpy as sc
import numpy as np
from scipy.stats import spearmanr
for s in ("0041323", "0041326"):
    a = sc.read_h5ad(f"data/processed/{s}/{s}_phenotyped.h5ad")
    um = a.obsm.get("X_umap_scvi", a.obsm["X_umap"])
    rx, _ = spearmanr(um[:, 0], a.obs["total_counts"])
    ry, _ = spearmanr(um[:, 1], a.obs["total_counts"])
    print(f"{s}: shape={a.shape}  leiden_scvi clusters={a.obs['leiden_scvi'].nunique()}")
    print(f"      Spearman(UMAP vs total_counts) = ({rx:.3f}, {ry:.3f})  (must be < 0.3)")
    print(f"      celltype distribution:")
    print(a.obs["celltype"].value_counts().to_string())
    print()
EOF
```

Expected: `leiden_scvi clusters` ≈ 15-20 (was 10-11 at resolution 0.5);
Spearman ρ both axes < 0.3; `celltype` distribution shows distinct
non-Acinar populations now that anti-acinar override is active.

Then visually check:
- `figures/{sample}/02_phenotyping/umap_qc_overlay.png` — total_counts /
  n_genes / cell_area panels should be uniform color (no gradient).
- `figures/{sample}/02_phenotyping/predicted_celltypes.png` — should
  render on the scVI UMAP with cleaner type boundaries than before.
- `figures/{sample}/02_phenotyping/dotplot_markers_by_cluster.png` —
  more clusters, each more cleanly enriched in one cell-type marker
  block.

## What to do if something fails

| Symptom | What it likely is | Fix |
|---|---|---|
| Spearman ρ > 0.3 | scVI didn't converge or under-trained | Bump `max_epochs=400` in cell 17, add `early_stopping_patience=40` |
| Single huge "Acinar" cluster after override | `lognorm` layer missing or override threshold too high | Check `available_markers` printout in cell 10 — confirm which markers are actually present |
| Memory OOM during scVI | Batch size too large | Drop `batch_size=4096` to 2048 in cell 17 |
| `numba.TypeError: readonly array` from squidpy | Numba+joblib readonly-array bug | Already worked around by `n_jobs=1` in cell 6 of 03 |
| nbconvert hangs mid-leiden | Leiden on 1.2M cells is slow | Wait — typically 30-50 min per leiden call. Use `py-spy` to confirm forward motion. |

If the pipeline truly stalls, kill it:
```bash
pkill -f "^/blue.*jupyter-nbconvert.*--execute"
pkill -f "ipykernel_launcher.*scratch"
pkill -f "^bash scripts/run_local_pipeline"
```

## Repo structure (post-cleanup)

```
notebooks/
├── 00_ingest_xenium.ipynb       # zarr ingest (rarely re-run)
├── 01_preprocessing_v2.ipynb    # OPTIONAL — stage 02 ingests from zarr directly
├── 02_phenotyping.ipynb         # PRIMARY — has anti-acinar override + leiden 1.0
├── 03_spatial_analysis.ipynb    # squidpy spatial stats (already-patched)
├── 04_group_comparisons.ipynb   # group-level DE
├── 05_tissue_comparisons.ipynb  # cross-tissue integration (NOT used in pancreas-only flow)
├── 06_xenium_phenocycler_integration.ipynb
└── legacy/                      # deprecated 01 variants — DO NOT USE
    ├── 01_preprocessing.ipynb
    └── 01_preprocessing_panc.ipynb

scripts/
├── run_local_pipeline.sh         # primary driver
├── slim_phenotyped.py            # post-process to shrink h5ad sizes
├── export_xenium_explorer_groups.py
├── slurm/                        # SLURM job scripts (unchanged)
└── legacy/                       # deprecated one-shot fix-ups — DO NOT USE
    ├── README.md                 # explains what each replaced
    ├── fix_umap_pacmap.py
    ├── fix_umap_qc.py
    ├── redo_umap_from_scvi.py
    └── ... (10 more)

data/
├── raw/{0041323,0041326}.zarr   # Xenium spatialdata zarr stores
└── processed/{0041323,0041326}/ # h5ads + per-sample CSVs

figures/{0041323,0041326}/{01_preprocessing,02_phenotyping,03_spatial_analysis}/
```

## Key architectural decisions (do not regress)

1. **Stage 02 ingests from zarr directly** (no stage 01 dependency).
   Implemented in cell 3 of `02_phenotyping.ipynb`. Stage 01's
   `_preprocessed.h5ad` (~100 GB) is not needed.

2. **scVI covariate is `cell_area` only.** Do NOT pass `total_counts`
   as a continuous covariate — the NB likelihood already models library
   size as a size factor; passing it again double-regresses and was
   the original library-size leak that PaCMAP was a workaround for.

3. **No HVG masking on PCA / scVI input.** The Xenium panel is curated
   (~5 k probes, all biologically chosen). HVG selection is computed
   for metadata only (`adata.var["highly_variable"]`) but PCA uses all
   panel genes.

4. **Squidpy graph stats are run on subsamples** (150k for
   co-occurrence, 200k for Moran's I, 100k for ligrec & Ripley).
   Centrality_scores is currently **skipped entirely** because
   `nx.closeness_centrality` on >100k nodes is days-to-weeks even
   subsampled.

5. **scanpy CPU UMAP, not rapids/cuML.** The cuML UMAP produced
   extreme outliers; scanpy CPU UMAP with `metric='cosine'`,
   `min_dist=0.5`, pinned `random_state=0` is reproducible and clean.

6. **`leiden` flavor is `igraph`, `n_iterations=-1`** — much faster
   than leidenalg on 1 M cells.

## Files modified today (2026-04-25)

- `notebooks/02_phenotyping.ipynb` — cells 15 + 17 patched (anti-acinar
  override, leiden res 1.0, plot deferred to scVI UMAP)
- `notebooks/03_spatial_analysis.ipynb` — cells 6 + 18 patched yesterday
  (already in saved spatial h5ads)
- `scripts/legacy/` and `notebooks/legacy/` — moved 13 stale scripts +
  2 legacy notebooks, with `README.md` in each explaining what
  replaced them
- `data/processed/legacy/` — **moved ~200 GB of pre-pivot artifacts
  out of `data/processed/{sample}/`**:
  - The broken pre-2026-04-24 `*_preprocessed.h5ad` (~100 GB each).
    These were made by the broken pipeline and would have been
    auto-loaded by cell 3 of `02_phenotyping.ipynb`, defeating the
    rerun. **Do NOT move them back.**
  - Stale `*_phenotyped_slim.h5ad` (~2.7 GB each, predate the current
    phenotyped.h5ad).
  - Stage-05 integration outputs (~10 GB) — explicitly out of scope for
    the per-sample workflow.
- `figures/legacy/` — moved stale `integrated_*` figure dirs and an
  orphan `notebooks/figures/` from past stage-05 attempts.
- `IMPLEMENTATION_SUMMARY.md`, `README.md`, `CLAUDE.md` — updated to
  reflect the post-cleanup structure
- `figures/{sample}/02_phenotyping/` — 4 stale figures regenerated
  on the scVI UMAP / leiden_scvi grouping; `final_celltypes_overview_grid.png`
  (PaCMAP leftover) deleted

After all moves, `data/processed/{sample}/` is the clean post-rerun
landing zone — only the current `_phenotyped.h5ad`, `_spatial_analysis.h5ad`,
and per-sample CSVs remain. Stage 02 cell 3 will see no `_preprocessed.h5ad`
and ingest directly from the raw zarr.
