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

## Pipeline (current, GPU-accelerated post-2026-04-25)

GPU annotations: **[GPU rsc]** uses `rapids_singlecell` (cuML kNN /
PCA / scale / spatial_autocorr); **[GPU cuML]** uses `cuml.cluster` /
`cuml.neighbors` directly; **[GPU torch]** uses scVI on B200 with TF32
Tensor Cores enabled in cell 1; **[CPU]** is intentionally CPU.

```
zarr (data/raw/{sample}.zarr)
  │
  ├─► notebooks/01_preprocessing_v2.ipynb  ← OPTIONAL (stage 02 can
  │       ingest zarr directly), but rerun for clean fresh outputs.
  │       Cell 1:    imports + B200 TF32 setup [GPU torch].
  │       Cell 18:   counts layer → HVG (seurat_v3, metadata only) [CPU]
  │                  → normalize_total + log1p + lognorm layer [CPU, fast]
  │                  → anndata_to_GPU + rsc.pp.scale [GPU rsc]
  │                  → rsc.pp.pca(n_comps=30) [GPU rsc]
  │                  → anndata_to_CPU(convert_all=True)
  │                  → rsc.pp.neighbors(k=15, cosine) [GPU rsc]
  │                  → sc.tl.umap(min_dist=0.5) [CPU, intentional]
  │                  → sc.tl.leiden(resolution=1.5, flavor='igraph') [CPU]
  │       → data/processed/{sample}/{sample}_preprocessed.h5ad
  │
  └─► notebooks/02_phenotyping.ipynb  ← PRIMARY
        Cell 1:    imports + B200 TF32 setup [GPU torch].
        Cell 3:    ingest zarr (or load preprocessed if present),
                   drop control probes, QC filter (min_counts=50,
                   min_genes=20, max_counts=98%-ile, min_cells=100),
                   counts layer preserved.
        Cell 5:    normalize_total + log1p + lognorm layer [CPU, fast].
                   Marker discovery only if leiden_1.5 already exists.
        Cell 10:   panel-aware pancreas marker dict (NO INS/GCG/SST/PPY).
                   LOW_CONFIDENCE_TYPES = {Alpha, Delta, Epsilon, Acinar}.
        Cell 13:   sc.tl.score_genes for each cell type [CPU, fast].
                   Plot deferred to cell 17 (so it lands on scVI UMAP).
        Cell 15:   argmax assignment → ANTI-ACINAR OVERRIDE
                   (per-cell mean of curated specific anti-acinar markers
                   in lognorm; reassign Acinar→alt if margin > 0)
                   → Indeterminate flagging for low-coverage types
                   (median - 1·MAD on celltype_confidence).
        Cell 17:   scVI training [GPU torch, TF32]
                   (n_latent=30, gene_likelihood='nb',
                    continuous_covariate_keys=['cell_area'], NO total_counts)
                   → rsc.pp.neighbors(use_rep='X_scvi', cosine) [GPU rsc]
                   → sc.tl.umap(min_dist=0.5, random_state=0) [CPU]
                   → sc.tl.leiden(resolution=1.0, flavor='igraph') [CPU]
                   → umap_qc_overlay.png (Spearman ρ < 0.3 gate)
                   → rsc.pp.neighbors(use_rep='X_pca') [GPU rsc] +
                     PCA-fallback UMAP (X_umap_pca, leiden_pca, both CPU).
                   → celltype_scores_umap.png + predicted_celltypes.png
                     rendered HERE on scVI UMAP.
        Cell 20:   per-cluster consensus annotation
                   (manual_celltype = leiden_scvi → majority predicted_celltype).
        Cell 22:   final celltypes overview, save to {sample}_phenotyped.h5ad.
        │
        └─► notebooks/03_spatial_analysis.ipynb
              Cell 1:    imports including cuKMeans, cuNearestNeighbors.
              Cell 2:    load + slim adata (drop X_pca, X_scvi, etc.).
              Cell 4:    sq.gr.spatial_neighbors(n_neighs=10, delaunay) [CPU, fast].
              Cell 6:    sq.gr.nhood_enrichment(n_perms=100, n_jobs=1) [CPU,
                         documented numba/joblib readonly-array workaround].
              Cell 8:    sq.gr.co_occurrence on 150k subsample [CPU].
              Cell 10:   rsc.gr.spatial_autocorr (Moran's I) [GPU rsc] —
                         attempts full 1.2M cells; 200k subsample retained
                         as fallback if OOM.
              Cell 12:   spatial niches: cuML NearestNeighbors (k=50, brute) [GPU cuML]
                         + cuML KMeans (n_clusters=12, scalable-k-means++) [GPU cuML]
                         + 2-pass numpy majority-vote smoothing [CPU, fast].
              Cell 14:   sq.gr.ligrec on 100k subsample, n_perms=50 [CPU].
              Cell 16:   sq.gr.ripley on 100k subsample [CPU].
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
   `centrality_scores` skipped entirely. (Moran's I via
   `rsc.gr.spatial_autocorr` runs full-1.2M on GPU now — see #7.)
4. **scanpy CPU UMAP, not rapids/cuML.** Rapids UMAP produced extreme
   outliers on this data. **Other ops (scale, PCA, neighbors, KMeans,
   Moran's I, score_genes) ARE on GPU via `rapids_singlecell` — only
   UMAP and leiden stay CPU.**
5. **`leiden` flavor=`igraph`, `n_iterations=-1`** — fast enough on
   CPU; we don't switch to cugraph leiden because `dask-cuda 24.12`
   has a hard incompat with the installed `dask 2026.1`.
6. **Pin random_state=0 everywhere** — PCA, neighbors, UMAP, leiden,
   scVI training. Reproducibility is non-negotiable.
7. **Driver-script env fix is required for GPU paths.**
   `scripts/run_local_pipeline.sh` exports `LD_PRELOAD` (conda
   `libstdc++.so.6` + `/apps/compilers/cuda/12.8.1/lib64/libcudart.so.12`)
   and adds `$CUDA_HOME/nvvm/lib64` to `LD_LIBRARY_PATH`. Without these,
   `cudf` 24.12 calls `pynvjitlink.patch_numba_linker()` and raises
   `RuntimeError` because `numba_cuda` 0.0.17.1 already provides those
   patches — every `import rapids_singlecell` fails and the notebooks
   silently fall back to CPU. The LD_PRELOAD pins runtime libcudart at
   12.8 (matching driver 12.8), so the MVC patch is skipped entirely.
8. **B200 Tensor Cores must be enabled for scVI.** Cell 1 of stages
   01/02/03 sets `torch.set_float32_matmul_precision('high')` and
   `torch.backends.cuda.matmul.allow_tf32 = True`. Without this scVI
   training takes ~30-50% longer on B200 (compute capability 10.0).
9. **`rapids_singlecell == 0.14.1`** — required for scanpy 1.12 API
   compat (the older 0.11.1 calls `_handle_mask_var` with the wrong
   signature for scanpy 1.12+).
10. **Upper QC filter is DENSITY-based (counts/cell_area), not absolute
    counts.** In FFPE Xenium, `cell_area` varies 5-10x and Spearman
    (total_counts, cell_area) ≈ 0.78 in this dataset, so an absolute-counts
    upper cap predominantly removes large *healthy* cells (acinar /
    ductal / large epithelial) — 88% of cells dropped by `counts > q98`
    have completely normal counts/area density. We use `counts/area > q98`
    to drop genuine outlier-density cells (likely doublets) instead.
    Lower bounds remain absolute (`min_counts=50`, `min_genes=20`)
    because those are sanity floors for noise/segmentation artifacts,
    not size-correlated. Implemented in stage 01 cell 14 and stage 02
    cell 3 (standalone branch).
11. **Lineage phenotyping (post-density-filter doublet refinement).**
    Implemented in stage 02 cell 22 (and as a post-hoc script
    `scripts/lineage_phenotype.py`). For each cell, count how many
    distinct markers from each canary lineage panel are detected at
    `lognorm ≥ log1p(3)` ("depth", sample-comparable). When ≥ 2
    mutually-exclusive groups are co-detected, the cell is resolved to
    the dominant lineage if `top_depth − second_depth ≥ 2 AND top_depth ≥ 2`;
    tagged `lineage_phenotyped` and `celltype_lineage` is set to the
    dominant lineage. ≥ 3 mutually-exclusive groups co-detected → tagged
    `doublet_suspected` (segmentation merge), kept as-is.
    `celltype` is NEVER overwritten by lineage refinement —
    `celltype_lineage` is the refined column for downstream consumers.
    Acinar canary panel was expanded from `[AMY1A]` to `[AMY1A, CUZD1]`
    so depth comparisons are fair across all lineages. The exact panels
    live in `scripts/verify_upper_tail_doublets.py::CANARY_PANELS`.

12. **Cell-type scoring panels MUST have ≥ 3 genes.** Single-gene
    "panels" (the historical `'Epsilon': ['GHRL']` entry in stage 02
    cell 10 was the only offender) are forbidden. Reason: the
    score_genes argmax under detection asymmetry biases toward the
    least-comprehensive panel — 0041326's lower per-cell detection let
    the 1-gene Epsilon score beat multi-gene panels, mislabeling 150K
    cells (~13% of the section) as Epsilon. Rare cell types with
    insufficient panel coverage are identified via **rule-based
    multi-criteria definitions** in `scripts/refine_rare_celltypes.py`,
    NOT via score_genes argmax. Currently rule-based: Gamma cells
    (AQP3 + ETV1 + pan-endocrine, NOT Beta/Delta) and Epsilon cells
    (GHRL + pan-endocrine, NOT Beta/Alpha/Delta/Gamma).

13. **Donor-sex confound.** Samples differ by donor sex: 0041323 is
    female, 0041326 is male. Y-chromosome genes (KDM5D, DDX3Y et al.)
    are at ~10% / ~7.7% in 0041326 vs ~0% in 0041323. The full
    sex-chromosome gene set is excluded from cell-type scoring AND
    from PCA/scVI input via the `panel_for_embedding` var-mask set by
    `scripts/finalize_panel.py`. scVI training adds `donor_sex` to
    `categorical_covariate_keys`. **Cross-sample comparison must be
    sex-aware**; do not pool male+female cells without controlling
    for the chromosome.

14. **Panel audit lives at `data/processed/panel_audit.csv`.** 106
    candidate exclusions across 12 categories (sperm/germline,
    photoreceptor, hepatocyte, cardiac/skeletal muscle, sex chromosome,
    bone, etc.) — genes biologically implausible in pancreas/LN/vessel
    tissue. Two boolean var-masks: `panel_for_scoring` (drops 100
    non-sex tissue-restricted genes) and `panel_for_embedding`
    (additionally drops 6 sex-chromosome genes). Set by
    `scripts/finalize_panel.py`.

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

### Run the whole pipeline (stages 01 → 02 → 03) from scratch
```bash
cd /blue/maigan/smith6jt/Xenium_Analysis
STAGES="01 02 03" nohup bash scripts/run_local_pipeline.sh \
    > logs/pipeline_driver_$(date +%Y%m%d_%H%M%S).log 2>&1 &
```
Total wall time on B200 GPU + 64 CPU node: ≈ 3.5 h for both samples
(stage 01 ≈ 30–40 min/sample, stage 02 ≈ 60–75 min/sample, stage 03
≈ 15–20 min/sample). Stage 02 alone on CPU was ~2 h 10 min/sample
pre-GPU; stage 03 was ~30 min/sample.

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

### Verifying GPU paths actually fired
```bash
# After a run finishes:
grep -E "TF32 enabled|rsc.pp|cuML KMeans|cuNearestNeighbors|GPU:" \
    notebooks/executed/0041323/02_0041323.ipynb | head -10
```
If you see only `sc.pp.*` and no `rsc.*` / `cuML` / `TF32` lines, the
LD_PRELOAD env fix in the driver script wasn't picked up and rapids
imports silently failed (notebook fell back to CPU). Re-source the
driver, then verify:
```bash
LD_PRELOAD=$CONDA_PREFIX/lib/libstdc++.so.6:/apps/compilers/cuda/12.8.1/lib64/libcudart.so.12 \
python -c "import rapids_singlecell as rsc; print(rsc.__version__)"
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
- DATA_README.md → expected data formats, h5ad contracts, headline CSVs
