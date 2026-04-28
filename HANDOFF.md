# Pickup handoff — 2026-04-25 (post-GPU-patch)

This document is the single source of truth for resuming work on the
Xenium analysis pipeline after disconnect. **Read this top-to-bottom
before doing anything else.**

## Current state

### What is on disk and ready to use

| Sample | Phenotyped h5ad (stage 02) | Spatial h5ad (stage 03) | Validated UMAP |
|---|---|---|---|
| 0041323 | `data/processed/0041323/0041323_phenotyped.h5ad` (59 GB, 2026-04-24 20:15) | `data/processed/0041323/0041323_spatial_analysis.h5ad` (10.4 GB, 2026-04-25 09:32) | ✓ scVI, ρ_leak={0.067, 0.185} |
| 0041326 | `data/processed/0041326/0041326_phenotyped.h5ad` (60 GB, 2026-04-24 22:15) | `data/processed/0041326/0041326_spatial_analysis.h5ad` (10.5 GB, 2026-04-25 09:50) | ✓ scVI, ρ_leak={0.125, 0.140} |

These are the **pre-GPU-patch CPU runs**. They're correct (library-size
leak gate passes) but stage 02 was at leiden res=0.5 and missing the
anti-acinar override. The current rerun is `STAGES="01 02 03"` from
scratch with GPU paths.

### What changed (and is now staged in the notebooks)

**1. Driver-script env fix (`scripts/run_local_pipeline.sh`)** — adds
`LD_PRELOAD` (conda `libstdc++.so.6` + `/apps/compilers/cuda/12.8.1/lib64/libcudart.so.12`)
and extends `LD_LIBRARY_PATH` to include `$CUDA_HOME/nvvm/lib64`. This
fixes the `cudf 24.12 → pynvjitlink → numba_cuda` MVC-patch conflict
that was making every `import rapids_singlecell` fail in the previous
runs. Without this fix, the GPU patches below silently no-op (notebooks
fall back to CPU). Driver vs runtime libcudart now both report 12.8 →
cudf skips the patch entirely.

**2. `rapids_singlecell` upgraded `0.11.1 → 0.14.1`** (single-package
upgrade) — needed for `_handle_mask_var` API compat with scanpy 1.12.

**3. `notebooks/02_phenotyping.ipynb`** (cells 1, 13, 15, 17):
- Cell 1: rsc/cupy/torch imports + B200 TF32 setup
  (`set_float32_matmul_precision('high')`, `allow_tf32=True`).
- Cell 17: both `sc.pp.neighbors` calls (scVI and PCA paths) replaced
  with `rsc.pp.neighbors` (cuML kNN, GPU). UMAP+leiden stay CPU.
- Cell 17: leiden resolution 0.5 → 1.0 (was producing ~10 mixed clusters).
- Cell 15: anti-acinar override (Beta IAPP/MAFA/NKX6-1/PDX1/SLC30A8;
  Ductal KRT19/CFTR/MUC1; etc. — see CLAUDE.md for full set).
- Cell 17: `predicted_celltypes.png` rendering moved to end (so it lands
  on scVI UMAP, not stale Xenium-builtin UMAP).

**4. `notebooks/01_preprocessing_v2.ipynb`** (cells 1, 18):
- Cell 1: rsc/cupy/torch imports + B200 TF32 setup.
- Cell 18: `sc.pp.scale → rsc.pp.scale`, `sc.pp.pca → rsc.pp.pca`,
  `sc.pp.neighbors → rsc.pp.neighbors`. Wraps GPU block with
  `rsc.get.anndata_to_GPU/CPU`. UMAP and leiden remain CPU.

**5b. QC upper-filter switched to density** (stage 01 cell 14 +
stage 02 cell 3 standalone branch): `counts > q98` → `counts/cell_area > q98`.
Empirical motivation: Spearman(total_counts, cell_area) ≈ 0.78 in this
dataset; 88% of cells dropped by absolute counts have normal density.
Density filter preserves large healthy acinar/ductal cells.

**5c. Lineage phenotyping** (stage 02 cell 22 + post-hoc
`scripts/lineage_phenotype.py`). New obs columns:
`lineage_status` ∈ {single_lineage, lineage_phenotyped, lineage_ambiguous,
doublet_suspected}; `celltype_lineage` (refined celltype, copy of
`celltype` unless reclassified); `lineage_dominant`, `lineage_top_depth`,
`lineage_second_depth`, `lineage_n_groups`. Catches AMY1A-overspill
mislabels (e.g. macrophages-with-engulfed-acinar argmax-scored to Acinar
get correctly resolved to Myeloid) by counting how many lineage-specific
markers each candidate lineage has at raw count ≥ 3, and reclassifying
the cell when one lineage's marker depth exceeds the runner-up by ≥ 2.
Acinar canary panel expanded from `[AMY1A]` to `[AMY1A, CUZD1]` so depth
comparisons are fair. **The original `celltype` column is preserved**;
downstream consumers can opt into `celltype_lineage`.

**5. `notebooks/03_spatial_analysis.ipynb`** (cells 1, 6, 10, 12, 18):
- Cell 1: rapids/cuML imports added.
- Cell 10: `sq.gr.spatial_autocorr → rsc.gr.spatial_autocorr` (GPU
  Moran's I). May allow running on full 1.2M cells; 200k subsample is
  retained as a graceful fallback if GPU OOMs.
- Cell 12: `MiniBatchKMeans → cuml.cluster.KMeans`,
  `sklearn NearestNeighbors → cuml.neighbors.NearestNeighbors` for the
  K=50 spatial graph + niche clustering. **Biggest stage-03 win.**
- Cell 6 (`nhood_enrichment(n_jobs=1, n_perms=100, seed=0)`) and
  cell 18 (defensive CSV exports) — unchanged from previous runs;
  preserved here because they're documented squidpy/numba workarounds.

Backups of pre-GPU notebooks are at `notebooks/01_preprocessing_v2.ipynb.bak.preplanB`,
`notebooks/02_phenotyping.ipynb.bak.preplanB`, and
`notebooks/03_spatial_analysis.ipynb.bak.preplanB`.

## How to resume — exactly

### 1. Submit the GPU node

```bash
sbatch /home/smith6jt/vscode_gpu.sh
```

`vscode_gpu.sh`: B200 GPU, 64 CPUs, 512 GB RAM, 72 h walltime. The new
GPU pipeline uses ~10–25 GB GPU memory peak (B200 has 192 GB) and
~150–200 GB RAM peak.

### 2. Connect VS Code and open a terminal in the tunnel

`cd /blue/maigan/smith6jt/Xenium_Analysis`

### 3. Run the full pipeline from scratch (recommended)

```bash
STAGES="01 02 03" nohup bash scripts/run_local_pipeline.sh \
    > logs/pipeline_driver_$(date +%Y%m%d_%H%M%S).log 2>&1 &
```

**Per-sample timing on the GPU pipeline (estimated):**
- Stage 01 preprocessing: ~30–40 min (ingest + HVG + PCA-on-GPU +
  neighbors-on-GPU + CPU UMAP + CPU leiden + squidpy decorative).
- Stage 02 phenotyping: ~60–75 min (cell-3 ingest + score_genes +
  anti-acinar override + scVI 200 epochs on B200 with TF32 +
  rsc.pp.neighbors twice + CPU UMAP twice + CPU leiden twice +
  figure rendering).
- Stage 03 spatial: ~15–20 min (squidpy spatial_neighbors +
  nhood_enrichment + co-occurrence + GPU Moran's I + GPU niche
  KMeans + ligrec + ripley).

**Total for both samples: ~3.5 h** (vs ~6 h on the previous CPU-only
runs of just 02+03).

### 4. Stage-only re-runs (after a notebook edit)

```bash
# Just phenotyping
STAGES="02" nohup bash scripts/run_local_pipeline.sh \
    > logs/pipeline_driver_$(date +%Y%m%d_%H%M%S).log 2>&1 &

# Just spatial (requires fresh phenotyped.h5ad upstream)
STAGES="03" nohup bash scripts/run_local_pipeline.sh \
    > logs/pipeline_driver_$(date +%Y%m%d_%H%M%S).log 2>&1 &
```

### 5. Verify GPU paths actually fired
The first thing to check after a run finishes — silent CPU fallback
is the most likely failure mode:
```bash
grep -E "TF32 enabled|rsc.pp|cuML KMeans|cuNearestNeighbors|GPU:" \
    notebooks/executed/0041323/02_0041323.ipynb | head -10
```
Expected output includes lines like `GPU: NVIDIA B200 (TF32 enabled, ...)`,
`rsc.pp.neighbors: 14.79s`, `cuML KMeans...`. If those are missing,
the LD_PRELOAD env fix didn't apply — re-source driver and verify with:
```bash
LD_PRELOAD=$CONDA_PREFIX/lib/libstdc++.so.6:/apps/compilers/cuda/12.8.1/lib64/libcudart.so.12 \
python -c "import rapids_singlecell as rsc; print(rsc.__version__)"
# should print: 0.14.1
```

### 6. Monitor

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

## Verification once `STAGES="01 02 03"` finishes

```bash
source /blue/maigan/smith6jt/miniforge3/etc/profile.d/conda.sh
conda activate xenium_analysis

python3 << 'EOF'
import scanpy as sc
import numpy as np
from scipy.stats import spearmanr
for s in ("0041323", "0041326"):
    # 1) preprocessed contract
    p = sc.read_h5ad(f"data/processed/{s}/{s}_preprocessed.h5ad")
    assert "leiden_1.5" in p.obs.columns, f"{s}: leiden_1.5 missing"
    assert "counts" in p.layers and "lognorm" in p.layers
    assert "X_pca" in p.obsm and "X_umap" in p.obsm
    print(f"{s} preprocessed: shape={p.shape}, leiden_1.5={p.obs['leiden_1.5'].nunique()} clusters")

    # 2) phenotyped: library-size leak gate
    a = sc.read_h5ad(f"data/processed/{s}/{s}_phenotyped.h5ad")
    um = a.obsm.get("X_umap_scvi", a.obsm["X_umap"])
    rx, _ = spearmanr(um[:, 0], a.obs["total_counts"])
    ry, _ = spearmanr(um[:, 1], a.obs["total_counts"])
    print(f"  phenotyped: leiden_scvi={a.obs['leiden_scvi'].nunique()} clusters")
    print(f"              Spearman(UMAP vs total_counts) = ({rx:.3f}, {ry:.3f})")
    assert max(abs(rx), abs(ry)) < 0.3, f"{s}: library-size leak"
    print("  celltype:")
    print(a.obs["celltype"].value_counts().to_string())

    # 3) spatial: niches + Moran I
    sp = sc.read_h5ad(f"data/processed/{s}/{s}_spatial_analysis.h5ad")
    assert "spatial_niche" in sp.obs.columns, f"{s}: spatial_niche missing"
    print(f"  spatial: niches={sp.obs['spatial_niche'].nunique()}")
    print()
print("ALL CHECKS PASSED")
EOF
```

Expected: `leiden_1.5` ≈ 24 clusters (preprocessed); `leiden_scvi`
≈ 15-20 at res=1.0; Spearman ρ both UMAP axes < 0.3; `celltype`
distribution shows distinct non-Acinar populations after anti-acinar
override; `spatial_niche` populated with ~12 niches.

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
| `RuntimeError: Cannot patch Numba: numba_cuda includes patches from pynvjitlink` at notebook start | LD_PRELOAD didn't apply (driver script wasn't sourced) | Re-run via `bash scripts/run_local_pipeline.sh` (don't bypass the driver). Verify `echo $LD_PRELOAD` includes `libcudart.so.12` and `libstdc++.so.6`. |
| Notebook output shows only `sc.pp.*` lines (no `rsc.*` / `cuML KMeans` / `TF32 enabled`) | rapids imports silently failed; CPU fallback active | Same as above — driver-script env wasn't picked up. |
| Spearman ρ > 0.3 | scVI didn't converge or under-trained | Bump `max_epochs=400` in cell 17, add `early_stopping_patience=40` |
| Single huge "Acinar" cluster after override | `lognorm` layer missing or override threshold too high | Check `available_markers` printout in cell 10 — confirm which markers are actually present |
| Memory OOM during scVI | Batch size too large | Drop `batch_size=4096` to 2048 in cell 17 |
| `cuML KMeans` GPU OOM | n_init too high or comp matrix too large | Drop `n_init=10` to 5 |
| `rsc.gr.spatial_autocorr` GPU OOM on full 1.2M | full Moran's I needs ~30 GB GPU mem | Re-add the 200k subsample (uncomment the `adata_mi = adata[sub_idx].copy()` block) |
| `numba.TypeError: readonly array` from squidpy | Numba+joblib readonly-array bug | Already worked around by `n_jobs=1` in cell 6 of 03 |
| nbconvert hangs mid-leiden | Leiden on 1.2M cells is slow | Wait — typically 5-15 min per leiden call (CPU igraph). Use `py-spy` to confirm forward motion. |

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

## 2026-04-27 — Stage 07 immune+islet proximity (added) + dual-seed islet ID

### What was added

1. **`notebooks/07_immune_islet_proximity.ipynb`** — pulls all immune
   cells + the indeterminate-with-immune-marker subset from both
   phenotyped h5ads, retrains scVI immune-only with `batch_key='sample'`
   to remove sample-of-origin variation, scores immune subtypes
   (T_helper/cytotoxic/reg/exhausted, NK, Macro_M1/M2/resident,
   Monocyte, DC, B_pan, B_plasma), identifies islets via DBSCAN on
   endocrine seed coordinates, computes per-cell distance-to-islet
   via cuML NN, bins (intra <25 μm, peri 25–50 μm, proximal 50–200 μm,
   distal >200 μm), and produces contingency χ² + per-islet
   infiltration scores normalized per 100 endocrine cells.

2. **Sample-comparable thresholds (lognorm).** Replaced all
   `raw_count >= N` thresholds with `lognorm >= log1p(N)` because
   0041326 has ~10% lower median total_counts and 298 panel genes
   drop by min_cells=100 (vs. 14 in 0041323). A static raw threshold
   was systematically more stringent in 0041326. Touched files:
   `scripts/lineage_phenotype.py`, `scripts/verify_upper_tail_doublets.py`,
   `scripts/show_canary_panels.py` (DETECTION RULE text), notebook 02
   cell 22 (lineage block), notebook 07 cells 5 (immune-marker
   recovery for indeterminate cells) + 15 (islet seed).

3. **Dual islet seeding.** Two parallel seed strategies:
   - **Strict hormone-marker seed** (`scripts`-internal, in notebook 07
     cell 15): `lognorm >= log1p(5)` of `[CHGA, INSM1, ISL1, NEUROD1,
     FEV]`. Output: `data/processed/islets_dbscan.csv` (110 islets in
     0041323, **2 islets in 0041326** — biological asymmetry, see below).
   - **Endo-label seed** (`scripts/islet_proximity_endolabel.py`):
     `celltype_lineage in {Endocrine, Beta, Alpha, Delta, Endocrine_pan}`
     (Epsilon excluded — 150k count in 0041326 is artefactual).
     Output: `data/processed/islets_dbscan_endolabel.csv` (265 islets in
     0041323, 318 islets in 0041326 — comparable across samples).

### The 0041326 islet asymmetry — what it means

With the **strict hormone-marker seed**, 0041326 yields only 2 islets
because hormone markers (CHGA/INSM1/ISL1/NEUROD1/FEV) are genuinely
under-expressed across this section. Per-100-endocrine numbers from the
strict CSV for 0041326 are **statistically meaningless** (every immune
cell in the section is NN-assigned to one of those 2 seeds).

With the **endo-label seed**, 0041326 recovers 318 islets — comparable
to 0041323's 265. The endo-label seed uses score-genes-argmax cell-type
assignments that don't require any single canonical marker to be highly
expressed, so a cell with merely above-cluster-average islet-program
expression still seeds. **Both samples now have islet-rate parity**
(7.30 vs 9.93 islets per 1000 endocrine-labeled cells).

### Headline biological findings (use endo-label CSVs for cross-sample comparison)

- **0041323 (less T1D progression?):** Most islets show NO immune
  infiltration; the small fraction that do are dominated by T_cytotoxic
  and T_reg. Median per-100-endocrine for every immune subtype = 0.
  Only DC shows a meaningful islet-zone enrichment (OR=2.58).

- **0041326 (more T1D progression?):** Pervasive low-level immune
  infiltration. Median per-100-endocrine: Monocyte=16, T_reg=12,
  Macro_M2=9, T_helper=4, B_plasma=3. Top-5 most-infiltrated islets
  carry 150–223 immune cells each, dominated by Macro_M2 (17–36),
  T_reg (22–36), and T_cytotoxic (8–18). Mild OR enrichments:
  B_plasma=1.36, T_cytotoxic=1.30, NK=1.15.

This asymmetry is consistent with 0041326 being further along the T1D
disease course — broader immune infiltration corresponds to greater
β-cell hormone marker depletion, which is exactly why the strict-seed
approach failed for that sample.

### CSVs produced

| File | Purpose |
|---|---|
| `data/processed/islets_dbscan.csv` | strict hormone-marker DBSCAN islets (use for high-stringency 0041323-only analyses) |
| `data/processed/islets_dbscan_endolabel.csv` | endo-label DBSCAN islets (use for cross-sample comparisons) |
| `data/processed/immune_proximity_summary.csv` | strict-seed contingency (subtype × bin × sample) |
| `data/processed/immune_proximity_summary_endolabel.csv` | endo-label contingency |
| `data/processed/islet_infiltration_per100endo.csv` | strict-seed per-islet immune counts + per100endo rates |
| `data/processed/islet_infiltration_per100endo_endolabel.csv` | endo-label per-islet immune counts + per100endo rates |

### Helper scripts

- `scripts/islet_proximity_endolabel.py` — generates the endo-label CSVs.
- `scripts/compare_islet_seeds.py` — side-by-side strict vs endo-label
  numerical summary (run after the above).
- `scripts/summarize_07.py` — dumps the strict-seed notebook 07
  summary (immune subtype distribution, per-islet stats, OR table for
  0041323).

### 2026-04-27 update — replaced strict-seed proximity tables with insulitis grades

The strict hormone-marker seed gave 0041326 only 2 islets, which made
its `immune_proximity_summary.csv` row uninterpretable (every immune
cell got binned to 'distal' by construction). User flagged the table
as "absolute garbage" and was right.

Replaced the headline tables with `scripts/insulitis_analysis.py` which:
- Uses the endo-label seed (`celltype_lineage` ∈ Endocrine/Beta/Alpha/Delta/
  Endocrine_pan) as the islet basis. Islet counts: 0041323=265, 0041326=318.
- Assigns each islet an `insulitis_grade` per Campbell-Thompson 2013
  nPOD criteria adapted to total CD45+ pool: `no_insulitis` (<3 immune
  within 50 μm), `peri_insulitis` (3-5), `insulitis` (≥6).
- Computes per-subtype `density_enrichment` = (zone fraction) /
  (distal fraction). >1 means concentrated near islets; <1 depleted.
  Sample-comparable, doesn't depend on per-islet sparsity.

Headline files now:
| File | Content |
|---|---|
| `data/processed/immune_proximity_summary.csv` | density enrichment per subtype × sample |
| `data/processed/islet_infiltration_per100endo.csv` | per-islet immune counts + insulitis_grade |
| `data/processed/islet_insulitis_grades.csv` | per-sample insulitis prevalence |
| `data/processed/islets_dbscan.csv` | endo-label seed islets (was strict-seed) |
| `*_strict.csv` | original strict-seed tables, archived |

Findings (**SUPERSEDED — see 2026-04-28 entry below**): 0041326 = 99.4%
insulitis (316/318 islets), 0041323 = 16.6% insulitis. Those numbers came
from an earlier immune-gating version (mean 42 immune cells/islet) that
no longer matches the current `*_immune_phenotyped.h5ad` (mean ~10
cells/islet). On the current immune set, the CT-3/6 absolute grade is
60.8% / 66.7% — much closer between samples — and the per-phenotype
rotation-null grades (the new headline) are 0–4% by subtype.

### 2026-04-28 update — per-phenotype dynamic insulitis grading

The single CT-3/6 absolute grade above flattens 0041326 to "everything is
insulitis". User asked for a size-, composition-, and immune-type-aware
replacement. `scripts/insulitis_analysis.py` rewritten to:

1. **Composition per islet** from per-cell `Beta_score` / `Alpha_score` /
   `Delta_score` (already present in `_phenotyped.h5ad` from stage 02 cell 13).
   Argmax with margin ≥ 0.1 else `Endo_unresolved`. `composition_class` ∈
   {Beta_rich, Alpha_rich, Delta_rich, Mixed} per 60% threshold.
2. **Per-phenotype rotation null** (1000 rotations of immune xy around the
   clustered-seed centroid) for each of 12 immune subtypes
   (`B_pan`, `B_plasma`, `DC`, `Macro_M1`, `Macro_M2`, `Macro_resident`,
   `Monocyte`, `NK`, `T_cytotoxic`, `T_exhausted`, `T_helper`, `T_reg`).
   Acinar/Schwann/Endothelial leakage dropped from immune count.
3. **Per-stratum thresholds**: `(sample × size_class)` with 5 size bins
   (`tiny ≤20`, `small 21-50`, `med 51-100`, `large 101-200`, `xl >200`),
   pooled per-100-endocrine null values, 95th / 99th percentile →
   `thresh_peri` / `thresh_insulitis`.
4. **Per-phenotype grade** per islet — independent {no, peri, insulitis}
   grade for each subtype, plus the CT-3/6 absolute grade preserved on
   `total_immune`. **Singular `insulitis_grade` column dropped.**
5. **Descriptive logistic regression** per (sample × subtype) with odds
   ratios + 95% CI, NO p-values (n=2 caveat).

#### New / changed CSVs (single writer = `insulitis_grade_absolute.py`)

| File | Schema |
|---|---|
| `islet_infiltration_per100endo.csv` | per-islet, 53 cols incl. `<S>_zone`, `<S>_per100endo`, `<S>_grade` for each of 12 subtypes, plus `composition_class`, `size_class`, `total_immune`, `total_per100endo`, `insulitis_grade_absolute` |
| `islet_insulitis_grades.csv` | long format: `sample, subtype, no_insulitis, peri_insulitis, insulitis, total_islets, pct_insulitis, pct_peri_or_insulitis` — 26 rows (2 samples × 13 = 12 subtypes + `absolute_total`) |
| `islet_insulitis_thresholds.csv` | NEW — `sample, size_class, subtype, n_islets, n_null_samples, thresh_peri, thresh_insulitis` (120 rows) |
| `islet_insulitis_regression.csv` | NEW — `sample, subtype, term, odds_ratio, ci_low, ci_high` (no p-values) |
| `*_preweighted.csv` | one-shot archive of pre-rewrite headline CSVs |

#### Findings (current run, 2026-04-28)

| Subtype | 0041323 % insulitis | 0041326 % insulitis |
|---|---:|---:|
| absolute_total (CT 3/6) | 60.8% (161/265) | 66.7% (212/318) |
| NK | 2.3% | 3.8% |
| Monocyte | 2.3% | 3.1% |
| T_helper | 2.3% | 2.2% |
| T_cytotoxic | 0.4% | 2.2% |
| DC | 1.5% | 1.6% |
| (others) | 0–1% | 0–1% |

The "99.4% / 16.6%" headlines from the 2026-04-27 entry above are **stale**
— they were computed from an earlier immune-gating version (mean 42
immune cells/islet) that no longer matches the current
`*_immune_phenotyped.h5ad` (mean ~10 cells/islet). The CT-3/6 absolute
grade on the current immune set is 60.8% / 66.7% — much closer between
samples than the prior headline suggested.

Per-phenotype prevalences are intentionally conservative — the
rotation-null 99th-percentile threshold flags only islets that are
statistically elevated *for that subtype × size_class stratum* relative
to the same sample's immune-density background. To loosen, drop the
threshold percentile in `insulitis_analysis.py` (currently 95/99 for
peri/insulitis; try 90/95 for a more permissive cut).

#### Composition findings

- 0041323: 64 Beta_rich islets (24%), 200 Mixed (76%), 1 Alpha_rich.
- 0041326: 6 Beta_rich (2%), 5 Alpha_rich (2%), 307 Mixed (96%) — β-cell
  identity loss consistent with more-advanced T1D progression.

#### Rotation-null caveat

Rotation around the clustered-seed centroid means rotated immune cells
can land outside the original tissue convex hull. Those positions
contribute zero to per-islet zone counts, biasing the null *low* and
making real-vs-null gaps look slightly larger than reality. A
tissue-mask-aware permutation would fix this; v1 does not implement it.

#### Duplicate-writer reconciliation

`scripts/rerun_immune_pipeline.py` step 6 now early-exits with a
redirect message; `scripts/build_notebook_07.py` writes its derivatives
to `*_legacy_notebook07.csv` so headline CSVs aren't silently
overwritten with the old schema.

#### Known follow-up

- Some per-(sample × subtype) regressions hit `LinAlgError: Singular
  matrix` because of perfect separation in small composition strata
  (e.g. 0041326 has only 5 Alpha_rich islets). Skipped silently; not a
  blocker.
- Per-phenotype thresholds may be too strict if the user's question is
  "where is *any* T-cell signal elevated"; loosen percentile or fold
  T_cytotoxic + T_exhausted into a single `T_effector` aggregate column
  if that question becomes important.

### 2026-04-28 (later) — Epsilon mislabel + full reanalysis pending

#### What's broken in the current saved h5ads

`scripts/audit_panel_exclusions.py` + manual review surfaced two structural
issues with the current `_phenotyped.h5ad` files:

1. **150,230 cells in 0041326 mislabeled `Epsilon`.** Cause: the stage 02
   cell 10 marker dict had `'Epsilon': ['GHRL']` — a 1-gene "panel".
   Under detection asymmetry (0041326 has lower per-cell detection),
   the score_genes argmax overweights single-gene panels relative to
   multi-gene panels. Result: ~13% of all cells in 0041326 got
   argmax-mislabeled Epsilon. 0041323 has 0 Epsilon (correct).
2. **Donor sex differs between samples.** KDM5D = 10.3% in 0041326 vs
   0.010% in 0041323; DDX3Y = 7.7% vs 0.017%. **0041326 is male,
   0041323 is female.** Y-chromosome variance is leaking into PCA, scVI,
   UMAP, leiden, and likely contaminating cell-type calls beyond
   Epsilon.

#### Resolution path — full reanalysis with the locked plan at `/home/smith6jt/.claude/plans/full-reanalysis-plan.md`

User-confirmed parameters:
- Margin gate: **0.2** (strict)
- Panel coverage threshold: **30%**
- Spatial coherence: **disabled for rare cell types** (only common
  clustering types get the check)
- scVI sex covariate: **YES**
- Cross-sample prevalence ratio: **reported only, not a gate** (n=2)
- Rare cell types (Epsilon, Gamma): identified via **rule-based
  multi-criteria definitions** (not score_genes argmax) — see
  `scripts/refine_rare_celltypes.py`

Helper scripts written 2026-04-28 (ready to use):
- `scripts/audit_panel_exclusions.py` — produces `data/processed/panel_audit.csv` (106 candidate exclusions)
- `scripts/finalize_panel.py` — adds `var['panel_for_scoring']`, `var['panel_for_embedding']`, `obs['donor_sex']`
- `scripts/refine_rare_celltypes.py` — post-stage-02 margin gate + panel coverage + rule-based Gamma/Epsilon + spatial coherence
- `scripts/verify_rare_celltypes.py` — per-(sample × type) verification suite

Notebook edits **still pending** (must be done before pipeline kickoff):
- `notebooks/02_phenotyping.ipynb` cell 10: drop `'Epsilon': ['GHRL']`; mask
  marker panels via `var['panel_for_scoring']`
- `notebooks/02_phenotyping.ipynb` cell 17: add `donor_sex` to scVI
  `categorical_covariate_keys`; restrict gene set to `panel_for_embedding`
- `notebooks/01_preprocessing_v2.ipynb` cell 18: HVG / PCA on
  `panel_for_embedding` only

Wall-time after kickoff: ~3.5–4 h on the B200 GPU node for stages 01 + 02 + 03,
followed by `refine_rare_celltypes.py` (~5 min), `insulitis_analysis.py`
(~30 s), `figs_07_immune.py` (~5 min), and `verify_rare_celltypes.py` (~10 min).
