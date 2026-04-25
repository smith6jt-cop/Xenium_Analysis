"""Re-compute UMAP + leiden from existing X_scvi latent space on GPU (rapids_singlecell).

Overwrites the spiky CPU UMAP with a GPU UMAP using tuned params:
  n_neighbors=30, min_dist=0.3, spread=1.2. cuGraph leiden replaces the CPU leiden.
For 1.2M-cell Xenium pancreas this runs in ~2-3 min on B200 vs ~60 min on CPU.

Usage:  XENIUM_SAMPLE_ID=0041323 python scripts/redo_umap_from_scvi.py
"""
import os
import gc
from pathlib import Path

os.environ.setdefault("NUMBA_CUDA_ENABLE_PYNVJITLINK", "1")
# rapids-singlecell 0.11 + our numba-cuda pins conflict on pynvjitlink patching
import pynvjitlink.patch as _pv
_pv.patch_numba_linker = lambda *a, **kw: None

import numpy as np
import pandas as pd
import scanpy as sc
import rapids_singlecell as rsc
import matplotlib.pyplot as plt

SAMPLE_ID = os.environ.get("XENIUM_SAMPLE_ID", "0041323")
PROCESSED = Path(f"/blue/maigan/smith6jt/Xenium_Analysis/data/processed/{SAMPLE_ID}")
FIGDIR = Path(f"/blue/maigan/smith6jt/Xenium_Analysis/figures/{SAMPLE_ID}/02_phenotyping")
FIGDIR.mkdir(parents=True, exist_ok=True)

h5 = PROCESSED / f"{SAMPLE_ID}_phenotyped.h5ad"
if not h5.exists():
    raise FileNotFoundError(h5)

print(f"Loading {h5} ...")
adata = sc.read_h5ad(h5)
print(f"  shape={adata.shape}, has X_scvi={'X_scvi' in adata.obsm}")
if "X_scvi" not in adata.obsm:
    raise RuntimeError("no X_scvi in obsm; cannot redo UMAP from scVI latent.")

# Slim heavy state to cut memory footprint before pushing to GPU
for k in ["neighbors", "scvi", "umap"]:
    adata.uns.pop(k, None)
for k in ["scvi_connectivities", "scvi_distances", "connectivities", "distances"]:
    if k in adata.obsp:
        del adata.obsp[k]
if "X_umap" in adata.obsm:
    del adata.obsm["X_umap"]
for lk in ["log_normalized", "scaled"]:
    if lk in adata.layers:
        del adata.layers[lk]
gc.collect()

print("Pushing adata to GPU ...")
rsc.get.anndata_to_GPU(adata)

print("rsc.pp.neighbors(k=50) ...")
rsc.pp.neighbors(adata, use_rep="X_scvi", n_neighbors=50, key_added="scvi")

print("rsc.tl.umap(min_dist=0.5, spread=2.0) ...")
rsc.tl.umap(adata, neighbors_key="scvi", min_dist=0.5, spread=2.0)

print("rsc.tl.leiden(resolution=0.5) ...")
rsc.tl.leiden(adata, resolution=0.5, neighbors_key="scvi", key_added="leiden_scvi")

print("Pulling adata back to CPU ...")
rsc.get.anndata_to_CPU(adata)

# --- Figures (CPU matplotlib on fresh X_umap) ---
print("Rendering figures ...")
sc.pl.umap(adata, color="celltype", legend_loc="on data", legend_fontsize=7,
           frameon=False, title=f"{SAMPLE_ID} celltypes (GPU UMAP)", show=False, size=2)
plt.savefig(FIGDIR / "final_celltypes_overview.png", dpi=300, bbox_inches="tight")
plt.close()

fig, axes = plt.subplots(2, 2, figsize=(16, 16))
sc.pl.umap(adata, color="celltype",          ax=axes[0, 0], show=False, title="Cell Types", size=2, frameon=False)
sc.pl.umap(adata, color="leiden_scvi",       ax=axes[0, 1], show=False, title="scVI Leiden", size=2, frameon=False)
sc.pl.umap(adata, color="total_counts",      ax=axes[1, 0], show=False, cmap="viridis", title="Total Counts", size=2, frameon=False)
sc.pl.umap(adata, color="n_genes_by_counts", ax=axes[1, 1], show=False, cmap="viridis", title="N Genes/cell", size=2, frameon=False)
plt.tight_layout()
plt.savefig(FIGDIR / "final_celltypes_overview_grid.png", dpi=300, bbox_inches="tight")
plt.close()

score_cols = [c for c in adata.obs.columns if c.endswith("_score")]
if score_cols:
    ncols = 3
    nrows = (len(score_cols) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axes = np.atleast_1d(axes).flatten()
    for i, s in enumerate(score_cols):
        sc.pl.umap(adata, color=s, ax=axes[i], show=False, cmap="viridis", title=s, size=1, frameon=False)
    for j in range(len(score_cols), len(axes)):
        axes[j].axis("off")
    plt.tight_layout()
    plt.savefig(FIGDIR / "celltype_scores_umap.png", dpi=300, bbox_inches="tight")
    plt.close()

if "predicted_celltype" in adata.obs.columns:
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    sc.pl.umap(adata, color="predicted_celltype", ax=axes[0], show=False, title="Predicted (score-based)", size=2, frameon=False)
    sc.pl.umap(adata, color="celltype", ax=axes[1], show=False, title="Final celltype", size=2, frameon=False)
    plt.tight_layout()
    plt.savefig(FIGDIR / "predicted_celltypes.png", dpi=300, bbox_inches="tight")
    plt.close()

print(f"Writing updated {h5} ...")
adata.write_h5ad(h5)
print("Done.")
