"""Clean UMAP with tight QC + standard scanpy PCA/UMAP + igraph Leiden + percentile-clipped axes.

Decisions based on empirical testing with this dataset:
  - Tight QC (total_counts>=100, n_genes_by_counts>=50) eliminates low-count spike cells
  - Standard sc.pp.pca + sc.pp.neighbors + sc.tl.umap produces robust embeddings;
    cuml UMAP was producing extreme outliers that collapsed the plot
  - `flavor='igraph', n_iterations=2, directed=False` for Leiden — orders of magnitude
    faster than leidenalg on 1M cells (per scanpy's own recommendation)
  - Figures rendered with 0.5-99.5 percentile axis clip to ignore UMAP outliers

Usage:  XENIUM_SAMPLE_ID=0041323 python scripts/fix_umap_qc.py
"""
import os
import gc
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

SAMPLE_ID = os.environ.get("XENIUM_SAMPLE_ID", "0041323")
ROOT = Path("/blue/maigan/smith6jt/Xenium_Analysis")
PROCESSED = ROOT / f"data/processed/{SAMPLE_ID}"
FIGDIR = ROOT / f"figures/{SAMPLE_ID}/02_phenotyping"
FIGDIR.mkdir(parents=True, exist_ok=True)

h5 = PROCESSED / f"{SAMPLE_ID}_phenotyped.h5ad"
if not h5.exists():
    raise FileNotFoundError(h5)

print(f"[{SAMPLE_ID}] Loading {h5} ...")
adata = sc.read_h5ad(h5)
print(f"[{SAMPLE_ID}] before QC: {adata.shape}")

# 1. Tight QC
MIN_COUNTS = 100
MIN_GENES  = 50
keep = (adata.obs["total_counts"] >= MIN_COUNTS) & (adata.obs["n_genes_by_counts"] >= MIN_GENES)
print(f"[{SAMPLE_ID}] QC drop: {int((~keep).sum()):,} cells")
adata = adata[keep].copy()
print(f"[{SAMPLE_ID}] after QC:  {adata.shape}")

# Clear stale embeddings
for k in list(adata.obsp.keys()): del adata.obsp[k]
for k in ["neighbors", "scvi", "umap", "pca", "_scvi_manager_uuid", "_scvi_uuid"]:
    adata.uns.pop(k, None)
for k in ["X_pca", "X_scvi", "X_umap"]:
    if k in adata.obsm: del adata.obsm[k]
gc.collect()

# 2. Normalize + HVG + scale
if "counts" not in adata.layers:
    raise RuntimeError("no 'counts' layer — cannot re-normalize")
print(f"[{SAMPLE_ID}] normalize + log1p ...")
adata.X = adata.layers["counts"].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers["log_normalized"] = adata.X.copy()

print(f"[{SAMPLE_ID}] HVG ...")
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")

print(f"[{SAMPLE_ID}] scale ...")
sc.pp.scale(adata, max_value=10)

# 3. Standard PCA + neighbors + UMAP
print(f"[{SAMPLE_ID}] sc.pp.pca (n_comps=50, HVG mask) ...")
sc.pp.pca(adata, n_comps=50, mask_var="highly_variable")

print(f"[{SAMPLE_ID}] sc.pp.neighbors (k=30) ...")
sc.pp.neighbors(adata, n_neighbors=30, use_rep="X_pca")

print(f"[{SAMPLE_ID}] sc.tl.umap (min_dist=0.5, spread=1.0) ...")
sc.tl.umap(adata, min_dist=0.5, spread=1.0, random_state=42)

print(f"[{SAMPLE_ID}] sc.tl.leiden (igraph flavor, resolution=0.5) ...")
sc.tl.leiden(
    adata, resolution=0.5, key_added="leiden_qc",
    flavor="igraph", n_iterations=2, directed=False,
)

# 4. Save
print(f"[{SAMPLE_ID}] writing {h5} ...")
if "counts" in adata.layers:
    adata.X = adata.layers["counts"]
adata.write_h5ad(h5)

# 5. Figures — clip to 0.5-99.5 percentile to hide outliers
print(f"[{SAMPLE_ID}] rendering figures (percentile-clipped axes) ...")
um = adata.obsm["X_umap"]
xlo, xhi = np.quantile(um[:, 0], [0.005, 0.995])
ylo, yhi = np.quantile(um[:, 1], [0.005, 0.995])
xpad = 0.05 * (xhi - xlo); ypad = 0.05 * (yhi - ylo)
xlim = (xlo - xpad, xhi + xpad)
ylim = (ylo - ypad, yhi + ypad)
print(f"[{SAMPLE_ID}] UMAP xlim={xlim}, ylim={ylim}")

def _render_grid(suffix):
    fig, axes = plt.subplots(2, 2, figsize=(16, 16))
    sc.pl.umap(adata, color="celltype",          ax=axes[0, 0], show=False, title="Cell Types",     size=2, frameon=False, legend_loc="on data", legend_fontsize=7)
    sc.pl.umap(adata, color="leiden_qc",         ax=axes[0, 1], show=False, title="Leiden (QC+PCA+UMAP)", size=2, frameon=False, legend_loc="on data", legend_fontsize=7)
    sc.pl.umap(adata, color="total_counts",      ax=axes[1, 0], show=False, cmap="viridis", title="Total Counts", size=2, frameon=False)
    sc.pl.umap(adata, color="n_genes_by_counts", ax=axes[1, 1], show=False, cmap="viridis", title="N Genes/cell", size=2, frameon=False)
    for ax in axes.flat:
        ax.set_xlim(xlim); ax.set_ylim(ylim)
    plt.tight_layout()
    plt.savefig(FIGDIR / f"final_celltypes_overview_grid{suffix}.png", dpi=300, bbox_inches="tight")
    plt.close()

_render_grid("")

# Celltype-only
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
sc.pl.umap(adata, color="celltype", legend_loc="on data", legend_fontsize=7,
           frameon=False, title=f"{SAMPLE_ID} — celltypes", show=False, size=2, ax=ax)
ax.set_xlim(xlim); ax.set_ylim(ylim)
plt.savefig(FIGDIR / "final_celltypes_overview.png", dpi=300, bbox_inches="tight")
plt.close()

# Celltype scores
score_cols = [c for c in adata.obs.columns if c.endswith("_score")]
if score_cols:
    ncols = 3
    nrows = (len(score_cols) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axes = np.atleast_1d(axes).flatten()
    for i, s in enumerate(score_cols):
        sc.pl.umap(adata, color=s, ax=axes[i], show=False, cmap="viridis",
                   title=s, size=1, frameon=False)
        axes[i].set_xlim(xlim); axes[i].set_ylim(ylim)
    for j in range(len(score_cols), len(axes)):
        axes[j].axis("off")
    plt.tight_layout()
    plt.savefig(FIGDIR / "celltype_scores_umap.png", dpi=300, bbox_inches="tight")
    plt.close()

if "predicted_celltype" in adata.obs.columns:
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    sc.pl.umap(adata, color="predicted_celltype", ax=axes[0], show=False, title="Predicted (score-based)", size=2, frameon=False)
    sc.pl.umap(adata, color="celltype", ax=axes[1], show=False, title="Final celltype", size=2, frameon=False)
    for ax in axes:
        ax.set_xlim(xlim); ax.set_ylim(ylim)
    plt.tight_layout()
    plt.savefig(FIGDIR / "predicted_celltypes.png", dpi=300, bbox_inches="tight")
    plt.close()

print(f"[{SAMPLE_ID}] Done.")
