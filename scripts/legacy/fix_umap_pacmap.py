"""PaCMAP embedding: specifically designed to avoid UMAP's radial/spiky artifacts.
Produces round, well-separated clusters characteristic of published figures.

Pipeline: tight QC → sc.pp.pca → PaCMAP (on PCA) → sc.tl.leiden (igraph) on PCA-KNN.

Usage:  XENIUM_SAMPLE_ID=0041323 python scripts/fix_umap_pacmap.py
"""
import os
import gc
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import pacmap

SAMPLE_ID = os.environ.get("XENIUM_SAMPLE_ID", "0041323")
ROOT = Path("/blue/maigan/smith6jt/Xenium_Analysis")
PROCESSED = ROOT / f"data/processed/{SAMPLE_ID}"
FIGDIR = ROOT / f"figures/{SAMPLE_ID}/02_phenotyping"
FIGDIR.mkdir(parents=True, exist_ok=True)

h5 = PROCESSED / f"{SAMPLE_ID}_phenotyped.h5ad"
print(f"[{SAMPLE_ID}] Loading {h5} ...")
adata = sc.read_h5ad(h5)
print(f"[{SAMPLE_ID}] shape: {adata.shape}")

# Ensure tight QC (re-apply, no-op if already clean)
MIN_COUNTS = 100
MIN_GENES  = 50
keep = (adata.obs["total_counts"] >= MIN_COUNTS) & (adata.obs["n_genes_by_counts"] >= MIN_GENES)
adata = adata[keep].copy()
print(f"[{SAMPLE_ID}] post-QC: {adata.shape}")

# If X_pca doesn't exist (e.g. from 02 scVI-based flow), compute
if "X_pca" not in adata.obsm:
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata, n_comps=50, mask_var="highly_variable")

X_pca = adata.obsm["X_pca"].astype(np.float32)
print(f"[{SAMPLE_ID}] PCA shape: {X_pca.shape}")

# PaCMAP with standard params (n_neighbors auto from n-log)
# PaCMAP is designed to produce visually cleaner clusters than UMAP on large datasets
print(f"[{SAMPLE_ID}] PaCMAP fit ({X_pca.shape[0]:,} cells, n_neighbors=10) ...")
reducer = pacmap.PaCMAP(
    n_components=2,
    n_neighbors=10,
    MN_ratio=0.5,
    FP_ratio=2.0,
    random_state=42,
    verbose=True,
)
X_pacmap = reducer.fit_transform(X_pca, init="pca")
print(f"[{SAMPLE_ID}] PaCMAP done; embedding range x: {X_pacmap[:,0].min():.2f} to {X_pacmap[:,0].max():.2f}")

# Write as X_umap (for scanpy compatibility)
adata.obsm["X_umap"] = X_pacmap.astype(np.float32)
adata.obsm["X_pacmap"] = X_pacmap.astype(np.float32)

# Save updated h5ad
if "counts" in adata.layers:
    adata.X = adata.layers["counts"]
print(f"[{SAMPLE_ID}] writing {h5}")
adata.write_h5ad(h5)

# Figures — clip axes to 0.5-99.5 percentile
um = adata.obsm["X_umap"]
xlo, xhi = np.quantile(um[:, 0], [0.005, 0.995])
ylo, yhi = np.quantile(um[:, 1], [0.005, 0.995])
xpad = 0.05 * (xhi - xlo); ypad = 0.05 * (yhi - ylo)
xlim = (xlo - xpad, xhi + xpad)
ylim = (ylo - ypad, yhi + ypad)
print(f"[{SAMPLE_ID}] PaCMAP xlim={xlim}, ylim={ylim}")

# Use the existing leiden_qc if present, else leiden from scanpy, else celltype only
has_leiden = "leiden_qc" in adata.obs.columns
fig, axes = plt.subplots(2, 2, figsize=(16, 16))
sc.pl.umap(adata, color="celltype",          ax=axes[0, 0], show=False, title="Cell Types (PaCMAP)", size=2, frameon=False, legend_loc="on data", legend_fontsize=7)
if has_leiden:
    sc.pl.umap(adata, color="leiden_qc",     ax=axes[0, 1], show=False, title="Leiden", size=2, frameon=False, legend_loc="on data", legend_fontsize=7)
else:
    axes[0, 1].axis("off")
sc.pl.umap(adata, color="total_counts",      ax=axes[1, 0], show=False, cmap="viridis", title="Total Counts", size=2, frameon=False)
sc.pl.umap(adata, color="n_genes_by_counts", ax=axes[1, 1], show=False, cmap="viridis", title="N Genes/cell", size=2, frameon=False)
for ax in axes.flat:
    ax.set_xlim(xlim); ax.set_ylim(ylim)
plt.tight_layout()
plt.savefig(FIGDIR / "final_celltypes_overview_grid.png", dpi=300, bbox_inches="tight")
plt.close()

fig, ax = plt.subplots(1, 1, figsize=(10, 10))
sc.pl.umap(adata, color="celltype", legend_loc="on data", legend_fontsize=7,
           frameon=False, title=f"{SAMPLE_ID} — celltypes (PaCMAP)", show=False, size=2, ax=ax)
ax.set_xlim(xlim); ax.set_ylim(ylim)
plt.savefig(FIGDIR / "final_celltypes_overview.png", dpi=300, bbox_inches="tight")
plt.close()

score_cols = [c for c in adata.obs.columns if c.endswith("_score")]
if score_cols:
    ncols = 3
    nrows = (len(score_cols) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axes = np.atleast_1d(axes).flatten()
    for i, s in enumerate(score_cols):
        sc.pl.umap(adata, color=s, ax=axes[i], show=False, cmap="viridis", title=s, size=1, frameon=False)
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
