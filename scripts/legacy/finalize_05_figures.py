"""Generate 05-integration figures from the saved integrated h5ad.

The 05 notebook runs UMAP on a 500k subsample (CPU UMAP on 2.5M cells is intractable
in this env). The resulting obsm['X_umap'] has NaN for non-subsampled cells, which
makes sc.pl.umap render poorly. This script plots only the subsample rows.

Usage:  python scripts/finalize_05_figures.py
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl

ROOT = Path("/blue/maigan/smith6jt/Xenium_Analysis")
SAMPLE_IDS = os.environ.get("XENIUM_SAMPLE_IDS", "0041323,0041326").split(",")
TAG = "_".join(SAMPLE_IDS)
FIGDIR = ROOT / f"figures/integrated_{TAG}"
FIGDIR.mkdir(parents=True, exist_ok=True)
H5 = ROOT / f"data/processed/integrated/{TAG}_integrated.h5ad"

print(f"Loading {H5} ...")
adata = sc.read_h5ad(H5)
print(f"  shape={adata.shape}")

# Composition figure
comp = pd.crosstab(
    adata.obs["sample_id"], adata.obs["celltype"], normalize="index"
) * 100
print("Composition:\n", comp)
comp.T.plot(kind="bar", figsize=(12, 6))
plt.ylabel("Percentage")
plt.xlabel("Cell Type")
plt.title("Cell Type Composition Across Samples")
plt.legend(title="Sample", bbox_to_anchor=(1.05, 1))
plt.tight_layout()
plt.savefig(FIGDIR / "composition_by_sample.png", dpi=300, bbox_inches="tight")
plt.close()
print("composition_by_sample.png saved")

# Subsample for plotting (use cells that have non-NaN UMAP)
if "X_umap" not in adata.obsm:
    print("No X_umap — skipping UMAP plots")
else:
    um = adata.obsm["X_umap"]
    finite = np.isfinite(um).all(axis=1)
    sub = adata[finite].copy()
    print(f"Plotting {sub.n_obs:,}/{adata.n_obs:,} cells with finite UMAP")

    fig, axes = plt.subplots(2, 2, figsize=(16, 16))
    sc.pl.umap(sub, color="sample_id",          ax=axes[0, 0], show=False, title="Sample",          size=2, frameon=False)
    sc.pl.umap(sub, color="celltype",           ax=axes[0, 1], show=False, title="Cell Type",       size=2, frameon=False)
    if "leiden_integrated" in sub.obs.columns:
        sc.pl.umap(sub, color="leiden_integrated", ax=axes[1, 0], show=False, title="scVI Leiden",    size=2, frameon=False)
    else:
        axes[1, 0].axis("off")
    sc.pl.umap(sub, color="n_genes_by_counts",  ax=axes[1, 1], show=False, cmap="viridis", title="N Genes", size=2, frameon=False)
    plt.tight_layout()
    plt.savefig(FIGDIR / "integrated_umap.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("integrated_umap.png saved")

# Spatial map colored by celltype + niche, split by sample
if "spatial" in adata.obsm:
    for sample in SAMPLE_IDS:
        mask = adata.obs["sample_id"] == sample
        if mask.sum() == 0:
            continue
        ad_s = adata[mask]
        rng = np.random.default_rng(0)
        n = ad_s.n_obs
        idx = rng.choice(n, size=min(200_000, n), replace=False)
        coords = np.asarray(ad_s.obsm["spatial"])[idx]
        celltype_codes = ad_s.obs["celltype"].astype("category").cat.codes.to_numpy()[idx]

        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        cmap = mpl.colormaps.get_cmap("tab20")
        sc_ = ax.scatter(coords[:, 0], coords[:, 1], c=celltype_codes, cmap=cmap, s=0.5, linewidths=0)
        ax.set_title(f"{sample} — celltype (spatial)")
        ax.set_aspect("equal", adjustable="datalim")
        ax.set_xlim(coords[:, 0].min(), coords[:, 0].max())
        ax.set_ylim(coords[:, 1].min(), coords[:, 1].max())
        plt.tight_layout()
        plt.savefig(FIGDIR / f"spatial_celltype_{sample}.png", dpi=300, bbox_inches="tight")
        plt.close()
        print(f"spatial_celltype_{sample}.png saved")

print("Done.")
