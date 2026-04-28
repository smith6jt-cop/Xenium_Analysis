"""Re-render umap_qc_overlay.png and final_celltypes_overview.png with
p95-clipped color scales and larger legend fonts.
"""
import time
from pathlib import Path

import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt

ROOT = Path("/blue/maigan/smith6jt/Xenium_Analysis")
SAMPLES = ("0041323", "0041326")


def render(sample):
    p = ROOT / f"data/processed/{sample}/{sample}_phenotyped.h5ad"
    fig_dir = ROOT / f"figures/{sample}/02_phenotyping"
    print(f"\n  --- {sample} ---")
    t = time.time()
    adata = sc.read_h5ad(p)
    print(f"    loaded {adata.shape} in {time.time()-t:.1f}s")

    # ---- umap_qc_overlay.png ----
    print(f"    rendering umap_qc_overlay.png ...")
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    for col, ax, title in [
        ("total_counts",      axes[0, 0], "Total counts"),
        ("n_genes_by_counts", axes[0, 1], "N genes/cell"),
        ("cell_area",         axes[1, 0], "Cell area"),
    ]:
        vmax = float(np.quantile(adata.obs[col], 0.95))
        vmin = float(adata.obs[col].min())
        sc.pl.umap(adata, color=col, ax=ax, show=False, cmap="viridis",
                   title=f"{title} (clipped p95={vmax:.0f})",
                   vmin=vmin, vmax=vmax)
    sc.pl.umap(adata, color="leiden_scvi", ax=axes[1, 1], show=False,
               title="leiden_scvi", legend_loc="on data",
               legend_fontsize=11, legend_fontoutline=2)
    plt.tight_layout()
    plt.savefig(fig_dir / "umap_qc_overlay.png", dpi=200, bbox_inches="tight")
    plt.close()

    # ---- final_celltypes_overview.png ----
    print(f"    rendering final_celltypes_overview.png ...")
    fig, axes = plt.subplots(2, 2, figsize=(15, 15))
    sc.pl.umap(adata, color="celltype", ax=axes[0, 0], show=False,
               title="Cell Types", legend_loc="on data",
               legend_fontsize=11, legend_fontoutline=2)
    sc.pl.umap(adata, color="leiden_scvi", ax=axes[0, 1], show=False,
               title="scVI Leiden", legend_loc="on data",
               legend_fontsize=11, legend_fontoutline=2)
    for col, ax, title in [
        ("total_counts",      axes[1, 0], "Total Counts"),
        ("n_genes_by_counts", axes[1, 1], "N Genes"),
    ]:
        vmax = float(np.quantile(adata.obs[col], 0.95))
        vmin = float(adata.obs[col].min())
        sc.pl.umap(adata, color=col, ax=ax, show=False, cmap="viridis",
                   title=f"{title} (clipped p95={vmax:.0f})",
                   vmin=vmin, vmax=vmax)
    plt.tight_layout()
    plt.savefig(fig_dir / "final_celltypes_overview.png", dpi=300, bbox_inches="tight")
    plt.close()
    del adata


for s in SAMPLES:
    render(s)
print("\nDONE")
