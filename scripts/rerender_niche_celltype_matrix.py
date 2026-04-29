"""Render a spatial-niche × cell-type matrix plot from the saved
{sample}_spatial_analysis.h5ad. Produces a 2-panel figure:
  left  : niche composition (rows = niches sum to 1)
  right : celltype distribution across niches (rows = celltypes sum to 1)

Saves to figures/{sample}/03_spatial_analysis/niche_celltype_matrix.png.

Usage:
    python scripts/rerender_niche_celltype_matrix.py                  # both
    python scripts/rerender_niche_celltype_matrix.py 0041323          # one
"""

from __future__ import annotations
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad

PROJECT = Path("/blue/maigan/smith6jt/Xenium_Analysis")
DEFAULT_SAMPLES = ["0041323", "0041326"]


def render(sample: str) -> None:
    h5ad = PROJECT / f"data/processed/{sample}/{sample}_spatial_analysis.h5ad"
    fig_dir = PROJECT / f"figures/{sample}/03_spatial_analysis"
    fig_dir.mkdir(parents=True, exist_ok=True)
    out_png = fig_dir / "niche_celltype_matrix.png"

    print(f"\n=== {sample} ===")
    t0 = time.time()
    a = ad.read_h5ad(h5ad)
    print(f"  loaded {a.shape} in {time.time()-t0:.1f}s")

    if "spatial_niche" not in a.obs.columns:
        raise SystemExit(f"  [{sample}] obs['spatial_niche'] missing")
    if "celltype" not in a.obs.columns:
        raise SystemExit(f"  [{sample}] obs['celltype'] missing")

    niche = a.obs["spatial_niche"].astype(str)
    celltype = a.obs["celltype"].astype(str)

    # Niche composition: each row is one niche; cells = fraction of that niche
    # made up by each cell type. Row sums = 1.
    comp = pd.crosstab(niche, celltype, normalize="index").fillna(0.0)
    comp = comp.reindex(sorted(comp.index, key=int))  # niches 0..11 in order

    # Celltype distribution: each row is one cell type; cells = fraction of
    # that cell type living in each niche. Row sums = 1.
    dist = pd.crosstab(niche, celltype, normalize="columns").fillna(0.0).T
    # Order columns niche-numerically; rows by total abundance (most common first)
    dist = dist.reindex(columns=sorted(dist.columns, key=int))
    abund = celltype.value_counts()
    dist = dist.reindex(index=[c for c in abund.index if c in dist.index])

    # Figure: two heatmaps side-by-side
    fig, axes = plt.subplots(
        1, 2, figsize=(22, 8),
        gridspec_kw={"width_ratios": [1.05, 1.0], "wspace": 0.25},
    )

    # --- Left: niche composition (rows = niches, columns = cell types) ---
    sns.heatmap(
        comp, cmap="viridis", ax=axes[0],
        annot=True, fmt=".2f", annot_kws={"size": 7},
        cbar_kws={"label": "fraction of niche", "shrink": 0.7},
        linewidths=0.4, linecolor="white",
    )
    axes[0].set_title("Niche composition  (rows sum = 1)", fontsize=12)
    axes[0].set_xlabel("cell type")
    axes[0].set_ylabel("spatial niche")
    axes[0].tick_params(axis="x", labelrotation=45, labelsize=9)
    for label in axes[0].get_xticklabels():
        label.set_horizontalalignment("right")
    axes[0].tick_params(axis="y", labelrotation=0, labelsize=10)

    # --- Right: celltype distribution across niches (rows = cell types) ---
    sns.heatmap(
        dist, cmap="magma", ax=axes[1],
        annot=True, fmt=".2f", annot_kws={"size": 7},
        cbar_kws={"label": "fraction of celltype", "shrink": 0.7},
        linewidths=0.4, linecolor="white",
    )
    axes[1].set_title("Cell-type distribution across niches  (rows sum = 1)",
                      fontsize=12)
    axes[1].set_xlabel("spatial niche")
    axes[1].set_ylabel("cell type")
    axes[1].tick_params(axis="x", labelrotation=0, labelsize=10)
    axes[1].tick_params(axis="y", labelrotation=0, labelsize=10)

    fig.suptitle(f"Spatial niche × cell type composition — {sample}",
                 fontsize=13, y=1.02)
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out_png}")
    print(f"  total: {time.time()-t0:.1f}s")


def main():
    samples = sys.argv[1:] if len(sys.argv) > 1 else DEFAULT_SAMPLES
    for s in samples:
        render(s)
    print("\nDONE")


if __name__ == "__main__":
    main()
