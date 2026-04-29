"""Re-run sq.gr.co_occurrence on the saved spatial_analysis h5ad (which has
post-fix-acinar celltype labels) and render with a cleaner layout that fixes:
  * y-tick labels overlapping the next subplot (sharey + 2-row grid)
  * verbose squidpy default titles like '$\\frac{p(exp|X)}{p(exp)}$'
    -> simplified to just the cluster name
  * one legend instead of repeating it on every subplot

Usage:
    python scripts/rerender_co_occurrence.py                # both samples
    python scripts/rerender_co_occurrence.py 0041323        # one sample
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
import scanpy as sc
import squidpy as sq

PROJECT = Path("/blue/maigan/smith6jt/Xenium_Analysis")
DEFAULT_SAMPLES = ["0041323", "0041326"]
N_SUB_CO = 150_000


def plot_co_occurrence_clean(adata, cluster_key="celltype", out_path=None,
                             ncols=4, panel_size=(3.6, 2.8)):
    """Custom co-occurrence renderer with clean axes and titles.

    Mirrors sq.pl.co_occurrence's data layout (out[idx, :, :].T per subplot,
    with one line per source cluster) but uses a 2D subplot grid with sharey
    so y-tick labels appear only on the leftmost column.
    """
    co = adata.uns[f"{cluster_key}_co_occurrence"]
    out = co["occ"]                 # (n_clusters, n_clusters, n_intervals)
    interval = co["interval"][1:]   # right edge of each distance bin

    cats = list(adata.obs[cluster_key].cat.categories)
    n = len(cats)
    nrows = (n + ncols - 1) // ncols

    # Use the per-cluster colors that scanpy/squidpy pickled in adata.uns,
    # falling back to tab20 if absent.
    palette_key = f"{cluster_key}_colors"
    if palette_key in adata.uns:
        palette = list(adata.uns[palette_key])
        if len(palette) < n:
            palette = list(plt.cm.tab20(np.linspace(0, 1, n)))
    else:
        palette = list(plt.cm.tab20(np.linspace(0, 1, n)))
    color_map = {cat: palette[i] for i, cat in enumerate(cats)}

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(panel_size[0] * ncols, panel_size[1] * nrows),
        sharey=True, sharex=True,
    )
    axes = np.atleast_1d(axes).flatten()

    for j, target in enumerate(cats):
        ax = axes[j]
        for i, source in enumerate(cats):
            ax.plot(interval, out[j, i, :], color=color_map[source], lw=1.2)
        ax.set_title(target, fontsize=10)
        ax.axhline(1.0, color="gray", linestyle="--", linewidth=0.6, alpha=0.6)
        ax.tick_params(labelsize=8)
        if j % ncols == 0:
            ax.set_ylabel("p(co-occur) / p(marginal)", fontsize=9)
        if j >= n - ncols:
            ax.set_xlabel("distance (μm)", fontsize=9)

    for j in range(n, len(axes)):
        axes[j].axis("off")

    # Single legend on the right
    handles = [plt.Line2D([0], [0], color=color_map[c], lw=2, label=c) for c in cats]
    fig.legend(handles=handles, loc="center right",
               bbox_to_anchor=(1.0, 0.5), frameon=False,
               title=cluster_key, fontsize=9, title_fontsize=10)

    fig.suptitle(f"Co-occurrence ratio vs distance, by {cluster_key}",
                 fontsize=12, y=0.99)
    fig.tight_layout(rect=(0, 0, 0.88, 0.97))

    if out_path is not None:
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
    return fig


def render(sample: str) -> None:
    h5ad = PROJECT / f"data/processed/{sample}/{sample}_spatial_analysis.h5ad"
    fig_dir = PROJECT / f"figures/{sample}/03_spatial_analysis"
    fig_dir.mkdir(parents=True, exist_ok=True)
    out_png = fig_dir / "co_occurrence.png"

    print(f"\n=== {sample} ===")
    t0 = time.time()
    a = sc.read_h5ad(h5ad)
    print(f"  loaded {a.shape} in {time.time()-t0:.1f}s")
    print(f"  celltype categories: {a.obs['celltype'].nunique()}")

    rng = np.random.default_rng(0)
    if a.n_obs > N_SUB_CO:
        sub_idx = rng.choice(a.n_obs, size=N_SUB_CO, replace=False)
        a_co = a[sub_idx].copy()
        print(f"  co_occurrence on {N_SUB_CO:,}/{a.n_obs:,} cells")
    else:
        a_co = a

    print("  running sq.gr.co_occurrence ...")
    t = time.time()
    sq.gr.co_occurrence(a_co, cluster_key="celltype")
    print(f"  done in {time.time()-t:.0f}s")

    plot_co_occurrence_clean(a_co, cluster_key="celltype", out_path=out_png)
    print(f"  wrote {out_png}")
    print(f"  total: {time.time()-t0:.0f}s")


def main():
    samples = sys.argv[1:] if len(sys.argv) > 1 else DEFAULT_SAMPLES
    for s in samples:
        render(s)
    print("\nDONE")


if __name__ == "__main__":
    main()
