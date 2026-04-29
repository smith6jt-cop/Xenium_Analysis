"""Re-render figures/{sample}/03_spatial_analysis/spatial_domains.png from
the saved {sample}_spatial_analysis.h5ad — picks up the new niche legend
without requiring a full re-run of stage 03.

Usage:
    python scripts/rerender_spatial_domains.py                 # both samples
    python scripts/rerender_spatial_domains.py 0041323         # one sample
"""

from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
import anndata as ad

PROJECT = Path("/blue/maigan/smith6jt/Xenium_Analysis")
DEFAULT_SAMPLES = ["0041323", "0041326"]
N_PLOT = 300_000  # match notebook subsample
K_COMP_LABEL = 50


def render(sample: str) -> None:
    h5ad = PROJECT / f"data/processed/{sample}/{sample}_spatial_analysis.h5ad"
    fig_dir = PROJECT / f"figures/{sample}/03_spatial_analysis"
    fig_dir.mkdir(parents=True, exist_ok=True)
    out_png = fig_dir / "spatial_domains.png"

    print(f"[{sample}] loading {h5ad.name}")
    a = ad.read_h5ad(h5ad)

    if "spatial" not in a.obsm:
        raise SystemExit(f"[{sample}] obsm['spatial'] missing")
    if "spatial_niche" not in a.obs.columns:
        raise SystemExit(f"[{sample}] obs['spatial_niche'] missing")

    coords = np.asarray(a.obsm["spatial"], dtype=np.float64)
    niche_str = a.obs["spatial_niche"].astype(str).to_numpy()
    smoothed = niche_str.astype(np.int32)
    n_niches = int(smoothed.max()) + 1
    n_cells = coords.shape[0]
    print(f"[{sample}] n_cells={n_cells:,}  n_niches={n_niches}")

    celltypes = a.obs["celltype"].astype("category")
    ct_cats = list(celltypes.cat.categories)
    sub_codes_full = celltypes.cat.codes.to_numpy()

    rng = np.random.default_rng(0)
    idx = rng.choice(n_cells, size=min(N_PLOT, n_cells), replace=False)

    tab20 = mpl.colormaps.get_cmap("tab20")
    niche_colors = [tab20(i % 20) for i in range(n_niches)]
    niche_cmap = ListedColormap(niche_colors)

    # 4-column layout: niche map | niche legend | celltype map | celltype legend
    fig, axes = plt.subplots(
        1, 4, figsize=(24, 8),
        gridspec_kw={
            "width_ratios": [1.0, 0.28, 1.0, 0.38],
            "wspace": 0.04,
        },
    )
    ax_niche, ax_leg, ax_ct, ax_ct_leg = axes
    ax_leg.axis("off")
    ax_ct_leg.axis("off")

    ax_niche.scatter(coords[idx, 0], coords[idx, 1],
                     c=smoothed[idx], cmap=niche_cmap,
                     vmin=-0.5, vmax=n_niches - 0.5,
                     s=0.5, linewidths=0)
    ax_niche.set_title(
        f"Spatial niches (smoothed, K_COMP={K_COMP_LABEL}, n_niches={n_niches})"
    )
    ax_niche.set_aspect("equal", adjustable="datalim")

    niche_sizes = np.bincount(smoothed, minlength=n_niches)
    niche_handles = [
        Line2D([0], [0], marker="o", linestyle="", markersize=9,
               markerfacecolor=niche_colors[i], markeredgecolor="none",
               label=f"Niche {i}  (n={niche_sizes[i]:,})")
        for i in range(n_niches)
    ]
    ax_leg.legend(handles=niche_handles, title="Spatial niche",
                  loc="center", frameon=False,
                  fontsize=9, title_fontsize=10,
                  handletextpad=0.6, labelspacing=1.1,
                  borderaxespad=0.0)

    # Cell-types panel: ListedColormap for stable color↔label mapping (so the
    # legend swatches match the scatter, mirroring the niche panel approach).
    n_ct = len(ct_cats)
    ct_colors = [tab20(i % 20) for i in range(n_ct)]
    ct_cmap = ListedColormap(ct_colors)

    sub_codes = sub_codes_full[idx]
    ax_ct.scatter(coords[idx, 0], coords[idx, 1],
                  c=sub_codes, cmap=ct_cmap,
                  vmin=-0.5, vmax=n_ct - 0.5,
                  s=0.5, linewidths=0)
    ax_ct.set_title("Cell types")
    ax_ct.set_aspect("equal", adjustable="datalim")
    ax_ct.tick_params(axis="y", left=False, labelleft=False)

    ct_sizes = np.bincount(sub_codes_full, minlength=n_ct)
    ct_handles = [
        Line2D([0], [0], marker="o", linestyle="", markersize=9,
               markerfacecolor=ct_colors[i], markeredgecolor="none",
               label=f"{ct_cats[i]}  (n={ct_sizes[i]:,})")
        for i in range(n_ct)
    ]
    ax_ct_leg.legend(handles=ct_handles, title="Cell type",
                     loc="center", frameon=False,
                     fontsize=9, title_fontsize=10,
                     handletextpad=0.6, labelspacing=1.1,
                     borderaxespad=0.0)

    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[{sample}] wrote {out_png}")


def main():
    samples = sys.argv[1:] if len(sys.argv) > 1 else DEFAULT_SAMPLES
    for s in samples:
        render(s)


if __name__ == "__main__":
    main()
