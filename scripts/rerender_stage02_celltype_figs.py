"""Refresh the celltype-dependent stage-02 UMAP figures from the saved
{sample}_phenotyped.h5ad. Needed because the originals were rendered before
scripts/fix_acinar_argmax.py rewrote the celltype labels — old figures show
B_cells inflated to ~464K cells, the bug from the [AMY1A]-only Acinar panel.

Refreshes:
    cluster_consensus_annotation.png   (manual_celltype on scVI UMAP)
    final_celltypes_overview.png       (celltype + leiden_scvi + counts/genes)
    predicted_celltypes.png            (predicted_celltype + confidence)
    scvi_vs_markers.png                (leiden_scvi vs predicted_celltype)
    celltype_scores_umap.png           (per-celltype *_score columns)

Usage:
    python scripts/rerender_stage02_celltype_figs.py                 # both
    python scripts/rerender_stage02_celltype_figs.py 0041323         # one
"""

from __future__ import annotations
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

PROJECT = Path("/blue/maigan/smith6jt/Xenium_Analysis")
DEFAULT_SAMPLES = ["0041323", "0041326"]


def render(sample: str) -> None:
    h5ad = PROJECT / f"data/processed/{sample}/{sample}_phenotyped.h5ad"
    fig_dir = PROJECT / f"figures/{sample}/02_phenotyping"
    fig_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n=== {sample} ===")
    t0 = time.time()
    a = sc.read_h5ad(h5ad)
    print(f"  loaded {a.shape} in {time.time()-t0:.1f}s")

    # Always plot on the scVI UMAP. Stage 02 stores X_umap == X_umap_scvi (cell 22),
    # but be defensive and prefer X_umap_scvi if present.
    if "X_umap_scvi" in a.obsm:
        a.obsm["X_umap"] = a.obsm["X_umap_scvi"]
    elif "X_umap" not in a.obsm:
        raise SystemExit(f"  [{sample}] no X_umap or X_umap_scvi in obsm — cannot render")

    # 1) cluster_consensus_annotation.png
    if "manual_celltype" in a.obs.columns:
        sc.pl.umap(a, color="manual_celltype", show=False, frameon=False, size=2)
        out = fig_dir / "cluster_consensus_annotation.png"
        plt.savefig(out, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"  wrote {out.name}  ({a.obs['manual_celltype'].nunique()} types)")

    # 2) final_celltypes_overview.png — 2x2 with celltype, leiden_scvi, counts, genes
    if "celltype" in a.obs.columns:
        fig, axes = plt.subplots(2, 2, figsize=(15, 15))
        sc.pl.umap(a, color="celltype", ax=axes[0, 0], show=False, size=2,
                   title="Cell Types", legend_loc="on data",
                   legend_fontsize=10, legend_fontoutline=2)
        sc.pl.umap(a, color="leiden_scvi", ax=axes[0, 1], show=False, size=2,
                   title="scVI Leiden", legend_loc="on data",
                   legend_fontsize=10, legend_fontoutline=2)
        for col, ax, title in [
            ("total_counts", axes[1, 0], "Total Counts"),
            ("n_genes_by_counts", axes[1, 1], "N Genes"),
        ]:
            if col in a.obs.columns:
                vmax = float(np.quantile(a.obs[col], 0.95))
                sc.pl.umap(a, color=col, ax=ax, show=False, cmap="viridis", size=2,
                           title=f"{title} (clipped p95={vmax:.0f})",
                           vmin=float(a.obs[col].min()), vmax=vmax)
            else:
                ax.axis("off")
        plt.tight_layout()
        out = fig_dir / "final_celltypes_overview.png"
        plt.savefig(out, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"  wrote {out.name}  ({a.obs['celltype'].nunique()} types)")

    # 3) predicted_celltypes.png — 1x2 with predicted_celltype + confidence
    if "predicted_celltype" in a.obs.columns:
        fig, axes = plt.subplots(1, 2, figsize=(15, 5))
        sc.pl.umap(a, color="predicted_celltype", ax=axes[0], show=False, size=2,
                   frameon=False)
        if "celltype_confidence" in a.obs.columns:
            sc.pl.umap(a, color="celltype_confidence", ax=axes[1], show=False,
                       cmap="viridis", size=2, frameon=False)
        else:
            axes[1].axis("off")
        plt.tight_layout()
        out = fig_dir / "predicted_celltypes.png"
        plt.savefig(out, dpi=200, bbox_inches="tight")
        plt.close()
        print(f"  wrote {out.name}")

    # 4) scvi_vs_markers.png — 1x2 with leiden_scvi + predicted_celltype
    if "predicted_celltype" in a.obs.columns and "leiden_scvi" in a.obs.columns:
        fig, axes = plt.subplots(1, 2, figsize=(15, 5))
        sc.pl.umap(a, color="leiden_scvi", ax=axes[0], show=False, size=2,
                   title="scVI Leiden")
        sc.pl.umap(a, color="predicted_celltype", ax=axes[1], show=False, size=2,
                   title="Marker-based prediction")
        plt.tight_layout()
        out = fig_dir / "scvi_vs_markers.png"
        plt.savefig(out, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"  wrote {out.name}")

    # 5) celltype_scores_umap.png — per-score UMAP grid
    score_cols = [c for c in a.obs.columns if c.endswith("_score")]
    if score_cols:
        n = len(score_cols)
        ncols = 3
        nrows = (n + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
        axes = np.atleast_1d(axes).flatten()
        for i, s in enumerate(score_cols):
            sc.pl.umap(a, color=s, ax=axes[i], show=False, cmap="viridis",
                       title=s, size=2, frameon=False)
        for j in range(n, len(axes)):
            axes[j].axis("off")
        plt.tight_layout()
        out = fig_dir / "celltype_scores_umap.png"
        plt.savefig(out, dpi=200, bbox_inches="tight")
        plt.close()
        print(f"  wrote {out.name}  ({len(score_cols)} score cols)")

    print(f"  total: {time.time()-t0:.1f}s")


def main():
    samples = sys.argv[1:] if len(sys.argv) > 1 else DEFAULT_SAMPLES
    for s in samples:
        render(s)
    print("\nDONE")


if __name__ == "__main__":
    main()
