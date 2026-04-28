"""Re-render marker dotplots/heatmaps with corrected separator positions
and color (white on dark cmap, black on light cmap).

For stage 01 (uses _preprocessed.h5ad, ~7 GB):
  - marker_heatmap.png   (sc.pl.matrixplot, cmap=magma → white separators)
  - marker_dotplot.png   (sc.pl.dotplot,    cmap=Reds  → black separators)

For stage 02 (uses _phenotyped.h5ad, ~8.5 GB):
  - marker_genes_heatmap.png       (sc.pl.rank_genes_groups_heatmap, viridis → white)
  - dotplot_markers_by_cluster.png (sc.pl.dotplot, Reds → black)
"""
from pathlib import Path
import time

import scanpy as sc
import matplotlib.pyplot as plt

ROOT = Path("/blue/maigan/smith6jt/Xenium_Analysis")
SAMPLES = ("0041323", "0041326")

PANCREAS_MARKERS = {
    "Beta":          ["IAPP", "MAFA", "NKX6-1", "PDX1", "SLC30A8", "ABCC8", "KCNJ11", "GLP1R", "PCSK1"],
    "Alpha":         ["ARX", "MAFB"],
    "Delta":         ["HHEX", "LEPR"],
    "Epsilon":       ["GHRL"],
    "Endocrine":     ["CHGA", "INSM1", "ISL1", "NEUROD1", "FEV", "PAX6", "SCG5"],
    "Acinar":        ["AMY1A"],
    "Ductal":        ["KRT19", "SOX9", "HNF1B", "MUC1", "CFTR"],
    "Stellate":      ["ACTA2", "PDGFRB", "RGS5", "PDGFRA"],
    "Endothelial":   ["PECAM1", "CDH5", "CLDN5", "PLVAP", "KDR", "ENG"],
    "T_cells":       ["CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "CD4", "PTPRC", "FOXP3"],
    "B_cells":       ["MS4A1", "CD79A", "CD79B", "CD19"],
    "Myeloid":       ["CD68", "CD163", "CSF1R", "MARCO", "CD14"],
    "Schwann":       ["MPZ", "SOX10", "NES"],
    "Proliferating": ["MKI67", "TOP2A", "PCNA"],
}


def add_separators(dp_or_mp, group_dict_or_n_per_group, n_groups=None,
                    color="black"):
    """Draw vertical lines BETWEEN gene groups on a scanpy DotPlot/MatrixPlot.

    Scanpy places genes at x=0.5, 1.5, 2.5, ... — boundaries between
    groups are at integer cumulative gene counts.
    """
    if dp_or_mp.ax_dict is None:
        dp_or_mp.make_figure()
    main_ax = dp_or_mp.ax_dict.get("mainplot_ax")
    if main_ax is None:
        return
    if isinstance(group_dict_or_n_per_group, dict):
        cum = 0
        items = list(group_dict_or_n_per_group.items())
        for i, (_, genes) in enumerate(items):
            cum += len(genes)
            if i < len(items) - 1:
                main_ax.axvline(cum, color=color, linewidth=1.6, alpha=0.95, zorder=10)
    else:
        n_per = int(group_dict_or_n_per_group)
        for k in range(1, int(n_groups)):
            main_ax.axvline(k * n_per, color=color, linewidth=1.6, alpha=0.95, zorder=10)


def stage_01(sample):
    p = ROOT / f"data/processed/{sample}/{sample}_preprocessed.h5ad"
    fig_dir = ROOT / f"figures/{sample}/01_preprocessing"
    print(f"\n  --- 01 / {sample} ---")
    t = time.time()
    adata = sc.read_h5ad(p)
    print(f"    loaded {adata.shape} in {time.time()-t:.1f}s")
    marker_genes = {ct: [g for g in genes if g in adata.var_names]
                     for ct, genes in PANCREAS_MARKERS.items()
                     if any(g in adata.var_names for g in genes)}

    # Matrixplot (magma → white separators)
    print(f"    rendering marker_heatmap.png (matrixplot, magma, white seps)")
    mp = sc.pl.matrixplot(
        adata, var_names=marker_genes, groupby="leiden_1.5",
        dendrogram=True, cmap="magma", figsize=(22, 10),
        vmin=0.0, vmax=1.0, standard_scale="var",
        swap_axes=False, var_group_rotation=90,
        colorbar_title="scaled mean\nexpression",
        show=False, return_fig=True,
    )
    add_separators(mp, marker_genes, color="white")
    plt.savefig(fig_dir / "marker_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close()

    # Dotplot (Reds → black separators)
    print(f"    rendering marker_dotplot.png (dotplot, Reds, black seps)")
    dp = sc.pl.dotplot(
        adata, var_names=marker_genes, groupby="leiden_1.5",
        dendrogram=True, cmap="Reds", figsize=(18, 9),
        vmin=0.0, vmax=1.0, standard_scale="var",
        swap_axes=False, var_group_rotation=90,
        colorbar_title="scaled mean\nexpression",
        size_title="fraction of cells\nin group (%)",
        show=False, return_fig=True,
    )
    add_separators(dp, marker_genes, color="black")
    plt.savefig(fig_dir / "marker_dotplot.png", dpi=300, bbox_inches="tight")
    plt.close()
    del adata


def stage_02(sample):
    p = ROOT / f"data/processed/{sample}/{sample}_phenotyped.h5ad"
    fig_dir = ROOT / f"figures/{sample}/02_phenotyping"
    print(f"\n  --- 02 / {sample} ---")
    t = time.time()
    adata = sc.read_h5ad(p)
    print(f"    loaded {adata.shape} in {time.time()-t:.1f}s")
    available_markers = {ct: [g for g in genes if g in adata.var_names]
                          for ct, genes in PANCREAS_MARKERS.items()
                          if any(g in adata.var_names for g in genes)}

    # Heatmap by leiden_1.5 (rank_genes_groups_heatmap, viridis → white)
    if "rank_genes_leiden" in adata.uns:
        print(f"    rendering marker_genes_heatmap.png (heatmap, viridis, white seps)")
        hm = sc.pl.rank_genes_groups_heatmap(
            adata, n_genes=5, key="rank_genes_leiden", show_gene_labels=True,
            figsize=(28, 28), cmap="viridis", vmin=-1, vmax=3,
            layer="lognorm", swap_axes=True, var_group_rotation=90, show=False,
            return_fig=True,
        )
        n_clusters = (
            adata.obs["leiden_1.5"].astype("category").cat.categories.size
            if "leiden_1.5" in adata.obs.columns else 0
        )
        if n_clusters >= 2 and hasattr(hm, "ax_dict"):
            add_separators(hm, 5, n_groups=n_clusters, color="white")
        plt.tight_layout()
        plt.savefig(fig_dir / "marker_genes_heatmap.png", dpi=300, bbox_inches="tight")
        plt.close()

    # Dotplot by leiden_1.5 (Reds → black)
    if "leiden_1.5" in adata.obs.columns:
        print(f"    rendering dotplot_markers_by_cluster.png (dotplot, Reds, black seps)")
        dp = sc.pl.dotplot(
            adata, available_markers, groupby="leiden_1.5",
            dendrogram=True, cmap="Reds", figsize=(18, 10),
            vmin=0.0, vmax=1.0, standard_scale="var",
            swap_axes=False, var_group_rotation=90, show=False,
            return_fig=True,
        )
        add_separators(dp, available_markers, color="black")
        plt.savefig(fig_dir / "dotplot_markers_by_cluster.png", dpi=300, bbox_inches="tight")
        plt.close()
    del adata


for s in SAMPLES:
    stage_01(s)
    stage_02(s)
print("\nDONE")
