"""Re-render dotplot_markers_by_cluster.png for each sample with a minimum
dot size so low-expression markers stay visible.

Mirrors notebook 02 cell 11 logic (groupby='leiden_1.5', dendrogram=True,
standard_scale='var') but adds smallest_dot=15 + mean_only_expressed=True.
smallest_dot=15 keeps low-expression dots visible while leaving clear
visual contrast with the largest dots (default ~200 area). The earlier
value of 40 was overcorrected — smallest dots looked nearly as large as
mid-range dots and washed out the size gradient.
Uses the same panel-aware pancreas_markers dict and panel_for_scoring mask
as the patched cell 10.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path("/blue/maigan/smith6jt/Xenium_Analysis")
SAMPLES = ("0041323", "0041326")

# Same panel-aware dict as the post-Acinar-fix notebook patch.
PANCREAS_MARKERS = {
    'Beta':          ['IAPP', 'MAFA', 'NKX6-1', 'PDX1', 'SLC30A8', 'ABCC8',
                       'KCNJ11', 'GLP1R', 'PCSK1'],
    'Alpha':         ['ARX', 'MAFB', 'GC', 'TM4SF4', 'TTR'],
    'Delta':         ['HHEX', 'LEPR', 'RBP4', 'BHLHE41', 'PCSK1'],
    'Endocrine':     ['CHGA', 'INSM1', 'ISL1', 'NEUROD1', 'FEV', 'PAX6', 'SCG5',
                       'SCG2', 'SCGN'],
    'Acinar':        ['AMY1A', 'CUZD1'],
    'Ductal':        ['KRT19', 'SOX9', 'HNF1B', 'MUC1', 'CFTR', 'KRT7'],
    'Stellate':      ['ACTA2', 'PDGFRB', 'RGS5', 'PDGFRA', 'DCN', 'COL1A1'],
    'Endothelial':   ['PECAM1', 'CDH5', 'CLDN5', 'PLVAP', 'KDR', 'ENG', 'VWF'],
    'T_cells':       ['CD3D', 'CD3E', 'CD3G', 'CD8A', 'CD8B', 'CD4', 'PTPRC',
                       'FOXP3'],
    'B_cells':       ['MS4A1', 'CD79A', 'CD79B', 'CD19'],
    'Myeloid':       ['CD68', 'CD163', 'CSF1R', 'MARCO', 'CD14'],
    'Schwann':       ['MPZ', 'SOX10', 'NES', 'S100B', 'PMP22', 'PLP1'],
    'Proliferating': ['MKI67', 'TOP2A', 'PCNA'],
}


def _add_group_separators(dp, group_dict, color="black"):
    if dp.ax_dict is None:
        dp.make_figure()
    main_ax = dp.ax_dict.get("mainplot_ax")
    if main_ax is None:
        return
    cum = 0
    items = list(group_dict.items())
    for i, (_, genes) in enumerate(items):
        cum += len(genes)
        if i < len(items) - 1:
            main_ax.axvline(cum, color=color, linewidth=1.6, alpha=0.95, zorder=10)


def render_one(sample):
    print(f"\n=== {sample} ===")
    p = ROOT / f"data/processed/{sample}/{sample}_phenotyped.h5ad"
    a = sc.read_h5ad(p)
    print(f"    shape: {a.shape}")

    panel_keep = (set(a.var_names[a.var["panel_for_scoring"].astype(bool)])
                    if "panel_for_scoring" in a.var.columns
                    else set(a.var_names))
    available = {ct: [g for g in genes if g in a.var_names and g in panel_keep]
                   for ct, genes in PANCREAS_MARKERS.items()}
    available = {k: v for k, v in available.items() if len(v) >= 2}
    print(f"    types in dotplot: {list(available.keys())}")

    if "leiden_1.5" not in a.obs.columns:
        print(f"    skipping — no leiden_1.5 column")
        return

    fig_dir = ROOT / f"figures/{sample}/02_phenotyping"
    fig_dir.mkdir(parents=True, exist_ok=True)

    dp = sc.pl.dotplot(
        a, available, groupby="leiden_1.5",
        dendrogram=True, cmap="Reds", figsize=(18, 10),
        vmin=0.0, vmax=1.0, standard_scale="var",
        swap_axes=False, var_group_rotation=90, show=False,
        return_fig=True,
        dot_min=0.0, dot_max=1.0,
        smallest_dot=15,
        mean_only_expressed=True,
    )
    _add_group_separators(dp, available, color="black")
    out = fig_dir / "dotplot_markers_by_cluster.png"
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"    wrote {out}")


def main():
    for s in SAMPLES:
        render_one(s)
    print("\nDONE")


if __name__ == "__main__":
    main()
