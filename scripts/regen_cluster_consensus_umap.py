"""Regenerate cluster_consensus_annotation.png as a two-panel
(general lineage | specific celltype) UMAP for both samples.

Reads the saved phenotyped h5ads (so Gamma/Epsilon rule-based refinements
from scripts/refine_rare_celltypes.py are picked up) and overwrites
figures/{sample}/02_phenotyping/cluster_consensus_annotation.png.

Mirrors the logic baked into notebooks/02_phenotyping.ipynb cell 20.
Run it after the pipeline if you want to refresh the figure without
re-executing the full notebook.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

REPO = Path(__file__).resolve().parent.parent
SAMPLES = ("0041323", "0041326")

LINEAGE_MAP = {
    "Alpha": "Endocrine", "Beta": "Endocrine", "Delta": "Endocrine",
    "Gamma": "Endocrine", "Epsilon": "Endocrine",
    "Endocrine": "Endocrine", "Endocrine_pan": "Endocrine",
    "Acinar": "Exocrine", "Ductal": "Exocrine",
    "T_cells": "Immune", "B_cells": "Immune", "Myeloid": "Immune",
    "Stellate": "Stromal",
    "Endothelial": "Vascular",
    "Schwann": "Neural",
    "Proliferating": "Cycling",
}

SPECIFIC_KEEP = {
    "Alpha", "Beta", "Delta", "Gamma", "Epsilon",
    "Acinar", "Ductal",
    "T_cells", "B_cells", "Myeloid",
    "Stellate", "Endothelial", "Schwann", "Proliferating",
}

EXOCRINE = {"Acinar", "Ductal"}

PALETTE = {
    "Alpha": "#d62728", "Beta": "#ff7f0e", "Delta": "#8c564b",
    "Gamma": "#fa8072", "Epsilon": "#e377c2", "Endocrine": "#b22222",
    "Acinar": "#7fbf3f", "Ductal": "#17becf", "Exocrine": "#2ca02c",
    "T_cells": "#1f3a93", "B_cells": "#87ceeb", "Myeloid": "#4682b4",
    "Immune": "#1f77b4",
    "Stellate": "#9467bd", "Stromal": "#6a3d9a",
    "Endothelial": "#e83e8c", "Vascular": "#e83e8c",
    "Schwann": "#808000", "Neural": "#808000",
    "Proliferating": "#7f7f7f", "Cycling": "#7f7f7f",
}


def _ordered(cats, counts):
    return sorted(cats, key=lambda c: -counts.get(c, 0))


def render_one(sample_id: str) -> None:
    h5 = REPO / "data" / "processed" / sample_id / f"{sample_id}_phenotyped.h5ad"
    fig_dir = REPO / "figures" / sample_id / "02_phenotyping"
    out = fig_dir / "cluster_consensus_annotation.png"

    if not h5.exists():
        print(f"[{sample_id}] missing {h5} — skipping")
        return
    fig_dir.mkdir(parents=True, exist_ok=True)

    print(f"[{sample_id}] loading {h5}")
    # backed='r' avoids a 60 GB RAM hit; we only need obs + obsm['X_umap'].
    adata = sc.read_h5ad(h5, backed="r")

    # The saved h5ad has `celltype` populated by cell 22 + post-hoc
    # scripts/refine_rare_celltypes.py (which is where Gamma/Epsilon and
    # Indeterminate_*_lowcoverage land). Fall back to manual_celltype if
    # an older h5ad without `celltype` is encountered.
    label_col = "celltype" if "celltype" in adata.obs.columns else "manual_celltype"
    print(f"[{sample_id}] using obs['{label_col}'] for annotation")

    if "X_umap" not in adata.obsm:
        # scVI UMAP is the canonical embedding for stage-02 figures.
        if "X_umap_scvi" in adata.obsm:
            umap = np.asarray(adata.obsm["X_umap_scvi"])
        else:
            raise KeyError(f"[{sample_id}] no X_umap or X_umap_scvi on h5ad")
    else:
        umap = np.asarray(adata.obsm["X_umap"])

    labels = pd.Series(adata.obs[label_col].astype(str).to_numpy(),
                       index=np.arange(adata.n_obs))
    keep_general = ~labels.str.startswith("Indeterminate")
    keep_specific = labels.isin(SPECIFIC_KEEP)
    keep_nonexo = keep_specific & ~labels.isin(EXOCRINE)

    # Build minimal in-memory AnnDatas for the two panels (no expression matrix).
    def _slim(mask, color_col, color_vals):
        import anndata as ad
        sub = ad.AnnData(
            X=np.zeros((int(mask.sum()), 1), dtype=np.float32),
            obs=pd.DataFrame({color_col: pd.Categorical(color_vals[mask.to_numpy()])},
                             index=np.arange(int(mask.sum())).astype(str)),
        )
        sub.obsm["X_umap"] = umap[mask.to_numpy()]
        return sub

    general_vals = labels.map(LINEAGE_MAP).fillna("Other")
    sub_g = _slim(keep_general, "_lineage", general_vals)
    sub_s = _slim(keep_specific, "_specific", labels)
    sub_n = _slim(keep_nonexo, "_specific", labels)

    g_counts = sub_g.obs["_lineage"].value_counts().to_dict()
    s_counts = sub_s.obs["_specific"].value_counts().to_dict()
    n_counts = sub_n.obs["_specific"].value_counts().to_dict()
    sub_g.obs["_lineage"] = sub_g.obs["_lineage"].cat.reorder_categories(
        _ordered(list(sub_g.obs["_lineage"].cat.categories), g_counts), ordered=False)
    sub_s.obs["_specific"] = sub_s.obs["_specific"].cat.reorder_categories(
        _ordered(list(sub_s.obs["_specific"].cat.categories), s_counts), ordered=False)
    sub_n.obs["_specific"] = sub_n.obs["_specific"].cat.reorder_categories(
        _ordered(list(sub_n.obs["_specific"].cat.categories), n_counts), ordered=False)

    print(f"[{sample_id}]   general panel cells: {sub_g.n_obs:,} "
          f"({len(g_counts)} classes: {list(sub_g.obs['_lineage'].cat.categories)})")
    print(f"[{sample_id}]   specific panel cells: {sub_s.n_obs:,} "
          f"({len(s_counts)} classes: {list(sub_s.obs['_specific'].cat.categories)})")
    print(f"[{sample_id}]   no-exocrine panel cells: {sub_n.n_obs:,} "
          f"({len(n_counts)} classes: {list(sub_n.obs['_specific'].cat.categories)})")

    fig, axes = plt.subplots(1, 3, figsize=(28, 9))
    sc.pl.umap(sub_g, color="_lineage", ax=axes[0], show=False,
               size=3, frameon=False, legend_loc="right margin",
               legend_fontsize=9, palette=PALETTE,
               title=f"General lineage (n={sub_g.n_obs:,})")
    sc.pl.umap(sub_s, color="_specific", ax=axes[1], show=False,
               size=3, frameon=False, legend_loc="right margin",
               legend_fontsize=9, palette=PALETTE,
               title=f"Specific celltype (n={sub_s.n_obs:,})")
    sc.pl.umap(sub_n, color="_specific", ax=axes[2], show=False,
               size=4, frameon=False, legend_loc="right margin",
               legend_fontsize=9, palette=PALETTE,
               title=f"Specific celltype, no exocrine (n={sub_n.n_obs:,})")
    for ax in axes:
        for coll in ax.collections:
            coll.set_rasterized(True)
    fig.suptitle(f"{sample_id} — cluster consensus annotation", y=0.99)
    plt.tight_layout()
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[{sample_id}] wrote {out}")

    adata.file.close()


def main() -> None:
    for s in SAMPLES:
        render_one(s)


if __name__ == "__main__":
    main()
