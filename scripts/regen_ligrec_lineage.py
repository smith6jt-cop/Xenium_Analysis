"""Re-run squidpy ligand-receptor analysis at the LINEAGE level (general
celltype categories: Endocrine, Exocrine, Immune, Stromal, Vascular, Neural,
Cycling) and render one dotplot per source lineage.

Why lineage rather than celltype? Subtype-level ligrec is unbalanced for this
panel — the rare endocrine sub-types (Beta=17 in 0041326, Gamma=2, Epsilon=3-9)
have too few cells for stable per-cluster mean expression / meaningful
permutation null, and dropping them by a 50-cell threshold means only Alpha
+ Delta + the 'Endocrine' general label survive on the endocrine side.
Aggregating to lineages eliminates that imbalance, gives 7 robust categories,
and shrinks the cluster-pair space ~15× (49 pairs vs 729) so plots are small
and per-source breakouts are readable.

Outputs (per sample):
  data/processed/{sample}/{sample}_ligrec_lineage_means.csv
  data/processed/{sample}/{sample}_ligrec_lineage_pvalues.csv
  data/processed/{sample}/{sample}_ligrec_lineage_metadata.csv
  figures/{sample}/03_spatial_analysis/ligand_receptor_dotplot_{lineage}.png  (one per source)
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq

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

N_SUB_LR = 100_000
N_PERMS = 50
TOP_N_PER_SOURCE = 40
MIN_MEAN = 0.3
PVAL_THRESH = 0.001
ALPHA = 1e-4
P_FLOOR = 1e-12
DPI = 130


def _build_lineage_subset(sample_id: str) -> sc.AnnData:
    """Load phenotyped h5ad, attach celltype_general, drop Indeterminate*,
    return a 100k-cell view ready for sq.gr.ligrec."""
    h5 = REPO / "data" / "processed" / sample_id / f"{sample_id}_phenotyped.h5ad"
    if not h5.exists():
        raise FileNotFoundError(h5)
    print(f"[{sample_id}] loading {h5}")
    adata = sc.read_h5ad(h5)
    print(f"[{sample_id}]   loaded {adata.shape}")

    ct = adata.obs["celltype"].astype(str)
    keep = ~ct.str.startswith("Indeterminate") & ct.isin(LINEAGE_MAP.keys())
    print(f"[{sample_id}]   after dropping Indeterminate*/unmapped: {int(keep.sum()):,}")
    adata = adata[keep.to_numpy()].copy()

    adata.obs["celltype_general"] = pd.Categorical(
        adata.obs["celltype"].astype(str).map(LINEAGE_MAP),
        categories=sorted(set(LINEAGE_MAP.values())),
    )
    print(f"[{sample_id}]   lineage breakdown:")
    for lin, n in adata.obs["celltype_general"].value_counts().items():
        print(f"      {lin:<10} {n:>9,}")

    if adata.n_obs > N_SUB_LR:
        rng = np.random.default_rng(1)
        sub_idx = rng.choice(adata.n_obs, size=N_SUB_LR, replace=False)
        adata = adata[sub_idx].copy()
        print(f"[{sample_id}]   subsampled to {adata.n_obs:,}")

    # squidpy ligrec needs lognorm in .X.
    if "lognorm" in adata.layers:
        adata.X = adata.layers["lognorm"].copy()
    elif "log_normalized" in adata.layers:
        adata.X = adata.layers["log_normalized"].copy()
    else:
        if "counts" in adata.layers:
            adata.X = adata.layers["counts"].copy()
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
    return adata


def _save_csvs(sample_id: str, lr: dict) -> None:
    proc = REPO / "data" / "processed" / sample_id
    for key in ("means", "pvalues", "metadata"):
        if key not in lr:
            continue
        df = lr[key]
        if not hasattr(df, "to_csv"):
            continue
        out = proc / f"{sample_id}_ligrec_lineage_{key}.csv"
        df.to_csv(out)
        print(f"[{sample_id}]   wrote {out.name}  ({df.shape})")


def _render_per_source(sample_id: str, lr: dict, lineages: list[str]) -> None:
    fig_dir = REPO / "figures" / sample_id / "03_spatial_analysis"
    fig_dir.mkdir(parents=True, exist_ok=True)

    means = lr["means"]
    pvals = lr["pvalues"]
    metadata = lr.get("metadata")

    pvals_clean = np.where(np.isfinite(pvals.to_numpy(dtype=np.float64)),
                           np.maximum(pvals.to_numpy(dtype=np.float64), P_FLOOR),
                           1.0)
    means_arr = np.where(np.isfinite(means.to_numpy(dtype=np.float64)),
                          means.to_numpy(dtype=np.float64), 0.0)

    written = []
    for src in lineages:
        col_mask = np.array([c1 == src for c1, c2 in means.columns])
        if not col_mask.any():
            continue
        m_src = means.loc[:, col_mask]
        p_src = pvals.loc[:, col_mask]
        means_src = means_arr[:, col_mask]
        pvals_src = pvals_clean[:, col_mask]

        sig_any = ((pvals_src < PVAL_THRESH) & (means_src >= MIN_MEAN)).any(axis=1)
        if sig_any.sum() == 0:
            print(f"[{sample_id}] {src}: no significant L-R pairs — skipping")
            continue
        m_src = m_src.loc[sig_any]
        p_src = p_src.loc[sig_any]
        score = (-np.log10(pvals_src[sig_any])).max(axis=1)
        order = np.argsort(-score)[:TOP_N_PER_SOURCE]
        m_top = m_src.iloc[order].copy()
        p_top = p_src.iloc[order].copy()

        res = {"means": m_top, "pvalues": p_top}
        if metadata is not None:
            res["metadata"] = metadata.reindex(m_top.index).copy()

        n_lr = m_top.shape[0]
        n_cp = m_top.shape[1]
        fig_w = max(6.0, 0.55 * n_cp + 4.0)
        fig_h = max(5.0, 0.18 * n_lr + 2.0)

        sq.pl.ligrec(
            res,
            means_range=(MIN_MEAN, np.inf),
            pvalue_threshold=PVAL_THRESH,
            alpha=ALPHA,
            remove_empty_interactions=True,
            swap_axes=True,
            title=f"{sample_id} — {src} → all targets  (top {n_lr} L-R interactions)",
            figsize=(fig_w, fig_h),
            dpi=DPI,
        )
        out = fig_dir / f"ligand_receptor_dotplot_{src}.png"
        plt.savefig(out, dpi=DPI, bbox_inches="tight")
        plt.close("all")
        size_mb = out.stat().st_size / 1e6
        print(f"[{sample_id}] {src}: wrote {out.name}  ({n_lr} pairs × {n_cp} targets, {size_mb:.2f} MB)")
        written.append(out)

    # Remove the legacy monolithic file if present so users don't open the
    # superseded subtype-level plot. Keep the per-source PNGs as the new index.
    legacy = fig_dir / "ligand_receptor_dotplot.png"
    if legacy.exists():
        legacy.unlink()
        print(f"[{sample_id}] removed legacy {legacy.name}")
    return written


def render_one(sample_id: str) -> None:
    adata = _build_lineage_subset(sample_id)

    print(f"[{sample_id}] running sq.gr.ligrec on celltype_general (n_perms={N_PERMS})")
    sq.gr.ligrec(
        adata,
        cluster_key="celltype_general",
        n_perms=N_PERMS,
        copy=False,
        use_raw=False,
    )
    lr = adata.uns["celltype_general_ligrec"]
    print(f"[{sample_id}]   means shape: {lr['means'].shape}, "
          f"pvalues shape: {lr['pvalues'].shape}")
    _save_csvs(sample_id, lr)

    lineages = list(adata.obs["celltype_general"].cat.categories)
    _render_per_source(sample_id, lr, lineages)


def main() -> None:
    for s in SAMPLES:
        render_one(s)


if __name__ == "__main__":
    main()
