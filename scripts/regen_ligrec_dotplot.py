"""Regenerate the ligand-receptor dotplot from the saved *_ligrec_*.csv.

Uses native sq.pl.ligrec (two-axis ligand→receptor labels, the squidpy
canonical style) but pre-filters the means/pvalues/metadata dict so the
output figure is bounded:

  1. Drop Indeterminate* clusters and clusters with < 50 cells.
  2. Filter L-R pairs to those with at least one (mean ≥ MIN_MEAN, p < PVAL).
  3. Keep the TOP_N L-R pairs by max(-log10(p)) across cluster pairs.
  4. Pass the filtered dict to sq.pl.ligrec with a bounded figsize + dpi.

Run after stage 03 finishes; doesn't require re-executing sq.gr.ligrec.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
SAMPLES = ("0041323", "0041326")

MIN_CELLS_PER_CLUSTER = 50
TOP_N_LR_PAIRS = 150
MIN_MEAN = 0.3
PVAL_THRESH = 0.001
ALPHA = 1e-4
P_FLOOR = 1e-12
DPI = 120


def _resolve_celltype_counts(sample_id: str) -> pd.Series | None:
    import scanpy as sc
    h5 = REPO / "data" / "processed" / sample_id / f"{sample_id}_spatial_analysis.h5ad"
    if not h5.exists():
        return None
    a = sc.read_h5ad(h5, backed="r")
    counts = a.obs["celltype"].value_counts()
    a.file.close()
    return counts


def render_one(sample_id: str) -> None:
    import squidpy as sq

    proc = REPO / "data" / "processed" / sample_id
    fig_dir = REPO / "figures" / sample_id / "03_spatial_analysis"
    means_csv = proc / f"{sample_id}_ligrec_means.csv"
    pvals_csv = proc / f"{sample_id}_ligrec_pvalues.csv"
    meta_csv = proc / f"{sample_id}_ligrec_metadata.csv"
    out = fig_dir / "ligand_receptor_dotplot.png"

    if not means_csv.exists() or not pvals_csv.exists():
        print(f"[{sample_id}] missing ligrec CSVs in {proc} — skipping")
        return
    fig_dir.mkdir(parents=True, exist_ok=True)

    print(f"[{sample_id}] loading {means_csv.name}")
    means = pd.read_csv(means_csv, header=[0, 1], index_col=[0, 1])
    print(f"[{sample_id}] loading {pvals_csv.name}")
    pvals = pd.read_csv(pvals_csv, header=[0, 1], index_col=[0, 1])
    metadata = None
    if meta_csv.exists():
        # metadata has 1-level columns and MultiIndex rows matching means.
        metadata = pd.read_csv(meta_csv, header=0, index_col=[0, 1])
        # Align to means index in case ordering differs.
        metadata = metadata.reindex(means.index)
    print(f"[{sample_id}] raw shape: means={means.shape}  pvalues={pvals.shape}")

    # --- 1) Drop Indeterminate* and tiny clusters from BOTH cluster axes. ---
    counts = _resolve_celltype_counts(sample_id)
    src_cats = sorted(set(means.columns.get_level_values(0)))
    keep = [c for c in src_cats
            if not c.startswith("Indeterminate")
            and (counts is None or counts.get(c, 0) >= MIN_CELLS_PER_CLUSTER)]
    if not keep:
        print(f"[{sample_id}] no clusters survive filter — skipping")
        return
    print(f"[{sample_id}] keeping {len(keep)} clusters: {keep}")

    keep_set = set(keep)
    col_mask = np.array([(c1 in keep_set) and (c2 in keep_set)
                         for c1, c2 in means.columns])
    means_f = means.loc[:, col_mask]
    pvals_f = pvals.loc[:, col_mask]
    print(f"[{sample_id}] cluster-pairs after filter: {means_f.shape[1]}")

    # --- 2) NaN/zero handling, then filter rows to those with >=1 sig hit. ---
    means_arr = np.where(np.isfinite(means_f.to_numpy(dtype=np.float64)),
                          means_f.to_numpy(dtype=np.float64), 0.0)
    pvals_clean = np.where(np.isfinite(pvals_f.to_numpy(dtype=np.float64)),
                            np.maximum(pvals_f.to_numpy(dtype=np.float64), P_FLOOR),
                            1.0)
    sig_any = ((pvals_clean < PVAL_THRESH) & (means_arr >= MIN_MEAN)).any(axis=1)
    means_f = means_f.loc[sig_any]
    pvals_f = pvals_f.loc[sig_any]
    pvals_clean = pvals_clean[sig_any]
    print(f"[{sample_id}] L-R pairs surviving threshold: {means_f.shape[0]}")

    if means_f.shape[0] == 0:
        print(f"[{sample_id}] no L-R pairs survive — relax thresholds")
        return

    # --- 3) Top-N L-R pairs by max(-log10(p)) across cluster-pairs. ---
    score = (-np.log10(pvals_clean)).max(axis=1)
    order = np.argsort(-score)[:TOP_N_LR_PAIRS]
    means_top = means_f.iloc[order].copy()
    pvals_top = pvals_f.iloc[order].copy()
    print(f"[{sample_id}] keeping top {means_top.shape[0]} L-R pairs "
          f"(score range {score[order].min():.2f}..{score[order].max():.2f})")

    # --- 4) Build the dict squidpy expects and call sq.pl.ligrec. ---
    res = {"means": means_top, "pvalues": pvals_top}
    if metadata is not None:
        meta_top = metadata.reindex(means_top.index).copy()
        res["metadata"] = meta_top

    n_lr = means_top.shape[0]
    n_cp = means_top.shape[1]
    # swap_axes=True puts cluster-pairs on x and L-R pairs on y (the
    # original layout the user preferred). Bounded figsize.
    fig_w = min(max(8.0, 0.18 * n_cp + 4.0), 22.0)
    fig_h = min(max(8.0, 0.18 * n_lr + 2.0), 30.0)

    sq.pl.ligrec(
        res,
        means_range=(MIN_MEAN, np.inf),
        pvalue_threshold=PVAL_THRESH,
        alpha=ALPHA,
        remove_empty_interactions=True,
        swap_axes=True,
        title=f"{sample_id} — top {n_lr} L-R interactions  "
              f"({len(keep)} celltypes, p<{PVAL_THRESH}, mean≥{MIN_MEAN})",
        figsize=(fig_w, fig_h),
        dpi=DPI,
    )
    plt.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close("all")
    print(f"[{sample_id}] wrote {out}  ({out.stat().st_size / 1e6:.2f} MB)")


def main() -> None:
    for s in SAMPLES:
        render_one(s)


if __name__ == "__main__":
    main()
