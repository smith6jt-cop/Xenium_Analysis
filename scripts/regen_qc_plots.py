"""Re-render stage-01 QC plots (qc_histo_*, qc_metrics_*) for both samples.

The shipped files were rendered before/after with two separate violin+scatter
calls that auto-scaled axes, which made the ~7%-cell delta invisible. New
design overlays before+after on the same axes with explicit cell counts and
color coding (grey before, blue after).

Loads raw zarr to recover the pre-filter distribution (the saved
`_preprocessed.h5ad` only has the post-filter cells).
"""
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import spatialdata as sd

ROOT = Path("/blue/maigan/smith6jt/Xenium_Analysis")
SAMPLES = ("0041323", "0041326")

# Filter params (match stage 01 cell 14)
MIN_GENES = 20
MIN_COUNTS = 50
MAX_COUNTS_QUANTILE = 0.98
MIN_CELLS = 100


def load_and_qc(sample: str):
    zarr_path = ROOT / f"data/raw/{sample}.zarr"
    print(f"  loading {zarr_path} ...")
    t = time.time()
    sdata = sd.read_zarr(zarr_path)
    adata = sdata.tables["table"]
    print(f"    raw: {adata.shape} ({time.time()-t:.1f}s)")

    control_mask = adata.var_names.str.match(
        r"^(NegControlProbe_|UnassignedCodeword_|NegControlCodeword_|antisense_|BLANK_)",
        case=False,
    )
    if control_mask.any():
        adata = adata[:, ~control_mask].copy()
        print(f"    dropped {int(control_mask.sum())} control probes")

    sc.pp.calculate_qc_metrics(adata, percent_top=(10, 20, 50, 150), inplace=True)
    return adata


def apply_filters(adata):
    sc.pp.filter_cells(adata, min_genes=MIN_GENES)
    sc.pp.filter_cells(adata, min_counts=MIN_COUNTS)
    density = adata.obs["total_counts"].astype(np.float64) / np.maximum(
        adata.obs["cell_area"].astype(np.float64), 1.0
    )
    cutoff = float(np.quantile(density, MAX_COUNTS_QUANTILE))
    adata = adata[(density <= cutoff).to_numpy()].copy()
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS)
    return adata, cutoff


def overlay_histograms(pre_obs, post_obs, fig_path, title_suffix=""):
    n_before = len(pre_obs)
    n_after = len(post_obs)
    n_dropped = n_before - n_after

    fig, axs = plt.subplots(1, 4, figsize=(25, 6))

    def panel(ax, before, after, title):
        sns.histplot(
            before, kde=False, ax=ax, color="#888888", alpha=0.45,
            element="step", linewidth=1.5, stat="count",
            label=f"before ({n_before:,})",
        )
        sns.histplot(
            after, kde=False, ax=ax, color="#1f77b4", alpha=0.75,
            element="step", linewidth=1.5, stat="count",
            label=f"after ({n_after:,})",
        )
        ax.set_title(title)
        ax.legend(loc="upper right")

    panel(axs[0], pre_obs["total_counts"], post_obs["total_counts"],
          "Total transcripts per cell")
    panel(axs[1], pre_obs["n_genes_by_counts"], post_obs["n_genes_by_counts"],
          "Unique transcripts per cell")
    panel(axs[2], pre_obs["cell_area"], post_obs["cell_area"],
          "Area of segmented cells")
    panel(axs[3],
          pre_obs["nucleus_area"] / pre_obs["cell_area"],
          post_obs["nucleus_area"] / post_obs["cell_area"],
          "Nucleus ratio")

    fig.suptitle(
        f"QC histograms: dropped {n_dropped:,} cells "
        f"({100.0*n_dropped/max(n_before,1):.1f}%) "
        f"by min/max filters + density q98{title_suffix}",
        y=1.02,
    )
    plt.tight_layout()
    plt.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"    wrote {fig_path}")


def overlay_metrics(pre_obs, post_obs, fig_path, kept_idx_set):
    n_before = len(pre_obs)
    n_after = len(post_obs)

    fig, axes = plt.subplots(1, 3, figsize=(25, 10))

    # Violin: total_counts before vs after
    df_tc = pd.concat([
        pd.DataFrame({"value": pre_obs["total_counts"].values, "phase": "before"}),
        pd.DataFrame({"value": post_obs["total_counts"].values, "phase": "after"}),
    ], ignore_index=True)
    sns.violinplot(
        data=df_tc, x="phase", y="value", ax=axes[0], inner="quartile",
        palette={"before": "#888888", "after": "#1f77b4"}, cut=0,
    )
    axes[0].set_title(f"Total counts per cell\nbefore={n_before:,}  after={n_after:,}")

    # Violin: n_genes before vs after
    df_ng = pd.concat([
        pd.DataFrame({"value": pre_obs["n_genes_by_counts"].values, "phase": "before"}),
        pd.DataFrame({"value": post_obs["n_genes_by_counts"].values, "phase": "after"}),
    ], ignore_index=True)
    sns.violinplot(
        data=df_ng, x="phase", y="value", ax=axes[1], inner="quartile",
        palette={"before": "#888888", "after": "#1f77b4"}, cut=0,
    )
    axes[1].set_title(f"Genes per cell\nbefore={n_before:,}  after={n_after:,}")

    # Scatter: pre-filter cells, kept (blue) vs dropped (red)
    rng = np.random.default_rng(0)
    idx = rng.choice(n_before, size=min(150_000, n_before), replace=False)
    sub = pre_obs.iloc[idx]
    sub_kept = np.array([str(c) in kept_idx_set for c in sub.index])
    axes[2].scatter(
        sub["total_counts"][sub_kept], sub["n_genes_by_counts"][sub_kept],
        s=3, c="#1f77b4", alpha=0.4, linewidths=0,
        label=f"kept ({sub_kept.sum():,})",
    )
    axes[2].scatter(
        sub["total_counts"][~sub_kept], sub["n_genes_by_counts"][~sub_kept],
        s=4, c="#d62728", alpha=0.7, linewidths=0,
        label=f"dropped ({(~sub_kept).sum():,})",
    )
    axes[2].set_xlabel("total_counts")
    axes[2].set_ylabel("n_genes_by_counts")
    axes[2].set_title(f"Counts vs genes (subsample {len(idx):,} cells)")
    axes[2].legend(loc="best", framealpha=0.9)

    plt.tight_layout()
    plt.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"    wrote {fig_path}")


def single_distribution_histograms(obs, fig_path):
    """Match the new cell-11 design — single 'before' distribution per panel."""
    n = len(obs)
    fig, axs = plt.subplots(1, 4, figsize=(25, 6))
    axs[0].set_title(f"Total transcripts per cell (n={n:,})")
    sns.histplot(obs["total_counts"], kde=False, ax=axs[0])
    axs[1].set_title(f"Unique transcripts per cell (n={n:,})")
    sns.histplot(obs["n_genes_by_counts"], kde=False, ax=axs[1])
    axs[2].set_title(f"Area of segmented cells (n={n:,})")
    sns.histplot(obs["cell_area"], kde=False, ax=axs[2])
    axs[3].set_title(f"Nucleus ratio (n={n:,})")
    sns.histplot(obs["nucleus_area"] / obs["cell_area"], kde=False, ax=axs[3])
    plt.tight_layout()
    plt.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"    wrote {fig_path}")


def single_metrics_violin(adata_obs, fig_path):
    """Match the original 'before' violin/scatter layout."""
    fig, axes = plt.subplots(1, 3, figsize=(25, 10))
    sns.violinplot(y=adata_obs["total_counts"].values, ax=axes[0], color="#888888", cut=0)
    axes[0].set_title(f"Total counts per cell (n={len(adata_obs):,})")
    sns.violinplot(y=adata_obs["n_genes_by_counts"].values, ax=axes[1], color="#888888", cut=0)
    axes[1].set_title(f"Genes per cell (n={len(adata_obs):,})")
    rng = np.random.default_rng(0)
    n = len(adata_obs)
    idx = rng.choice(n, size=min(150_000, n), replace=False)
    sub = adata_obs.iloc[idx]
    axes[2].scatter(sub["total_counts"], sub["n_genes_by_counts"],
                    s=3, c="#1f77b4", alpha=0.4, linewidths=0)
    axes[2].set_xlabel("total_counts"); axes[2].set_ylabel("n_genes_by_counts")
    axes[2].set_title(f"Counts vs genes (subsample {len(idx):,} cells)")
    plt.tight_layout()
    plt.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"    wrote {fig_path}")


for sample in SAMPLES:
    print(f"\n=== {sample} ===")
    fig_dir = ROOT / f"figures/{sample}/01_preprocessing"
    adata_pre = load_and_qc(sample)
    pre_obs = adata_pre.obs[
        ["total_counts", "n_genes_by_counts", "cell_area", "nucleus_area"]
    ].copy()
    n_before = len(pre_obs)

    adata_post, cutoff = apply_filters(adata_pre.copy())
    post_obs = adata_post.obs[
        ["total_counts", "n_genes_by_counts", "cell_area", "nucleus_area"]
    ].copy()
    n_after = len(post_obs)
    print(f"  filter: {n_before:,} -> {n_after:,} (cutoff = {cutoff:.3f} counts/area)")

    kept_idx_set = set(post_obs.index.astype(str))

    # Before-only plots (single distribution)
    single_distribution_histograms(pre_obs, fig_dir / "qc_histo_before_filtering.png")
    single_metrics_violin(pre_obs, fig_dir / "qc_metrics_before_filtering.png")

    # After plots are now overlays (before vs after)
    overlay_histograms(pre_obs, post_obs, fig_dir / "qc_histo_after_filtering.png")
    overlay_metrics(pre_obs, post_obs, fig_dir / "qc_metrics_after_filtering.png", kept_idx_set)

    del adata_pre, adata_post

print("\nALL QC plots regenerated for both samples")
