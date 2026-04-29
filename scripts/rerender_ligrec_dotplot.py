"""Re-render figures/{sample}/03_spatial_analysis/ligand_receptor_dotplot.png
from the saved {sample}_ligrec_{pvalues,means,metadata}.csv exports.

The saved _spatial_analysis.h5ad files no longer carry adata.uns['celltype_ligrec']
(stage 03 cell-18 strips it after CSV export), so we reconstruct the dict from
the three CSVs and pass it directly to sq.pl.ligrec, which accepts a
`Mapping[str, pd.DataFrame]` per its type signature.

Usage:
    python scripts/rerender_ligrec_dotplot.py                 # both samples
    python scripts/rerender_ligrec_dotplot.py 0041323         # one sample
"""

from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import squidpy as sq

PROJECT = Path("/blue/maigan/smith6jt/Xenium_Analysis")
DEFAULT_SAMPLES = ["0041323", "0041326"]


def load_ligrec_dict(sample: str) -> dict[str, pd.DataFrame]:
    proc = PROJECT / f"data/processed/{sample}"
    pv = pd.read_csv(proc / f"{sample}_ligrec_pvalues.csv",
                     header=[0, 1], index_col=[0, 1])
    mn = pd.read_csv(proc / f"{sample}_ligrec_means.csv",
                     header=[0, 1], index_col=[0, 1])
    md = pd.read_csv(proc / f"{sample}_ligrec_metadata.csv",
                     index_col=[0, 1])
    # Ensure means/pvalues share the row index of metadata (sq.pl.ligrec
    # joins on the (source, target) MultiIndex).
    pv.index.names = ["source", "target"]
    mn.index.names = ["source", "target"]
    md.index.names = ["source", "target"]
    pv.columns.names = ["cluster_1", "cluster_2"]
    mn.columns.names = ["cluster_1", "cluster_2"]
    return {"pvalues": pv, "means": mn, "metadata": md}


def render(sample: str) -> None:
    fig_dir = PROJECT / f"figures/{sample}/03_spatial_analysis"
    fig_dir.mkdir(parents=True, exist_ok=True)
    out_png = fig_dir / "ligand_receptor_dotplot.png"

    print(f"[{sample}] loading ligrec CSVs")
    d = load_ligrec_dict(sample)
    pv, mn = d["pvalues"], d["means"]
    print(f"[{sample}] pvalues shape={pv.shape}  means shape={mn.shape}")
    sig = (pv.values <= 0.001) & (mn.values >= 0.3)
    print(f"[{sample}] hits at means>=0.3 & p<=0.001: {int(sig.sum())} "
          f"(rows w/ ≥1 hit: {int(sig.any(axis=1).sum())})")

    sq.pl.ligrec(
        d,
        means_range=(0.3, np.inf),
        pvalue_threshold=0.001,
        alpha=1e-4,
        remove_empty_interactions=True,
        swap_axes=True,
        title="Significant ligand-receptor interactions",
    )
    fig = plt.gcf()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close("all")
    print(f"[{sample}] wrote {out_png}")


def main():
    samples = sys.argv[1:] if len(sys.argv) > 1 else DEFAULT_SAMPLES
    for s in samples:
        render(s)


if __name__ == "__main__":
    main()
