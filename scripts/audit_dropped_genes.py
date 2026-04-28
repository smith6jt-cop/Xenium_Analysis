"""Compare raw zarr panel vs preprocessed h5ad to identify dropped genes.

Stage 01 cell 14 drops genes via two steps:
1. control-probe regex (NegControlProbe_, UnassignedCodeword_,
   NegControlCodeword_, antisense_, BLANK_)
2. sc.pp.filter_genes(min_cells=100)

This script identifies which genes fell into which bucket per sample.
"""
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import spatialdata as sd

CONTROL_RE = re.compile(
    r"^(NegControlProbe_|UnassignedCodeword_|NegControlCodeword_|antisense_|BLANK_)",
    re.IGNORECASE,
)


def main(sample):
    print(f"\n{'='*72}\n=== {sample} ===\n{'='*72}")
    zarr = Path(f"data/raw/{sample}.zarr")
    pre = Path(f"data/processed/{sample}/{sample}_preprocessed.h5ad")

    sd_data = sd.read_zarr(zarr)
    raw_genes = list(sd_data.tables["table"].var_names)
    print(f"  raw zarr panel: {len(raw_genes):,} features")

    # Categorize raw features
    control = [g for g in raw_genes if CONTROL_RE.match(g)]
    real = [g for g in raw_genes if not CONTROL_RE.match(g)]
    print(f"  control probes: {len(control):,}")
    print(f"  real genes:     {len(real):,}")

    # Buckets among controls
    buckets = {}
    for g in control:
        prefix = re.match(r"^([A-Za-z]+_)", g)
        key = prefix.group(1) if prefix else "other"
        buckets[key] = buckets.get(key, 0) + 1
    if buckets:
        print(f"  control breakdown:")
        for k, v in sorted(buckets.items()):
            print(f"    {k:<24} {v:>4}")

    # Compare against preprocessed
    with h5py.File(pre, "r") as f:
        if "_index" in f["var"]:
            kept = list(f["var/_index"].asstr()[:])
        else:
            attr = f["var"].attrs.get("_index", "_index")
            kept = list(f[f"var/{attr}"].asstr()[:])
    kept_set = set(kept)
    real_set = set(real)
    print(f"\n  preprocessed kept genes: {len(kept):,}")

    # Genes dropped = raw - kept
    dropped = sorted(set(raw_genes) - kept_set)
    print(f"  dropped genes total: {len(dropped):,}")

    # Separate into control vs filter_genes(min_cells)
    dropped_control = [g for g in dropped if CONTROL_RE.match(g)]
    dropped_low_expr = [g for g in dropped if not CONTROL_RE.match(g)]
    print(f"    dropped as control probes:           {len(dropped_control):,}")
    print(f"    dropped by min_cells=100 (low expr): {len(dropped_low_expr):,}")

    if dropped_low_expr:
        # Look up n_cells_by_counts for the dropped genes from raw counts
        adata_raw = sd_data.tables["table"]
        # We need n_cells per gene in the FULL raw panel
        # Since adata_raw.X is the raw counts, count cells where each gene > 0
        X = adata_raw.X
        if hasattr(X, "toarray"):
            # sparse — use nnz per column
            from scipy.sparse import csc_matrix
            Xc = X.tocsc()
            n_cells_per_gene = np.diff(Xc.indptr)
        else:
            n_cells_per_gene = (np.asarray(X) > 0).sum(axis=0)
        gene_to_n = dict(zip(raw_genes, n_cells_per_gene))

        # Sort dropped low-expr by n_cells ascending
        dropped_with_n = sorted(
            [(g, int(gene_to_n.get(g, -1))) for g in dropped_low_expr],
            key=lambda x: x[1],
        )
        print(f"\n  dropped low-expr genes (sorted by n_cells; threshold = 100):")
        for g, n in dropped_with_n[:50]:
            print(f"    {g:<22} n_cells={n:>5}")
        if len(dropped_with_n) > 50:
            print(f"    ... {len(dropped_with_n) - 50} more")
        # Stats
        ns = [n for _, n in dropped_with_n]
        print(f"\n  dropped low-expr gene n_cells distribution:")
        print(f"    min={min(ns)}  med={int(np.median(ns))}  max={max(ns)}")
        print(f"    fraction with n_cells == 0: "
              f"{sum(1 for n in ns if n == 0)}/{len(ns)}")

    if dropped_control:
        print(f"\n  dropped control probe sample (first 20):")
        for g in dropped_control[:20]:
            print(f"    {g}")
        if len(dropped_control) > 20:
            print(f"    ... {len(dropped_control) - 20} more")


for s in ("0041323", "0041326"):
    main(s)
