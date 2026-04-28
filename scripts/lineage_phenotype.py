"""Lineage-phenotyping post-hoc on existing _phenotyped.h5ad / _spatial_analysis.h5ad.

For multi-lineage cells (>=2 mutually-exclusive groups detected at raw count >= 3),
count per-canary-lineage marker depth and resolve to the dominant lineage when its
depth exceeds the runner-up by a clear margin. Cells resolved this way are tagged
"lineage_phenotyped". Single-lineage cells are tagged "single_lineage".
3+ group cells are tagged "doublet_suspected" (kept as-is, no resolution).

New obs columns added (mirrored to _spatial_analysis.h5ad):
    lineage_status     : single_lineage | lineage_phenotyped | lineage_ambiguous | doublet_suspected
    celltype_lineage   : refined celltype (or copy of celltype if not phenotyped)
    lineage_dominant   : top-depth canary lineage among multi-lineage cells (NaN otherwise)
    lineage_top_depth  : n markers detected >= 3 in dominant lineage
    lineage_second_depth: n markers detected >= 3 in runner-up lineage
    lineage_n_groups   : n distinct mutually-exclusive groups detected

Reuses CANARY_PANELS and EXCLUSIVE_GROUP from
scripts/verify_upper_tail_doublets.py (after Acinar panel expanded to
[AMY1A, CUZD1]).
"""
import shutil
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, "scripts")
from verify_upper_tail_doublets import CANARY_PANELS, EXCLUSIVE_GROUP

ROOT = Path("/blue/maigan/smith6jt/Xenium_Analysis")
SAMPLES = ("0041323", "0041326")

# Detection threshold on the LOGNORM layer (log1p of normalize_total).
# normalize_total rescales each cell to the sample's median library size, so
# the same lognorm threshold maps to comparable normalized expression across
# samples regardless of overall sequencing depth. log(1+3) ≈ 1.386 preserves
# the prior "raw count >= 3" semantics for a median-depth cell.
LOG_DETECT = np.log1p(3.0)  # ≈ 1.386
DOMINANCE_MARGIN = 2  # top_depth - second_depth >= 2 to call "phenotyped"


def compute_lineage_columns(adata):
    """Return a DataFrame indexed by adata.obs_names with the new columns."""
    var_names = list(adata.var_names)
    gene_idx = {g: i for i, g in enumerate(var_names)}

    # Build dense canary submatrix
    all_canary = sorted({g for genes in CANARY_PANELS.values() for g in genes})
    present = [g for g in all_canary if g in gene_idx]
    cols = [gene_idx[g] for g in present]
    pos = {g: j for j, g in enumerate(present)}

    # Use lognorm layer for sample-comparable detection threshold.
    if "lognorm" in adata.layers:
        layer = adata.layers["lognorm"]
        layer_used = "lognorm"
    else:
        # Fallback for older h5ads: build on the fly from counts.
        print("  WARN: no 'lognorm' layer found, computing from counts on the fly")
        adata_tmp = adata.copy()
        adata_tmp.X = adata.layers["counts"].copy()
        sc.pp.normalize_total(adata_tmp)
        sc.pp.log1p(adata_tmp)
        layer = adata_tmp.X
        layer_used = "lognorm (computed)"
    print(f"  using detection layer: {layer_used} (threshold lognorm >= {LOG_DETECT:.3f})", flush=True)

    print(f"  building canary_dense {adata.n_obs:,} x {len(present)} ...", flush=True)
    t = time.time()
    if hasattr(layer, "toarray"):
        canary_dense = np.asarray(layer[:, cols].toarray())
    else:
        canary_dense = np.asarray(layer[:, cols])
    print(f"    {time.time()-t:.1f}s", flush=True)

    # Per-lineage depth: n markers in panel with lognorm >= LOG_DETECT
    print(f"  computing per-lineage depths...", flush=True)
    t = time.time()
    lineages = list(CANARY_PANELS.keys())
    depth_arr = np.zeros((adata.n_obs, len(lineages)), dtype=np.int8)
    sum_arr = np.zeros((adata.n_obs, len(lineages)), dtype=np.float32)
    for li, lin in enumerate(lineages):
        idxs = [pos[g] for g in CANARY_PANELS[lin] if g in pos]
        if not idxs:
            continue
        sub = canary_dense[:, idxs]
        depth_arr[:, li] = (sub >= LOG_DETECT).sum(axis=1).astype(np.int8)
        sum_arr[:, li] = sub.sum(axis=1)
    print(f"    {time.time()-t:.1f}s", flush=True)

    # n distinct exclusive groups detected (lineage detected if depth >= 1)
    print(f"  computing exclusive-group counts...", flush=True)
    t = time.time()
    detected = depth_arr >= 1  # (n_cells, n_lineages)
    # Map lineage -> group, group_idx
    group_names = sorted(set(EXCLUSIVE_GROUP[lin] for lin in lineages))
    group_to_idx = {g: i for i, g in enumerate(group_names)}
    lin_group_idx = np.array([group_to_idx[EXCLUSIVE_GROUP[lin]] for lin in lineages])
    # For each cell, mark which groups have any detected lineage
    group_hit = np.zeros((adata.n_obs, len(group_names)), dtype=bool)
    for li in range(len(lineages)):
        if not detected[:, li].any():
            continue
        gi = lin_group_idx[li]
        group_hit[:, gi] |= detected[:, li]
    n_groups = group_hit.sum(axis=1).astype(np.int8)
    print(f"    {time.time()-t:.1f}s", flush=True)

    # Identify dominant lineage per cell (argmax of depth, tie-break by sum)
    # Compose a (depth, sum) sort key using depth*1e6 + sum (sum < 1e6 for our scale)
    # composite = depth (primary) + scaled sum (tiebreak). depth maxes at ~7;
    # sum is lognorm sum (typically <50 per lineage). Scale sum to keep depth dominant.
    composite = depth_arr.astype(np.float64) * 1000.0 + sum_arr.astype(np.float64)
    top_idx = composite.argmax(axis=1)
    top_depth = depth_arr[np.arange(adata.n_obs), top_idx]

    # Second-best depth: zero out top, take next argmax depth
    composite2 = composite.copy()
    composite2[np.arange(adata.n_obs), top_idx] = -1
    second_idx = composite2.argmax(axis=1)
    second_depth = depth_arr[np.arange(adata.n_obs), second_idx]

    margin = top_depth.astype(np.int16) - second_depth.astype(np.int16)
    # Status assignment
    status = np.empty(adata.n_obs, dtype=object)
    is_doublet = n_groups >= 3
    is_multi = (n_groups == 2) & (~is_doublet)
    is_single = n_groups <= 1

    status[is_single] = "single_lineage"
    status[is_doublet] = "doublet_suspected"
    # Multi-lineage: phenotyped if margin>=2 AND top_depth>=2; else ambiguous
    can_resolve = is_multi & (margin >= DOMINANCE_MARGIN) & (top_depth >= 2)
    status[is_multi & can_resolve] = "lineage_phenotyped"
    status[is_multi & ~can_resolve] = "lineage_ambiguous"

    # Build celltype_lineage
    original_celltype = adata.obs["celltype"].astype(str).to_numpy()
    celltype_lineage = original_celltype.copy()
    # Only overwrite for lineage_phenotyped cells
    resolved_mask = status == "lineage_phenotyped"
    lineage_names_arr = np.array(lineages)
    celltype_lineage[resolved_mask] = lineage_names_arr[top_idx[resolved_mask]]

    # lineage_dominant column: only meaningful for non-single cells
    dominant = np.empty(adata.n_obs, dtype=object)
    dominant[:] = ""
    dominant[~is_single] = lineage_names_arr[top_idx[~is_single]]

    df = pd.DataFrame(
        {
            "lineage_status": pd.Categorical(
                status,
                categories=["single_lineage", "lineage_phenotyped",
                            "lineage_ambiguous", "doublet_suspected"],
            ),
            "celltype_lineage": pd.Categorical(celltype_lineage),
            "lineage_dominant": pd.Categorical(dominant),
            "lineage_top_depth": top_depth.astype(np.int8),
            "lineage_second_depth": second_depth.astype(np.int8),
            "lineage_n_groups": n_groups,
        },
        index=adata.obs_names.copy(),
    )
    return df


def apply_to_h5ad(h5ad_path, columns_df, backup=True):
    """Add the new obs columns to an existing h5ad and re-write."""
    h5ad_path = Path(h5ad_path)
    if not h5ad_path.exists():
        print(f"  SKIP (missing): {h5ad_path}")
        return

    if backup:
        bak = h5ad_path.with_suffix(h5ad_path.suffix + ".bak.prelineage")
        if not bak.exists():
            print(f"  backing up to {bak.name} ...", flush=True)
            t = time.time()
            shutil.copy2(h5ad_path, bak)
            print(f"    {time.time()-t:.1f}s", flush=True)

    print(f"  loading {h5ad_path.name} ...", flush=True)
    t = time.time()
    a = sc.read_h5ad(h5ad_path)
    print(f"    {a.shape} in {time.time()-t:.1f}s", flush=True)

    # Align by obs_names. Spatial h5ad is slimmed but should have the same n_obs.
    if a.n_obs != len(columns_df):
        print(f"  WARN n_obs mismatch: file has {a.n_obs:,} vs computed {len(columns_df):,}; "
              f"reindexing on obs_names")
    aligned = columns_df.reindex(a.obs_names)
    if aligned.isna().any().any():
        n_na = aligned.isna().any(axis=1).sum()
        print(f"  WARN {n_na:,} cells missing in computed table — falling back to single_lineage")
        # fill NaNs conservatively
        aligned["lineage_status"] = aligned["lineage_status"].cat.add_categories([]).fillna(
            "single_lineage"
        )
        aligned["celltype_lineage"] = aligned["celltype_lineage"].fillna(
            a.obs["celltype"].astype(str)
        )
        aligned["lineage_dominant"] = aligned["lineage_dominant"].fillna("")
        aligned["lineage_top_depth"] = aligned["lineage_top_depth"].fillna(0).astype(np.int8)
        aligned["lineage_second_depth"] = aligned["lineage_second_depth"].fillna(0).astype(np.int8)
        aligned["lineage_n_groups"] = aligned["lineage_n_groups"].fillna(0).astype(np.int8)

    for col in aligned.columns:
        a.obs[col] = aligned[col].values

    print(f"  writing {h5ad_path.name} ...", flush=True)
    t = time.time()
    a.write_h5ad(h5ad_path)
    print(f"    {time.time()-t:.1f}s", flush=True)


def main():
    for sample in SAMPLES:
        print(f"\n{'='*72}\n=== {sample} ===\n{'='*72}")
        pheno_path = ROOT / f"data/processed/{sample}/{sample}_phenotyped.h5ad"
        spatial_path = ROOT / f"data/processed/{sample}/{sample}_spatial_analysis.h5ad"

        # Load phenotyped (has counts layer)
        print(f"\n  loading {pheno_path.name} (for counts layer) ...", flush=True)
        t = time.time()
        a = sc.read_h5ad(pheno_path)
        print(f"    {a.shape} in {time.time()-t:.1f}s", flush=True)

        df = compute_lineage_columns(a)

        print(f"\n  lineage_status distribution:")
        print(df["lineage_status"].value_counts().to_string())
        diff_mask = a.obs["celltype"].astype(str).to_numpy() != df["celltype_lineage"].astype(str).to_numpy()
        print(f"\n  cells reclassified (celltype -> celltype_lineage): {int(diff_mask.sum()):,}")
        if diff_mask.any():
            tab = (
                pd.DataFrame({
                    "celltype": a.obs["celltype"].astype(str).to_numpy()[diff_mask],
                    "celltype_lineage": df["celltype_lineage"].astype(str).to_numpy()[diff_mask],
                })
                .groupby(["celltype", "celltype_lineage"], observed=True)
                .size()
                .sort_values(ascending=False)
            )
            print(tab.head(25).to_string())

        # Free memory before re-writing (the load was just for compute)
        del a

        # Apply to phenotyped (re-writes the file with new columns added)
        print(f"\n  --> updating {pheno_path.name}")
        apply_to_h5ad(pheno_path, df)

        # Mirror to spatial_analysis
        print(f"\n  --> updating {spatial_path.name}")
        apply_to_h5ad(spatial_path, df)


if __name__ == "__main__":
    main()
    print("\nDONE")
