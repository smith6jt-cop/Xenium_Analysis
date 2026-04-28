"""Rigorous verification of the density-based QC filter.

Steps:
1. Recompute pre-filter stats from raw zarr (drop control probes only).
2. Recompute filter cutoffs and apply min_genes / min_counts / density q98.
3. Compare to the saved _preprocessed.h5ad (must match cell count exactly).
4. Show n_genes distribution before/after filter; flag cells with n_genes > 1200.
5. Cross-reference high-n_genes cells with celltype labels from phenotyped.h5ad.
6. Sanity-check density filter: confirm the dropped cells are upper-tail density,
   not normal-density big cells.
"""
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import spatialdata as sd

ROOT = Path("/blue/maigan/smith6jt/Xenium_Analysis")
SAMPLES = ("0041323", "0041326")

MIN_GENES = 20
MIN_COUNTS = 50
MAX_COUNTS_QUANTILE = 0.98
MIN_CELLS = 100


def main(sample: str):
    print(f"\n{'='*70}\n=== {sample} ===\n{'='*70}")

    # Step 1: load raw zarr (the only source of pre-filter cells)
    zarr_path = ROOT / f"data/raw/{sample}.zarr"
    sdata = sd.read_zarr(zarr_path)
    adata = sdata.tables["table"]
    print(f"\nStep 1 — raw zarr: {adata.shape}")

    control_mask = adata.var_names.str.match(
        r"^(NegControlProbe_|UnassignedCodeword_|NegControlCodeword_|antisense_|BLANK_)",
        case=False,
    )
    if control_mask.any():
        adata = adata[:, ~control_mask].copy()
        print(f"  after dropping {int(control_mask.sum())} control probes: {adata.shape}")

    sc.pp.calculate_qc_metrics(adata, percent_top=(10, 20, 50, 150), inplace=True)

    n_pre = adata.n_obs
    n_genes_pre = adata.obs["n_genes_by_counts"].values.copy()
    total_counts_pre = adata.obs["total_counts"].values.copy()
    cell_area_pre = adata.obs["cell_area"].values.copy()
    cell_id_pre = adata.obs["cell_id"].astype(str).values.copy() if "cell_id" in adata.obs.columns else adata.obs.index.astype(str).values.copy()

    # Step 2: replicate filter logic exactly
    print(f"\nStep 2 — apply filter (min_genes={MIN_GENES}, min_counts={MIN_COUNTS}, "
          f"density>q{int(MAX_COUNTS_QUANTILE*100)})")
    n0 = adata.n_obs

    sc.pp.filter_cells(adata, min_genes=MIN_GENES)
    n1 = adata.n_obs
    print(f"  after min_genes={MIN_GENES}:    {n1:,}  (dropped {n0-n1:,})")

    sc.pp.filter_cells(adata, min_counts=MIN_COUNTS)
    n2 = adata.n_obs
    print(f"  after min_counts={MIN_COUNTS}:    {n2:,}  (dropped {n1-n2:,})")

    density = adata.obs["total_counts"].astype(np.float64).values / np.maximum(
        adata.obs["cell_area"].astype(np.float64).values, 1.0
    )
    cutoff = float(np.quantile(density, MAX_COUNTS_QUANTILE))
    print(f"  density q{int(MAX_COUNTS_QUANTILE*100)} cutoff: {cutoff:.4f} counts/area")

    keep = density <= cutoff
    n_drop_density = int((~keep).sum())
    adata = adata[keep].copy()
    n3 = adata.n_obs
    print(f"  after density filter: {n3:,}  (dropped {n_drop_density:,})")

    sc.pp.filter_genes(adata, min_cells=MIN_CELLS)
    n4 = adata.n_obs
    n_genes_kept = adata.n_vars
    print(f"  after filter_genes(min_cells={MIN_CELLS}): {n4:,} cells x {n_genes_kept:,} genes")

    # Step 3: compare to saved _preprocessed.h5ad
    pre_path = ROOT / f"data/processed/{sample}/{sample}_preprocessed.h5ad"
    with h5py.File(pre_path, "r") as f:
        n_disk = f["obs/total_counts"].shape[0]
        if "_index" in f["var"]:
            disk_var_names = f["var/_index"].asstr()[:]
        else:
            disk_var_names = f["var/" + f["var"].attrs.get("_index", "_index")].asstr()[:]

    print(f"\nStep 3 — disk verification:")
    print(f"  _preprocessed.h5ad cell count: {n_disk:,}")
    print(f"  recomputed cell count:         {n4:,}")
    if n_disk == n4:
        print(f"  EXACT MATCH ✓ filter logic is reproducible")
    else:
        delta = n_disk - n4
        print(f"  MISMATCH: diff={delta:,} (script may have a bug if delta != 0)")

    # Step 4: n_genes distribution before vs after
    n_genes_post = adata.obs["n_genes_by_counts"].values
    print(f"\nStep 4 — n_genes_by_counts distribution:")
    print(f"  panel size (post min_cells): {n_genes_kept:,} genes")
    print(f"  PRE filter (n={len(n_genes_pre):,}):")
    print(f"    mean={np.mean(n_genes_pre):.0f}  med={np.median(n_genes_pre):.0f}  "
          f"p50={np.quantile(n_genes_pre,0.50):.0f}  "
          f"p90={np.quantile(n_genes_pre,0.90):.0f}  "
          f"p98={np.quantile(n_genes_pre,0.98):.0f}  "
          f"p999={np.quantile(n_genes_pre,0.999):.0f}  "
          f"max={int(n_genes_pre.max())}")
    print(f"  POST filter (n={len(n_genes_post):,}):")
    print(f"    mean={np.mean(n_genes_post):.0f}  med={np.median(n_genes_post):.0f}  "
          f"p50={np.quantile(n_genes_post,0.50):.0f}  "
          f"p90={np.quantile(n_genes_post,0.90):.0f}  "
          f"p98={np.quantile(n_genes_post,0.98):.0f}  "
          f"p999={np.quantile(n_genes_post,0.999):.0f}  "
          f"max={int(n_genes_post.max())}")
    n_above_1200 = int((n_genes_post > 1200).sum())
    print(f"  cells with n_genes > 1200 in POST: {n_above_1200:,} "
          f"({100.0*n_above_1200/len(n_genes_post):.3f}%)")

    # Step 5: cross-ref with celltype from phenotyped h5ad
    pheno_path = ROOT / f"data/processed/{sample}/{sample}_phenotyped.h5ad"
    with h5py.File(pheno_path, "r") as f:
        pheno_n_genes = f["obs/n_genes_by_counts"][:]
        # celltype categorical
        ct_codes = f["obs/celltype/codes"][:]
        ct_cats = f["obs/celltype/categories"].asstr()[:]
        celltype = np.array(ct_cats)[ct_codes]

    print(f"\nStep 5 — high-n_genes cells by celltype (POST-pipeline):")
    print(f"  phenotyped n_obs: {len(pheno_n_genes):,}  (matches recomputed: "
          f"{'yes' if len(pheno_n_genes) == n4 else 'NO — pipeline n_obs drift'})")
    high_mask = pheno_n_genes > 1200
    print(f"  cells with n_genes > 1200 in phenotyped: {int(high_mask.sum()):,} "
          f"({100*high_mask.sum()/len(pheno_n_genes):.3f}%)")

    # Distribution of celltype within high-n_genes cells
    if high_mask.sum() > 0:
        df = pd.DataFrame({"celltype": celltype, "high": high_mask})
        ct_counts = df["celltype"].value_counts()
        ct_high_counts = df[df["high"]]["celltype"].value_counts()
        print(f"\n  celltype | n_high (>1200) | n_total | %_high_among_type | %_of_high_total")
        print(f"  {'-'*72}")
        rows = []
        for ct in ct_counts.index:
            nh = int(ct_high_counts.get(ct, 0))
            nt = int(ct_counts[ct])
            pct_within = 100.0 * nh / max(nt, 1)
            pct_of_high = 100.0 * nh / max(int(high_mask.sum()), 1)
            rows.append((ct, nh, nt, pct_within, pct_of_high))
        # sort by n_high desc
        rows.sort(key=lambda r: -r[1])
        for ct, nh, nt, pw, pa in rows:
            print(f"  {ct:<14} | {nh:>14,} | {nt:>7,} | {pw:>17.2f}% | {pa:>15.2f}%")

    # Step 6: sanity-check density filter
    print(f"\nStep 6 — density filter sanity check:")
    print(f"  q{int(MAX_COUNTS_QUANTILE*100)} of density = {cutoff:.3f}")
    print(f"  density of dropped cells:")
    dropped_density = density[~keep]
    print(f"    n={len(dropped_density):,}  min={dropped_density.min():.3f}  "
          f"med={np.median(dropped_density):.3f}  max={dropped_density.max():.3f}")
    print(f"  total_counts of dropped cells:")
    dropped_tc = total_counts_pre[adata.obs.index.size:][:0]  # placeholder
    # Re-derive: dropped cells are pre-filter cells with density > cutoff
    # But we already mutated adata, so use the pre-filter arrays
    pre_density = total_counts_pre / np.maximum(cell_area_pre, 1.0)
    pre_drop_mask = pre_density > cutoff
    print(f"    n={int(pre_drop_mask.sum()):,}  "
          f"mean total_counts={np.mean(total_counts_pre[pre_drop_mask]):.0f}  "
          f"max total_counts={int(total_counts_pre[pre_drop_mask].max())}")
    print(f"    mean cell_area={np.mean(cell_area_pre[pre_drop_mask]):.0f}  "
          f"mean n_genes={np.mean(n_genes_pre[pre_drop_mask]):.0f}")
    print(f"  KEPT cells stats:")
    pre_keep_mask = ~pre_drop_mask
    print(f"    mean total_counts={np.mean(total_counts_pre[pre_keep_mask]):.0f}  "
          f"mean cell_area={np.mean(cell_area_pre[pre_keep_mask]):.0f}  "
          f"mean n_genes={np.mean(n_genes_pre[pre_keep_mask]):.0f}")


for sample in SAMPLES:
    main(sample)

print("\n" + "="*70)
print("VERIFICATION COMPLETE")
print("="*70)
