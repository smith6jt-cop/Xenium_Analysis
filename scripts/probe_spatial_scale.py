"""Probe spatial coordinate scale to understand DBSCAN behavior."""
import numpy as np
import scanpy as sc
from scipy.spatial import cKDTree


for sample in ("0041323", "0041326"):
    print(f"\n=== {sample} ===")
    a = sc.read_h5ad(f"data/processed/{sample}/{sample}_phenotyped.h5ad")
    xy = np.asarray(a.obsm["spatial"], dtype=np.float64)
    print(f"  obsm['spatial'] dtype={a.obsm['spatial'].dtype}  "
          f"shape={a.obsm['spatial'].shape}")
    print(f"  range: x [{xy[:,0].min():.1f}, {xy[:,0].max():.1f}]  "
          f"y [{xy[:,1].min():.1f}, {xy[:,1].max():.1f}]")
    print(f"  span : x={xy[:,0].max()-xy[:,0].min():.0f}  "
          f"y={xy[:,1].max()-xy[:,1].min():.0f}")
    print(f"  mean cell_area: {a.obs['cell_area'].mean():.1f}  "
          f"sqrt={np.sqrt(a.obs['cell_area'].mean()):.1f}")

    rng = np.random.default_rng(0)
    idx = rng.choice(xy.shape[0], size=min(50_000, xy.shape[0]), replace=False)
    tree = cKDTree(xy[idx])
    d, _ = tree.query(xy[idx], k=2)
    nn_d = d[:, 1]
    print(f"  all-cell NN dist (50k subsample):")
    print(f"    p1={np.quantile(nn_d, 0.01):.2f}  p10={np.quantile(nn_d, 0.10):.2f}  "
          f"med={np.median(nn_d):.2f}  p90={np.quantile(nn_d, 0.90):.2f}")

    if "celltype_lineage" in a.obs.columns:
        endo_mask = a.obs["celltype_lineage"].astype(str).isin(
            {"Beta", "Alpha", "Delta", "Epsilon", "Endocrine", "Endocrine_pan"}
        ).values
        endo_xy = xy[endo_mask]
        if len(endo_xy) >= 100:
            tree_e = cKDTree(endo_xy)
            de, _ = tree_e.query(endo_xy, k=2)
            nn_e = de[:, 1]
            print(f"  endocrine-only NN dist (n={len(endo_xy):,}):")
            print(f"    p1={np.quantile(nn_e, 0.01):.2f}  p10={np.quantile(nn_e, 0.10):.2f}  "
                  f"med={np.median(nn_e):.2f}  p90={np.quantile(nn_e, 0.90):.2f}  "
                  f"max={nn_e.max():.1f}")
