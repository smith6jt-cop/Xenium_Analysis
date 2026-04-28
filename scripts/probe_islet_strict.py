"""Test stricter islet seed using specific Beta/Alpha hormone markers."""
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.cluster import DBSCAN

# Specific markers that are NOT pan-endocrine (avoid CHGA/INSM1 overspill)
SPECIFIC_HORMONE = ["IAPP", "MAFA", "PDX1", "NKX6-1", "SLC30A8", "ARX", "MAFB",
                    "GHRL", "HHEX"]

for sample in ("0041323", "0041326"):
    print(f"\n=== {sample} ===")
    a = sc.read_h5ad(f"data/processed/{sample}/{sample}_phenotyped.h5ad")
    present = [g for g in SPECIFIC_HORMONE if g in a.var_names]
    print(f"  specific markers: {present}")
    idxs = [a.var_names.get_loc(g) for g in present]
    sub = a.layers["counts"][:, idxs]
    if hasattr(sub, "toarray"):
        sub = sub.toarray()
    sub = np.asarray(sub)

    print(f"  cell counts at various raw-count detection thresholds:")
    for raw_thresh in (3, 5, 10):
        for n_thresh in (1, 2):
            mask = (sub >= raw_thresh).sum(axis=1) >= n_thresh
            print(f"    raw>={raw_thresh}, >=in_{n_thresh}_markers: {int(mask.sum()):,}")

    # Best seed: any cell with >= 1 specific marker at raw >= 5 (high expression)
    best_seed = (sub >= 5).any(axis=1)
    print(f"\n  using seed: any specific marker raw>=5 → {int(best_seed.sum()):,} cells")
    coords = np.asarray(a.obsm["spatial"], dtype=np.float32)[best_seed]

    print(f"  DBSCAN parameter sweep on this seed:")
    print(f"    {'eps':>5} {'min':>5} {'islets':>8} {'med':>6} {'p90':>6} {'max':>8} {'noise':>8}")
    for eps in (20, 30, 50, 80, 120):
        for ms in (5, 10, 20):
            db = DBSCAN(eps=float(eps), min_samples=int(ms), n_jobs=-1)
            labels = db.fit_predict(coords)
            n_islets = int((labels >= 0).max() + 1) if (labels >= 0).any() else 0
            n_noise = int((labels == -1).sum())
            sizes = pd.Series(labels[labels >= 0]).value_counts() if n_islets else pd.Series([], dtype=int)
            med = int(sizes.median()) if len(sizes) else 0
            p90 = int(sizes.quantile(0.90)) if len(sizes) else 0
            mx = int(sizes.max()) if len(sizes) else 0
            print(f"    {eps:>5} {ms:>5} {n_islets:>8,} {med:>6} {p90:>6} {mx:>8} {n_noise:>8,}")
