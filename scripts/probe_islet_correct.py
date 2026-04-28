"""Probe DBSCAN with the CORRECT n_islets calculation."""
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.cluster import DBSCAN

SPECIFIC_HORMONE = ["IAPP", "MAFA", "PDX1", "NKX6-1", "SLC30A8", "ARX", "MAFB",
                    "GHRL", "HHEX"]

for sample in ("0041323", "0041326"):
    print(f"\n=== {sample} ===")
    a = sc.read_h5ad(f"data/processed/{sample}/{sample}_phenotyped.h5ad")
    present = [g for g in SPECIFIC_HORMONE if g in a.var_names]
    idxs = [a.var_names.get_loc(g) for g in present]
    sub = a.layers["counts"][:, idxs]
    if hasattr(sub, "toarray"):
        sub = sub.toarray()
    sub = np.asarray(sub)

    # Seed: any specific hormone marker at raw >= 5
    seed = (sub >= 5).any(axis=1)
    coords = np.asarray(a.obsm["spatial"], dtype=np.float32)[seed]
    print(f"  seed cells (any hormone marker raw>=5): {coords.shape[0]:,}")

    print(f"  DBSCAN parameter sweep:")
    print(f"    {'eps':>5} {'min':>5} {'islets':>8} {'med':>6} {'p90':>6} {'max':>6} {'noise':>8}")
    for eps in (20, 30, 50, 80):
        for ms in (5, 10, 20):
            db = DBSCAN(eps=float(eps), min_samples=int(ms), n_jobs=-1)
            labels = db.fit_predict(coords)
            cluster_labels = labels[labels >= 0]
            n_islets = int(cluster_labels.max() + 1) if cluster_labels.size else 0
            n_noise = int((labels == -1).sum())
            sizes = pd.Series(cluster_labels).value_counts() if cluster_labels.size else pd.Series([], dtype=int)
            med = int(sizes.median()) if len(sizes) else 0
            p90 = int(sizes.quantile(0.90)) if len(sizes) else 0
            mx = int(sizes.max()) if len(sizes) else 0
            print(f"    {eps:>5} {ms:>5} {n_islets:>8,} {med:>6} {p90:>6} {mx:>6} {n_noise:>8,}")

    # Also try the broader endocrine label for reference
    if "celltype_lineage" in a.obs.columns:
        broad_mask = a.obs["celltype_lineage"].astype(str).isin(
            {"Beta", "Alpha", "Delta", "Epsilon", "Endocrine", "Endocrine_pan"}
        ).values
        coords_broad = np.asarray(a.obsm["spatial"], dtype=np.float32)[broad_mask]
        print(f"\n  reference: broad endocrine label ({coords_broad.shape[0]:,} cells)")
        for eps, ms in [(30, 20), (50, 20)]:
            db = DBSCAN(eps=float(eps), min_samples=int(ms), n_jobs=-1)
            labels = db.fit_predict(coords_broad)
            cluster_labels = labels[labels >= 0]
            n_islets = int(cluster_labels.max() + 1) if cluster_labels.size else 0
            n_noise = int((labels == -1).sum())
            print(f"    eps={eps} min={ms}: {n_islets:,} islets, {n_noise:,} noise")
