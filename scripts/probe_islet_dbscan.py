"""Probe DBSCAN parameters + multi-marker islet seed for both samples."""
import time
import numpy as np
import pandas as pd
import scanpy as sc
from cuml.cluster import DBSCAN as cuDBSCAN

ENDOCRINE_MARKERS = [
    "IAPP", "MAFA", "NKX6-1", "PDX1", "SLC30A8",
    "ARX", "MAFB",
    "CHGA", "INSM1", "ISL1", "NEUROD1", "FEV",
    "HHEX", "GHRL",
]


def main(sample, threshold_markers=(1, 2, 3)):
    print(f"\n=== {sample} ===")
    a = sc.read_h5ad(f"data/processed/{sample}/{sample}_phenotyped.h5ad")
    present = [g for g in ENDOCRINE_MARKERS if g in a.var_names]
    print(f"  markers on panel: {present}")
    idxs = [a.var_names.get_loc(g) for g in present]
    sub = a.layers["counts"][:, idxs]
    if hasattr(sub, "toarray"):
        sub = sub.toarray()
    sub = np.asarray(sub)
    detected = (sub >= 3).sum(axis=1)

    print(f"  cells with >=N endocrine markers at raw>=3:")
    for n in (1, 2, 3, 4, 5):
        mask = detected >= n
        print(f"    >= {n}: {int(mask.sum()):,} cells")

    coords_full = np.asarray(a.obsm["spatial"], dtype=np.float32)

    # Sweep param grid for threshold=2 (the most useful default)
    for thresh in threshold_markers:
        seed_mask = detected >= thresh
        if seed_mask.sum() < 100:
            continue
        coords = coords_full[seed_mask]
        print(f"\n  --- seed: >= {thresh} markers, {coords.shape[0]:,} cells ---")
        header = f"{'eps':>5} {'min':>5} {'n_islets':>10} {'med_size':>10} {'p90_size':>10} {'noise':>10}"
        print(f"    {header}")
        for eps in (20, 30, 50, 80, 120):
            for min_samples in (5, 10, 20):
                db = cuDBSCAN(eps=float(eps), min_samples=int(min_samples),
                              output_type="numpy")
                labels = db.fit_predict(coords)
                n_islets = int((labels >= 0).max() + 1) if (labels >= 0).any() else 0
                n_noise = int((labels == -1).sum())
                if n_islets > 0:
                    counts_per = pd.Series(labels[labels >= 0]).value_counts()
                    med = int(counts_per.median())
                    p90 = int(counts_per.quantile(0.90))
                else:
                    med = p90 = 0
                print(
                    f"    {eps:>5} {min_samples:>5} {n_islets:>10,} {med:>10} "
                    f"{p90:>10} {n_noise:>10,}"
                )


for s in ("0041323", "0041326"):
    main(s)
