"""Compare sklearn DBSCAN vs cuML DBSCAN on the same input."""
import time
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.cluster import DBSCAN
from cuml.cluster import DBSCAN as cuDBSCAN

a = sc.read_h5ad("data/processed/0041323/0041323_phenotyped.h5ad")
endo_mask = a.obs["celltype_lineage"].astype(str).isin(
    {"Beta", "Alpha", "Delta", "Epsilon", "Endocrine", "Endocrine_pan"}
).values
coords = np.asarray(a.obsm["spatial"], dtype=np.float32)[endo_mask]
print(f"endocrine coords shape: {coords.shape}")

print("\nsklearn DBSCAN(eps=30, min_samples=20):")
t = time.time()
db = DBSCAN(eps=30.0, min_samples=20, n_jobs=-1)
labels_sk = db.fit_predict(coords)
print(f"  {time.time()-t:.1f}s")
n_islets_sk = int((labels_sk >= 0).max() + 1) if (labels_sk >= 0).any() else 0
n_noise_sk = int((labels_sk == -1).sum())
sk_sizes = pd.Series(labels_sk[labels_sk >= 0]).value_counts() if n_islets_sk else pd.Series()
print(f"  islets: {n_islets_sk:,}  noise: {n_noise_sk:,}")
if len(sk_sizes):
    print(f"  size dist: med={int(sk_sizes.median())}  p90={int(sk_sizes.quantile(0.9))}  max={sk_sizes.max()}")

print("\ncuML DBSCAN(eps=30, min_samples=20, output_type='numpy'):")
t = time.time()
db_cu = cuDBSCAN(eps=30.0, min_samples=20, output_type="numpy")
labels_cu = db_cu.fit_predict(coords)
print(f"  {time.time()-t:.1f}s")
n_islets_cu = int((labels_cu >= 0).max() + 1) if (labels_cu >= 0).any() else 0
n_noise_cu = int((labels_cu == -1).sum())
cu_sizes = pd.Series(labels_cu[labels_cu >= 0]).value_counts() if n_islets_cu else pd.Series()
print(f"  islets: {n_islets_cu:,}  noise: {n_noise_cu:,}")
if len(cu_sizes):
    print(f"  size dist: med={int(cu_sizes.median())}  p90={int(cu_sizes.quantile(0.9))}  max={cu_sizes.max()}")
