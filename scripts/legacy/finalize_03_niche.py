"""Recompute spatial niches on the phenotyped h5ad with spatial-smoothed clustering
and save `<sample>_spatial_analysis.h5ad`.

Method (Banksy-style, simplified):
  1. For each cell, compute composition vector = fraction of each celltype among
     its K=50 nearest spatial neighbors.
  2. MiniBatchKMeans(n_clusters=12, n_init=20) on composition vectors → raw niches.
  3. Spatial smoothing: for each cell, replace its niche label with the mode of the
     niche labels among its 30 nearest spatial neighbors. Repeat 2 passes. This
     removes salt-and-pepper and produces contiguous tissue domains.

Usage:  XENIUM_SAMPLE_ID=0041323 python scripts/finalize_03_niche.py
"""
import os
import gc
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import MiniBatchKMeans

SAMPLE_ID = os.environ.get("XENIUM_SAMPLE_ID", "0041323")
PROCESSED = Path(f"/blue/maigan/smith6jt/Xenium_Analysis/data/processed/{SAMPLE_ID}")
FIGDIR = Path(f"/blue/maigan/smith6jt/Xenium_Analysis/figures/{SAMPLE_ID}/03_spatial_analysis")
FIGDIR.mkdir(parents=True, exist_ok=True)

phen = PROCESSED / f"{SAMPLE_ID}_phenotyped.h5ad"
if not phen.exists():
    raise FileNotFoundError(phen)

print(f"Loading {phen} ...")
adata = sc.read_h5ad(phen)
print(f"  shape={adata.shape}")

# Slim to reduce memory
for k in list(adata.uns.keys()):
    if k in {"_scvi_manager_uuid", "_scvi_uuid", "scvi",
             "leiden_1.5_nhood_enrichment", "leiden_1.5_centrality_scores",
             "dendrogram_leiden_1.5", "rank_genes_leiden",
             "hvg", "log1p", "neighbors", "pca", "umap", "spatial_neighbors"}:
        del adata.uns[k]
for k in list(adata.obsp.keys()):
    if k in {"connectivities", "distances", "scvi_connectivities", "scvi_distances"}:
        del adata.obsp[k]
for k in list(adata.layers.keys()):
    if k in {"log_normalized", "scaled"}:
        del adata.layers[k]
gc.collect()

coords = np.asarray(adata.obsm["spatial"], dtype=np.float64)
celltypes = adata.obs["celltype"].astype("category")
ct_cats = list(celltypes.cat.categories)
ct_codes = celltypes.cat.codes.to_numpy()
n_cells = coords.shape[0]

# --- 1. spatial k-NN for composition (K=50 for more context) ---
K_COMP = 50
print(f"Building spatial k-NN (k={K_COMP}) over {n_cells:,} cells ...")
nn = NearestNeighbors(n_neighbors=K_COMP, algorithm="kd_tree", n_jobs=-1)
nn.fit(coords)
nbrs = nn.kneighbors(coords, return_distance=False)
neighbor_codes = ct_codes[nbrs]

comp = np.zeros((n_cells, len(ct_cats)), dtype=np.float32)
for code in range(len(ct_cats)):
    comp[:, code] = (neighbor_codes == code).sum(axis=1)
comp /= K_COMP

# --- 2. MiniBatchKMeans with strong init ---
n_niches = 12
print(f"MiniBatchKMeans clustering into {n_niches} niches (n_init=20) ...")
km = MiniBatchKMeans(
    n_clusters=n_niches,
    random_state=42,
    batch_size=16384,
    n_init=20,
    max_iter=300,
    init='k-means++',
)
raw = km.fit_predict(comp).astype(np.int32)
print(f"Raw niche sizes: {np.bincount(raw).tolist()}")

# --- 3. Spatial smoothing: majority vote over K_SMOOTH nearest spatial neighbors ---
K_SMOOTH = 30
if K_SMOOTH > K_COMP:
    nn_smooth = NearestNeighbors(n_neighbors=K_SMOOTH, algorithm="kd_tree", n_jobs=-1)
    nn_smooth.fit(coords)
    sm_nbrs = nn_smooth.kneighbors(coords, return_distance=False)
else:
    sm_nbrs = nbrs[:, :K_SMOOTH]

def majority_vote(labels, nb_idx):
    """For each row i, return the most common label among nb_idx[i] labels."""
    out = labels.copy()
    # Vectorized majority via bincount per row would be expensive; do in chunks.
    n_labels = int(labels.max()) + 1
    chunk = 50_000
    for i in range(0, len(labels), chunk):
        idx = slice(i, min(i + chunk, len(labels)))
        nbr_labels = labels[nb_idx[idx]]  # shape (chunk, K_SMOOTH)
        # One-hot count per row, argmax
        oh = np.eye(n_labels, dtype=np.int32)[nbr_labels]  # (chunk, K_SMOOTH, n_labels)
        counts = oh.sum(axis=1)  # (chunk, n_labels)
        out[idx] = counts.argmax(axis=1)
    return out

smoothed = raw
for pass_i in range(2):
    new = majority_vote(smoothed, sm_nbrs)
    changed = (new != smoothed).sum()
    smoothed = new
    print(f"  smoothing pass {pass_i + 1}: {changed:,} cells changed label")

adata.obs["spatial_niche"] = pd.Categorical(smoothed.astype(str))
print(f"Final niche sizes:\n{adata.obs['spatial_niche'].value_counts().sort_index()}")

# --- 4. Plots ---
import matplotlib as mpl
import matplotlib.pyplot as plt

rng = np.random.default_rng(0)
idx = rng.choice(n_cells, size=min(300_000, n_cells), replace=False)

fig, axes = plt.subplots(1, 2, figsize=(18, 8))
axes[0].scatter(coords[idx, 0], coords[idx, 1],
                c=smoothed[idx], cmap="tab20", s=0.5, linewidths=0)
axes[0].set_title(f"Spatial niches (K_COMP={K_COMP}, smoothed, n_niches={n_niches})")
axes[0].set_xlim(coords[:, 0].min(), coords[:, 0].max())
axes[0].set_ylim(coords[:, 1].min(), coords[:, 1].max())
axes[0].set_aspect("equal", adjustable="datalim")

sub_codes = celltypes.iloc[idx].cat.codes.to_numpy()
cmap = mpl.colormaps.get_cmap("tab20").resampled(max(20, len(ct_cats)))
axes[1].scatter(coords[idx, 0], coords[idx, 1],
                c=sub_codes, cmap=cmap, s=0.5, linewidths=0)
axes[1].set_title("Cell types")
axes[1].set_xlim(coords[:, 0].min(), coords[:, 0].max())
axes[1].set_ylim(coords[:, 1].min(), coords[:, 1].max())
axes[1].set_aspect("equal", adjustable="datalim")

fig.savefig(FIGDIR / "spatial_domains.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print("spatial_domains.png saved")

# Composition heatmap
import seaborn as sns
niche_ct_matrix = pd.crosstab(
    pd.Series(smoothed, name="niche"),
    celltypes,
    normalize="index",
).fillna(0.0)
fig2, ax2 = plt.subplots(1, 1, figsize=(max(8, 0.8 * len(ct_cats)), 8))
sns.heatmap(niche_ct_matrix, cmap="viridis", ax=ax2, annot=True, fmt=".2f",
            cbar_kws={"label": "fraction"})
ax2.set_title("Niche composition (rows: niches, cols: celltype fractions)")
fig2.savefig(FIGDIR / "spatial_niche_composition.png", dpi=200, bbox_inches="tight")
plt.close(fig2)
print("spatial_niche_composition.png saved")

# Niche CSV
niche_csv = PROCESSED / f"{SAMPLE_ID}_spatial_niches.csv"
adata.obs[["celltype", "spatial_niche"]].to_csv(niche_csv)
print(f"niche CSV -> {niche_csv}")

# Sanitize uns (ligrec MultiIndex from earlier runs)
for k in list(adata.uns.keys()):
    v = adata.uns[k]
    if isinstance(v, dict) and any(not isinstance(ki, str) for ki in v.keys()):
        print(f"Dropping uns['{k}'] (non-str dict keys)")
        del adata.uns[k]

out = PROCESSED / f"{SAMPLE_ID}_spatial_analysis.h5ad"
print(f"Writing {out} ...")
adata.write_h5ad(out)
print("Done.")
