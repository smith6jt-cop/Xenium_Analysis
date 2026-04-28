"""Re-render 0041323 co_occurrence_and_spatial figure with both panels."""
import time
import sys
from pathlib import Path

import numpy as np
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt

SAMPLE = sys.argv[1] if len(sys.argv) > 1 else "0041323"
DATA = Path(f"data/processed/{SAMPLE}/{SAMPLE}_preprocessed.h5ad")
FIG_DIR = Path(f"figures/{SAMPLE}/01_preprocessing")

print(f"Loading {DATA} ...")
t = time.time()
adata = sc.read_h5ad(DATA)
print(f"  loaded {adata.shape} in {time.time()-t:.1f}s")
print(f"  leiden_1.5 in obs: {'leiden_1.5' in adata.obs.columns}")

N_SUB = 200_000
if adata.n_obs > N_SUB:
    rng = np.random.default_rng(0)
    idx = rng.choice(adata.n_obs, size=N_SUB, replace=False)
    adata_sub = adata[idx].copy()
else:
    adata_sub = adata
print(f"  co-occurrence subsample: {adata_sub.shape}")

t = time.time()
sq.gr.spatial_neighbors(adata_sub, coord_type="generic", n_neighs=10, delaunay=False)
print(f"  spatial_neighbors: {time.time()-t:.1f}s")

t = time.time()
sq.gr.co_occurrence(adata_sub, cluster_key="leiden_1.5", n_jobs=4, show_progress_bar=False)
print(f"  co_occurrence: {time.time()-t:.1f}s")

sizes = adata_sub.obs["leiden_1.5"].value_counts()
cluster_for_co = str(sizes.index[1]) if len(sizes) > 1 else str(sizes.index[0])
print(f"  highlighting cluster {cluster_for_co}")

# sq.pl.co_occurrence(..., ax=ax) is broken in squidpy 1.6.5 (duplicate `ax`
# kwarg from internal sns.lineplot). Pull the array directly from .uns and
# plot manually instead.
co_uns = adata_sub.uns["leiden_1.5_co_occurrence"]
occ = co_uns["occ"]            # shape: (n_clusters, n_clusters, n_intervals)
interval = co_uns["interval"]  # 1D distance bins (length n_intervals + 1)
cluster_categories = list(adata_sub.obs["leiden_1.5"].cat.categories)
cluster_idx = cluster_categories.index(cluster_for_co)

# Use bin midpoints as x-axis
mids = 0.5 * (interval[:-1] + interval[1:])
# occ shape is (n_clusters_a, n_clusters_b, n_intervals); slice for the focal cluster
focal = occ[cluster_idx, :, :]  # (n_clusters, n_intervals)

fig, axes = plt.subplots(1, 2, figsize=(20, 10))
import matplotlib as mpl
cmap = mpl.colormaps.get_cmap("tab20").resampled(max(20, len(cluster_categories)))
for j, other in enumerate(cluster_categories):
    axes[0].plot(mids, focal[j], label=other, color=cmap(j), linewidth=1.2, alpha=0.85)
axes[0].set_xlabel("distance"); axes[0].set_ylabel("p(other | focal cluster)")
axes[0].set_title(f"Co-occurrence: focal cluster {cluster_for_co} (subsample n={adata_sub.n_obs:,})")
axes[0].grid(True, alpha=0.3)
axes[0].legend(title="other cluster", ncol=2, fontsize=7, loc="best")

sq.pl.spatial_scatter(adata, color="leiden_1.5", shape=None, size=2, ax=axes[1])
axes[1].set_title("Spatial layout (full sample, leiden_1.5)")
plt.tight_layout()
plt.savefig(FIG_DIR / "co_occurrence_and_spatial.png", dpi=200, bbox_inches="tight")
print(f"  wrote {FIG_DIR}/co_occurrence_and_spatial.png")
