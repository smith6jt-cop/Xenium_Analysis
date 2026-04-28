"""End-to-end verification of all 6 pipeline outputs."""
import h5py
import numpy as np
from scipy.stats import spearmanr

print(f"{'sample':<10}{'stage':<10}{'n_obs':>10}{'key':>22}{'value':>12}")
print("-" * 70)

for s in ("0041323", "0041326"):
    # Preprocessed
    with h5py.File(f"data/processed/{s}/{s}_preprocessed.h5ad", "r") as f:
        n = f["obs/total_counts"].shape[0]
        cats = (
            f["obs/leiden_1.5/categories"].asstr()[:]
            if "leiden_1.5/categories" in f["obs"]
            else []
        )
        nc15 = len(set(cats))
        has_pca = "X_pca" in f["obsm"]
        has_umap = "X_umap" in f["obsm"]
        has_counts = "counts" in f["layers"]
        has_lognorm = "lognorm" in f["layers"]
    print(f"{s:<10}{'prep':<10}{n:>10,}{'leiden_1.5 nclu':>22}{nc15:>12}")
    print(
        f"{'':<10}{'':<10}{'':>10}{'layers OK':>22}"
        f"{str(has_counts and has_lognorm):>12}"
    )
    print(
        f"{'':<10}{'':<10}{'':>10}{'obsm OK':>22}"
        f"{str(has_pca and has_umap):>12}"
    )

    # Phenotyped
    with h5py.File(f"data/processed/{s}/{s}_phenotyped.h5ad", "r") as f:
        n = f["obs/total_counts"].shape[0]
        cats_scvi = (
            f["obs/leiden_scvi/categories"].asstr()[:]
            if "leiden_scvi/categories" in f["obs"]
            else []
        )
        nc_scvi = len(set(cats_scvi))
        um = (
            f["obsm/X_umap_scvi"][:]
            if "X_umap_scvi" in f["obsm"]
            else f["obsm/X_umap"][:]
        )
        tc = f["obs/total_counts"][:]
        rx, _ = spearmanr(um[:, 0], tc)
        ry, _ = spearmanr(um[:, 1], tc)
        ct_cats = (
            f["obs/celltype/categories"].asstr()[:]
            if "celltype/categories" in f["obs"]
            else []
        )
    print(f"{s:<10}{'pheno':<10}{n:>10,}{'leiden_scvi nclu':>22}{nc_scvi:>12}")
    print(f"{'':<10}{'':<10}{'':>10}{'Spearman x':>22}{rx:>12.3f}")
    print(f"{'':<10}{'':<10}{'':>10}{'Spearman y':>22}{ry:>12.3f}")
    print(f"{'':<10}{'':<10}{'':>10}{'celltype types':>22}{len(ct_cats):>12}")

    # Spatial
    with h5py.File(f"data/processed/{s}/{s}_spatial_analysis.h5ad", "r") as f:
        n = f["obs/total_counts"].shape[0]
        niche_cats = (
            f["obs/spatial_niche/categories"].asstr()[:]
            if "spatial_niche/categories" in f["obs"]
            else []
        )
        n_niches = len(set(niche_cats))
        has_moran = "moranI" in f["uns"]
    print(f"{s:<10}{'spatial':<10}{n:>10,}{'spatial_niche n':>22}{n_niches:>12}")
    print(f"{'':<10}{'':<10}{'':>10}{'moranI uns':>22}{str(has_moran):>12}")
    print()

print("ALL CHECKS COMPLETE")
