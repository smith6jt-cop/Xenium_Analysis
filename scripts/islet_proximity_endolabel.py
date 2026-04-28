"""Alternative islet identification using celltype_lineage label as seed.

The strict hormone-marker seed (lognorm >= log1p(5) of CHGA/INSM1/etc.)
recovers 110 islets in 0041323 but only 2 in 0041326, because 0041326
has genuinely lower hormone marker expression (median total_counts 309
vs 288, with 298 panel genes dropped by min_cells=100 vs 14 in 0041323).

This script seeds DBSCAN with cells already assigned an endocrine label
by score_genes argmax (celltype_lineage in {Endocrine, Beta, Alpha,
Delta, Endocrine_pan}). Epsilon is excluded — its 150k count in 0041326
is artefactual (likely a panel quirk; ghrelin cells are <1% of pancreas).

Outputs (parallel to the strict-seed CSVs, with _endolabel suffix):
  data/processed/islets_dbscan_endolabel.csv
  data/processed/immune_proximity_summary_endolabel.csv
  data/processed/islet_infiltration_per100endo_endolabel.csv
"""
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path("/blue/maigan/smith6jt/Xenium_Analysis")
SAMPLES = ("0041323", "0041326")

ENDO_LABELS = {"Endocrine", "Beta", "Alpha", "Delta", "Endocrine_pan"}
EPS_UM = 50.0
MIN_SAMPLES = 10
INTRA_UM = 25.0
PERI_UM = 50.0
PROXIMAL_UM = 200.0


def dbscan_islets(coords):
    try:
        from cuml.cluster import DBSCAN as cuDBSCAN
        import cupy as cp
        t = time.time()
        db = cuDBSCAN(eps=EPS_UM, min_samples=MIN_SAMPLES)
        labels = db.fit_predict(cp.asarray(coords, dtype=cp.float32))
        labels = labels.get() if hasattr(labels, "get") else np.asarray(labels)
        print(f"    cuML DBSCAN: {time.time()-t:.1f}s")
    except Exception as e:
        print(f"    cuML DBSCAN unavailable ({e}); falling back to sklearn")
        from sklearn.cluster import DBSCAN
        t = time.time()
        labels = DBSCAN(eps=EPS_UM, min_samples=MIN_SAMPLES, n_jobs=-1).fit_predict(coords)
        print(f"    sklearn DBSCAN: {time.time()-t:.1f}s")
    return np.asarray(labels)


def nn_distance(query_xy, seed_xy):
    try:
        from cuml.neighbors import NearestNeighbors as cuNN
        import cupy as cp
        t = time.time()
        nn = cuNN(n_neighbors=1, metric="euclidean")
        nn.fit(cp.asarray(seed_xy, dtype=cp.float32))
        d, idx = nn.kneighbors(cp.asarray(query_xy, dtype=cp.float32))
        d = d.get().ravel() if hasattr(d, "get") else np.asarray(d).ravel()
        idx = idx.get().ravel() if hasattr(idx, "get") else np.asarray(idx).ravel()
        print(f"    cuML NN: {time.time()-t:.1f}s")
    except Exception as e:
        print(f"    cuML NN unavailable ({e}); falling back to sklearn")
        from sklearn.neighbors import NearestNeighbors
        t = time.time()
        nn = NearestNeighbors(n_neighbors=1, metric="euclidean", n_jobs=-1)
        nn.fit(seed_xy)
        d, idx = nn.kneighbors(query_xy)
        d = d.ravel(); idx = idx.ravel()
        print(f"    sklearn NN: {time.time()-t:.1f}s")
    return d, idx


def process_sample(sample):
    print(f"\n{'='*72}\n=== {sample}\n{'='*72}")
    p = ROOT / f"data/processed/{sample}/{sample}_phenotyped.h5ad"
    print(f"  loading {p.name} (backed=r)...")
    t = time.time()
    a = sc.read_h5ad(p, backed="r")
    print(f"    loaded shape={a.shape} in {time.time()-t:.1f}s")
    coords_all = np.asarray(a.obsm["spatial"])
    ct_lineage = a.obs["celltype_lineage"].astype(str).values
    is_endo = np.isin(ct_lineage, list(ENDO_LABELS))
    n_endo = int(is_endo.sum())
    print(f"    endocrine-labeled cells: {n_endo:,}  "
          f"(excludes Epsilon)")
    print(f"    breakdown: ", {l: int((ct_lineage == l).sum())
                                  for l in sorted(ENDO_LABELS)
                                  if (ct_lineage == l).sum() > 0})
    seed_xy = coords_all[is_endo]
    seed_idx_in_full = np.where(is_endo)[0]

    print(f"  DBSCAN(eps={EPS_UM}, min_samples={MIN_SAMPLES}) on {n_endo:,} seeds...")
    labels = dbscan_islets(seed_xy)
    n_islets = int(labels[labels >= 0].max() + 1) if (labels >= 0).any() else 0
    n_clustered = int((labels >= 0).sum())
    print(f"    n_islets={n_islets}  clustered_seeds={n_clustered:,}  "
          f"noise={int((labels < 0).sum()):,}")

    if n_islets == 0:
        print(f"  WARNING: no islets found for {sample}; skipping")
        return None, None, None

    # Build per-islet records
    islet_rows = []
    for k in range(n_islets):
        mk = labels == k
        xy = seed_xy[mk]
        islet_id = f"{sample}_islet_{k:04d}"
        islet_rows.append({
            "sample": sample,
            "islet_id": islet_id,
            "centroid_x": float(xy[:, 0].mean()),
            "centroid_y": float(xy[:, 1].mean()),
            "n_seed_cells": int(mk.sum()),
            "x_min": float(xy[:, 0].min()),
            "x_max": float(xy[:, 0].max()),
            "y_min": float(xy[:, 1].min()),
            "y_max": float(xy[:, 1].max()),
            "approx_radius_um": float(np.sqrt(((xy - xy.mean(0))**2).sum(1)).max()),
        })
    islets_df = pd.DataFrame(islet_rows)
    print(f"    median seeds/islet: {islets_df['n_seed_cells'].median():.0f}  "
          f"max: {islets_df['n_seed_cells'].max()}")

    # Load immune labels for this sample (immune_phenotyped.h5ad has both
    # samples concatenated; subset to this sample).
    p_imm = ROOT / f"data/processed/{sample}/{sample}_immune_phenotyped.h5ad"
    print(f"  loading {p_imm.name}...")
    a.file.close()
    t = time.time()
    a_imm = sc.read_h5ad(p_imm)
    print(f"    loaded shape={a_imm.shape} in {time.time()-t:.1f}s")
    if "sample" in a_imm.obs.columns:
        a_imm = a_imm[a_imm.obs["sample"].astype(str) == sample].copy()
        print(f"    after sample filter: {a_imm.shape}")
    immune_xy = np.asarray(a_imm.obsm["spatial"])
    immune_subtype = a_imm.obs["immune_subtype"].astype(str).values

    # Restrict seeds to clustered (islet-member) seeds only
    in_islet_mask = labels >= 0
    seed_xy_clust = seed_xy[in_islet_mask]
    seed_label_clust = labels[in_islet_mask]

    print(f"  computing distance to nearest clustered seed for "
          f"{a_imm.n_obs:,} immune cells...")
    d, nn_idx = nn_distance(immune_xy, seed_xy_clust)
    nearest_islet_id = seed_label_clust[nn_idx]

    # Distance bins
    bins = np.full(len(d), "distal", dtype=object)
    bins[d < PROXIMAL_UM] = "proximal"
    bins[d < PERI_UM] = "intra_or_peri_islet"

    # Contingency: subtype × bin
    cont = pd.crosstab(immune_subtype, bins).reindex(columns=[
        "intra_or_peri_islet", "proximal", "distal"], fill_value=0)
    cont = cont.reset_index().rename(columns={"index": "immune_subtype"})
    cont.insert(0, "sample", sample)
    cont.columns = ["sample", "immune_subtype", "intra_or_peri_islet",
                    "proximal", "distal"]

    # Per-islet immune counts (cells assigned to islet AND within PERI_UM)
    in_peri = d < PERI_UM
    df_imm = pd.DataFrame({
        "subtype": immune_subtype[in_peri],
        "islet_id": [f"{sample}_islet_{k:04d}" for k in nearest_islet_id[in_peri]],
    })
    per_islet = (df_imm.groupby(["islet_id", "subtype"], observed=True)
                  .size().unstack("subtype", fill_value=0))
    # Ensure all subtypes present (use union from the contingency)
    all_subtypes = sorted(set(immune_subtype))
    for st in all_subtypes:
        if st not in per_islet.columns:
            per_islet[st] = 0
    per_islet = per_islet[all_subtypes].reset_index()

    # Merge n_endocrine from islets_df
    seeds_df = islets_df[["islet_id", "n_seed_cells"]].rename(
        columns={"n_seed_cells": "n_endocrine"})
    per_islet = per_islet.merge(seeds_df, on="islet_id", how="right").fillna(0)
    per_islet["sample"] = sample

    # per100endo
    for st in all_subtypes:
        per_islet[f"{st}_per100endo"] = (
            100.0 * per_islet[st] / per_islet["n_endocrine"].clip(lower=1))

    return islets_df, cont, per_islet


all_islets, all_cont, all_perislet = [], [], []
for s in SAMPLES:
    isl, cont, per = process_sample(s)
    if isl is not None:
        all_islets.append(isl)
        all_cont.append(cont)
        all_perislet.append(per)

islets_df = pd.concat(all_islets, ignore_index=True)
cont_df = pd.concat(all_cont, ignore_index=True)
perislet_df = pd.concat(all_perislet, ignore_index=True)

out_islets = ROOT / "data/processed/islets_dbscan_endolabel.csv"
out_cont = ROOT / "data/processed/immune_proximity_summary_endolabel.csv"
out_per = ROOT / "data/processed/islet_infiltration_per100endo_endolabel.csv"

islets_df.to_csv(out_islets, index=False)
cont_df.to_csv(out_cont, index=False)
perislet_df.to_csv(out_per, index=False)

print(f"\n{'='*72}\n=== SUMMARY\n{'='*72}")
print(f"  total islets:           {len(islets_df):,}")
for s in SAMPLES:
    n = (islets_df["sample"] == s).sum()
    print(f"    {s}: {n:,} islets")
print(f"  wrote: {out_islets}")
print(f"  wrote: {out_cont}")
print(f"  wrote: {out_per}")
print(f"\nContingency (subtype × bin × sample):")
print(cont_df.to_string(index=False))
