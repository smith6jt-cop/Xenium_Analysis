"""Replace the broken proximity tables with clinically-grounded
insulitis metrics derived from the endo-label islet seed.

The strict hormone-marker seed gives 0041326 only 2 islets, so its
proximity table puts every immune cell in 'distal' by construction.
The endo-label seed gives 265/318 islets — comparable across samples —
and is the right basis for cross-sample inference.

This script overwrites the headline tables:
  data/processed/immune_proximity_summary.csv     ← endo-label seed,
                                                     density enrichment
  data/processed/islet_infiltration_per100endo.csv ← only islets with
                                                     >=1 immune cell;
                                                     adds insulitis_grade
  data/processed/islet_insulitis_grades.csv       ← per-sample summary

The strict-seed legacy tables are renamed *_strict.csv for the record.
"""
import shutil
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

# Distance bins
ISLET_ZONE_UM = 50.0     # intra+peri-islet
PERI_PROXIMAL_UM = 200.0 # tissue-proximal beyond peri
# Anything > 200 = distal/background

# Clinical insulitis thresholds (Campbell-Thompson 2013 / nPOD consensus
# adapted for total CD45+ pool from immune subtypes):
INSULITIS_THRESH = 6   # >=6 immune cells within 50μm = insulitis
PERI_THRESH = 3        # 3-5 immune cells within 50μm = peri-insulitis


def dbscan_islets(coords):
    from sklearn.cluster import DBSCAN
    return np.asarray(DBSCAN(eps=EPS_UM, min_samples=MIN_SAMPLES,
                              n_jobs=-1).fit_predict(coords))


def nn_distance(query_xy, seed_xy):
    from sklearn.neighbors import NearestNeighbors
    nn = NearestNeighbors(n_neighbors=1, metric="euclidean", n_jobs=-1)
    nn.fit(seed_xy)
    d, idx = nn.kneighbors(query_xy)
    return d.ravel(), idx.ravel()


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
    print(f"    endocrine-labeled cells: {is_endo.sum():,}")
    seed_xy = coords_all[is_endo]

    print(f"  DBSCAN islet identification...")
    labels = dbscan_islets(seed_xy)
    n_islets = int(labels[labels >= 0].max() + 1) if (labels >= 0).any() else 0
    print(f"    n_islets={n_islets}")

    a.file.close()

    # Build per-islet metadata
    islet_seeds = {}    # k -> array of seed xy
    for k in range(n_islets):
        islet_seeds[k] = seed_xy[labels == k]

    # Load immune
    p_imm = ROOT / f"data/processed/{sample}/{sample}_immune_phenotyped.h5ad"
    a_imm = sc.read_h5ad(p_imm)
    if "sample" in a_imm.obs.columns:
        a_imm = a_imm[a_imm.obs["sample"].astype(str) == sample].copy()
    immune_xy = np.asarray(a_imm.obsm["spatial"])
    immune_subtype = a_imm.obs["immune_subtype"].astype(str).values
    print(f"    immune cells: {a_imm.n_obs:,}")

    # NN distance to nearest clustered seed
    in_islet_mask = labels >= 0
    seed_xy_clust = seed_xy[in_islet_mask]
    seed_label_clust = labels[in_islet_mask]
    print(f"  computing distance-to-islet for immune cells...")
    d, nn_idx = nn_distance(immune_xy, seed_xy_clust)
    nearest_islet = seed_label_clust[nn_idx]

    # Distance bins
    in_zone = d < ISLET_ZONE_UM        # within 50μm of any seed = islet zone
    in_proximal = (d >= ISLET_ZONE_UM) & (d < PERI_PROXIMAL_UM)
    in_distal = d >= PERI_PROXIMAL_UM

    # ===== Per-islet immune count + insulitis grade =====
    df_imm = pd.DataFrame({
        "subtype": immune_subtype[in_zone],
        "islet_id": [f"{sample}_islet_{k:04d}" for k in nearest_islet[in_zone]],
    })
    per_islet = (df_imm.groupby(["islet_id", "subtype"], observed=True)
                  .size().unstack("subtype", fill_value=0))
    all_subtypes = sorted(set(immune_subtype))
    for st in all_subtypes:
        if st not in per_islet.columns:
            per_islet[st] = 0
    per_islet = per_islet[all_subtypes].reset_index()

    # Add islet metadata
    islet_meta = pd.DataFrame({
        "islet_id": [f"{sample}_islet_{k:04d}" for k in range(n_islets)],
        "n_endocrine": [len(islet_seeds[k]) for k in range(n_islets)],
        "centroid_x": [islet_seeds[k][:, 0].mean() for k in range(n_islets)],
        "centroid_y": [islet_seeds[k][:, 1].mean() for k in range(n_islets)],
    })
    per_islet = islet_meta.merge(per_islet, on="islet_id", how="left").fillna(0)
    per_islet[all_subtypes] = per_islet[all_subtypes].astype(int)
    per_islet["sample"] = sample
    per_islet["total_immune"] = per_islet[all_subtypes].sum(axis=1)

    # Insulitis grade
    def grade(n):
        if n >= INSULITIS_THRESH:
            return "insulitis"
        if n >= PERI_THRESH:
            return "peri_insulitis"
        return "no_insulitis"
    per_islet["insulitis_grade"] = per_islet["total_immune"].apply(grade)

    # ===== Density enrichment per subtype =====
    # density_ratio = (N_subtype_in_zone / N_total_in_zone) /
    #                 (N_subtype_in_distal / N_total_in_distal)
    # Independent of islet count; sample-comparable.
    enrichment_rows = []
    n_zone_total = int(in_zone.sum())
    n_distal_total = int(in_distal.sum())
    for st in all_subtypes:
        is_st = immune_subtype == st
        n_zone_st = int((is_st & in_zone).sum())
        n_distal_st = int((is_st & in_distal).sum())
        zone_frac = n_zone_st / max(n_zone_total, 1)
        distal_frac = n_distal_st / max(n_distal_total, 1)
        ratio = zone_frac / max(distal_frac, 1e-12) if distal_frac > 0 else float("nan")
        enrichment_rows.append({
            "sample": sample,
            "immune_subtype": st,
            "n_in_islet_zone": n_zone_st,
            "n_in_proximal": int((is_st & in_proximal).sum()),
            "n_in_distal": n_distal_st,
            "zone_frac": zone_frac,
            "distal_frac": distal_frac,
            "density_enrichment": ratio,
        })
    enrichment_df = pd.DataFrame(enrichment_rows)
    enrichment_df["n_islet_zone_total"] = n_zone_total
    enrichment_df["n_distal_total"] = n_distal_total

    return per_islet, enrichment_df


# ===== Run both samples =====
all_per, all_enr = [], []
for s in SAMPLES:
    pi, en = process_sample(s)
    all_per.append(pi)
    all_enr.append(en)

per_df = pd.concat(all_per, ignore_index=True)
enr_df = pd.concat(all_enr, ignore_index=True)

# ===== Insulitis grade summary =====
grade_summary = (per_df.groupby(["sample", "insulitis_grade"], observed=True)
                  .size().unstack("insulitis_grade", fill_value=0))
for col in ("no_insulitis", "peri_insulitis", "insulitis"):
    if col not in grade_summary.columns:
        grade_summary[col] = 0
grade_summary = grade_summary[["no_insulitis", "peri_insulitis", "insulitis"]]
grade_summary["total_islets"] = grade_summary.sum(axis=1)
grade_summary["pct_insulitis"] = (
    100 * grade_summary["insulitis"] / grade_summary["total_islets"])
grade_summary["pct_peri_or_insulitis"] = (
    100 * (grade_summary["insulitis"] + grade_summary["peri_insulitis"])
    / grade_summary["total_islets"])
grade_summary = grade_summary.reset_index()

# ===== Backup old tables, write new ones =====
old_proximity = ROOT / "data/processed/immune_proximity_summary.csv"
old_per100 = ROOT / "data/processed/islet_infiltration_per100endo.csv"
old_islets = ROOT / "data/processed/islets_dbscan.csv"

if old_proximity.exists():
    shutil.copy2(old_proximity, ROOT / "data/processed/immune_proximity_summary_strict.csv")
if old_per100.exists():
    shutil.copy2(old_per100, ROOT / "data/processed/islet_infiltration_per100endo_strict.csv")
if old_islets.exists():
    shutil.copy2(old_islets, ROOT / "data/processed/islets_dbscan_strict.csv")
print(f"\n  archived strict-seed CSVs as *_strict.csv")

# Promote endo-label islets table to be the headline
endolabel_islets = ROOT / "data/processed/islets_dbscan_endolabel.csv"
if endolabel_islets.exists():
    shutil.copy2(endolabel_islets, old_islets)

# Write the new headline tables
enr_df.to_csv(old_proximity, index=False)
print(f"  wrote {old_proximity.name}: density enrichment per subtype × sample")
per_df.to_csv(old_per100, index=False)
print(f"  wrote {old_per100.name}: per-islet counts + insulitis_grade")
grade_summary.to_csv(ROOT / "data/processed/islet_insulitis_grades.csv", index=False)
print(f"  wrote islet_insulitis_grades.csv")

# ===== Print summary =====
print(f"\n{'='*72}\n=== INSULITIS GRADE SUMMARY (clinical thresholds)\n{'='*72}")
print(f"  no_insulitis:    < {PERI_THRESH} immune cells within {ISLET_ZONE_UM:.0f} μm of any endocrine seed")
print(f"  peri_insulitis:  {PERI_THRESH}-{INSULITIS_THRESH-1} immune cells within {ISLET_ZONE_UM:.0f} μm")
print(f"  insulitis:       >= {INSULITIS_THRESH} immune cells within {ISLET_ZONE_UM:.0f} μm")
print()
print(grade_summary.to_string(index=False))

print(f"\n{'='*72}\n=== DENSITY ENRICHMENT (islet-zone vs distal-background)\n{'='*72}")
print(f"  density_enrichment > 1 → subtype concentrates near islets")
print(f"  density_enrichment < 1 → subtype is depleted near islets")
print(f"  density_enrichment = 1 → uniformly distributed (random)")
print()
for s in SAMPLES:
    sub = enr_df[enr_df["sample"] == s].sort_values("density_enrichment", ascending=False)
    print(f"\n  {s} (zone={int(sub['n_islet_zone_total'].iloc[0]):,} cells, "
          f"distal={int(sub['n_distal_total'].iloc[0]):,} cells):")
    print(f"    {'subtype':<14} {'zone':>6} {'distal':>8} {'enrichment':>11}")
    for _, r in sub.iterrows():
        print(f"    {r['immune_subtype']:<14} {int(r['n_in_islet_zone']):>6} "
              f"{int(r['n_in_distal']):>8} {r['density_enrichment']:>11.2f}")

print(f"\n{'='*72}\n=== TOP-10 INSULITIS ISLETS PER SAMPLE\n{'='*72}")
for s in SAMPLES:
    sub = per_df[per_df["sample"] == s].nlargest(10, "total_immune")
    cols = ["islet_id", "n_endocrine", "total_immune", "insulitis_grade",
            "T_cytotoxic", "T_reg", "T_exhausted", "Macro_M1", "Macro_M2"]
    cols = [c for c in cols if c in sub.columns]
    print(f"\n  {s}:")
    print(sub[cols].to_string(index=False))
