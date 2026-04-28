"""Side-by-side comparison: strict hormone-marker seed vs. endo-label seed.

Strict seed (lognorm >= log1p(5) of CHGA/INSM1/ISL1/NEUROD1/FEV):
  data/processed/islets_dbscan.csv
Endo-label seed (celltype_lineage in {Endocrine,Beta,Alpha,Delta,Endocrine_pan}):
  data/processed/islets_dbscan_endolabel.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/blue/maigan/smith6jt/Xenium_Analysis")

isl_strict = pd.read_csv(ROOT / "data/processed/islets_dbscan.csv",
                          dtype={"sample": str})
isl_endo = pd.read_csv(ROOT / "data/processed/islets_dbscan_endolabel.csv",
                        dtype={"sample": str})
ct_strict = pd.read_csv(ROOT / "data/processed/immune_proximity_summary.csv",
                         dtype={"sample": str})
ct_endo = pd.read_csv(ROOT / "data/processed/immune_proximity_summary_endolabel.csv",
                       dtype={"sample": str})
per_strict = pd.read_csv(ROOT / "data/processed/islet_infiltration_per100endo.csv",
                          dtype={"sample": str})
per_endo = pd.read_csv(ROOT / "data/processed/islet_infiltration_per100endo_endolabel.csv",
                        dtype={"sample": str})

print("=" * 78)
print("ISLET SEEDING COMPARISON: strict hormone marker  vs.  endo-label")
print("=" * 78)

print("\n--- Islet counts ---")
print(f"  {'sample':<10} {'strict':>10} {'endo-label':>12} {'ratio':>8}")
for s in ("0041323", "0041326"):
    a = (isl_strict["sample"] == s).sum()
    b = (isl_endo["sample"] == s).sum()
    r = b / a if a > 0 else float("inf")
    print(f"  {s:<10} {a:>10,} {b:>12,} {r:>8.1f}x")

print("\n--- Median seeds / islet (size proxy) ---")
print(f"  {'sample':<10} {'strict':>12} {'endo-label':>14}")
for s in ("0041323", "0041326"):
    a = isl_strict[isl_strict["sample"] == s]["n_seed_cells"].median()
    b = isl_endo[isl_endo["sample"] == s]["n_seed_cells"].median()
    print(f"  {s:<10} {a:>12.0f} {b:>14.0f}")

print("\n--- Approx radius (μm) per islet, median ---")
print(f"  {'sample':<10} {'strict':>12} {'endo-label':>14}")
for s in ("0041323", "0041326"):
    a = isl_strict[isl_strict["sample"] == s]["approx_radius_um"].median()
    b = isl_endo[isl_endo["sample"] == s]["approx_radius_um"].median()
    print(f"  {s:<10} {a:>12.1f} {b:>14.1f}")

print("\n--- Endo-label contingency: subtype × bin × sample ---")
print(ct_endo.to_string(index=False))

print("\n--- Endo-label odds ratio (intra/peri vs distal), 0041323 ---")
sub = ct_endo[ct_endo["sample"] == "0041323"]
print(f"  {'subtype':<14} {'intra/peri':>11} {'distal':>10} {'OR':>8}")
for _, row in sub.iterrows():
    a = row["intra_or_peri_islet"]
    b = row["distal"]
    others_a = sub["intra_or_peri_islet"].sum() - a
    others_b = sub["distal"].sum() - b
    if a + others_a == 0 or b + others_b == 0:
        continue
    odds = (a * others_b) / max((b * others_a), 1)
    print(f"  {row['immune_subtype']:<14} {int(a):>11} {int(b):>10} {odds:>8.2f}")

print("\n--- Endo-label odds ratio (intra/peri vs distal), 0041326 ---")
sub = ct_endo[ct_endo["sample"] == "0041326"]
print(f"  {'subtype':<14} {'intra/peri':>11} {'distal':>10} {'OR':>8}")
for _, row in sub.iterrows():
    a = row["intra_or_peri_islet"]
    b = row["distal"]
    others_a = sub["intra_or_peri_islet"].sum() - a
    others_b = sub["distal"].sum() - b
    if a + others_a == 0 or b + others_b == 0:
        continue
    odds = (a * others_b) / max((b * others_a), 1)
    print(f"  {row['immune_subtype']:<14} {int(a):>11} {int(b):>10} {odds:>8.2f}")

print("\n--- Endo-label: median per-100-endocrine infiltration by sample ---")
per_cols = [c for c in per_endo.columns if c.endswith("_per100endo")]
med = per_endo.groupby("sample", observed=True)[per_cols].median().T
print(med.round(3).to_string())

print("\n--- Endo-label: top-5 most-infiltrated islets per sample ---")
non_subtype = {"n_endocrine", "sample", "islet_id"}
subtype_cols = [c for c in per_endo.columns
                  if c not in non_subtype and not c.endswith("_per100endo")]
per_endo["total_immune"] = per_endo[subtype_cols].sum(axis=1)
for s in ("0041323", "0041326"):
    sub = per_endo[per_endo["sample"] == s].nlargest(5, "total_immune")
    cols = ["islet_id", "n_endocrine", "total_immune"]
    for st in ("T_cytotoxic", "T_reg", "T_exhausted", "Macro_M1", "Macro_M2"):
        if st in sub.columns:
            cols.append(st)
    print(f"\n  {s}:")
    print(sub[cols].to_string(index=False))

print("\n--- Cross-sample sanity: islet recovery normalized to endocrine pool ---")
endo_pool = {"0041323": 36324, "0041326": 32026}
print(f"  {'sample':<10} {'endo cells':>12} {'islets':>9} {'islets/1000 endo':>18}")
for s in ("0041323", "0041326"):
    n = (isl_endo["sample"] == s).sum()
    pool = endo_pool[s]
    print(f"  {s:<10} {pool:>12,} {n:>9,} {1000*n/pool:>18.2f}")
