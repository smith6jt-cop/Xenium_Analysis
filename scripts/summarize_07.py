"""Concise summary of immune phenotyping + proximity results from notebook 07."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/blue/maigan/smith6jt/Xenium_Analysis")
NB = ROOT / "notebooks/executed/07/07_combined.ipynb"

print("=" * 72)
print("Immune phenotyping + islet proximity — summary")
print("=" * 72)

nb = json.load(NB.open())

print("\n=== Immune subtype distribution (per cell, both samples combined) ===")
for i, c in enumerate(nb["cells"]):
    if c["cell_type"] != "code":
        continue
    for o in c.get("outputs", []):
        if o.get("output_type") == "stream":
            txt = "".join(o.get("text", []))
            if "per-cell immune_subtype distribution" in txt:
                lines = txt.splitlines()
                for j, line in enumerate(lines):
                    if "per-cell immune_subtype" in line:
                        for k in range(j, min(j + 18, len(lines))):
                            print(f"  {lines[k]}")
                        break
                break

df = pd.read_csv(ROOT / "data/processed/islet_infiltration_per100endo.csv",
                  dtype={"sample": str})
print(f"\n=== Per-islet infiltration table ===")
print(f"  total islets: {len(df):,}")
sample_counts = df["sample"].value_counts().to_dict()
print(f"  per-sample islet count: {sample_counts}")

per_cols = [c for c in df.columns if c.endswith("_per100endo")]
print("\n  median infiltration per 100 endocrine seed cells (per islet):")
med = df.groupby("sample", observed=True)[per_cols].median().T
print(med.round(3).to_string())

# Total immune per islet (sum of all subtypes within proximal+ of that islet)
non_subtype_cols = {"n_endocrine", "sample", "islet_id"}
subtype_cols = [c for c in df.columns
                 if c not in non_subtype_cols and not c.endswith("_per100endo")]
df["total_immune"] = df[subtype_cols].sum(axis=1)
print("\n=== Top-5 most-infiltrated islets per sample ===")
for s in df["sample"].unique():
    sub = df[df["sample"] == s].nlargest(5, "total_immune")
    cols_to_show = ["islet_id", "n_endocrine", "total_immune"]
    for st in ("T_cytotoxic", "T_reg", "T_exhausted", "Macro_M1", "Macro_M2"):
        if st in sub.columns:
            cols_to_show.append(st)
    print(f"\n  {s}:")
    print(sub[cols_to_show].to_string(index=False))

print("\n=== Contingency: subtype × distance_bin × sample ===")
ct = pd.read_csv(ROOT / "data/processed/immune_proximity_summary.csv",
                  dtype={"sample": str})
print(ct.to_string(index=False))

print("\n=== Per-subtype intra/peri vs distal odds ratio (0041323 only — better statistics) ===")
sub_323 = ct[ct["sample"] == "0041323"]
print(f"  {'subtype':<14} {'intra/peri':>11} {'distal':>10} {'OR':>8}")
for _, row in sub_323.iterrows():
    a = row["intra_or_peri_islet"]
    b = row["distal"]
    others_a = sub_323["intra_or_peri_islet"].sum() - a
    others_b = sub_323["distal"].sum() - b
    if a + others_a == 0 or b + others_b == 0:
        continue
    odds = (a * others_b) / max((b * others_a), 1)
    print(f"  {row['immune_subtype']:<14} {int(a):>11} {int(b):>10} {odds:>8.2f}")
