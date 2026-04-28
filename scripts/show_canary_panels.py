"""Print canary panels + exclusive groups + per-sample availability."""
import sys
from pathlib import Path

import h5py

sys.path.insert(0, "scripts")
from verify_upper_tail_doublets import CANARY_PANELS, EXCLUSIVE_GROUP

print("=" * 72)
print("CURRENT CANARY PANELS (used in scripts/verify_upper_tail_doublets.py)")
print("=" * 72)
print()
print(f"{'Lineage':<16} {'Exclusive group':<18} Markers")
print("-" * 72)
for lin, genes in CANARY_PANELS.items():
    grp = EXCLUSIVE_GROUP[lin]
    print(f"{lin:<16} {grp:<18} {genes}")

print()
print("=" * 72)
print("MUTUAL-EXCLUSION GROUPS")
print("=" * 72)
print("Lineages in the SAME group CAN co-express (legitimate biology).")
print("Lineages in DIFFERENT groups should NOT co-express (= doublet signal).")
print()
groups = {}
for lin, g in EXCLUSIVE_GROUP.items():
    groups.setdefault(g, []).append(lin)
for g, lins in groups.items():
    print(f"  {g:<22} = {lins}")

print()
print("=" * 72)
print("DETECTION RULE")
print("=" * 72)
print("A lineage is 'detected' in a cell if ANY of its markers has lognorm >= log1p(3) ≈ 1.386")
print("(sample-comparable: lognorm = log1p(normalize_total(counts)), per-sample median-rescaled).")
print("A cell is flagged as a doublet if it triggers >= 3 distinct exclusive groups.")

print()
print("=" * 72)
print("PER-SAMPLE AVAILABILITY (markers that survived min_cells=100 filter)")
print("=" * 72)
for sample in ("0041323", "0041326"):
    path = Path(f"data/processed/{sample}/{sample}_preprocessed.h5ad")
    with h5py.File(path, "r") as f:
        if "_index" in f["var"]:
            var_names = f["var/_index"].asstr()[:]
        else:
            attr = f["var"].attrs.get("_index", "_index")
            var_names = f[f"var/{attr}"].asstr()[:]
    panel_set = set(var_names)
    print(f"\n  {sample} (panel size {len(var_names):,} genes):")
    for lin, genes in CANARY_PANELS.items():
        present = [g for g in genes if g in panel_set]
        missing = [g for g in genes if g not in panel_set]
        ok = "OK" if present else "EMPTY!"
        print(f"    {lin:<16} {ok:<8} present={present}  missing={missing}")
