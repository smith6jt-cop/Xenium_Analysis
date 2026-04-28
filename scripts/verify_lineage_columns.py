"""Verify all new lineage columns landed in both h5ads for both samples."""
import h5py

NEW_COLS = (
    "lineage_status",
    "celltype_lineage",
    "lineage_dominant",
    "lineage_top_depth",
    "lineage_second_depth",
    "lineage_n_groups",
)

for s in ("0041323", "0041326"):
    print(f"\n=== {s} ===")
    for fname in (
        f"data/processed/{s}/{s}_phenotyped.h5ad",
        f"data/processed/{s}/{s}_spatial_analysis.h5ad",
    ):
        with h5py.File(fname, "r") as f:
            print(f"  {fname.split('/')[-1]}:")
            for col in NEW_COLS:
                exists = col in f["obs"]
                marker = "present" if exists else "MISSING"
                print(f"    {col:<22} {marker}")
            # Show a few sample values from celltype_lineage vs celltype
            if "celltype_lineage" in f["obs"]:
                if "categories" in f["obs/celltype_lineage"]:
                    cats = f["obs/celltype_lineage/categories"].asstr()[:]
                    codes = f["obs/celltype_lineage/codes"][:]
                    import numpy as np
                    vc = np.bincount(codes, minlength=len(cats))
                    pairs = sorted(zip(cats, vc), key=lambda kv: -kv[1])
                    print(f"    celltype_lineage value_counts:")
                    for c, n in pairs:
                        print(f"      {c:<22} {n:>10,}")
            if "lineage_status" in f["obs"]:
                if "categories" in f["obs/lineage_status"]:
                    cats = f["obs/lineage_status/categories"].asstr()[:]
                    codes = f["obs/lineage_status/codes"][:]
                    import numpy as np
                    vc = np.bincount(codes, minlength=len(cats))
                    print(f"    lineage_status value_counts:")
                    for c, n in zip(cats, vc):
                        print(f"      {c:<22} {n:>10,}")
