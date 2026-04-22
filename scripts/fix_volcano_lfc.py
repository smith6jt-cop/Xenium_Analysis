"""
Fix volcano plot LFC thresholds and pseudocount in three THYHDL notebooks.

Root causes of empty volcano plots (zero significant genes):
1. Pseudocount of 1 in the custom LFC formula compresses fold changes ~50%
   when mean expression values are near 1 (typical for Xenium post-normalization).
   Using 1e-2 gives <2% compression on expressed genes while bounding edge cases
   (zero-expression genes) to ~|LFC| <= 7.
2. lfc_thresh=2.0 is too stringent for Xenium targeted panels (~280 genes)
   in thymus tissue with relatively homogeneous T-cell populations.

Fixes:
A. DEG cell: Change pseudocount from 1 to 1e-2
B. Volcano cell: Change lfc_thresh from 2.0 to 0.5 (standard for Xenium data)
"""

import json
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

NOTEBOOKS = [
    (REPO / "notebooks" / "01_preprocessing_v2_THYHDL065.ipynb", 33, 37),
    (REPO / "notebooks" / "01_preprocessing_v2_THYHDL073.ipynb", 34, 38),
    (REPO / "notebooks" / "01_preprocessing_v2_THYHDL172.ipynb", 34, 38),
]

DEG_REPLACEMENTS = [
    (
        "# log2 fold change with pseudocount of 1\n",
        "# log2 fold change with pseudocount of 1e-2 (<2% compression on expressed\n"
        "    # genes; pseudocount=1 compresses LFC ~50% when mean expression is near 1)\n",
    ),
    (
        "np.log2(mean_in + 1) - np.log2(mean_out + 1)",
        "np.log2(mean_in + 1e-2) - np.log2(mean_out + 1e-2)",
    ),
    (
        "Recomputed logfoldchanges from normalized counts (pseudocount=1)",
        "Recomputed logfoldchanges from normalized counts (pseudocount=1e-2)",
    ),
]

VOLCANO_REPLACEMENTS = [
    ("lfc_thresh = 2.0", "lfc_thresh = 0.5"),
]


def apply_replacements(src, replacements, cell_idx, nb_name):
    for old, new in replacements:
        assert old in src, (
            f"Could not find {old!r} in cell {cell_idx} of {nb_name}"
        )
        src = src.replace(old, new)
    return src


def src_to_lines(src):
    """Convert a single string back to the list-of-lines format notebooks use."""
    lines = src.split("\n")
    return [
        line + "\n" if i < len(lines) - 1 else line
        for i, line in enumerate(lines)
    ]


def main():
    for nb_path, deg_idx, volcano_idx in NOTEBOOKS:
        print(f"\nEditing {nb_path.name} ...")

        with open(nb_path) as f:
            nb = json.load(f)

        # --- Fix A: DEG cell pseudocount ---
        deg_cell = nb["cells"][deg_idx]
        deg_src = "".join(deg_cell["source"])
        deg_src = apply_replacements(deg_src, DEG_REPLACEMENTS, deg_idx, nb_path.name)
        deg_cell["source"] = src_to_lines(deg_src)
        deg_cell["outputs"] = []
        deg_cell["execution_count"] = None
        print(f"  Cell {deg_idx} (DEG): pseudocount 1 -> 1e-2")

        # --- Fix B: Volcano cell threshold ---
        vol_cell = nb["cells"][volcano_idx]
        vol_src = "".join(vol_cell["source"])
        vol_src = apply_replacements(vol_src, VOLCANO_REPLACEMENTS, volcano_idx, nb_path.name)
        vol_cell["source"] = src_to_lines(vol_src)
        vol_cell["outputs"] = []
        vol_cell["execution_count"] = None
        print(f"  Cell {volcano_idx} (volcano): lfc_thresh 2.0 -> 0.5")

        with open(nb_path, "w") as f:
            json.dump(nb, f, indent=1)
            f.write("\n")

    print("\nDone. Re-run DEG and volcano cells in each notebook to regenerate.")


if __name__ == "__main__":
    main()
