"""
Update preprocessing parameters in the three THYHDL notebooks:
1. normalize_total: add target_sum=None
2. highly_variable_genes: n_top_genes=300
3. Remove scaling entirely (sc.pp.scale + scaled layer)
4. neighbors: add metric='cosine'
5. Downstream cells: layer='scaled' -> layer='log_normalized'
"""

import json
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

# (path, cell_18_leiden_key, dotplot_cell, matrixplot_cell)
NOTEBOOKS = [
    (REPO / "notebooks" / "01_preprocessing_v2_THYHDL065.ipynb", "leiden_1.5", 32, 34),
    (REPO / "notebooks" / "01_preprocessing_v2_THYHDL073.ipynb", "leiden_1.0", 33, 35),
    (REPO / "notebooks" / "01_preprocessing_v2_THYHDL172.ipynb", "leiden_1.5", 33, 35),
]


def replace_in_cell(cell, old, new):
    """Replace a string in a cell's source, returning True if found."""
    src = "".join(cell["source"])
    if old not in src:
        return False
    src = src.replace(old, new)
    # Convert back to line list
    lines = src.split("\n")
    result = []
    for i, line in enumerate(lines):
        if i < len(lines) - 1:
            result.append(line + "\n")
        else:
            result.append(line)
    cell["source"] = result
    return True


def remove_lines_containing(cell, patterns):
    """Remove source lines that contain any of the given patterns."""
    cell["source"] = [
        line for line in cell["source"]
        if not any(p in line for p in patterns)
    ]


def main():
    for nb_path, leiden_key, dotplot_idx, matrixplot_idx in NOTEBOOKS:
        print(f"\nEditing {nb_path.name} ...")

        with open(nb_path) as f:
            nb = json.load(f)

        # --- Cell 18: preprocessing pipeline ---
        cell18 = nb["cells"][18]
        assert cell18["cell_type"] == "code"

        # 1. normalize_total: add target_sum=None
        ok = replace_in_cell(
            cell18,
            "sc.pp.normalize_total(adata, exclude_highly_expressed=False, inplace=True)",
            "sc.pp.normalize_total(adata, target_sum=None, exclude_highly_expressed=False, inplace=True)",
        )
        print(f"  normalize_total target_sum=None: {'OK' if ok else 'SKIPPED (already changed?)'}")

        # 2. highly_variable_genes: n_top_genes=300
        src = "".join(cell18["source"])
        for old_n in ["n_top_genes=500", "n_top_genes=2000"]:
            if old_n in src:
                ok = replace_in_cell(cell18, old_n, "n_top_genes=300")
                print(f"  n_top_genes {old_n} -> 300: {'OK' if ok else 'FAILED'}")
                break
        else:
            print("  n_top_genes: SKIPPED (already 300?)")

        # 3. Remove scaling
        before = len(cell18["source"])
        remove_lines_containing(cell18, [
            "sc.pp.scale(",
            "adata.layers['scaled']",
        ])
        removed = before - len(cell18["source"])
        print(f"  Remove scaling: removed {removed} lines")

        # 4. neighbors: add metric='cosine'
        ok = replace_in_cell(
            cell18,
            "sc.pp.neighbors(adata)",
            "sc.pp.neighbors(adata, metric='cosine')",
        )
        print(f"  neighbors metric='cosine': {'OK' if ok else 'SKIPPED (already changed?)'}")

        cell18["outputs"] = []
        cell18["execution_count"] = None

        # --- Downstream cells: layer='scaled' -> layer='log_normalized' ---
        for ci in [dotplot_idx, matrixplot_idx]:
            cell = nb["cells"][ci]
            ok = replace_in_cell(cell, "layer='scaled'", "layer='log_normalized'")
            print(f"  Cell {ci} layer='scaled' -> 'log_normalized': {'OK' if ok else 'SKIPPED'}")
            cell["outputs"] = []
            cell["execution_count"] = None

        with open(nb_path, "w") as f:
            json.dump(nb, f, indent=1)
            f.write("\n")

    print("\nAll notebooks updated.")


if __name__ == "__main__":
    main()
