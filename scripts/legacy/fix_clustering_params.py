"""
Fix clustering after preprocessing changes:
1. Reduce PCA n_comps to 20 (300 HVGs don't support 50 PCs)
2. Lower Leiden resolution to 0.5
3. Update leiden key references throughout each notebook
"""

import json
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

# (path, old_leiden_key, new_leiden_key)
NOTEBOOKS = [
    (REPO / "notebooks" / "01_preprocessing_v2_THYHDL065.ipynb", "leiden_1.5", "leiden_0.5"),
    (REPO / "notebooks" / "01_preprocessing_v2_THYHDL073.ipynb", "leiden_1.0", "leiden_0.5"),
    (REPO / "notebooks" / "01_preprocessing_v2_THYHDL172.ipynb", "leiden_1.5", "leiden_0.5"),
]

NEW_RESOLUTION = 0.5
NEW_N_COMPS = 20


def main():
    for nb_path, old_key, new_key in NOTEBOOKS:
        print(f"\nEditing {nb_path.name} ({old_key} -> {new_key}) ...")

        with open(nb_path) as f:
            nb = json.load(f)

        # --- Cell 18: PCA n_comps and leiden resolution ---
        cell18 = nb["cells"][18]
        src = "".join(cell18["source"])

        # Add n_comps to PCA
        src = src.replace(
            "sc.pp.pca(adata, svd_solver='arpack', use_highly_variable=True)",
            f"sc.pp.pca(adata, n_comps={NEW_N_COMPS}, svd_solver='arpack', use_highly_variable=True)",
        )

        # Update leiden call
        old_leiden = f"sc.tl.leiden(adata, resolution={old_key.split('_')[1]}, key_added='{old_key}')"
        new_leiden = f"sc.tl.leiden(adata, resolution={NEW_RESOLUTION}, key_added='{new_key}')"
        src = src.replace(old_leiden, new_leiden)

        # Convert back to source lines
        lines = src.split("\n")
        cell18["source"] = [
            line + "\n" if i < len(lines) - 1 else line
            for i, line in enumerate(lines)
        ]
        cell18["outputs"] = []
        cell18["execution_count"] = None

        # --- All cells: replace old leiden key with new ---
        replacements = 0
        for i, cell in enumerate(nb["cells"]):
            if i == 18:
                continue  # already handled
            new_source = []
            for line in cell["source"]:
                if old_key in line:
                    new_source.append(line.replace(old_key, new_key))
                    replacements += 1
                else:
                    new_source.append(line)
            cell["source"] = new_source

        print(f"  Cell 18: PCA n_comps={NEW_N_COMPS}, resolution={NEW_RESOLUTION}")
        print(f"  Replaced {replacements} leiden key references in other cells")

        with open(nb_path, "w") as f:
            json.dump(nb, f, indent=1)
            f.write("\n")

    print("\nAll notebooks updated.")


if __name__ == "__main__":
    main()
