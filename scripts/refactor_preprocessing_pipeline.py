"""
Refactor preprocessing notebooks to:
1. Separate normalization tracks (HVG vs clustering)
2. Add sc.pp.scale() on HVG subset before PCA
3. Change target_sum=None to target_sum=1e4
4. Add adaptive parameters based on dataset size
5. Replace hardcoded leiden keys with LEIDEN_KEY variable
6. Update DE cell to use adaptive parameters
"""

import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

NOTEBOOKS = [
    {
        "path": REPO / "notebooks" / "01_preprocessing_v2_THYHDL065.ipynb",
        "old_leiden_key": "leiden_0.5",
        "de_cell": 33,
        "subsample_cell": 26,
        "summary_cell": 46,
    },
    {
        "path": REPO / "notebooks" / "01_preprocessing_v2_THYHDL073.ipynb",
        "old_leiden_key": "leiden_0.5",
        "de_cell": 34,
        "subsample_cell": 26,
        "summary_cell": 47,
    },
    {
        "path": REPO / "notebooks" / "01_preprocessing_v2_THYHDL172.ipynb",
        "old_leiden_key": "leiden_0.2",
        "de_cell": 34,
        "subsample_cell": 26,
        "summary_cell": 47,
    },
]


def src_to_lines(text: str) -> list:
    """Convert multi-line string to notebook source format (list of lines)."""
    lines = text.split("\n")
    return [line + "\n" if i < len(lines) - 1 else line
            for i, line in enumerate(lines)]


def make_code_cell(source_str: str) -> dict:
    """Create a new code cell dict from source string."""
    return {
        "cell_type": "code",
        "execution_count": None,
        "id": "",
        "metadata": {},
        "outputs": [],
        "source": src_to_lines(source_str),
    }


def clear_cell(cell: dict) -> None:
    """Clear outputs and execution count."""
    cell["outputs"] = []
    cell["execution_count"] = None


# ── New cell contents ────────────────────────────────────────────────────

CELL_18_NORM_HVG = """\
# --- Track A: Save counts, normalize, log-transform, select HVGs ---
adata.layers['counts'] = adata.X.copy()

# Normalize to 10,000 counts per cell (standard CPM-like normalization)
sc.pp.normalize_total(adata, target_sum=1e4, exclude_highly_expressed=False, inplace=True)
sc.pp.log1p(adata)
adata.layers['log_normalized'] = adata.X.copy()

print(f"Dataset: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
print(f"Layers: {list(adata.layers.keys())}")

# HVG selection on raw counts (seurat_v3 requires counts, not normalized)
n_top_genes = min(300, max(100, int(adata.n_vars * 0.065)))
sc.pp.highly_variable_genes(adata, layer='counts', n_top_genes=n_top_genes,
                            flavor='seurat_v3', subset=False)
print(f'Highly variable genes: {adata.var.highly_variable.sum()} / {adata.n_vars} '
      f'({adata.var.highly_variable.sum()/adata.n_vars*100:.1f}%)')\
"""

CELL_19_ADAPTIVE_PARAMS = """\
# --- Adaptive parameters based on dataset size ---
n_cells = adata.n_obs
n_hvgs = int(adata.var.highly_variable.sum())

# PCA: n_comps capped by n_hvgs - 1
n_comps = min(50, n_hvgs - 1)

# Leiden resolution: scale with sqrt(n_cells), base 0.5 at ~80K cells
LEIDEN_RESOLUTION = round(0.5 * np.sqrt(n_cells / 80000), 2)
LEIDEN_RESOLUTION = max(0.1, min(1.5, LEIDEN_RESOLUTION))
LEIDEN_KEY = f'leiden_{LEIDEN_RESOLUTION}'

# DE subsampling: lower floor for smaller datasets
MAX_PER_CLUSTER = 200
MIN_PER_CLUSTER = 25 if n_cells < 50000 else 50

# Co-occurrence subsample: larger datasets get smaller fractions
SUBSAMPLE_FRACTION = min(0.5, max(0.1, 40000 / n_cells))

print(f"=== Adaptive Parameters for {n_cells:,} cells ===")
print(f"  PCA components: {n_comps}")
print(f"  Leiden resolution: {LEIDEN_RESOLUTION} (key: '{LEIDEN_KEY}')")
print(f"  DE subsampling: max {MAX_PER_CLUSTER}, min {MIN_PER_CLUSTER} per cluster")
print(f"  Co-occurrence subsample: {SUBSAMPLE_FRACTION:.2f}")\
"""

CELL_20_PCA = """\
# --- Track B: Scale HVG subset -> PCA (preserves .X as log-normalized) ---
adata_hvg = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.pp.pca(adata_hvg, n_comps=n_comps, svd_solver='arpack')

# Copy PCA results back to full object
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca'].copy()
adata.uns['pca'] = adata_hvg.uns['pca'].copy()
adata.varm['PCs'] = np.zeros((adata.n_vars, n_comps))
adata.varm['PCs'][adata.var.highly_variable.values] = adata_hvg.varm['PCs'].copy()

del adata_hvg
gc.collect()

print(f"PCA: {n_comps} components on {adata.var.highly_variable.sum()} scaled HVGs")
print(f".X unchanged: log-normalized data")\
"""

CELL_21_CLUSTERING = """\
# --- Clustering pipeline ---
sc.pp.neighbors(adata, metric='cosine', n_pcs=min(n_comps, 20))
sc.tl.umap(adata, min_dist=0.06, spread=2.2)
sc.tl.leiden(adata, resolution=LEIDEN_RESOLUTION, key_added=LEIDEN_KEY)

n_clusters = adata.obs[LEIDEN_KEY].nunique()
cluster_sizes = adata.obs[LEIDEN_KEY].value_counts().sort_index()
print(f"\\nClustering: {n_clusters} clusters at resolution {LEIDEN_RESOLUTION}")
for cl, count in cluster_sizes.items():
    flag = " <-- small cluster" if count < MIN_PER_CLUSTER else ""
    print(f"  Cluster {cl}: {count:,}{flag}")\
"""


def replace_leiden_in_line(line: str, old_key: str) -> str:
    """Replace hardcoded leiden key with LEIDEN_KEY variable in a source line.

    Handles patterns like:
      'leiden_0.5'  -> LEIDEN_KEY
      "leiden_0.5"  -> LEIDEN_KEY
      ["leiden_0.5"] -> [LEIDEN_KEY]
      f'...{adata_de.obs["leiden_0.5"]...' -> f'...{adata_de.obs[LEIDEN_KEY]...'
    """
    if old_key not in line:
        return line

    # Pattern: f-string with ["leiden_0.5"] inside braces - replace quotes within
    # e.g. f'...{adata_de.obs["leiden_0.5"]...}'
    line = line.replace(f'["{old_key}"]', '[LEIDEN_KEY]')
    line = line.replace(f"['{old_key}']", '[LEIDEN_KEY]')

    # Pattern: keyword=["leiden_0.5"] -> keyword=[LEIDEN_KEY]
    line = line.replace(f'["{old_key}"]', '[LEIDEN_KEY]')
    line = line.replace(f"['{old_key}']", '[LEIDEN_KEY]')

    # Pattern: keyword='leiden_0.5' -> keyword=LEIDEN_KEY
    line = line.replace(f"'{old_key}'", 'LEIDEN_KEY')
    line = line.replace(f'"{old_key}"', 'LEIDEN_KEY')

    return line


def update_de_cell(cell: dict, old_key: str) -> None:
    """Update the DE cell to use adaptive parameters and LEIDEN_KEY."""
    src = "".join(cell["source"])

    # Remove the .X restoration line (no longer needed)
    src = src.replace(
        "# Restore log-normalized data to .X for DE\n"
        "adata.X = adata.layers['log_normalized'].copy()\n\n",
        ""
    )

    # Remove hardcoded MAX/MIN_PER_CLUSTER definitions (use from Cell 19)
    src = src.replace("MAX_PER_CLUSTER = 200\n", "")
    src = src.replace("MIN_PER_CLUSTER = 50\n", "")

    # Remove the now-redundant numpy import (already in Cell 1)
    src = src.replace("import numpy as np\n", "")

    # Replace leiden key references
    lines = src.split("\n")
    new_lines = [replace_leiden_in_line(l, old_key) for l in lines]
    src = "\n".join(new_lines)

    cell["source"] = src_to_lines(src)
    clear_cell(cell)


def update_subsample_cell(cell: dict) -> None:
    """Update subsample cell to use SUBSAMPLE_FRACTION."""
    src = "".join(cell["source"])
    src = src.replace("fraction=0.5", "fraction=SUBSAMPLE_FRACTION")
    cell["source"] = src_to_lines(src)
    clear_cell(cell)


def add_gc_import(cell: dict) -> None:
    """Add 'import gc' to the imports cell if not present."""
    src = "".join(cell["source"])
    if "import gc" not in src:
        # Add after the last import line
        lines = src.split("\n")
        last_import_idx = 0
        for i, line in enumerate(lines):
            if line.startswith("import ") or line.startswith("from "):
                last_import_idx = i
        lines.insert(last_import_idx + 1, "import gc")
        cell["source"] = src_to_lines("\n".join(lines))
        clear_cell(cell)


def main():
    for config in NOTEBOOKS:
        nb_path = config["path"]
        old_key = config["old_leiden_key"]
        de_cell_idx = config["de_cell"]
        subsample_cell_idx = config["subsample_cell"]
        summary_cell_idx = config["summary_cell"]

        print(f"\n{'='*60}")
        print(f"Editing {nb_path.name} (old key: {old_key})")
        print(f"{'='*60}")

        with open(nb_path) as f:
            nb = json.load(f)

        n_cells_before = len(nb["cells"])
        print(f"  Cells before: {n_cells_before}")

        # --- Step 1: Add gc import to Cell 1 ---
        add_gc_import(nb["cells"][1])
        print("  Added 'import gc' to imports cell")

        # --- Step 2: Replace Cell 18 with normalization + HVG ---
        nb["cells"][18] = make_code_cell(CELL_18_NORM_HVG)
        print("  Replaced Cell 18: normalization + HVG selection")

        # --- Step 3: Insert 3 new cells after Cell 18 ---
        new_cells = [
            make_code_cell(CELL_19_ADAPTIVE_PARAMS),
            make_code_cell(CELL_20_PCA),
            make_code_cell(CELL_21_CLUSTERING),
        ]
        for i, cell in enumerate(new_cells):
            nb["cells"].insert(19 + i, cell)
        print("  Inserted 3 new cells (19: params, 20: PCA, 21: clustering)")

        # Cell indices have shifted by +3 for everything after 18
        de_cell_idx_new = de_cell_idx + 3
        subsample_cell_idx_new = subsample_cell_idx + 3
        summary_cell_idx_new = summary_cell_idx + 3

        # --- Step 4: Update DE cell ---
        update_de_cell(nb["cells"][de_cell_idx_new], old_key)
        print(f"  Updated DE cell (now at index {de_cell_idx_new})")

        # --- Step 5: Update subsample cell ---
        update_subsample_cell(nb["cells"][subsample_cell_idx_new])
        print(f"  Updated subsample cell (now at index {subsample_cell_idx_new})")

        # --- Step 6: Replace leiden key in ALL downstream cells (after Cell 21) ---
        replacements = 0
        for i in range(22, len(nb["cells"])):
            if i == de_cell_idx_new:
                continue  # already handled
            if i == subsample_cell_idx_new:
                continue  # already handled
            new_source = []
            for line in nb["cells"][i]["source"]:
                new_line = replace_leiden_in_line(line, old_key)
                if new_line != line:
                    replacements += 1
                new_source.append(new_line)
            if new_source != nb["cells"][i]["source"]:
                nb["cells"][i]["source"] = new_source
                clear_cell(nb["cells"][i])
            nb["cells"][i]["source"] = new_source
        print(f"  Replaced {replacements} leiden key references in downstream cells")

        n_cells_after = len(nb["cells"])
        print(f"  Cells after: {n_cells_after} (+{n_cells_after - n_cells_before})")

        # --- Write back ---
        with open(nb_path, "w") as f:
            json.dump(nb, f, indent=1)
            f.write("\n")

        print(f"  Saved: {nb_path.name}")

    print(f"\n{'='*60}")
    print("All notebooks updated.")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
