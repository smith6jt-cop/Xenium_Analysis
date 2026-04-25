"""
Fix logfoldchanges and volcano plot display in the three THYHDL notebooks.

Part A: Replace the fold change recomputation in DEG cells to use normalized
counts (expm1 per cell, then mean) instead of raw counts. This is the
Seurat-style approach and avoids both Jensen's inequality and compression
from low raw counts.

Part B: Update volcano plots — LFC threshold to ±2, cap -log10(pval) at 50.
"""

import json
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

# (path, deg_cell, volcano_cell)
NOTEBOOKS = [
    (REPO / "notebooks" / "01_preprocessing_v2_THYHDL065.ipynb", 33, 37),
    (REPO / "notebooks" / "01_preprocessing_v2_THYHDL073.ipynb", 34, 38),
    (REPO / "notebooks" / "01_preprocessing_v2_THYHDL172.ipynb", 34, 38),
]

# Old block to find and replace in the DEG cell
OLD_LFC_BLOCK = """# Recompute logfoldchanges from raw counts to avoid Jensen's inequality artifact
# (scanpy uses expm1(mean(log1p(x))) which produces extreme values for sparse genes)
from scipy.sparse import issparse
counts_de = adata_de.layers['counts']
clusters_de = adata_de.obs['leiden_0.5'].astype(str)
gene_names = adata_de.var_names
rgg = adata.uns['rank_genes_groups']

for group_name in rgg['names'].dtype.names:
    mask_in = (clusters_de == group_name).values
    mask_out = ~mask_in

    if issparse(counts_de):
        mean_in = np.asarray(counts_de[mask_in].mean(axis=0)).flatten()
        mean_out = np.asarray(counts_de[mask_out].mean(axis=0)).flatten()
    else:
        mean_in = counts_de[mask_in].mean(axis=0)
        mean_out = counts_de[mask_out].mean(axis=0)

    # log2 fold change with pseudocount of 1
    lfc = np.log2(mean_in + 1) - np.log2(mean_out + 1)

    # Map gene names in ranked order to the new fold changes
    gene_to_lfc = dict(zip(gene_names, lfc))
    ranked_genes = rgg['names'][group_name]
    rgg['logfoldchanges'][group_name] = np.array(
        [gene_to_lfc[g] for g in ranked_genes], dtype='float32'
    )

print('Recomputed logfoldchanges from raw counts (pseudocount=1)')"""

NEW_LFC_BLOCK = """# Recompute logfoldchanges from normalized counts (Seurat-style)
# Back-transform log-normalized data per cell via expm1, then average per group.
# This avoids Jensen's inequality (scanpy's expm1-of-mean) and compression from
# low raw counts.
from scipy.sparse import issparse
log_norm = adata_de.layers['log_normalized']
if issparse(log_norm):
    norm_counts = log_norm.expm1()
else:
    norm_counts = np.expm1(log_norm)

clusters_de = adata_de.obs['leiden_0.5'].astype(str)
gene_names = adata_de.var_names
rgg = adata.uns['rank_genes_groups']

for group_name in rgg['names'].dtype.names:
    mask_in = (clusters_de == group_name).values
    mask_out = ~mask_in

    mean_in = np.asarray(norm_counts[mask_in].mean(axis=0)).flatten()
    mean_out = np.asarray(norm_counts[mask_out].mean(axis=0)).flatten()

    # log2 fold change with pseudocount of 1
    lfc = np.log2(mean_in + 1) - np.log2(mean_out + 1)

    # Map gene names in ranked order to the new fold changes
    gene_to_lfc = dict(zip(gene_names, lfc))
    ranked_genes = rgg['names'][group_name]
    rgg['logfoldchanges'][group_name] = np.array(
        [gene_to_lfc[g] for g in ranked_genes], dtype='float32'
    )

print('Recomputed logfoldchanges from normalized counts (pseudocount=1)')"""


def src_to_lines(src):
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

        # --- Part A: Fix DEG fold change computation ---
        deg_cell = nb["cells"][deg_idx]
        deg_src = "".join(deg_cell["source"])

        assert OLD_LFC_BLOCK in deg_src, (
            f"Could not find old LFC block in cell {deg_idx} of {nb_path.name}"
        )
        deg_src = deg_src.replace(OLD_LFC_BLOCK, NEW_LFC_BLOCK)
        deg_cell["source"] = src_to_lines(deg_src)
        deg_cell["outputs"] = []
        deg_cell["execution_count"] = None
        print(f"  Cell {deg_idx} (DEG): replaced fold change computation")

        # --- Part B: Fix volcano plot thresholds ---
        vol_cell = nb["cells"][volcano_idx]
        vol_src = "".join(vol_cell["source"])

        # LFC threshold: 1.0 -> 2.0
        vol_src = vol_src.replace("lfc_thresh = 1.0", "lfc_thresh = 2.0")

        # P-value cap: 1e-300 -> 1e-50
        vol_src = vol_src.replace(
            "df['pvals_adj'].clip(lower=1e-300)",
            "df['pvals_adj'].clip(lower=1e-50)",
        )

        vol_cell["source"] = src_to_lines(vol_src)
        vol_cell["outputs"] = []
        vol_cell["execution_count"] = None
        print(f"  Cell {volcano_idx} (volcano): lfc_thresh=2.0, pval cap=1e-50")

        with open(nb_path, "w") as f:
            json.dump(nb, f, indent=1)
            f.write("\n")

    print("\nAll notebooks updated.")


if __name__ == "__main__":
    main()
