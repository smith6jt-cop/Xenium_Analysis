"""
Comprehensive fix for DEG statistics in the three THYHDL notebooks.

Root causes of extreme p-values (-log10 reaching 300):
1. Subsampling targets ~10k cells. With 6-10 clusters, this means 1000-1700
   cells per cluster. Wilcoxon p-values scale exponentially with sample size —
   at N=1000 per group, even moderate differences give -log10(p) > 60.
2. The Wilcoxon rank-sum test has power ~ sqrt(n). At n=1000+, it flags
   noise-level differences as significant.
3. No minimum expression filter — genes expressed in very few cells contribute
   unreliable results.

Fixes applied:
A. DEG cell: Cap subsampling at 200 cells per cluster (keeps -log10(p) in
   the 10-40 range for genuine markers). Keep the Wilcoxon test (appropriate
   for non-normal count data) but at a sample size where it produces
   interpretable p-values. Add pts=True for downstream filtering.
   Keep the Seurat-style fold change recomputation from normalized counts.
B. Volcano cell: Revert p-value cap (1e-50 -> 1e-300). With proper sample
   sizes, p-values will naturally stay in a reasonable range. Update threshold
   to ±2.
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

NEW_DEG_CELL = """\
# Restore log-normalized data to .X for DE
adata.X = adata.layers['log_normalized'].copy()

# Subsample for DE: max 200 cells per cluster to keep p-values interpretable.
# Wilcoxon p-values scale exponentially with N — at 1000+ cells per group,
# even noise-level differences yield p < 1e-100. Capping at 200 keeps
# -log10(pval) in the 10-40 range for genuine markers.
import numpy as np
np.random.seed(42)

MAX_PER_CLUSTER = 200
MIN_PER_CLUSTER = 50
sampled_indices = []
for group in adata.obs['leiden_0.5'].cat.categories:
    idx = adata.obs.index[adata.obs['leiden_0.5'] == group]
    n_sample = min(MAX_PER_CLUSTER, len(idx))
    n_sample = max(MIN_PER_CLUSTER, n_sample)
    n_sample = min(n_sample, len(idx))  # don't exceed cluster size
    sampled_indices.extend(np.random.choice(idx, size=n_sample, replace=False))

adata_de = adata[sampled_indices].copy()
print(f'Subsampled from {adata.n_obs} to {adata_de.n_obs} cells '
      f'({adata_de.n_obs/adata.n_obs*100:.1f}%) for DE')
print(f'Cluster sizes: {adata_de.obs["leiden_0.5"].value_counts().to_dict()}')

sc.tl.rank_genes_groups(adata_de, groupby='leiden_0.5', method='wilcoxon',
                        use_raw=False, pts=True)
adata.uns['rank_genes_groups'] = adata_de.uns['rank_genes_groups']

# Recompute logfoldchanges from normalized counts (Seurat-style).
# Scanpy uses expm1(mean(log1p(x))) which suffers from Jensen's inequality.
# Correct approach: expm1 per cell first, then average per group.
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

    gene_to_lfc = dict(zip(gene_names, lfc))
    ranked_genes = rgg['names'][group_name]
    rgg['logfoldchanges'][group_name] = np.array(
        [gene_to_lfc[g] for g in ranked_genes], dtype='float32'
    )

print('Recomputed logfoldchanges from normalized counts (pseudocount=1)')
del adata_de"""


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

        # --- Part A: Replace entire DEG cell ---
        deg_cell = nb["cells"][deg_idx]
        deg_cell["source"] = src_to_lines(NEW_DEG_CELL)
        deg_cell["outputs"] = []
        deg_cell["execution_count"] = None
        print(f"  Cell {deg_idx} (DEG): replaced with capped subsampling (200/cluster)")

        # --- Part B: Fix volcano plot ---
        vol_cell = nb["cells"][volcano_idx]
        vol_src = "".join(vol_cell["source"])

        # Ensure lfc_thresh = 2.0
        vol_src = vol_src.replace("lfc_thresh = 1.0", "lfc_thresh = 2.0")

        # Revert p-value cap back to 1e-300 (no artificial cap needed now)
        vol_src = vol_src.replace(
            "df['pvals_adj'].clip(lower=1e-50)",
            "df['pvals_adj'].clip(lower=1e-300)",
        )

        vol_cell["source"] = src_to_lines(vol_src)
        vol_cell["outputs"] = []
        vol_cell["execution_count"] = None
        print(f"  Cell {volcano_idx} (volcano): lfc_thresh=2.0, reverted p-value cap")

        with open(nb_path, "w") as f:
            json.dump(nb, f, indent=1)
            f.write("\n")

    print("\nAll notebooks updated.")


if __name__ == "__main__":
    main()
