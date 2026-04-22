"""
Edit notebook 04_group_comparisons.ipynb to add subsample-based parameter
tuning cells for Harmony/UMAP, and update the full integration cell to use
tuned parameters.

Inserts 3 new cells between 'integration-header' and 'integration-code',
then updates 'integration-code' to reference BEST_PARAMS.
"""

import json
import sys
from pathlib import Path

NOTEBOOK_PATH = Path(__file__).resolve().parent.parent / "notebooks" / "04_group_comparisons.ipynb"


def make_source(text: str) -> list[str]:
    """Convert a multi-line string to notebook source format (list of lines with \\n)."""
    lines = text.split("\n")
    result = []
    for i, line in enumerate(lines):
        if i < len(lines) - 1:
            result.append(line + "\n")
        else:
            result.append(line)  # last line: no trailing newline
    return result


def make_markdown_cell(source_str: str, cell_id: str) -> dict:
    return {
        "cell_type": "markdown",
        "id": cell_id,
        "metadata": {},
        "source": make_source(source_str),
    }


def make_code_cell(source_str: str, cell_id: str) -> dict:
    return {
        "cell_type": "code",
        "execution_count": None,
        "id": cell_id,
        "metadata": {},
        "outputs": [],
        "source": make_source(source_str),
    }


# ---------------------------------------------------------------------------
# Cell contents
# ---------------------------------------------------------------------------

MARKDOWN_SUBSAMPLE_HEADER = """\
## 1a. Subsample-Based Parameter Search

Subsample ~75k cells proportionally from each sample, then sweep Harmony + neighbor + UMAP
parameters to find a combination that produces a clean, well-separated UMAP without spikes.
Run the full-scale integration (Section 1b) only after confirming good parameters here."""

SUBSAMPLE_SWEEP_CODE = """\
import gc
import time
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from scipy.sparse.csgraph import connected_components

# ── Configuration ──────────────────────────────────────────────────
SUBSAMPLE_N = 75_000
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

# ── Stratified subsample ──────────────────────────────────────────
sample_counts = adata.obs['sample_id'].value_counts()
sample_fracs  = sample_counts / sample_counts.sum()
sample_n      = (sample_fracs * SUBSAMPLE_N).round().astype(int)
sample_n.iloc[0] += SUBSAMPLE_N - sample_n.sum()  # adjust rounding

idx = []
for sid, n in sample_n.items():
    pool = adata.obs.index[adata.obs['sample_id'] == sid]
    idx.extend(np.random.choice(pool, size=min(n, len(pool)), replace=False))

adata_sub = adata[idx].copy()
print(f"Subsampled {adata_sub.n_obs:,} cells from {adata.n_obs:,}")
print(adata_sub.obs['sample_id'].value_counts())

# ── Preprocess subsample ──────────────────────────────────────────
adata_sub.X = adata_sub.layers['counts'].copy()
sc.pp.highly_variable_genes(
    adata_sub, n_top_genes=2000, flavor='seurat_v3',
    batch_key='sample_id', subset=False
)
sc.pp.normalize_total(adata_sub)
sc.pp.log1p(adata_sub)

adata_hvg_sub = adata_sub[:, adata_sub.var.highly_variable].copy()
sc.pp.scale(adata_hvg_sub, max_value=10)
sc.pp.pca(adata_hvg_sub, n_comps=50, svd_solver='arpack')
adata_sub.obsm['X_pca'] = adata_hvg_sub.obsm['X_pca'].copy()
adata_sub.uns['pca'] = adata_hvg_sub.uns['pca'].copy()
del adata_hvg_sub; gc.collect()

# ── PCA scree plot ────────────────────────────────────────────────
pca_var = adata_sub.uns['pca']['variance_ratio']
cumvar  = np.cumsum(pca_var)

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
axes[0].plot(range(1, len(pca_var)+1), pca_var, 'o-', markersize=3)
axes[0].set_xlabel('PC'); axes[0].set_ylabel('Variance ratio')
axes[0].set_title('PCA Scree Plot')
axes[0].axhline(0.01, ls='--', c='red', alpha=0.5)

axes[1].plot(range(1, len(cumvar)+1), cumvar, 'o-', markersize=3)
axes[1].set_xlabel('PC'); axes[1].set_ylabel('Cumulative variance')
axes[1].set_title('Cumulative Variance Explained')
for thresh in [0.80, 0.90, 0.95]:
    n_pc = int(np.searchsorted(cumvar, thresh)) + 1
    axes[1].axhline(thresh, ls='--', alpha=0.3)
    axes[1].annotate(f'{thresh:.0%} at PC {n_pc}', (n_pc, thresh), fontsize=8, ha='left')

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'param_sweep_pca_scree.png', dpi=150, bbox_inches='tight')
plt.show()
print(f"Variance explained — 20 PCs: {cumvar[19]:.1%}, 30 PCs: {cumvar[29]:.1%}, 50 PCs: {cumvar[49]:.1%}")

# ── Phase 1: Harmony parameter sweep ─────────────────────────────
harmony_configs = [
    {'max_iter_harmony': 10, 'theta': 1.0},   # current defaults
    {'max_iter_harmony': 20, 'theta': 2.0},   # stronger correction
]

print("\\n" + "=" * 60)
print("PHASE 1: Harmony parameter sweep (fixed nn=50, pcs=10, md=0.3)")
print("=" * 60)

harmony_results = {}
for hp in harmony_configs:
    label = f"theta={hp['theta']}, max_iter={hp['max_iter_harmony']}"
    print(f"\\n  Testing: {label}")
    t0 = time.time()

    adata_test = adata_sub.copy()
    sc.external.pp.harmony_integrate(adata_test, key='sample_id', **hp)
    sc.pp.neighbors(adata_test, use_rep='X_pca_harmony', n_neighbors=50, n_pcs=10)

    n_comp, _ = connected_components(adata_test.obsp['connectivities'], directed=False)
    print(f"    Connected components: {n_comp}")

    sc.tl.umap(adata_test, min_dist=0.3, spread=1.0, random_state=42)
    elapsed = time.time() - t0
    print(f"    Time: {elapsed:.1f}s")

    harmony_results[label] = adata_test

fig, axes = plt.subplots(len(harmony_results), 2, figsize=(14, 6 * len(harmony_results)))
if len(harmony_results) == 1:
    axes = axes[np.newaxis, :]
for i, (label, ad) in enumerate(harmony_results.items()):
    sc.pl.umap(ad, color='sample_id', ax=axes[i, 0], show=False, size=1)
    axes[i, 0].set_title(f'{label} — by sample')
    sc.pl.umap(ad, color='total_counts', ax=axes[i, 1], show=False, size=1)
    axes[i, 1].set_title(f'{label} — total counts')
plt.tight_layout()
plt.savefig(FIGURES_DIR / 'param_sweep_harmony.png', dpi=150, bbox_inches='tight')
plt.show()

del harmony_results; gc.collect()

# ── Phase 2: Neighbor + UMAP sweep ───────────────────────────────
# SET THESE based on Phase 1 results (defaults assume theta=2.0 wins)
BEST_HARMONY = {'max_iter_harmony': 20, 'theta': 2.0}

print("\\n" + "=" * 60)
print("PHASE 2: Neighbor + UMAP parameter sweep")
print("=" * 60)

adata_harm = adata_sub.copy()
sc.external.pp.harmony_integrate(adata_harm, key='sample_id', **BEST_HARMONY)

neighbor_configs = [
    {'n_neighbors': 15, 'n_pcs': 10},
    {'n_neighbors': 30, 'n_pcs': 10},
    {'n_neighbors': 50, 'n_pcs': 10},
    {'n_neighbors': 15, 'n_pcs': 15},
    {'n_neighbors': 30, 'n_pcs': 15},
    {'n_neighbors': 50, 'n_pcs': 15},
]
umap_configs = [
    {'min_dist': 0.3, 'spread': 1.0},
    {'min_dist': 0.5, 'spread': 1.0},
]

results = {}
for np_cfg in neighbor_configs:
    for up_cfg in umap_configs:
        label = (f"nn={np_cfg['n_neighbors']}, pcs={np_cfg['n_pcs']}, "
                 f"md={up_cfg['min_dist']}")
        print(f"  {label}", end=" ... ")
        t0 = time.time()

        adata_test = adata_harm.copy()
        sc.pp.neighbors(adata_test, use_rep='X_pca_harmony', **np_cfg)

        n_comp, _ = connected_components(adata_test.obsp['connectivities'], directed=False)
        sc.tl.umap(adata_test, random_state=42, **up_cfg)
        elapsed = time.time() - t0
        print(f"{elapsed:.1f}s, components={n_comp}")

        results[label] = adata_test

# Plot grid
n_combos = len(results)
ncols = 4
nrows = (n_combos + ncols - 1) // ncols
fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 5 * nrows))
axes_flat = axes.flatten()
for i, (label, ad) in enumerate(results.items()):
    sc.pl.umap(ad, color='sample_id', ax=axes_flat[i], show=False,
               size=0.5, title=label, legend_loc='none')
for j in range(i + 1, len(axes_flat)):
    axes_flat[j].axis('off')
plt.tight_layout()
plt.savefig(FIGURES_DIR / 'param_sweep_neighbors_umap.png', dpi=150, bbox_inches='tight')
plt.show()

del results, adata_harm, adata_sub; gc.collect()
print("\\nParameter sweep complete. Review plots and set BEST_PARAMS below.")\
"""

BEST_PARAMS_CODE = """\
# ── BEST PARAMETERS (update after reviewing sweep plots) ─────────
BEST_PARAMS = {
    # Harmony
    'max_iter_harmony': 20,
    'theta': 2.0,
    # Neighbors
    'n_neighbors': 50,
    'n_pcs': 10,       # elbow at ~10 PCs; beyond that is noise
    # UMAP
    'min_dist': 0.3,
    'spread': 1.0,
    # Leiden
    'resolution': 1.0,
}
print("Selected parameters for full-scale run:")
for k, v in BEST_PARAMS.items():
    print(f"  {k}: {v}")\
"""

UPDATED_INTEGRATION_CODE = """\
import gc
from scipy.sparse.csgraph import connected_components

# --- 1b. Full-scale Joint Integration (using tuned parameters) ---
adata.X = adata.layers['counts'].copy()

# HVG selection per batch (uses raw counts with seurat_v3)
sc.pp.highly_variable_genes(
    adata, n_top_genes=2000, flavor='seurat_v3',
    batch_key='sample_id', subset=False
)
print(f'Highly variable genes: {adata.var.highly_variable.sum()}')

# Normalize and log-transform
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.layers['log_normalized'] = adata.X.copy()

# PCA on HVGs only
adata_hvg = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.pp.pca(adata_hvg, n_comps=max(BEST_PARAMS['n_pcs'], 50), svd_solver='arpack')
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca'].copy()
adata.uns['pca'] = adata_hvg.uns['pca'].copy()
del adata_hvg
gc.collect()

# Harmony batch correction with tuned parameters
print('Running Harmony integration...')
sc.external.pp.harmony_integrate(
    adata, key='sample_id',
    max_iter_harmony=BEST_PARAMS['max_iter_harmony'],
    theta=BEST_PARAMS['theta'],
)
print('Harmony done.')

# Neighbors on corrected embedding
sc.pp.neighbors(
    adata, use_rep='X_pca_harmony',
    n_neighbors=BEST_PARAMS['n_neighbors'],
    n_pcs=BEST_PARAMS['n_pcs'],
)

# Connectivity check
n_components, _ = connected_components(adata.obsp['connectivities'], directed=False)
print(f'Neighbor graph connected components: {n_components}')
if n_components > 1:
    print('WARNING: graph has disconnected components — UMAP may show spikes.')
    print('Consider increasing n_neighbors or n_pcs.')

# UMAP and clustering
sc.tl.umap(
    adata,
    min_dist=BEST_PARAMS['min_dist'],
    spread=BEST_PARAMS['spread'],
    random_state=42,
)
sc.tl.leiden(adata, resolution=BEST_PARAMS['resolution'], key_added=CLUSTER_KEY)

print(f'\\nJoint clusters: {adata.obs[CLUSTER_KEY].nunique()}')
print(f'Cells per cluster:')
print(adata.obs[CLUSTER_KEY].value_counts().sort_index())\
"""


def main():
    with open(NOTEBOOK_PATH) as f:
        nb = json.load(f)

    # Find insertion point: between integration-header and integration-code
    header_idx = None
    code_idx = None
    for i, cell in enumerate(nb["cells"]):
        cid = cell.get("id", "")
        if cid == "integration-header":
            header_idx = i
        elif cid == "integration-code":
            code_idx = i

    if header_idx is None or code_idx is None:
        print("ERROR: Could not find integration-header or integration-code cells", file=sys.stderr)
        sys.exit(1)

    # Check if we already inserted (idempotency)
    for cell in nb["cells"]:
        if cell.get("id") == "subsample-header":
            print("Subsample cells already present — removing old ones first.")
            nb["cells"] = [c for c in nb["cells"] if c.get("id") not in
                           ("subsample-header", "subsample-sweep", "best-params")]
            # Re-find indices after removal
            for i, cell in enumerate(nb["cells"]):
                cid = cell.get("id", "")
                if cid == "integration-header":
                    header_idx = i
                elif cid == "integration-code":
                    code_idx = i
            break

    insert_at = header_idx + 1  # right after integration-header

    new_cells = [
        make_markdown_cell(MARKDOWN_SUBSAMPLE_HEADER, "subsample-header"),
        make_code_cell(SUBSAMPLE_SWEEP_CODE, "subsample-sweep"),
        make_code_cell(BEST_PARAMS_CODE, "best-params"),
    ]

    for i, cell in enumerate(new_cells):
        nb["cells"].insert(insert_at + i, cell)

    # Update integration-code cell (now shifted by 3)
    new_code_idx = code_idx + len(new_cells)
    nb["cells"][new_code_idx]["source"] = make_source(UPDATED_INTEGRATION_CODE)
    nb["cells"][new_code_idx]["outputs"] = []
    nb["cells"][new_code_idx]["execution_count"] = None

    with open(NOTEBOOK_PATH, "w") as f:
        json.dump(nb, f, indent=1)

    print(f"Notebook updated: {NOTEBOOK_PATH}")
    print(f"  Inserted 3 cells at index {insert_at}")
    print(f"  Updated integration-code cell at index {new_code_idx}")


if __name__ == "__main__":
    main()
