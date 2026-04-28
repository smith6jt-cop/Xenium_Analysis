"""One-shot notebook patcher for the 2026-04-28 reanalysis.

Edits THREE cells across two notebooks:
  notebooks/01_preprocessing_v2.ipynb
    cell 18 — apply panel_for_embedding mask to HVG/scale/PCA/neighbors

  notebooks/02_phenotyping.ipynb
    cell 10 — drop 'Epsilon': ['GHRL'] (1-gene panel that bit us);
               apply panel_for_scoring mask; expand Alpha/Delta panels
               where additional markers are available; warn if any panel
               drops below 2 genes after exclusion.
    cell 17 — add categorical_covariate_keys=['donor_sex'] to scVI;
               train on panel_for_embedding view to drop Y-chromosome
               variance from the latent.

Backups saved as *.bak.preEpsilonFix next to each notebook.

This script is idempotent: re-running detects the patched form and is a no-op.
"""
import json
import shutil
import sys
from pathlib import Path

ROOT = Path("/blue/maigan/smith6jt/Xenium_Analysis")
NB01 = ROOT / "notebooks/01_preprocessing_v2.ipynb"
NB02 = ROOT / "notebooks/02_phenotyping.ipynb"
BACKUP_SUFFIX = ".bak.preEpsilonFix"


def _src(cell):
    return "".join(cell.get("source", []))


def _set_src(cell, text):
    """Assign source as a list of lines with trailing \n (notebook convention)."""
    lines = text.splitlines(keepends=True)
    if lines and not lines[-1].endswith("\n"):
        lines[-1] += "\n"
    cell["source"] = lines


def patch_nb01_cell18(cell):
    src = _src(cell)
    if "panel_for_embedding" in src:
        print("    [01 cell 18] already patched, skipping")
        return False
    new = """# Counts layer is required as scVI / DE input. lognorm is the post-norm scale.
# HVG / scale / PCA / neighbors run on the panel_for_embedding subset
# (drops the 100 tissue-restricted + 6 sex-chromosome genes flagged in
# data/processed/panel_audit.csv). Sex-chromosome exclusion is critical:
# 0041326 is male (KDM5D ~10%), 0041323 is female (KDM5D ~0%) —
# Y-chromosome genes would otherwise drive the dominant PC and embedding
# axis, conflating donor-sex with biology.

# --- Set panel masks + donor_sex inline from panel_audit.csv if missing.
# This makes the cell self-contained so a fresh-from-zarr run still gets
# the audit-aware feature set.
if 'panel_for_embedding' not in adata.var.columns:
    import pandas as _pd
    _audit = _pd.read_csv('data/processed/panel_audit.csv')
    _nonsex = set(_audit.loc[(_audit['exclude_recommended']=='yes') &
                              (_audit['category']!='sex_chromosome'), 'gene'])
    _sex_chrom = set(_audit.loc[_audit['category']=='sex_chromosome', 'gene'])
    adata.var['panel_for_scoring'] = ~adata.var_names.isin(_nonsex)
    adata.var['panel_for_embedding'] = (~adata.var_names.isin(_nonsex)
                                          & ~adata.var_names.isin(_sex_chrom))
    print(f"panel_for_scoring: {int(adata.var['panel_for_scoring'].sum()):,} kept")
    print(f"panel_for_embedding: {int(adata.var['panel_for_embedding'].sum()):,} kept")

if 'donor_sex' not in adata.obs.columns:
    _LOG_DETECT = float(np.log1p(1.0))
    _sex_genes = [g for g in ('KDM5D', 'DDX3Y') if g in adata.var_names]
    if _sex_genes:
        _layer = adata.layers['lognorm'] if 'lognorm' in adata.layers else adata.X
        _n_max = 0
        for _g in _sex_genes:
            _i = adata.var_names.get_loc(_g)
            _col = _layer[:, _i]
            if hasattr(_col, 'toarray'):
                _col = _col.toarray()
            _col = np.asarray(_col).ravel()
            _n_max = max(_n_max, int((_col >= _LOG_DETECT).sum()))
        _sex = 'male' if _n_max / max(adata.n_obs, 1) >= 0.02 else 'female'
    else:
        _sex = 'unknown'
    import pandas as _pd
    adata.obs['donor_sex'] = _pd.Categorical([_sex]*adata.n_obs,
                                               categories=['female','male','unknown'])
    print(f"donor_sex tag: {_sex}")

adata.layers['counts'] = adata.X.copy()

# HVG metadata only — required for the HVG_score column the QC report uses.
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3', layer='counts')
print(f"HVG: {adata.var['highly_variable'].sum():,} of {adata.n_vars:,} flagged (metadata only).")

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.layers['lognorm'] = adata.X.copy()
print("normalize_total + log1p applied; layers['lognorm'] saved.")

# Panel mask — set by scripts/finalize_panel.py before this notebook runs.
# Fall back gracefully if the user has skipped that step (e.g. dev runs).
if 'panel_for_embedding' in adata.var.columns:
    embed_mask = adata.var['panel_for_embedding'].astype(bool).to_numpy()
    print(f"panel_for_embedding: {int(embed_mask.sum()):,} / {adata.n_vars:,} genes used for PCA/scVI.")
else:
    print("WARNING: var['panel_for_embedding'] not set — running PCA on full panel. "
          "Run scripts/finalize_panel.py against this h5ad before re-running for clean reproduction.")
    embed_mask = None

# Build a view restricted to the embedding-mask genes for scale + PCA.
# We compute on the subset and copy results back so adata keeps full panel.
if embed_mask is not None:
    adata_emb = adata[:, embed_mask].copy()
else:
    adata_emb = adata

# GPU scale (zero_center=False so sparse stays sparse)
print("rsc.pp.scale (GPU)...")
rsc.get.anndata_to_GPU(adata_emb)
rsc.pp.scale(adata_emb, zero_center=False, max_value=10)

# GPU PCA
print("rsc.pp.pca (GPU)...")
rsc.pp.pca(adata_emb, n_comps=30, random_state=0)

# Copy results back to the full-panel adata
rsc.get.anndata_to_CPU(adata_emb, convert_all=True)
adata.obsm['X_pca'] = adata_emb.obsm['X_pca']
adata.uns['pca'] = adata_emb.uns.get('pca', {})
del adata_emb

# Sanity: PC1 should not be dominated by Y-chromosome variance.
# (If panel_for_embedding excluded Y genes correctly, KDM5D/DDX3Y won't
# even be in adata_emb. This check is for the unguarded fallback case.)
from scipy.stats import spearmanr as _spearmanr
for sex_g in ('KDM5D', 'DDX3Y'):
    if sex_g in adata.var_names:
        i = adata.var_names.get_loc(sex_g)
        col = adata.layers['lognorm'][:, i]
        if hasattr(col, 'toarray'):
            col = col.toarray().ravel()
        else:
            col = np.asarray(col).ravel()
        r1, _ = _spearmanr(adata.obsm['X_pca'][:, 0], col)
        if abs(r1) > 0.1:
            print(f"  WARNING: Spearman(PC1, {sex_g}) = {r1:.3f} > 0.1 — sex-chromosome leak into PCA.")
        else:
            print(f"  Spearman(PC1, {sex_g}) = {r1:.3f}  OK")

# GPU neighbors on PCA
print("rsc.pp.neighbors (GPU)...")
rsc.pp.neighbors(adata, n_neighbors=15, use_rep='X_pca', metric='cosine', random_state=0)

# CPU UMAP (rapids UMAP outliers — see CLAUDE.md decision #4)
print("sc.tl.umap (CPU)...")
sc.tl.umap(adata, min_dist=0.5, random_state=0)

# CPU leiden, igraph flavor
print("sc.tl.leiden (CPU, igraph) ...")
sc.tl.leiden(adata, resolution=1.5, flavor='igraph', n_iterations=-1,
             directed=False, random_state=0, key_added='leiden_1.5')
print(f"leiden_1.5: {adata.obs['leiden_1.5'].nunique()} clusters.")
"""
    _set_src(cell, new)
    print("    [01 cell 18] patched")
    return True


def patch_nb02_cell10(cell):
    src = _src(cell)
    if "Rule-based Gamma" in src or "Epsilon identified post-hoc" in src:
        print("    [02 cell 10] already patched, skipping")
        return False
    new = """# Panel-aware pancreas markers for the T1D-focused hAtlas_v1.1+100 panel.
# Absent from panel: INS, GCG, SST, PPY, IRX2, PRSS1/2, CTRB1/2, CELA1/2A/3A,
# REG1A/3A. Beta uses surrogates; Alpha/Delta use multi-marker panels.
#
# RARE TYPES (Epsilon, Gamma) are NOT scored via score_genes argmax —
# they're identified post-hoc in scripts/refine_rare_celltypes.py via
# multi-criteria rules (e.g. pan-endocrine + AQP3+ETV1 + NOT Beta/Delta
# = Gamma). Score-genes argmax with weak/single-gene panels mislabeled
# 150K cells in 0041326 as Epsilon under the previous design (the bug
# that triggered this rebuild).
#
# Acinar reduces to AMY1A + AMY2A + AMY2B; the anti-acinar override
# in cell 15 routes mis-Acinar calls based on stronger alt-lineage
# markers.
pancreas_markers = {
    'Beta':          ['IAPP', 'MAFA', 'NKX6-1', 'PDX1', 'SLC30A8', 'ABCC8',
                      'KCNJ11', 'GLP1R', 'PCSK1'],
    'Alpha':         ['ARX', 'MAFB', 'GC', 'TM4SF4', 'TTR'],
    'Delta':         ['HHEX', 'LEPR', 'RBP4', 'BHLHE41', 'PCSK1'],
    # Epsilon identified post-hoc via rule-based multi-criteria — DO NOT add
    # back to score_genes argmax. See scripts/refine_rare_celltypes.py.
    # Gamma identified post-hoc via rule-based multi-criteria — DO NOT add
    # back to score_genes argmax (PPY is not on the panel; AQP3 + ETV1 +
    # pan-endocrine is the multi-marker definition).
    'Endocrine':     ['CHGA', 'INSM1', 'ISL1', 'NEUROD1', 'FEV', 'PAX6', 'SCG5',
                      'SCG2', 'SCGN'],
    'Acinar':        ['AMY1A', 'AMY2A', 'AMY2B'],
    'Ductal':        ['KRT19', 'SOX9', 'HNF1B', 'MUC1', 'CFTR', 'KRT7'],
    'Stellate':      ['ACTA2', 'PDGFRB', 'RGS5', 'PDGFRA', 'DCN', 'COL1A1'],
    'Endothelial':   ['PECAM1', 'CDH5', 'CLDN5', 'PLVAP', 'KDR', 'ENG', 'VWF'],
    'T_cells':       ['CD3D', 'CD3E', 'CD3G', 'CD8A', 'CD8B', 'CD4', 'PTPRC',
                      'FOXP3'],
    'B_cells':       ['MS4A1', 'CD79A', 'CD79B', 'CD19'],
    'Myeloid':       ['CD68', 'CD163', 'CSF1R', 'MARCO', 'CD14'],
    'Schwann':       ['MPZ', 'SOX10', 'NES', 'S100B', 'PMP22', 'PLP1'],
    'Proliferating': ['MKI67', 'TOP2A', 'PCNA'],
}

# Low-confidence types (panel coverage thin). Used to downgrade ambiguous
# assignments via the post-hoc panel-coverage check in
# scripts/refine_rare_celltypes.py. Epsilon removed (no longer score_genes-target).
LOW_CONFIDENCE_TYPES = {'Alpha', 'Delta', 'Acinar'}

# Apply the panel-audit exclusion mask. Set by scripts/finalize_panel.py.
# Drops 100 non-sex tissue-restricted genes + 6 sex-chromosome genes from
# the per-type panels.
if 'panel_for_scoring' in adata.var.columns:
    panel_keep = set(adata.var_names[adata.var['panel_for_scoring'].astype(bool)])
    print(f"panel_for_scoring: keeping {len(panel_keep):,} of {adata.n_vars:,} genes")
else:
    print("WARNING: var['panel_for_scoring'] not set — running with full panel. "
          "Run scripts/finalize_panel.py before re-running.")
    panel_keep = set(adata.var_names)

available_markers = {}
for ct, genes in pancreas_markers.items():
    kept = [g for g in genes if g in adata.var_names and g in panel_keep]
    if not kept:
        print(f"  {ct}: NO markers on panel; dropped from scoring")
        continue
    if len(kept) < 2:
        # CLAUDE.md decision #12: panels with <2 genes too noisy for argmax.
        # The single-gene panel was the Epsilon bug.
        print(f"  {ct}: only {len(kept)} marker(s) ({kept}); dropped from scoring "
              f"to prevent single-gene mislabel bug")
        continue
    available_markers[ct] = kept
    if len(kept) < 3:
        print(f"  WARNING: {ct} has only {len(kept)} markers ({kept}); flagging "
              f"as LOW_CONFIDENCE_TYPES even if not currently in the set")
        LOW_CONFIDENCE_TYPES.add(ct)

print(f"\\nCell types scored via score_genes argmax: {list(available_markers.keys())}")
print(f"LOW_CONFIDENCE_TYPES (downgrade-candidates after argmax): {sorted(LOW_CONFIDENCE_TYPES)}")
print("Epsilon + Gamma identified post-hoc via rule-based multi-criteria — "
      "see scripts/refine_rare_celltypes.py.")
"""
    _set_src(cell, new)
    print("    [02 cell 10] patched")
    return True


def patch_nb02_cell17(cell):
    src = _src(cell)
    if "donor_sex" in src and "panel_for_embedding" in src:
        print("    [02 cell 17] already patched, skipping")
        return False
    # Replace just the scVI setup_anndata call + the model construction.
    # We keep the rest of the cell intact (UMAP + leiden + QC overlay + PCA fallback).
    old_setup = """scvi.model.SCVI.setup_anndata(
    adata,
    layer='counts',
    batch_key=None,
    continuous_covariate_keys=['cell_area'],
)
vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood='nb', dropout_rate=0.1)"""
    new_setup = """# Train scVI on the panel_for_embedding subset only — drops Y-chromosome
# genes (KDM5D, DDX3Y et al.) so the donor-sex difference between samples
# (0041326 male, 0041323 female) doesn't dominate the latent. Also drops
# the 100 tissue-restricted genes flagged in panel_audit.csv that have
# no biological role in pancreas / LN / vessel context.
# donor_sex obs-column is set by scripts/finalize_panel.py from
# KDM5D + DDX3Y detection rates.
if 'panel_for_embedding' in adata.var.columns:
    embed_mask = adata.var['panel_for_embedding'].astype(bool).to_numpy()
    print(f"scVI training on {int(embed_mask.sum()):,} of {adata.n_vars:,} panel_for_embedding genes")
    adata_scvi = adata[:, embed_mask].copy()
else:
    print("WARNING: var['panel_for_embedding'] not set; training scVI on full panel "
          "(sex-chromosome leak likely). Run scripts/finalize_panel.py first.")
    adata_scvi = adata

cat_covariates = []
if 'donor_sex' in adata_scvi.obs.columns:
    cat_covariates = ['donor_sex']
    print(f"scVI categorical_covariate_keys: {cat_covariates}")

scvi.model.SCVI.setup_anndata(
    adata_scvi,
    layer='counts',
    batch_key=None,
    continuous_covariate_keys=['cell_area'],
    categorical_covariate_keys=cat_covariates if cat_covariates else None,
)
vae = scvi.model.SCVI(adata_scvi, n_layers=2, n_latent=30, gene_likelihood='nb', dropout_rate=0.1)"""
    if old_setup not in src:
        print("    [02 cell 17] expected scVI setup block not found — manual review needed")
        return False
    new = src.replace(old_setup, new_setup)
    # Also replace the get_latent_representation line so it reads from adata_scvi
    new = new.replace(
        "adata.obsm['X_scvi'] = vae.get_latent_representation()",
        "adata.obsm['X_scvi'] = vae.get_latent_representation()  # row-aligned with adata\ndel adata_scvi"
    )
    _set_src(cell, new)
    print("    [02 cell 17] patched")
    return True


def patch_notebook(path, cell_index, patch_fn, label):
    print(f"\n=== {path.name} (cell {cell_index} — {label}) ===")
    backup = path.with_suffix(path.suffix + BACKUP_SUFFIX)
    if not backup.exists():
        shutil.copy2(path, backup)
        print(f"    backup → {backup.name}")
    with open(path) as f:
        nb = json.load(f)
    if cell_index >= len(nb["cells"]):
        print(f"    cell {cell_index} out of range (notebook has {len(nb['cells'])}); skip")
        return
    if patch_fn(nb["cells"][cell_index]):
        with open(path, "w") as f:
            json.dump(nb, f, indent=1, ensure_ascii=False)
            f.write("\n")
        print(f"    {path.name} written")


def main():
    patch_notebook(NB01, 18, patch_nb01_cell18, "preprocessing scale/PCA/neighbors")
    patch_notebook(NB02, 10, patch_nb02_cell10, "marker dict + panel mask")
    patch_notebook(NB02, 17, patch_nb02_cell17, "scVI w/ donor_sex covariate + panel mask")
    print("\nNotebook patches applied. Verify by re-running notebook cells.")


if __name__ == "__main__":
    main()
