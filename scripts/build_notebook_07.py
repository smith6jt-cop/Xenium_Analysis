"""Build notebooks/07_immune_islet_proximity.ipynb.

The notebook subsets immune cells across both samples, retrains scVI on the
immune-only subset (with sample as batch covariate), assigns immune subtypes
via canary panels, identifies islets via DBSCAN on endocrine spatial coords,
and analyzes immune subtype distribution by distance to islet.
"""
import json
from pathlib import Path

NB = Path("/blue/maigan/smith6jt/Xenium_Analysis/notebooks/07_immune_islet_proximity.ipynb")


def md(text):
    return {"cell_type": "markdown", "metadata": {}, "source": text.splitlines(keepends=True)}


def code(text):
    return {
        "cell_type": "code",
        "metadata": {},
        "execution_count": None,
        "outputs": [],
        "source": text.splitlines(keepends=True),
    }


cells = []

cells.append(md("""# 07 — Immune phenotyping + islet proximity (T1D-focused)

Subsets all immune cells (T_cells, B_cells, Myeloid, NK candidates) across
both samples, retrains scVI on the immune-only subset (with sample as a
batch covariate so cross-sample subtypes align), assigns immune subtypes
via curated marker panels, identifies islets via DBSCAN on endocrine
spatial coordinates, and analyzes immune subtype distribution by distance
to islet.

**Inputs**: `_phenotyped.h5ad` for both samples (must already have
`celltype_lineage` from the lineage_phenotyping step in `02_phenotyping.ipynb`
or `scripts/lineage_phenotype.py`).

**Outputs**:
- `data/processed/{sample}/{sample}_immune_phenotyped.h5ad` — per-sample
  immune subset with subtype labels and distance-to-islet annotations
- `data/processed/immune_proximity_summary.csv` — cross-sample contingency
- `figures/{sample}/07_immune/*.png` — per-sample figures
- `figures/combined_07_immune/*.png` — combined-sample figures
"""))

cells.append(code("""# Imports + B200 TF32 setup (mirrors stage 02 cell 1)
import os
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

import rapids_singlecell as rsc
import cupy as cp
from cuml.cluster import DBSCAN as cuDBSCAN
from cuml.neighbors import NearestNeighbors as cuNN

if torch.cuda.is_available():
    torch.set_float32_matmul_precision('high')
    torch.backends.cuda.matmul.allow_tf32 = True
    torch.backends.cudnn.allow_tf32 = True
    print(f"GPU: {torch.cuda.get_device_name(0)} (TF32 enabled)")

sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white', frameon=False)

print(f"scanpy {sc.__version__}  scvi {scvi.__version__}  rsc {rsc.__version__}")
"""))

cells.append(md("""## 1. Load both samples + concat with sample batch key"""))

cells.append(code("""SAMPLES = ("0041323", "0041326")
DATA_ROOT = Path(os.environ.get("XENIUM_DATA_ROOT", "../data/processed"))
FIG_ROOT = Path(os.environ.get("XENIUM_FIG_ROOT", "../figures"))

adatas = []
for s in SAMPLES:
    p = DATA_ROOT / s / f"{s}_phenotyped.h5ad"
    print(f"loading {p} ...")
    a = sc.read_h5ad(p)
    a.obs["sample"] = s
    a.obs_names = [f"{s}__{n}" for n in a.obs_names]
    adatas.append(a)
    print(f"  {a.shape}")

# Concatenate. Use only shared genes (var intersection).
adata_full = ad.concat(adatas, join="inner", merge="same", index_unique=None,
                       label="sample", keys=SAMPLES)
adata_full.obs["sample"] = adata_full.obs["sample"].astype("category")
print(f"\\nconcatenated: {adata_full.shape}")
print(f"sample distribution:")
print(adata_full.obs["sample"].value_counts().to_string())

del adatas
"""))

cells.append(md("""## 2. Subset to immune cells

Includes:
- `celltype` ∈ {T_cells, B_cells, Myeloid}
- `celltype_lineage` ∈ {T_cells, B_cells, Myeloid} (recovers
  lineage_phenotyped reclassifications, e.g. the 194 Acinar→Myeloid
  macrophages-with-engulfed-acinar from 0041323)
- Cells in any Indeterminate bucket with raw counts ≥ 3 of CD3E/CD79A/CD68
  (canonical immune canary markers)
"""))

cells.append(code("""IMMUNE_LABELS = {"T_cells", "B_cells", "Myeloid"}

# Build the subset mask
ct = adata_full.obs["celltype"].astype(str)
ct_lin = adata_full.obs["celltype_lineage"].astype(str) if "celltype_lineage" in adata_full.obs.columns else ct
indet_mask = ct.str.startswith("Indeterminate")

immune_by_label = ct.isin(IMMUNE_LABELS) | ct_lin.isin(IMMUNE_LABELS)

# Recover indeterminate cells with strong immune marker signal
canary_immune = {"T_cells": ["CD3E", "CD3G", "PTPRC", "CD2"],
                 "B_cells": ["MS4A1", "CD79A", "CD79B"],
                 "Myeloid": ["CD68", "CD163", "CSF1R"]}
present_canary = [g for genes in canary_immune.values() for g in genes if g in adata_full.var_names]
canary_idx = [adata_full.var_names.get_loc(g) for g in present_canary]
counts_canary = adata_full.layers["counts"][:, canary_idx]
if hasattr(counts_canary, "toarray"):
    counts_canary = counts_canary.toarray()
counts_canary = np.asarray(counts_canary)
immune_marker_hit = (counts_canary >= 3).any(axis=1)
recovered = indet_mask.values & immune_marker_hit

immune_mask = (immune_by_label.values | recovered)
adata = adata_full[immune_mask].copy()
print(f"immune subset: {adata.shape}")
print(f"  by label:        {int(immune_by_label.sum()):,}")
print(f"  recovered indet: {int(recovered.sum()):,}")
print(f"per-sample:")
print(adata.obs["sample"].value_counts().to_string())
print(f"per-celltype (original label):")
print(adata.obs["celltype"].value_counts().to_string())

del adata_full
"""))

cells.append(md("""## 3. HVG + GPU scale + GPU PCA on immune subset

The immune subset has fewer cells (≈170k) so HVG selection is more
informative when restricted to this set. Restoring counts → normalize →
log → seurat_v3 HVG → keep all panel genes (no subsetting), scale on GPU,
PCA on GPU.
"""))

cells.append(code("""# Restore raw counts as .X for HVG/normalize
adata.X = adata.layers["counts"].copy()
print("normalize_total + log1p ...")
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.layers["lognorm"] = adata.X.copy()

# HVG (seurat_v3 needs raw counts via layer arg)
print("HVG seurat_v3 (immune-set) ...")
sc.pp.highly_variable_genes(adata, layer="counts", n_top_genes=2000,
                             flavor="seurat_v3", subset=False)
print(f"  HVG count: {int(adata.var['highly_variable'].sum())}")

# GPU scale + PCA
print("anndata_to_GPU + rsc.pp.scale + rsc.pp.pca ...")
rsc.get.anndata_to_GPU(adata)
rsc.pp.scale(adata, zero_center=False, max_value=10)
rsc.pp.pca(adata, n_comps=50, random_state=0)
rsc.get.anndata_to_CPU(adata, convert_all=True)
adata.X = adata.layers["lognorm"].copy()
print(f"  X_pca shape: {adata.obsm['X_pca'].shape}")
"""))

cells.append(md("""## 4. scVI retrain on immune subset

Uses `sample` as a batch covariate so cross-sample subtypes align (immune
subtypes share markers across samples; the latent space should not segregate
by sample). `n_latent=20` is enough for an immune-only set this size.
"""))

cells.append(code("""scvi.model.SCVI.setup_anndata(
    adata,
    layer='counts',
    batch_key='sample',
    continuous_covariate_keys=['cell_area'],
)

vae = scvi.model.SCVI(adata, n_layers=2, n_latent=20,
                      gene_likelihood='nb', dropout_rate=0.1)
print("training scVI on immune subset ...")
vae.train(
    max_epochs=200,
    early_stopping=True,
    early_stopping_patience=20,
    accelerator=('gpu' if torch.cuda.is_available() else 'cpu'),
    devices=1,
    batch_size=4096,
    datasplitter_kwargs={"num_workers": 8, "pin_memory": True,
                         "persistent_workers": True},
)
adata.obsm['X_scvi_immune'] = vae.get_latent_representation()
print(f"X_scvi_immune shape: {adata.obsm['X_scvi_immune'].shape}")
"""))

cells.append(md("""## 5. Neighbors + UMAP + leiden on immune scVI latent"""))

cells.append(code("""print("rsc.pp.neighbors on X_scvi_immune ...")
rsc.pp.neighbors(adata, use_rep='X_scvi_immune', n_neighbors=30,
                 metric='cosine', random_state=0, key_added='immune')

print("sc.tl.umap (CPU) ...")
sc.tl.umap(adata, neighbors_key='immune', min_dist=0.3, spread=1.0,
           random_state=0)

print("sc.tl.leiden res=0.8 (CPU igraph) ...")
sc.tl.leiden(adata, resolution=0.8, neighbors_key='immune',
             key_added='leiden_immune', flavor='igraph',
             n_iterations=-1, directed=False, random_state=0)
print(f"  {adata.obs['leiden_immune'].nunique()} immune clusters")
"""))

cells.append(md("""## 6. Score immune subtype markers + assign subtypes

Subtype panels chosen from probed availability:
- **T_helper**: CD4, ICOS, IL7R, CCR7
- **T_cytotoxic**: CD8A, CD8B, GZMA, GZMB, GZMK, PRF1, NKG7
- **T_reg**: FOXP3, IL2RA, CTLA4, IKZF2
- **T_exhausted**: PDCD1, LAG3, HAVCR2, TIGIT, TOX
- **NK**: NCAM1, NKG7, KLRD1, KLRF1, FCGR3A
- **Macro_M1**: TNF, IL1B, IL6, NOS2, CD86
- **Macro_M2**: CD163, MRC1, MERTK
- **Macro_resident**: MARCO, TIMD4, LYVE1
- **Monocyte**: CD14, FCN1
- **DC**: ITGAX, CD1C, FCER1A, CLEC9A, FLT3, LAMP3
- **B_pan**: MS4A1, CD79A, CD79B, CD19
- **B_plasma**: SDC1, MZB1, PRDM1
"""))

cells.append(code("""SUBTYPE_PANELS = {
    "T_helper":      ["CD4", "ICOS", "IL7R", "CCR7"],
    "T_cytotoxic":   ["CD8A", "CD8B", "GZMA", "GZMB", "GZMK", "PRF1", "NKG7"],
    "T_reg":         ["FOXP3", "IL2RA", "CTLA4", "IKZF2"],
    "T_exhausted":   ["PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX"],
    "NK":            ["NCAM1", "NKG7", "KLRD1", "KLRF1", "FCGR3A"],
    "Macro_M1":      ["TNF", "IL1B", "IL6", "NOS2", "CD86"],
    "Macro_M2":      ["CD163", "MRC1", "MERTK"],
    "Macro_resident":["MARCO", "TIMD4", "LYVE1"],
    "Monocyte":      ["CD14", "FCN1"],
    "DC":            ["ITGAX", "CD1C", "FCER1A", "CLEC9A", "FLT3", "LAMP3"],
    "B_pan":         ["MS4A1", "CD79A", "CD79B", "CD19"],
    "B_plasma":      ["SDC1", "MZB1", "PRDM1"],
}

# Score each subtype using sc.tl.score_genes on the lognorm layer
for st, genes in SUBTYPE_PANELS.items():
    present = [g for g in genes if g in adata.var_names]
    if not present:
        print(f"  {st}: no markers on panel — skipped")
        continue
    sc.tl.score_genes(adata, present, score_name=f"{st}_score", use_raw=False)
    print(f"  scored {st}: {present}")

# Argmax across subtype scores
score_cols = [f"{st}_score" for st in SUBTYPE_PANELS if f"{st}_score" in adata.obs.columns]
adata.obs["immune_subtype"] = (
    adata.obs[score_cols].idxmax(axis=1).str.replace("_score", "")
)

# Confidence margin: top - second
top2 = np.sort(adata.obs[score_cols].values, axis=1)[:, -2:]
adata.obs["immune_subtype_margin"] = top2[:, 1] - top2[:, 0]

# Per-cluster consensus (majority vote within leiden_immune)
consensus = (adata.obs.groupby('leiden_immune', observed=True)['immune_subtype']
             .agg(lambda s: s.value_counts().idxmax()).to_dict())
adata.obs['immune_cluster_consensus'] = (
    adata.obs['leiden_immune'].map(consensus).astype('category')
)

print("\\nper-cell immune_subtype distribution:")
print(adata.obs['immune_subtype'].value_counts().to_string())
print("\\nper-cluster consensus:")
print(adata.obs['immune_cluster_consensus'].value_counts().to_string())
"""))

cells.append(md("""## 7. Identify islets via DBSCAN on endocrine spatial coordinates

Per sample: take the cells (from the FULL phenotyped h5ad, not the immune
subset) where `celltype_lineage` ∈ {Beta, Alpha, Delta, Epsilon, Endocrine,
Endocrine_pan}. Cluster their `obsm['spatial']` coords with DBSCAN
(`eps=30 μm`, `min_samples=20`) → each connected cluster = an islet. Compute
per-islet centroid, area (convex hull), and a unique `islet_id` (per sample).
"""))

cells.append(code("""ISLET_LABELS = {"Beta", "Alpha", "Delta", "Epsilon", "Endocrine", "Endocrine_pan"}
DBSCAN_EPS_UM = 30.0
DBSCAN_MIN_SAMPLES = 20

islets_by_sample = {}  # sample_id -> DataFrame[islet_id, centroid_x, centroid_y, n_cells]

for s in SAMPLES:
    p = DATA_ROOT / s / f"{s}_phenotyped.h5ad"
    print(f"\\n=== {s} ===")
    a_full = sc.read_h5ad(p)
    ct_lin = a_full.obs.get("celltype_lineage", a_full.obs["celltype"]).astype(str)
    islet_mask = ct_lin.isin(ISLET_LABELS).values
    coords = np.asarray(a_full.obsm["spatial"], dtype=np.float32)[islet_mask]
    print(f"  endocrine cells for DBSCAN: {coords.shape[0]:,}")

    if coords.shape[0] == 0:
        print(f"  no endocrine cells — skipping islets")
        continue

    db = cuDBSCAN(eps=DBSCAN_EPS_UM, min_samples=DBSCAN_MIN_SAMPLES,
                  output_type='numpy')
    labels = db.fit_predict(coords)
    n_islets = int((labels >= 0).max() + 1) if (labels >= 0).any() else 0
    n_noise = int((labels == -1).sum())
    print(f"  islets identified: {n_islets:,}  (noise: {n_noise:,} endocrine cells not in any islet)")

    rows = []
    for islet_id in range(n_islets):
        mem = labels == islet_id
        c = coords[mem]
        rows.append({
            "sample": s,
            "islet_id": f"{s}_islet_{islet_id:04d}",
            "centroid_x": float(c[:, 0].mean()),
            "centroid_y": float(c[:, 1].mean()),
            "n_cells": int(mem.sum()),
            "x_min": float(c[:, 0].min()),
            "x_max": float(c[:, 0].max()),
            "y_min": float(c[:, 1].min()),
            "y_max": float(c[:, 1].max()),
            "approx_radius_um": float(0.5 * (c[:, 0].max() - c[:, 0].min()
                                              + c[:, 1].max() - c[:, 1].min()) / 2),
        })
    islets_df = pd.DataFrame(rows)
    islets_by_sample[s] = islets_df
    print(f"  islet size dist: median={islets_df['n_cells'].median():.0f}  "
          f"p90={islets_df['n_cells'].quantile(0.90):.0f}  "
          f"max={islets_df['n_cells'].max()}")
    del a_full

# Save combined islet table
all_islets = pd.concat(islets_by_sample.values(), ignore_index=True)
all_islets.to_csv(DATA_ROOT / "islets_dbscan.csv", index=False)
print(f"\\ntotal islets across both samples: {len(all_islets):,}")
"""))

cells.append(md("""## 8. Distance to nearest islet (per immune cell, per sample)

Uses cuML NearestNeighbors on the islet centroids (per sample) — for each
immune cell, finds the nearest islet centroid and the distance in μm.
Subsequently bins the distance into discrete categories.
"""))

cells.append(code("""DIST_BINS_UM = [0, 50, 200, np.inf]
DIST_LABELS = ["intra_or_peri_islet", "proximal", "distal"]

immune_dist = []
for s in SAMPLES:
    if s not in islets_by_sample or len(islets_by_sample[s]) == 0:
        continue
    sub = adata[adata.obs["sample"] == s].copy()
    if sub.n_obs == 0:
        continue
    immune_xy = np.asarray(sub.obsm["spatial"], dtype=np.float32)
    islet_xy = islets_by_sample[s][["centroid_x", "centroid_y"]].values.astype(np.float32)

    nn = cuNN(n_neighbors=1, algorithm="brute", output_type="numpy")
    nn.fit(islet_xy)
    dist, idx = nn.kneighbors(immune_xy, return_distance=True)
    nearest_islet = islets_by_sample[s].iloc[idx[:, 0].astype(int)]["islet_id"].values
    sub.obs["dist_to_islet_um"] = dist[:, 0]
    sub.obs["nearest_islet_id"] = nearest_islet
    sub.obs["distance_bin"] = pd.cut(
        dist[:, 0], bins=DIST_BINS_UM, labels=DIST_LABELS, include_lowest=True,
    )
    immune_dist.append(sub.obs[["sample", "immune_subtype",
                                 "immune_cluster_consensus",
                                 "dist_to_islet_um", "nearest_islet_id",
                                 "distance_bin"]].copy())
    # Mirror back into combined adata
    adata.obs.loc[sub.obs.index, "dist_to_islet_um"] = dist[:, 0]
    adata.obs.loc[sub.obs.index, "nearest_islet_id"] = nearest_islet
    adata.obs.loc[sub.obs.index, "distance_bin"] = sub.obs["distance_bin"].values
    print(f"{s}: {sub.n_obs:,} immune cells, distance dist:")
    print(f"  {sub.obs['distance_bin'].value_counts().to_string()}")

adata.obs["distance_bin"] = pd.Categorical(adata.obs["distance_bin"],
                                            categories=DIST_LABELS)
"""))

cells.append(md("""## 9. Contingency analysis + statistical tests"""))

cells.append(code("""from scipy.stats import chi2_contingency, fisher_exact

# Build contingency table: subtype × distance_bin × sample
ct = pd.crosstab(
    [adata.obs["sample"], adata.obs["immune_subtype"]],
    adata.obs["distance_bin"],
)
print("contingency: (sample, subtype) × distance_bin")
print(ct.to_string())

# Save full contingency to CSV — write to _legacy variant so we don't
# overwrite the per-phenotype-aware headline owned by
# scripts/insulitis_analysis.py (post-2026-04-28 rewrite).
ct.to_csv(DATA_ROOT / "immune_proximity_summary_legacy_notebook07.csv")

# Per-sample chi-square: is immune subtype independent of distance bin?
print("\\nper-sample chi-square (subtype × distance_bin):")
for s in SAMPLES:
    sub = adata.obs[adata.obs["sample"] == s]
    table = pd.crosstab(sub["immune_subtype"], sub["distance_bin"])
    chi2, pval, dof, _ = chi2_contingency(table.values)
    print(f"  {s}: chi2={chi2:.1f}  dof={dof}  p={pval:.2e}  shape={table.shape}")

# T1D-relevant per-subtype enrichment: which subtypes are over-represented
# intra/peri-islet vs distal?
print("\\nper-subtype intra/peri-islet enrichment (intra_or_peri vs distal):")
for s in SAMPLES:
    sub = adata.obs[adata.obs["sample"] == s]
    print(f"\\n  {s}:")
    print(f"  {'subtype':<14} {'intra/peri':>11} {'distal':>10} {'odds_ratio':>12} {'p':>10}")
    for st in sub["immune_subtype"].astype("category").cat.categories:
        a = ((sub["immune_subtype"] == st) & (sub["distance_bin"] == "intra_or_peri_islet")).sum()
        b = ((sub["immune_subtype"] == st) & (sub["distance_bin"] == "distal")).sum()
        c = ((sub["immune_subtype"] != st) & (sub["distance_bin"] == "intra_or_peri_islet")).sum()
        d = ((sub["immune_subtype"] != st) & (sub["distance_bin"] == "distal")).sum()
        if a + c == 0 or b + d == 0:
            continue
        odds, pval = fisher_exact([[a, b], [c, d]])
        print(f"  {st:<14} {int(a):>11,} {int(b):>10,} {odds:>12.2f} {pval:>10.2e}")
"""))

cells.append(md("""## 10. Per-islet infiltration score

For each islet: count immune cells of each subtype within 50 μm; normalize
by islet size (n_endocrine_cells, a proxy for islet area). Useful T1D
metrics: CD8 cytotoxic / islet, Treg / islet, M2-vs-M1 macrophage ratio.
"""))

cells.append(code("""infiltration_rows = []
for s in SAMPLES:
    if s not in islets_by_sample:
        continue
    immune_sub = adata[
        (adata.obs["sample"] == s) &
        (adata.obs["distance_bin"].isin(["intra_or_peri_islet", "proximal"]))
    ].copy()
    print(f"\\n{s}: {immune_sub.n_obs:,} immune cells within proximal+ of islets")

    # Group by nearest_islet_id × subtype
    counts = (
        immune_sub.obs.groupby(["nearest_islet_id", "immune_subtype"], observed=True)
        .size().unstack(fill_value=0)
    )
    islets_meta = islets_by_sample[s].set_index("islet_id")
    counts = counts.reindex(islets_meta.index).fillna(0).astype(int)
    counts["n_endocrine"] = islets_meta["n_cells"]
    counts["sample"] = s
    counts["islet_id"] = counts.index
    infiltration_rows.append(counts.reset_index(drop=True))

infiltration = pd.concat(infiltration_rows, ignore_index=True) if infiltration_rows else pd.DataFrame()
if not infiltration.empty:
    # Normalize per islet
    for st in [c for c in infiltration.columns
               if c not in {"n_endocrine", "sample", "islet_id"}]:
        infiltration[f"{st}_per100endo"] = (
            100 * infiltration[st] / np.maximum(infiltration["n_endocrine"], 1)
        )

    # Per-sample summary stats
    print("\\nper-sample infiltration summary (median per islet):")
    summary_cols = [c for c in infiltration.columns if c.endswith("_per100endo")]
    print(
        infiltration.groupby("sample", observed=True)[summary_cols]
        .median()
        .round(3)
        .to_string()
    )

    # Write to _legacy variant so we don't overwrite the per-phenotype headline
    # owned by scripts/insulitis_analysis.py (post-2026-04-28 rewrite).
    infiltration.to_csv(DATA_ROOT / "islet_infiltration_per100endo_legacy_notebook07.csv", index=False)
    print(f"\\nsaved per-islet infiltration to islet_infiltration_per100endo_legacy_notebook07.csv")
"""))

cells.append(md("""## 11. Figures"""))

cells.append(code("""# 11a-11b. Combined-sample UMAP + all-subtype distance violin are now
# produced by scripts/figs_07_immune.py::fig13_combined_immune_umap and
# ::fig14_all_subtype_distance_violin, reading from the current
# _immune_phenotyped.h5ad files. Notebook 07's prior outputs to
# figures/combined_07_immune/ predated the lognorm sample-comparable
# threshold fix and were moved to figures/legacy/combined_07_immune_pre_lognorm_fix/.
print("[deprecated] combined-sample immune UMAP + per-subtype distance violin "
      "are now produced by scripts/figs_07_immune.py (fig13/fig14). "
      "Run that script instead of the notebook for these plots.")

# 11c. Subtype × distance_bin heatmap (per-sample)
for s in SAMPLES:
    sub_obs = adata.obs[adata.obs["sample"] == s]
    if len(sub_obs) == 0:
        continue
    table = pd.crosstab(sub_obs["immune_subtype"], sub_obs["distance_bin"],
                        normalize="index")
    fig, ax = plt.subplots(figsize=(8, max(6, 0.4*len(table))))
    sns.heatmap(table, annot=True, fmt=".2f", cmap="rocket_r", ax=ax,
                cbar_kws={"label": "row-normalized fraction"})
    ax.set_title(f"{s}: immune subtype × distance bin")
    sample_fig_dir = FIG_ROOT / s / "07_immune"
    sample_fig_dir.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(sample_fig_dir / "subtype_distance_heatmap.png",
                dpi=200, bbox_inches="tight")
    plt.show()
"""))

cells.append(code("""# 11d. Per-sample spatial overlay: immune subtypes + islet outlines
import matplotlib as mpl
import matplotlib.patches as patches

for s in SAMPLES:
    if s not in islets_by_sample:
        continue
    sub = adata[adata.obs["sample"] == s].copy()
    if sub.n_obs == 0:
        continue
    xy = np.asarray(sub.obsm["spatial"], dtype=np.float32)

    # Subsample for plot speed
    rng = np.random.default_rng(0)
    n_plot = min(80_000, sub.n_obs)
    idx = rng.choice(sub.n_obs, size=n_plot, replace=False)

    fig, ax = plt.subplots(figsize=(14, 14))
    cats = sub.obs["immune_subtype"].astype("category").iloc[idx]
    codes = cats.cat.codes.to_numpy()
    cmap = mpl.colormaps.get_cmap("tab20").resampled(max(20, len(cats.cat.categories)))
    sc_obj = ax.scatter(xy[idx, 0], xy[idx, 1], c=codes, cmap=cmap,
                        s=2, linewidths=0, alpha=0.7)

    # Overlay islet centroids as hollow circles
    islets_meta = islets_by_sample[s]
    for _, row in islets_meta.iterrows():
        # Approximate islet radius
        radius = max(15, row["approx_radius_um"])
        circ = patches.Circle((row["centroid_x"], row["centroid_y"]),
                              radius=radius, fill=False, edgecolor="red",
                              linewidth=0.6, alpha=0.7)
        ax.add_patch(circ)

    # Legend for subtypes (manually built)
    handles = [patches.Patch(color=cmap(i / max(1, len(cats.cat.categories)-1)),
                              label=str(cat))
               for i, cat in enumerate(cats.cat.categories)]
    ax.legend(handles=handles, loc="lower right", fontsize=8, framealpha=0.9,
              title="subtype")
    ax.set_aspect("equal")
    ax.set_title(f"{s}: immune cells (subsampled {n_plot:,}) + islet outlines (red)")
    sample_fig_dir = FIG_ROOT / s / "07_immune"
    sample_fig_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(sample_fig_dir / "spatial_immune_islets.png",
                dpi=200, bbox_inches="tight")
    plt.show()
"""))

cells.append(md("""## 12. Save annotated immune subset"""))

cells.append(code("""for s in SAMPLES:
    sub = adata[adata.obs["sample"] == s].copy()
    out = DATA_ROOT / s / f"{s}_immune_phenotyped.h5ad"
    sub.write_h5ad(out)
    print(f"  wrote {out}: {sub.shape}")
print("\\nDONE")
"""))


nb_obj = {
    "cells": cells,
    "metadata": {
        "kernelspec": {
            "display_name": "Python 3 (xenium_analysis)",
            "language": "python",
            "name": "python3",
        },
        "language_info": {"name": "python", "version": "3.12"},
    },
    "nbformat": 4,
    "nbformat_minor": 5,
}

NB.parent.mkdir(parents=True, exist_ok=True)
with NB.open("w") as f:
    json.dump(nb_obj, f, indent=1)

print(f"wrote {NB} with {len(cells)} cells")
