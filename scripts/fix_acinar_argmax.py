"""Hotfix for celltype assignment after the cell-10 Acinar-panel oversight.

Problem
-------
notebooks/02_phenotyping.ipynb cell 10 (post-patch) listed Acinar markers as
['AMY1A', 'AMY2A', 'AMY2B']. Only AMY1A + CUZD1 are actually on the panel
(AMY2A/B, CTRB1/2, CPA1/2, CELA*, PRSS1/2, REG*, GP2 — all absent).

The patched cell 10 enforces a ≥ 2-gene panel rule, so Acinar collapsed to
[AMY1A] and was DROPPED from score_genes argmax. With ~60% of pancreas cells
being acinar (~700K cells), they argmax-distributed across whichever score
was highest — bumping B_cells to 464K (vs ~5K real B cells in pancreas).

Fix
---
This script does NOT re-run scVI/UMAP/leiden (those are correct). It
re-scores Acinar with the right panel, re-argmaxes celltype, re-applies the
anti-acinar override, and updates lineage labels. ~5 min/sample.

Inputs:
    data/processed/{sample}/{sample}_phenotyped.h5ad     (post-stage-02)
    data/processed/{sample}/{sample}_spatial_analysis.h5ad  (post-stage-03)

Outputs (in-place):
    Updated celltype, predicted_celltype, celltype_lineage in BOTH h5ads.
    Adds Acinar_score column. Backs up to *.h5ad.bak.preAcinarFix once.
"""
import shutil
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path("/blue/maigan/smith6jt/Xenium_Analysis")
SAMPLES = ("0041323", "0041326")

# Acinar markers actually present on the panel (verified via panel_audit.csv).
ACINAR_PANEL = ["AMY1A", "CUZD1"]

# Anti-acinar override panels — same as in notebook 02 cell 15.
ANTI_ACINAR_MARKERS = {
    "Beta":        ["IAPP", "MAFA", "NKX6-1", "PDX1", "SLC30A8"],
    "Endocrine":   ["CHGA", "INSM1", "ISL1", "NEUROD1", "FEV"],
    "Ductal":      ["KRT19", "CFTR", "MUC1"],
    "Stellate":    ["ACTA2", "PDGFRB", "RGS5"],
    "Endothelial": ["PECAM1", "CDH5", "CLDN5", "PLVAP"],
    "T_cells":     ["CD3D", "CD3E", "CD8A", "PTPRC"],
    "B_cells":     ["MS4A1", "CD79A"],
    "Myeloid":     ["CD68", "CD163", "MARCO"],
    "Schwann":     ["MPZ", "SOX10"],
}

# Margin threshold for argmax confidence (matches user-locked decision)
MARGIN_THRESH = 0.2


def per_cell_marker_mean(adata, gene_list, layer="lognorm"):
    present = [g for g in gene_list if g in adata.var_names]
    if not present:
        return None
    X = adata[:, present].layers[layer]
    if hasattr(X, "toarray"):
        return np.asarray(X.mean(axis=1)).ravel()
    return np.asarray(X.mean(axis=1)).ravel()


def fix_one(sample):
    print(f"\n{'='*72}\n=== {sample}\n{'='*72}")
    p_pheno = ROOT / f"data/processed/{sample}/{sample}_phenotyped.h5ad"
    p_spatial = ROOT / f"data/processed/{sample}/{sample}_spatial_analysis.h5ad"

    for p in (p_pheno, p_spatial):
        bak = p.with_suffix(p.suffix + ".bak.preAcinarFix")
        if not bak.exists():
            print(f"  backing up {p.name} -> {bak.name}")
            shutil.copy2(p, bak)

    print(f"  loading {p_pheno.name}...")
    t = time.time()
    a = sc.read_h5ad(p_pheno)
    print(f"    loaded {a.shape} in {time.time()-t:.1f}s")

    # ----- Score Acinar with the right panel -----
    available = [g for g in ACINAR_PANEL if g in a.var_names]
    if not available:
        raise SystemExit(f"  no Acinar markers on panel for {sample}!")
    print(f"  scoring Acinar with {available} on lognorm layer...")
    if "lognorm" in a.layers:
        original_X = a.X
        a.X = a.layers["lognorm"]
        sc.tl.score_genes(a, gene_list=available, score_name="Acinar_score",
                            use_raw=False, random_state=0)
        a.X = original_X
    else:
        sc.tl.score_genes(a, gene_list=available, score_name="Acinar_score",
                            use_raw=False, random_state=0)

    # ----- Re-argmax across all _score columns -----
    score_cols = [c for c in a.obs.columns
                    if c.endswith("_score") and not c.startswith("Indeterminate")]
    print(f"  re-argmax across {len(score_cols)} score columns: {score_cols}")
    M = np.stack([a.obs[c].values.astype(np.float32) for c in score_cols], axis=1)
    sort_idx = np.argsort(-M, axis=1)
    top = M[np.arange(len(M)), sort_idx[:, 0]]
    second = M[np.arange(len(M)), sort_idx[:, 1]]
    margin = top - second
    type_names = [c.replace("_score", "") for c in score_cols]
    top_type = np.array([type_names[i] for i in sort_idx[:, 0]])

    # Apply margin gate (0.2 — user-locked)
    new_pred = top_type.copy()
    a.obs["predicted_celltype"] = pd.Categorical(top_type)
    a.obs["celltype_confidence"] = margin

    print(f"  pre-override predicted_celltype:")
    print(pd.Series(new_pred).value_counts().head(10).to_string())

    # ----- Anti-acinar override -----
    n_pre_acinar = int((new_pred == "Acinar").sum())
    print(f"\n  initial Acinar argmax calls: {n_pre_acinar:,}")

    alt_mean = {}
    alt_thr = {}
    for alt_ct, marker_list in ANTI_ACINAR_MARKERS.items():
        m = per_cell_marker_mean(a, marker_list, layer="lognorm")
        if m is None or m.size == 0:
            continue
        alt_mean[alt_ct] = m
        # Threshold: median over cells initially called alt_ct
        alt_cells = new_pred == alt_ct
        if alt_cells.sum() >= 50:
            alt_thr[alt_ct] = float(np.median(m[alt_cells]))
        else:
            alt_thr[alt_ct] = float(np.percentile(m, 90))
        print(f"    {alt_ct}: thr={alt_thr[alt_ct]:.3f}, "
              f"alt-cells={int(alt_cells.sum()):,}")

    # For cells argmax=Acinar, compute the alt-mean strengths and reassign
    # if any alt set's mean > its threshold.
    is_acinar = new_pred == "Acinar"
    n_reassigned = 0
    for alt_ct in ANTI_ACINAR_MARKERS:
        if alt_ct not in alt_mean:
            continue
        candidates = is_acinar & (alt_mean[alt_ct] > alt_thr[alt_ct])
        new_pred[candidates] = alt_ct
        is_acinar = new_pred == "Acinar"  # recompute after reassign
        n_reassigned += int(candidates.sum())
    print(f"  anti-acinar override: {n_reassigned:,} cells reassigned away from Acinar")

    # Apply margin gate AFTER override
    needs_indeterminate = (margin < MARGIN_THRESH) & (new_pred == top_type)
    n_indeterminate = int(needs_indeterminate.sum())
    final = new_pred.copy()
    final[needs_indeterminate] = ["Indeterminate_" + t for t in new_pred[needs_indeterminate]]
    print(f"  margin gate: {n_indeterminate:,} cells flagged as Indeterminate (margin < {MARGIN_THRESH})")

    a.obs["celltype"] = pd.Categorical(final)
    if "manual_celltype" in a.obs.columns:
        a.obs["manual_celltype"] = pd.Categorical(final)
    if "celltype_lineage" in a.obs.columns:
        # celltype_lineage starts as celltype unless lineage refinement reassigned
        # For simplicity, reset to celltype here. refine_rare_celltypes.py will
        # apply Gamma/Epsilon rules afterward.
        a.obs["celltype_lineage"] = pd.Categorical(final)

    # ----- Print final distribution -----
    print(f"\n  final celltype distribution:")
    print(pd.Series(final).value_counts().head(15).to_string())

    # ----- Save phenotyped h5ad -----
    print(f"\n  saving {p_pheno.name}...")
    t = time.time()
    a.write(p_pheno, compression="gzip")
    print(f"    saved in {time.time()-t:.1f}s")

    # ----- Mirror celltype + Acinar_score to spatial h5ad -----
    print(f"  syncing celltype to {p_spatial.name}...")
    t = time.time()
    sa = sc.read_h5ad(p_spatial)
    # Align by obs_names (both should have same row order from same source)
    if (sa.obs_names != a.obs_names).any():
        # Reindex if not aligned
        common_obs = a.obs.loc[sa.obs_names]
        sa.obs["celltype"] = common_obs["celltype"].values
        if "celltype_lineage" in common_obs.columns:
            sa.obs["celltype_lineage"] = common_obs["celltype_lineage"].values
        sa.obs["Acinar_score"] = common_obs["Acinar_score"].values
        sa.obs["predicted_celltype"] = common_obs["predicted_celltype"].values
    else:
        sa.obs["celltype"] = a.obs["celltype"].values
        if "celltype_lineage" in a.obs.columns:
            sa.obs["celltype_lineage"] = a.obs["celltype_lineage"].values
        sa.obs["Acinar_score"] = a.obs["Acinar_score"].values
        sa.obs["predicted_celltype"] = a.obs["predicted_celltype"].values
    sa.write(p_spatial, compression="gzip")
    print(f"    saved in {time.time()-t:.1f}s")


def main():
    for s in SAMPLES:
        fix_one(s)
    print(f"\n{'='*72}\n=== DONE\n{'='*72}")


if __name__ == "__main__":
    main()
