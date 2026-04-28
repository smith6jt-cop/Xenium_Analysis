"""For cells with n_genes > 1200, check for mutually-exclusive lineage co-expression.

If a single cell expresses high signal for >1 lineage that should be exclusive
(e.g., epithelial Ductal + lymphoid T-cell), the cell is a segmentation merge
(doublet) and the QC filter failed for that cell.

Lineage signals are computed as the per-cell mean of (lognorm) expression
across a small curated specific marker set per lineage. Specificity > sensitivity
here — we only want markers that should NEVER co-occur in the same lineage.
"""
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path("/blue/maigan/smith6jt/Xenium_Analysis")

# Tight, mutually-exclusive lineage panels. Drop pan-endocrine markers (CHGA,
# INSM1) because those legitimately appear in all endocrine subtypes.
LINEAGE_MARKERS = {
    "Endocrine_Beta":  ["IAPP", "MAFA", "NKX6-1", "PDX1", "SLC30A8", "ABCC8", "KCNJ11"],
    "Endocrine_Alpha": ["ARX", "MAFB"],
    "Endocrine_Delta": ["HHEX", "LEPR"],
    "Endocrine_Epsilon": ["GHRL"],
    "Acinar":          ["AMY1A"],
    "Ductal":          ["KRT19", "CFTR", "MUC1", "SOX9", "HNF1B"],
    "Stellate":        ["ACTA2", "PDGFRB", "RGS5", "PDGFRA"],
    "Endothelial":     ["PECAM1", "CDH5", "CLDN5", "PLVAP", "KDR", "ENG"],
    "T_cells":         ["CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "CD4", "PTPRC", "FOXP3"],
    "B_cells":         ["MS4A1", "CD79A", "CD79B", "CD19"],
    "Myeloid":         ["CD68", "CD163", "CSF1R", "MARCO", "CD14"],
    "Schwann":         ["MPZ", "SOX10", "NES"],
    "Proliferating":   ["MKI67", "TOP2A", "PCNA"],
}

# Some lineages can co-express (e.g., Endocrine subtypes share CHGA-family
# co-expression). Build groups that should be mutually exclusive.
EXCLUSIVE_GROUPS = {
    "Endocrine": ["Endocrine_Beta", "Endocrine_Alpha", "Endocrine_Delta", "Endocrine_Epsilon"],
    "Epithelial_exocrine": ["Acinar", "Ductal"],
    "Stromal": ["Stellate", "Endothelial"],
    "Lymphoid": ["T_cells", "B_cells"],
    "Myeloid": ["Myeloid"],
    "Neural": ["Schwann"],
    "Proliferating": ["Proliferating"],
}
GROUP_OF = {ln: g for g, lns in EXCLUSIVE_GROUPS.items() for ln in lns}


def lineage_signal(X, var_names, gene_idx_lookup, marker_dict):
    """Return DataFrame (n_cells x n_lineages) with mean lognorm expression."""
    out = {}
    for ln, genes in marker_dict.items():
        present = [gene_idx_lookup[g] for g in genes if g in gene_idx_lookup]
        if not present:
            out[ln] = np.zeros(X.shape[0])
            continue
        if hasattr(X, "toarray"):
            sub = np.asarray(X[:, present].mean(axis=1)).ravel()
        else:
            sub = X[:, present].mean(axis=1)
        out[ln] = sub
    return pd.DataFrame(out)


def evaluate_cell(row, threshold=0.4):
    """Per-cell: list lineages above threshold; flag if multiple GROUPS hit."""
    above = row[row > threshold].sort_values(ascending=False)
    groups = sorted(set(GROUP_OF[ln] for ln in above.index if ln in GROUP_OF))
    return above, groups


def main(sample, threshold=0.4):
    print(f"\n{'='*72}\n=== {sample} (purity threshold mean lognorm > {threshold}) ===\n{'='*72}")
    adata = sc.read_h5ad(ROOT / f"data/processed/{sample}/{sample}_phenotyped.h5ad")

    high_mask = adata.obs["n_genes_by_counts"].values > 1200
    n_high = int(high_mask.sum())
    print(f"  cells with n_genes > 1200: {n_high}")

    if n_high == 0:
        return

    sub = adata[high_mask].copy()
    if "lognorm" in sub.layers:
        X = sub.layers["lognorm"]
        layer_name = "lognorm"
    else:
        X = sub.X
        layer_name = "X (assumed log-normalized)"
    var_names = list(sub.var_names)
    gene_idx_lookup = {g: i for i, g in enumerate(var_names)}

    # Filter marker dict to only present genes
    present_markers = {
        ln: [g for g in genes if g in gene_idx_lookup]
        for ln, genes in LINEAGE_MARKERS.items()
    }
    print(f"  markers per lineage (post-panel filter):")
    for ln, gs in present_markers.items():
        print(f"    {ln:<22}: {gs}")

    df = lineage_signal(X, var_names, gene_idx_lookup, present_markers)
    df.index = sub.obs.index

    print(f"\n  per-cell evaluation (mean lognorm of each lineage panel):")
    n_clean = 0
    n_doublet = 0
    for cid, row in df.iterrows():
        celltype = sub.obs.loc[cid, "celltype"]
        n_genes = int(sub.obs.loc[cid, "n_genes_by_counts"])
        cell_area = float(sub.obs.loc[cid, "cell_area"])
        total_counts = int(sub.obs.loc[cid, "total_counts"])
        above, groups = evaluate_cell(row, threshold)
        is_doublet = len(groups) > 1
        if is_doublet:
            n_doublet += 1
        else:
            n_clean += 1
        flag = "DOUBLET" if is_doublet else "OK"
        print(f"\n  cell {cid}  celltype={celltype}  n_genes={n_genes}  "
              f"cell_area={cell_area:.0f}  total_counts={total_counts}  --> {flag}")
        if len(above):
            print(f"    lineages above threshold: {dict(above.round(3))}")
        else:
            print(f"    NO lineage above threshold "
                  f"(top-3: {dict(row.sort_values(ascending=False).head(3).round(3))})")
        if is_doublet:
            print(f"    ⚠ EXCLUSIVE-GROUP COLLISION: {groups}")

    print(f"\n  SUMMARY: {n_clean} consistent / {n_doublet} doublet-like  out of {n_high}")


for sample in ("0041323", "0041326"):
    main(sample, threshold=0.4)
