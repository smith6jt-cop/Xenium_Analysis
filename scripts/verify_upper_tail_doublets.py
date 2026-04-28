"""Doublet vs real-cell validation for the upper tail of n_genes.

Approach: take cells in the top 1% by n_genes_by_counts. For each, extract
LOGNORM (layers["lognorm"]) values of canary genes that are highly lineage-specific
and biologically mutually exclusive. A real big healthy cell should detect
at most ONE lineage group above background. A segmentation-merge doublet
will detect multiple.

Detection rule per canary gene: lognorm >= log1p(3) ≈ 1.386 — sample-comparable
because lognorm = log1p(normalize_total(counts)) rescales each cell to the
sample's median library size before log-transformation. Equivalent to
"raw count >= 3 in a median-depth cell" but consistent across samples with
different overall sequencing depth.

Lineage groups (mutually exclusive — co-detection is a doublet signal):
  - Acinar:        AMY1A, CUZD1
  - Ductal:        KRT19, CFTR, MUC1
  - Beta:          IAPP, MAFA, NKX6-1, PDX1, SLC30A8
  - Alpha:         ARX, MAFB
  - Endocrine_pan: CHGA, INSM1, ISL1, NEUROD1, FEV
  - Stellate:      ACTA2, PDGFRB, RGS5
  - Endothelial:   PECAM1, CDH5, CLDN5, PLVAP
  - T_cells:       CD3D, CD3E, CD8A, PTPRC
  - B_cells:       MS4A1, CD79A
  - Myeloid:       CD68, CD163, MARCO
  - Schwann:       MPZ, SOX10

Beta/Alpha/Endocrine_pan are NOT mutually exclusive (all islet endocrine).
Ductal/Acinar are NOT mutually exclusive (both pancreatic exocrine epithelium).
Everything else IS mutually exclusive with everything else.
"""
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path("/blue/maigan/smith6jt/Xenium_Analysis")

CANARY_PANELS = {
    "Acinar":        ["AMY1A", "CUZD1"],
    "Ductal":        ["KRT19", "CFTR", "MUC1"],
    "Beta":          ["IAPP", "MAFA", "NKX6-1", "PDX1", "SLC30A8"],
    "Alpha":         ["ARX", "MAFB"],
    "Endocrine_pan": ["CHGA", "INSM1", "ISL1", "NEUROD1", "FEV"],
    "Stellate":      ["ACTA2", "PDGFRB", "RGS5"],
    "Endothelial":   ["PECAM1", "CDH5", "CLDN5", "PLVAP"],
    "T_cells":       ["CD3D", "CD3E", "CD8A", "PTPRC"],
    "B_cells":       ["MS4A1", "CD79A"],
    "Myeloid":       ["CD68", "CD163", "MARCO"],
    "Schwann":       ["MPZ", "SOX10"],
}

# Mutual-exclusion groups — co-detection across distinct groups = doublet.
# Endocrine_pan / Beta / Alpha share an islet identity (same group).
# Ductal / Acinar share an exocrine epithelial identity (same group).
EXCLUSIVE_GROUP = {
    "Acinar": "Exocrine_epi",
    "Ductal": "Exocrine_epi",
    "Beta": "Endocrine_islet",
    "Alpha": "Endocrine_islet",
    "Endocrine_pan": "Endocrine_islet",
    "Stellate": "Stromal",
    "Endothelial": "Endothelial",  # endothelium != stellate; keep separate
    "T_cells": "Lymphoid",
    "B_cells": "Lymphoid",
    "Myeloid": "Myeloid",
    "Schwann": "Neural",
}

# Detection threshold on the LOGNORM layer (sample-comparable).
# log1p(3) ≈ 1.386 maps to "raw count >= 3 in a median-depth cell" but works
# uniformly across samples regardless of sequencing depth.
LOG_DETECT = np.log1p(3.0)  # ≈ 1.386
MIN_DETECT_TRANSCRIPTS = LOG_DETECT  # alias kept for backward-compat in the call sites
LOWER_TAIL_QUANTILE = 0.99  # top 1%


def load(sample):
    return sc.read_h5ad(ROOT / f"data/processed/{sample}/{sample}_phenotyped.h5ad")


def lineage_detected(raw_counts_row, gene_idx_lookup, marker_dict, min_detect):
    """Return dict {lineage: bool} — True iff any marker has raw count >= min_detect."""
    detected = {}
    for lin, genes in marker_dict.items():
        idxs = [gene_idx_lookup[g] for g in genes if g in gene_idx_lookup]
        if not idxs:
            detected[lin] = False
            continue
        vals = raw_counts_row[idxs]
        if hasattr(vals, "toarray"):
            vals = np.asarray(vals.toarray()).ravel()
        detected[lin] = bool(np.any(vals >= min_detect))
    return detected


def lineage_max_count(raw_counts_row, gene_idx_lookup, marker_dict):
    """Return dict {lineage: (gene, max_raw_count)}."""
    out = {}
    for lin, genes in marker_dict.items():
        best = (None, 0)
        for g in genes:
            if g not in gene_idx_lookup:
                continue
            v = raw_counts_row[gene_idx_lookup[g]]
            if hasattr(v, "toarray"):
                v = float(np.asarray(v.toarray()).ravel()[0])
            else:
                v = float(v)
            if v > best[1]:
                best = (g, v)
        out[lin] = best
    return out


def main(sample, min_detect=MIN_DETECT_TRANSCRIPTS, q=LOWER_TAIL_QUANTILE):
    print(f"\n{'='*72}\n=== {sample}  top {(1-q)*100:.0f}%  min_detect={min_detect} ===\n{'='*72}")
    adata = load(sample)
    n = adata.n_obs

    # Top tail
    n_genes = adata.obs["n_genes_by_counts"].values
    cutoff = float(np.quantile(n_genes, q))
    mask = n_genes > cutoff
    n_high = int(mask.sum())
    print(f"  n_genes p{int(q*100)} cutoff: {cutoff:.0f}  → {n_high:,} cells in upper tail "
          f"({100*n_high/n:.2f}% of {n:,})")

    sub = adata[mask].copy()
    var_names = list(sub.var_names)
    gene_idx_lookup = {g: i for i, g in enumerate(var_names)}

    # Build a small dense canary matrix: n_high cells x N_canary genes.
    all_canary_genes = sorted({g for genes in CANARY_PANELS.values() for g in genes})
    present_canary = [g for g in all_canary_genes if g in gene_idx_lookup]
    canary_idx = [gene_idx_lookup[g] for g in present_canary]
    # Use lognorm layer for sample-comparable detection (compute on the fly
    # if missing from the h5ad).
    if "lognorm" in sub.layers:
        counts = sub.layers["lognorm"]
    else:
        sub.X = sub.layers["counts"].copy()
        sc.pp.normalize_total(sub)
        sc.pp.log1p(sub)
        counts = sub.X
    if hasattr(counts, "toarray"):
        canary_dense = np.asarray(counts[:, canary_idx].toarray())
    else:
        canary_dense = np.asarray(counts[:, canary_idx])
    canary_pos = {g: j for j, g in enumerate(present_canary)}
    print(f"  built canary_dense: shape={canary_dense.shape}  (cells x canary_genes)")

    print(f"  evaluating {n_high:,} cells against {len(CANARY_PANELS)} canary lineages...")

    # Per-cell detection
    def cell_lineage_detected(canary_row, marker_dict, min_detect):
        out = {}
        for lin, genes in marker_dict.items():
            present = [canary_pos[g] for g in genes if g in canary_pos]
            if not present:
                out[lin] = False
                continue
            out[lin] = bool(np.any(canary_row[present] >= min_detect))
        return out

    rows = []
    for i in range(sub.n_obs):
        det = cell_lineage_detected(canary_dense[i], CANARY_PANELS, min_detect)
        groups_hit = sorted({EXCLUSIVE_GROUP[lin] for lin, d in det.items() if d})
        rows.append({
            "cell_id": sub.obs.index[i],
            "celltype": sub.obs["celltype"].iloc[i],
            "n_genes": int(sub.obs["n_genes_by_counts"].iloc[i]),
            "total_counts": int(sub.obs["total_counts"].iloc[i]),
            "cell_area": float(sub.obs["cell_area"].iloc[i]),
            "lineages_detected": [lin for lin, d in det.items() if d],
            "exclusive_groups_hit": groups_hit,
            "n_groups_hit": len(groups_hit),
        })
    df = pd.DataFrame(rows)

    print(f"\n  distribution of n_groups_hit (mutually-exclusive groups co-detected):")
    for k, v in df["n_groups_hit"].value_counts().sort_index().items():
        print(f"    {k} group(s): {v:,}  ({100*v/n_high:.2f}%)")

    print(f"\n  n_groups_hit by assigned celltype (top {n_high:,} cells):")
    pivot = df.pivot_table(
        index="celltype", columns="n_groups_hit", aggfunc="size", fill_value=0
    )
    print(pivot.to_string())

    # Doublet candidates: 3+ exclusive groups (more than just exocrine/endocrine pair)
    doublet_strict = df[df["n_groups_hit"] >= 3]
    doublet_loose = df[df["n_groups_hit"] >= 2]
    print(f"\n  STRICT doublets (>=3 mutually-exclusive groups): {len(doublet_strict):,}")
    print(f"  LOOSE  doublets (>=2 mutually-exclusive groups): {len(doublet_loose):,}")

    # Sample a few examples per category for inspection
    if len(doublet_strict) > 0:
        print(f"\n  --- 5 strict doublet examples ---")
        for _, r in doublet_strict.sample(min(5, len(doublet_strict)), random_state=0).iterrows():
            print(f"    cell {r.cell_id}  ct={r.celltype}  n_genes={r.n_genes}  "
                  f"counts={r.total_counts}  area={r.cell_area:.0f}")
            print(f"      groups={r.exclusive_groups_hit}  lineages={r.lineages_detected}")

    # For all upper-tail cells, look at LINEAGE GROUPS distribution to characterize
    print(f"\n  cells where the assigned celltype's lineage WAS detected:")
    # For each cell, check if its assigned celltype name maps to a CANARY_PANELS key (or close)
    label_to_canary = {
        "Acinar": "Acinar", "Ductal": "Ductal",
        "Beta": "Beta", "Alpha": "Alpha", "Delta": "Endocrine_pan",
        "Endocrine": "Endocrine_pan", "Epsilon": "Endocrine_pan",
        "Stellate": "Stellate", "Endothelial": "Endothelial",
        "T_cells": "T_cells", "B_cells": "B_cells",
        "Myeloid": "Myeloid", "Schwann": "Schwann",
    }
    df["expected_lineage"] = df["celltype"].map(label_to_canary)
    df["expected_detected"] = df.apply(
        lambda r: r["expected_lineage"] in r["lineages_detected"]
        if isinstance(r["expected_lineage"], str) else False,
        axis=1,
    )
    pct = 100 * df["expected_detected"].sum() / max(len(df), 1)
    print(f"    {df['expected_detected'].sum():,} / {len(df):,}  ({pct:.1f}%)")

    return df


for sample in ("0041323", "0041326"):
    main(sample)

print("\n" + "="*72)
print("DOUBLET VALIDATION COMPLETE")
print("="*72)
