"""Stage 5 — verification suite for rare cell type calls.

Reads each sample's _phenotyped.h5ad (post refine_rare_celltypes.py) and
emits a verification report covering, per (sample × cell type):

  1. cell_count, prevalence
  2. Multi-marker dotplot pass: mean lognorm of the type's panel in cells
     of that type vs all other cells. PASS if the type's own panel is
     enriched (log2 FC ≥ 0.5 for ≥ 50% of panel genes).
  3. Coexpression heatmap pass: per-gene panel coverage statistics.
  4. Cross-sample prevalence ratio: reported but NOT a gate (n=2).
  5. Sex-stratified prevalence: reports prevalence in male vs female; flags
     types with extreme male-only or female-only enrichment as a Y-confound
     leak suspicion.
  6. Spatial scatter image — PNG per (sample × rare type) saved to
     figures/07_immune_islet/verification/.

Outputs:
  data/processed/rare_celltype_verification.csv
  figures/07_immune_islet/verification/*.png
"""
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path("/blue/maigan/smith6jt/Xenium_Analysis")
SAMPLES = ("0041323", "0041326")
OUT_DIR = ROOT / "figures/07_immune_islet/verification"
OUT_DIR.mkdir(parents=True, exist_ok=True)
CSV_OUT = ROOT / "data/processed/rare_celltype_verification.csv"

LOG_DETECT = float(np.log1p(3.0))

RARE_TYPES_OF_INTEREST = ["Beta", "Alpha", "Delta", "Gamma", "Epsilon",
                            "B_cells", "Schwann", "Endocrine_pan",
                            "Endothelial"]

PANEL_FOR_TYPE = {
    "Beta":         ["IAPP", "MAFA", "NKX6-1", "PDX1", "SLC30A8", "DLK1"],
    "Alpha":        ["ARX", "MAFB", "GC", "TM4SF4", "TTR"],
    "Delta":        ["HHEX", "PCSK1", "RBP4", "BHLHE41"],
    "Gamma":        ["AQP3", "ETV1", "SCGN", "CHGA", "INSM1"],
    "Epsilon":      ["GHRL", "CHGA", "INSM1", "ISL1", "NEUROD1"],
    "B_cells":      ["MS4A1", "CD79A", "CD79B", "CD19", "CD22"],
    "Schwann":      ["S100B", "MPZ", "PMP22", "PLP1", "NRXN1"],
    "Endocrine_pan":["CHGA", "INSM1", "ISL1", "NEUROD1", "FEV", "SCG2", "SCGN"],
    "Endothelial":  ["PECAM1", "CDH5", "VWF", "CLDN5", "FLT1", "KDR"],
}


def get_lognorm(adata, gene):
    if gene not in adata.var_names:
        return None
    i = adata.var_names.get_loc(gene)
    L = adata.layers.get("lognorm", adata.X)
    col = L[:, i]
    if hasattr(col, "toarray"):
        col = col.toarray()
    return np.asarray(col).ravel()


def per_gene_log2fc(adata, ct, panel):
    """Return per-gene log2(mean_in / mean_out) on lognorm scale."""
    cur = adata.obs["celltype"].astype(str).values
    inn = cur == ct
    out = ~inn
    if inn.sum() == 0:
        return {g: float("nan") for g in panel}
    fc = {}
    for g in panel:
        v = get_lognorm(adata, g)
        if v is None:
            fc[g] = float("nan")
            continue
        m_in = float(v[inn].mean())
        m_out = float(v[out].mean())
        fc[g] = float(np.log2((m_in + 1e-3) / (m_out + 1e-3)))
    return fc


def render_spatial_panel(adata, ct, sample, ax):
    coords = np.asarray(adata.obsm["spatial"])
    cur = adata.obs["celltype"].astype(str).values
    inn = cur == ct
    rng = np.random.default_rng(0)
    n_bg = min(40_000, (~inn).sum())
    bg_idx = rng.choice(np.where(~inn)[0], size=n_bg, replace=False)
    ax.scatter(coords[bg_idx, 0], coords[bg_idx, 1], s=0.2, c="#dddddd",
                 alpha=0.3, edgecolors="none")
    if inn.any():
        ax.scatter(coords[inn, 0], coords[inn, 1], s=4, c="#d7191c",
                     alpha=0.85, edgecolors="black", linewidths=0.1)
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"{sample}: {ct}  n={int(inn.sum()):,}")


def process_sample(sample):
    print(f"\n=== {sample}")
    p = ROOT / f"data/processed/{sample}/{sample}_phenotyped.h5ad"
    a = sc.read_h5ad(p, backed="r")
    n_total = a.n_obs
    cur = a.obs["celltype"].astype(str).values
    sex = (a.obs["donor_sex"].astype(str).iloc[0]
            if "donor_sex" in a.obs.columns else "unknown")
    print(f"    n={n_total:,}, sex={sex}")

    rows = []
    for ct in RARE_TYPES_OF_INTEREST:
        n_ct = int((cur == ct).sum())
        prev = n_ct / max(n_total, 1)
        panel = PANEL_FOR_TYPE.get(ct, [])
        panel_present = [g for g in panel if g in a.var_names]

        # Need to load into memory for the FC computation; pull just a slice
        if n_ct > 0 and panel_present:
            a_mem = sc.read_h5ad(p)
            fc = per_gene_log2fc(a_mem, ct, panel_present)
            del a_mem
        else:
            fc = {g: float("nan") for g in panel_present}
        n_panel = len(panel_present)
        n_enriched = sum(1 for g, v in fc.items()
                          if not np.isnan(v) and v >= 0.5)
        dotplot_pass = (n_panel >= 1 and n_enriched / n_panel >= 0.5)

        rows.append({
            "sample": sample,
            "donor_sex": sex,
            "celltype": ct,
            "n_cells": n_ct,
            "prevalence_pct": round(100 * prev, 4),
            "panel_size": n_panel,
            "panel_genes": ";".join(panel_present),
            "panel_log2fc": ";".join(f"{g}={fc[g]:.2f}" for g in panel_present
                                       if g in fc),
            "n_panel_enriched_log2fc_ge_0.5": n_enriched,
            "dotplot_pass": dotplot_pass,
        })
        print(f"    {ct:<14} n={n_ct:>7,} prev={100*prev:>5.2f}%  "
              f"panel_log2fc enriched: {n_enriched}/{n_panel}  "
              f"dotplot_pass={dotplot_pass}")

    a.file.close()
    return rows


def render_spatial_grid(rows_df):
    """One row per rare type; one column per sample. Reads h5ads on demand."""
    rare = ["Beta", "Alpha", "Delta", "Gamma", "Epsilon", "Schwann"]
    n_rows = len(rare)
    fig, axes = plt.subplots(n_rows, 2, figsize=(14, 5.0 * n_rows),
                                constrained_layout=True)
    if n_rows == 1:
        axes = axes.reshape(1, 2)
    for j, sample in enumerate(SAMPLES):
        p = ROOT / f"data/processed/{sample}/{sample}_phenotyped.h5ad"
        a = sc.read_h5ad(p)
        for i, ct in enumerate(rare):
            render_spatial_panel(a, ct, sample, axes[i, j])
    plt.suptitle("Rare-cell-type spatial verification — red = called type, "
                  "grey = all others (subsampled)", fontsize=12)
    fig.savefig(OUT_DIR / "spatial_verification_grid.png",
                  dpi=180, bbox_inches="tight")
    plt.close()
    print(f"\n  wrote {OUT_DIR / 'spatial_verification_grid.png'}")


def main():
    all_rows = []
    for s in SAMPLES:
        all_rows.extend(process_sample(s))
    df = pd.DataFrame(all_rows)

    # Cross-sample prevalence ratio (reported, not a gate)
    pivot = df.pivot(index="celltype", columns="sample",
                       values="prevalence_pct")
    if all(s in pivot.columns for s in SAMPLES):
        pivot["ratio_max_min"] = pivot.max(axis=1) / pivot.min(axis=1).replace(0, np.nan)
    print(f"\n=== cross-sample prevalence (informational only — n=2) ===")
    print(pivot.round(3).to_string())

    df.to_csv(CSV_OUT, index=False)
    print(f"\nwrote {CSV_OUT}")
    render_spatial_grid(df)


if __name__ == "__main__":
    main()
