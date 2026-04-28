"""Stage 2.5 — refine rare-cell-type calls AFTER stage 02 phenotyping.

Three passes over each sample's _phenotyped.h5ad:
    A. MARGIN GATE — re-classify cells with insufficient top-vs-runner-up
       gap as Indeterminate_<top_type>.
    B. PANEL COVERAGE — re-classify cells where < COVERAGE_THRESH of the
       assigned type's panel is detected, as Indeterminate_<type>_lowcoverage.
    C. RULE-BASED rare-cell-type identification — for cells currently labeled
       any endocrine type (or Indeterminate_endocrine), apply explicit
       multi-criteria rules to identify Gamma and Epsilon cells. Multi-marker
       co-expression — cannot suffer the 1-gene argmax mislabel.
    D. SPATIAL COHERENCE — for COMMON clustering types only (NOT rare types),
       re-classify isolated cells as Indeterminate_<type>_isolated.

Run after notebook 02 has produced data/processed/{s}/{s}_phenotyped.h5ad and
finalize_panel.py has stamped the var flags + donor_sex tag.

Outputs:
    Updates _phenotyped.h5ad in-place with refined celltype/celltype_lineage.
    Adds obs columns for audit:
        margin_top_minus_second, panel_coverage_<type>,
        epsilon_marker_count, gamma_marker_count, pan_endocrine_marker_count,
        spatial_coherence_pass.
"""
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial import cKDTree

ROOT = Path("/blue/maigan/smith6jt/Xenium_Analysis")
SAMPLES = ("0041323", "0041326")

# ---- User-locked parameters (full-reanalysis-plan.md, 2026-04-28) ----
MARGIN_THRESH = 0.2
COVERAGE_THRESH = 0.30
LOG_DETECT = float(np.log1p(3.0))           # ≈ 1.386 (sample-comparable)
LOG_HIGH = float(np.log1p(5.0))             # ≈ 1.792 (high-confidence)
SPATIAL_K = 50
SPATIAL_MIN_NEIGHBORS = 2
COMMON_CLUSTERING_TYPES = {"Beta", "Alpha", "Acinar", "Ductal", "Stellate",
                            "Schwann"}
RARE_GLOBAL_THRESH = 0.05  # types below this prevalence skip spatial coherence

# ---- Multi-marker definitions ----
PAN_ENDOCRINE = ("CHGA", "INSM1", "ISL1", "NEUROD1", "FEV", "SCG2", "SCGN")
GAMMA_MARKERS = ("AQP3", "ETV1")
GAMMA_HIGH_AQP3 = "AQP3"
EPSILON_MARKER = "GHRL"
NOT_BETA_MARKERS = ("IAPP", "MAFA")
NOT_DELTA_MARKERS = ("HHEX",)
NOT_EPSILON_MARKER = "GHRL"
NOT_GAMMA_MARKERS = ("AQP3", "ETV1")

ENDO_LIKE = {"Beta", "Alpha", "Delta", "Endocrine", "Endocrine_pan",
              "Indeterminate_Beta", "Indeterminate_Alpha", "Indeterminate_Delta",
              "Indeterminate_Endocrine", "Indeterminate_endocrine"}


def get_lognorm_col(adata, gene):
    """Return (n_cells,) lognorm vector for a gene, or None if absent."""
    if gene not in adata.var_names:
        return None
    i = adata.var_names.get_loc(gene)
    L = adata.layers.get("lognorm", adata.X)
    col = L[:, i]
    if hasattr(col, "toarray"):
        col = col.toarray()
    return np.asarray(col).ravel()


def build_lognorm_matrix(adata, gene_list):
    """Build (n_cells, n_genes) dense matrix for a small gene set."""
    cols = []
    kept = []
    for g in gene_list:
        v = get_lognorm_col(adata, g)
        if v is not None:
            cols.append(v)
            kept.append(g)
    if not cols:
        return np.zeros((adata.n_obs, 0), dtype=np.float32), kept
    return np.stack(cols, axis=1).astype(np.float32), kept


def apply_margin_gate(adata):
    print(f"  pass A: margin gate (threshold {MARGIN_THRESH})")
    score_cols = [c for c in adata.obs.columns if c.endswith("_score")
                   and c not in {"highly_variable_score"}]
    if not score_cols:
        print("    no _score columns found; skipping margin gate")
        return 0
    M = np.stack([adata.obs[c].values.astype(np.float32) for c in score_cols],
                  axis=1)
    sort_idx = np.argsort(-M, axis=1)
    top = M[np.arange(len(M)), sort_idx[:, 0]]
    second = M[np.arange(len(M)), sort_idx[:, 1]]
    margin = top - second
    adata.obs["margin_top_minus_second"] = margin
    type_names = [c.replace("_score", "") for c in score_cols]
    top_type = np.array([type_names[i] for i in sort_idx[:, 0]])
    current = adata.obs["celltype"].astype(str).values.copy()
    # Only reclassify cells that match an existing _score-based assignment
    # (skip Indeterminate_* etc. which were already routed in cell 15)
    not_indeterminate = ~pd.Series(current).str.startswith("Indeterminate").values
    needs_gate = not_indeterminate & (margin < MARGIN_THRESH) & (current == top_type)
    n_reclass = int(needs_gate.sum())
    print(f"    {n_reclass:,} cells fall below margin threshold")
    new_lab = current.copy()
    new_lab[needs_gate] = ["Indeterminate_" + t for t in top_type[needs_gate]]
    adata.obs["celltype"] = pd.Categorical(new_lab)
    return n_reclass


def apply_panel_coverage(adata, marker_panels):
    print(f"  pass B: panel coverage (threshold {COVERAGE_THRESH:.0%})")
    n_reclass = 0
    current = adata.obs["celltype"].astype(str).values.copy()
    for ct, genes in marker_panels.items():
        mask = current == ct
        if mask.sum() == 0:
            continue
        L, kept = build_lognorm_matrix(adata, genes)
        if not kept:
            continue
        sub = L[mask]
        coverage = (sub >= LOG_DETECT).mean(axis=1)
        adata.obs.loc[mask, f"panel_coverage_{ct}"] = coverage
        below = mask.copy()
        below_idx = np.where(mask)[0][coverage < COVERAGE_THRESH]
        if len(below_idx) == 0:
            continue
        new_label = f"Indeterminate_{ct}_lowcoverage"
        for i in below_idx:
            current[i] = new_label
        n_reclass += len(below_idx)
        print(f"    {ct}: {len(below_idx):,} cells below {COVERAGE_THRESH:.0%} "
              f"coverage of {len(kept)}-gene panel")
    adata.obs["celltype"] = pd.Categorical(current)
    return n_reclass


def apply_rare_celltype_rules(adata):
    print(f"  pass C: rule-based Gamma + Epsilon identification")
    pan = build_lognorm_matrix(adata, PAN_ENDOCRINE)[0]
    pan_count = (pan >= LOG_DETECT).sum(axis=1) if pan.shape[1] > 0 else np.zeros(adata.n_obs, dtype=int)
    adata.obs["pan_endocrine_marker_count"] = pan_count.astype(np.int16)

    aqp3 = get_lognorm_col(adata, "AQP3")
    etv1 = get_lognorm_col(adata, "ETV1")
    ghrl = get_lognorm_col(adata, "GHRL")
    iapp = get_lognorm_col(adata, "IAPP")
    mafa = get_lognorm_col(adata, "MAFA")
    hhex = get_lognorm_col(adata, "HHEX")
    n = adata.n_obs

    def _zero_if_missing(v):
        return v if v is not None else np.zeros(n, dtype=np.float32)

    aqp3 = _zero_if_missing(aqp3)
    etv1 = _zero_if_missing(etv1)
    ghrl = _zero_if_missing(ghrl)
    iapp = _zero_if_missing(iapp)
    mafa = _zero_if_missing(mafa)
    hhex = _zero_if_missing(hhex)

    not_beta = (iapp < LOG_DETECT) & (mafa < LOG_DETECT)
    not_delta = (hhex < LOG_DETECT)

    gamma_markers_count = ((aqp3 >= LOG_DETECT).astype(int)
                            + (etv1 >= LOG_DETECT).astype(int))
    adata.obs["gamma_marker_count"] = gamma_markers_count.astype(np.int16)

    not_gamma_for_epsilon = (aqp3 < LOG_HIGH) | (etv1 < LOG_DETECT)

    epsilon_markers_count = (ghrl >= LOG_DETECT).astype(int)
    adata.obs["epsilon_marker_count"] = epsilon_markers_count.astype(np.int16)
    not_epsilon_for_gamma = (ghrl < LOG_DETECT)

    is_gamma = (
        (pan_count >= 2)
        & ((gamma_markers_count >= 2) | (aqp3 >= LOG_HIGH))
        & not_beta
        & not_delta
        & not_epsilon_for_gamma
    )
    is_epsilon = (
        (pan_count >= 2)
        & (ghrl >= LOG_DETECT)
        & not_beta
        & not_delta
        & not_gamma_for_epsilon
    )

    # Restrict to cells that the upstream argmax tagged endocrine-related.
    # Avoids hijacking acinar/ductal/immune cells with stray AQP3/ETV1.
    current = adata.obs["celltype"].astype(str).values
    endo_mask = pd.Series(current).isin(ENDO_LIKE).values
    is_gamma &= endo_mask
    is_epsilon &= endo_mask
    overlap = is_gamma & is_epsilon
    if overlap.any():
        # Resolve overlap: AQP3 strength wins over GHRL alone
        is_epsilon[overlap & (aqp3[overlap] >= LOG_HIGH)] = False
        is_gamma[overlap & ~(aqp3[overlap] >= LOG_HIGH)] = False

    new_lab = current.copy()
    new_lab[is_gamma] = "Gamma"
    new_lab[is_epsilon] = "Epsilon"
    adata.obs["celltype"] = pd.Categorical(new_lab)

    # Mirror to celltype_lineage where appropriate
    if "celltype_lineage" in adata.obs.columns:
        cl = adata.obs["celltype_lineage"].astype(str).values.copy()
        cl[is_gamma] = "Gamma"
        cl[is_epsilon] = "Epsilon"
        adata.obs["celltype_lineage"] = pd.Categorical(cl)

    n_g = int(is_gamma.sum())
    n_e = int(is_epsilon.sum())
    print(f"    Gamma: {n_g:,} cells")
    print(f"    Epsilon: {n_e:,} cells")
    return n_g, n_e


def apply_spatial_coherence(adata):
    print(f"  pass D: spatial coherence (k={SPATIAL_K}, "
          f"min_same_neighbors={SPATIAL_MIN_NEIGHBORS}, "
          f"common-clustering types only)")
    if "spatial" not in adata.obsm:
        print("    no obsm['spatial']; skipping")
        return 0
    coords = np.asarray(adata.obsm["spatial"])
    current = adata.obs["celltype"].astype(str).values.copy()
    type_counts = pd.Series(current).value_counts()
    n_total = adata.n_obs

    tree = cKDTree(coords)
    _, idx = tree.query(coords, k=SPATIAL_K + 1)
    idx = idx[:, 1:]  # drop self

    n_reclass = 0
    pass_flag = np.ones(n_total, dtype=bool)
    for ct in COMMON_CLUSTERING_TYPES:
        prevalence = type_counts.get(ct, 0) / max(n_total, 1)
        if prevalence < RARE_GLOBAL_THRESH:
            print(f"    {ct}: prevalence {prevalence:.4f} < {RARE_GLOBAL_THRESH:.0%} "
                  f"→ SKIP (rare in this sample)")
            continue
        mask = current == ct
        if mask.sum() == 0:
            continue
        same_type_neighbors = (current[idx[mask]] == ct).sum(axis=1)
        below = same_type_neighbors < SPATIAL_MIN_NEIGHBORS
        n_below = int(below.sum())
        if n_below > 0:
            below_idx = np.where(mask)[0][below]
            for i in below_idx:
                current[i] = f"Indeterminate_{ct}_isolated"
                pass_flag[i] = False
            n_reclass += n_below
            print(f"    {ct}: {n_below:,} cells reclassified isolated")
    adata.obs["celltype"] = pd.Categorical(current)
    adata.obs["spatial_coherence_pass"] = pass_flag
    return n_reclass


def process_sample(sample):
    print(f"\n{'='*72}\n=== {sample}\n{'='*72}")
    p = ROOT / f"data/processed/{sample}/{sample}_phenotyped.h5ad"
    print(f"  loading {p.name}...")
    t = time.time()
    a = sc.read_h5ad(p)
    print(f"    loaded {a.shape} in {time.time()-t:.1f}s")

    # Marker panels for panel-coverage check.
    # Pulled from notebook 02 cell 10 conventions, restricted to genes on panel.
    raw_panels = {
        "Beta":         ["IAPP", "MAFA", "NKX6-1", "PDX1", "SLC30A8", "DLK1"],
        "Alpha":        ["ARX", "MAFB", "GC", "TM4SF4", "TTR"],
        "Delta":        ["HHEX", "PCSK1", "RBP4", "BHLHE41"],
        "Acinar":       ["CTRB1", "CPA1", "CPA2", "CELA1", "PRSS1", "PRSS2",
                          "AMY2A", "AMY2B", "REG1A"],
        "Ductal":       ["KRT19", "KRT7", "CFTR", "MUC1", "SOX9"],
        "Stellate":     ["DCN", "PDGFRA", "PDGFRB", "COL1A1", "COL1A2", "COL3A1"],
        "Endothelial":  ["PECAM1", "CDH5", "VWF", "CLDN5", "FLT1", "KDR"],
        "T_cells":      ["CD3E", "CD3G", "CD2", "CD5", "PTPRC"],
        "B_cells":      ["MS4A1", "CD79A", "CD79B", "CD19", "CD22"],
        "Myeloid":      ["CD68", "CD163", "CD14", "ITGAM", "FCGR3A"],
        "Schwann":      ["S100B", "MPZ", "PMP22", "PLP1", "NRXN1"],
        "Endocrine":    list(PAN_ENDOCRINE),
        "Endocrine_pan":list(PAN_ENDOCRINE),
    }
    panels = {ct: [g for g in genes if g in a.var_names]
                for ct, genes in raw_panels.items()}

    n_a = apply_margin_gate(a)
    n_b = apply_panel_coverage(a, panels)
    n_g, n_e = apply_rare_celltype_rules(a)
    n_d = apply_spatial_coherence(a)

    # Final celltype tally
    print(f"\n  final celltype tally:")
    print(a.obs["celltype"].astype(str).value_counts().head(20).to_string())

    print(f"\n  saving updated {p.name}...")
    t = time.time()
    a.write(p, compression="gzip")
    print(f"    saved in {time.time()-t:.1f}s")

    return {
        "sample": sample,
        "margin_reclassified": n_a,
        "coverage_reclassified": n_b,
        "gamma_calls": n_g,
        "epsilon_calls": n_e,
        "spatial_isolated": n_d,
    }


def main():
    rows = []
    for s in SAMPLES:
        rows.append(process_sample(s))
    print(f"\n{'='*72}\n=== SUMMARY\n{'='*72}")
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
