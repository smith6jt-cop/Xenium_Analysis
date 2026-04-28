"""Comprehensive figure suite for the immune+islet analysis.

Outputs to figures/07_immune_islet/. Recomputes dist_to_islet using the
endo-label seed (consistent with the new headline tables) and writes
the updated columns back to data/processed/{sample}/{sample}_immune_phenotyped.h5ad.

Figures (designed to be reviewer-grade):

  01_insulitis_grade_summary.png   stacked bars per sample
  02_density_enrichment_heatmap.png subtype × sample
  03_distance_cdf.png              CDF of distance-to-islet per subtype
  04_islet_size_vs_immune.png      n_endocrine vs total_immune per islet
  05_top_insulitis_islets.png      spatial zoom-in of top 6 insulitis islets/sample
  06_immune_marker_dotplot.png     all-subtype phenotyping dotplot
  07_tcell_marker_dotplot.png      T-cell canonical-marker dotplot
  08_tcell_scores_umap.png         T-cell UMAP × cytotoxic/exhausted/reg/activation
  09_tcell_distance_violin.png     distance-to-islet violin per subtype
  10_tcell_de_islet_vs_distal.png  T-cell DE: in-islet vs distal lognorm logFC
  11_spatial_overview.png          tissue map with islets + immune cells
  12_islet_immune_composition.png  stacked composition per top islets
"""
import time
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore", category=FutureWarning)

ROOT = Path("/blue/maigan/smith6jt/Xenium_Analysis")
OUT = ROOT / "figures/07_immune_islet"
OUT.mkdir(parents=True, exist_ok=True)
SAMPLES = ("0041323", "0041326")
ENDO_LABELS = {"Endocrine", "Beta", "Alpha", "Delta", "Endocrine_pan"}
EPS_UM = 50.0
MIN_SAMPLES = 10
ISLET_ZONE_UM = 50.0
PROXIMAL_UM = 200.0

# ===== Marker panels =====
TCELL_MARKERS = {
    "Pan-T":       ["CD3E", "CD3G", "CD2", "CD5", "PTPRC"],
    "CD4":         ["CD4", "IL7R", "CCR7"],
    "CD8":         ["CD8A", "CD8B"],
    "Cytotoxic":   ["GZMA", "GZMB", "GZMK", "GZMH", "PRF1", "NKG7", "IFNG"],
    "Regulatory":  ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "ENTPD1"],
    "Exhausted":   ["PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX", "EOMES"],
    "Activated":   ["MKI67", "TNFRSF9", "CD69", "ICOS"],
    "Memory":      ["SELL", "IL7R"],
    "Resident":    ["ITGAE", "ITGA1", "CXCR6", "RUNX3"],
    "Tfh":         ["CXCR5", "BCL6"],
}

ALL_IMMUNE_MARKERS = {
    "T-cells":   ["CD3E", "CD3G", "CD2", "CD5", "PTPRC", "CD4", "CD8A", "CD8B"],
    "T_reg":     ["FOXP3", "IL2RA", "CTLA4", "IKZF2"],
    "T_cyto":    ["GZMB", "GZMK", "PRF1", "IFNG", "NKG7"],
    "T_exh":     ["PDCD1", "HAVCR2", "LAG3", "TOX"],
    "NK":        ["NCAM1", "KLRD1", "KLRF1", "NCR1"],
    "B":         ["MS4A1", "CD79A", "CD79B", "CD19", "CD22"],
    "Plasma":    ["XBP1", "PRDM1", "MZB1"],
    "Mono/Mac":  ["CD68", "CD163", "MARCO", "CD14", "FCGR3A", "ITGAM"],
    "M1":        ["TNF", "IL6", "IL1B", "CXCL10", "CXCL9"],
    "M2":        ["MRC1", "CCL18", "CCL22"],
    "DC":        ["ITGAX", "CLEC9A", "CD1C", "LAMP3"],
}

TCELL_SCORE_GENES = {
    "Cytotoxicity":  ["GZMA", "GZMB", "GZMK", "GZMH", "PRF1", "NKG7", "IFNG"],
    "Exhaustion":    ["PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX", "EOMES"],
    "Regulation":    ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "ENTPD1"],
    "Activation":    ["MKI67", "TNFRSF9", "CD69", "ICOS"],
    "Tissue_residency": ["ITGAE", "ITGA1", "CXCR6", "RUNX3"],
}

# Curated immune/T-cell-relevant gene set for the in-islet vs distal DE.
# Without this restriction, the top hits get polluted by spatial-
# contamination genes (CHGA from neighboring β-cells, AMY1A/CUZD1 from
# acinar transcripts diffusing into distal-region segments). Drawing
# only from this list keeps the DE biologically interpretable as a
# T-cell program shift.
TCELL_DE_GENES = sorted(set(sum([
    # TCR / lineage
    ["CD3E", "CD3G", "CD2", "CD5", "CD4", "CD8A", "CD8B", "PTPRC",
     "ZAP70", "LCK", "FYN", "ITK", "CD7", "TRAC", "TRBC2"],
    # Cytotoxic
    ["GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "NKG7", "IFNG", "TNF",
     "FASLG", "FAS"],
    # Regulatory
    ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "ENTPD1", "TIGIT", "TGFB1",
     "IL10"],
    # Exhaustion
    ["PDCD1", "HAVCR2", "LAG3", "TOX", "EOMES", "TBX21", "BTLA"],
    # Costimulatory / activation
    ["ICOS", "TNFRSF9", "TNFRSF18", "CD28", "CD27", "CD40LG", "CD69"],
    # Chemokines / receptors
    ["CXCL9", "CXCL10", "CXCL11", "CXCL12", "CCL5", "CCL19", "CCL21",
     "IL16", "CCR4", "CCR5", "CCR6", "CCR7", "CCR9", "CXCR3", "CXCR4",
     "CXCR5", "CXCR6", "S1PR1"],
    # Cytokines
    ["IL4", "IL7", "IL10", "IL15", "IL17A", "IL21", "IFNG", "TNF",
     "TGFB1", "TGFB2", "TGFB3"],
    # IFN response / antigen presentation
    ["ISG15", "MX1", "IFIT1", "IFIT2", "IFIT3", "IRF1", "IRF3", "IRF7",
     "IRF8", "STAT1", "STAT3", "OAS1", "OAS3", "B2M"],
    # Adhesion / migration
    ["ITGAE", "ITGA1", "ITGA4", "ITGAL", "ITGB2", "ITGB7", "SELL",
     "CD44", "ITGAM"],
    # Survival / proliferation
    ["BCL2", "BCL2L1", "MKI67", "MCM2", "CCNB1", "TOP2A"],
    # Memory
    ["IL7R", "GPR183"],
    # Th polarization
    ["GATA3", "RORC", "BCL6"],
], [])))


def get_seed_xy(sample):
    """Return endocrine seed coords (clustered + cluster labels) for one sample."""
    p = ROOT / f"data/processed/{sample}/{sample}_phenotyped.h5ad"
    print(f"    loading seeds from {p.name}...")
    a = sc.read_h5ad(p, backed="r")
    coords = np.asarray(a.obsm["spatial"])
    ct = a.obs["celltype_lineage"].astype(str).values
    seed = coords[np.isin(ct, list(ENDO_LABELS))]
    a.file.close()
    from sklearn.cluster import DBSCAN
    labels = DBSCAN(eps=EPS_UM, min_samples=MIN_SAMPLES,
                     n_jobs=-1).fit_predict(seed)
    keep = labels >= 0
    return seed[keep], labels[keep]


def recompute_distances():
    """Recompute dist_to_islet_um + nearest_islet_id with endo-label seed,
    write back to immune h5ads."""
    from sklearn.neighbors import NearestNeighbors
    print("\n=== Recomputing dist_to_islet using endo-label seed ===")
    for s in SAMPLES:
        print(f"  {s}:")
        seed_xy, seed_lab = get_seed_xy(s)
        print(f"    {seed_xy.shape[0]:,} clustered seeds, "
              f"{seed_lab.max()+1} islets")
        p_imm = ROOT / f"data/processed/{s}/{s}_immune_phenotyped.h5ad"
        a = sc.read_h5ad(p_imm)
        if "sample" in a.obs.columns:
            mask = a.obs["sample"].astype(str) == s
            sub = a[mask].copy()
        else:
            sub = a
        nn = NearestNeighbors(n_neighbors=1, n_jobs=-1).fit(seed_xy)
        d, idx = nn.kneighbors(np.asarray(sub.obsm["spatial"]))
        d = d.ravel(); idx = idx.ravel()
        islet_id = [f"{s}_islet_{seed_lab[k]:04d}" for k in idx]
        bins = np.full(len(d), "distal", dtype=object)
        bins[d < PROXIMAL_UM] = "proximal"
        bins[d < ISLET_ZONE_UM] = "intra_or_peri"
        sub.obs["dist_to_islet_um"] = d
        sub.obs["nearest_islet_id"] = pd.Categorical(islet_id)
        sub.obs["distance_bin"] = pd.Categorical(
            bins, categories=["intra_or_peri", "proximal", "distal"])
        sub.write(p_imm, compression="gzip")
        print(f"    wrote updated {p_imm.name}")


def load_immune():
    """Load both immune h5ads, return concatenated."""
    parts = []
    for s in SAMPLES:
        p = ROOT / f"data/processed/{s}/{s}_immune_phenotyped.h5ad"
        a = sc.read_h5ad(p)
        if "sample" in a.obs.columns:
            a = a[a.obs["sample"].astype(str) == s].copy()
        else:
            a.obs["sample"] = s
        parts.append(a)
    if len(parts) > 1:
        full = parts[0].concatenate(parts[1:], batch_key=None,
                                       index_unique=None,
                                       join="inner")
    else:
        full = parts[0]
    return full


# ===== Figure helpers =====
def setup_fig(nrows, ncols, w, h):
    fig, axes = plt.subplots(nrows, ncols, figsize=(w, h),
                                constrained_layout=True)
    return fig, axes


# ===== 01: Insulitis grade summary =====
def fig01_insulitis_grades():
    print("\n=== fig 01: insulitis grade summary ===")
    g = pd.read_csv(ROOT / "data/processed/islet_insulitis_grades.csv",
                     dtype={"sample": str})
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5),
                                       constrained_layout=True)
    cols = ["no_insulitis", "peri_insulitis", "insulitis"]
    colors = ["#cccccc", "#fdae61", "#d7191c"]
    bottom = np.zeros(len(g))
    x = np.arange(len(g))
    for c, col in zip(cols, colors):
        ax1.bar(x, g[c].values, bottom=bottom, color=col, label=c,
                  edgecolor="black", linewidth=0.5)
        bottom += g[c].values
    ax1.set_xticks(x)
    ax1.set_xticklabels(g["sample"])
    ax1.set_ylabel("Number of islets")
    ax1.set_title("Insulitis grade per sample (counts)")
    ax1.legend(loc="upper left", framealpha=0.95)
    # Right: percentages
    pct = g[cols].div(g[cols].sum(axis=1), axis=0) * 100
    bottom = np.zeros(len(g))
    for c, col in zip(cols, colors):
        ax2.bar(x, pct[c].values, bottom=bottom, color=col, label=c,
                 edgecolor="black", linewidth=0.5)
        bottom += pct[c].values
    ax2.set_xticks(x)
    ax2.set_xticklabels(g["sample"])
    ax2.set_ylabel("% of islets")
    ax2.set_title("Insulitis grade per sample (percent)")
    for i, row in g.iterrows():
        ax2.text(i, 105, f"n={int(row['total_islets'])}",
                  ha="center", fontsize=10)
    ax2.set_ylim(0, 115)
    plt.suptitle("Clinical insulitis grading (Campbell-Thompson 2013, "
                  "≥6 immune cells within 50 μm = insulitis)",
                  fontsize=11)
    fig.savefig(OUT / "01_insulitis_grade_summary.png",
                  dpi=200, bbox_inches="tight")
    plt.close()


# ===== 02: Density enrichment heatmap =====
def fig02_density_enrichment():
    print("=== fig 02: density enrichment heatmap ===")
    e = pd.read_csv(ROOT / "data/processed/immune_proximity_summary.csv",
                     dtype={"sample": str})
    pivot = e.pivot(index="immune_subtype", columns="sample",
                       values="density_enrichment")
    # Order subtypes by mean enrichment
    pivot = pivot.loc[pivot.mean(axis=1).sort_values(ascending=False).index]
    fig, ax = plt.subplots(figsize=(7, 7), constrained_layout=True)
    log_pivot = np.log2(pivot.clip(lower=0.01))
    vmax = max(abs(log_pivot.min().min()), abs(log_pivot.max().max()))
    im = ax.imshow(log_pivot.values, aspect="auto", cmap="RdBu_r",
                    vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns)
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index)
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            v = pivot.values[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                     color="black" if abs(np.log2(max(v, .01))) < vmax*0.6 else "white",
                     fontsize=9)
    cb = plt.colorbar(im, ax=ax, label="log₂(density enrichment)\n"
                                          "(islet-zone vs distal)")
    ax.set_title("Density enrichment per immune subtype × sample\n"
                  "values shown are linear-scale enrichment ratios")
    fig.savefig(OUT / "02_density_enrichment_heatmap.png",
                  dpi=200, bbox_inches="tight")
    plt.close()


# ===== 03: Distance CDF =====
def fig03_distance_cdf(adata):
    print("=== fig 03: distance CDF per subtype ===")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True,
                                constrained_layout=True)
    subtypes = sorted(adata.obs["immune_subtype"].unique())
    cmap = plt.get_cmap("tab20")
    for ax, s in zip(axes, SAMPLES):
        m = adata.obs["sample"].astype(str) == s
        for i, st in enumerate(subtypes):
            mm = m & (adata.obs["immune_subtype"].astype(str) == st)
            d = adata.obs.loc[mm, "dist_to_islet_um"].values
            if len(d) == 0:
                continue
            d_sort = np.sort(d)
            cdf = np.arange(1, len(d)+1) / len(d)
            ax.plot(d_sort, cdf, label=f"{st} (n={len(d):,})",
                     color=cmap(i % 20), lw=1.6)
        ax.axvline(ISLET_ZONE_UM, ls="--", color="grey", alpha=0.7)
        ax.axvline(PROXIMAL_UM, ls=":", color="grey", alpha=0.7)
        ax.set_xlim(0, 600)
        ax.set_xlabel("Distance to nearest endocrine seed (μm)")
        ax.set_ylabel("Cumulative fraction")
        ax.set_title(f"{s}")
        ax.grid(True, alpha=0.3)
    axes[1].legend(bbox_to_anchor=(1.02, 1), loc="upper left",
                     fontsize=8, framealpha=0.95)
    plt.suptitle("Per-subtype CDF of distance to nearest islet — "
                  "left-shift = subtype concentrates near islets",
                  fontsize=11)
    fig.savefig(OUT / "03_distance_cdf.png",
                  dpi=200, bbox_inches="tight")
    plt.close()


# ===== 04: islet size vs immune burden =====
def fig04_islet_size_immune():
    print("=== fig 04: islet size vs immune burden ===")
    df = pd.read_csv(ROOT / "data/processed/islet_infiltration_per100endo.csv",
                      dtype={"sample": str})
    fig, axes = plt.subplots(1, 2, figsize=(14, 6),
                                constrained_layout=True)
    grade_color = {"no_insulitis": "#cccccc",
                     "peri_insulitis": "#fdae61",
                     "insulitis": "#d7191c"}
    for ax, s in zip(axes, SAMPLES):
        sub = df[df["sample"] == s]
        for grade, color in grade_color.items():
            mm = sub["insulitis_grade"] == grade
            ax.scatter(sub.loc[mm, "n_endocrine"],
                         sub.loc[mm, "total_immune"],
                         s=20, c=color, alpha=0.7,
                         label=f"{grade} (n={int(mm.sum())})",
                         edgecolors="black", linewidths=0.3)
        ax.set_xscale("log")
        ax.set_yscale("symlog", linthresh=1)
        ax.set_xlabel("Endocrine cells in islet (n)")
        ax.set_ylabel("Immune cells within 50 μm")
        ax.set_title(f"{s} — n_islets={len(sub)}")
        ax.legend(loc="upper left", framealpha=0.95, fontsize=9)
        ax.grid(True, alpha=0.3)
    plt.suptitle("Per-islet immune burden vs islet size",
                  fontsize=11)
    fig.savefig(OUT / "04_islet_size_vs_immune.png",
                  dpi=200, bbox_inches="tight")
    plt.close()


# ===== 05: Top insulitis islets spatial zoom =====
def fig05_top_insulitis(adata):
    print("=== fig 05: top insulitis islets zoom ===")
    df = pd.read_csv(ROOT / "data/processed/islet_infiltration_per100endo.csv",
                      dtype={"sample": str})
    n_show = 6
    fig, axes = plt.subplots(2, n_show, figsize=(n_show*4.0, 9),
                                constrained_layout=True)
    for row_idx, s in enumerate(SAMPLES):
        seed_xy, seed_lab = get_seed_xy(s)
        sub = df[df["sample"] == s].nlargest(n_show, "total_immune")
        # Subset immune of this sample
        m_imm = adata.obs["sample"].astype(str) == s
        imm_xy = adata.obsm["spatial"][m_imm]
        imm_subtype = adata.obs.loc[m_imm, "immune_subtype"].values
        for j, (_, isl) in enumerate(sub.iterrows()):
            ax = axes[row_idx, j]
            cx, cy = isl["centroid_x"], isl["centroid_y"]
            R = 150
            # Endocrine seeds in window
            d_seed = np.sqrt((seed_xy[:, 0]-cx)**2 + (seed_xy[:, 1]-cy)**2)
            mseed = d_seed < R
            ax.scatter(seed_xy[mseed, 0], seed_xy[mseed, 1],
                         s=8, c="#ffd966", edgecolors="#806600",
                         linewidths=0.2, label="endocrine", zorder=2)
            # Immune cells in window, colored by subtype family
            d_imm = np.sqrt((imm_xy[:, 0]-cx)**2 + (imm_xy[:, 1]-cy)**2)
            mimm = d_imm < R
            family = np.array([
                "T_cytotoxic" if x == "T_cytotoxic"
                else "T_reg" if x == "T_reg"
                else "T_other" if x.startswith("T_")
                else "B/Plasma" if x.startswith("B_")
                else "Macro_M2" if x == "Macro_M2"
                else "Macro_M1" if x == "Macro_M1"
                else "NK" if x == "NK"
                else "Other" for x in imm_subtype])
            fam_color = {"T_cytotoxic": "#d7191c",
                          "T_reg": "#fdae61",
                          "T_other": "#fee08b",
                          "B/Plasma": "#abdda4",
                          "Macro_M1": "#2b83ba",
                          "Macro_M2": "#5e3c99",
                          "NK": "#e34a33",
                          "Other": "#999999"}
            for fam, c in fam_color.items():
                mm = mimm & (family == fam)
                if mm.sum() == 0:
                    continue
                ax.scatter(imm_xy[mm, 0], imm_xy[mm, 1],
                             s=14, c=c, edgecolors="black",
                             linewidths=0.2,
                             label=f"{fam} ({int(mm.sum())})",
                             zorder=3)
            ax.set_xlim(cx-R, cx+R)
            ax.set_ylim(cy-R, cy+R)
            ax.set_aspect("equal")
            ax.set_xticks([]); ax.set_yticks([])
            ax.set_title(f"islet {isl['islet_id'].split('_')[-1]}\n"
                          f"endo={int(isl['n_endocrine'])}  "
                          f"imm={int(isl['total_immune'])}",
                          fontsize=11)
            if j == 0:
                ax.set_ylabel(s, fontsize=13, fontweight="bold")
        # Legend on the last column
        axes[row_idx, -1].legend(bbox_to_anchor=(1.02, 1),
                                       loc="upper left", fontsize=9)
    plt.suptitle("Top insulitis islets — 300 μm window",
                  fontsize=11)
    fig.savefig(OUT / "05_top_insulitis_islets.png",
                  dpi=200, bbox_inches="tight")
    plt.close()


# ===== 06 / 07: Marker dotplots =====
def fig06_immune_marker_dotplot(adata):
    print("=== fig 06: all-immune marker dotplot ===")
    panel = set(adata.var_names)
    use = {k: [g for g in v if g in panel] for k, v in ALL_IMMUNE_MARKERS.items()}
    use = {k: v for k, v in use.items() if v}
    flat = [g for v in use.values() for g in v]
    sc.settings.set_figure_params(dpi=200)
    dp = sc.pl.dotplot(adata, var_names=use, groupby="immune_subtype",
                          use_raw=False, layer="lognorm",
                          standard_scale="var",
                          return_fig=True,
                          figsize=(max(13, len(flat)*0.27), 7))
    dp.make_figure()
    fig = plt.gcf()
    mainplot_ax = dp.ax_dict.get("mainplot_ax")
    if mainplot_ax is not None:
        cum = 0
        keys = list(use.keys())
        for k in keys[:-1]:
            cum += len(use[k])
            mainplot_ax.axvline(cum, color="black", lw=0.9, alpha=0.7,
                                  zorder=10)
    fig.text(0.5, 1.10,
             "Immune marker expression by subtype "
             "(both samples combined, lognorm, scaled per gene)",
             ha="center", fontsize=11, transform=fig.transFigure)
    fig.savefig(OUT / "06_immune_marker_dotplot.png",
                  dpi=200, bbox_inches="tight")
    plt.close()


def fig07_tcell_marker_dotplot(adata):
    print("=== fig 07: T-cell marker dotplot ===")
    tmask = adata.obs["immune_subtype"].astype(str).str.startswith("T_")
    tadata = adata[tmask].copy()
    panel = set(tadata.var_names)
    use = {k: [g for g in v if g in panel] for k, v in TCELL_MARKERS.items()}
    use = {k: v for k, v in use.items() if v}
    sc.settings.set_figure_params(dpi=200)
    dp = sc.pl.dotplot(tadata, var_names=use,
                          groupby="immune_subtype",
                          use_raw=False, layer="lognorm",
                          standard_scale="var",
                          return_fig=True,
                          figsize=(16, 6))
    dp.make_figure()
    fig = plt.gcf()
    mainplot_ax = dp.ax_dict.get("mainplot_ax")
    if mainplot_ax is not None:
        cum = 0
        keys = list(use.keys())
        for k in keys[:-1]:
            cum += len(use[k])
            mainplot_ax.axvline(cum, color="black", lw=0.9, alpha=0.7,
                                  zorder=10)
    fig.text(0.5, 1.30,
             "T-cell marker expression by subtype "
             "(canonical T-cell programs, lognorm, scaled per gene)",
             ha="center", fontsize=12, transform=fig.transFigure)
    fig.savefig(OUT / "07_tcell_marker_dotplot.png",
                  dpi=200, bbox_inches="tight")
    plt.close()


# ===== 08: T-cell scores UMAP =====
def fig08_tcell_scores(adata):
    print("=== fig 08: T-cell scores UMAP ===")
    tmask = adata.obs["immune_subtype"].astype(str).str.startswith("T_")
    tad = adata[tmask].copy()
    print(f"    T-cells: {tad.n_obs:,}")
    # Quality gate: drop the noisy tail that produces UMAP scatter.
    # Low-confidence T-cells (mostly indeterminate-with-T-marker pickups)
    # do not have coherent T-cell programs and embed as spray.
    if "T_cells_score" in tad.obs.columns:
        thresh = float(np.quantile(tad.obs["T_cells_score"], 0.05))
        keep = tad.obs["T_cells_score"] > thresh
        print(f"    keeping {int(keep.sum()):,} after T_cells_score > "
              f"{thresh:.2f} filter ({100*keep.mean():.1f}%)")
        tad = tad[keep].copy()
    if "n_genes_by_counts" in tad.obs.columns:
        ng_thresh = max(20, int(np.quantile(tad.obs["n_genes_by_counts"], 0.02)))
        keep = tad.obs["n_genes_by_counts"] >= ng_thresh
        print(f"    keeping {int(keep.sum()):,} with n_genes >= {ng_thresh}")
        tad = tad[keep].copy()
    panel = set(tad.var_names)
    if "lognorm" in tad.layers:
        tad.X = tad.layers["lognorm"]
    for label, genes in TCELL_SCORE_GENES.items():
        present = [g for g in genes if g in panel]
        if not present:
            continue
        sc.tl.score_genes(tad, gene_list=present,
                            score_name=f"score_{label}",
                            use_raw=False, random_state=0)
    # Recompute neighbors + UMAP using the existing scVI-immune latent
    # restricted to T-cells. Building a T-cell-only graph eliminates the
    # scatter caused by T-cells that originally embedded near non-T
    # neighbors in the full-immune UMAP.
    if "X_scvi_immune" in tad.obsm:
        rep = "X_scvi_immune"
    elif "X_scvi" in tad.obsm:
        rep = "X_scvi"
    else:
        sc.pp.pca(tad, n_comps=20, random_state=0)
        rep = "X_pca"
    print(f"    recomputing neighbors + UMAP on T-cells using {rep} (k=20, cosine)...")
    sc.pp.neighbors(tad, use_rep=rep, n_neighbors=20, metric="cosine",
                     random_state=0)
    sc.tl.umap(tad, min_dist=0.5, random_state=0)
    umap_key = "X_umap"
    score_cols = [f"score_{l}" for l in TCELL_SCORE_GENES.keys()
                    if f"score_{l}" in tad.obs.columns]
    fig, axes = plt.subplots(2, 3, figsize=(18, 11),
                                constrained_layout=True)
    axes = axes.ravel()
    # First panel: subtype with legend on side, larger font
    sc.pl.embedding(tad, basis=umap_key, color="immune_subtype",
                       ax=axes[0], show=False, size=10,
                       legend_loc="right margin", legend_fontsize=11,
                       title="T-cell subtype")
    for i, col in enumerate(score_cols):
        ax = axes[i+1]
        v = tad.obs[col].values
        vmax = float(np.quantile(v, 0.99))
        vmin = float(np.quantile(v, 0.01))
        sc.pl.embedding(tad, basis=umap_key, color=col,
                           ax=ax, show=False, size=10,
                           cmap="viridis", vmin=vmin, vmax=vmax,
                           title=col.replace("score_", ""))
    for j in range(len(score_cols)+1, 6):
        axes[j].axis("off")
    plt.suptitle("T-cell programs scored on the immune UMAP "
                  "(both samples combined)", fontsize=13)
    fig.savefig(OUT / "08_tcell_scores_umap.png",
                  dpi=200, bbox_inches="tight")
    plt.close()
    return tad   # return scored tad for re-use in fig09/10


# ===== 09: T-cell distance violin =====
def fig09_tcell_distance(adata):
    print("=== fig 09: T-cell distance violin ===")
    tmask = adata.obs["immune_subtype"].astype(str).str.startswith("T_")
    df = adata.obs[tmask][["immune_subtype", "dist_to_islet_um", "sample"]].copy()
    df["dist_to_islet_um"] = df["dist_to_islet_um"].astype(float)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5),
                                sharey=True, constrained_layout=True)
    order = sorted(df["immune_subtype"].unique())
    for ax, s in zip(axes, SAMPLES):
        sub = df[df["sample"].astype(str) == s]
        data = [np.log10(sub.loc[sub["immune_subtype"] == k,
                                    "dist_to_islet_um"].values + 1)
                for k in order]
        parts = ax.violinplot(data, showmeans=False, showmedians=True,
                                 widths=0.8)
        for pc in parts["bodies"]:
            pc.set_facecolor("#fdae61")
            pc.set_alpha(0.7)
            pc.set_edgecolor("black")
        ax.axhline(np.log10(ISLET_ZONE_UM+1), ls="--", color="red",
                     alpha=0.6, label=f"{int(ISLET_ZONE_UM)} μm (islet zone)")
        ax.axhline(np.log10(PROXIMAL_UM+1), ls=":", color="grey",
                     alpha=0.6, label=f"{int(PROXIMAL_UM)} μm (proximal)")
        ax.set_xticks(range(1, len(order)+1))
        ax.set_xticklabels(order, rotation=30, ha="right")
        ax.set_ylabel("log₁₀(distance to islet + 1) μm")
        ax.set_title(f"{s}  (n_T-cells = {len(sub):,})")
        ax.legend(loc="upper left", fontsize=8)
        ax.grid(True, alpha=0.3, axis="y")
    plt.suptitle("Distance to nearest islet by T-cell subtype",
                  fontsize=11)
    fig.savefig(OUT / "09_tcell_distance_violin.png",
                  dpi=200, bbox_inches="tight")
    plt.close()


# ===== 10: T-cell DE in-islet vs distal =====
def fig10_tcell_de(adata):
    print("=== fig 10: T-cell in-islet vs distal DE ===")
    tmask = adata.obs["immune_subtype"].astype(str).str.startswith("T_")
    tad = adata[tmask].copy()
    if "lognorm" in tad.layers:
        L = tad.layers["lognorm"]
    else:
        L = tad.X
    if hasattr(L, "toarray"):
        L = L.toarray()
    L = np.asarray(L)
    bins = tad.obs["distance_bin"].astype(str).values
    in_islet = bins == "intra_or_peri"
    distal = bins == "distal"
    print(f"    in_islet: {in_islet.sum():,}  distal: {distal.sum():,}")
    if in_islet.sum() < 30:
        print("    too few in-islet T-cells; skipping fig10")
        return
    # Restrict the DE to a curated immune/T-cell gene set so the top
    # hits reflect T-cell biology, not spatial contamination from
    # neighboring acinar (AMY1A/CUZD1) or endocrine (CHGA) cells.
    var_lookup = {g: i for i, g in enumerate(tad.var_names)}
    curated_idx = np.array([var_lookup[g] for g in TCELL_DE_GENES
                                if g in var_lookup])
    curated_names = [tad.var_names[i] for i in curated_idx]
    print(f"    DE restricted to {len(curated_idx)} curated T-cell/immune genes")

    subtypes = sorted(set(tad.obs["immune_subtype"].astype(str).unique()))
    fig, axes = plt.subplots(1, len(subtypes), figsize=(4*len(subtypes), 8.5),
                                sharey=False, constrained_layout=True)
    if len(subtypes) == 1:
        axes = [axes]
    for ax, st in zip(axes, subtypes):
        m_st = tad.obs["immune_subtype"].astype(str).values == st
        m_in = m_st & in_islet
        m_di = m_st & distal
        if m_in.sum() < 10 or m_di.sum() < 10:
            ax.text(0.5, 0.5, f"{st}\n(too few cells)",
                     transform=ax.transAxes, ha="center", va="center")
            ax.axis("off")
            continue
        # Slice to curated genes only
        L_in = L[m_in][:, curated_idx]
        L_di = L[m_di][:, curated_idx]
        mean_in = L_in.mean(axis=0)
        mean_di = L_di.mean(axis=0)
        delta = mean_in - mean_di
        # Drop genes with no expression at all in either bin
        keep = (mean_in + mean_di) > 0.02
        idx_local = np.where(keep)[0]
        delta_k = delta[idx_local]
        order_local = idx_local[np.argsort(delta_k)]
        top_up = order_local[-15:][::-1]
        top_dn = order_local[:5]
        plot_idx = np.concatenate([top_up, top_dn])
        plot_genes = [curated_names[i] for i in plot_idx]
        plot_delta = delta[plot_idx]
        colors = ["#d7191c" if x > 0 else "#2b83ba" for x in plot_delta]
        ax.barh(range(len(plot_genes)), plot_delta, color=colors,
                  edgecolor="black", linewidth=0.4)
        ax.set_yticks(range(len(plot_genes)))
        ax.set_yticklabels(plot_genes, fontsize=9)
        ax.invert_yaxis()
        ax.axvline(0, color="black", lw=0.5)
        ax.set_xlabel("Δ mean lognorm\n(in-islet − distal)")
        ax.set_title(f"{st}\nin={int(m_in.sum()):,} di={int(m_di.sum()):,}",
                       fontsize=10)
        ax.grid(True, alpha=0.3, axis="x")
    plt.suptitle("T-cell program shift at islet zone vs distal — curated "
                  f"immune-relevant genes only ({len(curated_idx)} genes)\n"
                  "red = up at islet, blue = down at islet",
                  fontsize=11)
    fig.savefig(OUT / "10_tcell_de_islet_vs_distal.png",
                  dpi=200, bbox_inches="tight")
    plt.close()


# ===== 11: Spatial overview =====
def fig11_spatial_overview(adata):
    print("=== fig 11: spatial overview ===")
    fig, axes = plt.subplots(1, 2, figsize=(20, 11),
                                constrained_layout=True)
    for ax, s in zip(axes, SAMPLES):
        seed_xy, seed_lab = get_seed_xy(s)
        m_imm = adata.obs["sample"].astype(str) == s
        imm_xy = adata.obsm["spatial"][m_imm]
        imm_subtype = adata.obs.loc[m_imm, "immune_subtype"].astype(str).values
        # Endocrine drawn FIRST and LARGER so islets are visible
        ax.scatter(seed_xy[:, 0], seed_xy[:, 1],
                     s=2.0, c="#ffb300", alpha=0.7, edgecolors="none",
                     label=f"endocrine seeds ({len(seed_xy):,})", zorder=2)
        family = np.array([
            "T-cell" if x.startswith("T_")
            else "B/Plasma" if x.startswith("B_")
            else "Macrophage" if x.startswith("Macro_")
            else "NK/DC/Mono" for x in imm_subtype])
        fam_color = {"T-cell": "#d7191c",
                      "B/Plasma": "#5e3c99",
                      "Macrophage": "#2b83ba",
                      "NK/DC/Mono": "#fdae61"}
        for fam, c in fam_color.items():
            mm = family == fam
            if mm.sum() == 0:
                continue
            ax.scatter(imm_xy[mm, 0], imm_xy[mm, 1],
                         s=0.3, c=c, alpha=0.45,
                         edgecolors="none",
                         label=f"{fam} ({int(mm.sum()):,})", zorder=1)
        # Overlay islet centroids as crosshairs
        n_islets = int(seed_lab.max() + 1) if (seed_lab >= 0).any() else 0
        if n_islets > 0:
            cents = np.array([seed_xy[seed_lab == k].mean(axis=0)
                                  for k in range(n_islets)])
            ax.scatter(cents[:, 0], cents[:, 1],
                         s=40, c="none", edgecolors="black",
                         linewidths=0.6, marker="o", zorder=3,
                         label=f"islet centroids ({n_islets})")
        ax.set_aspect("equal")
        ax.set_title(f"{s}", fontsize=14, fontweight="bold")
        ax.set_xticks([]); ax.set_yticks([])
        ax.legend(loc="upper right", framealpha=0.95, markerscale=4,
                    fontsize=10)
    plt.suptitle("Tissue-wide spatial distribution: endocrine seeds, "
                  "immune cells, and islet centroids", fontsize=13)
    fig.savefig(OUT / "11_spatial_overview.png",
                  dpi=200, bbox_inches="tight")
    plt.close()


# ===== 12: Per-islet immune composition =====
def fig12_islet_composition():
    print("=== fig 12: top-islet composition stacked bars ===")
    df = pd.read_csv(ROOT / "data/processed/islet_infiltration_per100endo.csv",
                      dtype={"sample": str})
    non_subtype = {"islet_id", "n_endocrine", "centroid_x", "centroid_y",
                     "sample", "total_immune", "insulitis_grade"}
    subtype_cols = [c for c in df.columns if c not in non_subtype]
    n_show = 12
    fig, axes = plt.subplots(2, 1, figsize=(14, 10),
                                constrained_layout=True)
    cmap = plt.get_cmap("tab20")
    for ax, s in zip(axes, SAMPLES):
        top = df[df["sample"] == s].nlargest(n_show, "total_immune").copy()
        if len(top) == 0 or top["total_immune"].max() == 0:
            ax.text(0.5, 0.5, f"{s}: no infiltrated islets",
                     transform=ax.transAxes, ha="center")
            continue
        # Normalize to fractions for stacking
        ids = top["islet_id"].str.split("_").str[-1].values
        bottom = np.zeros(len(top))
        x = np.arange(len(top))
        for i, st in enumerate(subtype_cols):
            vals = top[st].values
            ax.bar(x, vals, bottom=bottom, label=st, color=cmap(i % 20),
                     edgecolor="black", linewidth=0.3)
            bottom += vals
        ax.set_xticks(x)
        ax.set_xticklabels([f"{i}\n(endo={int(top.iloc[k]['n_endocrine'])})"
                              for k, i in enumerate(ids)], fontsize=9)
        ax.set_ylabel("Immune cells within 50 μm")
        ax.set_title(f"{s} — top {n_show} most-infiltrated islets")
        ax.legend(loc="upper right", bbox_to_anchor=(1.18, 1.0),
                    fontsize=8, framealpha=0.95)
        ax.grid(True, alpha=0.3, axis="y")
    plt.suptitle("Per-islet immune subtype composition (top infiltrated)",
                  fontsize=11)
    fig.savefig(OUT / "12_islet_immune_composition.png",
                  dpi=200, bbox_inches="tight")
    plt.close()


# ===== Run =====
if __name__ == "__main__":
    t0 = time.time()
    recompute_distances()
    print(f"\n=== loading concatenated immune adata ===")
    a = load_immune()
    print(f"    shape: {a.shape}")
    print(f"    samples: {a.obs['sample'].astype(str).value_counts().to_dict()}")

    fig01_insulitis_grades()
    fig02_density_enrichment()
    fig03_distance_cdf(a)
    fig04_islet_size_immune()
    fig05_top_insulitis(a)
    fig06_immune_marker_dotplot(a)
    fig07_tcell_marker_dotplot(a)
    fig08_tcell_scores(a)
    fig09_tcell_distance(a)
    fig10_tcell_de(a)
    fig11_spatial_overview(a)
    fig12_islet_composition()

    print(f"\n=== DONE in {time.time()-t0:.0f}s ===")
    print(f"Figures saved to: {OUT}")
    for p in sorted(OUT.glob("*.png")):
        sz = p.stat().st_size / 1024
        print(f"  {p.name}  ({sz:.0f} KB)")
