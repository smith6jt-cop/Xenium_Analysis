"""End-to-end re-analysis of immune cells from the gated subset.

Reads:
  data/processed/immune_combined_gated.h5ad   ← from regate_immune.py

Does:
  1. Subtype scoring (score_genes) for all 12 fine subtypes restricted to
     biologically appropriate parent gates:
        T_cells parent  → T_helper, T_cytotoxic, T_reg, T_exhausted
        NK parent       → NK
        Myeloid parent  → Macro_M1, Macro_M2, Macro_resident, Monocyte
        DC parent       → DC
        B_cells parent  → B_pan, B_plasma
  2. Retrains scVI on the gated immune subset with batch_key='sample'
     (clean sample-batch correction over true immune cells only).
  3. Recomputes neighbors + UMAP on the new scVI embedding.
  4. Recomputes per-immune-cell distance to nearest endo-label seed.
  5. Computes insulitis grades + density enrichment + per-islet
     infiltration tables, OVERWRITING the headline CSVs:
        immune_proximity_summary.csv      (density enrichment)
        islet_infiltration_per100endo.csv (per-islet counts)
        islet_insulitis_grades.csv        (per-sample summary)
  6. Saves a final immune h5ad with all subtype scores + UMAP +
     distance columns:
        data/processed/immune_combined_phenotyped.h5ad

Outputs are sample-comparable: no static raw-count thresholds anywhere.
"""
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path("/blue/maigan/smith6jt/Xenium_Analysis")
SAMPLES = ("0041323", "0041326")
ENDO_LABELS = {"Endocrine", "Beta", "Alpha", "Delta", "Endocrine_pan"}
EPS_UM = 50.0
MIN_SAMPLES = 10
ISLET_ZONE_UM = 50.0
PROXIMAL_UM = 200.0
INSULITIS_THRESH = 6
PERI_THRESH = 3

# ---- Subtype panels (restricted to parent-gate cells) ----
SUBTYPE_PANELS = {
    "T_cells": {
        "T_helper":     ["CD4", "IL7R", "CCR7", "CD2"],
        "T_cytotoxic":  ["CD8A", "CD8B", "GZMA", "GZMB", "GZMK", "GZMH",
                            "PRF1", "NKG7", "IFNG"],
        "T_reg":        ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "ENTPD1"],
        "T_exhausted":  ["PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX", "EOMES"],
    },
    "NK": {
        "NK":           ["NCAM1", "KLRD1", "KLRF1", "NCR1", "NCR3", "FCGR3A"],
    },
    "Myeloid": {
        "Macro_M1":     ["TNF", "IL6", "IL1B", "CXCL10", "CXCL9", "CD86"],
        "Macro_M2":     ["MRC1", "CD163", "CCL18", "CCL22", "VSIG4"],
        "Macro_resident":["CD68", "MARCO", "VSIG4", "CSF1R"],
        "Monocyte":     ["CD14", "FCGR3A", "S100A8", "S100A9", "VCAN"],
    },
    "DC": {
        "DC":           ["CD1C", "CLEC9A", "CD1A", "LAMP3", "CLEC10A",
                            "FCER1A"],
    },
    "B_cells": {
        "B_pan":        ["MS4A1", "CD79A", "CD79B", "CD19", "CD22"],
        "B_plasma":     ["XBP1", "PRDM1", "MZB1", "CD27"],
    },
}


def get_seed_xy(sample):
    """Return endocrine clustered seed coords + cluster labels."""
    p = ROOT / f"data/processed/{sample}/{sample}_phenotyped.h5ad"
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


def main():
    t0 = time.time()
    print("=" * 72)
    print("STEP 1: Load gated immune adata")
    print("=" * 72, flush=True)
    p = ROOT / "data/processed/immune_combined_gated.h5ad"
    a = sc.read_h5ad(p)
    print(f"  shape: {a.shape}", flush=True)
    print(f"  per-sample: {a.obs['sample'].astype(str).value_counts().to_dict()}",
          flush=True)
    print(f"  per-gate: {a.obs['celltype_gated'].value_counts().to_dict()}",
          flush=True)
    if "lognorm" in a.layers:
        a.X = a.layers["lognorm"]
    panel = set(a.var_names)

    # ---- STEP 2: score subtypes within each parent gate ----
    print(f"\n{'='*72}\nSTEP 2: Score immune subtypes within parent gates\n{'='*72}",
          flush=True)
    a.obs["immune_subtype"] = pd.Categorical(a.obs["celltype_gated"].astype(str))
    for parent, subtypes in SUBTYPE_PANELS.items():
        parent_mask = a.obs["celltype_gated"].astype(str).values == parent
        n_parent = int(parent_mask.sum())
        if n_parent == 0:
            continue
        print(f"  {parent}  (n={n_parent:,}):", flush=True)
        sub = a[parent_mask].copy()
        for st_name, genes in subtypes.items():
            present = [g for g in genes if g in panel]
            if not present:
                continue
            sc.tl.score_genes(sub, gene_list=present,
                                score_name=f"{st_name}_score",
                                use_raw=False, random_state=0)
        # Argmax-assignment within parent
        score_cols = [f"{k}_score" for k in subtypes if f"{k}_score" in sub.obs.columns]
        if not score_cols:
            print(f"    no score_cols, keeping {parent} label")
            continue
        scores = sub.obs[score_cols].values
        best_idx = scores.argmax(axis=1)
        best_st = [score_cols[i].replace("_score", "") for i in best_idx]
        # Apply back to a.obs
        new_immune_subtype = a.obs["immune_subtype"].astype(str).copy().values
        new_immune_subtype[parent_mask] = best_st
        a.obs["immune_subtype"] = pd.Categorical(new_immune_subtype)
        # Show distribution
        vc = pd.Series(best_st).value_counts()
        print(f"    {vc.to_dict()}", flush=True)
        # Copy individual score columns back too
        for col in score_cols:
            a.obs.loc[parent_mask, col] = sub.obs[col].values
            a.obs[col] = a.obs[col].fillna(0.0)

    # ---- STEP 3: scVI retraining on the cleaned subset ----
    print(f"\n{'='*72}\nSTEP 3: Retrain scVI on gated immune (batch=sample)\n{'='*72}",
          flush=True)
    try:
        import scvi
        import torch
        torch.set_float32_matmul_precision("high")
        torch.backends.cuda.matmul.allow_tf32 = True
        a.layers["counts"] = a.layers.get("counts", a.X)
        scvi.model.SCVI.setup_anndata(
            a, layer="counts", batch_key="sample",
            continuous_covariate_keys=["cell_area"]
            if "cell_area" in a.obs.columns else None)
        model = scvi.model.SCVI(a, n_latent=20, gene_likelihood="nb")
        print(f"  training scVI...", flush=True)
        t = time.time()
        model.train(max_epochs=100, batch_size=512, early_stopping=True,
                     accelerator="gpu", devices=1)
        print(f"    {time.time()-t:.0f}s", flush=True)
        a.obsm["X_scvi_immune"] = model.get_latent_representation()
        rep = "X_scvi_immune"
    except Exception as e:
        print(f"  scVI failed ({e}); falling back to PCA on lognorm",
              flush=True)
        sc.pp.pca(a, n_comps=20, random_state=0)
        rep = "X_pca"

    # ---- STEP 4: neighbors + UMAP on cleaned embedding ----
    print(f"\n{'='*72}\nSTEP 4: neighbors + UMAP\n{'='*72}", flush=True)
    sc.pp.neighbors(a, use_rep=rep, n_neighbors=15, metric="cosine",
                      random_state=0)
    sc.tl.umap(a, min_dist=0.5, random_state=0)
    sc.tl.leiden(a, resolution=0.8, flavor="igraph", n_iterations=-1,
                   random_state=0, key_added="leiden_immune")
    print(f"  leiden_immune distribution: {a.obs['leiden_immune'].value_counts().head(10).to_dict()}",
          flush=True)

    # ---- STEP 5: distance to islet + insulitis grades + tables ----
    print(f"\n{'='*72}\nSTEP 5: distance + insulitis + tables\n{'='*72}",
          flush=True)
    from sklearn.neighbors import NearestNeighbors
    a.obs["dist_to_islet_um"] = np.nan
    a.obs["nearest_islet_id"] = ""
    all_per_islet = []
    all_enrichment = []
    grade_rows = []
    for s in SAMPLES:
        m = a.obs["sample"].astype(str) == s
        if m.sum() == 0:
            continue
        seed_xy, seed_lab = get_seed_xy(s)
        n_islets = int(seed_lab.max() + 1) if (seed_lab >= 0).any() else 0
        print(f"  {s}: {n_islets} islets, {m.sum():,} immune cells", flush=True)

        # NN
        imm_xy = np.asarray(a[m].obsm["spatial"])
        nn = NearestNeighbors(n_neighbors=1, n_jobs=-1).fit(seed_xy)
        d, idx = nn.kneighbors(imm_xy); d = d.ravel(); idx = idx.ravel()
        nearest = seed_lab[idx]
        a.obs.loc[m, "dist_to_islet_um"] = d
        a.obs.loc[m, "nearest_islet_id"] = [
            f"{s}_islet_{k:04d}" for k in nearest]

        in_zone = d < ISLET_ZONE_UM
        in_proximal = (d >= ISLET_ZONE_UM) & (d < PROXIMAL_UM)
        in_distal = d >= PROXIMAL_UM

        immune_subtype = a.obs.loc[m, "immune_subtype"].astype(str).values
        all_subtypes = sorted(set(immune_subtype))

        # Per-islet counts + insulitis grade
        df_imm = pd.DataFrame({
            "subtype": immune_subtype[in_zone],
            "islet_id": [f"{s}_islet_{k:04d}" for k in nearest[in_zone]],
        })
        per_islet = (df_imm.groupby(["islet_id", "subtype"], observed=True)
                      .size().unstack("subtype", fill_value=0))
        for st in all_subtypes:
            if st not in per_islet.columns:
                per_islet[st] = 0
        per_islet = per_islet[all_subtypes].reset_index()
        islet_meta = pd.DataFrame({
            "islet_id": [f"{s}_islet_{k:04d}" for k in range(n_islets)],
            "n_endocrine": [int((seed_lab == k).sum()) for k in range(n_islets)],
            "centroid_x": [seed_xy[seed_lab == k, 0].mean() for k in range(n_islets)],
            "centroid_y": [seed_xy[seed_lab == k, 1].mean() for k in range(n_islets)],
        })
        per_islet = islet_meta.merge(per_islet, on="islet_id", how="left").fillna(0)
        per_islet[all_subtypes] = per_islet[all_subtypes].astype(int)
        per_islet["sample"] = s
        per_islet["total_immune"] = per_islet[all_subtypes].sum(axis=1)
        per_islet["insulitis_grade"] = per_islet["total_immune"].apply(
            lambda n: "insulitis" if n >= INSULITIS_THRESH
            else ("peri_insulitis" if n >= PERI_THRESH else "no_insulitis"))
        all_per_islet.append(per_islet)

        # Density enrichment
        n_zone_total = int(in_zone.sum())
        n_distal_total = int(in_distal.sum())
        for st in all_subtypes:
            is_st = immune_subtype == st
            n_zone_st = int((is_st & in_zone).sum())
            n_distal_st = int((is_st & in_distal).sum())
            zone_frac = n_zone_st / max(n_zone_total, 1)
            distal_frac = n_distal_st / max(n_distal_total, 1)
            ratio = (zone_frac / distal_frac) if distal_frac > 0 else float("nan")
            all_enrichment.append({
                "sample": s,
                "immune_subtype": st,
                "n_in_islet_zone": n_zone_st,
                "n_in_proximal": int((is_st & in_proximal).sum()),
                "n_in_distal": n_distal_st,
                "zone_frac": zone_frac,
                "distal_frac": distal_frac,
                "density_enrichment": ratio,
                "n_islet_zone_total": n_zone_total,
                "n_distal_total": n_distal_total,
            })
        # Grade summary
        gc = per_islet["insulitis_grade"].value_counts().to_dict()
        for k in ("no_insulitis", "peri_insulitis", "insulitis"):
            gc.setdefault(k, 0)
        grade_rows.append({"sample": s,
                            "no_insulitis": gc["no_insulitis"],
                            "peri_insulitis": gc["peri_insulitis"],
                            "insulitis": gc["insulitis"],
                            "total_islets": int(n_islets)})

    # Distance bin column
    bins = np.full(a.n_obs, "distal", dtype=object)
    d_all = a.obs["dist_to_islet_um"].fillna(np.inf).values
    bins[d_all < PROXIMAL_UM] = "proximal"
    bins[d_all < ISLET_ZONE_UM] = "intra_or_peri"
    a.obs["distance_bin"] = pd.Categorical(
        bins, categories=["intra_or_peri", "proximal", "distal"])

    per_df = pd.concat(all_per_islet, ignore_index=True)
    enr_df = pd.DataFrame(all_enrichment)
    grade_df = pd.DataFrame(grade_rows)
    grade_df["pct_insulitis"] = 100 * grade_df["insulitis"] / grade_df["total_islets"]
    grade_df["pct_peri_or_insulitis"] = 100 * (grade_df["insulitis"]
                                                  + grade_df["peri_insulitis"]
                                                  ) / grade_df["total_islets"]

    # ---- Headline CSVs are now owned by scripts/insulitis_analysis.py ----
    # The per-phenotype dynamic-grading rewrite (2026-04-28) moved the three
    # headline tables (immune_proximity_summary.csv,
    # islet_infiltration_per100endo.csv, islet_insulitis_grades.csv) to a
    # single writer. Running this block would silently revert the schema.
    print(f"\n{'='*72}\nSTEP 6: skipped — headline CSV writes redirected\n{'='*72}",
            flush=True)
    print("  Headline tables (immune_proximity_summary.csv, "
            "islet_infiltration_per100endo.csv, islet_insulitis_grades.csv)\n"
            "  are now owned by scripts/insulitis_analysis.py "
            "(per-phenotype dynamic grading).", flush=True)
    print("  Run: python scripts/insulitis_analysis.py for headline tables", flush=True)
    # Continue to immune h5ad saves — those are still owned by this script.

    # ---- Save final immune h5ad ----
    p_out = ROOT / "data/processed/immune_combined_phenotyped.h5ad"
    print(f"\n  writing {p_out.name}...", flush=True)
    a.write(p_out, compression="gzip")
    print(f"    done ({p_out.stat().st_size/1e6:.0f} MB)", flush=True)

    # ---- Per-sample immune h5ads (small, for figs script) ----
    for s in SAMPLES:
        m = a.obs["sample"].astype(str) == s
        if m.sum() == 0:
            continue
        sub = a[m].copy()
        p_s = ROOT / f"data/processed/{s}/{s}_immune_phenotyped.h5ad"
        sub.write(p_s, compression="gzip")
        print(f"    wrote {p_s.name} ({p_s.stat().st_size/1e6:.0f} MB)",
              flush=True)

    print(f"\n=== DONE in {time.time()-t0:.0f}s ===", flush=True)


if __name__ == "__main__":
    main()
