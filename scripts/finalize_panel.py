"""Stage 0 — finalize panel masks + donor_sex tag on each phenotyped/raw h5ad.

Reads:
    data/processed/panel_audit.csv        — 106 candidate exclusions
    data/processed/{sample}/{sample}_phenotyped.h5ad  (or _preprocessed.h5ad)

Writes back, in-place:
    adata.var['panel_for_scoring']    bool — True for genes safe in cell-type
                                       scoring (drops 100 non-sex tissue-restricted)
    adata.var['panel_for_embedding']  bool — additionally drops 6 sex-chromosome
                                       genes; for HVG/PCA/scVI input
    adata.var['exclude_category']     str — audit category (or '' for kept genes)
    adata.obs['donor_sex']            str — 'male' / 'female' from KDM5D+DDX3Y

This script is idempotent: re-running it just overwrites the flags. It writes
to whichever h5ad file is passed; the same logic applies to raw zarr ingestion
output and to phenotyped outputs.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path("/blue/maigan/smith6jt/Xenium_Analysis")
SAMPLES = ("0041323", "0041326")
AUDIT_CSV = ROOT / "data/processed/panel_audit.csv"

# Sex-call thresholds (lognorm scale)
# KDM5D + DDX3Y are reliably 7-10% in male, near-zero in female.
SEX_GENE_THRESH_LOGNORM = np.log1p(1.0)
SEX_DETECTION_FRAC = 0.02  # ≥ 2% of cells expressing → male


def annotate_var_flags(adata, audit_df):
    """Add panel_for_scoring, panel_for_embedding, exclude_category to var."""
    audit_lookup = audit_df.set_index("gene")["category"].to_dict()
    cat_col = adata.var_names.map(lambda g: audit_lookup.get(g, ""))
    adata.var["exclude_category"] = pd.Categorical(cat_col)
    sex_chrom_genes = set(audit_df[audit_df["category"] == "sex_chromosome"]["gene"])
    nonsex_excluded = set(audit_df[audit_df["exclude_recommended"] == "yes"]["gene"]) - sex_chrom_genes

    adata.var["panel_for_scoring"] = ~adata.var_names.isin(nonsex_excluded)
    adata.var["panel_for_embedding"] = (~adata.var_names.isin(nonsex_excluded)
                                          & ~adata.var_names.isin(sex_chrom_genes))


def call_donor_sex(adata):
    """Return 'male' or 'female' from KDM5D + DDX3Y detection."""
    layer = adata.layers["lognorm"] if "lognorm" in adata.layers else adata.X
    sex_genes = [g for g in ("KDM5D", "DDX3Y") if g in adata.var_names]
    if not sex_genes:
        return "unknown"
    n_cells = adata.n_obs
    n_detected_max = 0
    for g in sex_genes:
        i = adata.var_names.get_loc(g)
        col = layer[:, i]
        if hasattr(col, "toarray"):
            col = col.toarray().ravel()
        col = np.asarray(col).ravel()
        n_det = int((col >= SEX_GENE_THRESH_LOGNORM).sum())
        n_detected_max = max(n_detected_max, n_det)
    frac = n_detected_max / max(n_cells, 1)
    return "male" if frac >= SEX_DETECTION_FRAC else "female"


def process_one(p, audit_df, dry_run=False):
    print(f"\n=== {p.name} ===")
    a = sc.read_h5ad(p)
    print(f"    shape: {a.shape}")
    annotate_var_flags(a, audit_df)
    n_score = int(a.var["panel_for_scoring"].sum())
    n_embed = int(a.var["panel_for_embedding"].sum())
    print(f"    panel_for_scoring: {n_score:,} / {a.n_vars:,}")
    print(f"    panel_for_embedding: {n_embed:,} / {a.n_vars:,}")
    sex = call_donor_sex(a)
    a.obs["donor_sex"] = pd.Categorical([sex] * a.n_obs,
                                          categories=["female", "male", "unknown"])
    print(f"    donor_sex: {sex}")
    if dry_run:
        print("    DRY RUN — not saving")
        return
    a.write(p, compression="gzip")
    print(f"    saved")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paths", nargs="+", default=None,
                       help="h5ad paths to annotate (default: both samples' "
                            "_phenotyped.h5ad)")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    audit_df = pd.read_csv(AUDIT_CSV)
    print(f"loaded audit: {len(audit_df):,} genes, "
          f"{(audit_df['exclude_recommended']=='yes').sum()} flagged")

    if args.paths:
        paths = [Path(p) for p in args.paths]
    else:
        paths = [ROOT / f"data/processed/{s}/{s}_phenotyped.h5ad" for s in SAMPLES]

    for p in paths:
        if not p.exists():
            print(f"\n=== {p.name}: NOT FOUND, skipping")
            continue
        process_one(p, audit_df, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
