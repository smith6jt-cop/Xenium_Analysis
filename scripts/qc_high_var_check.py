"""Identify genes with normalized variance > 10 in 0041323_preprocessed."""
import h5py
import numpy as np

with h5py.File("data/processed/0041323/0041323_preprocessed.h5ad", "r") as f:
    var_norm = f["var/variances_norm"][:]
    means = f["var/means"][:]
    var_raw = f["var/variances"][:]
    is_hv = f["var/highly_variable"][:].astype(bool)
    rank = f["var/highly_variable_rank"][:]
    grp = f["var"]
    if "_index" in grp:
        idx_name = grp["_index"].asstr()[:]
    else:
        idx_attr = grp.attrs.get("_index", "_index")
        idx_name = grp[idx_attr].asstr()[:]
    n_cells = f["var/n_cells_by_counts"][:] if "n_cells_by_counts" in grp else None
    pct_drop = f["var/pct_dropout_by_counts"][:] if "pct_dropout_by_counts" in grp else None

mask = var_norm > 10
idxs = np.where(mask)[0]
order = idxs[np.argsort(-var_norm[idxs])]

print(f"Genes with variances_norm > 10: {int(mask.sum())}")
print()
header = f"{'gene':<22}{'var_norm':>10}{'mean':>10}{'var_raw':>12}{'n_cells':>10}{'pct_drop':>10}{'HV':>5}{'rank':>8}"
print(header)
print("-" * len(header))
for i in order:
    nc = int(n_cells[i]) if n_cells is not None else -1
    pd_ = float(pct_drop[i]) if pct_drop is not None else -1.0
    rk = int(rank[i]) if not np.isnan(rank[i]) else -1
    print(
        f"{idx_name[i]:<22}{var_norm[i]:>10.3f}{means[i]:>10.4f}"
        f"{var_raw[i]:>12.4f}{nc:>10}{pd_:>10.2f}{int(is_hv[i]):>5}{rk:>8}"
    )

print()
print("Top 20 by RAW variance (for context):")
print(header.replace("HV", "  ").replace("rank", "    "))
print("-" * len(header))
top_raw = np.argsort(-var_raw)[:20]
for i in top_raw:
    nc = int(n_cells[i]) if n_cells is not None else -1
    pd_ = float(pct_drop[i]) if pct_drop is not None else -1.0
    rk = int(rank[i]) if not np.isnan(rank[i]) else -1
    print(
        f"{idx_name[i]:<22}{var_norm[i]:>10.3f}{means[i]:>10.4f}"
        f"{var_raw[i]:>12.4f}{nc:>10}{pd_:>10.2f}{int(is_hv[i]):>5}{rk:>8}"
    )

print()
print("Distribution of variances_norm:")
for q in (0.5, 0.75, 0.9, 0.95, 0.98, 0.99, 0.995, 1.0):
    print(f"  q{q:.3f}: {np.quantile(var_norm, q):.3f}")
print(f"  count > 1:  {int((var_norm > 1).sum())}")
print(f"  count > 5:  {int((var_norm > 5).sum())}")
print(f"  count > 10: {int((var_norm > 10).sum())}")
