"""Sample-comparable immune cell identification via marker gating.

The stage 02 score_genes-argmax phenotyping is unbalanced for immune
cell types: 0041323 had 0 Myeloid (false negatives), and 0041326 had
127k Myeloid of which 80% failed a basic marker gate.

This script applies sample-comparable rule-based gating using the
sample-agnostic lognorm layer:

  Myeloid:  >= 1 of [CD68, CD163, MARCO, CD14, CSF1R, VSIG4, FCGR1A]
            (ITGAM/ITGAX/FCGR3A dropped — also expressed by NK)
  T_cells:  >= 1 of [CD3E, CD3G, CD2, CD8A, CD8B, CD4, CD7, LCK, ZAP70]
  B_cells:  >= 1 of [MS4A1, CD79A, CD79B, CD19, CD22]
  NK:       >= 1 of [NCAM1, NCR1, NCR3, KLRD1, KLRF1, KLRB1] AND
                  not passing T_cells gate (NK is CD3-)
  DC:       >= 1 of [CD1C, CD1A, CLEC9A, LAMP3, CLEC10A, FCER1A]

All gates require lognorm >= log1p(3) ≈ 1.39 (sample-comparable).

When a cell passes multiple gates (rare), the gate with the highest
mean lognorm of present markers wins.

Outputs (does NOT modify the giant phenotyped h5ad files):
  data/processed/{sample}/{sample}_celltype_gated.csv
      cell_id, celltype_gated, gate_pass_*  (one row per cell)
  data/processed/{sample}/{sample}_immune_gated.h5ad
      gate-passing cells only, ready for downstream analysis
  data/processed/immune_combined_gated.h5ad
      both samples concatenated
"""
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path("/blue/maigan/smith6jt/Xenium_Analysis")
SAMPLES = ("0041323", "0041326")

LOG_DETECT = float(np.log1p(3.0))  # ≈ 1.386

# Order is precedence for multi-gate cells: more-specific gates first.
GATES = {
    "B_cells": ["MS4A1", "CD79A", "CD79B", "CD19", "CD22"],
    "T_cells": ["CD3E", "CD3G", "CD2", "CD8A", "CD8B", "CD4", "CD7",
                  "LCK", "ZAP70"],
    "NK":      ["NCAM1", "NCR1", "NCR3", "KLRD1", "KLRF1", "KLRB1"],
    "DC":      ["CD1C", "CD1A", "CLEC9A", "LAMP3", "CLEC10A", "FCER1A"],
    "Myeloid": ["CD68", "CD163", "MARCO", "CD14", "CSF1R", "VSIG4",
                  "FCGR1A"],
}


def gate_sample(sample):
    print(f"\n{'='*72}\n=== {sample}\n{'='*72}", flush=True)
    p = ROOT / f"data/processed/{sample}/{sample}_phenotyped.h5ad"
    print(f"  loading {p.name}...", flush=True)
    t = time.time()
    a = sc.read_h5ad(p)
    print(f"    shape={a.shape}  ({time.time()-t:.1f}s)", flush=True)
    panel = set(a.var_names)
    L = a.layers["lognorm"]
    if hasattr(L, "tocsr"):
        L = L.tocsr()

    n = a.n_obs
    gate_pass = {}
    gate_score = {}
    for name, genes in GATES.items():
        present = [g for g in genes if g in panel]
        idx = [list(a.var_names).index(g) for g in present]
        sub = L[:, idx]
        if hasattr(sub, "toarray"):
            sub = sub.toarray()
        passes = (sub >= LOG_DETECT).any(axis=1)
        score = np.where(sub >= LOG_DETECT, sub, 0.0).mean(axis=1)
        gate_pass[name] = passes
        gate_score[name] = score
        print(f"  {name:8s} gate: {int(passes.sum()):>7,} cells pass "
              f"({100*passes.mean():5.2f}%)  markers={present}", flush=True)

    # NK is CD3- by definition: enforce
    if "T_cells" in gate_pass and "NK" in gate_pass:
        nk_minus_t = int((gate_pass["NK"] & gate_pass["T_cells"]).sum())
        gate_pass["NK"] = gate_pass["NK"] & ~gate_pass["T_cells"]
        print(f"  NK ∩ T_cells overlap (NKT-like): {nk_minus_t:,}  "
              f"→ assigned to T_cells (CD3+)", flush=True)

    # Resolve multi-pass: assign to gate with highest score
    any_pass = np.zeros(n, dtype=bool)
    for name in GATES:
        any_pass |= gate_pass[name]
    print(f"  any-immune-gate-pass: {int(any_pass.sum()):,} "
          f"({100*any_pass.mean():.2f}%)", flush=True)

    orig_ct = a.obs["celltype"].astype(str).values
    gated = orig_ct.copy()
    # Demote false positives
    immune_orig_set = set(["Myeloid", "T_cells", "B_cells"])
    orig_immune = np.array([c in immune_orig_set for c in orig_ct])
    false_pos = orig_immune & ~any_pass
    print(f"  false positives demoted (orig immune label, no gate): "
          f"{int(false_pos.sum()):,}", flush=True)
    gated[false_pos] = "Indeterminate"

    if any_pass.sum() > 0:
        score_stack = np.stack([gate_score[k] for k in GATES.keys()], axis=1)
        best = np.argmax(score_stack, axis=1)
        gate_keys = list(GATES.keys())
        for i_gate, k in enumerate(gate_keys):
            mask = any_pass & (best == i_gate) & gate_pass[k]
            gated[mask] = k
            print(f"    assigned {int(mask.sum()):>7,} cells to {k}",
                  flush=True)

    # Sanity: total assigned
    counts = pd.Series(gated).value_counts()
    print(f"  gated celltype distribution:\n{counts.head(15).to_string()}",
          flush=True)

    # === Write a small CSV of gating columns, NEVER overwrite phenotyped ===
    df_csv = pd.DataFrame({
        "cell_id": a.obs.index.values,
        "celltype": orig_ct,
        "celltype_gated": gated,
        "gate_pass_myeloid": gate_pass["Myeloid"],
        "gate_pass_t_cells": gate_pass["T_cells"],
        "gate_pass_b_cells": gate_pass["B_cells"],
        "gate_pass_nk":      gate_pass["NK"],
        "gate_pass_dc":      gate_pass["DC"],
        "gate_score_myeloid": gate_score["Myeloid"],
        "gate_score_t_cells": gate_score["T_cells"],
        "gate_score_b_cells": gate_score["B_cells"],
        "gate_score_nk":      gate_score["NK"],
        "gate_score_dc":      gate_score["DC"],
    })
    p_csv = ROOT / f"data/processed/{sample}/{sample}_celltype_gated.csv"
    df_csv.to_csv(p_csv, index=False)
    print(f"  wrote {p_csv.name}  ({p_csv.stat().st_size/1e6:.1f} MB)",
          flush=True)

    # === Immune-only h5ad ===
    immune_mask = any_pass
    print(f"  total gate-passing immune: {int(immune_mask.sum()):,}",
          flush=True)
    imm = a[immune_mask].copy()
    imm.obs["celltype_gated"] = pd.Categorical(gated[immune_mask])
    imm.obs["sample"] = sample
    p_imm = ROOT / f"data/processed/{sample}/{sample}_immune_gated.h5ad"
    print(f"  writing {p_imm.name}...", flush=True)
    imm.write(p_imm, compression="gzip")
    print(f"    done ({p_imm.stat().st_size/1e6:.0f} MB)", flush=True)
    return imm


if __name__ == "__main__":
    immune_parts = []
    for s in SAMPLES:
        imm = gate_sample(s)
        immune_parts.append(imm)

    print(f"\n{'='*72}\n=== Concatenating both samples\n{'='*72}",
          flush=True)
    combined = ad.concat(immune_parts, axis=0, join="inner",
                          index_unique="-", merge="first")
    print(f"  combined immune shape: {combined.shape}", flush=True)
    for s in SAMPLES:
        n = (combined.obs["sample"].astype(str) == s).sum()
        print(f"    {s}: {n:,} cells", flush=True)
    n0, n1 = ((combined.obs["sample"] == s).sum() for s in SAMPLES)
    print(f"  ratio: {n1/max(n0,1):.2f}x  (was 9.21x with score_genes-only)",
          flush=True)

    print(f"\n  composition by gate:", flush=True)
    for k in GATES:
        n = (combined.obs["celltype_gated"].astype(str) == k).sum()
        print(f"    {k}: {n:,}", flush=True)

    p_out = ROOT / "data/processed/immune_combined_gated.h5ad"
    print(f"  writing {p_out.name}...", flush=True)
    combined.write(p_out, compression="gzip")
    print(f"    done ({p_out.stat().st_size/1e6:.0f} MB)", flush=True)
    print(f"\n=== DONE ===", flush=True)
