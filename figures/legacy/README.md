# Legacy figures — DO NOT USE

These figures are from prior experimental pipeline branches that are no
longer part of the current workflow:

- `integrated_0041323_0041326/` — outputs from `notebooks/05_tissue_comparisons.ipynb`
  cross-sample integration. The current per-sample workflow does NOT do
  integration (decided 2026-04-24). The notebook still exists in
  `notebooks/` for reference but is not part of the default driver.

- `notebooks_figures_orphan/` — same integration figures, written
  with a relative path that landed under `notebooks/figures/`. Moved
  here so the working tree stays clean.

- `combined_07_immune_pre_lognorm_fix/` — combined-sample immune UMAP
  + per-subtype distance violin produced by `notebooks/07_immune_islet_proximity.ipynb`
  on 2026-04-27 12:12, **before** the lognorm sample-comparable
  threshold fix landed in the immune phenotyping (h5ads regenerated
  2026-04-27 14:03). Replaced 2026-04-28 by `figures/07_immune_islet/13_combined_immune_umap.png`
  and `14_all_subtype_distance_violin.png` rendered from the
  current `_immune_phenotyped.h5ad` files.

If you need integration in the future, regenerate by running stage 05
explicitly: `STAGES="05" bash scripts/run_local_pipeline.sh`.
