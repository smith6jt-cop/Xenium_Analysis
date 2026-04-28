"""Render the 12-figure immune+islet suite from the cleaned, gated data.
Skips the dist-to-islet recompute (already done by rerun_immune_pipeline.py).
"""
import sys
sys.path.insert(0, "scripts")
from figs_07_immune import (
    fig01_insulitis_grades, fig02_density_enrichment, fig03_distance_cdf,
    fig04_islet_size_immune, fig05_top_insulitis,
    fig06_immune_marker_dotplot, fig07_tcell_marker_dotplot,
    fig08_tcell_scores, fig09_tcell_distance, fig10_tcell_de,
    fig11_spatial_overview, fig12_islet_composition,
    load_immune,
)
import time

t0 = time.time()
print("=== loading concatenated immune adata ===")
a = load_immune()
print(f"    shape: {a.shape}")
print(f"    samples: {a.obs['sample'].astype(str).value_counts().to_dict()}")
print(f"    subtypes: {a.obs['immune_subtype'].value_counts().to_dict()}")

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
