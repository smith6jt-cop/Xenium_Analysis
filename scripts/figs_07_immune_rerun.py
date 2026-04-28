"""Re-render only the figures that had layout/legend issues.
Skips dist_to_islet recompute (already done) and the panels that were OK."""
import sys
sys.path.insert(0, "scripts")
from figs_07_immune import (
    fig05_top_insulitis, fig06_immune_marker_dotplot,
    fig07_tcell_marker_dotplot, fig08_tcell_scores,
    fig11_spatial_overview, load_immune,
)

print("=== loading concatenated immune adata ===")
a = load_immune()
print(f"    shape: {a.shape}")

fig05_top_insulitis(a)
fig06_immune_marker_dotplot(a)
fig07_tcell_marker_dotplot(a)
fig08_tcell_scores(a)
fig11_spatial_overview(a)
print("=== DONE ===")
