"""Probe how sc.pl.dotplot / sc.pl.matrixplot lay out their axes."""
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt

# Small synthetic dataset
adata = sc.datasets.pbmc3k_processed()
adata = adata[:500].copy()

avail = set(adata.var_names)
candidate = {
    "T": ["CD3D", "CD3E", "CD8A"],
    "B": ["MS4A1", "CD79A", "CD79B"],
    "Mono": ["CD14", "LYZ", "FCGR3A"],
}
groups = {k: [g for g in v if g in avail] for k, v in candidate.items()}
groups = {k: v for k, v in groups.items() if v}
print(f"groups (available): {groups}")

print("=== sc.pl.dotplot return_fig=True ===")
dp = sc.pl.dotplot(adata, var_names=groups, groupby="louvain",
                    show=False, return_fig=True)
print(f"  ax_dict before make_figure: {dp.ax_dict}")
dp.make_figure()
print(f"  ax_dict after make_figure : {list(dp.ax_dict.keys()) if dp.ax_dict else 'still None'}")
fig = dp.fig
print(f"  fig.axes count: {len(fig.axes)}")
for i, ax in enumerate(fig.axes):
    pos = ax.get_position()
    xt = ax.get_xticks()
    xtl = [t.get_text() for t in ax.get_xticklabels()]
    print(f"    axes[{i}] x=[{pos.x0:.3f}-{pos.x1:.3f}] y=[{pos.y0:.3f}-{pos.y1:.3f}]  "
          f"xticks={list(xt)[:5]}  labels={xtl[:5]}")

mainax = dp.ax_dict["mainplot_ax"]
print(f"\nmainplot_ax xlim = {mainax.get_xlim()}")
print(f"mainplot_ax xticks = {mainax.get_xticks()}")
print(f"mainplot_ax xticklabels = {[t.get_text() for t in mainax.get_xticklabels()]}")
plt.close("all")

print("\n=== sc.pl.matrixplot return_fig=True ===")
mp = sc.pl.matrixplot(adata, var_names=groups, groupby="louvain",
                       cmap="magma", show=False, return_fig=True)
print(f"  ax_dict before make_figure: {mp.ax_dict}")
mp.make_figure()
print(f"  ax_dict after make_figure : {list(mp.ax_dict.keys()) if mp.ax_dict else 'still None'}")
mainax_mp = mp.ax_dict["mainplot_ax"]
print(f"  mainplot_ax xlim = {mainax_mp.get_xlim()}")
print(f"  mainplot_ax xticks = {mainax_mp.get_xticks()}")
plt.close("all")
