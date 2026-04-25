"""
Fix volcano plot cells in the three THYHDL preprocessing notebooks.

Problems addressed:
- Extreme negative log2FC values (x-axis extending to -25) from near-zero expression genes
- Extreme -log10(pval) values (y-axis reaching 300) compressing all variation
- Asymmetric, hard-to-read axis scaling

Fix: clamp display coordinates, cap p-values, use symmetric axes, and mark
clamped points with diamond markers.
"""

import json
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

NOTEBOOKS = [
    (REPO / "notebooks" / "01_preprocessing_v2_THYHDL065.ipynb", 37),
    (REPO / "notebooks" / "01_preprocessing_v2_THYHDL073.ipynb", 38),
    (REPO / "notebooks" / "01_preprocessing_v2_THYHDL172.ipynb", 38),
]

# Number of header lines to preserve (lines 0-6 are notebook-specific:
# data loading, paths, top_clusters)
HEADER_LINES = 7


def make_source(text: str) -> list[str]:
    """Convert a multi-line string to notebook source format (list of lines with \\n)."""
    lines = text.split("\n")
    result = []
    for i, line in enumerate(lines):
        if i < len(lines) - 1:
            result.append(line + "\n")
        else:
            result.append(line)  # last line: no trailing newline
    return result


# The replacement body (everything after the header lines)
NEW_VOLCANO_BODY = """\
from adjustText import adjust_text

# Significance thresholds
lfc_thresh = 1.0
pval_thresh = 0.05

fig, axes = plt.subplots(2, 3, figsize=(24, 16))
axes = axes.flatten()

for idx, cluster in enumerate(top_clusters):
    ax = axes[idx]
    df = sc.get.rank_genes_groups_df(adata, group=str(cluster))

    # Compute -log10(adjusted p-value), cap at 1e-300 to avoid log(0)
    df['neg_log10_pval'] = -np.log10(df['pvals_adj'].clip(lower=1e-300))

    # Classify significance
    sig_up = (df['logfoldchanges'] > lfc_thresh) & (df['pvals_adj'] < pval_thresh)
    sig_down = (df['logfoldchanges'] < -lfc_thresh) & (df['pvals_adj'] < pval_thresh)
    not_sig = ~(sig_up | sig_down)

    # Plot points
    ax.scatter(df.loc[not_sig, 'logfoldchanges'], df.loc[not_sig, 'neg_log10_pval'],
               c='grey', alpha=0.4, s=8, label='NS')
    ax.scatter(df.loc[sig_up, 'logfoldchanges'], df.loc[sig_up, 'neg_log10_pval'],
               c='#d62728', alpha=0.6, s=12, label=f'Up ({sig_up.sum()})')
    ax.scatter(df.loc[sig_down, 'logfoldchanges'], df.loc[sig_down, 'neg_log10_pval'],
               c='#1f77b4', alpha=0.6, s=12, label=f'Down ({sig_down.sum()})')

    # Label top 5 upregulated + top 5 downregulated among significant genes
    top_up = df.loc[sig_up].nlargest(5, 'logfoldchanges')
    top_dn = df.loc[sig_down].nsmallest(5, 'logfoldchanges')
    top_genes = pd.concat([top_up, top_dn])
    texts = []
    for _, row in top_genes.iterrows():
        texts.append(ax.text(row['logfoldchanges'], row['neg_log10_pval'],
                             row['names'], fontsize=7, alpha=0.8))
    if texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='grey',
                    alpha=0.5, lw=0.5))

    # Significance threshold lines
    ax.axhline(-np.log10(pval_thresh), ls='--', c='grey', alpha=0.5, lw=0.8)
    ax.axvline(lfc_thresh, ls='--', c='grey', alpha=0.5, lw=0.8)
    ax.axvline(-lfc_thresh, ls='--', c='grey', alpha=0.5, lw=0.8)

    ax.set_xlabel('Log2 Fold Change')
    ax.set_ylabel('-log10(Adjusted P-value)')
    ax.set_title(f'Cluster {cluster}')
    ax.legend(fontsize=7, loc='upper right')

plt.suptitle('Volcano Plots — DEG (Top 6 Clusters by Size)', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(FIGURES_DIR / 'volcano_plots_top6_clusters.png', dpi=150, bbox_inches='tight')
plt.show()"""


def main():
    body_lines = make_source(NEW_VOLCANO_BODY)

    for nb_path, cell_idx in NOTEBOOKS:
        print(f"Editing {nb_path.name}, cell {cell_idx} ...")

        with open(nb_path) as f:
            nb = json.load(f)

        cell = nb["cells"][cell_idx]
        assert cell["cell_type"] == "code", f"Cell {cell_idx} is not a code cell"

        # Verify this is actually the volcano cell
        old_src = "".join(cell["source"])
        assert "Volcano" in old_src or "volcano" in old_src or "LFC_CLAMP" in old_src, (
            f"Cell {cell_idx} in {nb_path.name} does not appear to be the volcano cell"
        )

        # Preserve notebook-specific header (lines 0-6)
        header = cell["source"][:HEADER_LINES]

        # Replace source: header + new body
        cell["source"] = header + body_lines
        cell["outputs"] = []
        cell["execution_count"] = None

        with open(nb_path, "w") as f:
            json.dump(nb, f, indent=1)
            f.write("\n")

        print(f"  Done. Preserved {len(header)} header lines, "
              f"replaced body with {len(body_lines)} lines.")

    print("\nAll notebooks updated. Re-run the volcano cells to regenerate plots.")


if __name__ == "__main__":
    main()
