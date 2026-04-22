"""
Structural tests for the preprocessing pipeline refactor.
Verifies that the refactoring script correctly modified all three notebooks.
"""

import json
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent

NOTEBOOKS = [
    {
        "path": REPO / "notebooks" / "01_preprocessing_v2_THYHDL065.ipynb",
        "expected_cells": 51,  # was 48 + 3 new
        "old_leiden_key": "leiden_0.5",
    },
    {
        "path": REPO / "notebooks" / "01_preprocessing_v2_THYHDL073.ipynb",
        "expected_cells": 52,  # was 49 + 3 new
        "old_leiden_key": "leiden_0.5",
    },
    {
        "path": REPO / "notebooks" / "01_preprocessing_v2_THYHDL172.ipynb",
        "expected_cells": 52,  # was 49 + 3 new
        "old_leiden_key": "leiden_0.2",
    },
]


def load_notebook(path):
    with open(path) as f:
        return json.load(f)


def cell_source(nb, idx):
    return "".join(nb["cells"][idx]["source"])


@pytest.fixture(params=NOTEBOOKS, ids=lambda c: c["path"].name)
def notebook_config(request):
    config = request.param
    config["nb"] = load_notebook(config["path"])
    return config


class TestNotebookStructure:
    def test_valid_json(self, notebook_config):
        """Notebook loads as valid JSON."""
        assert notebook_config["nb"] is not None

    def test_cell_count(self, notebook_config):
        """Notebook has correct number of cells after insertion."""
        actual = len(notebook_config["nb"]["cells"])
        expected = notebook_config["expected_cells"]
        assert actual == expected, f"Expected {expected} cells, got {actual}"

    def test_gc_import(self, notebook_config):
        """Imports cell includes 'import gc'."""
        src = cell_source(notebook_config["nb"], 1)
        assert "import gc" in src


class TestCell18NormHVG:
    def test_target_sum_1e4(self, notebook_config):
        """Cell 18 uses target_sum=1e4."""
        src = cell_source(notebook_config["nb"], 18)
        assert "target_sum=1e4" in src

    def test_no_target_sum_none(self, notebook_config):
        """Cell 18 does NOT use target_sum=None."""
        src = cell_source(notebook_config["nb"], 18)
        assert "target_sum=None" not in src

    def test_hvg_on_counts_layer(self, notebook_config):
        """Cell 18 runs HVG on raw counts layer."""
        src = cell_source(notebook_config["nb"], 18)
        assert "layer='counts'" in src

    def test_adaptive_n_top_genes(self, notebook_config):
        """Cell 18 computes adaptive n_top_genes."""
        src = cell_source(notebook_config["nb"], 18)
        assert "n_top_genes=n_top_genes" in src
        assert "n_top_genes = min(300" in src

    def test_no_pca_in_cell18(self, notebook_config):
        """Cell 18 does NOT contain PCA (moved to Cell 20)."""
        src = cell_source(notebook_config["nb"], 18)
        assert "sc.pp.pca" not in src

    def test_no_umap_in_cell18(self, notebook_config):
        """Cell 18 does NOT contain UMAP (moved to Cell 21)."""
        src = cell_source(notebook_config["nb"], 18)
        assert "sc.tl.umap" not in src


class TestCell19AdaptiveParams:
    def test_leiden_resolution(self, notebook_config):
        """Cell 19 defines LEIDEN_RESOLUTION."""
        src = cell_source(notebook_config["nb"], 19)
        assert "LEIDEN_RESOLUTION" in src

    def test_leiden_key(self, notebook_config):
        """Cell 19 defines LEIDEN_KEY."""
        src = cell_source(notebook_config["nb"], 19)
        assert "LEIDEN_KEY" in src

    def test_min_per_cluster(self, notebook_config):
        """Cell 19 defines MIN_PER_CLUSTER adaptively."""
        src = cell_source(notebook_config["nb"], 19)
        assert "MIN_PER_CLUSTER" in src
        assert "n_cells < 50000" in src

    def test_subsample_fraction(self, notebook_config):
        """Cell 19 defines SUBSAMPLE_FRACTION."""
        src = cell_source(notebook_config["nb"], 19)
        assert "SUBSAMPLE_FRACTION" in src


class TestCell20PCA:
    def test_scale_on_hvg_subset(self, notebook_config):
        """Cell 20 scales HVG subset, not full adata."""
        src = cell_source(notebook_config["nb"], 20)
        assert "sc.pp.scale(adata_hvg" in src

    def test_pca_on_subset(self, notebook_config):
        """Cell 20 runs PCA on HVG subset."""
        src = cell_source(notebook_config["nb"], 20)
        assert "sc.pp.pca(adata_hvg" in src

    def test_copy_back_embeddings(self, notebook_config):
        """Cell 20 copies PCA embeddings back to full adata."""
        src = cell_source(notebook_config["nb"], 20)
        assert "adata.obsm['X_pca']" in src
        assert "adata.uns['pca']" in src

    def test_cleanup(self, notebook_config):
        """Cell 20 deletes temporary HVG copy."""
        src = cell_source(notebook_config["nb"], 20)
        assert "del adata_hvg" in src
        assert "gc.collect()" in src


class TestCell21Clustering:
    def test_leiden_uses_variable(self, notebook_config):
        """Cell 21 uses LEIDEN_RESOLUTION and LEIDEN_KEY variables."""
        src = cell_source(notebook_config["nb"], 21)
        assert "resolution=LEIDEN_RESOLUTION" in src
        assert "key_added=LEIDEN_KEY" in src

    def test_umap(self, notebook_config):
        """Cell 21 contains UMAP."""
        src = cell_source(notebook_config["nb"], 21)
        assert "sc.tl.umap" in src

    def test_cluster_size_warnings(self, notebook_config):
        """Cell 21 warns about small clusters."""
        src = cell_source(notebook_config["nb"], 21)
        assert "small cluster" in src


class TestDownstreamCells:
    def test_no_hardcoded_leiden_key(self, notebook_config):
        """No hardcoded leiden key strings remain after Cell 19."""
        old_key = notebook_config["old_leiden_key"]
        nb = notebook_config["nb"]
        for i in range(20, len(nb["cells"])):
            src = cell_source(nb, i)
            # Allow it in output cells (which we may not have cleared) and markdown
            if nb["cells"][i]["cell_type"] != "code":
                continue
            # Check source lines only (not outputs)
            for line in nb["cells"][i]["source"]:
                # Skip comment lines
                stripped = line.strip()
                if stripped.startswith("#"):
                    continue
                assert f"'{old_key}'" not in line and f'"{old_key}"' not in line, \
                    f"Found hardcoded '{old_key}' in cell {i}: {line.strip()}"

    def test_de_cell_no_x_restore(self, notebook_config):
        """DE cell does not restore .X from log_normalized layer."""
        nb = notebook_config["nb"]
        # DE cell is now at original_idx + 3
        for i in range(22, len(nb["cells"])):
            src = cell_source(nb, i)
            if "rank_genes_groups" in src and "wilcoxon" in src:
                assert "adata.X = adata.layers" not in src, \
                    f"DE cell {i} still restores .X from layers"
                break

    def test_de_cell_no_hardcoded_de_params(self, notebook_config):
        """DE cell does not redefine MAX_PER_CLUSTER or MIN_PER_CLUSTER."""
        nb = notebook_config["nb"]
        for i in range(22, len(nb["cells"])):
            src = cell_source(nb, i)
            if "rank_genes_groups" in src and "wilcoxon" in src:
                assert "MAX_PER_CLUSTER = 200" not in src, \
                    f"DE cell {i} still defines MAX_PER_CLUSTER"
                assert "MIN_PER_CLUSTER = 50" not in src, \
                    f"DE cell {i} still defines MIN_PER_CLUSTER"
                break

    def test_subsample_uses_variable(self, notebook_config):
        """Subsample cell uses SUBSAMPLE_FRACTION instead of 0.5."""
        nb = notebook_config["nb"]
        for i in range(22, len(nb["cells"])):
            src = cell_source(nb, i)
            if "subsample" in src.lower() and "sdata.tables" in src:
                assert "SUBSAMPLE_FRACTION" in src, \
                    f"Subsample cell {i} doesn't use SUBSAMPLE_FRACTION"
                break

    def test_summary_uses_leiden_key(self, notebook_config):
        """Summary cell uses LEIDEN_KEY variable."""
        nb = notebook_config["nb"]
        for i in range(len(nb["cells"]) - 5, len(nb["cells"])):
            src = cell_source(nb, i)
            if "Number of clusters" in src:
                assert "LEIDEN_KEY" in src, \
                    f"Summary cell {i} doesn't use LEIDEN_KEY"
                break
