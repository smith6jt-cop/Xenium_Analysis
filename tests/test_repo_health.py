"""
Repo-health tests.

These run on every PR and validate that the repo is in a shape that will
actually work on our two target platforms:

- A Linux workstation with the ``xenium_analysis`` conda env.
- The HiPerGator cluster (SLURM + conda).

They are cheap: no heavy imports, no real data required.  The goal is to
catch regressions in things that are annoying to debug remotely:

- ``environment.yml`` parses and pins a supported Python.
- SLURM scripts use LF line endings, start with a shebang, and reference the
  correct conda env.
- Every notebook is valid JSON with a known kernel.
- The ``utils`` package imports its core submodules cleanly.
- ``config.ini`` parses and defines the sections the rest of the code reads.
"""

from __future__ import annotations

import configparser
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent

SLURM_DIR = REPO / "scripts" / "slurm"
NOTEBOOKS_DIR = REPO / "notebooks"
ENV_FILE = REPO / "environment.yml"
CONFIG_FILE = REPO / "config.ini"


# ---------------------------------------------------------------------------
# environment.yml
# ---------------------------------------------------------------------------

def _parse_env_yaml(path: Path) -> dict:
    """Hand-rolled minimal YAML parser for the simple environment.yml layout.

    We avoid pulling in PyYAML for CI so the test file has no runtime
    dependency beyond stdlib + pytest.  The parser handles:
      name: <value>
      channels: [- items]
      dependencies: mixed list of strings and nested pip mapping.
    """
    text = path.read_text().splitlines()
    result: dict = {"channels": [], "dependencies": [], "pip": []}
    section = None
    in_pip = False
    for raw in text:
        line = raw.rstrip()
        if not line or line.lstrip().startswith("#"):
            continue
        if line.startswith("name:"):
            result["name"] = line.split(":", 1)[1].strip()
            section = None
            continue
        if line.startswith("channels:"):
            section = "channels"
            continue
        if line.startswith("dependencies:"):
            section = "dependencies"
            continue
        stripped = line.lstrip()
        if stripped.startswith("- ") and section == "channels":
            result["channels"].append(stripped[2:].strip())
        elif stripped == "- pip:" and section == "dependencies":
            in_pip = True
        elif in_pip and stripped.startswith("- "):
            # Four-space indent level means still inside pip: block.
            if line.startswith("    - "):
                result["pip"].append(line[6:].strip())
            else:
                in_pip = False
                result["dependencies"].append(stripped[2:].strip())
        elif stripped.startswith("- ") and section == "dependencies":
            in_pip = False
            result["dependencies"].append(stripped[2:].strip())
    return result


class TestEnvironmentYaml:
    def test_file_exists(self):
        assert ENV_FILE.exists(), "environment.yml missing"

    def test_parses(self):
        env = _parse_env_yaml(ENV_FILE)
        assert env["name"] == "xenium_analysis"
        assert env["channels"], "no conda channels declared"

    def test_python_pinned(self):
        env = _parse_env_yaml(ENV_FILE)
        py_entries = [d for d in env["dependencies"] if d.startswith("python")]
        assert py_entries, "environment.yml must pin a python version"
        # Must be pinned with '=' (conda syntax) and to 3.10+
        spec = py_entries[0]
        assert "=" in spec, f"python not pinned: {spec}"
        version = spec.split("=", 1)[1]
        major, minor = version.split(".")[:2]
        assert int(major) == 3 and int(minor) >= 10, f"python {version} too old"

    def test_conda_forge_listed(self):
        env = _parse_env_yaml(ENV_FILE)
        assert "conda-forge" in env["channels"]

    def test_core_deps_present(self):
        env = _parse_env_yaml(ENV_FILE)
        deps = {d.split("=")[0].split("<")[0].split(">")[0].strip() for d in env["dependencies"]}
        for pkg in ("scanpy", "squidpy", "anndata", "zarr"):
            assert pkg in deps, f"{pkg} missing from environment.yml"


# ---------------------------------------------------------------------------
# SLURM scripts
# ---------------------------------------------------------------------------

SLURM_SCRIPTS = sorted(SLURM_DIR.glob("*.sh")) if SLURM_DIR.exists() else []

# Not every shell script under ``scripts/slurm/`` is itself a ``sbatch`` job
# -- some are orchestration / submitter scripts that only call ``sbatch``
# under the hood (``run_all_samples.sh``).  The ``#SBATCH`` header marks
# the ones that should satisfy the full set of batch-job invariants.
def _is_batch_job(script: Path) -> bool:
    text = script.read_text(errors="replace")
    return "#SBATCH" in text


BATCH_JOB_SCRIPTS = [s for s in SLURM_SCRIPTS if _is_batch_job(s)]


@pytest.mark.parametrize("script", SLURM_SCRIPTS, ids=lambda p: p.name)
class TestShellHygiene:
    """Basic shell hygiene that applies to any ``.sh`` under scripts/slurm/."""

    def test_shebang(self, script: Path):
        first = script.read_bytes().splitlines()[0]
        assert first.startswith(b"#!"), f"{script.name} missing shebang"
        assert b"bash" in first, f"{script.name} should use bash"

    def test_lf_line_endings(self, script: Path):
        """Scripts edited on Windows can pick up CRLF which SLURM rejects."""
        content = script.read_bytes()
        assert b"\r\n" not in content, f"{script.name} has CRLF line endings"


@pytest.mark.parametrize("script", BATCH_JOB_SCRIPTS, ids=lambda p: p.name)
class TestSlurmBatchJobs:
    """Invariants that only apply to scripts with ``#SBATCH`` headers.

    Orchestration-only scripts (no ``#SBATCH``) still get the basic shell
    hygiene checks via :class:`TestShellHygiene`.
    """

    def test_has_job_name(self, script: Path):
        text = script.read_text()
        assert "#SBATCH --job-name=" in text, f"{script.name} missing --job-name"

    def test_has_resource_directives(self, script: Path):
        text = script.read_text()
        for directive in ("--time=", "--mem", "--cpus-per-task"):
            assert directive in text, f"{script.name} missing {directive}"

    def test_activates_conda_env(self, script: Path):
        text = script.read_text()
        assert "xenium_analysis" in text, (
            f"{script.name} does not reference the xenium_analysis conda env"
        )

    def test_outputs_to_logs_dir(self, script: Path):
        """Output files should land in a dedicated logs/ dir (git-ignored)."""
        text = script.read_text()
        assert "--output=logs/" in text, f"{script.name} should write to logs/"


def test_slurm_dir_has_scripts():
    assert SLURM_SCRIPTS, "expected at least one SLURM script under scripts/slurm/"


def test_slurm_dir_has_batch_jobs():
    assert BATCH_JOB_SCRIPTS, (
        "expected at least one ``#SBATCH``-annotated script under scripts/slurm/"
    )


# ---------------------------------------------------------------------------
# Notebooks
# ---------------------------------------------------------------------------

NOTEBOOKS = sorted(NOTEBOOKS_DIR.glob("*.ipynb")) if NOTEBOOKS_DIR.exists() else []


@pytest.mark.parametrize("nb_path", NOTEBOOKS, ids=lambda p: p.name)
class TestNotebookValidity:
    def test_valid_json(self, nb_path: Path):
        with nb_path.open() as f:
            nb = json.load(f)
        assert "cells" in nb
        assert "metadata" in nb

    def test_has_kernelspec(self, nb_path: Path):
        with nb_path.open() as f:
            nb = json.load(f)
        kernelspec = nb.get("metadata", {}).get("kernelspec", {})
        assert kernelspec, f"{nb_path.name} missing kernelspec"
        # We expect a Python kernel; the exact name varies by host.
        assert "python" in kernelspec.get("language", "").lower()

    def test_nbformat_version(self, nb_path: Path):
        with nb_path.open() as f:
            nb = json.load(f)
        # 4.x is what JupyterLab writes and what nbconvert on HiPerGator needs.
        assert nb.get("nbformat") == 4


def test_notebooks_dir_has_pipeline_notebooks():
    names = {p.name for p in NOTEBOOKS}
    required = {
        "01_preprocessing.ipynb",
        "02_phenotyping.ipynb",
        "03_spatial_analysis.ipynb",
        "04_group_comparisons.ipynb",
        "05_tissue_comparisons.ipynb",
    }
    missing = required - names
    assert not missing, f"pipeline notebooks missing: {sorted(missing)}"


# ---------------------------------------------------------------------------
# utils package imports
# ---------------------------------------------------------------------------

class TestUtilsImports:
    def test_xenium_explorer_export_imports_clean(self):
        # This module intentionally has no heavy deps so it must always
        # import without scanpy being installed.
        import importlib

        mod = importlib.import_module("utils.xenium_explorer_export")
        for name in (
            "export_for_xenium_explorer",
            "export_groups_to_csv",
            "export_groups_to_zarr",
            "generate_color_palette",
        ):
            assert hasattr(mod, name), f"utils.xenium_explorer_export.{name} missing"


# ---------------------------------------------------------------------------
# config.ini
# ---------------------------------------------------------------------------

class TestConfigIni:
    def test_parses(self):
        cfg = configparser.ConfigParser(inline_comment_prefixes=("#",))
        cfg.read(CONFIG_FILE)
        for section in (
            "general", "paths", "samples", "preprocessing",
            "clustering", "phenotyping", "spatial",
            "differential_expression", "integration", "hipergator",
        ):
            assert cfg.has_section(section), f"config.ini missing [{section}]"

    def test_preprocessing_thresholds_are_numeric(self):
        cfg = configparser.ConfigParser(inline_comment_prefixes=("#",))
        cfg.read(CONFIG_FILE)
        for key in ("min_genes_per_cell", "min_cells_per_gene",
                    "min_counts_per_cell"):
            cfg.getint("preprocessing", key)
        for key in ("max_mt_percent",):
            cfg.getfloat("preprocessing", key)
