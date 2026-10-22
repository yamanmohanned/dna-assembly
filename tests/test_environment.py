"""Phase 0 smoke tests: the toolchain promised in docs/PLAN.md is present and the layout exists."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent

REQUIRED_MODULES = [
    "numpy",
    "numba",
    "pandas",
    "scipy",
    "Bio",
    "plotly",
    "matplotlib",
    "streamlit",
    "yaml",
    "pptx",
    "hypothesis",
]

REQUIRED_DIRS = [
    "genomeasm",
    "app",
    "configs",
    "data/raw",
    "data/processed",
    "results",
    "reports",
    "tests",
    "docs",
    "docs/phases",
]


def test_python_version_is_supported() -> None:
    assert sys.version_info >= (3, 11)


@pytest.mark.parametrize("name", REQUIRED_MODULES)
def test_required_module_imports(name: str) -> None:
    importlib.import_module(name)


def test_numba_compiles_and_matches_numpy() -> None:
    import numba
    import numpy as np

    @numba.njit(cache=False)
    def total(a: np.ndarray) -> int:
        s = 0
        for i in range(a.shape[0]):
            s += a[i]
        return s

    a = np.arange(1000, dtype=np.int64)
    assert total(a) == int(a.sum())


@pytest.mark.parametrize("rel", REQUIRED_DIRS)
def test_repository_layout(rel: str) -> None:
    assert (ROOT / rel).is_dir(), f"missing directory {rel}"


@pytest.mark.parametrize(
    "rel",
    ["docs/PLAN.md", "docs/DECISIONS.md", "docs/MASTER_PROMPT.md", "CLAUDE.md"],
)
def test_project_documents_exist(rel: str) -> None:
    assert (ROOT / rel).is_file(), f"missing document {rel}"
