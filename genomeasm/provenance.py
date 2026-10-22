"""Provenance: what every result file records so that it can be traced and reproduced (rule R3).

A result file is only trustworthy if we know which settings, which code version, which random
seed and which machine produced it. ``collect`` gathers exactly that.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
from collections.abc import Mapping
from datetime import UTC, datetime
from importlib import metadata
from pathlib import Path
from typing import Any

LIBRARIES = ("numpy", "scipy", "pandas", "numba", "biopython", "matplotlib", "plotly", "PyYAML")


def config_hash(config: Mapping[str, Any]) -> str:
    """SHA-256 of the configuration in canonical JSON form (sorted keys, no spaces)."""
    canonical = json.dumps(config, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def _git(repo_dir: Path, *args: str) -> str | None:
    try:
        done = subprocess.run(
            ["git", *args], cwd=repo_dir, capture_output=True, text=True, check=False, timeout=30
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return done.stdout.strip() if done.returncode == 0 else None


GENERATED_PATHS = ("results", "reports")


def git_info(repo_dir: Path) -> dict[str, Any]:
    """Current commit and whether the code has uncommitted changes.

    Generated outputs (``results/`` and ``reports/``) are ignored when deciding ``dirty``, so
    results written after committing the code are recorded as coming from a clean commit.
    """
    commit = _git(repo_dir, "rev-parse", "HEAD")
    excludes = [f":(exclude){path}" for path in GENERATED_PATHS]
    status = _git(repo_dir, "status", "--porcelain", "--", ".", *excludes)
    return {"commit": commit, "dirty": None if status is None else bool(status)}


def library_versions() -> dict[str, str | None]:
    versions: dict[str, str | None] = {}
    for name in LIBRARIES:
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = None
    return versions


def machine_info() -> dict[str, Any]:
    return {
        "platform": platform.platform(),
        "processor": platform.processor(),
        "python": platform.python_version(),
        "cpu_count": os.cpu_count(),
    }


def collect(
    config: Mapping[str, Any], *, seed: int | None, repo_dir: Path | None = None
) -> dict[str, Any]:
    """Provenance block: config hash, git commit, seed, library versions, machine, timestamp."""
    return {
        "config_hash": config_hash(config),
        "git": git_info(repo_dir or Path.cwd()),
        "seed": seed,
        "library_versions": library_versions(),
        "machine": machine_info(),
        "timestamp_utc": datetime.now(UTC).isoformat(timespec="seconds"),
    }
