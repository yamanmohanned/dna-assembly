"""Everything the Reproducibility page shows: configuration, seeds, libraries, commit, checksums."""

from __future__ import annotations

from typing import Any

import yaml

from genomeasm import provenance
from genomeasm.appkit.paths import CONFIG, RAW, RESULTS, ROOT
from genomeasm.io import read_sums, sha256_file
from genomeasm.solvers.lkh import lkh_binary


def config_info() -> dict[str, Any]:
    """The Phase 5 configuration: its file hash, its canonical hash and the seed scheme."""
    cfg = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    return {
        "path": str(CONFIG.relative_to(ROOT)).replace("\\", "/"),
        "file_sha256": sha256_file(CONFIG),
        "config_hash": provenance.config_hash(cfg),
        "master_seed": cfg["seed"],
        "seeds": cfg["seeds"],
        "budget_seconds": cfg["budget_seconds"],
        "tuning_trials": cfg["tuning"]["trials"],
        "seed_scheme": (
            "read-set and solver seeds = SeedSequence([master seed, stage, instance index, seed "
            "index]).generate_state(2); tuning seeds start at 1000, evaluation seeds are 0 to 29"
        ),
    }


def environment() -> dict[str, Any]:
    """Machine, libraries, git state and solver availability of this very process."""
    binary = lkh_binary()
    return {
        "machine": provenance.machine_info(),
        "libraries": provenance.library_versions(),
        "git": provenance.git_info(ROOT),
        "lkh3_binary": None if binary is None else str(binary),
    }


def verify_checksums() -> list[dict[str, Any]]:
    """Compare every data file with ``data/raw/SHA256SUMS`` (the same check as ``load_genome``)."""
    rows = []
    for name, expected in sorted(read_sums(RAW).items()):
        path = RAW / name
        actual = sha256_file(path) if path.is_file() else None
        rows.append(
            {"file": name, "expected": expected, "actual": actual, "ok": actual == expected}
        )
    return rows


def result_files() -> list[dict[str, Any]]:
    """Fingerprints of the committed Phase 5 result files."""
    rows = []
    for path in sorted(p for p in RESULTS.rglob("*") if p.is_file()):
        rows.append(
            {
                "file": str(path.relative_to(ROOT)).replace("\\", "/"),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    return rows
