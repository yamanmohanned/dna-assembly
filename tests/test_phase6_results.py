"""Checks on the committed timing of the application (results/phase6) and PHASE_6.md claims."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pytest

from genomeasm.appkit import lab

ROOT = Path(__file__).resolve().parents[1]
TIMING = ROOT / "results" / "phase6" / "interactive_timing.json"

pytestmark = pytest.mark.skipif(not TIMING.exists(), reason="Phase 6 timing not generated")


@pytest.fixture(scope="module")
def timing() -> dict[str, Any]:
    loaded = json.loads(TIMING.read_text(encoding="utf-8"))
    assert isinstance(loaded, dict)
    return loaded


def test_every_measured_interactive_path_ends_before_the_thirty_second_rule(
    timing: dict[str, Any],
) -> None:
    assert timing["limit_seconds"] == lab.MAX_TOTAL_SECONDS == 30.0
    totals = {name: s["total"] for name, s in timing["scenarios"].items()}
    assert timing["worst_total_seconds"] == pytest.approx(max(totals.values()))
    assert timing["ok"] is True and timing["worst_total_seconds"] < 30.0
    for name, scenario in timing["scenarios"].items():
        parts = {k: v for k, v in scenario.items() if k in ("build", "solve", "render", "simulate")}
        assert scenario["total"] >= sum(parts.values()) - 1e-6, name


def test_the_measurement_used_the_limits_the_interface_enforces(timing: dict[str, Any]) -> None:
    limits = timing["interface_limits"]
    assert limits == {
        "max_reads": lab.MAX_READS,
        "max_window": lab.MAX_WINDOW,
        "max_solver_seconds": lab.MAX_SECONDS,
    }
    # the worst instance really is the largest one the interface allows
    solver_run = timing["scenarios"]["solver_lab_alns_20s"]
    assert 0.95 * lab.MAX_READS <= solver_run["reads"] <= lab.MAX_READS


def test_the_solvers_searched_for_the_whole_budget_and_were_not_starved_by_the_charts(
    timing: dict[str, Any],
) -> None:
    for name in ("solver_lab_alns_20s", "solver_lab_tabu_cl_20s"):
        scenario = timing["scenarios"][name]
        assert 19.9 <= scenario["solve"] <= 21.0, name
        assert scenario["iterations"] > 5000, name  # the first live version managed about 10
        # the live redraws cost about a tenth of the run
        assert scenario["live_redraws"] * scenario["redraw_seconds_each"] < 0.2 * scenario["solve"]
    assembly = timing["scenarios"]["assembly_viewer_alns_15s"]
    assert 15.0 <= assembly["solve"] <= 15.5


def test_lkh3_is_timed_or_reported_unavailable_never_substituted(timing: dict[str, Any]) -> None:
    lkh = timing["scenarios"]["solver_lab_lkh3_15s"]
    assert ("unavailable" in lkh) != ("solve" in lkh)
    if "solve" in lkh:
        assert lkh["solve"] <= 17.0 + 3.0  # 15 s search plus the overhead of its calibration run


def test_the_timing_has_provenance_of_a_clean_commit(timing: dict[str, Any]) -> None:
    git = timing["provenance"]["git"]
    assert re.fullmatch(r"[0-9a-f]{40}", git["commit"]) and git["dirty"] is False
    assert timing["provenance"]["machine"]["cpu_count"] >= 1
