"""The committed Phase 4 results: solver correctness validation."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
import yaml

from genomeasm import provenance

ROOT = Path(__file__).resolve().parent.parent
CONFIG_PATH = ROOT / "configs" / "phase4_correctness.yaml"
RESULTS_PATH = ROOT / "results" / "phase4" / "correctness.json"
POLICIES = ["random", "roundrobin", "roulette", "epsilon_greedy", "ucb1"]


@pytest.fixture(scope="module")
def committed() -> dict[str, Any]:
    return json.loads(RESULTS_PATH.read_text(encoding="utf-8"))


def test_every_check_passed(committed: dict[str, Any]) -> None:
    assert committed["all_passed"] is True and committed["n_failed"] == 0
    assert committed["n_checks"] == len(committed["checks"]) >= 10
    assert all(c["passed"] for c in committed["checks"])


def test_results_match_the_config_and_have_full_provenance(committed: dict[str, Any]) -> None:
    cfg = yaml.safe_load(CONFIG_PATH.read_text(encoding="utf-8"))
    assert committed["config"] == cfg, "configs/phase4_correctness.yaml changed: rerun"
    prov = committed["provenance"]
    assert prov["config_hash"] == provenance.config_hash(cfg) and prov["seed"] == cfg["seed"]
    assert len(prov["git"]["commit"]) == 40 and prov["git"]["dirty"] is False
    assert committed["lkh_available"] is True  # LKH-3 was really run, never replaced


def test_mandatory_checks_are_present(committed: dict[str, Any]) -> None:
    ids = {c["id"] for c in committed["checks"]}
    expected = {
        "delta.equals_full_recomputation",
        "exact.brute_force_equals_held_karp",
        "solvers.valid_permutations",
        "reproducibility.same_seed_same_result",
        "preview.valid_and_verified",
        "optimality.tabu",
        "optimality.lkh3",
        *(f"optimality.alns_{p}" for p in POLICIES),
    }
    assert expected <= ids


def test_delta_check_covered_many_moves_without_a_mismatch(committed: dict[str, Any]) -> None:
    assert committed["delta"]["mismatches"] == 0
    assert committed["delta"]["moves_checked"] >= 100_000


def test_small_instance_table_is_complete_and_heuristics_reach_the_optimum(
    committed: dict[str, Any],
) -> None:
    small = committed["small"]
    assert len(small["instances"]) == 5 * 12 + 12
    solvers = small["solvers"]
    for name in ("greedy", "nearest_neighbour", "tabu", "lkh3", *(f"alns_{p}" for p in POLICIES)):
        assert name in solvers
    assert sum(name.startswith("alns_fixed_") for name in solvers) == 6
    required = committed["config"]["tolerances"]["min_success_rate"]
    for name in ("tabu", "lkh3", *(f"alns_{p}" for p in POLICIES)):
        assert solvers[name]["success_rate"] >= required, name
    for stats in solvers.values():  # the gap is documented for every solver
        assert 0.0 <= stats["mean_gap_percent"] <= stats["max_gap_percent"] <= 100.0


def test_local_search_never_ends_below_its_greedy_start(committed: dict[str, Any]) -> None:
    preview = committed["preview"]["solvers"]
    assert committed["preview"]["instance"]["n_reads"] > 400
    for name in ("tabu", *(f"alns_{p}" for p in POLICIES)):
        assert preview[name]["fitness"] >= preview["greedy"]["fitness"], name
    assert "lkh3" in preview and preview["lkh3"]["info"]["lkh_reported_cost"] is not None
