"""Checks on the committed Phase 5 results (results/phase5) and the claims of PHASE_5.md.

Structural tests follow from the design (one result per solver and read set, budgets, separation of
tuning and evaluation). The ``finding`` tests are regression tests of statements made in the report:
they are not theorems, and if the experiments are regenerated and one fails, the report must change.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from genomeasm.experiments import analysis, design, plan, runner

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "phase5"
FIGURES = ROOT / "reports" / "figures"

pytestmark = pytest.mark.skipif(
    not (RESULTS / "analysis" / "tests.csv").exists(), reason="Phase 5 results not generated"
)

CFG = design.load_config(ROOT / "configs" / "phase5_experiments.yaml")
ALL_SOLVERS = [*plan.MAIN_ALGORITHMS]
TIMED = [a for a in ALL_SOLVERS if a not in ("greedy", "nearest_neighbour")]


@pytest.fixture(scope="module")
def records() -> list[dict[str, Any]]:
    return runner.read_records(RESULTS / "main_runs.jsonl.gz")


@pytest.fixture(scope="module")
def runs(records: list[dict[str, Any]]) -> pd.DataFrame:
    return analysis.load_runs(records)[0]


@pytest.fixture(scope="module")
def tests_table() -> pd.DataFrame:
    return pd.read_csv(RESULTS / "analysis" / "tests.csv")


@pytest.fixture(scope="module")
def summary() -> dict[str, Any]:
    return json.loads((RESULTS / "tuning_summary.json").read_text(encoding="utf-8"))


def _pooled(tests: pd.DataFrame, a: str, b: str) -> pd.Series:
    rows = tests[(tests["scope"] == "pooled") & (tests["a"] == a) & (tests["b"] == b)]
    assert len(rows) == 1, (a, b)
    return rows.iloc[0]


# ---- structure ----------------------------------------------------------------------------------


def test_every_solver_has_exactly_one_result_per_read_set(runs: pd.DataFrame) -> None:
    n_sets = len(CFG["evaluation_instances"]) * int(CFG["seeds"]["evaluation"])
    assert n_sets == 300
    assert set(runs["algorithm"]) == set(ALL_SOLVERS)
    for name, part in runs.groupby("algorithm"):
        keys = set(zip(part["instance"], part["seed_index"], strict=True))
        assert len(part) == len(keys) == n_sets, name
    ids = {s["id"] for s in CFG["evaluation_instances"]}
    assert set(runs["instance"]) == ids
    assert set(runs["seed_index"]) == set(range(int(CFG["seeds"]["evaluation"])))


def test_the_error_records_of_the_wsl_crash_are_all_repaired(records: list[dict[str, Any]]) -> None:
    errors = [r for r in records if r["kind"] == "error"]
    have = {
        (r["instance"], r["seed_index"])
        for r in records
        if r["kind"] == "run" and r["requested"] == "lkh3"
    }
    assert {e["algorithm"] for e in errors} == {"lkh3"} and len(errors) > 0
    assert all((e["instance"], e["seed_index"]) in have for e in errors)  # decision D-044


def test_timed_solvers_stop_at_the_budget(records: list[dict[str, Any]]) -> None:
    budget = float(CFG["budget_seconds"])
    seconds = pd.DataFrame(
        [(r["requested"], r["seconds"]) for r in records if r["kind"] == "run"],
        columns=["algorithm", "seconds"],
    )
    for name, part in seconds.groupby("algorithm"):
        if name in ("greedy", "nearest_neighbour"):
            assert part["seconds"].max() < 0.5, name
        elif name == "lkh3":
            # LKH-3's time is set by a calibration run, so it is not exact (limitation, D-041)
            assert abs(part["seconds"].median() - budget) < 0.15
            assert part["seconds"].max() < 1.1 * budget
        else:
            assert budget <= part["seconds"].min() and part["seconds"].max() < budget + 0.1, name


def test_tuning_and_evaluation_read_sets_are_separate(summary: dict[str, Any]) -> None:
    tuning = [
        r for r in runner.read_records(RESULTS / "tuning_runs.jsonl.gz") if r.get("kind") == "run"
    ]
    assert {r["instance"] for r in tuning} == {s["id"] for s in CFG["tuning_instances"]}
    assert min(r["seed_index"] for r in tuning) >= 1000
    assert {r["requested"] for r in tuning} == set(plan.ALGORITHMS_TUNED)
    for name in plan.ALGORITHMS_TUNED:
        n_trials = len(plan.trial_list(CFG, name))
        assert {r["trial"] for r in tuning if r["requested"] == name} == set(range(n_trials))


def test_best_known_is_the_best_of_all_solvers_and_gaps_are_not_negative(
    runs: pd.DataFrame,
) -> None:
    assert (runs["gap"] >= 0).all()
    best = runs.groupby("task")["fitness"].transform("max")
    assert (runs["best_known"] == best).all()
    assert (runs[runs["algorithm"] == "lkh3"]["gap"] < 1e-9).mean() > 0.2  # it often is the best


def test_local_searches_never_return_less_than_their_greedy_start(runs: pd.DataFrame) -> None:
    greedy = runs[runs["algorithm"] == "greedy"].set_index("task")["fitness"]
    for name in ("tabu", "tabu_cl"):  # both start from the greedy solution and keep the best
        part = runs[runs["algorithm"] == name].set_index("task")["fitness"]
        assert (part >= greedy.loc[part.index]).all(), name


def test_tuning_summary_is_consistent(summary: dict[str, Any]) -> None:
    assert set(summary["algorithms"]) == set(plan.ALGORITHMS_TUNED)
    for name, entry in summary["algorithms"].items():
        regrets = [t["mean_regret"] for t in entry["trials"]]
        assert entry["best_mean_regret"] == pytest.approx(min(regrets))
        assert entry["params"] == plan.trial_list(CFG, name)[entry["best_trial"]]
        assert len(regrets) == len(plan.trial_list(CFG, name))
        if entry["default_mean_regret"] is not None:
            assert entry["best_mean_regret"] <= entry["default_mean_regret"] + 1e-12
    alns = summary["alns_regret_by_policy"]
    assert summary["best_alns_policy"] == min(alns, key=lambda a: alns[a])
    assert summary["provenance_added"]["later"]["git"]["dirty"] is False  # tabu_cl, D-045
    assert summary["provenance"]["git"]["dirty"] is False


# ---- statistics ---------------------------------------------------------------------------------


def test_stored_tests_can_be_recomputed_from_the_raw_runs(
    runs: pd.DataFrame, tests_table: pd.DataFrame
) -> None:
    pivot = runs.pivot_table(index="task", columns="algorithm", values="fitness")
    best = runs.groupby("task")["best_known"].first().loc[pivot.index]
    rel = pivot.div(best, axis=0)
    for a, b in (("alns_ucb1", "alns_random"), ("alns_roulette", "alns_random")):
        row = _pooled(tests_table, a, b)
        d = rel[a] - rel[b]
        assert row["n"] == len(d)
        assert row["p"] == pytest.approx(stats.wilcoxon(d, zero_method="wilcox").pvalue)
        assert row["a12"] == pytest.approx(analysis.a12(rel[a].to_numpy(), rel[b].to_numpy()))
        assert row["median_diff_pp"] == pytest.approx(100 * np.median(d))
        assert row["wins"] + row["ties"] + row["losses"] == len(d)
    assert (tests_table["p_holm"] >= tests_table["p"] - 1e-15).all()
    assert (tests_table["p_holm"] <= 1.0).all()


def test_the_stored_claim_follows_from_the_stored_tests(tests_table: pd.DataFrame) -> None:
    info = json.loads((RESULTS / "analysis" / "analysis.json").read_text(encoding="utf-8"))
    assert info["claim_bandit_beats_random"] == json.loads(
        json.dumps(analysis.claim_verdict(tests_table, float(CFG["statistics"]["alpha"])))
    )
    assert info["read_sets"] == 300
    assert info["best_alns_policy_from_tuning"] == "alns_random"


# ---- findings reported in docs/phases/PHASE_5.md (regression tests) -----------------------------


def test_finding_ucb1_does_not_beat_random_and_epsilon_greedy_is_slightly_worse(
    tests_table: pd.DataFrame,
) -> None:
    ucb = _pooled(tests_table, "alns_ucb1", "alns_random")
    eps = _pooled(tests_table, "alns_epsilon_greedy", "alns_random")
    assert ucb["p_holm"] > 0.05 and 0.45 < ucb["a12"] < 0.55
    assert (
        eps["p_holm"] < 0.05 and eps["a12"] < 0.5
    )  # random is better (a small effect: A12 above 0.4)
    assert eps["a12"] > 0.4


def test_finding_roundrobin_and_roulette_are_slightly_better_than_random(
    tests_table: pd.DataFrame,
) -> None:
    for name in ("alns_roundrobin", "alns_roulette"):
        row = _pooled(tests_table, name, "alns_random")
        assert row["p_holm"] < 0.05 and 0.5 < row["a12"] < 0.6


def test_finding_lkh3_is_never_beaten_and_greedy_is_far_behind(
    runs: pd.DataFrame, tests_table: pd.DataFrame
) -> None:
    for name in ("alns_random", "alns_ucb1", "tabu_cl", "tabu"):
        pivot = runs.pivot_table(index="task", columns="algorithm", values="fitness")
        assert (pivot[name] <= pivot["lkh3"]).all(), name
    med = runs.groupby("algorithm")["gap"].median()
    assert med["lkh3"] == 0.0 and med["greedy"] > 0.3 and med["nearest_neighbour"] > 2.0
    alns = [med[a] for a in med.index if a.startswith("alns_") and "fixed" not in a]
    assert max(alns) < med["tabu_cl"] < med["tabu"] < med["greedy"]


def test_finding_candidate_lists_help_but_alns_is_still_better(tests_table: pd.DataFrame) -> None:
    better = _pooled(tests_table, "tabu_cl", "tabu")
    assert better["p_holm"] < 0.001 and better["wins"] > 5 * better["losses"]
    alns = _pooled(tests_table, "alns_random", "tabu_cl")
    assert alns["p_holm"] < 0.001 and alns["wins"] > alns["losses"]


def test_finding_a_mixture_of_operators_beats_the_random_and_greedy_single_operators(
    tests_table: pd.DataFrame,
) -> None:
    for name in ("alns_fixed_random+greedy", "alns_fixed_random+regret"):
        row = _pooled(tests_table, "alns_random", name)
        assert row["p_holm"] < 0.001 and row["median_diff_pp"] > 0


def test_finding_biological_quality_is_similar_for_all_solvers_and_has_no_misjoins(
    runs: pd.DataFrame,
) -> None:
    by = runs.groupby("algorithm")
    assert (by["bio_misjoins"].max() == 0).all()
    assert (by["bio_genome_fraction"].median() > 0.99).all()
    assert (by["bio_identity"].median() > 0.99).all()
    assert (
        by["bio_n_contigs"].median().min() > 20
    )  # the read errors fragment every solution (D-031)


# ---- scaling and figures ------------------------------------------------------------------------


def test_scaling_runs_are_complete_and_lkh3_total_time_grows_with_n() -> None:
    recs = runner.read_records(RESULTS / "scaling_runs.jsonl.gz")
    runs = [r for r in recs if r["kind"] == "run"]
    n_sets = len(CFG["scaling_instances"]) * int(CFG["seeds"]["scaling"])
    for name in CFG["scaling_solvers"]:
        assert sum(r["requested"] == name for r in runs) == n_sets, name
    table = pd.read_csv(RESULTS / "analysis" / "scaling.csv").sort_values("n")
    assert table["n"].is_monotonic_increasing and len(table) == len(CFG["scaling_instances"])
    assert table["lkh3_seconds"].is_monotonic_increasing
    assert table["build_seconds"].is_monotonic_increasing


@pytest.mark.parametrize(
    "stem",
    [
        "fig04_convergence",
        "fig05_final_fitness",
        "fig06_gap_to_lkh3",
        "fig07_ablation",
        "fig08_arm_probability",
        "fig09_operator_heatmap",
        "fig10_runtime_vs_n",
        "fig11_pareto",
    ],
)
def test_figures_exist_in_every_format_with_caption_and_matching_config_hash(stem: str) -> None:
    for suffix in (".svg", ".pdf", ".png"):
        assert (FIGURES / f"{stem}{suffix}").stat().st_size > 1000
    caption = json.loads((FIGURES / f"{stem}.caption.json").read_text(encoding="utf-8"))
    info = json.loads((RESULTS / "analysis" / "analysis.json").read_text(encoding="utf-8"))
    assert len(caption["caption"]) > 40
    assert caption["config_hash"].startswith(info["provenance"]["config_hash"][:12])
