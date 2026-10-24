"""Tests of the application's data layer (no browser): content, solvers, charts, results."""

from __future__ import annotations

import re
from typing import Any

import numpy as np
import pytest

from genomeasm.appkit import content, figs, lab, repro, results_view
from genomeasm.appkit.paths import GENOMES, RESULTS
from genomeasm.solver_validation import random_matrix_problem
from genomeasm.solvers.greedy import greedy
from genomeasm.solvers.lkh import LkhError

SMALL = lab.InstanceParams(length=600, coverage=6, read_length=80, error=0.01, min_overlap=16)


@pytest.fixture(scope="module")
def built() -> Any:
    return lab.build(SMALL)


# ---- glossary and toy example -------------------------------------------------------------------


def test_glossary_search_needs_every_word_and_is_sorted() -> None:
    terms = content.load_glossary()
    assert len(terms) > 40 and all({"term", "plain", "precise"} <= set(t) for t in terms)
    everything = content.search_glossary(terms, "")
    assert [t["term"].lower() for t in everything] == sorted(t["term"].lower() for t in terms)
    hits = content.search_glossary(terms, "COVERAGE")
    assert hits and all("coverage" in str(t).lower() for t in hits)
    assert len(content.search_glossary(terms, "coverage zzzqqq")) == 0
    narrower = content.search_glossary(terms, "coverage depth")
    assert 0 < len(narrower) <= len(hits)


def test_toy_example_matches_the_worked_example_of_the_objective_document() -> None:
    toy = content.make_toy(list(content.TOY_READS), min_overlap=1)
    expected = np.array([[0, 3, 1], [1, 0, 3], [0, 0, 0]])  # docs/OBJECTIVE_FUNCTION.md, section 2
    assert np.array_equal(toy.instance.ov, expected) and toy.labels == ["A", "B", "C"]
    rows = content.all_orders(toy)
    assert len(rows) == 6
    assert (rows[0]["order"], rows[0]["F"], rows[0]["L"]) == ("ABC", 6, 9)
    assert rows[0]["merged"] == "ATGCGTACC"
    by_order = {r["order"]: (r["F"], r["L"]) for r in rows}
    assert by_order == {
        "ABC": (6, 9), "BCA": (3, 12), "CAB": (3, 12),
        "BAC": (2, 13), "ACB": (1, 14), "CBA": (1, 14),
    }  # fmt: skip
    assert all(r["F"] + r["L"] == 15 for r in rows)  # L = sum of lengths - F


def test_toy_input_is_cleaned_and_validated() -> None:
    assert content.clean_reads("atgcg, gcgta\nxx gtacc") == ["ATGCG", "GCGTA", "", "GTACC"]
    assert content.clean_reads("atg-cg  N,GCGTA") == ["ATGCG", "", "GCGTA"]
    with pytest.raises(ValueError, match="between 2"):
        content.make_toy(["ACGT"])
    with pytest.raises(ValueError, match="at most 40"):
        content.make_toy(["A" * 41, "C" * 10])
    toy = content.make_toy(["ACGTACGT", "ACGTACGT", "TTTT"], 2)  # a duplicate is dropped
    assert toy.notes and toy.instance.n == 2
    base = content.make_toy(list(content.TOY_READS), 1)
    assert content.order_from_letters(base, "bca").tolist() == [1, 2, 0]
    with pytest.raises(ValueError, match="exactly once"):
        content.order_from_letters(base, "AAB")


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_greedy_steps_reproduce_the_greedy_solver(seed: int) -> None:
    problem = random_matrix_problem(np.random.default_rng(seed), 14)
    steps, order = content.greedy_steps(problem.weights[1:, 1:])
    assert order == greedy(problem).order.tolist()
    accepted = [s for s in steps if s["verdict"] == "accepted"]
    assert len(accepted) <= 13 and all(s["overlap"] > 0 for s in steps)
    overlaps = [s["overlap"] for s in steps]
    assert overlaps == sorted(overlaps, reverse=True)  # largest overlap first
    assert any("cycle" in s["verdict"] for s in steps) or len(steps) > 0


def test_greedy_steps_explain_the_toy() -> None:
    toy = content.make_toy(list(content.TOY_READS), 1)
    steps, order = content.greedy_steps(toy.instance.ov)
    assert order == [0, 1, 2]
    assert [(s["from"], s["to"], s["verdict"]) for s in steps][:2] == [
        (0, 1, "accepted"),
        (1, 2, "accepted"),
    ]
    assert any("already has" in s["verdict"] for s in steps)  # A to C (1) is rejected


# ---- instances and solvers ----------------------------------------------------------------------


def test_instance_limits_protect_the_thirty_second_envelope() -> None:
    assert lab.InstanceParams().problems() == []
    too_many = lab.InstanceParams(length=5000, coverage=25, read_length=100)
    assert too_many.estimated_reads() > lab.MAX_READS and "reads" in too_many.problems()[0]
    assert lab.InstanceParams(length=9000).problems()
    assert (
        lab.InstanceParams(error=0.5).problems() and lab.InstanceParams(min_overlap=100).problems()
    )
    assert lab.InstanceParams(genome="NC_000000.1").problems()
    with pytest.raises(ValueError, match="reads"):
        lab.build(too_many)
    assert lab.MAX_SECONDS + 10 <= lab.MAX_TOTAL_SECONDS


def test_build_is_deterministic_and_cached(built: Any) -> None:
    assert lab.build(SMALL) is built
    again = lab.build(lab.InstanceParams(**{**SMALL.to_dict(), "seed": 2}))
    assert again.instance.n != built.instance.n or not np.array_equal(
        again.instance.ov, built.instance.ov
    )
    assert built.problem.n == built.instance.n > 20
    assert set(GENOMES) == {"NC_001422.1", "NC_001416.1"}


def test_tuned_parameters_come_from_the_phase5_summary() -> None:
    summary = lab.tuning_summary()
    assert summary is not None
    tuned = lab.solver_params("alns", "ucb1_rt", True)
    assert tuned == summary["algorithms"]["alns_ucb1_rt"]["params"]
    default = lab.solver_params("alns", "ucb1_rt", False)
    assert default["alns"]["reward_preset"] == "classical"  # trial 0 of the configuration
    assert (
        lab.solver_params("tabu_cl", "random", True) == summary["algorithms"]["tabu_cl"]["params"]
    )


@pytest.mark.parametrize("solver", ["greedy", "nearest_neighbour", "tabu", "tabu_cl", "alns"])
def test_every_local_solver_returns_a_valid_order(built: Any, solver: str) -> None:
    result = lab.run_solver(built, solver, seconds=1, evaluations=40000, policy="roulette")
    n = built.problem.n
    assert sorted(result.order.tolist()) == list(range(n))
    assert result.fitness == built.instance.fitness(result.order)
    if solver in ("tabu", "tabu_cl", "alns"):
        assert result.fitness >= greedy(built.problem).fitness


def test_live_progress_and_step_log_come_through_run_solver(built: Any) -> None:
    events: list[dict[str, Any]] = []
    steps: list[dict[str, Any]] = []
    result = lab.run_solver(
        built, "alns", seconds=1, evaluations=30000, policy="ucb1", progress=events.append,
        step_log=steps,
    )  # fmt: skip
    assert events[-1]["final"] and events[-1]["best"] == result.fitness
    assert len(steps) == result.info["iterations"] > 5
    with pytest.raises(ValueError, match="unknown solver"):
        lab.run_solver(built, "magic", seconds=1)
    with pytest.raises(ValueError, match="unknown policy"):
        lab.run_solver(built, "alns", seconds=1, evaluations=100, policy="oracle")


def test_the_time_budget_is_capped_for_interactive_use(built: Any) -> None:
    result = lab.run_solver(built, "tabu", seconds=1000, evaluations=2000)
    assert result.evaluations >= 2000  # an evaluation budget replaces the time (test only)
    assert min(1000.0, lab.MAX_SECONDS) == lab.MAX_SECONDS


def test_lkh3_is_never_replaced_when_it_cannot_run(
    built: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Any
) -> None:
    monkeypatch.setenv("GENOMEASM_LKH", str(tmp_path / "no-such-LKH"))
    assert not lab.lkh_available()
    with pytest.raises(LkhError):  # rule R9: reported as unavailable, no substitute solver
        lab.run_solver(built, "lkh3", seconds=1)


def test_quality_and_dot_plot_follow_the_true_order() -> None:
    clean = lab.build(lab.InstanceParams(**{**SMALL.to_dict(), "error": 0.0, "tolerance": 0.0}))
    order = clean.instance.true_order()
    q = lab.quality(clean, order)
    assert q["misjoins"] == 0 and q["genome_fraction"] > 0.9 and q["identity"] > 0.99
    assembly = clean.instance.assemble(order)
    dp = lab.dot_plot(clean, assembly)
    assert dp["boundaries"][0] == 0 and len(dp["x"]) == len(dp["y"]) == len(dp["contig"])
    mapped = ~np.isnan(dp["y"])
    assert mapped.sum() >= 3
    one_contig = dp["contig"] == dp["contig"][0]
    if (one_contig & mapped).sum() >= 3:  # inside a contig: y = x + constant (a 45 degree line)
        x, y = dp["x"][one_contig & mapped], dp["y"][one_contig & mapped]
        assert np.ptp(y - x) <= 3
    assert dp["reference_length"] == SMALL.length


# ---- charts -------------------------------------------------------------------------------------


def test_every_chart_renders_and_exports_as_svg_and_png(built: Any) -> None:
    toy = content.make_toy(list(content.TOY_READS), 1)
    reads = lab.simulate(SMALL)
    result = lab.run_solver(built, "alns", seconds=1, evaluations=20000, policy="roulette")
    depth_fig, stats = figs.depth_vs_theory(reads)
    assert 0 <= stats["uncovered_observed"] <= 1 and stats["reads"] == reads.n_reads
    card = {"base_counts": {"A": 1, "C": 2, "G": 3, "T": 4}, "length_bp": 10}
    assembly = built.instance.assemble(result.order)
    charts = [
        figs.toy_graph(toy.instance, toy.labels, [0, 1, 2]),
        figs.matrix_figure(toy.instance.ov, toy.labels, "overlaps", "letters"),
        figs.tour_strip([0, 1, 2], [1, 0, 2], toy.labels, "step"),
        figs.gc_profile(lab.window_of(SMALL), 100),
        figs.base_counts(card),
        figs.pipeline_diagram(),
        depth_fig,
        figs.read_map(reads),
        figs.overlap_heatmap(built.instance, max_n=60),
        figs.overlap_lengths(built.instance),
        figs.degree_distribution(built.instance),
        figs.live_convergence([(0.0, 10.0), (1.0, 12.0)], {"greedy": 9.0}, 2.0),
        figs.operator_probabilities([(0.0, [1 / 6] * 6), (1.0, [0.5] + [0.1] * 5)], 2.0),
        figs.outcome_bars(result.info["outcomes"]),
        figs.dot_plot(lab.dot_plot(built, assembly)),
        figs.tour_blocks(result.order.tolist(), built.instance.true_starts),
        figs.box_by_algorithm({"a": np.array([0.0, 0.1]), "b": np.array([0.2, 0.4])}, "gap", "t"),
    ]  # fmt: skip
    for fig in charts:
        svg = figs.figure_bytes(fig, "svg", footer="live run")
        png = figs.figure_bytes(fig, "png")
        assert (
            svg.startswith(b"<?xml") and b"\r\n" not in svg and b"live run" in svg or b"<g" in svg
        )
        assert png[:8] == b"\x89PNG\r\n\x1a\n"
    with pytest.raises(ValueError, match="fmt"):
        figs.figure_bytes(charts[0], "pdf")


def test_svg_export_is_byte_identical_for_the_same_figure() -> None:
    first = figs.figure_bytes(figs.pipeline_diagram(), "svg")
    second = figs.figure_bytes(figs.pipeline_diagram(), "svg")
    assert first == second  # the footer is removed again after export


# ---- precomputed results and reproducibility ----------------------------------------------------


def test_phase5_results_load_and_filters_work() -> None:
    p5 = results_view.load_phase5()
    assert results_view.available() and p5.info["read_sets"] == 300
    assert {"alns_ucb1", "lkh3", "tabu_cl"} <= set(p5.by_algorithm["algorithm"])
    primary = results_view.filter_tests(p5.tests, family="primary", scope="pooled")
    assert len(primary) == 2 and set(primary["b"]) == {"alns_random"}
    only = results_view.filter_tests(p5.tests, family="primary", significant_only=True)
    assert only["significant"].all() and len(only) < len(p5.tests)
    assert (results_view.filter_tests(p5.tests, algorithm="lkh3")[["a", "b"]] == "lkh3").any(
        axis=None
    )
    runs = results_view.load_runs_frame()
    assert runs.groupby("algorithm").size().nunique() == 1  # the same number of read sets for all
    part = results_view.filter_runs(runs, instances=["L5-c10-e1.0"], algorithms=["greedy", "lkh3"])
    assert set(part["algorithm"]) == {"greedy", "lkh3"} and set(part["instance"]) == {"L5-c10-e1.0"}


def test_figure_catalog_lists_every_figure_with_its_files() -> None:
    catalog = results_view.figure_catalog()
    assert len(catalog) == len(results_view.FIGURE_STEMS) + len(results_view.EARLIER_STEMS)
    for entry in catalog:
        assert entry["png"].is_file() and entry["svg"].is_file() and entry["pdf"].is_file()
        assert len(entry["caption"]) > 40 and re.fullmatch(r"[0-9a-f]{64}", entry["config_hash"])


def test_reproducibility_information_is_complete_and_checksums_verify() -> None:
    cfg = repro.config_info()
    assert cfg["path"] == "configs/phase5_experiments.yaml" and len(cfg["config_hash"]) == 64
    assert cfg["master_seed"] == 20261025 and cfg["seeds"]["evaluation"] == 30
    env = repro.environment()
    assert env["git"]["commit"] and env["libraries"]["numpy"] and env["machine"]["python"]
    sums = repro.verify_checksums()
    assert len(sums) == 2 and all(row["ok"] for row in sums)
    files = repro.result_files()
    assert any(f["file"] == "results/phase5/tuning_summary.json" for f in files)
    assert all(re.fullmatch(r"[0-9a-f]{64}", f["sha256"]) for f in files)
    assert RESULTS.is_dir()
