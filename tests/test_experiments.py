"""Phase 5 machinery: biological metrics, statistics, tuning plan, tasks and one tiny real task."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

from genomeasm.assembly import Assembly, Contig
from genomeasm.experiments import analysis, design, plan, runner
from genomeasm.experiments.bio import ReferenceIndex, assembly_quality, n50, order_quality

ROOT = Path(__file__).resolve().parent.parent
CFG = design.load_config(ROOT / "configs" / "phase5_experiments.yaml")
REF = "".join(np.random.default_rng(3).choice(list("ACGT"), size=1500))


def contig(seq: str) -> Contig:
    return Contig(sequence=seq, reads=[0], offsets=[0])


# ---- biological metrics ------------------------------------------------------------------------


def test_n50_hand_examples() -> None:
    assert n50([]) == 0 and n50([10]) == 10
    assert n50([2, 2, 2, 3, 3, 4, 8]) == 4  # sorted 8,4,3,3,2,2,2: 8+4 = 12 >= 12 (half of 24)
    assert n50([100, 100]) == 100


def test_reference_index_locates_exact_and_mutated_windows() -> None:
    index = ReferenceIndex(REF)
    window = REF[400:500]
    assert index.locate(window) == (400, 100, 100)
    mutated = window[:50] + ("A" if window[50] != "A" else "C") + window[51:]
    diagonal, matches, aligned = index.locate(mutated)  # type: ignore[misc]
    assert (diagonal, matches, aligned) == (400, 99, 100)
    assert index.locate("ACGT" * 25) is None or index.locate("ACGT" * 25)[1] < 100  # type: ignore[index]


def test_perfect_contig_has_full_identity_and_a_misjoin_is_detected() -> None:
    index = ReferenceIndex(REF)
    good = assembly_quality(Assembly([contig(REF[100:700])]), index)
    assert good["identity"] == 1.0 and good["misjoins"] == 0 and good["unmapped_fraction"] == 0
    assert good["genome_fraction"] == pytest.approx(600 / len(REF))
    assert good["n_contigs"] == 1 and good["n50"] == 600 and good["assembled_length"] == 600
    jumpy = assembly_quality(Assembly([contig(REF[100:300] + REF[900:1100])]), index)
    assert jumpy["misjoins"] == 1  # the contig jumps from position 300 to 900
    junk = assembly_quality(
        Assembly([contig("".join(np.random.default_rng(9).choice(list("ACGT"), 300)))]), index
    )
    assert junk["genome_fraction"] == 0.0 and junk["unmapped_fraction"] > 0.5


def test_order_quality_of_perfect_reversed_and_shuffled_orders() -> None:
    starts = np.arange(0, 1000, 10, dtype=np.int64)
    order = np.arange(starts.size)
    perfect = order_quality(order, starts, 100, [list(order)])
    assert perfect["spearman"] == pytest.approx(1.0)
    assert perfect["spearman_contigs"] == pytest.approx(1.0)
    assert perfect["adjacent_correct"] == 1.0
    reversed_ = order_quality(order[::-1], starts, 100, [list(order[::-1])])
    assert reversed_["spearman"] == pytest.approx(-1.0)
    halves = [list(order[50:]), list(order[:50])]  # two correct contigs placed in the wrong order
    swapped = order_quality(np.concatenate(halves), starts, 100, halves)
    assert swapped["spearman"] < 0 and swapped["spearman_contigs"] == pytest.approx(1.0)


# ---- statistics ---------------------------------------------------------------------------------


def test_holm_and_effect_sizes() -> None:
    assert analysis.holm([0.01, 0.04, 0.03]) == pytest.approx([0.03, 0.06, 0.06])
    assert analysis.holm([0.5]) == [0.5] and analysis.holm([0.0001, 0.9])[0] == pytest.approx(
        0.0002
    )
    assert analysis.a12(np.array([2.0, 3.0]), np.array([1.0, 2.0])) == pytest.approx(0.875)
    assert analysis.a12(np.array([1.0]), np.array([1.0])) == 0.5


def test_paired_test_detects_a_shift_and_not_a_tie() -> None:
    rng = np.random.default_rng(1)
    base = 0.9 + 0.05 * rng.random(40)
    shifted = analysis.paired_test(base + 0.01, base, 2000, 0)
    assert shifted["p"] < 1e-6 and shifted["median_diff_pp"] == pytest.approx(1.0)
    assert shifted["ci_low_pp"] == pytest.approx(1.0) == pytest.approx(shifted["ci_high_pp"])
    assert shifted["wins"] == 40
    tie = analysis.paired_test(base, base, 2000, 0)
    assert tie["p"] == 1.0 and tie["ties"] == 40 and tie["a12"] == 0.5


def _fake_records(n_sets: int = 8) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for s in range(n_sets):
        records.append(
            {"kind": "instance", "stage": "main", "instance": "X", "seed_index": s, "n": 10,
             "n_input": 11, "edges": 5, "build_seconds": 0.1, "F_true": 90,
             "bio_true": {"identity": 1.0}}
        )  # fmt: skip
        for name, value in (
            ("alns_random", 100),
            ("alns_ucb1", 101),
            ("lkh3", 102),
            ("greedy", 95),
        ):
            records.append(
                {"kind": "run", "stage": "main", "instance": "X", "seed_index": s, "n": 10,
                 "fitness": value + (s % 2), "evaluations": 5, "seconds": 1.0, "merged_length": 1,
                 "requested": name, "bio": {"identity": 1.0}, "info": {},
                 "trace": [[0.0, 1, 90], [0.5, 3, 95]]}
            )  # fmt: skip
    return records


def test_load_runs_and_gaps() -> None:
    runs, inst = analysis.load_runs(_fake_records())
    assert len(runs) == 32 and len(inst) == 8
    row = runs[(runs["algorithm"] == "alns_ucb1") & (runs["seed_index"] == 0)].iloc[0]
    assert row["best_known"] == 102 and row["gap"] == pytest.approx(100 * 1 / 102)
    assert row["gap_lkh"] == pytest.approx(100 * 1 / 102)
    assert runs[runs["algorithm"] == "lkh3"]["gap"].eq(0).all()


def test_tests_are_holm_corrected_within_families_and_the_verdict_follows_the_rules() -> None:
    runs, _ = analysis.load_runs(_fake_records())
    cfg = {**CFG, "statistics": {**CFG["statistics"], "bootstrap_resamples": 500,
           "primary": [["alns_ucb1", "alns_random"]], "secondary_vs_random": [], "vs_roulette": [],
           "single_operator_vs_random": False,
           "main_pairs": [["alns_best", "greedy"], ["lkh3", "alns_best"]]}}  # fmt: skip
    tests = analysis.run_tests(runs, cfg, "alns_ucb1")
    pooled = tests[tests["scope"] == "pooled"]
    main = pooled[pooled["family"] == "main"]
    assert set(zip(main["a"], main["b"], strict=True)) == {
        ("alns_ucb1", "greedy"),
        ("lkh3", "alns_ucb1"),
    }
    assert (main["p_holm"] >= main["p"]).all() and main["p_holm"].le(1).all()
    verdict = analysis.claim_verdict(tests, 0.05)
    assert verdict[0]["policy"] == "alns_ucb1" and verdict[0]["median_diff_pp"] > 0
    assert verdict[0]["verdict"] == "supported"
    flipped = tests.copy()
    flipped["median_diff_pp"] = -1.0
    assert "rejected" in analysis.claim_verdict(flipped, 0.05)[0]["verdict"]
    flipped["significant"] = False
    assert "not supported" in analysis.claim_verdict(flipped, 0.05)[0]["verdict"]


def test_convergence_steps_and_time_to_target() -> None:
    trace = [[0.1, 1, 10], [0.5, 4, 12], [2.0, 9, 15]]
    grid = np.array([0.05, 0.1, 0.4, 0.5, 3.0])
    out = analysis._step(trace, 0, grid)
    assert np.isnan(out[0]) and out[1:].tolist() == [10, 10, 12, 15]
    assert analysis.time_to_target(trace, 12) == 0.5 and analysis.time_to_target(trace, 99) is None
    assert analysis.time_to_target(trace, 15, offset=1.5) == 3.5


# ---- planning -----------------------------------------------------------------------------------


def test_trial_zero_is_the_default_and_sampling_is_deterministic() -> None:
    for algorithm in plan.TUNED:
        first, again = plan.trial_list(CFG, algorithm), plan.trial_list(CFG, algorithm)
        assert first == again and len(first) == CFG["tuning"]["trials"]
        assert first[0] == plan.default_params(CFG, algorithm)
        assert len({json.dumps(t, sort_keys=True) for t in first}) > 15  # real variety
    space = CFG["tuning"]["space"]["alns"]
    for trial in plan.trial_list(CFG, "alns_ucb1")[1:]:
        assert all(trial["alns"][k] in space[k] for k in space)
    assert len(plan.trial_list(CFG, "lkh3")) == len(CFG["tuning"]["lkh_variants"])
    assert plan.trial_list(CFG, "alns_random") != plan.trial_list(CFG, "alns_ucb1")


def test_seeds_are_independent_of_task_order_and_differ_between_stages() -> None:
    a = design.seeds_for(CFG, "main", 3, 7)
    assert a == design.seeds_for(CFG, "main", 3, 7)
    others = {design.seeds_for(CFG, "main", 3, 8), design.seeds_for(CFG, "main", 4, 7),
              design.seeds_for(CFG, "tuning", 3, 7),
              design.seeds_for(CFG, "scaling", 3, 7)}  # fmt: skip
    assert a not in others and len(others) == 4


def test_tuning_selection_picks_the_smallest_mean_regret() -> None:
    cfg = {**CFG, "tuning": {**CFG["tuning"], "trials": 3}}
    records = []
    for algorithm in plan.ALGORITHMS_TUNED:
        for trial in range(len(plan.trial_list(cfg, algorithm))):
            for s in range(2):
                bonus = {"alns_roulette": 5, "alns_ucb1": 3}.get(algorithm, 0)
                fitness = 1000 - (trial != 1) * 10 + bonus
                records.append({"kind": "run", "instance": "T", "seed_index": s,
                                "requested": algorithm, "trial": trial,
                                "fitness": fitness})  # fmt: skip
    summary = plan.select_tuned(cfg, records)
    assert all(summary["algorithms"][a]["best_trial"] == 1 for a in plan.ALGORITHMS_TUNED)
    assert summary["best_alns_policy"] == "alns_roulette"
    assert plan.tuned_params(summary, "alns_fixed_random+greedy")["policy"] == {}


def test_main_tasks_have_the_designed_shape() -> None:
    summary = plan.select_tuned(CFG, _tuning_records())
    tasks = plan.main_tasks(CFG, summary, seeds=2)
    assert len(tasks) == 2 * len(CFG["evaluation_instances"])
    names = [r["algorithm"] for r in tasks[0]["runs"]]
    assert sorted(names) == sorted(plan.MAIN_ALGORITHMS) and names[:2] == [
        "greedy",
        "nearest_neighbour",
    ]
    assert names != [r["algorithm"] for r in tasks[1]["runs"]]  # the timed order is shuffled
    assert len({runner.task_key(t) for t in tasks}) == len(tasks)
    assert plan.main_tasks(CFG, summary, seeds=2)[0]["runs"] == tasks[0]["runs"]  # reproducible
    scaling = plan.scaling_tasks(CFG, summary, seeds=1)
    assert len(scaling) == 4 and scaling[0]["lkh_search_seconds"] == CFG["lkh_scaling_seconds"]
    assert len(plan.tuning_tasks(CFG, seeds=3)) == 3 * 6 * len(plan.ALGORITHMS_TUNED)
    tuning = plan.tuning_tasks(CFG, seeds=1, trials=2)[0]
    assert (
        tuning["spec"]["genome"] == "NC_001422.1"
        and tuning["spec"]["start"] + tuning["spec"]["length"] <= 2500
    )


def _tuning_records() -> list[dict[str, Any]]:
    cfg = CFG
    return [
        {
            "kind": "run",
            "instance": "T",
            "seed_index": 0,
            "requested": a,
            "trial": t,
            "fitness": 100 - t,
        }
        for a in plan.ALGORITHMS_TUNED
        for t in range(len(plan.trial_list(cfg, a)))
    ]


def test_evaluation_instances_do_not_overlap_the_tuning_region_or_seeds() -> None:
    tuning = {(s["genome"], s["start"], s["start"] + s["length"]) for s in CFG["tuning_instances"]}
    assert tuning == {("NC_001422.1", 0, 2500)}
    for s in CFG["evaluation_instances"]:
        if s["genome"] == "NC_001422.1":
            assert s["start"] >= 2500  # phiX174 evaluation windows start after the tuning window
    assert {s["genome"] for s in CFG["evaluation_instances"]} == {"NC_001416.1", "NC_001422.1"}
    assert all(seed >= 1000 for seed in range(1000, 1003)) and CFG["seeds"]["evaluation"] == 30


# ---- parameters and one tiny real task ----------------------------------------------------------


def test_parameter_conversion() -> None:
    kw = design.tabu_kwargs(
        {"tenure_factor": 2.0, "moves_factor": 10, "restart_factor": 5, "kick_fraction": 0.1}, 400
    )
    assert kw == {"tenure": 40, "moves_per_iteration": 4000, "restart_after": 2000, "kick": 40}
    cl = design.tabu_kwargs(plan.default_params(CFG, "tabu_cl")["tabu"], 400)
    assert cl["candidate_fraction"] == 0.5 and cl["candidates"] == 6
    assert "candidate_fraction" not in design.tabu_kwargs(
        plan.default_params(CFG, "tabu")["tabu"], 9
    )
    params = plan.default_params(CFG, "alns_ucb1")
    akw = design.alns_kwargs(CFG, params)
    assert akw["destroy_fraction"] == (0.10, 0.40) and akw["rewards"] == (1.0, 0.2727, 0.3939, 0.0)


def test_one_tiny_real_task_runs_end_to_end() -> None:
    summary = plan.select_tuned(CFG, _tuning_records())
    spec = CFG["evaluation_instances"][0]
    names = ["greedy", "tabu", "tabu_cl", "alns_ucb1", "alns_fixed_random+greedy"]
    task = plan.main_tasks(CFG, summary, seeds=1, specs=[spec], algorithms=names)[0]
    task["seconds"] = 0.4
    records = runner.run_task(task)
    kinds = [r["kind"] for r in records]
    assert kinds.count("instance") == 1 and kinds.count("run") == 5 and "error" not in kinds
    instance = records[0]
    assert instance["n"] > 100 and instance["bio_true"]["misjoins"] == 0
    runs = {r["requested"]: r for r in records if r["kind"] == "run"}
    assert runs["greedy"]["seconds"] < 0.4 and runs["tabu"]["fitness"] >= runs["greedy"]["fitness"]
    assert runs["tabu_cl"]["fitness"] >= runs["greedy"]["fitness"]
    assert (
        runs["alns_ucb1"]["info"]["outcomes"] and runs["alns_ucb1"]["bio"]["genome_fraction"] > 0.9
    )
    assert runs["alns_ucb1"]["trace"][0][0] <= runs["alns_ucb1"]["trace"][-1][0]


def test_finished_tasks_are_not_run_again(tmp_path: Path) -> None:
    summary = plan.select_tuned(CFG, _tuning_records())
    tasks = plan.main_tasks(CFG, summary, seeds=1, specs=CFG["evaluation_instances"][:2])
    work = tmp_path / "work.jsonl"
    work.write_text(
        "".join(json.dumps({"kind": "task_done", "key": runner.task_key(t)}) + "\n" for t in tasks)
    )
    said: list[str] = []
    runner.execute(tasks, work, workers=1, say=said.append)
    assert "0 to run" in said[0]
    runner.compress(work, tmp_path / "out.jsonl.gz")
    assert runner.read_records(tmp_path / "out.jsonl.gz") == []  # bookkeeping lines are hidden


def test_markdown_table_helper() -> None:
    text = analysis.markdown(pd.DataFrame({"a": [1, 2], "b": [0.5, 1.25]}))
    assert text.splitlines()[0] == "| a | b |" and "| 2 | 1.250 |" in text


# ---- algorithms added after the first run (tabu_cl, decision D-045) ----------------------------


def test_adding_tabu_cl_does_not_change_the_other_tuning_trials() -> None:
    assert (
        plan.ALGORITHMS_TUNED[-1] == "tabu_cl"
    )  # appended: earlier random-search seeds are intact
    for algorithm in plan.ALGORITHMS_TUNED[:-1]:
        assert plan.ALGORITHMS_TUNED.index(algorithm) < plan.ALGORITHMS_TUNED.index("tabu_cl")
    space = CFG["tuning"]["space"]["tabu_cl"]
    trials = plan.trial_list(CFG, "tabu_cl")
    assert trials[0] == plan.default_params(CFG, "tabu_cl")
    for trial in trials[1:]:
        assert all(trial["tabu"][k] in space[k] for k in space)
    assert plan.trial_list(CFG, "tabu") != trials


def test_select_tuned_keeps_the_entries_that_were_frozen_before() -> None:
    cfg = {**CFG, "tuning": {**CFG["tuning"], "trials": 3}}

    def records(best: int) -> list[dict[str, Any]]:
        return [
            {"kind": "run", "instance": "T", "seed_index": 0, "requested": a, "trial": t,
             "fitness": 100 - 10 * (t != best)}
            for a in plan.ALGORITHMS_TUNED
            for t in range(3)
        ]  # fmt: skip

    first = plan.select_tuned(cfg, records(best=1))
    assert (
        plan.select_tuned(cfg, records(best=2), frozen=first)["algorithms"] == first["algorithms"]
    )
    only_old = {"algorithms": {k: v for k, v in first["algorithms"].items() if k != "tabu_cl"}}
    mixed = plan.select_tuned(cfg, records(best=2), frozen=only_old)
    assert mixed["algorithms"]["tabu_cl"]["best_trial"] == 2
    assert mixed["algorithms"]["alns_ucb1"] == first["algorithms"]["alns_ucb1"]
    assert mixed["algorithms"]["tabu"]["best_trial"] == 1


def test_missing_tasks_cover_only_the_read_sets_without_a_result() -> None:
    summary = plan.select_tuned(CFG, _tuning_records())
    n_sets = len(CFG["evaluation_instances"]) * CFG["seeds"]["evaluation"]
    none = plan.missing_tasks(CFG, summary, [], "tabu_cl")
    assert len(none) == n_sets and {t["group"] for t in none} == {"add_tabu_cl"}
    assert all(len(t["runs"]) == 1 and t["runs"][0]["algorithm"] == "tabu_cl" for t in none)
    assert all(t["emit_instance"] is False and "lkh_search_seconds" not in t for t in none)
    spec = CFG["evaluation_instances"][0]["id"]
    have = [{"kind": "run", "requested": "tabu_cl", "stage": "main", "instance": spec,
             "seed_index": s} for s in range(CFG["seeds"]["evaluation"])]  # fmt: skip
    rest = plan.missing_tasks(CFG, summary, have, "tabu_cl")
    assert len(rest) == n_sets - CFG["seeds"]["evaluation"]
    assert not {t["spec"]["id"] for t in rest} & {spec}
    other_stage = [{**r, "stage": "scaling"} for r in have]
    assert len(plan.missing_tasks(CFG, summary, other_stage, "tabu_cl")) == n_sets
    scaling = plan.missing_tasks(CFG, summary, [], "lkh3", stage="scaling")
    assert {t["lkh_search_seconds"] for t in scaling} == {CFG["lkh_scaling_seconds"]}
    assert {t["group"] for t in scaling} == {"lkh_rerun"}
    assert len(plan.missing_lkh_tasks(CFG, summary, [])) == n_sets
    assert "tabu_cl" in CFG["scaling_solvers"]
