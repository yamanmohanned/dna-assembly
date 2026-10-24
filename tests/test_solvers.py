"""Solvers: delta evaluation, exact validators, baselines, tabu, ALNS and LKH-3."""

from __future__ import annotations

import itertools
from pathlib import Path

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from genomeasm import objective
from genomeasm.instance import AssemblyInstance
from genomeasm.io import Genome, load_genome
from genomeasm.objective import is_permutation
from genomeasm.simulate import ReadLengthModel, SimulationConfig, simulate_reads
from genomeasm.solver_validation import random_matrix_problem
from genomeasm.solvers.alns import ARM_NAMES, alns
from genomeasm.solvers.common import Budget, Problem, Recorder, order_of, tour_of
from genomeasm.solvers.exact import brute_force, held_karp
from genomeasm.solvers.greedy import greedy, nearest_neighbour
from genomeasm.solvers.lkh import (
    LkhError,
    lkh,
    lkh_binary,
    parse_tour_file,
    to_wsl_path,
    verify_tour,
)
from genomeasm.solvers.moves import (
    apply_insert,
    apply_swap,
    insert_delta,
    swap_delta,
    tour_fitness,
)
from genomeasm.solvers.tabu import _candidate_lists, tabu_search

ROOT = Path(__file__).resolve().parent.parent
TOY = ["ATGCG", "GCGTA", "GTACC"]
POLICIES = ["random", "roundrobin", "roulette", "epsilon_greedy", "ucb1", "ucb1_rt"]


def toy_problem() -> Problem:
    return Problem.from_instance(AssemblyInstance.from_sequences(TOY, min_overlap=1))


def valid(result_order: np.ndarray, n: int) -> bool:
    return is_permutation(result_order, n)


# ---- problem and budget -----------------------------------------------------------------------


def test_problem_agrees_with_the_objective_module() -> None:
    inst = AssemblyInstance.from_sequences(TOY, min_overlap=1)
    problem = Problem.from_instance(inst)
    assert problem.n == 3 and problem.total_length == 15
    assert np.array_equal(problem.cost_matrix(), inst.cost())
    for order in itertools.permutations(range(3)):
        assert problem.fitness(list(order)) == inst.fitness(order)
    assert np.all(problem.weights[0] == 0) and np.all(problem.weights[:, 0] == 0)


def test_budget_needs_a_limit_and_reports_use() -> None:
    with pytest.raises(ValueError):
        Budget()
    b = Budget(seconds=10, evaluations=100)
    assert not b.used_up(1.0, 50) and b.used_up(11.0, 50) and b.used_up(1.0, 100)
    assert b.fraction(5.0, 20) == 0.5 and b.fraction(1.0, 100) == 1.0 and b.fraction(99, 1) == 1.0


def test_recorder_trace_only_records_improvements() -> None:
    recorder = Recorder(toy_problem())
    for value, evals in [(3, 1), (3, 2), (5, 3), (4, 4), (6, 5)]:
        recorder.note(value, evals)
    assert [v for _, _, v in recorder.trace] == [3, 5, 6]


# ---- delta evaluation: every move against a full recomputation ----------------------------------


def test_every_swap_and_insert_delta_equals_a_full_recomputation_exhaustively() -> None:
    rng = np.random.default_rng(1)
    for n in (4, 6):
        problem = random_matrix_problem(rng, n)
        w = problem.weights
        t = tour_of(rng.permutation(n).astype(np.int64))
        for i, j in itertools.combinations(range(1, n + 1), 2):
            moved = t.copy()
            apply_swap(moved, i, j)
            assert swap_delta(t, w, i, j) == tour_fitness(moved, w) - tour_fitness(t, w)
            assert is_permutation(moved[1:] - 1, n) and moved[0] == 0
        for i in range(1, n + 1):
            for j in range(0, n + 1):
                if j in (i, i - 1):
                    continue
                moved = t.copy()
                apply_insert(moved, i, j)
                assert insert_delta(t, w, i, j) == tour_fitness(moved, w) - tour_fitness(t, w)
                assert is_permutation(moved[1:] - 1, n) and moved[0] == 0
                # the node really ended up right after the old neighbour t[j]
                position = int(np.flatnonzero(moved == t[i])[0])
                assert moved[position - 1] == t[j]


@settings(deadline=None, max_examples=60)
@given(seed=st.integers(min_value=0, max_value=10**6), n=st.integers(min_value=3, max_value=30))
def test_random_move_sequences_keep_the_running_fitness_exact(seed: int, n: int) -> None:
    rng = np.random.default_rng(seed)
    problem = random_matrix_problem(rng, n)
    w = problem.weights
    t = tour_of(rng.permutation(n).astype(np.int64))
    running = tour_fitness(t, w)
    for _ in range(40):
        if rng.random() < 0.5 and n >= 2:
            i, j = sorted(rng.choice(np.arange(1, n + 1), size=2, replace=False).tolist())
            running += swap_delta(t, w, i, j)
            apply_swap(t, i, j)
        else:
            i = int(rng.integers(1, n + 1))
            j = int(rng.choice([k for k in range(n + 1) if k not in (i, i - 1)]))
            running += insert_delta(t, w, i, j)
            apply_insert(t, i, j)
        assert running == tour_fitness(t, w)


# ---- exact validators ---------------------------------------------------------------------------


def test_toy_optimum() -> None:
    problem = toy_problem()
    for solver in (brute_force, held_karp):
        result = solver(problem)
        assert result.fitness == 6 and result.order.tolist() == [0, 1, 2]
        assert result.merged_length == 9


@settings(deadline=None, max_examples=40)
@given(seed=st.integers(min_value=0, max_value=10**6), n=st.integers(min_value=2, max_value=8))
def test_brute_force_and_held_karp_agree(seed: int, n: int) -> None:
    problem = random_matrix_problem(np.random.default_rng(seed), n)
    a, b = brute_force(problem), held_karp(problem)
    assert a.fitness == b.fitness
    assert valid(a.order, n) and valid(b.order, n)
    assert problem.fitness(b.order) == b.fitness


def test_exact_solvers_refuse_instances_that_are_too_large() -> None:
    rng = np.random.default_rng(0)
    with pytest.raises(ValueError, match="limited"):
        brute_force(random_matrix_problem(rng, 10))
    with pytest.raises(ValueError, match="limited"):
        held_karp(random_matrix_problem(rng, 16))
    assert held_karp(random_matrix_problem(rng, 1)).order.tolist() == [0]


# ---- greedy and nearest neighbour ---------------------------------------------------------------


def test_greedy_and_nearest_neighbour_on_the_toy_example() -> None:
    problem = toy_problem()
    assert greedy(problem).order.tolist() == [0, 1, 2] and greedy(problem).fitness == 6
    nn = nearest_neighbour(problem)
    assert valid(nn.order, 3) and nn.fitness == 6


@settings(deadline=None, max_examples=40)
@given(seed=st.integers(min_value=0, max_value=10**6), n=st.integers(min_value=1, max_value=40))
def test_baselines_return_valid_deterministic_orders(seed: int, n: int) -> None:
    problem = random_matrix_problem(np.random.default_rng(seed), n)
    for solver in (greedy, nearest_neighbour):
        a, b = solver(problem), solver(problem)
        assert valid(a.order, n) and np.array_equal(a.order, b.order)
        assert a.fitness == problem.fitness(a.order)


def test_greedy_never_closes_a_cycle_and_takes_the_largest_overlap_first() -> None:
    ov = np.zeros((3, 3), dtype=np.int32)
    ov[0, 1], ov[1, 2], ov[2, 0] = 9, 8, 7  # a 3-cycle: the last edge must be refused
    result = greedy(Problem.from_matrix(ov, np.full(3, 20)))
    assert result.order.tolist() == [0, 1, 2] and result.fitness == 17


def test_greedy_matches_the_true_order_on_error_free_real_reads() -> None:
    """Observed on phiX174 (regression test, not a theorem): without errors greedy is optimal."""
    genome = load_genome("NC_001422.1", ROOT / "data" / "raw").genome
    reads = simulate_reads(
        genome, SimulationConfig(coverage=5, read_length=ReadLengthModel.fixed(100), seed=1)
    )
    inst = AssemblyInstance.from_readset(reads)
    result = greedy(Problem.from_instance(inst))
    assert result.fitness == inst.fitness(inst.true_order())
    assert result.fitness == inst.fitness(result.order)


# ---- tabu search --------------------------------------------------------------------------------


@settings(deadline=None, max_examples=15)
@given(seed=st.integers(min_value=0, max_value=10**6), n=st.integers(min_value=3, max_value=12))
def test_tabu_returns_valid_solutions_no_worse_than_greedy(seed: int, n: int) -> None:
    problem = random_matrix_problem(np.random.default_rng(seed), n)
    result = tabu_search(problem, Budget(evaluations=2000), seed=seed)
    assert valid(result.order, n)
    assert result.fitness == problem.fitness(result.order) >= greedy(problem).fitness


def test_tabu_reaches_the_optimum_on_small_instances() -> None:
    rng = np.random.default_rng(5)
    for k in range(10):
        problem = random_matrix_problem(rng, int(rng.integers(6, 13)))
        result = tabu_search(problem, Budget(evaluations=30000), seed=k)
        assert result.fitness == held_karp(problem).fitness


def test_tabu_is_reproducible_per_seed_with_an_evaluation_budget() -> None:
    problem = random_matrix_problem(np.random.default_rng(9), 30)
    a = tabu_search(problem, Budget(evaluations=40000), seed=3)
    b = tabu_search(problem, Budget(evaluations=40000), seed=3)
    c = tabu_search(problem, Budget(evaluations=40000), seed=4)
    assert np.array_equal(a.order, b.order) and a.fitness == b.fitness
    assert a.evaluations == b.evaluations >= 40000
    assert valid(c.order, 30)


def test_tabu_options_and_restarts() -> None:
    problem = random_matrix_problem(np.random.default_rng(2), 40)
    sampled = tabu_search(
        problem, Budget(evaluations=60000), seed=1, moves_per_iteration=200, init="random"
    )
    assert valid(sampled.order, 40)
    restarted = tabu_search(problem, Budget(evaluations=60000), seed=1, restart_after=5, kick=3)
    assert restarted.info["restarts"] > 0 and valid(restarted.order, 40)
    with pytest.raises(ValueError, match="init"):
        tabu_search(problem, Budget(evaluations=100), init="best")
    assert valid(
        tabu_search(
            random_matrix_problem(np.random.default_rng(0), 2), Budget(evaluations=10)
        ).order,
        2,
    )


def test_tabu_candidate_lists_hold_the_best_partners_by_overlap() -> None:
    problem = random_matrix_problem(np.random.default_rng(11), 25)
    w = np.ascontiguousarray(problem.weights)
    pred, succ = _candidate_lists(w, 4)
    for x in range(1, 26):
        for lists, column in ((pred, w[:, x]), (succ, w[x, :])):
            col = column.astype(np.int64)
            col[0] = col[x] = 0
            chosen = [int(a) for a in lists[x] if a >= 0]
            assert x not in chosen and 0 not in chosen and len(set(chosen)) == len(chosen)
            # the values of the chosen partners are the largest ones, in decreasing order
            expected = sorted(col[col > 0], reverse=True)[:4]
            assert [int(col[a]) for a in chosen] == [int(v) for v in expected]


def test_tabu_without_candidate_moves_keeps_its_random_stream() -> None:
    problem = random_matrix_problem(np.random.default_rng(4), 40)
    plain = tabu_search(problem, Budget(evaluations=30000), seed=2, moves_per_iteration=200)
    zero = tabu_search(
        problem, Budget(evaluations=30000), seed=2, moves_per_iteration=200, candidate_fraction=0.0
    )
    assert np.array_equal(plain.order, zero.order) and plain.evaluations == zero.evaluations


@pytest.mark.parametrize("fraction", [0.25, 0.9, 1.0])
def test_tabu_candidate_moves_are_valid_consistent_and_reproducible(fraction: float) -> None:
    # the position table that the candidate moves rely on is checked through the final fitness:
    # a stale position would make the incremental fitness differ from the recomputed one
    rng = np.random.default_rng(8)
    for n in (20, 60):
        problem = random_matrix_problem(rng, n)
        a, b = (
            tabu_search(
                problem,
                Budget(evaluations=60000),
                seed=1,
                moves_per_iteration=5 * n,
                restart_after=30,
                kick=3,
                candidate_fraction=fraction,
            )
            for _ in range(2)
        )
        assert valid(a.order, n)
        assert a.fitness == problem.fitness(a.order) >= greedy(problem).fitness
        assert np.array_equal(a.order, b.order) and a.evaluations == b.evaluations
        assert a.info["restarts"] > 0


def test_tabu_candidate_moves_work_on_real_reads() -> None:
    full = load_genome("NC_001422.1", ROOT / "data" / "raw").genome
    window = Genome.from_sequence("window", full.sequence[:1500])
    reads = simulate_reads(
        window, SimulationConfig(coverage=8, read_length=ReadLengthModel.fixed(100), seed=2)
    )
    inst = AssemblyInstance.from_readset(reads)
    problem = Problem.from_instance(inst)
    result = tabu_search(
        problem,
        Budget(evaluations=150000),
        seed=1,
        moves_per_iteration=3 * problem.n,
        candidate_fraction=0.9,
    )
    assert valid(result.order, problem.n)
    assert result.fitness == inst.fitness(result.order) >= greedy(problem).fitness


# ---- ALNS ---------------------------------------------------------------------------------------


@pytest.mark.parametrize("policy", POLICIES)
def test_alns_reaches_the_optimum_on_small_instances(policy: str) -> None:
    rng = np.random.default_rng(6)
    for k in range(6):
        problem = random_matrix_problem(rng, int(rng.integers(6, 11)))
        result = alns(problem, Budget(evaluations=150000), seed=k, policy=policy)
        assert valid(result.order, problem.n)
        assert result.fitness == held_karp(problem).fitness


@pytest.mark.parametrize("arm", range(len(ARM_NAMES)))
def test_every_single_operator_pair_alone_gives_valid_solutions(arm: int) -> None:
    problem = random_matrix_problem(np.random.default_rng(arm), 25)
    result = alns(problem, Budget(evaluations=20000), seed=1, policy="fixed", fixed_arm=arm)
    assert valid(result.order, 25) and result.fitness == problem.fitness(result.order)
    used = np.array(result.info["outcomes"]).sum(axis=1)
    assert used[arm] == used.sum() > 0  # only the chosen operator pair was used


@pytest.mark.parametrize("policy", POLICIES)
def test_alns_is_reproducible_per_seed_with_an_evaluation_budget(policy: str) -> None:
    problem = random_matrix_problem(np.random.default_rng(4), 30)
    a = alns(problem, Budget(evaluations=60000), seed=8, policy=policy)
    b = alns(problem, Budget(evaluations=60000), seed=8, policy=policy)
    assert np.array_equal(a.order, b.order) and a.fitness == b.fitness
    assert a.evaluations == b.evaluations and np.array_equal(
        a.info["arm_history"], b.info["arm_history"]
    )


def test_alns_bookkeeping_and_temperature_calibration() -> None:
    problem = random_matrix_problem(np.random.default_rng(7), 35)
    result = alns(problem, Budget(evaluations=80000), seed=2, policy="ucb1", trace_every=10)
    outcomes = np.array(result.info["outcomes"])
    assert outcomes.sum() == result.info["iterations"] == len(result.info["arm_history"])
    assert result.info["t0"] > 0 and result.evaluations >= 80000
    assert len(result.info["probability_trace"]) == result.info["iterations"] // 10
    assert all(abs(sum(p) - 1) < 1e-9 for p in result.info["probability_trace"])
    assert [f for _, _, f in result.trace] == sorted(f for _, _, f in result.trace)


def test_alns_validates_its_inputs() -> None:
    problem = random_matrix_problem(np.random.default_rng(0), 12)
    with pytest.raises(ValueError, match="unknown policy"):
        alns(problem, Budget(evaluations=100), policy="magic")
    assert valid(
        alns(random_matrix_problem(np.random.default_rng(0), 3), Budget(evaluations=10)).order, 3
    )


# ---- LKH-3 (needs the build in tools/; skipped when it is absent) -------------------------------


def test_wsl_path_conversion() -> None:
    assert to_wsl_path(Path("E:/a b/c.txt")).startswith("/mnt/e/")
    assert to_wsl_path(Path("E:/a b/c.txt")).endswith("/a b/c.txt")


def test_tour_file_parsing_and_verification_reject_bad_tours() -> None:
    problem = toy_problem()  # nodes: dummy plus 3 reads; optimal order A, B, C costs 9
    text = "COMMENT : Length = 9\nTYPE : TOUR\nDIMENSION : 4\nTOUR_SECTION\n1\n2\n3\n4\n-1\nEOF\n"
    ids, cost = parse_tour_file(text)
    assert ids.tolist() == [1, 2, 3, 4] and cost == 9
    assert verify_tour(problem, ids, cost).tolist() == [0, 1, 2]
    with pytest.raises(LkhError, match="Hamiltonian"):
        verify_tour(problem, np.array([1, 2, 2, 4]), 9)
    with pytest.raises(LkhError, match="Hamiltonian"):
        verify_tour(problem, np.array([1, 2, 3]), 9)
    with pytest.raises(LkhError, match="reports cost"):
        verify_tour(problem, ids, 8)
    with pytest.raises(LkhError, match="TOUR_SECTION"):
        parse_tour_file("no tour here")


def test_missing_lkh_is_reported_and_never_replaced(tmp_path: Path) -> None:
    with pytest.raises(LkhError, match="never replaced"):
        lkh(toy_problem(), seconds=1, binary=tmp_path / "does-not-exist")


@pytest.mark.skipif(lkh_binary() is None, reason="LKH-3 is not built in tools/")
def test_lkh_finds_the_optimum_of_small_instances_and_its_tour_verifies() -> None:
    rng = np.random.default_rng(12)
    for k in range(3):
        problem = random_matrix_problem(rng, int(rng.integers(8, 13)))
        result = lkh(problem, seconds=1, seed=k)
        assert valid(result.order, problem.n)
        assert result.fitness == held_karp(problem).fitness
        tour = objective.tour_from_order(result.order)
        assert objective.cycle_cost(problem.cost_matrix(), tour) == result.merged_length
    assert order_of(tour_of(np.arange(3))).tolist() == [0, 1, 2]
