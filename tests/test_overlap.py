"""Overlaps: the toy example, equivalence with the naive definition, tolerance, preprocessing."""

from __future__ import annotations

import itertools
from pathlib import Path

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from genomeasm import objective
from genomeasm.instance import AssemblyInstance
from genomeasm.io import load_genome
from genomeasm.overlap import (
    expected_spurious_pairs,
    is_substring_free,
    make_substring_free,
    mismatch_limit,
    overlap_matrix,
    overlap_matrix_naive,
    overlap_matrix_tolerant,
    recommend_min_overlap,
    truth_agreement,
)
from genomeasm.simulate import ReadLengthModel, SimulationConfig, simulate_reads

ROOT = Path(__file__).resolve().parent.parent

TOY = ["ATGCG", "GCGTA", "GTACC"]  # reads A, B, C of the master prompt
TOY_OV = [[0, 3, 1], [1, 0, 3], [0, 0, 0]]  # [[-,3,1],[1,-,3],[0,0,-]] with "-" stored as 0


# ---- (a) the toy example from the master prompt ---------------------------------------------


def test_toy_overlap_matrix_matches_the_master_prompt() -> None:
    assert overlap_matrix(TOY, 1).tolist() == TOY_OV
    assert overlap_matrix_naive(TOY, 1).tolist() == TOY_OV
    assert overlap_matrix_tolerant(TOY, 1, 0.0).tolist() == TOY_OV


def test_toy_overlaps_are_directional() -> None:
    ov = overlap_matrix(TOY, 1)
    assert ov[0, 1] == 3 and ov[1, 0] == 1  # A then B saves 3, B then A saves only 1
    assert ov[2, 0] == 0 and ov[0, 2] == 1


def test_toy_best_order_is_a_b_c_with_merged_string() -> None:
    inst = AssemblyInstance.from_sequences(TOY, min_overlap=1)
    scores = {p: inst.fitness(p) for p in itertools.permutations(range(3))}
    best = max(scores, key=lambda p: scores[p])
    assert best == (0, 1, 2) and scores[best] == 6
    assert inst.merged_length(best) == 9 == sum(map(len, TOY)) - 6
    assert inst.assemble(best).sequences == ["ATGCGTACC"]
    # the order with the smallest ATSP tour cost is the same one
    cost = inst.cost()
    tours = {p: objective.cycle_cost(cost, objective.tour_from_order(p)) for p in scores}
    assert min(tours, key=lambda p: tours[p]) == best and tours[best] == 9


# ---- (d) the fast code equals the naive definition ----------------------------------------------


@settings(deadline=None, max_examples=150)
@given(
    reads=st.lists(st.text(alphabet="ACGT", min_size=2, max_size=22), min_size=1, max_size=14),
    tau=st.integers(min_value=1, max_value=5),
)
def test_fast_overlap_equals_naive_on_random_reads(reads: list[str], tau: int) -> None:
    naive = overlap_matrix_naive(reads, tau)
    assert np.array_equal(overlap_matrix(reads, tau), naive)
    assert np.array_equal(overlap_matrix_tolerant(reads, tau, 0.0), naive)


@settings(deadline=None, max_examples=60)
@given(
    genome=st.text(alphabet="ACGT", min_size=40, max_size=120),
    data=st.data(),
    tau=st.integers(min_value=2, max_value=8),
)
def test_fast_overlap_equals_naive_on_reads_cut_from_one_genome(
    genome: str, data: st.DataObject, tau: int
) -> None:
    reads = []
    for _ in range(data.draw(st.integers(min_value=2, max_value=12))):
        length = data.draw(st.integers(min_value=6, max_value=30))
        start = data.draw(st.integers(min_value=0, max_value=len(genome) - length))
        reads.append(genome[start : start + length])
    naive = overlap_matrix_naive(reads, tau)
    assert np.array_equal(overlap_matrix(reads, tau), naive)
    assert np.array_equal(overlap_matrix_tolerant(reads, tau, 0.0), naive)


def test_overlap_respects_the_length_bounds() -> None:
    reads = ["ACGTACGT", "ACGT", "TTACGT"]
    ov = overlap_matrix(reads, 2)
    for i, j in itertools.permutations(range(3), 2):
        assert ov[i, j] == 0 or 2 <= ov[i, j] < min(len(reads[i]), len(reads[j]))
    # an overlap equal to the shorter read would be containment, which is excluded by k < min
    assert overlap_matrix_naive(["AAAA", "AAAAAA"], 1)[0, 1] == 3


def test_min_overlap_must_be_positive() -> None:
    for func in (overlap_matrix, overlap_matrix_naive):
        with pytest.raises(ValueError, match="at least 1"):
            func(TOY, 0)
    with pytest.raises(ValueError):
        overlap_matrix_tolerant(TOY, 1, 1.0)


# ---- error-tolerant overlaps --------------------------------------------------------------------


@pytest.mark.parametrize(
    ("k", "tolerance", "limit"),
    [(8, 0.0, 0), (20, 0.05, 1), (19, 0.05, 0), (100, 0.29, 29), (100, 0.07, 7), (10, 0.15, 1)],
)
def test_mismatch_limit_is_floor_of_tolerance_times_k(k: int, tolerance: float, limit: int) -> None:
    assert mismatch_limit(k, tolerance) == limit  # 0.29 * 100 is 28.999999999999996 in floats


def test_tolerant_overlap_accepts_a_mismatch_that_exact_overlap_rejects() -> None:
    a, b = (
        "TTTTACGTACGT",
        "ACGTACGAAAAA",
    )  # the last 8 letters of a differ from b's first 8 at 1 place
    assert overlap_matrix_naive([a, b], 5)[0, 1] == 0
    assert overlap_matrix_tolerant([a, b], 5, 0.15)[0, 1] == 8  # floor(0.15 * 8) = 1 mismatch
    assert overlap_matrix_tolerant([a, b], 5, 0.10)[0, 1] == 0  # floor(0.10 * 8) = 0


def test_tolerant_overlap_is_never_shorter_than_the_exact_one() -> None:
    rng = np.random.default_rng(5)
    genome = "".join(rng.choice(list("ACGT"), size=300))
    reads = [genome[s : s + 40] for s in rng.integers(0, 260, size=25)]
    exact = overlap_matrix(reads, 8)
    tolerant = overlap_matrix_tolerant(reads, 8, 0.1)
    assert np.all(tolerant >= exact)


# ---- substring-free preprocessing ---------------------------------------------------------------


def test_duplicates_and_contained_reads_are_removed_and_counted() -> None:
    reads = ["ACGTAC", "CGTA", "ACGTAC", "TTTTTT", "GTAC", "AAAACGTACGG"]
    pre = make_substring_free(reads)
    assert pre.n_input == 6 and pre.n_duplicates == 1  # the second ACGTAC
    # CGTA and GTAC are inside ACGTAC; ACGTAC is inside AAAACGTACGG; TTTTTT stays
    assert pre.reads == ["TTTTTT", "AAAACGTACGG"]
    assert pre.kept.tolist() == [3, 5]
    assert pre.n_contained == 3
    assert is_substring_free(pre.reads)


def test_preprocessing_of_empty_and_single_inputs() -> None:
    assert make_substring_free([]).reads == []
    single = make_substring_free(["ACGT"])
    assert single.reads == ["ACGT"] and single.kept.tolist() == [0]


def test_equal_length_reads_are_only_removed_as_duplicates() -> None:
    pre = make_substring_free(["ACGT", "CGTA", "ACGT", "GTAC"])
    assert pre.reads == ["ACGT", "CGTA", "GTAC"] and pre.n_contained == 0 and pre.n_duplicates == 1


@settings(deadline=None, max_examples=150)
@given(reads=st.lists(st.text(alphabet="AC", min_size=1, max_size=9), min_size=0, max_size=14))
def test_preprocessing_property(reads: list[str]) -> None:
    pre = make_substring_free(reads)
    assert is_substring_free(pre.reads)
    assert [reads[i] for i in pre.kept] == pre.reads
    assert pre.n_input == len(reads)
    assert pre.n_input - pre.n_duplicates - pre.n_contained == len(pre.reads)
    for read in reads:  # nothing is lost: every read is kept or sits inside a kept read
        assert any(read in kept for kept in pre.reads)
    assert list(pre.kept) == sorted(pre.kept)


# ---- choosing tau -------------------------------------------------------------------------------


@pytest.mark.parametrize(("n", "tau"), [(1, 8), (2, 8), (539, 10), (1000, 11), (10000, 14)])
def test_recommended_min_overlap_values(n: int, tau: int) -> None:
    assert recommend_min_overlap(n) == tau


@given(n=st.integers(min_value=2, max_value=200_000))
def test_recommended_overlap_meets_the_budget_and_is_minimal(n: int) -> None:
    tau = recommend_min_overlap(n, minimum=1)
    assert expected_spurious_pairs(n, tau) <= 1.0 + 1e-9
    if tau > 1:
        assert expected_spurious_pairs(n, tau - 1) > 1.0


def test_recommended_overlap_grows_with_the_number_of_reads() -> None:
    values = [recommend_min_overlap(n, minimum=1) for n in (10, 100, 1000, 10_000, 100_000)]
    assert values == sorted(values) and values[0] < values[-1]


# ---- comparison with the ground truth of simulated reads ----------------------------------------


def test_truth_agreement_hand_example() -> None:
    starts = np.array([0, 50, 100], dtype=np.int64)
    spans = np.array([100, 100, 100], dtype=np.int64)
    ov = np.zeros((3, 3), dtype=np.int32)
    ov[0, 1] = 50  # the true overlap
    ov[1, 2] = 60  # a longer chance match hides the true overlap of 50
    ov[0, 2] = 12  # reads 0 and 2 do not overlap at all
    counts = truth_agreement(ov, starts, spans, 10)
    assert (counts["true_pairs"], counts["recovered"], counts["exact"], counts["longer"]) == (
        2,
        2,
        1,
        1,
    )
    assert (counts["edges"], counts["spurious"]) == (3, 1)
    assert counts["hist_true"][50] == 1 and counts["hist_true"][60] == 1
    assert counts["hist_other"][12] == 1


def test_overlap_graph_of_real_error_free_reads_finds_every_true_overlap() -> None:
    genome = load_genome("NC_001422.1", ROOT / "data" / "raw").genome
    reads = simulate_reads(
        genome, SimulationConfig(coverage=8, read_length=ReadLengthModel.fixed(100), seed=31)
    )
    inst = AssemblyInstance.from_readset(reads)
    assert inst.true_starts is not None and inst.true_spans is not None
    truth = truth_agreement(inst.ov, inst.true_starts, inst.true_spans, inst.min_overlap)
    assert truth["true_pairs"] > 1000
    assert truth["recovered"] == truth["true_pairs"]  # recall 1 for error-free reads
    assert truth["exact"] >= 0.999 * truth["true_pairs"]
    assert truth["spurious"] <= 5
