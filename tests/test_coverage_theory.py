"""Coverage measurement and theory: hand examples, brute-force enumeration and invariants."""

from __future__ import annotations

import itertools

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from genomeasm import coverage as cov


def arr(values: list[int]) -> np.ndarray:
    return np.array(values, dtype=np.int64)


# ---- measurement ------------------------------------------------------------------------


def test_depth_hand_example() -> None:
    depth = cov.coverage_depth(10, arr([0, 3, 3]), arr([4, 3, 2]))
    assert depth.tolist() == [1, 1, 1, 3, 2, 1, 0, 0, 0, 0]
    assert cov.uncovered_fraction(depth) == 0.4
    assert cov.depth_histogram(depth, 5).tolist() == pytest.approx([0.4, 0.4, 0.1, 0.1, 0.0])


@pytest.mark.parametrize(
    ("starts", "spans"),
    [([-1], [3]), ([8], [3]), ([0], [0]), ([0, 1], [3])],
)
def test_depth_rejects_reads_outside_the_genome(starts: list[int], spans: list[int]) -> None:
    with pytest.raises(ValueError):
        cov.coverage_depth(10, arr(starts), arr(spans))


@settings(deadline=None, max_examples=80)
@given(
    genome_length=st.integers(min_value=20, max_value=300),
    data=st.data(),
)
def test_depth_sums_to_total_span(genome_length: int, data: st.DataObject) -> None:
    n = data.draw(st.integers(min_value=1, max_value=40))
    spans = [data.draw(st.integers(min_value=1, max_value=genome_length)) for _ in range(n)]
    starts = [data.draw(st.integers(min_value=0, max_value=genome_length - s)) for s in spans]
    depth = cov.coverage_depth(genome_length, arr(starts), arr(spans))
    assert depth.sum() == sum(spans)
    assert np.min(depth) >= 0


def test_contigs_hand_examples() -> None:
    starts, spans = arr([0, 10, 20]), arr([15, 15, 15])  # neighbours overlap by 5 bases
    assert cov.count_contigs(starts, spans, 5) == 1
    assert cov.count_contigs(starts, spans, 6) == 3
    # a read contained in an earlier one overlaps by its own length
    starts, spans = arr([0, 20]), arr([100, 10])
    assert cov.count_contigs(starts, spans, 10) == 1
    assert cov.count_contigs(starts, spans, 11) == 2
    assert cov.count_contigs(arr([]), arr([]), 5) == 0


def test_contig_count_does_not_depend_on_input_order() -> None:
    rng = np.random.default_rng(3)
    starts = rng.integers(0, 900, size=60)
    spans = np.full(60, 100, dtype=np.int64)
    expected = cov.count_contigs(starts, spans, 30)
    for _ in range(5):
        order = rng.permutation(60)
        assert cov.count_contigs(starts[order], spans[order], 30) == expected


def test_contig_threshold_must_be_positive() -> None:
    with pytest.raises(ValueError, match="at least 1"):
        cov.count_contigs(arr([0]), arr([10]), 0)


# ---- theory against brute-force enumeration -----------------------------------------------


def test_exact_uncovered_fraction_equals_full_enumeration_fixed_length() -> None:
    g, ell, n = 6, 3, 2
    starts_all = range(g - ell + 1)
    total = 0.0
    outcomes = list(itertools.product(starts_all, repeat=n))
    for combo in outcomes:
        depth = cov.coverage_depth(g, arr(list(combo)), arr([ell] * n))
        total += cov.uncovered_fraction(depth)
    brute = total / len(outcomes)
    formula = cov.expected_uncovered_fraction(g, n, arr([ell]), np.array([1.0]))
    assert formula == pytest.approx(brute, abs=1e-15)


def test_exact_uncovered_fraction_equals_full_enumeration_mixed_lengths() -> None:
    g, n = 6, 2
    support, pmf = arr([2, 3]), np.array([0.5, 0.5])
    reads = [(ell, s, 0.5 / (g - ell + 1)) for ell in (2, 3) for s in range(g - ell + 1)]
    brute = 0.0
    for (l1, s1, p1), (l2, s2, p2) in itertools.product(reads, repeat=n):
        depth = cov.coverage_depth(g, arr([s1, s2]), arr([l1, l2]))
        brute += p1 * p2 * cov.uncovered_fraction(depth)
    assert cov.expected_uncovered_fraction(g, n, support, pmf) == pytest.approx(brute, abs=1e-15)


def test_exact_depth_histogram_equals_full_enumeration() -> None:
    g, ell, n, bins = 7, 3, 3, 6
    outcomes = list(itertools.product(range(g - ell + 1), repeat=n))
    brute = np.zeros(bins)
    for combo in outcomes:
        depth = cov.coverage_depth(g, arr(list(combo)), arr([ell] * n))
        brute += cov.depth_histogram(depth, bins)
    brute /= len(outcomes)
    formula = cov.expected_depth_histogram(g, n, arr([ell]), np.array([1.0]), bins)
    assert formula == pytest.approx(brute, abs=1e-14)


def test_exact_contig_formula_equals_full_enumeration() -> None:
    ell, theta, n, n_starts = 5, 4, 3, 4  # max gap x = l - theta = 1; genome length 8
    outcomes = list(itertools.product(range(n_starts), repeat=n))
    total = sum(cov.count_contigs(arr(list(c)), arr([ell] * n), theta) for c in outcomes)
    assert cov.expected_contigs_exact(n_starts, n, ell - theta) == pytest.approx(
        total / len(outcomes), abs=1e-12
    )


def test_exact_contig_formula_for_one_read_and_zero_gap() -> None:
    assert cov.expected_contigs_exact(100, 1, 20) == pytest.approx(1.0)
    with pytest.raises(ValueError):
        cov.expected_contigs_exact(10, 3, -1)


# ---- invariants of the exact model ----------------------------------------------------------


@pytest.mark.parametrize(("g", "ell", "n"), [(5386, 100, 269), (2000, 50, 400), (1000, 100, 30)])
def test_expected_histogram_has_total_mass_one_and_mean_depth_c(g: int, ell: int, n: int) -> None:
    bins = n + 1  # depth can never exceed the number of reads
    hist = cov.expected_depth_histogram(g, n, arr([ell]), np.array([1.0]), bins)
    assert hist.sum() == pytest.approx(1.0, abs=1e-12)
    assert np.dot(np.arange(bins), hist) == pytest.approx(n * ell / g, rel=1e-12)


def test_cover_probability_is_one_over_s_at_the_ends_and_flat_inside() -> None:
    g, ell = 1000, 100
    q = cov.start_cover_probability(g, arr([ell]), np.array([1.0]))
    s = g - ell + 1
    assert q[0] == pytest.approx(1 / s)
    assert q[-1] == pytest.approx(1 / s)
    assert np.allclose(q[ell - 1 : g - ell + 1], ell / s)
    assert q.sum() == pytest.approx(ell)  # a read covers l positions on average


def test_lander_waterman_functions() -> None:
    assert cov.lw_uncovered_fraction(0.0) == 1.0
    assert cov.lw_uncovered_fraction(3.0) == pytest.approx(np.exp(-3))
    assert cov.lw_depth_pmf(4.0, 80).sum() == pytest.approx(1.0)
    assert cov.lw_depth_pmf(4.0, 5)[0] == pytest.approx(np.exp(-4))
    assert cov.lw_expected_contigs(100, 5.0, 0, 100) == pytest.approx(100 * np.exp(-5))
    assert cov.nominal_coverage(269, 100, 5386) == pytest.approx(269 * 100 / 5386)


def test_exact_contigs_approach_lander_waterman_for_a_long_genome() -> None:
    ell, theta, g = 100, 30, 2_000_000
    n = round(5 * g / ell)
    exact = cov.expected_contigs_exact(g - ell + 1, n, ell - theta)
    lw = cov.lw_expected_contigs(n, n * ell / g, theta, ell)
    # Reads that start at the same position are one occupied position, so the exact count is
    # lower by a factor (1 - exp(-lambda)) / lambda with lambda = N/S = c/l, about 2.5% here.
    lam = n / (g - ell + 1)
    assert exact == pytest.approx(lw * (1 - np.exp(-lam)) / lam, rel=1e-3)
    assert exact < lw


# ---- bin grouping -----------------------------------------------------------------------------


def test_bin_groups_are_contiguous_complete_and_large_enough() -> None:
    g = 5386
    exact = cov.expected_depth_histogram(g, 269, arr([100]), np.array([1.0]), 50)
    groups = cov.group_bins(exact, g, 100.0)
    assert groups[0][0] == 0 and groups[-1][1] == 50
    assert all(a[1] == b[0] for a, b in itertools.pairwise(groups))
    for a, b in groups:
        assert exact[a:b].sum() * g >= 100.0
    assert cov.grouped_fractions(exact, groups).sum() == pytest.approx(1.0)


def test_grouped_fractions_put_overflow_into_the_last_group() -> None:
    hist = np.array([0.5, 0.3, 0.1])  # 0.1 of the positions are deeper than the last bin
    grouped = cov.grouped_fractions(hist, [(0, 1), (1, 3)])
    assert grouped.tolist() == pytest.approx([0.5, 0.5])
    runs = cov.grouped_fractions(np.vstack([hist, hist]), [(0, 1), (1, 3)])
    assert runs.shape == (2, 2)


def test_single_group_when_nothing_reaches_the_minimum() -> None:
    assert cov.group_bins(np.array([0.001, 0.001]), 100, 100.0) == [(0, 2)]
