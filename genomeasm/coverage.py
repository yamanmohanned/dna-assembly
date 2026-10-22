"""Coverage measurement and coverage theory (Lander-Waterman and the exact finite-genome model).

Plain-language story: imagine throwing N short strips of paper at random onto a long strip
(the genome). Coverage is how many layers of paper lie on top of each point. Some points get
none by bad luck. This module counts the layers and predicts, with mathematics, how many points
should be bare, so that the simulator can be checked against theory.

Definitions (genome length G, N reads, read length l, one strand, reads lie fully inside the
genome, start positions uniform; decision D-017):

* nominal coverage        c = N * l / G
* depth at position x     number of reads that cover x
* Lander-Waterman (LW)    approximation: depth ~ Poisson(c), uncovered fraction ~ exp(-c).
  It ignores the genome ends and treats the genome as very long compared with a read.
* exact finite-genome     with S = G - l + 1 possible start positions, a read covers position x
  with probability q(x) = w(x) / S, where w(x) is the number of starts whose read contains x:
  w(x) = min(x, S-1) - max(0, x-l+1) + 1. Reads are independent, so the depth at x is
  Binomial(N, q(x)) and P(x uncovered) = (1 - q(x))**N. This is exact for this simulator.
  Near the two ends q(x) is smaller than in the middle, which makes the uncovered fraction
  larger than exp(-c): that edge effect is real and is reported, not hidden.

Contigs (derived here, verified by Monte Carlo in tests; decision D-021): reads are in the same
"true contig" if consecutive reads (sorted by start) overlap by at least theta bases. For fixed
l this means consecutive occupied start positions differ by at most x = l - theta. A contig
starts at an occupied start position u when no read starts in the x positions before u. With
m_u = min(x, u) and independent uniform starts,

    P(u occupied and window empty) = (1 - m_u/S)**N - (1 - (m_u+1)/S)**N

so E[#contigs] = sum over u of that probability (exact). For long genomes this is close to the
Lander-Waterman value N * exp(-c * (1 - theta/l)).
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from scipy.stats import binom, poisson

IntArray = NDArray[np.int64]
FloatArray = NDArray[np.float64]


# --------------------------------------------------------------------------------------
# Measurement on a set of read positions
# --------------------------------------------------------------------------------------


def coverage_depth(genome_length: int, starts: IntArray, spans: IntArray) -> IntArray:
    """Depth at every genome position: how many reads cover it (difference array + cumsum)."""
    if starts.shape != spans.shape:
        raise ValueError("starts and spans must have the same shape")
    ends = starts + spans
    if starts.size and (np.min(starts) < 0 or np.max(ends) > genome_length or np.min(spans) < 1):
        raise ValueError("every read must lie inside the genome and have a positive span")
    delta = np.bincount(starts, minlength=genome_length + 1) - np.bincount(
        ends, minlength=genome_length + 1
    )
    return np.cumsum(delta[:genome_length]).astype(np.int64)


def uncovered_fraction(depth: IntArray) -> float:
    """Share of genome positions that no read covers."""
    return float(np.mean(depth == 0))


def depth_histogram(depth: IntArray, n_bins: int) -> FloatArray:
    """Fraction of positions having depth 0, 1, ..., n_bins-1 (deeper positions are not listed)."""
    counts = np.bincount(depth, minlength=n_bins)[:n_bins]
    return counts.astype(np.float64) / depth.size


def group_bins(
    expected: FloatArray, genome_length: int, min_positions: float
) -> list[tuple[int, int]]:
    """Merge neighbouring depth bins until each group expects >= ``min_positions`` positions.

    Bins are scanned from depth 0 upwards and a group is closed as soon as its expected number of
    positions (probability x genome length) reaches the minimum; a leftover at the top joins the
    last group. This is the usual way to avoid bins that almost never occur, which have no
    measurable spread and would make a goodness-of-fit comparison meaningless.
    """
    groups: list[tuple[int, int]] = []
    start, accumulated = 0, 0.0
    for i, mass in enumerate(expected):
        accumulated += float(mass) * genome_length
        if accumulated >= min_positions:
            groups.append((start, i + 1))
            start, accumulated = i + 1, 0.0
    if start < len(expected):
        if groups:
            groups[-1] = (groups[-1][0], len(expected))
        else:
            groups.append((0, len(expected)))
    return groups


def grouped_fractions(histogram: FloatArray, groups: list[tuple[int, int]]) -> FloatArray:
    """Sum depth-histogram bins into groups; depths beyond the last bin join the last group.

    ``histogram`` has shape ``(bins,)`` or ``(runs, bins)``; the result has one column per group.
    """
    hist = np.atleast_2d(histogram)
    overflow = 1.0 - hist.sum(axis=1)
    out = np.stack([hist[:, a:b].sum(axis=1) for a, b in groups], axis=1)
    out[:, -1] += overflow
    return out if histogram.ndim == 2 else out[0]


def count_contigs(starts: IntArray, spans: IntArray, min_overlap: int) -> int:
    """Number of true contigs: groups of reads chained by reference overlaps >= min_overlap."""
    if min_overlap < 1:
        raise ValueError("min_overlap must be at least 1")
    if starts.size == 0:
        return 0
    order = np.argsort(starts, kind="stable")
    s = starts[order]
    e = (starts + spans)[order]
    reach = np.maximum.accumulate(e)[:-1]  # furthest end among the earlier reads
    overlap_with_best_earlier = np.minimum(reach, e[1:]) - s[1:]
    return 1 + int(np.sum(overlap_with_best_earlier < min_overlap))


# --------------------------------------------------------------------------------------
# Theory
# --------------------------------------------------------------------------------------


def nominal_coverage(n_reads: int, mean_read_length: float, genome_length: int) -> float:
    """c = N * l / G."""
    return n_reads * mean_read_length / genome_length


def start_cover_probability(
    genome_length: int, length_support: IntArray, length_pmf: FloatArray
) -> FloatArray:
    """q(x): probability that one random read covers position x, for every x in 0..G-1."""
    x = np.arange(genome_length, dtype=np.int64)
    q = np.zeros(genome_length, dtype=np.float64)
    for length, probability in zip(length_support.tolist(), length_pmf.tolist(), strict=True):
        n_starts = genome_length - length + 1
        w = np.minimum(x, n_starts - 1) - np.maximum(0, x - length + 1) + 1
        q += probability * w / n_starts
    return q


def expected_uncovered_fraction(
    genome_length: int, n_reads: int, length_support: IntArray, length_pmf: FloatArray
) -> float:
    """Exact expected uncovered fraction for this simulator (edges and finite size included)."""
    q = start_cover_probability(genome_length, length_support, length_pmf)
    return float(np.mean((1.0 - q) ** n_reads))


def expected_depth_histogram(
    genome_length: int,
    n_reads: int,
    length_support: IntArray,
    length_pmf: FloatArray,
    n_bins: int,
) -> FloatArray:
    """Exact expected fraction of positions at each depth 0..n_bins-1."""
    q = start_cover_probability(genome_length, length_support, length_pmf)
    depths = np.arange(n_bins)[:, None]
    return np.asarray(binom.pmf(depths, n_reads, q[None, :]).mean(axis=1), dtype=np.float64)


def lw_uncovered_fraction(coverage: float) -> float:
    """Lander-Waterman: exp(-c)."""
    return float(np.exp(-coverage))


def lw_depth_pmf(coverage: float, n_bins: int) -> FloatArray:
    """Lander-Waterman: Poisson(c) probabilities for depth 0..n_bins-1."""
    return np.asarray(poisson.pmf(np.arange(n_bins), coverage), dtype=np.float64)


def expected_contigs_exact(n_starts: int, n_reads: int, max_gap: int) -> float:
    """Exact E[#true contigs] for fixed read length (see module docstring for the derivation).

    ``n_starts`` = G - l + 1, ``max_gap`` = l - theta.
    """
    if max_gap < 0:
        raise ValueError("max_gap must be non-negative")
    u = np.arange(n_starts, dtype=np.float64)
    m = np.minimum(float(max_gap), u)
    terms = (1.0 - m / n_starts) ** n_reads - (1.0 - (m + 1.0) / n_starts) ** n_reads
    return float(terms.sum())


def lw_expected_contigs(n_reads: int, coverage: float, min_overlap: int, read_length: int) -> float:
    """Lander-Waterman approximation N * exp(-c * (1 - theta / l))."""
    return float(n_reads * np.exp(-coverage * (1.0 - min_overlap / read_length)))
