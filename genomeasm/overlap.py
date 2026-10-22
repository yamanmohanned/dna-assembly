"""Overlaps between reads: the overlap graph that turns assembly into a travelling-salesman problem.

Plain-language story: two reads can be glued when the end of one equals the beginning of the other,
like two jigsaw pieces whose edges match. The longer the matching edge, the more sure we are that
they were neighbours in the genome. The overlap of read i followed by read j is the length of the
longest such matching edge. Gluing a chain of reads saves the overlap each time, so the best
assembly is the order of reads that saves the most letters in total.

Definitions (master prompt, "Mathematical specification"):

* Reads are first made *substring-free*: duplicates are dropped and a read that is contained in
  another read is dropped (it adds no new letters).
* ``ov(i, j)`` is the largest ``k`` with ``tau <= k < min(len_i, len_j)`` such that the last ``k``
  letters of read ``i`` equal the first ``k`` letters of read ``j``; otherwise 0. It is
  *directional*: ``ov(i, j) != ov(j, i)`` in general. ``ov(i, i)`` is defined as 0.
* Error-tolerant overlap: the largest ``k`` in the same range such that the number of positions
  where the suffix of ``i`` and the prefix of ``j`` differ (Hamming distance, no insertions or
  deletions) is at most ``floor(tolerance * k + 1e-9)``. The **score is the length k itself**:
  mismatches are allowed but not charged. The 1e-9 only protects against floating-point error
  such as ``0.29 * 100 = 28.999999999999996``. With ``tolerance = 0`` it is the exact overlap.
* ``tau`` is the minimum overlap. Too small and unrelated reads match by chance; too large and
  true neighbours are missed. ``recommend_min_overlap`` chooses it from a union bound (below).

Spurious overlaps in a random-sequence model (decision D-027): for two unrelated reads, the chance
that some ``k >= tau`` gives a match is at most ``sum_{k >= tau} 4**-k < (4/3) * 4**-tau``. Over
``n * (n - 1)`` ordered pairs, the expected number of spurious overlaps is therefore below
``n * (n - 1) * (4/3) * 4**-tau``. Real genomes contain repeats, so this is a guide, not a
guarantee; the measured numbers are in results/phase3/.
"""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numba
import numpy as np
from numpy.typing import NDArray

IntMatrix = NDArray[np.int32]
IntArray = NDArray[np.int64]

_ROUNDING_GUARD = 1e-9


def mismatch_limit(k: int, tolerance: float) -> int:
    """Largest number of mismatches allowed in an overlap of length ``k``: floor(tolerance * k)."""
    return int(math.floor(tolerance * k + _ROUNDING_GUARD))


# --------------------------------------------------------------------------------------
# Making the reads substring-free
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class Preprocessed:
    """Reads after removing duplicates and contained reads, with the bookkeeping."""

    reads: list[str]  # the kept reads, in their original order
    kept: IntArray  # index of each kept read in the original list
    n_input: int
    n_duplicates: int  # identical to an earlier read
    n_contained: int  # a (shorter) distinct read that is a substring of another read


def make_substring_free(sequences: Sequence[str]) -> Preprocessed:
    """Drop duplicate reads, then reads that are substrings of another read.

    Of several identical reads the first is kept. A read contained in a longer read is dropped;
    chains (a inside b inside c) leave only c. The result has no read inside another.
    """
    n_input = len(sequences)
    if n_input == 0:
        return Preprocessed([], np.zeros(0, dtype=np.int64), 0, 0, 0)
    first_seen: set[str] = set()
    unique: list[int] = []
    for i, sequence in enumerate(sequences):
        if sequence not in first_seen:
            first_seen.add(sequence)
            unique.append(i)
    n_duplicates = n_input - len(unique)

    texts = [sequences[i] for i in unique]
    shortest = min(len(t) for t in texts)
    # Only a read longer than the shortest one can contain another distinct read. Index every
    # window of that shortest length: if r is inside s then r[:shortest] is a window of s.
    windows: dict[str, list[int]] = defaultdict(list)
    for idx, text in enumerate(texts):
        if len(text) > shortest:
            for window in {text[p : p + shortest] for p in range(len(text) - shortest + 1)}:
                windows[window].append(idx)
    keep: list[int] = []
    for idx, text in enumerate(texts):
        contained = any(
            other != idx and len(texts[other]) > len(text) and text in texts[other]
            for other in windows.get(text[:shortest], ())
        )
        if not contained:
            keep.append(idx)
    kept = np.array([unique[i] for i in keep], dtype=np.int64)
    return Preprocessed(
        reads=[texts[i] for i in keep],
        kept=kept,
        n_input=n_input,
        n_duplicates=n_duplicates,
        n_contained=len(unique) - len(keep),
    )


def is_substring_free(reads: Sequence[str]) -> bool:
    """True if no read is identical to, or contained in, another read (slow reference check)."""
    for i, a in enumerate(reads):
        for j, b in enumerate(reads):
            if i != j and a in b:
                return False
    return True


# --------------------------------------------------------------------------------------
# Exact overlaps
# --------------------------------------------------------------------------------------


def overlap_matrix_naive(reads: Sequence[str], min_overlap: int) -> IntMatrix:
    """Reference implementation straight from the definition (quadratic, only for testing)."""
    _check_tau(min_overlap)
    n = len(reads)
    ov = np.zeros((n, n), dtype=np.int32)
    for i, a in enumerate(reads):
        for j, b in enumerate(reads):
            if i == j:
                continue
            for k in range(min(len(a), len(b)) - 1, min_overlap - 1, -1):
                if a[len(a) - k :] == b[:k]:
                    ov[i, j] = k
                    break
    return ov


def overlap_matrix(reads: Sequence[str], min_overlap: int) -> IntMatrix:
    """Exact suffix-prefix overlaps of all ordered pairs, using a prefix index.

    Every read ``j`` is filed under its first ``tau`` letters. For each read ``i`` and each suffix
    of length ``k >= tau`` (longest first) only the reads filed under the first ``tau`` letters of
    that suffix are compared, so the cost is about ``n * length`` look-ups instead of ``n**2``.
    """
    _check_tau(min_overlap)
    tau = min_overlap
    n = len(reads)
    ov = np.zeros((n, n), dtype=np.int32)
    by_prefix: dict[str, list[int]] = defaultdict(list)
    for j, read in enumerate(reads):
        if len(read) > tau:  # needs k >= tau and k < len, so len >= tau + 1
            by_prefix[read[:tau]].append(j)
    for i, a in enumerate(reads):
        length = len(a)
        for start in range(1, length - tau + 1):  # suffix a[start:], longest first
            candidates = by_prefix.get(a[start : start + tau])
            if not candidates:
                continue
            k = length - start
            suffix = a[start:]
            for j in candidates:
                if j != i and ov[i, j] == 0 and k < len(reads[j]) and reads[j].startswith(suffix):
                    ov[i, j] = k
    return ov


# --------------------------------------------------------------------------------------
# Error-tolerant overlaps
# --------------------------------------------------------------------------------------


@numba.njit(parallel=True, cache=True)
def _tolerant_kernel(
    flat: NDArray[np.uint8],
    starts: NDArray[np.int64],
    lengths: NDArray[np.int64],
    tau: int,
    tolerance: float,
    out: NDArray[np.int32],
) -> None:
    n = lengths.shape[0]
    for i in numba.prange(n):
        li = lengths[i]
        oi = starts[i]
        for j in range(n):
            if i == j:
                continue
            lj = lengths[j]
            oj = starts[j]
            kmax = min(li, lj) - 1
            for k in range(kmax, tau - 1, -1):
                limit = int(math.floor(tolerance * k + 1e-9))
                mismatches = 0
                base = oi + li - k
                ok = True
                for p in range(k):
                    if flat[base + p] != flat[oj + p]:
                        mismatches += 1
                        if mismatches > limit:
                            ok = False
                            break
                if ok:
                    out[i, j] = k
                    break


def overlap_matrix_tolerant(reads: Sequence[str], min_overlap: int, tolerance: float) -> IntMatrix:
    """Error-tolerant overlaps (Hamming distance at most floor(tolerance * k)), exact by search.

    Every ordered pair and every length is tested, with early exit on too many mismatches, so no
    overlap that meets the definition is missed (unlike seed heuristics). The cost is quadratic in
    the number of reads; it is meant for the noise experiments, with a few thousand reads at most.
    """
    _check_tau(min_overlap)
    if not 0.0 <= tolerance < 1.0:
        raise ValueError("tolerance must lie in [0, 1)")
    n = len(reads)
    out = np.zeros((n, n), dtype=np.int32)
    if n == 0:
        return out
    lengths = np.array([len(r) for r in reads], dtype=np.int64)
    starts = np.concatenate([[0], np.cumsum(lengths)[:-1]]).astype(np.int64)
    flat = np.frombuffer("".join(reads).encode("ascii"), dtype=np.uint8).copy()
    _tolerant_kernel(flat, starts, lengths, min_overlap, float(tolerance), out)
    return out


# --------------------------------------------------------------------------------------
# Choosing tau
# --------------------------------------------------------------------------------------


def expected_spurious_pairs(n_reads: int, min_overlap: int) -> float:
    """Upper bound on the expected number of chance overlaps among unrelated random reads."""
    return n_reads * (n_reads - 1) * (4.0 / 3.0) * 4.0**-min_overlap


def recommend_min_overlap(
    n_reads: int, *, spurious_budget: float = 1.0, minimum: int = 8, maximum: int | None = None
) -> int:
    """Smallest tau for which fewer than ``spurious_budget`` chance overlaps are expected.

    A *smaller* tau keeps more true overlaps (fewer contigs), so the smallest tau that meets the
    budget is chosen. ``minimum`` guards against tiny instances.
    """
    if n_reads < 2:
        return minimum
    bound = n_reads * (n_reads - 1) * (4.0 / 3.0) / spurious_budget
    tau = max(minimum, math.ceil(math.log(bound, 4)))
    return tau if maximum is None else min(tau, maximum)


def _check_tau(min_overlap: int) -> None:
    if min_overlap < 1:
        raise ValueError("min_overlap must be at least 1")


# --------------------------------------------------------------------------------------
# Comparing an overlap matrix with the ground truth of simulated reads
# --------------------------------------------------------------------------------------


def truth_agreement(
    ov: IntMatrix, starts: IntArray, spans: IntArray, min_overlap: int, *, chunk: int = 400
) -> dict[str, Any]:
    """Compare overlaps with the true positions of the reads (simulator ground truth).

    For reads ``i`` and ``j`` that start at ``s_i < s_j`` and end at ``e_i`` (exclusive), the true
    directed overlap is ``t = e_i - s_j`` when positive and the read ``j`` is not inside ``i``.
    Counted over ordered pairs:

    * ``true_pairs``: pairs with ``t >= tau`` and ``t < min(span_i, span_j)``.
    * ``recovered``: true pairs with ``ov >= t`` (the true overlap, or a longer chance one).
    * ``exact``: true pairs with ``ov == t``.
    * ``longer``: true pairs with ``ov > t`` (a longer chance match hides the true overlap).
    * ``edges``: pairs with ``ov > 0``; ``spurious``: edges whose pair has no true overlap of at
      least ``tau`` (``t < tau`` or not a true neighbour).
    * ``hist_true`` / ``hist_other``: how many edges of each overlap length belong to true pairs
      and to other pairs (arrays indexed by overlap length).
    For error-free reads ``recovered == true_pairs`` must hold, because the true overlap is always
    a valid suffix-prefix match.
    """
    n = ov.shape[0]
    ends = starts + spans
    counts: dict[str, Any] = {
        "true_pairs": 0,
        "recovered": 0,
        "exact": 0,
        "longer": 0,
        "edges": 0,
        "spurious": 0,
    }
    width = int(np.max(ov)) + 1 if ov.size else 1
    hist_true = np.zeros(width, dtype=np.int64)  # overlap lengths of edges that are true pairs
    hist_other = np.zeros(width, dtype=np.int64)  # overlap lengths of all other edges
    for lo in range(0, n, chunk):
        hi = min(n, lo + chunk)
        t = ends[lo:hi, None] - starts[None, :]
        valid_geometry = (starts[None, :] > starts[lo:hi, None]) & (t > 0)
        limit = np.minimum(spans[lo:hi, None], spans[None, :])
        true_pair = valid_geometry & (t >= min_overlap) & (t < limit)
        block = ov[lo:hi]
        counts["true_pairs"] += int(true_pair.sum())
        counts["recovered"] += int((true_pair & (block >= t)).sum())
        counts["exact"] += int((true_pair & (block == t)).sum())
        counts["longer"] += int((true_pair & (block > t)).sum())
        edge = block > 0
        counts["edges"] += int(edge.sum())
        counts["spurious"] += int((edge & ~true_pair).sum())
        hist_true += np.bincount(block[edge & true_pair], minlength=width)
        hist_other += np.bincount(block[edge & ~true_pair], minlength=width)
    counts["hist_true"] = hist_true
    counts["hist_other"] = hist_other
    return counts
