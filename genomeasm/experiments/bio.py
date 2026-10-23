"""Biological quality of an assembly: does the order of the reads rebuild the real genome?

Plain-language story: a high objective value F only says the reads were glued with long overlaps.
Here we check the result against the answer key. Each assembled piece (contig) is cut into
windows of 100 letters; every window is located in the real genome by shared 15-letter words, and
we count how many letters agree. From that we get the share of the genome that was rebuilt, how
correct the letters are, how many pieces there are, how long the typical piece is (N50) and whether
a piece jumps between distant places of the genome (a misjoin).

Definitions (all on the simulated window of the genome, one strand, no indels in the reads):

* ``identity``: matching letters / aligned letters over all windows that could be located (at least
  3 shared words on one diagonal).
* ``genome_fraction``: share of genome positions covered by a located window with identity >= 0.9.
* ``n_contigs``, ``n50`` (length L such that contigs of length >= L hold half of all assembled
  letters), ``assembled_length`` (= L(pi)), ``largest_contig``.
* ``misjoins``: number of neighbouring windows inside one contig that are located more than 10
  letters apart (the contig jumps), ``unmapped_fraction``: share of windows not located.
* ``spearman``: rank correlation between the position of a read in the order and its true start.
* ``adjacent_correct``: share of neighbouring reads in the order whose true starts differ by at most
  one read length (true neighbours).
"""

from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Sequence
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.stats import spearmanr

from genomeasm.assembly import Assembly

K = 15
WINDOW = 100
MIN_VOTES = 3
MIN_WINDOW_IDENTITY = 0.9
JUMP = 10

IntArray = NDArray[np.int64]


class ReferenceIndex:
    """All 15-letter words of the reference (the real genome window) with their positions."""

    def __init__(self, sequence: str) -> None:
        self.sequence = sequence
        self.codes = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)
        self.words: dict[str, list[int]] = defaultdict(list)
        for p in range(len(sequence) - K + 1):
            self.words[sequence[p : p + K]].append(p)

    def locate(self, window: str) -> tuple[int, int, int] | None:
        """(diagonal, matches, aligned) of the window on the reference, or None if not located."""
        votes: Counter[int] = Counter()
        for offset in range(len(window) - K + 1):
            for position in self.words.get(window[offset : offset + K], ()):
                votes[position - offset] += 1
        if not votes:
            return None
        diagonal, count = votes.most_common(1)[0]
        if count < MIN_VOTES:
            return None
        codes = np.frombuffer(window.encode("ascii"), dtype=np.uint8)
        first = max(0, -diagonal)
        last = min(len(window), len(self.codes) - diagonal)
        if last <= first:
            return None
        matches = int(np.sum(codes[first:last] == self.codes[diagonal + first : diagonal + last]))
        return diagonal, matches, last - first


def n50(lengths: Sequence[int]) -> int:
    """Length L such that contigs of length at least L contain half of all assembled letters."""
    if not lengths:
        return 0
    total, running = sum(lengths), 0
    for length in sorted(lengths, reverse=True):
        running += length
        if 2 * running >= total:
            return int(length)
    return 0


def assembly_quality(assembly: Assembly, reference: ReferenceIndex) -> dict[str, Any]:
    """Genome fraction, identity, contig statistics and misjoins of an assembly."""
    matches = aligned = windows = unmapped = misjoins = 0
    covered = np.zeros(len(reference.sequence), dtype=bool)
    for contig in assembly.contigs:
        sequence = contig.sequence
        previous: int | None = None
        for start in range(0, len(sequence), WINDOW):
            window = sequence[start : start + WINDOW]
            if len(window) < 30:
                continue
            windows += 1
            hit = reference.locate(window)
            if hit is None:
                unmapped += 1
                previous = None
                continue
            diagonal, m, a = hit
            matches, aligned = matches + m, aligned + a
            shift = diagonal - start  # the contig's own diagonal (constant if nothing jumps)
            if previous is not None and abs(shift - previous) > JUMP:
                misjoins += 1
            previous = shift
            if a >= 0.5 * len(window) and m >= MIN_WINDOW_IDENTITY * a:
                lo = max(0, diagonal + max(0, -diagonal))
                covered[lo : min(len(covered), diagonal + len(window))] = True
    lengths = [len(c.sequence) for c in assembly.contigs]
    return {
        "identity": matches / aligned if aligned else 0.0,
        "genome_fraction": float(covered.mean()),
        "n_contigs": len(lengths),
        "n50": n50(lengths),
        "assembled_length": int(sum(lengths)),
        "largest_contig": int(max(lengths)) if lengths else 0,
        "misjoins": misjoins,
        "unmapped_fraction": unmapped / windows if windows else 1.0,
    }


def _rho(starts: IntArray) -> float:
    if starts.size < 3 or np.all(starts == starts[0]):
        return 1.0
    return float(spearmanr(np.arange(starts.size), starts).statistic)


def order_quality(
    order: IntArray, true_starts: IntArray, read_length: int, contigs: list[list[int]]
) -> dict[str, float]:
    """Rank correlations with the true positions and the share of correct neighbours.

    ``spearman`` uses the order as the solver returned it. Contigs are placed in an arbitrary
    order (the edges between them carry no overlap), which lowers it, so ``spearman_contigs`` is
    also given: the same correlation after sorting the contigs by the median true start of
    their reads, which measures only how well reads are ordered inside the contigs.
    """
    starts = true_starts[order]
    adjacent = float(np.mean(np.abs(np.diff(starts)) <= read_length)) if order.size > 1 else 1.0
    arranged = sorted(contigs, key=lambda reads: float(np.median(true_starts[reads])))
    inside = true_starts[np.array([r for reads in arranged for r in reads], dtype=np.int64)]
    return {
        "spearman": _rho(starts),
        "spearman_contigs": _rho(inside),
        "adjacent_correct": adjacent,
    }
