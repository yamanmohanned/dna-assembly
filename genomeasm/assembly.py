"""From an order of reads to assembled sequence: layout and consensus.

Plain-language story: once the reads are in order, slide each one along so that it sits on top of
the previous one by exactly their overlap. That picture is the *layout*. Where several reads lie on
top of each other, a vote decides each letter: the most common letter wins (the *consensus*), which
repairs isolated typing mistakes. Where two neighbours do not overlap at all, the chain breaks and
a new *contig* (a gap-free piece of the assembly) begins.

Precise definitions:

* Layout of an order ``pi``: read ``pi_1`` starts at 0 and
  ``offset(pi_{k+1}) = offset(pi_k) + len(pi_k) - ov(pi_k, pi_{k+1})``.
* A contig ends where ``ov(pi_k, pi_{k+1}) = 0``.
* Consensus: for each column of a contig, count the A, C, G, T of all reads covering it and take
  the most frequent; **ties go to the first of A, C, G, T**. Only reads of the same contig vote.
* Total length of all contigs equals ``L(pi)`` exactly (consensus never changes a length).
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from genomeasm.objective import as_order

IntArray = NDArray[np.int64]
IntMatrix = NDArray[np.int32]

_ENCODE = np.zeros(256, dtype=np.uint8)
for _code, _letter in enumerate("ACGT"):
    _ENCODE[ord(_letter)] = _code
_DECODE = np.frombuffer(b"ACGT", dtype=np.uint8)


@dataclass(frozen=True)
class Contig:
    """One gap-free piece of the assembly."""

    sequence: str
    reads: list[int]  # read indices in layout order
    offsets: list[int]  # start column of each read inside the contig


@dataclass(frozen=True)
class Assembly:
    """The contigs produced from one order of reads."""

    contigs: list[Contig]

    @property
    def n_contigs(self) -> int:
        return len(self.contigs)

    @property
    def total_length(self) -> int:
        return sum(len(c.sequence) for c in self.contigs)

    @property
    def sequences(self) -> list[str]:
        return [c.sequence for c in self.contigs]


def merged_string(order: Sequence[int] | IntArray, reads: Sequence[str], ov: IntMatrix) -> str:
    """Glue the reads in order, appending only the letters not already covered by the overlap.

    This is the plain merge used to check ``L(pi)``; it assumes the overlaps are exact. Where the
    overlap is 0 the next read is simply appended.
    """
    arr = as_order(order, len(reads))
    if arr.size == 0:
        return ""
    parts = [reads[int(arr[0])]]
    for a, b in zip(arr[:-1], arr[1:], strict=True):
        parts.append(reads[int(b)][int(ov[a, b]) :])
    return "".join(parts)


def layout(
    order: Sequence[int] | IntArray, lengths: IntArray, ov: IntMatrix
) -> tuple[list[list[int]], list[list[int]]]:
    """Split the order into contigs; returns (reads per contig, offsets per contig)."""
    arr = as_order(order, ov.shape[0])
    groups: list[list[int]] = []
    offsets: list[list[int]] = []
    for position, read in enumerate(arr.tolist()):
        previous = int(arr[position - 1]) if position else -1
        if position == 0 or ov[previous, read] == 0:
            groups.append([read])
            offsets.append([0])
        else:
            groups[-1].append(read)
            offsets[-1].append(offsets[-1][-1] + int(lengths[previous]) - int(ov[previous, read]))
    return groups, offsets


def consensus(reads: Sequence[str], members: Sequence[int], offsets: Sequence[int]) -> str:
    """Majority vote per column over one contig's reads (ties: A before C before G before T)."""
    width = max(off + len(reads[m]) for m, off in zip(members, offsets, strict=True))
    counts = np.zeros((4, width), dtype=np.int32)
    for member, offset in zip(members, offsets, strict=True):
        codes = _ENCODE[np.frombuffer(reads[member].encode("ascii"), dtype=np.uint8)]
        counts[codes, offset + np.arange(codes.size)] += 1
    return _DECODE[np.argmax(counts, axis=0)].tobytes().decode("ascii")


def assemble(
    order: Sequence[int] | IntArray, reads: Sequence[str], lengths: IntArray, ov: IntMatrix
) -> Assembly:
    """Layout and consensus of an order of reads."""
    groups, offsets = layout(order, lengths, ov)
    contigs = [
        Contig(sequence=consensus(reads, members, offs), reads=members, offsets=offs)
        for members, offs in zip(groups, offsets, strict=True)
    ]
    return Assembly(contigs=contigs)
