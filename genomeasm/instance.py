"""An assembly instance: substring-free reads, their overlaps, and everything derived from them.

This is the single object that the solvers (Phase 4) and the experiments (Phase 5) work with. It
joins the pieces built in this phase: preprocessing (``overlap.make_substring_free``), the overlap
matrix, the objective ``F`` / ``L``, the ATSP cost matrix with the dummy node, and layout +
consensus. When the reads come from the simulator it also keeps their ground truth, so the *true
ordering* (reads sorted by their real start position) and its fitness are available as a
reference.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

from genomeasm import objective
from genomeasm.assembly import Assembly, assemble
from genomeasm.overlap import (
    make_substring_free,
    overlap_matrix,
    overlap_matrix_tolerant,
    recommend_min_overlap,
)
from genomeasm.simulate import ReadSet

IntArray = NDArray[np.int64]
IntMatrix = NDArray[np.int32]


@dataclass(frozen=True)
class AssemblyInstance:
    """Substring-free reads with their overlap matrix, plus optional ground truth."""

    reads: list[str]
    lengths: IntArray
    ov: IntMatrix
    min_overlap: int
    tolerance: float
    kept: IntArray  # index of each read in the list the instance was built from
    n_input: int
    n_duplicates: int
    n_contained: int
    true_starts: IntArray | None = None  # ground truth of the kept reads (simulated data only)
    true_spans: IntArray | None = None

    @classmethod
    def from_sequences(
        cls,
        sequences: Sequence[str],
        *,
        min_overlap: int | str = "auto",
        tolerance: float = 0.0,
        starts: IntArray | None = None,
        spans: IntArray | None = None,
    ) -> AssemblyInstance:
        """Build an instance from raw read strings (optionally with their true positions)."""
        pre = make_substring_free(sequences)
        tau = recommend_min_overlap(len(pre.reads)) if min_overlap == "auto" else int(min_overlap)
        if tolerance > 0:
            ov = overlap_matrix_tolerant(pre.reads, tau, tolerance)
        else:
            ov = overlap_matrix(pre.reads, tau)
        return cls(
            reads=pre.reads,
            lengths=np.array([len(r) for r in pre.reads], dtype=np.int64),
            ov=ov,
            min_overlap=tau,
            tolerance=tolerance,
            kept=pre.kept,
            n_input=pre.n_input,
            n_duplicates=pre.n_duplicates,
            n_contained=pre.n_contained,
            true_starts=None if starts is None else starts[pre.kept],
            true_spans=None if spans is None else spans[pre.kept],
        )

    @classmethod
    def from_readset(
        cls, reads: ReadSet, *, min_overlap: int | str = "auto", tolerance: float = 0.0
    ) -> AssemblyInstance:
        """Build an instance from simulated reads, keeping their ground truth."""
        return cls.from_sequences(
            reads.sequences,
            min_overlap=min_overlap,
            tolerance=tolerance,
            starts=reads.starts,
            spans=reads.spans,
        )

    @property
    def n(self) -> int:
        return len(self.reads)

    def cost(self) -> IntArray:
        """The (n+1) x (n+1) ATSP cost matrix with the dummy node 0."""
        return objective.cost_matrix(self.ov, self.lengths)

    def fitness(self, order: Sequence[int] | IntArray) -> int:
        """F(pi), total overlap (maximise)."""
        return objective.fitness(order, self.ov)

    def merged_length(self, order: Sequence[int] | IntArray) -> int:
        """L(pi) = sum of lengths minus F(pi) (minimise)."""
        return objective.merged_length(order, self.lengths, self.ov)

    def assemble(self, order: Sequence[int] | IntArray) -> Assembly:
        """Layout and consensus for an order."""
        return assemble(order, self.reads, self.lengths, self.ov)

    def true_order(self) -> IntArray:
        """Reads sorted by their real start position (needs simulated ground truth)."""
        if self.true_starts is None:
            raise ValueError("this instance has no ground truth")
        return np.argsort(self.true_starts, kind="stable").astype(np.int64)

    def write_tsplib(self, path: Path, name: str) -> None:
        """Write the ATSP instance for LKH-3 or any other TSPLIB solver."""
        comment = (
            f"DNA fragment assembly: {self.n} reads + dummy node 0, "
            f"min_overlap={self.min_overlap}, tolerance={self.tolerance}"
        )
        objective.write_tsplib_atsp(path, self.cost(), name=name, comment=comment)

    def summary(self) -> dict[str, Any]:
        """Counts that every report of an instance needs."""
        return {
            "n_input_reads": self.n_input,
            "n_duplicates_removed": self.n_duplicates,
            "n_contained_removed": self.n_contained,
            "n_reads": self.n,
            "min_overlap": self.min_overlap,
            "tolerance": self.tolerance,
            "total_read_length": int(self.lengths.sum()),
            "n_overlap_edges": int((self.ov > 0).sum()),
        }
