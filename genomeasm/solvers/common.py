"""Shared types of all solvers: the problem, the budget and the result.

Every solver works on the *cycle form* of the problem: nodes 0..n where node 0 is the dummy and
read i is node i + 1. ``Problem.weights`` is the (n+1) x (n+1) matrix W with W[a, b] = ov(a, b)
for reads and 0 in row 0 and column 0. A tour is an array ``t`` of length n + 1 with ``t[0] = 0``;
its fitness is ``F(t) = sum_p W[t[p], t[(p + 1) % (n + 1)]]``, which equals the objective F of the
order ``t[1:] - 1`` (the dummy edges add 0). The matching ATSP cost of the cycle is
``total_length - F``.

Budgets (rule R6): all heuristics get the same wall-clock time and every solver also reports the
number of *evaluations*, defined as the number of times the objective change of one candidate
move was computed (one delta evaluation, one insertion-gain evaluation, or one full recomputation
counts 1). LKH-3 does not report evaluations.
"""

from __future__ import annotations

import time
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from genomeasm.instance import AssemblyInstance

IntArray = NDArray[np.int64]


@dataclass(frozen=True)
class Problem:
    """An assembly instance in the cycle form used by the solvers."""

    weights: NDArray[np.int32]
    lengths: IntArray  # read lengths, needed for the ATSP cost matrix (LKH-3)

    @property
    def n(self) -> int:
        return int(self.weights.shape[0]) - 1

    @property
    def total_length(self) -> int:
        return int(self.lengths.sum())

    @classmethod
    def from_instance(cls, inst: AssemblyInstance) -> Problem:
        return cls.from_matrix(inst.ov, inst.lengths)

    @classmethod
    def from_matrix(cls, ov: NDArray[np.int32], lengths: IntArray) -> Problem:
        """Build a problem from any overlap matrix (also random ones, for testing the solvers)."""
        weights = np.zeros((ov.shape[0] + 1, ov.shape[0] + 1), dtype=np.int32)
        weights[1:, 1:] = ov
        np.fill_diagonal(weights, 0)
        return cls(weights=weights, lengths=np.asarray(lengths, dtype=np.int64))

    def cost_matrix(self) -> IntArray:
        """The ATSP cost matrix of docs/OBJECTIVE_FUNCTION.md (as ``objective.cost_matrix``)."""
        cost = np.zeros_like(self.weights, dtype=np.int64)
        cost[0, 1:] = self.lengths
        cost[1:, 1:] = self.lengths[None, :] - self.weights[1:, 1:]
        np.fill_diagonal(cost, 0)
        return cost

    def fitness(self, order: IntArray | list[int]) -> int:
        """F of an order of reads 0..n-1 (the same value as ``objective.fitness``)."""
        t = tour_of(np.asarray(order, dtype=np.int64))
        return int(self.weights[t, np.roll(t, -1)].sum(dtype=np.int64))


def tour_of(order: IntArray) -> IntArray:
    """Tour array [0, order + 1]."""
    return np.concatenate([[0], np.asarray(order, dtype=np.int64) + 1])


def order_of(tour: IntArray) -> IntArray:
    """Order of reads 0..n-1 from a tour that starts at the dummy node."""
    return np.asarray(tour[1:], dtype=np.int64) - 1


@dataclass(frozen=True)
class Budget:
    """Stop when the time or the number of evaluations is used up (at least one must be set)."""

    seconds: float | None = None
    evaluations: int | None = None

    def __post_init__(self) -> None:
        if self.seconds is None and self.evaluations is None:
            raise ValueError("a budget needs seconds, evaluations or both")

    def used_up(self, elapsed: float, evaluations: int) -> bool:
        return (self.seconds is not None and elapsed >= self.seconds) or (
            self.evaluations is not None and evaluations >= self.evaluations
        )

    def fraction(self, elapsed: float, evaluations: int) -> float:
        """How much of the budget is used, between 0 and 1 (the larger of the two measures)."""
        parts = []
        if self.seconds is not None:
            parts.append(elapsed / self.seconds)
        if self.evaluations is not None:
            parts.append(evaluations / self.evaluations)
        return float(min(1.0, max(parts)))


@dataclass
class SolveResult:
    """What a solver returns."""

    solver: str
    order: IntArray
    fitness: int
    merged_length: int
    evaluations: int
    seconds: float
    trace: list[tuple[float, int, int]] = field(default_factory=list)  # (seconds, evals, best F)
    info: dict[str, Any] = field(default_factory=dict)


class Recorder:
    """Keeps the best-so-far trace of a run (for convergence curves)."""

    def __init__(self, problem: Problem, clock: Callable[[], float] = time.perf_counter) -> None:
        self.problem = problem
        self.clock = clock
        self.start = clock()
        self.trace: list[tuple[float, int, int]] = []
        self.best = -1

    def elapsed(self) -> float:
        return self.clock() - self.start

    def note(self, best_fitness: int, evaluations: int) -> None:
        if best_fitness > self.best:
            self.best = best_fitness
            self.trace.append((self.elapsed(), evaluations, best_fitness))

    def result(
        self, solver: str, tour: IntArray, evaluations: int, info: dict[str, Any] | None = None
    ) -> SolveResult:
        order = order_of(tour)
        fit = self.problem.fitness(order)
        self.note(fit, evaluations)
        return SolveResult(
            solver=solver,
            order=order,
            fitness=fit,
            merged_length=self.problem.total_length - fit,
            evaluations=evaluations,
            seconds=self.elapsed(),
            trace=self.trace,
            info=info or {},
        )
