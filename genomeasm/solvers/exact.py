"""Exact solvers used only to validate the heuristics: brute force (n <= 9) and Held-Karp (n <= 15).

Held-Karp (dynamic programming over subsets; the citation is TODO-VERIFY) computes, for every
subset S of reads and every read j in S, the best fitness of a chain that visits exactly S and ends
in j: ``best[S, j] = max_i best[S - {j}, i] + ov(i, j)``. The answer is the best entry for the full
set. It needs ``O(2^n n^2)`` time and ``O(2^n n)`` memory, hence n <= 15 here.
"""

from __future__ import annotations

import itertools

import numba
import numpy as np
from numpy.typing import NDArray

from genomeasm.solvers.common import Problem, Recorder, SolveResult, tour_of

MAX_BRUTE_FORCE = 9
MAX_HELD_KARP = 15


def _result(problem: Problem, order: NDArray[np.int64], evaluations: int, name: str) -> SolveResult:
    recorder = Recorder(problem)
    return recorder.result(name, tour_of(order), evaluations)


def brute_force(problem: Problem) -> SolveResult:
    """Try every order of the reads (n <= 9) and keep the best fitness."""
    n = problem.n
    if n > MAX_BRUTE_FORCE:
        raise ValueError(f"brute force is limited to {MAX_BRUTE_FORCE} reads")
    w = problem.weights[1:, 1:].astype(np.int64)
    best, best_order, evaluations = -1, tuple(range(n)), 0
    for perm in itertools.permutations(range(n)):
        value = sum(int(w[a, b]) for a, b in itertools.pairwise(perm))
        evaluations += 1
        if value > best:
            best, best_order = value, perm
    return _result(problem, np.array(best_order, dtype=np.int64), evaluations, "brute_force")


@numba.njit(cache=True)
def _held_karp(w: NDArray[np.int64], n: int) -> tuple[NDArray[np.int64], NDArray[np.int8]]:
    full = 1 << n
    best = np.full((full, n), -1, dtype=np.int64)
    parent = np.full((full, n), -1, dtype=np.int8)
    for j in range(n):
        best[1 << j, j] = 0
    for mask in range(1, full):
        for j in range(n):
            if not (mask >> j) & 1 or best[mask, j] < 0:
                continue
            for k in range(n):
                if (mask >> k) & 1:
                    continue
                new_mask = mask | (1 << k)
                value = best[mask, j] + w[j, k]
                if value > best[new_mask, k]:
                    best[new_mask, k] = value
                    parent[new_mask, k] = j
    return best, parent


def held_karp(problem: Problem) -> SolveResult:
    """Optimal order by dynamic programming (n <= 15)."""
    n = problem.n
    if n > MAX_HELD_KARP:
        raise ValueError(f"Held-Karp is limited to {MAX_HELD_KARP} reads")
    if n == 1:
        return _result(problem, np.zeros(1, dtype=np.int64), 1, "held_karp")
    w = problem.weights[1:, 1:].astype(np.int64)
    best, parent = _held_karp(w, n)
    full = (1 << n) - 1
    last = int(np.argmax(best[full]))
    order = [last]
    mask = full
    while parent[mask, order[-1]] >= 0:
        prev = int(parent[mask, order[-1]])
        mask ^= 1 << order[-1]
        order.append(prev)
    order.reverse()
    return _result(problem, np.array(order, dtype=np.int64), (1 << n) * n * n, "held_karp")
