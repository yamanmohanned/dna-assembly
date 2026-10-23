"""Constructive baselines: global greedy edge merging and nearest neighbour.

Greedy (the classical assembly heuristic): sort all overlaps from largest to smallest and accept an
overlap ``i -> j`` when read ``i`` still has no successor, read ``j`` still has no predecessor and
the edge does not close a cycle. The result is a set of chains; they are joined in the order of
their first read (joins carry no overlap, otherwise greedy would have taken them). Ties between
equal overlaps are broken by the smaller ``(i, j)``, so the result is deterministic.

Nearest neighbour: start at the read with the smallest best incoming overlap (the most likely
start of the genome), then repeatedly go to the unused read with the largest overlap from the
current one (smallest index on ties); when no unused read overlaps, continue with the unused read
of smallest index.

Both are deterministic and finish in a fraction of a second, so they do not use their budget.
"""

from __future__ import annotations

import numba
import numpy as np
from numpy.typing import NDArray

from genomeasm.solvers.common import Problem, Recorder, SolveResult, tour_of


@numba.njit(cache=True)
def _greedy_chains(
    src: NDArray[np.int64], dst: NDArray[np.int64], n: int
) -> tuple[NDArray[np.int64], NDArray[np.int64]]:
    succ = np.full(n, -1, dtype=np.int64)
    pred = np.full(n, -1, dtype=np.int64)
    head = np.arange(n, dtype=np.int64)  # first read of the chain containing each read (at tails)
    tail = np.arange(n, dtype=np.int64)  # last read of the chain containing each read (at heads)
    for e in range(src.shape[0]):
        i, j = src[e], dst[e]
        if succ[i] != -1 or pred[j] != -1:
            continue
        if head[i] == j:  # j is the head of the chain that ends in i: this would close a cycle
            continue
        h, t = head[i], tail[j]
        succ[i] = j
        pred[j] = i
        tail[h] = t
        head[t] = h
    return succ, pred


def greedy(problem: Problem) -> SolveResult:
    """Global greedy edge merging by largest overlap."""
    recorder = Recorder(problem)
    n = problem.n
    ov = problem.weights[1:, 1:]
    src, dst = np.nonzero(ov > 0)
    keys = np.lexsort((dst, src, -ov[src, dst].astype(np.int64)))  # largest first, then (i, j)
    succ, pred = _greedy_chains(src[keys].astype(np.int64), dst[keys].astype(np.int64), n)
    order: list[int] = []
    for start in range(n):
        if pred[start] == -1:  # a chain begins here
            node = start
            while node != -1:
                order.append(node)
                node = int(succ[node])
    return recorder.result("greedy", tour_of(np.array(order, dtype=np.int64)), int(src.size))


@numba.njit(cache=True)
def _nearest_neighbour(ov: NDArray[np.int32]) -> NDArray[np.int64]:
    n = ov.shape[0]
    used = np.zeros(n, dtype=np.bool_)
    order = np.empty(n, dtype=np.int64)
    best_in = np.zeros(n, dtype=np.int64)
    for j in range(n):
        for i in range(n):
            if i != j and ov[i, j] > best_in[j]:
                best_in[j] = ov[i, j]
    cur = 0
    for j in range(n):
        if best_in[j] < best_in[cur]:
            cur = j
    for step in range(n):
        order[step] = cur
        used[cur] = True
        nxt, best = -1, 0
        for j in range(n):
            if not used[j] and ov[cur, j] > best:
                nxt, best = j, ov[cur, j]
        if nxt == -1:
            for j in range(n):
                if not used[j]:
                    nxt = j
                    break
        cur = nxt
    return order


def nearest_neighbour(problem: Problem) -> SolveResult:
    """Nearest-neighbour construction (secondary baseline)."""
    recorder = Recorder(problem)
    n = problem.n
    order = _nearest_neighbour(np.ascontiguousarray(problem.weights[1:, 1:]))
    return recorder.result("nearest_neighbour", tour_of(order), n * n)
