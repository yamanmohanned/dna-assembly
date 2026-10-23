"""Constant-time move evaluation on a tour (swap and insert), plus the move application.

A tour ``t`` has ``t[0] = 0`` (the dummy) and the reads at positions 1..n; positions wrap around
(``t[n + 1]`` is ``t[0]``). With weights ``W`` the fitness is the sum of ``W`` over the cycle edges,
so a move that replaces a few edges changes the fitness by (new edges) - (old edges): O(1).

* ``swap_delta(t, W, i, j)``  exchange the nodes at positions ``i < j`` (both in 1..n).
* ``insert_delta(t, W, i, j)``  take the node at position ``i`` out and put it back right after
  position ``j`` (``j`` in 0..n, ``j != i`` and ``j != i - 1``; ``j = 0`` means first).

Every delta is tested against a full recomputation (tests/test_solvers.py).
"""

from __future__ import annotations

import numba
import numpy as np
from numpy.typing import NDArray

IntArray = NDArray[np.int64]


@numba.njit(cache=True)
def tour_fitness(t: IntArray, w: NDArray[np.int32]) -> int:
    n1 = t.shape[0]
    total = 0
    for p in range(n1):
        total += w[t[p], t[(p + 1) % n1]]
    return total


@numba.njit(cache=True)
def swap_delta(t: IntArray, w: NDArray[np.int32], i: int, j: int) -> int:
    n1 = t.shape[0]
    a, x, b = t[i - 1], t[i], t[(i + 1) % n1]
    c, y, d = t[j - 1], t[j], t[(j + 1) % n1]
    if j == i + 1:
        old = w[a, x] + w[x, y] + w[y, d]
        new = w[a, y] + w[y, x] + w[x, d]
    else:
        old = w[a, x] + w[x, b] + w[c, y] + w[y, d]
        new = w[a, y] + w[y, b] + w[c, x] + w[x, d]
    return new - old


@numba.njit(cache=True)
def insert_delta(t: IntArray, w: NDArray[np.int32], i: int, j: int) -> int:
    n1 = t.shape[0]
    p, x, q = t[i - 1], t[i], t[(i + 1) % n1]
    a, b = t[j], t[(j + 1) % n1]
    return (w[p, q] - w[p, x] - w[x, q]) + (w[a, x] + w[x, b] - w[a, b])


@numba.njit(cache=True)
def apply_swap(t: IntArray, i: int, j: int) -> None:
    t[i], t[j] = t[j], t[i]


@numba.njit(cache=True)
def apply_insert(t: IntArray, i: int, j: int) -> None:
    x = t[i]
    if j > i:
        for k in range(i, j):
            t[k] = t[k + 1]
        t[j] = x
    else:  # j < i - 1: the node ends up at position j + 1
        for k in range(i, j + 1, -1):
            t[k] = t[k - 1]
        t[j + 1] = x
