"""Tabu search on the permutation of reads (swap and insert neighbourhoods, O(1) delta evaluation).

Plain-language story: stand on a solution, look at many small changes (swap two reads, or move one
read somewhere else), step to the best one even if it is worse than where you stand, but forbid
moving the same reads again for a while (the *tabu list*) so you do not walk in circles. A forbidden
move is still allowed if it would give a new overall best (*aspiration*). If nothing improves for a
long time, shake the best solution and start again (*restart*).

Parameters (all exposed; tuned in Phase 5 on separate tuning instances):

* ``tenure``: base tabu tenure; each moved read is tabu for a random time in
  ``[ceil(0.5 * tenure), ceil(1.5 * tenure)]`` iterations. Default ``max(5, round(sqrt(n)))``.
* ``moves_per_iteration``: how many random moves are evaluated per iteration; 0 means the full
  neighbourhood (all swaps and all inserts). Default 0 for ``n <= 60``, else ``20 * n``.
* ``restart_after``: iterations without a new best before a restart. Default ``max(200, 10 * n)``.
* ``kick``: number of random swaps applied to the best tour at a restart. Default
  ``max(3, round(0.05 * n))``.
* ``init``: ``"greedy"`` (default) or ``"random"``.
* ``candidate_fraction`` and ``candidates``: with probability ``candidate_fraction`` a sampled
  move is a *candidate-list insert* (a read is put right after one of its ``candidates`` best
  predecessors or right before one of its best successors, by overlap) instead of a uniformly
  random move. Default 0 (the plain random-sampling tabu search, "tabu" in the experiments);
  "tabu_cl" uses them.

Evaluations: one per delta computed (full recomputations at start and restart count 1 each).
"""

from __future__ import annotations

import math

import numba
import numpy as np
from numpy.typing import NDArray

from genomeasm.solvers.common import Budget, IntArray, Problem, Recorder, SolveResult, tour_of
from genomeasm.solvers.greedy import greedy
from genomeasm.solvers.moves import (
    apply_insert,
    apply_swap,
    insert_delta,
    swap_delta,
    tour_fitness,
)


@numba.njit(cache=True)
def seed_numba(seed: int) -> None:
    np.random.seed(seed)


@numba.njit(cache=True)
def _tabu_chunk(
    t: IntArray,
    w: NDArray[np.int32],
    cur: int,
    best_t: IntArray,
    best: int,
    tabu_until: IntArray,
    it0: int,
    n_iter: int,
    tenure_lo: int,
    tenure_hi: int,
    moves: int,
    stall: int,
    restart_after: int,
    kick: int,
    pos: IntArray,
    cand_pred: IntArray,
    cand_succ: IntArray,
    cand_prob: float,
    ncand: int,
) -> tuple[int, int, int, int, int, int]:
    """Run ``n_iter`` iterations. Returns (cur, best, evaluations, stall, restarts, last_it)."""
    n = t.shape[0] - 1
    evals = 0
    restarts = 0
    for step in range(n_iter):
        it = it0 + step
        chosen_kind, ci, cj, cdelta = -1, 0, 0, -(1 << 60)
        fallback_kind, fi, fj, fdelta = -1, 0, 0, -(1 << 60)
        full = moves == 0
        total = n * (n - 1) // 2 + n * (n - 1) if full else moves
        for m in range(total):
            if full:
                # enumerate: first all swaps (i < j), then all inserts (i, j)
                if m < n * (n - 1) // 2:
                    kind = 0
                    k = m
                    i = 1
                    while k >= n - i:
                        k -= n - i
                        i += 1
                    j = i + 1 + k
                else:
                    kind = 1
                    k = m - n * (n - 1) // 2
                    i = 1 + k // (n - 1)
                    j = k % (n - 1)  # j runs over 0..n without i - 1 and i
                    if j >= i - 1:
                        j += 2
            elif cand_prob > 0.0 and np.random.random() < cand_prob:
                # candidate-list insert: put a random read next to one of its best partners
                x0 = np.random.randint(1, n + 1)
                r = np.random.randint(0, ncand)
                kind = 1
                i = pos[x0]
                if np.random.random() < 0.5:
                    a0 = cand_pred[x0, r]  # right after a read that overlaps into x0
                    if a0 < 0:
                        continue
                    j = pos[a0]
                else:
                    b0 = cand_succ[x0, r]  # right before a read that x0 overlaps into
                    if b0 < 0:
                        continue
                    j = pos[b0] - 1
                if j == i or j == i - 1:
                    continue
            else:
                kind = np.random.randint(0, 2)
                i = np.random.randint(1, n + 1)
                if kind == 0:
                    j = np.random.randint(1, n + 1)
                    if j == i:
                        continue
                    if j < i:
                        i, j = j, i
                else:
                    j = np.random.randint(0, n + 1)
                    if j == i or j == i - 1:
                        continue
            if kind == 0:
                delta = swap_delta(t, w, i, j)
                x, y = t[i], t[j]
                is_tabu = tabu_until[x] > it or tabu_until[y] > it
            else:
                delta = insert_delta(t, w, i, j)
                x, y = t[i], 0
                is_tabu = tabu_until[x] > it
            evals += 1
            if is_tabu and cur + delta <= best:
                if delta > fdelta:
                    fallback_kind, fi, fj, fdelta = kind, i, j, delta
                continue
            if delta > cdelta:
                chosen_kind, ci, cj, cdelta = kind, i, j, delta
        if chosen_kind == -1:  # every move is tabu: take the best tabu one anyway
            chosen_kind, ci, cj, cdelta = fallback_kind, fi, fj, fdelta
        if chosen_kind == -1:
            continue
        if chosen_kind == 0:
            x, y = t[ci], t[cj]
            apply_swap(t, ci, cj)
            pos[t[ci]] = ci
            pos[t[cj]] = cj
            tabu_until[y] = it + np.random.randint(tenure_lo, tenure_hi + 1)
        else:
            x = t[ci]
            apply_insert(t, ci, cj)
            for p in range(min(ci, cj), min(max(ci, cj) + 2, n + 1)):
                pos[t[p]] = p
        tabu_until[x] = it + np.random.randint(tenure_lo, tenure_hi + 1)
        cur += cdelta
        if cur > best:
            best = cur
            best_t[:] = t
            stall = 0
        else:
            stall += 1
        if stall >= restart_after:  # restart from the best tour, shaken by random swaps
            t[:] = best_t
            for _ in range(kick):
                a = np.random.randint(1, n + 1)
                b = np.random.randint(1, n + 1)
                if a != b:
                    apply_swap(t, a, b)
            for p in range(n + 1):
                pos[t[p]] = p
            cur = tour_fitness(t, w)
            evals += 1
            tabu_until[:] = 0
            stall = 0
            restarts += 1
    return cur, best, evals, stall, restarts, it0 + n_iter


def _candidate_lists(w: NDArray[np.int32], k: int) -> tuple[IntArray, IntArray]:
    """For every read its ``k`` best predecessors and successors by overlap (-1 = none)."""
    n1 = w.shape[0]
    pred = np.full((n1, k), -1, dtype=np.int64)
    succ = np.full((n1, k), -1, dtype=np.int64)
    for x in range(1, n1):
        for lists, column in ((pred, w[:, x]), (succ, w[x, :])):
            col = column.astype(np.int64)
            col[0] = col[x] = 0
            best = np.argsort(-col, kind="stable")[:k]
            lists[x, : best.size] = np.where(col[best] > 0, best, -1)
    return pred, succ


def tabu_search(
    problem: Problem,
    budget: Budget,
    seed: int = 0,
    *,
    tenure: int | None = None,
    moves_per_iteration: int | None = None,
    restart_after: int | None = None,
    kick: int | None = None,
    init: str = "greedy",
    chunk: int = 20,
    candidate_fraction: float = 0.0,
    candidates: int = 6,
) -> SolveResult:
    """Tabu search; see the module docstring for the parameters."""
    recorder = Recorder(problem)
    n = problem.n
    w = np.ascontiguousarray(problem.weights)
    if n < 3:
        from genomeasm.solvers.exact import brute_force

        return brute_force(problem)
    seed_numba(seed)
    rng = np.random.default_rng(seed)
    if init == "greedy":
        t = tour_of(greedy(problem).order)
    elif init == "random":
        t = tour_of(rng.permutation(n).astype(np.int64))
    else:
        raise ValueError("init must be 'greedy' or 'random'")
    tenure = tenure if tenure is not None else max(5, round(math.sqrt(n)))
    lo, hi = max(1, math.ceil(0.5 * tenure)), max(1, math.ceil(1.5 * tenure))
    moves = moves_per_iteration if moves_per_iteration is not None else (0 if n <= 60 else 20 * n)
    restart_after = restart_after if restart_after is not None else max(200, 10 * n)
    kick = kick if kick is not None else max(3, round(0.05 * n))
    cur = int(tour_fitness(t, w))
    best, best_t = cur, t.copy()
    tabu_until = np.zeros(n + 1, dtype=np.int64)
    pos = np.empty(n + 1, dtype=np.int64)
    pos[t] = np.arange(n + 1)
    k = max(1, candidates)
    cand_pred, cand_succ = _candidate_lists(w, k)
    evals, it, stall, restarts = 1, 0, 0, 0
    recorder.note(best, evals)
    while not budget.used_up(recorder.elapsed(), evals):
        cur, best, e, stall, r, it = _tabu_chunk(
            t,
            w,
            cur,
            best_t,
            best,
            tabu_until,
            it,
            chunk,
            lo,
            hi,
            moves,
            stall,
            restart_after,
            kick,
            pos,
            cand_pred,
            cand_succ,
            float(candidate_fraction),
            k,
        )
        evals += int(e)
        restarts += int(r)
        recorder.note(int(best), evals)
    info = {"iterations": int(it), "restarts": restarts, "tenure": tenure, "moves": moves}
    return recorder.result("tabu", best_t, evals, info)
