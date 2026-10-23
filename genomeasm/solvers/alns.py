"""Adaptive Large Neighbourhood Search (ALNS) with bandit-based operator selection.

Plain-language story: each iteration breaks the current solution (*destroy*: take some reads out)
and rebuilds it (*repair*: put them back, hopefully better). Accept the new solution if better,
or sometimes if it is a little worse (*simulated annealing*: the allowed worsening shrinks as the
search ages). Which destroy and repair tool to use is chosen by a *policy* (``genomeasm.bandit``).

Operators (an *arm* is a pair ``(destroy, repair)``; arm index = ``2 * destroy + repair``,
names in ``ARM_NAMES``):

* destroy 0 ``random``: remove ``q`` random reads.
* destroy 1 ``worst``: remove poorly connected reads: the contribution of read ``x`` between
  ``p`` and ``s`` is ``W[p,x] + W[x,s] - W[p,s]``; reads are drawn from the list sorted by
  ascending contribution at index ``floor(u ** p_worst * m)`` (Ropke and Pisinger's randomisation).
* destroy 2 ``segment``: remove ``q`` consecutive reads.
* repair 0 ``greedy``: repeatedly insert the removed read with the best insertion gain at its best
  place (regret-1).
* repair 1 ``regret``: repeatedly insert the removed read with the largest *regret* (sum over
  ``m = 2..k`` of best gain minus ``m``-th best gain), at its best place (ties: larger best gain).
  ``k = regret_k`` (default 2).
  Both use cached best insertion places that are updated only where the new edges can matter.
  Insertion gains used for *choosing* places carry a random jitter of +-``noise * max(W)`` (Ropke
  and Pisinger's noise idea): without it a small destroy-repair move mostly rebuilds the same tour
  (about 90% of the moves on small instances in the first run, decision D-034). The fitness is
  always recomputed exactly.

Acceptance: improving solutions always; a worse one with probability ``exp(delta / T)``. The
temperature falls geometrically from ``T0`` to ``T0 * final_ratio`` as the budget is used. **T0 is
calibrated** from the search itself: up to ``calibration_moves`` destroy-repair moves (stopping
when ``calibration_fraction`` = 10% of the budget is used) are tried from the start solution and
``T0 = median(worsening) / ln(1 / accept_probability0)``, so that a median-sized worsening is
accepted with probability ``accept_probability0`` (default 0.5) at the start. The calibration moves
count in the budget.

Rewards (normalised, Ropke and Pisinger's 33 / 9 / 13 / 0 divided by 33, parameter ``rewards``): new
global best 1.0, better than current 0.273, accepted although worse 0.394, anything else 0. Moves
that keep the fitness (a *lateral* move to another tour, or the same tour again) and rejected
moves earn 0: with many zero-weight edges lateral moves are frequent and would otherwise teach a
policy to chase moves that make no progress (found in the first run, decision D-034).

Evaluations: one per insertion-gain computed plus one full recomputation per iteration.
"""

from __future__ import annotations

import math
from typing import Any

import numba
import numpy as np
from numpy.typing import NDArray

from genomeasm.bandit import make_policy
from genomeasm.solvers.common import Budget, IntArray, Problem, Recorder, SolveResult, tour_of
from genomeasm.solvers.greedy import greedy
from genomeasm.solvers.moves import tour_fitness
from genomeasm.solvers.tabu import seed_numba

DESTROY = ("random", "worst", "segment")
REPAIR = ("greedy", "regret")
ARM_NAMES = tuple(f"{d}+{r}" for d in DESTROY for r in REPAIR)  # arm = 2 * destroy + repair
OUTCOMES = ("best", "improved", "accepted", "lateral", "unchanged", "rejected")
# rewards for: new global best, better than current, accepted although worse, anything else
DEFAULT_REWARDS = (1.0, 9.0 / 33.0, 13.0 / 33.0, 0.0)
_NEG = -(1 << 40)


@numba.njit(cache=True)
def _destroy_random(t: IntArray, q: int) -> IntArray:
    n = t.shape[0] - 1
    pos = np.random.permutation(n)[:q] + 1
    return t[pos].copy()


@numba.njit(cache=True)
def _destroy_segment(t: IntArray, q: int) -> IntArray:
    n = t.shape[0] - 1
    s = np.random.randint(1, n - q + 2)
    return t[s : s + q].copy()


@numba.njit(cache=True)
def _destroy_worst(t: IntArray, w: NDArray[np.int32], q: int, p_worst: float) -> IntArray:
    n1 = t.shape[0]
    n = n1 - 1
    contribution = np.empty(n, dtype=np.int64)
    for p in range(1, n + 1):
        a, x, b = t[p - 1], t[p], t[(p + 1) % n1]
        contribution[p - 1] = w[a, x] + w[x, b] - w[a, b]
    avail = np.argsort(contribution) + 1  # positions, smallest contribution first
    removed = np.empty(q, dtype=np.int64)
    m = n
    for r in range(q):
        idx = int(np.random.random() ** p_worst * m)
        removed[r] = t[avail[idx]]
        for s in range(idx, m - 1):
            avail[s] = avail[s + 1]
        m -= 1
    return removed


@numba.njit(cache=True)
def _jitter(noise: int) -> int:
    if noise == 0:
        return 0
    return np.random.randint(-noise, noise + 1)


@numba.njit(cache=True)
def _offer(g_row: IntArray, p_row: IntArray, gain: int, pred: int) -> None:
    k = g_row.shape[0]
    if gain <= g_row[k - 1]:
        return
    pos = k - 1
    while pos > 0 and g_row[pos - 1] < gain:
        g_row[pos] = g_row[pos - 1]
        p_row[pos] = p_row[pos - 1]
        pos -= 1
    g_row[pos] = gain
    p_row[pos] = pred


@numba.njit(cache=True)
def _scan(
    succ: IntArray, w: NDArray[np.int32], y: int, g_row: IntArray, p_row: IntArray, noise: int
) -> int:
    g_row[:] = _NEG
    p_row[:] = -1
    a = 0
    count = 0
    while True:
        b = succ[a]
        _offer(g_row, p_row, w[a, y] + w[y, b] - w[a, b] + _jitter(noise), a)
        count += 1
        a = b
        if a == 0:
            break
    return count


@numba.njit(cache=True)
def _repair(succ: IntArray, w: NDArray[np.int32], removed: IntArray, k: int, noise: int) -> int:
    """Insert all removed reads (regret-k; k = 1 is greedy insertion). Returns evaluations."""
    q = removed.shape[0]
    gains = np.full((q, k), _NEG, dtype=np.int64)
    preds = np.full((q, k), -1, dtype=np.int64)
    placed = np.zeros(q, dtype=np.bool_)
    evals = 0
    for r in range(q):
        evals += _scan(succ, w, removed[r], gains[r], preds[r], noise)
    for _ in range(q):
        sel, sel_regret, sel_gain = -1, _NEG, _NEG
        for r in range(q):
            if placed[r]:
                continue
            g0 = gains[r, 0]
            regret = 0
            for m in range(1, k):
                if preds[r, m] >= 0:
                    regret += g0 - gains[r, m]
            if sel == -1 or regret > sel_regret or (regret == sel_regret and g0 > sel_gain):
                sel, sel_regret, sel_gain = r, regret, g0
        x = removed[sel]
        a = preds[sel, 0]
        b = succ[a]
        succ[x] = b
        succ[a] = x
        placed[sel] = True
        for r in range(q):
            if placed[r]:
                continue
            y = removed[r]
            broken = False
            for m in range(k):
                if preds[r, m] == a:
                    broken = True
            if broken:
                evals += _scan(succ, w, y, gains[r], preds[r], noise)
            else:
                _offer(gains[r], preds[r], w[a, y] + w[y, x] - w[a, x] + _jitter(noise), a)
                _offer(gains[r], preds[r], w[x, y] + w[y, b] - w[x, b] + _jitter(noise), x)
                evals += 2
    return evals


@numba.njit(cache=True)
def alns_step(
    t: IntArray,
    w: NDArray[np.int32],
    destroy: int,
    repair_k: int,
    q: int,
    p_worst: float,
    noise: int,
    out: IntArray,
) -> tuple[int, int]:
    """One destroy-repair move from tour ``t``; the new tour is written to ``out``."""
    n1 = t.shape[0]
    if destroy == 0:
        removed = _destroy_random(t, q)
    elif destroy == 1:
        removed = _destroy_worst(t, w, q, p_worst)
    else:
        removed = _destroy_segment(t, q)
    mark = np.zeros(n1, dtype=np.bool_)
    for v in removed:
        mark[v] = True
    succ = np.zeros(n1, dtype=np.int64)
    prev = t[0]
    for p in range(1, n1):
        v = t[p]
        if not mark[v]:
            succ[prev] = v
            prev = v
    succ[prev] = t[0]
    evals = _repair(succ, w, removed, repair_k, noise)
    node = 0
    for p in range(n1):
        out[p] = node
        node = succ[node]
    return tour_fitness(out, w), evals + 1


def alns(
    problem: Problem,
    budget: Budget,
    seed: int = 0,
    *,
    policy: str = "ucb1",
    policy_params: dict[str, Any] | None = None,
    fixed_arm: int | None = None,
    destroy_fraction: tuple[float, float] = (0.10, 0.40),
    q_cap: int = 30,
    p_worst: float = 3.0,
    regret_k: int = 2,
    noise: float = 0.3,
    rewards: tuple[float, float, float, float] = DEFAULT_REWARDS,
    accept_probability0: float = 0.5,
    final_ratio: float = 1e-3,
    calibration_moves: int = 100,
    calibration_fraction: float = 0.1,
    init: str = "greedy",
    trace_every: int = 100,
) -> SolveResult:
    """ALNS with the given operator-selection policy (``fixed_arm``: one operator pair alone)."""
    recorder = Recorder(problem)
    n = problem.n
    if n < 4:
        from genomeasm.solvers.exact import brute_force

        return brute_force(problem)
    w = np.ascontiguousarray(problem.weights)
    seed_numba(seed)
    rng = np.random.default_rng(seed)
    cur = tour_of(greedy(problem).order if init == "greedy" else rng.permutation(n))
    cur_f = int(tour_fitness(cur, w))
    best_t, best_f = cur.copy(), cur_f
    cand = np.empty_like(cur)
    # q reads are removed per move, drawn from [q_lo, q_hi]: fractions of n, at most q_cap
    # (and q_hi at least 3, so that tiny instances still get a real neighbourhood)
    q_lo = max(1, min(q_cap, round(destroy_fraction[0] * n)))
    q_hi = min(n, max(q_lo, min(q_cap, max(3, round(destroy_fraction[1] * n)))))
    noise_scale = round(noise * int(np.max(w)))  # insertion gains get +-noise_scale of jitter
    n_arms = len(ARM_NAMES)
    if policy == "fixed":
        policy_params = {**(policy_params or {}), "arm": 0 if fixed_arm is None else fixed_arm}
    pol = make_policy(policy, n_arms, rng, **(policy_params or {}))
    r_best, r_improved, r_accepted, r_other = rewards
    evals, iterations = 1, 0
    recorder.note(best_f, evals)

    def move(arm: int) -> tuple[int, int]:
        destroy, repair = divmod(arm, 2)
        q = int(rng.integers(q_lo, q_hi + 1))
        return alns_step(
            cur, w, destroy, 1 if repair == 0 else regret_k, q, p_worst, noise_scale, cand
        )

    # calibrate the start temperature from the moves a random operator makes at the start solution
    worsenings: list[int] = []
    for _ in range(calibration_moves):
        new_f, e = move(int(rng.integers(n_arms)))
        evals += int(e)
        if new_f < cur_f:
            worsenings.append(cur_f - int(new_f))
        elif new_f > best_f:
            best_f, best_t[:] = int(new_f), cand
        if budget.fraction(recorder.elapsed(), evals) >= calibration_fraction:
            break
    t0 = float(np.median(worsenings)) / math.log(1.0 / accept_probability0) if worsenings else 1.0
    t0 = max(t0, 1e-9)

    outcomes = np.zeros((n_arms, len(OUTCOMES)), dtype=np.int64)
    arm_history: list[int] = []
    probability_trace: list[list[float]] = []
    while not budget.used_up(recorder.elapsed(), evals):
        arm = pol.select()
        new_f, e = move(arm)
        new_f, evals, iterations = int(new_f), evals + int(e), iterations + 1
        delta = new_f - cur_f
        temperature = t0 * final_ratio ** budget.fraction(recorder.elapsed(), evals)
        if new_f > best_f:
            outcome, reward, accept = 0, r_best, True
        elif delta > 0:
            outcome, reward, accept = 1, r_improved, True
        elif delta == 0:  # same fitness: a lateral move to another tour, or the same tour again
            same = np.array_equal(cand, cur)
            outcome, reward, accept = (4 if same else 3), r_other, True
        elif rng.random() < math.exp(delta / temperature):
            outcome, reward, accept = 2, r_accepted, True
        else:
            outcome, reward, accept = 5, r_other, False
        outcomes[arm, outcome] += 1
        pol.update(arm, reward)
        arm_history.append(arm)
        if iterations % trace_every == 0:
            probability_trace.append([float(x) for x in pol.probabilities()])
        if accept:
            cur, cand = cand, cur
            cur_f = new_f
            if new_f > best_f:
                best_f, best_t[:] = new_f, cur
                recorder.note(best_f, evals)
    info = {
        "policy": policy,
        "fixed_arm": fixed_arm,
        "iterations": iterations,
        "t0": t0,
        "q_range": (q_lo, q_hi),
        "arm_names": ARM_NAMES,
        "outcome_names": OUTCOMES,
        "outcomes": outcomes.tolist(),
        "arm_history": np.array(arm_history, dtype=np.uint8),
        "probability_trace": probability_trace,
        "trace_every": trace_every,
    }
    return recorder.result(
        f"alns_{policy}" + ("" if fixed_arm is None else f"_{fixed_arm}"), best_t, evals, info
    )
