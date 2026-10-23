"""Compile (or load from the numba cache) every solver kernel before anything is timed.

The first call of a numba function compiles it, which can take many seconds. Experiments call
``warmup()`` once per process so that wall-clock budgets measure the search, not the compiler.
"""

from __future__ import annotations

import numpy as np

from genomeasm.solvers.alns import alns
from genomeasm.solvers.common import Budget, Problem
from genomeasm.solvers.exact import held_karp
from genomeasm.solvers.greedy import greedy, nearest_neighbour
from genomeasm.solvers.tabu import tabu_search


def warmup() -> None:
    rng = np.random.default_rng(0)
    ov = rng.integers(0, 6, size=(12, 12)).astype(np.int32)
    np.fill_diagonal(ov, 0)
    problem = Problem.from_matrix(ov, np.full(12, 20, dtype=np.int64))
    greedy(problem)
    nearest_neighbour(problem)
    held_karp(problem)
    tabu_search(problem, Budget(evaluations=3000), seed=0)
    tabu_search(problem, Budget(evaluations=3000), seed=0, moves_per_iteration=20)
    alns(problem, Budget(evaluations=20000), seed=0, policy="roundrobin", calibration_moves=6)
