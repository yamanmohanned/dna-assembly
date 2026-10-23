"""Correctness validation of the solvers (Phase 4 exit gate).

Plain-language story: before we compare solvers we make sure each one is *right*. The O(1) shortcut
for the size of a move must give the same answer as recomputing everything; the two exact solvers
must agree; every solver must return a genuine order of the reads; and on instances small enough to
be solved exactly, the heuristics must find the best order (or we report how far off they are).
Pass rules are in ``configs/phase4_correctness.yaml``; results go to
``results/phase4/correctness.json``.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np

from genomeasm import provenance
from genomeasm.checks import Checks, mean_and_se
from genomeasm.checks import add_check as _add
from genomeasm.instance import AssemblyInstance
from genomeasm.io import Genome
from genomeasm.objective import is_permutation
from genomeasm.simulate import ReadLengthModel, SimulationConfig, simulate_reads
from genomeasm.solvers.alns import ARM_NAMES, alns
from genomeasm.solvers.common import Budget, Problem, SolveResult, tour_of
from genomeasm.solvers.exact import brute_force, held_karp
from genomeasm.solvers.greedy import greedy, nearest_neighbour
from genomeasm.solvers.lkh import LkhError, lkh, lkh_binary
from genomeasm.solvers.moves import (
    apply_insert,
    apply_swap,
    insert_delta,
    swap_delta,
    tour_fitness,
)
from genomeasm.solvers.tabu import tabu_search
from genomeasm.solvers.warmup import warmup

Runner = Callable[[Problem, int], SolveResult]


def random_matrix_problem(rng: np.random.Generator, n: int) -> Problem:
    """A random overlap matrix (synthetic: used only to test algorithms, never for results)."""
    ov = rng.integers(1, 30, size=(n, n)).astype(np.int32)
    ov[rng.random((n, n)) < 0.5] = 0
    np.fill_diagonal(ov, 0)
    return Problem.from_matrix(ov, np.full(n, 40, dtype=np.int64))


def real_window_problem(genome: Genome, rng: np.random.Generator, cfg: dict[str, Any]) -> Problem:
    """A small instance: a window of the real genome, simulated reads, tolerant overlaps."""
    w = cfg["window"]
    while True:
        length = int(rng.integers(w["length_min"], w["length_max"] + 1))
        start = int(rng.integers(0, genome.length - length + 1))
        piece = Genome.from_sequence(genome.accession, genome.sequence[start : start + length])
        reads = simulate_reads(
            piece,
            SimulationConfig(
                coverage=w["coverage"],
                read_length=ReadLengthModel.fixed(w["read_length"]),
                substitution_rate=w["substitution_rate"],
                seed=int(rng.integers(2**31)),
            ),
        )
        inst = AssemblyInstance.from_readset(
            reads, min_overlap=w["min_overlap"], tolerance=w["tolerance"]
        )
        if 8 <= inst.n <= 15:
            return Problem.from_instance(inst)


def _delta_block(
    cfg: dict[str, Any], problems: list[Problem], rng: np.random.Generator, checks: Checks
) -> dict[str, Any]:
    moves = int(cfg["delta"]["moves_per_instance"])
    checked = mismatches = 0
    for problem in problems:
        n, w = problem.n, problem.weights
        t = tour_of(rng.permutation(n).astype(np.int64))
        for _ in range(moves):
            if rng.random() < 0.5:
                i, j = sorted(rng.choice(np.arange(1, n + 1), size=2, replace=False).tolist())
                before = tour_fitness(t, w)
                predicted = swap_delta(t, w, i, j)
                apply_swap(t, i, j)
            else:
                i = int(rng.integers(1, n + 1))
                j = int(rng.choice([k for k in range(n + 1) if k not in (i, i - 1)]))
                before = tour_fitness(t, w)
                predicted = insert_delta(t, w, i, j)
                apply_insert(t, i, j)
            actual = tour_fitness(t, w) - before
            checked += 1
            mismatches += int(actual != predicted)
            assert is_permutation(t[1:] - 1, n) and t[0] == 0
    _add(
        checks,
        "delta.equals_full_recomputation",
        "every swap and insert delta equals the difference of two full recomputations",
        mismatches == 0,
        moves_checked=checked,
        mismatches=mismatches,
    )
    return {"moves_checked": checked, "mismatches": mismatches}


def _exact_block(cfg: dict[str, Any], rng: np.random.Generator, checks: Checks) -> dict[str, Any]:
    block = cfg["exact_agreement"]
    disagreements = 0
    for _ in range(int(block["instances"])):
        n = int(rng.integers(block["n_min"], block["n_max"] + 1))
        problem = random_matrix_problem(rng, n)
        a, b = brute_force(problem), held_karp(problem)
        disagreements += int(a.fitness != b.fitness)
    _add(
        checks,
        "exact.brute_force_equals_held_karp",
        "brute force and Held-Karp find the same optimal fitness",
        disagreements == 0,
        instances=block["instances"],
        disagreements=disagreements,
    )
    return {"instances": block["instances"], "disagreements": disagreements}


def _solver_table(cfg: dict[str, Any]) -> dict[str, Runner]:
    budget = Budget(evaluations=int(cfg["small"]["evaluation_budget"]))
    fixed_budget = Budget(evaluations=int(cfg["small"]["fixed_arm_budget"]))

    def run_greedy(p: Problem, s: int) -> SolveResult:
        return greedy(p)

    def run_nn(p: Problem, s: int) -> SolveResult:
        return nearest_neighbour(p)

    def run_tabu(p: Problem, s: int) -> SolveResult:
        return tabu_search(p, budget, seed=s)

    def run_alns(policy: str, arm: int | None = None) -> Runner:
        def run(p: Problem, s: int) -> SolveResult:
            used = budget if arm is None else fixed_budget
            return alns(p, used, seed=s, policy=policy, fixed_arm=arm)

        return run

    table: dict[str, Runner] = {
        "greedy": run_greedy,
        "nearest_neighbour": run_nn,
        "tabu": run_tabu,
    }
    for policy in cfg["policies"]:
        table[f"alns_{policy}"] = run_alns(policy)
    for arm, name in enumerate(ARM_NAMES):
        table[f"alns_fixed_{name}"] = run_alns("fixed", arm)
    return table


def _small_block(
    cfg: dict[str, Any],
    genome: Genome,
    rng: np.random.Generator,
    checks: Checks,
    say: Callable[[str], None],
) -> dict[str, Any]:
    block = cfg["small"]
    problems = [
        (f"random_matrix_n{n}", random_matrix_problem(rng, n))
        for n in block["sizes"]
        for _ in range(int(block["per_size"]))
    ] + [
        ("real_window", real_window_problem(genome, rng, block))
        for _ in range(int(block["real_windows"]))
    ]
    optimum = [held_karp(p).fitness for _, p in problems]
    table = _solver_table(cfg)
    stats: dict[str, dict[str, Any]] = {}
    all_valid = True
    for name, runner in table.items():
        say(f"small instances: {name}")
        gaps, hits = [], 0
        for k, ((_, problem), best) in enumerate(zip(problems, optimum, strict=True)):
            result = runner(problem, k)
            all_valid &= is_permutation(result.order, problem.n)
            all_valid &= result.fitness == problem.fitness(result.order)
            gap = 100.0 * (best - result.fitness) / best if best else 0.0
            gaps.append(gap)
            hits += int(result.fitness == best)
        stats[name] = {
            "runs": len(problems),
            "optimal": hits,
            "success_rate": hits / len(problems),
            "mean_gap_percent": float(np.mean(gaps)),
            "max_gap_percent": float(np.max(gaps)),
        }
    if lkh_binary() is not None:
        say("small instances: lkh3")
        gaps, hits, valid = [], 0, True
        subset = list(zip(problems, optimum, strict=True))[
            :: max(1, len(problems) // int(block["lkh_instances"]))
        ]
        for k, ((_, problem), best) in enumerate(subset[: int(block["lkh_instances"])]):
            result = lkh(problem, seconds=float(block["lkh_seconds"]), seed=k)
            valid &= is_permutation(result.order, problem.n)
            gaps.append(100.0 * (best - result.fitness) / best if best else 0.0)
            hits += int(result.fitness == best)
        all_valid &= valid
        stats["lkh3"] = {
            "runs": len(gaps),
            "optimal": hits,
            "success_rate": hits / len(gaps),
            "mean_gap_percent": float(np.mean(gaps)),
            "max_gap_percent": float(np.max(gaps)),
        }
    _add(
        checks,
        "solvers.valid_permutations",
        "every returned solution is a valid permutation and its reported fitness is correct",
        all_valid,
        solvers=len(stats),
        instances=len(problems),
    )
    required = float(cfg["tolerances"]["min_success_rate"])
    for name, s in stats.items():
        if name in ("greedy", "nearest_neighbour") or name.startswith("alns_fixed_"):
            continue
        _add(
            checks,
            f"optimality.{name}",
            "finds the exact optimum on at least the required share of small instances",
            s["success_rate"] >= required,
            success_rate=s["success_rate"],
            required=required,
            runs=s["runs"],
            max_gap_percent=s["max_gap_percent"],
        )
    return {"instances": [name for name, _ in problems], "solvers": stats}


def _reproducibility_block(
    cfg: dict[str, Any], rng: np.random.Generator, checks: Checks
) -> dict[str, Any]:
    problem = random_matrix_problem(rng, 25)
    budget = Budget(evaluations=30000)
    runs: dict[str, tuple[SolveResult, SolveResult]] = {
        "tabu": (tabu_search(problem, budget, seed=7), tabu_search(problem, budget, seed=7)),
    }
    for policy in cfg["policies"]:
        runs[f"alns_{policy}"] = (
            alns(problem, budget, seed=7, policy=policy),
            alns(problem, budget, seed=7, policy=policy),
        )
    same = {
        k: bool(np.array_equal(a.order, b.order) and a.fitness == b.fitness)
        for k, (a, b) in runs.items()
    }
    _add(
        checks,
        "reproducibility.same_seed_same_result",
        "with an evaluation budget, the same seed gives the same order and fitness",
        all(same.values()),
        solvers=sorted(same),
    )
    return {"same_seed_identical": same}


def _preview_block(
    cfg: dict[str, Any], genome: Genome, checks: Checks, say: Callable[[str], None]
) -> dict[str, Any]:
    p = cfg["preview"]
    reads = simulate_reads(
        genome,
        SimulationConfig(
            coverage=p["coverage"],
            read_length=ReadLengthModel.fixed(p["read_length"]),
            substitution_rate=p["substitution_rate"],
            seed=p["read_seed"],
        ),
    )
    inst = AssemblyInstance.from_readset(
        reads, min_overlap=p["min_overlap"], tolerance=p["tolerance"]
    )
    problem = Problem.from_instance(inst)
    budget = Budget(seconds=float(p["seconds"]))
    seed = int(p["solver_seed"])
    results: dict[str, SolveResult] = {
        "greedy": greedy(problem),
        "nearest_neighbour": nearest_neighbour(problem),
    }
    say("preview: tabu")
    results["tabu"] = tabu_search(problem, budget, seed=seed)
    for policy in cfg["policies"]:
        say(f"preview: alns {policy}")
        results[f"alns_{policy}"] = alns(problem, budget, seed=seed, policy=policy)
    if lkh_binary() is not None:
        say("preview: lkh3")
        results["lkh3"] = lkh(problem, seconds=float(p["seconds"]), seed=seed)
    valid = all(
        is_permutation(r.order, inst.n) and r.fitness == inst.fitness(r.order)
        for r in results.values()
    )
    _add(
        checks,
        "preview.valid_and_verified",
        "on the realistic instance every solution is a permutation and its fitness is re-verified",
        valid,
        n=inst.n,
    )
    return {
        "instance": {**inst.summary(), "true_order_fitness": inst.fitness(inst.true_order())},
        "note": "one seed, one instance: a smoke comparison, not a result (statistics: Phase 5)",
        "solvers": {
            name: {
                "fitness": r.fitness,
                "merged_length": r.merged_length,
                "evaluations": r.evaluations,
                "seconds": r.seconds,
                "trace_points": len(r.trace),
                "info": {
                    k: v
                    for k, v in r.info.items()
                    if k in ("iterations", "restarts", "t0", "lkh_reported_cost")
                },
            }
            for name, r in results.items()
        },
    }


def run_solver_validation(
    cfg: dict[str, Any],
    genomes: dict[str, Genome],
    *,
    repo_dir: Path | None = None,
    progress: Callable[[str], None] | None = None,
) -> dict[str, Any]:
    """Run every Phase 4 check and return the result dictionary."""
    say = progress or (lambda message: None)
    warmup()
    genome = genomes[cfg["genome"]]
    streams = np.random.SeedSequence(int(cfg["seed"])).spawn(5)
    rngs = [np.random.default_rng(s) for s in streams]
    checks: Checks = []
    d = cfg["delta"]
    delta_problems = [
        random_matrix_problem(rngs[0], int(d["synthetic_n"]))
        for _ in range(int(d["synthetic_instances"]))
    ]
    delta_problems += [real_window_problem(genome, rngs[0], cfg["small"]) for _ in range(3)]
    say("delta evaluation")
    delta = _delta_block(cfg, delta_problems, rngs[0], checks)
    say("exact solvers")
    exact = _exact_block(cfg, rngs[1], checks)
    small = _small_block(cfg, genome, rngs[2], checks, say)
    reproducibility = _reproducibility_block(cfg, rngs[3], checks)
    preview = _preview_block(cfg, genome, checks, say)
    return {
        "name": cfg["name"],
        "config": cfg,
        "provenance": provenance.collect(cfg, seed=int(cfg["seed"]), repo_dir=repo_dir),
        "lkh_available": lkh_binary() is not None,
        "delta": delta,
        "exact": exact,
        "small": small,
        "reproducibility": reproducibility,
        "preview": preview,
        "checks": checks,
        "n_checks": len(checks),
        "n_failed": sum(not c["passed"] for c in checks),
        "all_passed": all(c["passed"] for c in checks),
    }


__all__ = ["LkhError", "mean_and_se", "run_solver_validation"]
