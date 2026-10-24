"""Experiment design: instances, solver parameters and the function that runs one solver.

An *instance* is a window of a real genome, reads simulated from it (substitution errors) and the
tolerant overlap graph of those reads. Every (instance, seed) pair is a *read set*; all solvers of
an experiment run on the same read set (paired design). Seeds are derived from the experiment seed
with ``numpy.random.SeedSequence`` so that nothing depends on the order in which tasks run.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from genomeasm.experiments.bio import ReferenceIndex, assembly_quality, order_quality
from genomeasm.instance import AssemblyInstance
from genomeasm.io import Genome, load_genome
from genomeasm.simulate import ReadLengthModel, SimulationConfig, simulate_reads
from genomeasm.solvers.alns import ARM_NAMES, alns
from genomeasm.solvers.common import Budget, Problem, SolveResult
from genomeasm.solvers.greedy import greedy, nearest_neighbour
from genomeasm.solvers.lkh import lkh
from genomeasm.solvers.tabu import tabu_search

STAGES = {"tuning": 1, "main": 2, "scaling": 3}
FIXED_ALGORITHMS = [f"alns_fixed_{name}" for name in ARM_NAMES]


def load_config(path: Path) -> dict[str, Any]:
    return yaml.safe_load(Path(path).read_text(encoding="utf-8"))


@dataclass(frozen=True)
class InstanceSpec:
    id: str
    genome: str
    start: int
    length: int
    coverage: float
    error: float
    read_length: int
    tolerance: float
    min_overlap: int

    @classmethod
    def from_config(cls, cfg: dict[str, Any], entry: dict[str, Any]) -> InstanceSpec:
        return cls(
            id=entry["id"],
            genome=entry["genome"],
            start=int(entry["start"]),
            length=int(entry["length"]),
            coverage=float(entry["coverage"]),
            error=float(entry["error"]),
            read_length=int(entry["read_length"]),
            tolerance=float(entry.get("tolerance", cfg["overlap"]["tolerance"])),
            min_overlap=int(entry.get("min_overlap", cfg["overlap"]["min_overlap"])),
        )


def seeds_for(cfg: dict[str, Any], stage: str, spec_index: int, seed_index: int) -> tuple[int, int]:
    """(read seed, solver seed) of a read set, independent of the order tasks are run in."""
    state = np.random.SeedSequence(
        [int(cfg["seed"]), STAGES[stage], spec_index, seed_index]
    ).generate_state(2)
    return int(state[0]) % 2**31, int(state[1]) % 2**31


@dataclass
class Built:
    spec: InstanceSpec
    instance: AssemblyInstance
    problem: Problem
    window: str
    reference: ReferenceIndex
    build_seconds: float


_GENOMES: dict[str, Genome] = {}


def genome(accession: str, raw_dir: Path = Path("data/raw")) -> Genome:
    if accession not in _GENOMES:
        _GENOMES[accession] = load_genome(accession, raw_dir).genome
    return _GENOMES[accession]


def build(spec: InstanceSpec, read_seed: int) -> Built:
    """Simulate the reads of one read set and build its overlap graph."""
    started = time.perf_counter()
    window = genome(spec.genome).sequence[spec.start : spec.start + spec.length]
    reads = simulate_reads(
        Genome.from_sequence(spec.genome, window),
        SimulationConfig(
            coverage=spec.coverage,
            read_length=ReadLengthModel.fixed(spec.read_length),
            substitution_rate=spec.error,
            seed=read_seed,
        ),
    )
    inst = AssemblyInstance.from_readset(
        reads, min_overlap=spec.min_overlap, tolerance=spec.tolerance
    )
    return Built(
        spec, inst, Problem.from_instance(inst), window, ReferenceIndex(window),
        time.perf_counter() - started,
    )  # fmt: skip


# ---- parameters -> solver arguments ------------------------------------------------------------


def tabu_kwargs(params: dict[str, Any], n: int) -> dict[str, Any]:
    kwargs = {
        "tenure": max(3, round(params["tenure_factor"] * n**0.5)),
        "moves_per_iteration": round(params["moves_factor"] * n),
        "restart_after": max(50, round(params["restart_factor"] * n)),
        "kick": max(2, round(params["kick_fraction"] * n)),
    }
    if "candidate_fraction" in params:
        kwargs["candidate_fraction"] = float(params["candidate_fraction"])
        kwargs["candidates"] = int(params["candidates"])
    return kwargs


def alns_kwargs(cfg: dict[str, Any], params: dict[str, Any]) -> dict[str, Any]:
    shared = params["alns"]
    return {
        "destroy_fraction": (shared["destroy_lo"], max(shared["destroy_hi"], shared["destroy_lo"])),
        "q_cap": int(shared["q_cap"]),
        "noise": float(shared["noise"]),
        "final_ratio": float(shared["final_ratio"]),
        "accept_probability0": float(shared["accept_probability0"]),
        "regret_k": int(shared["regret_k"]),
        "p_worst": float(shared["p_worst"]),
        "rewards": tuple(cfg["reward_presets"][shared["reward_preset"]]),
    }


def run_algorithm(
    name: str,
    built: Built,
    cfg: dict[str, Any],
    params: dict[str, Any],
    seconds: float,
    seed: int,
    lkh_search_seconds: float | None = None,
) -> SolveResult:
    """Run one solver with its parameters. ``params`` has keys tabu / alns / policy / lkh."""
    p, n = built.problem, built.problem.n
    budget = Budget(seconds=seconds)
    if name == "greedy":
        return greedy(p)
    if name == "nearest_neighbour":
        return nearest_neighbour(p)
    if name in ("tabu", "tabu_cl"):
        return tabu_search(p, budget, seed=seed, **tabu_kwargs(params["tabu"], n))
    if name.startswith("alns_fixed_"):
        arm = ARM_NAMES.index(name.removeprefix("alns_fixed_"))
        return alns(p, budget, seed=seed, policy="fixed", fixed_arm=arm, **alns_kwargs(cfg, params))
    if name.startswith("alns_"):
        policy = name.removeprefix("alns_")
        return alns(
            p, budget, seed=seed, policy=policy, policy_params=params["policy"].get(policy, {}),
            **alns_kwargs(cfg, params),
        )  # fmt: skip
    if name == "lkh3":
        extra = params.get("lkh", {})
        search = lkh_search_seconds
        if search is None:  # equal wall-clock: subtract the preprocessing overhead of this instance
            overhead = lkh(p, seconds=0.2, seed=seed, extra_parameters=extra).seconds - 0.2
            search = max(1.0, seconds - overhead)
        result = lkh(p, seconds=search, seed=seed, extra_parameters=extra)
        result.info["search_seconds"] = search
        return result
    raise ValueError(f"unknown algorithm {name!r}")


def thin_trace(trace: list[tuple[float, int, int]], keep: int = 120) -> list[list[float]]:
    """At most ``keep`` points of the best-so-far trace, always including the first and last."""
    if len(trace) <= keep:
        return [[round(t, 4), e, f] for t, e, f in trace]
    idx = np.unique(np.linspace(0, len(trace) - 1, keep).astype(int))
    return [[round(trace[i][0], 4), trace[i][1], trace[i][2]] for i in idx]


def bio_metrics(built: Built, order: Any) -> dict[str, Any]:
    """Biological quality of the assembly made from an order of the reads."""
    inst = built.instance
    assert inst.true_starts is not None
    assembly = inst.assemble(order)
    return {
        **assembly_quality(assembly, built.reference),
        **order_quality(
            order, inst.true_starts, built.spec.read_length, [c.reads for c in assembly.contigs]
        ),
    }


def evaluate(built: Built, result: SolveResult, with_bio: bool = True) -> dict[str, Any]:
    """The record of one solver run: fitness, cost, trace and biological quality."""
    record: dict[str, Any] = {
        "algorithm": result.solver,
        "fitness": int(result.fitness),
        "merged_length": int(result.merged_length),
        "evaluations": int(result.evaluations),
        "seconds": round(result.seconds, 4),
        "trace": thin_trace(result.trace),
    }
    if with_bio:
        record["bio"] = bio_metrics(built, result.order)
    keep = ("iterations", "restarts", "t0", "search_seconds", "outcomes", "arm_names")
    info = {k: v for k, v in result.info.items() if k in keep}
    if "probability_trace" in result.info:
        pt = result.info["probability_trace"]
        step = max(1, len(pt) // 100)
        info["probability_trace"] = pt[::step]
        info["trace_every"] = result.info["trace_every"] * step
    record["info"] = info
    return record
