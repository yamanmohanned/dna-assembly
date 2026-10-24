"""Instances and solver runs for the interactive pages (simulator, overlap graph, solver lab).

Every interactive run must finish in under 30 seconds on the presenting laptop (master prompt,
Phase 6). The limits below keep the size of an instance and the solver budget inside that
envelope; ``tests/test_appkit.py`` measures the worst case instead of assuming it. Larger runs are
shown from the precomputed results of Phase 5 and are labelled as such.
"""

from __future__ import annotations

import json
import subprocess
import time
from collections.abc import Callable
from dataclasses import asdict, dataclass
from functools import lru_cache
from typing import Any

import numpy as np

from genomeasm.appkit.paths import CONFIG, GENOMES, RAW, RESULTS
from genomeasm.assembly import Assembly
from genomeasm.experiments import design, plan
from genomeasm.experiments.bio import WINDOW, ReferenceIndex
from genomeasm.instance import AssemblyInstance
from genomeasm.io import Genome
from genomeasm.simulate import ReadLengthModel, ReadSet, SimulationConfig, simulate_reads
from genomeasm.solvers.alns import alns
from genomeasm.solvers.common import Budget, IntArray, Problem, SolveResult
from genomeasm.solvers.greedy import greedy, nearest_neighbour
from genomeasm.solvers.lkh import WSL_DISTRO, lkh_binary
from genomeasm.solvers.tabu import tabu_search

# ---- limits that keep an interactive run under 30 seconds ----------------------------------------

MAX_READS = 1000  # the overlap graph of 1,000 reads takes about 3 s; 2,000 reads about 10 s
MAX_WINDOW = 5000
MAX_SECONDS = 20.0
MAX_TOTAL_SECONDS = 30.0

SOLVERS = {
    "greedy": "Greedy edge merging",
    "nearest_neighbour": "Nearest neighbour",
    "tabu": "Tabu search (swap + insert)",
    "tabu_cl": "Tabu search + candidate lists (post hoc)",
    "alns": "ALNS with an operator-selection policy",
    "lkh3": "LKH-3 (reference solver)",
}
POLICIES = {
    "random": "Random",
    "roundrobin": "Round robin",
    "roulette": "Roulette wheel (classical adaptive)",
    "epsilon_greedy": "Epsilon-greedy bandit",
    "ucb1": "UCB1 bandit (as pre-registered)",
    "ucb1_rt": "UCB1 bandit, random ties (post-hoc fix)",
}


@dataclass(frozen=True)
class InstanceParams:
    """Everything that defines a simulated read set and its overlap graph."""

    genome: str = "NC_001422.1"
    start: int = 0
    length: int = 3000
    coverage: float = 10.0
    read_length: int = 100
    error: float = 0.02
    seed: int = 1
    tolerance: float = 0.05
    min_overlap: int = 20

    def estimated_reads(self) -> int:
        return round(self.coverage * self.length / self.read_length)

    def problems(self) -> list[str]:
        """Reasons why this instance is not allowed (empty list: fine)."""
        issues = []
        if self.genome not in GENOMES:
            issues.append(f"unknown genome {self.genome}")
        if not 2 <= self.read_length <= self.length:
            issues.append("the read length must be between 2 and the window length")
        if self.estimated_reads() > MAX_READS:
            issues.append(
                f"about {self.estimated_reads()} reads: more than the {MAX_READS} allowed for a "
                "live run (the overlap graph would take too long); lower the coverage or window"
            )
        if self.length > MAX_WINDOW:
            issues.append(f"the window is longer than {MAX_WINDOW} letters")
        if not 0 <= self.error <= 0.2:
            issues.append("the error rate must lie between 0 and 0.2")
        if not 1 <= self.min_overlap < self.read_length:
            issues.append("the minimum overlap must lie between 1 and the read length - 1")
        return issues

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@lru_cache(maxsize=4)
def genome_sequence(accession: str) -> str:
    return design.genome(accession, RAW).sequence


def window_of(params: InstanceParams) -> str:
    return genome_sequence(params.genome)[params.start : params.start + params.length]


@lru_cache(maxsize=16)
def simulate(params: InstanceParams) -> ReadSet:
    """Simulate the reads of a window of a real genome (substitution errors, one strand)."""
    window = window_of(params)
    return simulate_reads(
        Genome.from_sequence(params.genome, window),
        SimulationConfig(
            coverage=params.coverage,
            read_length=ReadLengthModel.fixed(params.read_length),
            substitution_rate=params.error,
            seed=params.seed,
        ),
    )


@lru_cache(maxsize=8)
def build(params: InstanceParams) -> design.Built:
    """Reads, overlap graph and reference index of an instance (cached per parameter set)."""
    issues = params.problems()
    if issues:
        raise ValueError("; ".join(issues))
    started = time.perf_counter()
    reads = simulate(params)
    window = window_of(params)
    inst = AssemblyInstance.from_readset(
        reads, min_overlap=params.min_overlap, tolerance=params.tolerance
    )
    spec = design.InstanceSpec(
        id="live",
        genome=params.genome,
        start=params.start,
        length=params.length,
        coverage=params.coverage,
        error=params.error,
        read_length=params.read_length,
        tolerance=params.tolerance,
        min_overlap=params.min_overlap,
    )
    return design.Built(
        spec, inst, Problem.from_instance(inst), window, ReferenceIndex(window),
        time.perf_counter() - started,
    )  # fmt: skip


def tiny_instance(seed: int) -> design.Built:
    """A tiny exact-overlap instance with 5 to 8 reads (for the step-through mode).

    Windows of phiX174 are tried from ``seed`` upwards until the substring-free read set has 5 to 8
    reads, so the same seed always gives the same instance.
    """
    for k in range(400):
        params = InstanceParams(
            genome="NC_001422.1", start=(37 * (seed + k)) % 4800, length=240, coverage=2.6,
            read_length=80, error=0.0, seed=seed + k, tolerance=0.0, min_overlap=10,
        )  # fmt: skip
        built = build(params)
        if 5 <= built.instance.n <= 8:
            return built
    raise RuntimeError("no tiny instance found")  # pragma: no cover


# ---- solvers -------------------------------------------------------------------------------------


@lru_cache(maxsize=1)
def config() -> dict[str, Any]:
    return design.load_config(CONFIG)


@lru_cache(maxsize=1)
def tuning_summary() -> dict[str, Any] | None:
    path = RESULTS / "tuning_summary.json"
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else None


def solver_params(solver: str, policy: str, tuned: bool) -> dict[str, Any]:
    """Parameters of a solver: the tuned ones from Phase 5 or the defaults of the configuration."""
    name = {"alns": f"alns_{policy}"}.get(solver, solver)
    summary = tuning_summary()
    if tuned and summary is not None and name in summary["algorithms"]:
        return dict(summary["algorithms"][name]["params"])
    return plan.default_params(config(), name)


def lkh_available() -> bool:
    return lkh_binary() is not None


_LKH_STATE: dict[str, Any] = {"time": 0.0, "value": None}
LKH_FIX = (
    "LKH-3 is a Linux program that runs inside Ubuntu on WSL. How to fix: (1) open PowerShell and "
    f"run `wsl -d {WSL_DISTRO} -- echo ok`; (2) if it prints an error such as E_UNEXPECTED, "
    "restart "
    "the computer (or run `Restart-Service LxssManager` as administrator); (3) press the button "
    "to check again."
)
LKH_READY = (
    "LKH-3 is ready: Ubuntu (WSL) answers. Its total time includes file writing and preprocessing."
)


def lkh_status(force: bool = False, ttl: float = 60.0) -> tuple[bool, str]:
    """Can LKH-3 run right now? The program must exist and WSL must answer (cached for a minute)."""
    now = time.time()
    cached = _LKH_STATE["value"]
    if not force and cached is not None and now - _LKH_STATE["time"] < ttl:
        return (bool(cached[0]), str(cached[1]))
    if lkh_binary() is None:
        value = (
            False,
            "The LKH-3 program was not found (tools/LKH3-master/LKH): see tools/README.md.",
        )
    else:
        try:
            command = ["wsl", "-d", WSL_DISTRO, "--", "echo", "ok"]
            done = subprocess.run(command, capture_output=True, timeout=25, check=False)
            answer = done.stdout.decode("utf-8", "ignore").replace(chr(0), "").strip()
            if done.returncode == 0 and answer == "ok":
                value = (True, LKH_READY)
            else:
                value = (False, f"WSL does not work right now. {LKH_FIX}")
        except (OSError, subprocess.SubprocessError):
            value = (False, f"WSL did not answer within 25 seconds. {LKH_FIX}")
    _LKH_STATE.update(time=now, value=value)
    return value


def evaluations_text(solver: str, count: int) -> str:
    """LKH-3 does not count evaluations: say so instead of showing a misleading 0."""
    return "not reported" if solver == "lkh3" else f"{count:,}"


def run_solver(
    built: design.Built,
    solver: str,
    *,
    seconds: float,
    seed: int = 1,
    policy: str = "ucb1_rt",
    tuned: bool = True,
    progress: Callable[[dict[str, Any]], None] | None = None,
    step_log: list[dict[str, Any]] | None = None,
    evaluations: int | None = None,
) -> SolveResult:
    """Run one solver on a built instance. ``evaluations`` replaces the time budget (tests)."""
    if solver not in SOLVERS:
        raise ValueError(f"unknown solver {solver!r}")
    if solver == "alns" and policy not in POLICIES:
        raise ValueError(f"unknown policy {policy!r}")
    seconds = min(float(seconds), MAX_SECONDS)
    problem, n = built.problem, built.problem.n
    if solver == "greedy":
        return greedy(problem)
    if solver == "nearest_neighbour":
        return nearest_neighbour(problem)
    budget = Budget(evaluations=evaluations) if evaluations is not None else Budget(seconds=seconds)
    params = solver_params(solver, policy, tuned)
    if solver in ("tabu", "tabu_cl"):
        return tabu_search(
            problem, budget, seed=seed, progress=progress, **design.tabu_kwargs(params["tabu"], n)
        )
    if solver == "alns":
        return alns(
            problem,
            budget,
            seed=seed,
            policy=policy,
            policy_params=params["policy"].get(policy, {}),
            progress=progress,
            step_log=step_log,
            **design.alns_kwargs(config(), params),
        )
    # LKH-3: equal wall-clock needs a calibration run; it can fail (e.g. WSL down): the caller
    # reports it as unavailable, it is never replaced by another solver (rule R9)
    return design.run_algorithm(
        "lkh3", built, config(), params, min(seconds, 15.0), seed, lkh_search_seconds=None
    )


# ---- quality of a result -------------------------------------------------------------------------


def quality(built: design.Built, order: IntArray) -> dict[str, Any]:
    """Biological quality of the assembly made from ``order`` (identity, contigs, N50, ...)."""
    return design.bio_metrics(built, order)


def dot_plot(built: design.Built, assembly: Assembly) -> dict[str, Any]:
    """Where every 100-letter window of the assembly lies on the reference.

    Returns arrays x (position in the assembled text, contigs written one after the other),
    y (position on the reference, NaN if the window could not be located), contig (index of the
    contig of every window) and ``boundaries`` (start of every contig in x).
    """
    xs: list[int] = []
    ys: list[float] = []
    contig_of: list[int] = []
    boundaries: list[int] = []
    offset = 0
    for k, contig in enumerate(assembly.contigs):
        boundaries.append(offset)
        sequence = contig.sequence
        for start in range(0, len(sequence), WINDOW):
            window = sequence[start : start + WINDOW]
            if len(window) < 30:
                continue
            hit = built.reference.locate(window)
            xs.append(offset + start)
            ys.append(float("nan") if hit is None else float(hit[0]))
            contig_of.append(k)
        offset += len(sequence)
    return {
        "x": np.array(xs, dtype=float),
        "y": np.array(ys, dtype=float),
        "contig": np.array(contig_of, dtype=int),
        "boundaries": boundaries,
        "assembled_length": offset,
        "reference_length": len(built.window),
    }
