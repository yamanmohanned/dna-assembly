"""LKH-3 reference solver, run through Ubuntu (WSL) and verified independently.

LKH-3 (Keld Helsgaun's implementation of the Lin-Kernighan-Helsgaun heuristic) is built from the
source in ``tools/`` (version and provenance in docs/DECISIONS.md D-033). This module writes the
ATSP instance as a TSPLIB file (``TYPE: ATSP``, ``EDGE_WEIGHT_FORMAT: FULL_MATRIX``), writes a
parameter file, runs ``LKH`` and reads the tour back. The returned tour is **never trusted**: it
must be a Hamiltonian cycle over all nodes and its cost, recomputed here from the cost matrix, must
equal both the fitness-based length and the cost that LKH-3 reports.

Budget: ``TIME_LIMIT`` is the search time per run; LKH-3 also spends time before the search
(preprocessing, ascent). Both the total wall-clock time and LKH-3's own search time are reported;
experiments compare the total (rule R6). LKH-3 does not report a number of evaluations.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

import numpy as np

from genomeasm import objective
from genomeasm.solvers.common import IntArray, Problem, SolveResult, order_of

DEFAULT_LKH = Path("tools/LKH3-master/LKH")
WSL_DISTRO = "Ubuntu-24.04"


class LkhError(RuntimeError):
    """LKH-3 is missing, failed, or returned a tour that does not check out."""


def lkh_binary(path: Path | None = None) -> Path | None:
    """Path of the LKH-3 executable, or None if it has not been built."""
    candidate = Path(path) if path else Path(os.environ.get("GENOMEASM_LKH", str(DEFAULT_LKH)))
    return candidate.resolve() if candidate.is_file() else None


def to_wsl_path(path: Path) -> str:
    """Windows path to the matching /mnt/<drive>/... path inside WSL."""
    p = Path(path).resolve()
    drive = p.drive.rstrip(":").lower()
    rest = "/".join(p.parts[1:])
    return f"/mnt/{drive}/{rest}"


def parse_tour_file(text: str) -> tuple[IntArray, int | None]:
    """(1-based node ids in tour order, cost from the ``Length`` comment if present)."""
    cost_match = re.search(r"Length\s*=\s*(-?\d+)", text)
    section = text.split("TOUR_SECTION", 1)
    if len(section) != 2:
        raise LkhError("no TOUR_SECTION in the tour file")
    ids: list[int] = []
    for token in section[1].split():
        value = int(token)
        if value == -1:
            break
        ids.append(value)
    return np.array(ids, dtype=np.int64), int(cost_match.group(1)) if cost_match else None


def verify_tour(problem: Problem, node_ids: IntArray, reported_cost: int | None) -> IntArray:
    """Check LKH-3's tour and return the order of reads (raises ``LkhError`` otherwise)."""
    n = problem.n
    tour = node_ids - 1  # TSPLIB nodes are 1-based; node 1 is our dummy node 0
    if tour.size != n + 1 or not np.array_equal(np.sort(tour), np.arange(n + 1)):
        raise LkhError(f"LKH-3 tour is not a Hamiltonian cycle over {n + 1} nodes")
    cost = problem.cost_matrix()
    cycle = objective.cycle_cost(cost, tour)
    order = objective.order_from_tour(tour, n)
    fitness_length = problem.total_length - problem.fitness(order)
    if cycle != fitness_length:
        raise LkhError(f"cycle cost {cycle} differs from sum(len) - F = {fitness_length}")
    if reported_cost is not None and reported_cost != cycle:
        raise LkhError(f"LKH-3 reports cost {reported_cost}, the tour costs {cycle}")
    return order


def lkh(
    problem: Problem,
    seconds: float,
    seed: int = 0,
    *,
    binary: Path | None = None,
    keep_dir: Path | None = None,
    extra_parameters: dict[str, str] | None = None,
) -> SolveResult:
    """Run LKH-3 on the instance for about ``seconds`` of search time and verify the result."""
    exe = lkh_binary(binary)
    if exe is None:
        raise LkhError(
            "LKH-3 is not built (expected tools/LKH3-master/LKH); it is never replaced by "
            "another solver"
        )
    workdir = Path(keep_dir) if keep_dir else Path(tempfile.mkdtemp(prefix="genomeasm_lkh_"))
    workdir.mkdir(parents=True, exist_ok=True)
    try:
        objective.write_tsplib_atsp(
            workdir / "inst.atsp", problem.cost_matrix(), name="assembly", comment="genomeasm"
        )
        parameters = {
            "PROBLEM_FILE": "inst.atsp",
            "TOUR_FILE": "inst.tour",
            "RUNS": "1",
            "MAX_TRIALS": "1000000000",
            "TIME_LIMIT": f"{seconds:g}",
            "SEED": str(seed),
            "TRACE_LEVEL": "1",
            **(extra_parameters or {}),
        }
        (workdir / "inst.par").write_text(
            "".join(f"{k} = {v}\n" for k, v in parameters.items()), encoding="ascii", newline="\n"
        )
        command = f'cd "{to_wsl_path(workdir)}" && "{to_wsl_path(exe)}" inst.par'
        started = time.perf_counter()
        done = subprocess.run(
            ["wsl", "-d", WSL_DISTRO, "--", "bash", "-c", command],
            capture_output=True,
            text=True,
            timeout=max(60.0, 10 * seconds + 60),
            check=False,
        )
        wall = time.perf_counter() - started
        tour_file = workdir / "inst.tour"
        if done.returncode != 0 or not tour_file.is_file():
            raise LkhError(f"LKH-3 failed (exit {done.returncode}): {done.stdout[-400:]}")
        node_ids, reported = parse_tour_file(tour_file.read_text(encoding="ascii"))
        order = verify_tour(problem, node_ids, reported)
        fit = problem.fitness(order)
        trace = [
            (float(t), 0, problem.total_length - int(cost))
            for cost, t in re.findall(r"Cost\s*=\s*(-?\d+).*?Time\s*=\s*([\d.]+)", done.stdout)
        ]
        return SolveResult(
            solver="lkh3",
            order=order,
            fitness=fit,
            merged_length=problem.total_length - fit,
            evaluations=0,
            seconds=wall,
            trace=trace,
            info={"lkh_reported_cost": reported, "stdout_tail": done.stdout[-300:]},
        )
    finally:
        if keep_dir is None:
            shutil.rmtree(workdir, ignore_errors=True)


__all__ = ["LkhError", "lkh", "lkh_binary", "order_of", "parse_tour_file", "verify_tour"]
