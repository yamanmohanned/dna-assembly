"""Measure the worst interactive runs of the application (the 30-second rule, measured not assumed).

Each scenario walks through what a page does after a click: build the instance, run the solver
with the progress callback and live redraws, compute the quality, draw the final charts. The
worst parameter values the interface allows are used. Numba kernels are compiled before timing,
as the application does at start-up.
"""

from __future__ import annotations

import time
from typing import Any

from genomeasm.appkit import figs, lab
from genomeasm.solvers.warmup import warmup

LIMIT = lab.MAX_TOTAL_SECONDS


def _worst_instance(tolerance: float = 0.05, min_overlap: int = 20) -> lab.InstanceParams:
    return lab.InstanceParams(
        genome="NC_001416.1", start=10000, length=lab.MAX_WINDOW, coverage=20.0, read_length=100,
        error=0.03, seed=1, tolerance=tolerance, min_overlap=min_overlap,
    )  # fmt: skip


def _timed(steps: dict[str, Any]) -> dict[str, float]:
    out: dict[str, float] = {}
    for name, step in steps.items():
        started = time.perf_counter()
        step()
        out[name] = round(time.perf_counter() - started, 3)
    out["total"] = round(sum(out.values()), 3)
    return out


def scenario_solver_lab(solver: str, policy: str = "ucb1_rt") -> dict[str, Any]:
    """Solver Lab after a click: build (cold), run for the maximum budget with live redraws."""
    params = _worst_instance()
    lab.build.cache_clear()
    lab.simulate.cache_clear()
    state: dict[str, Any] = {}
    redraws: list[float] = []

    def on_progress(event: dict[str, Any]) -> None:
        started = time.perf_counter()
        history = [(event["seconds"], event.get("probabilities", [1 / 6] * 6))]
        figs.live_convergence([(event["seconds"], event["best"])], {"greedy": 0.0}, 20.0)
        figs.operator_probabilities(history * 2, 20.0)
        redraws.append(time.perf_counter() - started)

    def build() -> None:
        state["built"] = lab.build(params)

    def run() -> None:
        state["result"] = lab.run_solver(
            state["built"], solver, seconds=lab.MAX_SECONDS, policy=policy, progress=on_progress
        )

    def quality() -> None:
        state["quality"] = lab.quality(state["built"], state["result"].order)

    def render() -> None:
        result = state["result"]
        figs.figure_bytes(figs.live_convergence([(0.0, 1.0)], {}, 20.0), "png")
        if result.info.get("outcomes"):
            figs.outcome_bars(result.info["outcomes"])

    timings = _timed({"build": build, "solve": run, "quality": quality, "render": render})
    timings["reads"] = float(state["built"].problem.n)
    timings["live_redraws"] = float(len(redraws))
    timings["redraw_seconds_each"] = round(sum(redraws) / max(len(redraws), 1), 3)
    return timings


def scenario_assembly_viewer(solver: str = "alns") -> dict[str, Any]:
    """Assembly Viewer: build, solve for 15 s, assemble, quality, dot plot."""
    params = _worst_instance()
    state: dict[str, Any] = {}

    def build() -> None:
        state["built"] = lab.build(params)

    def run() -> None:
        state["result"] = lab.run_solver(state["built"], solver, seconds=15.0)

    def assemble() -> None:
        state["assembly"] = state["built"].instance.assemble(state["result"].order)
        state["quality"] = lab.quality(state["built"], state["result"].order)
        state["dot"] = lab.dot_plot(state["built"], state["assembly"])

    def render() -> None:
        figs.figure_bytes(figs.dot_plot(state["dot"]), "png")

    return _timed(
        {"build": build, "solve": run, "assemble_and_quality": assemble, "render": render}
    )


def scenario_overlap_graph() -> dict[str, Any]:
    """Overlap Graph at the loosest settings (tolerance 10%, tau 5) with all charts."""
    params = _worst_instance(tolerance=0.10, min_overlap=5)
    state: dict[str, Any] = {}

    def build() -> None:
        state["built"] = lab.build(params)

    def render() -> None:
        inst = state["built"].instance
        for fig in (
            figs.overlap_heatmap(inst, 300),
            figs.overlap_lengths(inst),
            figs.degree_distribution(inst),
        ):
            figs.figure_bytes(fig, "png")

    return _timed({"build": build, "render": render})


def scenario_simulator() -> dict[str, Any]:
    """Simulator: simulate and draw both charts for the largest allowed set."""
    params = _worst_instance()
    state: dict[str, Any] = {}

    def simulate() -> None:
        lab.simulate.cache_clear()
        state["reads"] = lab.simulate(params)

    def render() -> None:
        figs.figure_bytes(figs.depth_vs_theory(state["reads"])[0], "png")
        figs.figure_bytes(figs.read_map(state["reads"]), "png")

    return _timed({"simulate": simulate, "render": render})


def measure() -> dict[str, Any]:
    """All worst-case scenarios; ``ok`` is true if every total stays under the 30-second rule."""
    warmup()
    lab.build(lab.InstanceParams(length=600, coverage=4))  # imports and caches are warm
    scenarios: dict[str, Any] = {
        "solver_lab_alns_20s": scenario_solver_lab("alns"),
        "solver_lab_tabu_cl_20s": scenario_solver_lab("tabu_cl"),
        "assembly_viewer_alns_15s": scenario_assembly_viewer(),
        "overlap_graph_loosest": scenario_overlap_graph(),
        "simulator_largest": scenario_simulator(),
    }
    worst = max(s["total"] for s in scenarios.values())
    return {
        "limit_seconds": LIMIT,
        "worst_total_seconds": worst,
        "ok": worst < LIMIT,
        "scenarios": scenarios,
        "interface_limits": {
            "max_reads": lab.MAX_READS,
            "max_window": lab.MAX_WINDOW,
            "max_solver_seconds": lab.MAX_SECONDS,
        },
    }
