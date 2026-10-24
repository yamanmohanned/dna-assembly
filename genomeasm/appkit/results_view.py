"""Loading the precomputed results of Phase 5 for the Experiments page (read-only)."""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd

from genomeasm.appkit.paths import FIGURES, RESULTS
from genomeasm.experiments import analysis, runner

FIGURE_STEMS = {
    "fig04_convergence": "Convergence of every solver (mean gap against time and evaluations)",
    "fig05_final_fitness": "Final fitness of every solver after the same budget",
    "fig06_gap_to_lkh3": "Gap to LKH-3 against coverage and error rate",
    "fig07_ablation": "Ablation: operator-selection policies against random selection",
    "fig08_arm_probability": "How the policies spread their choices over the operator pairs",
    "fig09_operator_heatmap": "Success of each destroy + repair operator pair",
    "fig10_runtime_vs_n": "Runtime against instance size",
    "fig11_pareto": "Quality against time (Pareto front)",
}
EARLIER_STEMS = {
    "fig01a_coverage_histogram": "Coverage histogram against theory (Phase 2)",
    "fig02_overlap_length_distribution": "Overlap-length distribution (Phase 3)",
    "fig03_overlap_heatmap": "Overlap-matrix heatmap (Phase 3)",
}


@dataclass
class Phase5:
    summary: dict[str, Any]
    info: dict[str, Any]
    by_algorithm: pd.DataFrame
    by_instance: pd.DataFrame
    tests: pd.DataFrame
    scaling: pd.DataFrame
    true_reference: pd.DataFrame
    stages: dict[str, Any]


def available() -> bool:
    return (RESULTS / "analysis" / "tests.csv").exists()


@lru_cache(maxsize=1)
def load_phase5() -> Phase5:
    """The tables and summaries of Phase 5 (raises FileNotFoundError if they are missing)."""
    a = RESULTS / "analysis"

    def read_json(path: Path) -> dict[str, Any]:
        loaded = json.loads(path.read_text(encoding="utf-8"))
        assert isinstance(loaded, dict)
        return loaded

    return Phase5(
        summary=read_json(RESULTS / "tuning_summary.json"),
        info=read_json(a / "analysis.json"),
        by_algorithm=pd.read_csv(a / "by_algorithm.csv"),
        by_instance=pd.read_csv(a / "by_instance.csv"),
        tests=pd.read_csv(a / "tests.csv"),
        scaling=pd.read_csv(a / "scaling.csv"),
        true_reference=pd.read_csv(a / "true_order_reference.csv"),
        stages=read_json(RESULTS / "stages.json"),
    )


@lru_cache(maxsize=1)
def load_runs_frame() -> pd.DataFrame:
    """Every main run (one row per solver and read set) with its gap and biological metrics."""
    records = runner.read_records(RESULTS / "main_runs.jsonl.gz")
    runs, _ = analysis.load_runs(records)
    return runs


def figure_catalog() -> list[dict[str, Any]]:
    """The figures of Phases 2, 3 and 5 that exist on disk, with caption and file names."""
    entries = []
    for stems in (FIGURE_STEMS, EARLIER_STEMS):
        for stem, title in stems.items():
            caption_file = FIGURES / f"{stem}.caption.json"
            if not caption_file.exists():
                continue
            caption = json.loads(caption_file.read_text(encoding="utf-8"))
            entries.append(
                {
                    "stem": stem,
                    "title": title,
                    "caption": caption["caption"],
                    "source": caption["source"],
                    "config_hash": caption["config_hash"],
                    "png": FIGURES / f"{stem}.png",
                    "svg": FIGURES / f"{stem}.svg",
                    "pdf": FIGURES / f"{stem}.pdf",
                }
            )
    return entries


def filter_tests(
    tests: pd.DataFrame,
    *,
    family: str | None = None,
    scope: str | None = None,
    significant_only: bool = False,
    algorithm: str | None = None,
) -> pd.DataFrame:
    """Rows of the statistics table that match every given filter."""
    out = tests
    if family:
        out = out[out["family"] == family]
    if scope:
        out = out[out["scope"] == scope]
    if significant_only:
        out = out[out["significant"]]
    if algorithm:
        out = out[(out["a"] == algorithm) | (out["b"] == algorithm)]
    return out.reset_index(drop=True)


def filter_runs(
    runs: pd.DataFrame, *, instances: list[str] | None = None, algorithms: list[str] | None = None
) -> pd.DataFrame:
    out = runs
    if instances:
        out = out[out["instance"].isin(instances)]
    if algorithms:
        out = out[out["algorithm"].isin(algorithms)]
    return out
