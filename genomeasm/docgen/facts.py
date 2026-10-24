"""The numbers behind the report, the slides and the oral-exam sheet, read from ``results/``."""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import cached_property
from typing import Any

import pandas as pd

from genomeasm.appkit import lab, results_view
from genomeasm.appkit.paths import RESULTS, ROOT
from genomeasm.appkit.results_view import Phase5

NAMES = {
    "greedy": "Greedy",
    "nearest_neighbour": "Nearest neighbour",
    "tabu": "Tabu search",
    "tabu_cl": "Tabu search + candidate lists (post hoc)",
    "alns_random": "ALNS, random choice",
    "alns_roundrobin": "ALNS, round robin",
    "alns_roulette": "ALNS, roulette wheel",
    "alns_epsilon_greedy": "ALNS, epsilon-greedy",
    "alns_ucb1": "ALNS, UCB1",
    "alns_ucb1_rt": "ALNS, UCB1 with random ties (post hoc)",
    "lkh3": "LKH-3",
}


def name(algorithm: str) -> str:
    return NAMES.get(algorithm, algorithm.replace("alns_fixed_", "single operator: "))


def p_text(value: float) -> str:
    return (
        "p < 0.0001" if value < 1e-4 else f"p = {value:.4f}" if value < 0.1 else f"p = {value:.2f}"
    )


@dataclass
class Facts:
    p5: Phase5
    runs: pd.DataFrame
    config: dict[str, Any]
    timing: dict[str, Any]

    @classmethod
    def load(cls) -> Facts:
        timing = json.loads(
            (ROOT / "results" / "phase6" / "interactive_timing.json").read_text("utf-8")
        )
        return cls(results_view.load_phase5(), results_view.load_runs_frame(), lab.config(), timing)

    # ---- sizes -----------------------------------------------------------------------------------
    @property
    def n_sets(self) -> int:
        return int(self.p5.info["read_sets"])

    @property
    def n_settings(self) -> int:
        return len(self.config["evaluation_instances"])

    @property
    def seeds(self) -> int:
        return int(self.config["seeds"]["evaluation"])

    @property
    def budget(self) -> float:
        return float(self.config["budget_seconds"])

    @property
    def n_algorithms(self) -> int:
        return int(self.runs["algorithm"].nunique())

    # ---- results ---------------------------------------------------------------------------------
    @cached_property
    def by_algorithm(self) -> pd.DataFrame:
        return self.p5.by_algorithm.set_index("algorithm")

    def mean_gap(self, algorithm: str) -> float:
        return float(self.by_algorithm.loc[algorithm, "mean_gap"])

    def median_gap(self, algorithm: str) -> float:
        return float(self.by_algorithm.loc[algorithm, "median"])

    def test(self, a: str, b: str, scope: str = "pooled") -> pd.Series:
        t = self.p5.tests
        rows = t[(t["scope"] == scope) & (t["a"] == a) & (t["b"] == b)]
        if rows.empty:
            raise KeyError((a, b, scope))
        return rows.iloc[0]

    @cached_property
    def lkh_never_beaten(self) -> bool:
        pivot = self.runs.pivot_table(index="task", columns="algorithm", values="fitness")
        return bool(pivot.drop(columns="lkh3").le(pivot["lkh3"], axis=0).all().all())

    @cached_property
    def greedy_equals_best(self) -> float:
        pivot = self.runs.pivot_table(index="task", columns="algorithm", values="fitness")
        best = self.runs.groupby("task")["best_known"].first().loc[pivot.index]
        return float((pivot["greedy"] == best).mean())

    @cached_property
    def tabu_equals_greedy(self) -> int:
        pivot = self.runs.pivot_table(index="task", columns="algorithm", values="fitness")
        return int((pivot["tabu"] == pivot["greedy"]).sum())

    @cached_property
    def arm_zero_share(self) -> float:
        """Share of the uses of the first operator pair under the original UCB1."""
        from genomeasm.experiments import runner

        records = runner.read_records(RESULTS / "main_runs.jsonl.gz")
        total = [0.0] * 6
        for r in records:
            if r["kind"] == "run" and r["requested"] == "alns_ucb1":
                for k, row in enumerate(r["info"]["outcomes"]):
                    total[k] += sum(row)
        return total[0] / sum(total)

    @cached_property
    def lkh_seconds(self) -> list[float]:
        from genomeasm.experiments import runner

        records = runner.read_records(RESULTS / "main_runs.jsonl.gz")
        return [r["seconds"] for r in records if r["kind"] == "run" and r["requested"] == "lkh3"]

    @property
    def lkh_over(self) -> float:
        """Share of the LKH-3 runs that took more than 1% longer than the budget."""
        return sum(s > 1.01 * self.budget for s in self.lkh_seconds) / len(self.lkh_seconds)

    @property
    def lkh_max(self) -> float:
        return max(self.lkh_seconds)

    def bio(self, algorithm: str, field: str) -> float:
        return float(self.by_algorithm.loc[algorithm, f"bio_{field}"])

    @property
    def verdict_ucb(self) -> dict[str, Any]:
        return dict(self.p5.info["claim_bandit_beats_random"][0])

    @property
    def verdict_eps(self) -> dict[str, Any]:
        return dict(self.p5.info["claim_bandit_beats_random"][1])

    @property
    def tuned(self) -> dict[str, Any]:
        return dict(self.p5.summary["algorithms"])

    @property
    def worst_live_seconds(self) -> float:
        return float(self.timing["worst_total_seconds"])
