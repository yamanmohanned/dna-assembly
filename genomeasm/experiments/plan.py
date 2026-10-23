"""Planning of the Phase 5 experiments: tuning tasks, tuned parameters, main and scaling tasks."""

from __future__ import annotations

import copy
from typing import Any

import numpy as np
import pandas as pd

from genomeasm.experiments.design import FIXED_ALGORITHMS, InstanceSpec, seeds_for  # noqa: F401
from genomeasm.experiments.runner import read_records  # noqa: F401  (re-exported for the CLI)

TUNED = [
    "tabu",
    *(f"alns_{p}" for p in ("random", "roundrobin", "roulette", "epsilon_greedy", "ucb1")),
]
ALGORITHMS_TUNED = [*TUNED, "lkh3"]


def default_params(cfg: dict[str, Any], algorithm: str) -> dict[str, Any]:
    d = cfg["tuning"]["defaults"]
    if algorithm == "tabu":
        return {"tabu": dict(d["tabu"])}
    if algorithm == "lkh3":
        return {"lkh": {}}
    policy = algorithm.removeprefix("alns_")
    return {"alns": dict(d["alns"]), "policy": {policy: dict(d["policy"][policy])}}


def sample_params(cfg: dict[str, Any], algorithm: str, index: int, trial: int) -> dict[str, Any]:
    """Random configuration of the search space (trial 0 is the default configuration)."""
    if trial == 0:
        return default_params(cfg, algorithm)
    rng = np.random.default_rng(np.random.SeedSequence([int(cfg["seed"]), 99, index, trial]))
    space = cfg["tuning"]["space"]

    def pick(options: list[Any]) -> Any:
        return options[int(rng.integers(len(options)))]

    if algorithm == "tabu":
        return {"tabu": {k: pick(v) for k, v in space["tabu"].items()}}
    policy = algorithm.removeprefix("alns_")
    return {
        "alns": {k: pick(v) for k, v in space["alns"].items()},
        "policy": {policy: {k: pick(v) for k, v in space["policy"].get(policy, {}).items()}},
    }


def trial_list(cfg: dict[str, Any], algorithm: str) -> list[dict[str, Any]]:
    if algorithm == "lkh3":
        return [{"lkh": dict(v)} for v in cfg["tuning"]["lkh_variants"].values()]
    index = ALGORITHMS_TUNED.index(algorithm)
    return [sample_params(cfg, algorithm, index, t) for t in range(int(cfg["tuning"]["trials"]))]


def _task(
    cfg: dict[str, Any],
    stage: str,
    spec_index: int,
    spec: dict[str, Any],
    seed_index: int,
    group: str,
    runs: list[dict[str, Any]],
    seconds: float,
    **extra: Any,
) -> dict[str, Any]:
    read_seed, solver_seed = seeds_for(cfg, stage, spec_index, seed_index)
    return {
        "stage": stage,
        "spec": spec,
        "seed_index": seed_index,
        "read_seed": read_seed,
        "solver_seed": solver_seed,
        "seconds": seconds,
        "runs": runs,
        "group": group,
        "config": cfg,
        **extra,
    }


def tuning_tasks(
    cfg: dict[str, Any], seeds: int | None = None, trials: int | None = None
) -> list[dict[str, Any]]:
    cfg = copy.deepcopy(cfg)
    if trials is not None:
        cfg["tuning"]["trials"] = trials
    tasks = []
    for algorithm in ALGORITHMS_TUNED:
        trial_params = trial_list(cfg, algorithm)
        runs = [
            {"algorithm": algorithm, "trial": t, "params": p} for t, p in enumerate(trial_params)
        ]
        for si, spec in enumerate(cfg["tuning_instances"]):
            for s in range(seeds if seeds is not None else int(cfg["seeds"]["tuning"])):
                tasks.append(
                    _task(
                        cfg,
                        "tuning",
                        si,
                        spec,
                        1000 + s,
                        algorithm,
                        runs,
                        float(cfg["tuning_budget_seconds"]),
                        bio=False,
                    )  # fmt: skip
                )
    return tasks


def regret_table(records: list[dict[str, Any]]) -> pd.DataFrame:
    """Mean regret to the best-known fitness (per read set) of every algorithm and trial."""
    df = pd.DataFrame([r for r in records if r.get("kind") == "run"])
    df["task"] = df["instance"] + "#" + df["seed_index"].astype(str)
    best = df.groupby("task")["fitness"].transform("max")
    df["regret"] = 100.0 * (best - df["fitness"]) / best
    table = df.groupby(["requested", "trial"])["regret"].agg(["mean", "std", "count"]).reset_index()
    return table.rename(columns={"requested": "algorithm"})


def select_tuned(cfg: dict[str, Any], records: list[dict[str, Any]]) -> dict[str, Any]:
    """The best configuration of every algorithm and the best ALNS policy (tuning data only)."""
    table = regret_table(records)
    summary: dict[str, Any] = {"algorithms": {}}
    for algorithm in ALGORITHMS_TUNED:
        rows = table[table["algorithm"] == algorithm].sort_values(["mean", "trial"])
        trials = trial_list(cfg, algorithm)
        best_trial = int(rows.iloc[0]["trial"])
        default_row = rows[rows["trial"] == 0]
        summary["algorithms"][algorithm] = {
            "best_trial": best_trial,
            "best_mean_regret": float(rows.iloc[0]["mean"]),
            "default_mean_regret": float(default_row.iloc[0]["mean"]) if len(default_row) else None,
            "params": trials[best_trial],
            "trials": [
                {
                    "trial": int(r.trial),
                    "mean_regret": float(r.mean),
                    "params": trials[int(r.trial)],
                }
                for r in table[table["algorithm"] == algorithm].itertuples()
            ],
        }
    alns_rows = {a: summary["algorithms"][a]["best_mean_regret"] for a in TUNED if a != "tabu"}
    summary["best_alns_policy"] = min(alns_rows, key=lambda a: alns_rows[a])
    summary["alns_regret_by_policy"] = alns_rows
    return summary


def tuned_params(summary: dict[str, Any], algorithm: str) -> dict[str, Any]:
    """Parameters for an algorithm in the main experiment (single operators reuse alns_random)."""
    if algorithm in FIXED_ALGORITHMS:
        return {"alns": summary["algorithms"]["alns_random"]["params"]["alns"], "policy": {}}
    return summary["algorithms"][algorithm]["params"]


MAIN_ALGORITHMS = [
    "greedy",
    "nearest_neighbour",
    "tabu",
    *(f"alns_{p}" for p in ("random", "roundrobin", "roulette", "epsilon_greedy", "ucb1")),
    *FIXED_ALGORITHMS,
    "lkh3",
]


def _params_for(summary: dict[str, Any], algorithm: str) -> dict[str, Any]:
    if algorithm in ("greedy", "nearest_neighbour"):
        return {}
    return tuned_params(summary, algorithm)


def main_tasks(
    cfg: dict[str, Any],
    summary: dict[str, Any],
    stage: str = "main",
    seeds: int | None = None,
    specs: list[dict[str, Any]] | None = None,
    algorithms: list[str] | None = None,
    **extra: Any,
) -> list[dict[str, Any]]:
    """One task per read set; the order of the timed solvers is shuffled per task."""
    if specs is None:
        specs = cfg["evaluation_instances"]
    names = algorithms or MAIN_ALGORITHMS
    n_seeds = seeds if seeds is not None else int(cfg["seeds"]["evaluation"])
    tasks = []
    for si, spec in enumerate(specs):
        for s in range(n_seeds):
            rng = np.random.default_rng(np.random.SeedSequence([int(cfg["seed"]), 7, si, s]))
            quick = [a for a in names if a in ("greedy", "nearest_neighbour")]
            timed = [a for a in names if a not in quick]
            order = quick + [timed[i] for i in rng.permutation(len(timed))]
            runs = [{"algorithm": a, "params": _params_for(summary, a)} for a in order]
            tasks.append(
                _task(cfg, stage, si, spec, s, "all", runs, float(cfg["budget_seconds"]), **extra)
            )
    return tasks


def missing_lkh_tasks(
    cfg: dict[str, Any], summary: dict[str, Any], records: list[dict[str, Any]], stage: str = "main"
) -> list[dict[str, Any]]:
    """Tasks that run only LKH-3 on the read sets where it has no result (e.g. WSL crashed).

    The read sets are rebuilt from the same seeds, so they are identical to the original ones.
    """
    have = {
        (r["instance"], r["seed_index"])
        for r in records
        if r.get("kind") == "run" and r.get("requested") == "lkh3" and r["stage"] == stage
    }
    specs = cfg["evaluation_instances"] if stage == "main" else cfg["scaling_instances"]
    tasks = []
    for si, spec in enumerate(specs):
        for s in range(int(cfg["seeds"]["evaluation" if stage == "main" else "scaling"])):
            if (spec["id"], s) not in have:
                run = {"algorithm": "lkh3", "params": _params_for(summary, "lkh3")}
                extra: dict[str, Any] = {"emit_instance": False}
                if stage == "scaling":
                    extra["lkh_search_seconds"] = float(cfg["lkh_scaling_seconds"])
                tasks.append(
                    _task(
                        cfg,
                        stage,
                        si,
                        spec,
                        s,
                        "lkh_rerun",
                        [run],
                        float(cfg["budget_seconds"]),
                        **extra,
                    )
                )
    return tasks


def scaling_tasks(
    cfg: dict[str, Any], summary: dict[str, Any], seeds: int | None = None
) -> list[dict[str, Any]]:
    return main_tasks(
        cfg,
        summary,
        stage="scaling",
        seeds=seeds if seeds is not None else int(cfg["seeds"]["scaling"]),
        specs=cfg["scaling_instances"],
        algorithms=list(cfg["scaling_solvers"]),
        lkh_search_seconds=float(cfg["lkh_scaling_seconds"]),
    )
