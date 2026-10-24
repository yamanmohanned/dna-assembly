"""Statistics and tables of the Phase 5 experiments (rule R7), computed only from the raw runs.

Fitness is compared as the *relative gap to the best-known fitness of the same read set*:
``gap = 100 * (best - F) / best``, where ``best`` is the largest fitness any solver reached on that
read set (so gap 0 means "as good as the best known"). Read sets are the units of replication:
solvers are compared **paired** on the same read set. Tests (two-sided, ``alpha`` in the config):

* Wilcoxon signed-rank test on the paired differences of relative fitness (zero differences are
  discarded; if every difference is zero the p-value is 1),
* effect sizes: Vargha-Delaney ``A12`` = P(A > B) + 0.5 P(A = B) over all pairs of observations,
  Cliff's ``delta = 2 * A12 - 1``, and the paired win rate P(d > 0) + 0.5 P(d = 0),
* bootstrap 95% confidence interval of the **median paired difference** (percentile method,
  ``bootstrap_resamples`` resamples, fixed seed),
* Holm correction within every family of comparisons (family and scope are columns of the table).

A positive ``median_diff_pp`` means the first solver of the pair is better (smaller gap) by that
many percentage points of fitness.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from genomeasm.experiments.plan import FIXED_ALGORITHMS, MAIN_ALGORITHMS

BIO_FIELDS = [
    "identity", "genome_fraction", "n_contigs", "n50", "assembled_length", "largest_contig",
    "misjoins", "unmapped_fraction", "spearman", "spearman_contigs", "adjacent_correct",
]  # fmt: skip


RUN_KEYS = ("stage", "instance", "seed_index", "n", "fitness", "evaluations", "seconds")
INSTANCE_KEYS = (
    "stage",
    "instance",
    "seed_index",
    "n",
    "n_input",
    "edges",
    "build_seconds",
    "F_true",
)


def load_runs(records: list[dict[str, Any]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(runs, instances) tables with biological metrics flattened into columns."""
    runs, inst = [], []
    for r in records:
        if r.get("kind") == "run":
            row = {k: r[k] for k in RUN_KEYS}
            row["algorithm"] = r["requested"]
            row.update({f"bio_{k}": v for k, v in r["bio"].items()})
            row["info"], row["trace"] = r["info"], r["trace"]
            runs.append(row)
        elif r.get("kind") == "instance":
            row = {k: r[k] for k in INSTANCE_KEYS}
            row.update({f"bio_{k}": v for k, v in r["bio_true"].items()})
            inst.append(row)
    runs_df, inst_df = pd.DataFrame(runs), pd.DataFrame(inst)
    for df in (runs_df, inst_df):
        df["task"] = df["instance"] + "#" + df["seed_index"].astype(str)
    best = runs_df.groupby("task")["fitness"].transform("max")
    runs_df["best_known"] = best
    runs_df["gap"] = 100.0 * (best - runs_df["fitness"]) / best
    runs_df["rel"] = runs_df["fitness"] / best
    lkh = runs_df[runs_df["algorithm"] == "lkh3"].set_index("task")["fitness"]
    runs_df["gap_lkh"] = (
        100.0 * (runs_df["task"].map(lkh) - runs_df["fitness"]) / runs_df["task"].map(lkh)
    )
    return runs_df, inst_df


def holm(p_values: list[float]) -> list[float]:
    """Holm-Bonferroni adjusted p-values (step-down, monotone)."""
    order = np.argsort(p_values)
    adjusted = np.empty(len(p_values))
    running = 0.0
    for rank, idx in enumerate(order):
        running = max(running, min(1.0, (len(p_values) - rank) * p_values[idx]))
        adjusted[idx] = running
    return [float(x) for x in adjusted]


def a12(x: np.ndarray, y: np.ndarray) -> float:
    """Vargha-Delaney A12: probability that a value of x exceeds one of y (ties count half)."""
    greater = (x[:, None] > y[None, :]).mean()
    equal = (x[:, None] == y[None, :]).mean()
    return float(greater + 0.5 * equal)


def paired_test(a: np.ndarray, b: np.ndarray, resamples: int, seed: int) -> dict[str, float]:
    """Paired comparison of two samples of relative fitness (higher is better)."""
    d = a - b
    if np.any(d != 0):
        w = stats.wilcoxon(d, zero_method="wilcox", alternative="two-sided")
        statistic, p = float(w.statistic), float(w.pvalue)
    else:
        statistic, p = 0.0, 1.0
    rng = np.random.default_rng(seed)
    boots = np.median(d[rng.integers(0, d.size, size=(resamples, d.size))], axis=1)
    lo, hi = np.percentile(boots, [2.5, 97.5])
    ratio = a12(a, b)
    return {
        "n": int(d.size),
        "median_diff_pp": float(100 * np.median(d)),
        "ci_low_pp": float(100 * lo),
        "ci_high_pp": float(100 * hi),
        "wilcoxon_W": statistic,
        "p": p,
        "a12": ratio,
        "cliffs_delta": 2 * ratio - 1,
        "paired_win_rate": float(np.mean(d > 0) + 0.5 * np.mean(d == 0)),
        "wins": int(np.sum(d > 0)),
        "ties": int(np.sum(d == 0)),
        "losses": int(np.sum(d < 0)),
    }


def families(cfg: dict[str, Any], best_alns: str) -> dict[str, list[tuple[str, str]]]:
    s = cfg["statistics"]
    main = [[best_alns if x == "alns_best" else x for x in pair] for pair in s["main_pairs"]]
    out = {
        "primary": [tuple(p) for p in s["primary"]],
        "secondary_vs_random": [(a, "alns_random") for a in s["secondary_vs_random"]],
        "bandits_vs_roulette": [(a, "alns_roulette") for a in s["vs_roulette"]],
        "main": [tuple(p) for p in main],
    }
    if s.get("single_operator_vs_random"):
        out["random_vs_single_operator"] = [("alns_random", a) for a in FIXED_ALGORITHMS]
    return out


def run_tests(runs: pd.DataFrame, cfg: dict[str, Any], best_alns: str) -> pd.DataFrame:
    """All paired tests, in every scope (pooled and per instance), Holm-corrected per family."""
    s = cfg["statistics"]
    wide = runs.pivot_table(index=["instance", "seed_index"], columns="algorithm", values="rel")
    rows = []
    scopes = ["pooled", *sorted(runs["instance"].unique())]
    for scope in scopes:
        part = wide if scope == "pooled" else wide.xs(scope, level="instance")
        for family, pairs in families(cfg, best_alns).items():
            block = []
            for a, b in pairs:
                if a in part and b in part:
                    block.append(
                        {
                            "scope": scope,
                            "family": family,
                            "a": a,
                            "b": b,
                            **paired_test(
                                part[a].to_numpy(),
                                part[b].to_numpy(),
                                int(s["bootstrap_resamples"]),
                                int(cfg["seed"]),
                            ),
                        }  # fmt: skip
                    )
            for row, p_adj in zip(block, holm([r["p"] for r in block]), strict=True):
                row["p_holm"] = p_adj
                row["significant"] = bool(p_adj < s["alpha"])
            rows.extend(block)
    return pd.DataFrame(rows)


def summarise(runs: pd.DataFrame, by: list[str]) -> pd.DataFrame:
    """Median and quartiles of the gap and the biological metrics per group."""
    g = runs.groupby(by)
    out = g["gap"].agg(
        n="count", median="median", q1=lambda x: x.quantile(0.25), q3=lambda x: x.quantile(0.75)
    )
    out["mean_fitness"] = g["fitness"].mean()
    out["median_seconds"] = g["seconds"].median()
    out["median_evaluations"] = g["evaluations"].median()
    for field in (
        "genome_fraction",
        "identity",
        "n_contigs",
        "n50",
        "misjoins",
        "adjacent_correct",
    ):
        out[f"bio_{field}"] = g[f"bio_{field}"].median()
    return out.reset_index()


def claim_verdict(tests: pd.DataFrame, alpha: float) -> list[dict[str, Any]]:
    """The pre-registered claim: bandit selection beats random selection (pooled, Holm)."""
    rows = []
    for r in tests[(tests["scope"] == "pooled") & (tests["family"] == "primary")].itertuples():
        zero = "; the median paired difference is 0, the effect is small"
        if r.significant and r.median_diff_pp > 0:
            verdict = "supported"
        elif r.significant and r.median_diff_pp < 0:
            verdict = "rejected (random is significantly better)"
        elif r.significant and r.a12 > 0.5:  # median difference 0: the effect size decides
            verdict = "supported" + zero
        elif r.significant and r.a12 < 0.5:
            verdict = "rejected (random is significantly better)" + zero
        elif r.significant:
            verdict = "significant but without a direction (A12 = 0.5)"
        else:
            verdict = "not supported (no significant difference)"
        rows.append({"policy": r.a, "baseline": r.b, "median_diff_pp": r.median_diff_pp,
                     "ci": [r.ci_low_pp, r.ci_high_pp], "p_holm": r.p_holm, "a12": r.a12,
                     "alpha": alpha, "verdict": verdict})  # fmt: skip
    return rows


# ---- convergence ---------------------------------------------------------------------------------


def _step(trace: list[list[float]], axis: int, grid: np.ndarray, offset: float = 0.0) -> np.ndarray:
    """Best-so-far value of a trace on a grid (NaN before the first point)."""
    xs = np.array([p[axis] for p in trace], dtype=float) + offset
    ys = np.array([p[2] for p in trace], dtype=float)
    idx = np.searchsorted(xs, grid, side="right") - 1
    return np.where(idx >= 0, ys[np.clip(idx, 0, None)], np.nan)


def convergence(runs: pd.DataFrame, time_grid: np.ndarray, eval_grid: np.ndarray) -> dict[str, Any]:
    """Median and quartiles over read sets of the relative fitness along time and evaluations."""
    out: dict[str, Any] = {
        "time": time_grid.tolist(),
        "evaluations": eval_grid.tolist(),
        "curves": {},
    }
    for name in MAIN_ALGORITHMS:
        part = runs[runs["algorithm"] == name]
        by_time, by_eval = [], []
        for r in part.itertuples():
            offset = 0.0
            if name == "lkh3":  # LKH-3's own clock starts after its preprocessing
                offset = float(r.seconds) - float(r.info.get("search_seconds", r.seconds))
            by_time.append(_step(r.trace, 0, time_grid, offset) / r.best_known)
            by_eval.append(_step(r.trace, 1, eval_grid) / r.best_known)
        entry: dict[str, Any] = {}
        for key, rows in (("time", by_time), ("evaluations", by_eval)):
            arr = np.array(rows)
            with np.errstate(all="ignore"):
                entry[key] = {
                    "median": np.nanmedian(arr, axis=0).tolist(),
                    "q1": np.nanpercentile(arr, 25, axis=0).tolist(),
                    "q3": np.nanpercentile(arr, 75, axis=0).tolist(),
                }
        out["curves"][name] = entry
    return out


def time_to_target(trace: list[list[float]], target: float, offset: float = 0.0) -> float | None:
    """First time at which the best-so-far fitness reaches ``target`` (None if never)."""
    for t, _, f in trace:
        if f >= target:
            return float(t) + offset
    return None


def scaling_table(runs: pd.DataFrame, inst: pd.DataFrame, budget: float) -> pd.DataFrame:
    """Runtime against the number of reads: build, greedy, LKH-3 wall-clock, time to near-best."""
    rows = []
    for instance, part in runs.groupby("instance"):
        row: dict[str, Any] = {"instance": instance, "n": float(part["n"].median())}
        row["build_seconds"] = float(inst[inst["instance"] == instance]["build_seconds"].median())
        for name in ("greedy", "lkh3"):
            row[f"{name}_seconds"] = float(part[part["algorithm"] == name]["seconds"].median())
        for name in ("tabu", "alns_random", "alns_ucb1"):
            times, censored = [], 0
            for r in part[part["algorithm"] == name].itertuples():
                t = time_to_target(r.trace, (1 - 0.001) * r.best_known)
                censored += t is None
                times.append(budget if t is None else t)
            row[f"{name}_time_to_0.1pct"] = float(np.median(times))
            row[f"{name}_censored"] = censored
            row[f"{name}_final_gap"] = float(part[part["algorithm"] == name]["gap"].median())
        rows.append(row)
    return pd.DataFrame(rows).sort_values("n")


def markdown(df: pd.DataFrame, floatfmt: str = ".3f") -> str:
    """A pipe table without extra dependencies."""
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for row in df.itertuples(index=False):
        cells = [format(v, floatfmt) if isinstance(v, float) else str(v) for v in row]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines) + "\n"


def write_tables(tables: dict[str, pd.DataFrame], out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, df in tables.items():
        df.to_csv(out_dir / f"{name}.csv", index=False, lineterminator="\n")
