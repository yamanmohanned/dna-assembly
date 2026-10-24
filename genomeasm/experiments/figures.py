"""Figures 4 to 11 of the required list, drawn only from the raw results of Phase 5."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.ticker import NullFormatter

from genomeasm.experiments.plan import MAIN_ALGORITHMS
from genomeasm.plots import OKABE_ITO as C
from genomeasm.plots import apply_style, save_figure
from genomeasm.solvers.alns import ARM_NAMES, OUTCOMES

COLORS = {
    "greedy": "#777777", "nearest_neighbour": "#BBBBBB", "tabu": C["blue"], "tabu_cl": "#0B3C5D",
    "alns_random": C["vermillion"], "alns_roundrobin": C["sky"], "alns_roulette": C["orange"],
    "alns_epsilon_greedy": C["green"], "alns_ucb1": C["purple"], "alns_ucb1_rt": "#882255", "lkh3": C["black"],
}  # fmt: skip
LABELS = {
    "greedy": "Greedy", "nearest_neighbour": "Nearest neighbour", "tabu": "Tabu search",
    "tabu_cl": "Tabu search + candidate lists",
    "alns_random": "ALNS random", "alns_roundrobin": "ALNS round robin",
    "alns_roulette": "ALNS roulette wheel", "alns_epsilon_greedy": "ALNS epsilon-greedy",
    "alns_ucb1": "ALNS UCB1", "alns_ucb1_rt": "ALNS UCB1 (random ties)", "lkh3": "LKH-3",
}  # fmt: skip
MAIN = list(COLORS)


def _label(name: str) -> str:
    return LABELS.get(name, name.replace("alns_fixed_", "single: "))


def _stars(p: float) -> str:
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"


def _gap_curve(curve: dict[str, Any], floor: float = 1e-3) -> np.ndarray:
    """Mean gap (%) of a convergence curve; NaN where fewer than 90% of the read sets have data."""
    gap = 100 * (1 - np.array(curve["mean"]))
    gap = np.where(np.array(curve["valid"]) >= 0.9, gap, np.nan)
    return np.maximum(gap, floor)


def fig_convergence(conv: dict[str, Any], n_sets: int) -> tuple[Figure, str]:
    apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    for ax, key, grid, xlabel in (
        (axes[0], "time", np.array(conv["time"]), "Wall-clock time (s, log scale)"),
        (axes[1], "evaluations", np.array(conv["evaluations"]), "Evaluations (log scale)"),
    ):
        for name in MAIN:
            if name == "nearest_neighbour" or (name == "lkh3" and key == "evaluations"):
                continue
            ax.plot(grid, _gap_curve(conv["curves"][name][key]), color=COLORS[name], lw=2.0,
                    ls="--" if name == "lkh3" else "-", label=_label(name))  # fmt: skip
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Mean gap to the best known (%, log scale)")
    axes[0].legend(fontsize=8.5, loc="lower left")
    axes[0].set_title("(a) against time")
    axes[1].set_title("(b) against evaluations")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    return fig, (
        f"Convergence over {n_sets} read sets (10 instance settings x 30 seeds): mean gap of the best "
        "solution so far to the best fitness any solver reached on the same read set. A curve is drawn "
        "where at least 90% of the read sets already have a solution. LKH-3's curve starts after its "
        "preprocessing, has no evaluation count and sits on the 0.001% floor (its gap is 0 by "
        "construction of the best-known value). Nearest neighbour is left out (gap about 5%)."
    )


def fig_final(runs: pd.DataFrame) -> tuple[Figure, str]:
    apply_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5), gridspec_kw={"width_ratios": [1, 1.25]})
    names = [a for a in MAIN if a in set(runs["algorithm"])]
    data = [runs[runs["algorithm"] == a]["gap"].to_numpy() for a in names]
    box = ax1.boxplot(data, vert=False, patch_artist=True, showfliers=True, widths=0.6,
                      flierprops={"markersize": 2.5, "alpha": 0.4})  # fmt: skip
    for patch, a in zip(box["boxes"], names, strict=True):
        patch.set_facecolor(COLORS[a])
        patch.set_alpha(0.6)
    ax1.set_yticks(range(1, len(names) + 1), [_label(a) for a in names])
    ax1.invert_yaxis()
    ax1.set_xlabel("Gap to the best known fitness (%, lower is better)")
    ax1.set_title(f"(a) all read sets pooled (N = {len(data[0])} per solver)")
    shown = ["greedy", "tabu", "tabu_cl", "alns_random", "alns_ucb1", "lkh3"]
    med = runs.groupby(["instance", "algorithm"])["gap"].median().unstack()
    order = list(runs.drop_duplicates("instance")["instance"])
    for a in shown:
        ax2.plot(range(len(order)), med.loc[order, a], "o-", color=COLORS[a], label=_label(a))
    ax2.set_xticks(range(len(order)), order, rotation=45, ha="right", fontsize=8.5)
    ax2.set_ylabel("Median gap to the best known (%)")
    ax2.set_title("(b) median per instance setting")
    ax2.legend(fontsize=8.5)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    return fig, (
        "Final fitness of all solvers after the same wall-clock budget. (a) Boxes: quartiles, whiskers "
        "1.5 x IQR, points: outliers. (b) Median gap per instance setting. Gap = 100 x (best known "
        "- F) / best known on the same read set."
    )


def fig_gap_lkh(runs: pd.DataFrame) -> tuple[Figure, str]:
    apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
    sets = (
        (axes[0], ["L5-c5-e1.0", "L5-c10-e1.0", "L5-c20-e1.0"], [5, 10, 20], "Coverage (x)"),
        (axes[1], ["L5-c10-e0.5", "L5-c10-e1.0", "L5-c10-e2.0"], [0.5, 1.0, 2.0], "Substitution error rate (%)"),
    )  # fmt: skip
    for ax, ids, xs, xlabel in sets:
        for k, a in enumerate(["greedy", "tabu", "tabu_cl", "alns_random", "alns_ucb1"]):
            q = [runs[(runs["instance"] == i) & (runs["algorithm"] == a)]["gap_lkh"] for i in ids]
            med = np.array([x.median() for x in q])
            lo = med - np.array([x.quantile(0.25) for x in q])
            hi = np.array([x.quantile(0.75) for x in q]) - med
            ax.errorbar(np.array(xs) * (1 + 0.02 * k), med, yerr=[lo, hi], fmt="o-", capsize=3,
                        color=COLORS[a], label=_label(a))  # fmt: skip
        ax.axhline(0, color="#444444", lw=0.8)
        ax.set_xlabel(xlabel)
        ax.set_xscale("log")
        ax.set_xticks(xs, [str(x) for x in xs])
        ax.xaxis.set_minor_formatter(NullFormatter())
    axes[0].set_ylabel("Gap to LKH-3 (%; positive = below LKH-3)")
    axes[0].set_title("(a) against coverage (error 1%)")
    axes[1].set_title("(b) against error rate (coverage 10x)")
    axes[0].legend(fontsize=9)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    return fig, (
        "Gap of each solver to LKH-3 on the same read set, lambda window of 5,000 bases, 30 read sets "
        "per point (median and interquartile range). Gap = 100 x (F of LKH-3 - F) / F of LKH-3: "
        "negative values mean the solver beat LKH-3 within the same total wall-clock time."
    )


def fig_ablation(tests: pd.DataFrame) -> tuple[Figure, str]:
    apply_style()
    pooled = tests[tests["scope"] == "pooled"]
    rows = []
    for r in pooled.itertuples():
        if r.family in ("primary", "secondary_vs_random", "post_hoc") and r.b == "alns_random":
            rows.append((r.a, r.median_diff_pp, r.ci_low_pp, r.ci_high_pp, r.p_holm, r.family))
        elif r.family == "random_vs_single_operator":  # shown as: single operator minus random
            rows.append((r.b, -r.median_diff_pp, -r.ci_high_pp, -r.ci_low_pp, r.p_holm, r.family))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2), gridspec_kw={"width_ratios": [1.3, 1]})
    for y, (_, m, lo, hi, p, fam) in enumerate(rows):
        colour = C["blue"] if fam != "random_vs_single_operator" else C["vermillion"]
        colour = COLORS["alns_ucb1_rt"] if fam == "post_hoc" else colour
        ax1.errorbar(m, y, xerr=[[m - lo], [hi - m]], fmt="o", color=colour, capsize=3)
        ax1.text(hi + 0.0005, y, _stars(p), va="center", fontsize=10)
    ax1.axvline(0, color="#444444", lw=0.8)
    ax1.set_yticks(range(len(rows)), [_label(r[0]) for r in rows])
    ax1.invert_yaxis()
    ax1.set_xlabel("Median fitness difference to ALNS random (percentage points; > 0 = better)")
    ax1.set_title("(a) pooled over all read sets, with 95% bootstrap CI")
    for pol, colour in (
        ("alns_ucb1", COLORS["alns_ucb1"]),
        ("alns_epsilon_greedy", COLORS["alns_epsilon_greedy"]),
        ("alns_ucb1_rt", COLORS["alns_ucb1_rt"]),
    ):
        part = tests[
            (tests["a"] == pol) & (tests["b"] == "alns_random") & (tests["scope"] != "pooled")
        ]
        ax2.plot(
            part["median_diff_pp"].to_numpy(),
            range(len(part)),
            "o",
            color=colour,
            label=_label(pol),
        )
    ids = tests[(tests["scope"] != "pooled")]["scope"].drop_duplicates().tolist()
    ax2.set_yticks(range(len(ids)), ids, fontsize=8.5)
    ax2.invert_yaxis()
    ax2.axvline(0, color="#444444", lw=0.8)
    ax2.set_xlabel("Median difference to ALNS random (pp)")
    ax2.set_title("(b) per instance setting (30 read sets each)")
    ax2.legend(fontsize=9)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    return fig, (
        "Ablation: operator-selection policies against random selection, paired by read set. Significance "
        "marks come from two-sided Wilcoxon signed-rank tests with Holm correction within each family "
        "(*** p<0.001, ** p<0.01, * p<0.05, ns not significant). Red rows: a single operator pair "
        "used alone minus the random mixture. UCB1 with random ties is a post-hoc variant (D-046)."
    )


def fig_arm_probability(runs: pd.DataFrame) -> tuple[Figure, str]:
    apply_style()
    fig, axes = plt.subplots(1, 4, figsize=(17, 4.4), sharey=True)
    palette = [C["blue"], C["sky"], C["green"], C["orange"], C["vermillion"], C["purple"]]
    grid = np.linspace(0, 1, 60)
    policies = ["alns_roulette", "alns_epsilon_greedy", "alns_ucb1", "alns_ucb1_rt"]
    for ax, policy in zip(axes, policies, strict=True):
        stack = []
        for r in runs[runs["algorithm"] == policy].itertuples():
            trace = np.array(r.info["probability_trace"])
            if trace.shape[0] >= 2:
                x = np.linspace(0, 1, trace.shape[0])
                stack.append(
                    np.stack([np.interp(grid, x, trace[:, a]) for a in range(trace.shape[1])])
                )
        mean = np.mean(stack, axis=0)
        ax.stackplot(grid, mean, colors=palette, labels=list(ARM_NAMES), alpha=0.9)
        ax.set_title(f"{_label(policy)} (mean of {len(stack)} runs)", fontsize=11)
        ax.set_xlabel("Fraction of the run")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
    axes[0].set_ylabel("Probability of choosing the operator pair")
    axes[3].legend(
        fontsize=8, loc="center left", bbox_to_anchor=(1.0, 0.5), title="destroy + repair"
    )
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    return fig, (
        "How the learning policies spread their choices over the six operator pairs during a run (mean "
        "over all read sets). Roulette: its weights; epsilon-greedy: epsilon/6 plus the greedy arm; UCB1: "
        "the arm it would pick next (one-hot), so the mean shows how often each arm leads. UCB1 with "
        "lowest-index ties collapses onto the first arm once rewards are all zero (D-046); with random "
        "ties (rightmost) it does not."
    )


def fig_operator_heatmap(runs: pd.DataFrame) -> tuple[Figure, str]:
    apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
    for ax, policy in zip(axes, ["alns_random", "alns_ucb1"], strict=True):
        total = np.zeros((len(ARM_NAMES), len(OUTCOMES)))
        for r in runs[runs["algorithm"] == policy].itertuples():
            total += np.array(r.info["outcomes"])
        share = total / np.maximum(total.sum(axis=1, keepdims=True), 1)
        image = ax.imshow(100 * share, cmap="viridis", aspect="auto", vmin=0, vmax=100)
        for i in range(share.shape[0]):
            for j in range(share.shape[1]):
                ax.text(j, i, f"{100 * share[i, j]:.1f}", ha="center", va="center", fontsize=8,
                        color="white" if share[i, j] < 0.5 else "black")  # fmt: skip
        ax.set_xticks(range(len(OUTCOMES)), OUTCOMES, rotation=30, ha="right")
        ax.set_yticks(
            range(len(ARM_NAMES)),
            [f"{a} ({int(t):,})" for a, t in zip(ARM_NAMES, total.sum(axis=1), strict=True)],
            fontsize=8.5,
        )
        ax.set_title(_label(policy))
        ax.grid(False)
    fig.colorbar(image, ax=axes, label="% of the uses of the operator pair")
    return fig, (
        "What each destroy + repair operator pair achieved, as a share of its uses (all read sets pooled; "
        "the number of uses is in brackets). With random selection every pair is used equally, so the rows "
        "show the operators' intrinsic usefulness; with UCB1 the rows show where the bandit put its effort."
    )


def fig_runtime(table: pd.DataFrame, budget: float) -> tuple[Figure, str]:
    apply_style()
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    n = table["n"].to_numpy()
    series = [
        ("build_seconds", "Overlap graph construction", C["green"], "o"),
        ("greedy_seconds", "Greedy", COLORS["greedy"], "s"),
        ("lkh3_seconds", "LKH-3 total wall-clock (1 s search)", COLORS["lkh3"], "^"),
        ("tabu_time_to_0.1pct", "Tabu: time to within 0.1% of best known", COLORS["tabu"], "D"),
        ("tabu_cl_time_to_0.1pct", "Tabu + candidate lists: same", COLORS["tabu_cl"], "P"),
        ("alns_ucb1_time_to_0.1pct", "ALNS UCB1: same", COLORS["alns_ucb1"], "v"),
    ]
    for col, label, colour, marker in series:
        ax.plot(n, np.maximum(table[col].to_numpy(), 1e-3), marker + "-", color=colour, label=label)
    ax.axhline(budget, color="#999999", ls=":", lw=1)
    ax.text(
        n[0], budget * 1.08, f"budget {budget:g} s (values at the budget are censored)", fontsize=8
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(n, [f"{int(v)}" for v in n])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel("Number of reads n (log scale)")
    ax.set_ylabel("Seconds (log scale)")
    ax.legend(fontsize=8)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    return fig, (
        "Runtime against instance size (lambda windows of 2,500 to 20,000 bases, coverage 10x, 1% errors, "
        "10 read sets per size, medians). LKH-3's preprocessing makes its total time grow much faster "
        "than the other quantities."
    )


def fig_pareto(conv: dict[str, Any]) -> tuple[Figure, str]:
    apply_style()
    fig, ax = plt.subplots(figsize=(8, 5.2))
    grid = np.array(conv["time"])
    points = []
    for name in MAIN:
        gap = _gap_curve(conv["curves"][name]["time"])
        ok = ~np.isnan(gap)
        ax.plot(
            grid[ok], np.maximum(gap[ok], 1e-3), "o-", ms=3, color=COLORS[name], label=_label(name)
        )
        points += [(t, g, name) for t, g in zip(grid[ok], gap[ok], strict=True)]
    front, best = [], np.inf
    for t, g, _ in sorted(points):
        if g < best - 1e-12:
            best = g
            front.append((t, max(g, 1e-3)))
    ax.step(
        *zip(*front, strict=True),
        where="post",
        color="black",
        lw=2.2,
        alpha=0.5,
        label="Pareto front",
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Wall-clock time (s, log scale)")
    ax.set_ylabel("Mean gap to the best known (%, log scale)")
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    return fig, (
        "Quality against time: the mean gap of every solver at every time of the budget. The black "
        "step line is the Pareto front: no solver reaches a smaller gap in less time. Gaps below 0.001% "
        "are drawn at 0.001%."
    )


def make_figures(
    runs: pd.DataFrame, tests: pd.DataFrame, conv: dict[str, Any], scaling: pd.DataFrame,
    budget: float, out_dir: Path, source: str, config_hash: str,
) -> list[Path]:  # fmt: skip
    """Draw and save figures 4 to 11."""
    n_sets = int(runs.drop_duplicates(["instance", "seed_index"]).shape[0])
    jobs = [
        ("fig04_convergence", fig_convergence(conv, n_sets)),
        ("fig05_final_fitness", fig_final(runs)),
        ("fig06_gap_to_lkh3", fig_gap_lkh(runs)),
        ("fig07_ablation", fig_ablation(tests)),
        ("fig08_arm_probability", fig_arm_probability(runs)),
        ("fig09_operator_heatmap", fig_operator_heatmap(runs)),
        ("fig10_runtime_vs_n", fig_runtime(scaling, budget)),
        ("fig11_pareto", fig_pareto(conv)),
    ]
    paths: list[Path] = []
    for stem, (fig, caption) in jobs:
        paths += save_figure(
            fig, stem, out_dir, caption=caption, source=source, config_hash=config_hash
        )
    return paths


__all__ = ["MAIN_ALGORITHMS", "make_figures"]
