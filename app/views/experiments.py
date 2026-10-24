"""Page 6, Experiments: the precomputed results of Phase 5 with filters, charts and statistics."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd  # noqa: E402
import streamlit as st  # noqa: E402
import ui  # noqa: E402

from genomeasm.appkit import figs, results_view  # noqa: E402

st.title("Experiments: what the comparison found")
ui.precomputed_badge("Precomputed (Phase 5): 300 read sets, 8 s per solver, not run live")
if not results_view.available():
    st.error("The Phase 5 results are missing: run `genomeasm experiments analyze` first.")
    st.stop()
p5 = results_view.load_phase5()

LABELS = {
    "greedy": "Greedy", "nearest_neighbour": "Nearest neighbour", "tabu": "Tabu",
    "tabu_cl": "Tabu + candidate lists (post hoc)", "alns_random": "ALNS random",
    "alns_roundrobin": "ALNS round robin", "alns_roulette": "ALNS roulette wheel",
    "alns_epsilon_greedy": "ALNS epsilon-greedy", "alns_ucb1": "ALNS UCB1",
    "alns_ucb1_rt": "ALNS UCB1 random ties (post hoc)", "lkh3": "LKH-3",
}  # fmt: skip


def label(name: str) -> str:
    return LABELS.get(name, name.replace("alns_fixed_", "single operator: "))


verdict_tab, charts_tab, stats_tab, compare_tab, scaling_tab, tuning_tab = st.tabs(
    ["Verdict", "Charts", "Statistics", "Compare solvers", "Scaling", "Tuning"]
)

with verdict_tab:
    st.subheader("Does bandit selection beat random selection?")
    st.markdown(
        "The claim, fixed before the experiments: *bandit-based operator selection finds better "
        "assemblies than random operator selection under an equal budget.* Result of the two "
        "pre-registered tests (all 300 read sets pooled, Holm-corrected):"
    )
    primary = pd.DataFrame(p5.info["claim_bandit_beats_random"])
    primary["policy"] = primary["policy"].map(label)
    primary["baseline"] = primary["baseline"].map(label)
    ui.table(primary[["policy", "baseline", "p_holm", "a12", "verdict"]])
    st.error(
        "Not supported: UCB1 is not different from random, epsilon-greedy is slightly worse.",
        icon=":material/close:",
    )
    st.markdown(
        "**A12** is the chance that the first method beats the second on a read set (0.5 = no difference)."
    )
    st.markdown(
        "**Post hoc, exploratory** (added after the first results; the original UCB1 had a tie-breaking defect, decision D-046):"
    )
    post = pd.DataFrame(p5.info["post_hoc_ucb1_tie_breaking"])
    post["policy"] = post["policy"].map(label)
    post["baseline"] = post["baseline"].map(label)
    ui.table(post[["policy", "baseline", "p_holm", "a12", "verdict"]])
    st.markdown(
        "**Final quality of every solver** (gap to the best solution found, in %; lower is better)"
    )
    summary = p5.by_algorithm[
        ["algorithm", "n", "mean_gap", "median", "q1", "q3", "median_seconds"]
    ].copy()
    summary["algorithm"] = summary["algorithm"].map(label)
    summary.columns = [
        "solver",
        "read sets",
        "mean gap %",
        "median gap %",
        "Q1",
        "Q3",
        "median seconds",
    ]
    ui.table(summary.sort_values("mean gap %"))
    st.caption(
        "LKH-3 was never beaten on any read set. 'Best known' is the best solution any of our solvers found, not a proven optimum."
    )

with charts_tab:
    catalog = results_view.figure_catalog()
    choice = st.selectbox(
        "Figure",
        [e["stem"] for e in catalog],
        format_func=lambda s: next(e["title"] for e in catalog if e["stem"] == s),
    )
    entry = next(e for e in catalog if e["stem"] == choice)
    st.image(str(entry["png"]), width="stretch")
    st.caption(entry["caption"])
    st.caption(f"Source: {entry['source']}, configuration hash {entry['config_hash'][:12]}")
    cols = st.columns([1, 1, 1, 4])
    for column, key, mime in zip(
        cols, ("svg", "png", "pdf"), ("image/svg+xml", "image/png", "application/pdf"), strict=False
    ):
        column.download_button(
            key.upper(), entry[key].read_bytes(), f"{choice}.{key}", mime,
            key=f"exp-{choice}-{key}", icon=":material/download:", on_click="ignore",
        )  # fmt: skip

with stats_tab:
    st.markdown(
        "Every paired test of Phase 5. **Mean/median difference** is in percentage points of fitness "
        "(positive: the first solver is better). **p (Holm)** is corrected within each family of tests."
    )
    families = sorted(p5.tests["family"].unique())
    instance_scopes = sorted(s for s in p5.tests["scope"].unique() if s != "pooled")
    f1, f2, f3, f4 = st.columns([2, 2, 2, 1])
    family = f1.selectbox("Family", ["all", *families])
    scope = f2.selectbox(
        "Scope",
        ["pooled", *instance_scopes],
        help="pooled: all 300 read sets; or one instance setting",
    )
    algorithm = f3.selectbox(
        "Solver involved", ["any", *sorted(set(p5.tests["a"]) | set(p5.tests["b"]))]
    )
    only = f4.toggle("Significant only", False)
    table = results_view.filter_tests(
        p5.tests,
        family=None if family == "all" else family,
        scope=scope,
        significant_only=only,
        algorithm=None if algorithm == "any" else algorithm,
    )
    st.write(f"{len(table)} tests")
    shown = table[
        [
            "family",
            "a",
            "b",
            "mean_diff_pp",
            "median_diff_pp",
            "p",
            "p_holm",
            "a12",
            "wins",
            "ties",
            "losses",
            "significant",
        ]
    ]
    ui.table(
        shown,
        column_config={
            "p": st.column_config.NumberColumn(format="%.2e"),
            "p_holm": st.column_config.NumberColumn(format="%.2e"),
        },
    )
    st.download_button(
        "Download this table (CSV)", shown.to_csv(index=False), "phase5_tests.csv", "text/csv",
        icon=":material/download:", on_click="ignore",
    )  # fmt: skip

with compare_tab:
    runs = results_view.load_runs_frame()
    st.markdown(
        "Pick instance settings and solvers; the chart shows the distribution over the read sets."
    )
    instances = sorted(runs["instance"].unique())
    k1, k2, k3 = st.columns([2, 3, 2])
    chosen_instances = k1.multiselect("Instance settings", instances, default=instances)
    names = sorted(
        runs["algorithm"].unique(),
        key=lambda a: p5.by_algorithm.set_index("algorithm").loc[a, "mean_gap"],
    )
    default = [
        a
        for a in (
            "greedy",
            "tabu_cl",
            "alns_random",
            "alns_roulette",
            "alns_ucb1",
            "alns_ucb1_rt",
            "lkh3",
        )
        if a in names
    ]
    chosen = k2.multiselect("Solvers", names, default=default, format_func=label)
    metric = k3.selectbox(
        "Measure",
        [
            "gap",
            "fitness",
            "seconds",
            "bio_genome_fraction",
            "bio_identity",
            "bio_n_contigs",
            "bio_n50",
        ],
    )
    part = results_view.filter_runs(runs, instances=chosen_instances, algorithms=chosen)
    if part.empty:
        st.info("Select at least one instance setting and one solver.")
    else:
        data = figs.runs_frame_to_series(part, chosen, metric)
        ui.show_figure(
            figs.box_by_algorithm(
                {label(a): v for a, v in data.items()},
                {"gap": "Gap to the best known (%)"}.get(metric, metric),
                f"{len(chosen_instances)} instance setting(s)",
            ),
            "compare_solvers",
            caption=f"{part.drop_duplicates(['instance', 'seed_index']).shape[0]} read sets; boxes show quartiles, whiskers 1.5 x IQR.",
            footer="Source: results/phase5/main_runs.jsonl.gz",
        )
        agg = part.groupby("algorithm")[metric].agg(["count", "mean", "median"]).reset_index()
        agg["algorithm"] = agg["algorithm"].map(label)
        ui.table(agg.rename(columns={"count": "n"}))

with scaling_tab:
    st.markdown(
        "Runtime against the number of reads (lambda windows of 2,500 to 20,000 letters, coverage 10x, 10 read sets each)."
    )
    ui.table(p5.scaling)
    png = results_view.figure_catalog()
    entry = next(e for e in png if e["stem"] == "fig10_runtime_vs_n")
    st.image(str(entry["png"]), width="stretch")
    st.caption(entry["caption"])

with tuning_tab:
    st.markdown(
        "Hyper-parameters were tuned by random search on **separate** instances (phiX174, positions 0 to 2,500), "
        "24 trials per algorithm. Mean regret = percent below the best fitness any trial reached."
    )
    rows = [
        {
            "algorithm": label(a),
            "default regret %": v["default_mean_regret"],
            "tuned regret %": v["best_mean_regret"],
            "best trial": v["best_trial"],
            "tuned parameters": str(v["params"]),
        }
        for a, v in p5.summary["algorithms"].items()
    ]
    ui.table(pd.DataFrame(rows))
    st.caption(
        f"Best ALNS policy on the tuning data: {label(p5.summary['best_alns_policy'])}. The best of 24 noisy trials looks better than it is (winner's curse)."
    )
