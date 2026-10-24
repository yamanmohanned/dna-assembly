"""Page 5, Solver Lab: run a solver with live charts, or step through a tiny instance."""

from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd  # noqa: E402
import streamlit as st  # noqa: E402
import ui  # noqa: E402

from genomeasm.appkit import content, figs, lab  # noqa: E402
from genomeasm.solvers.alns import ARM_NAMES, alns  # noqa: E402
from genomeasm.solvers.common import Budget, Problem  # noqa: E402
from genomeasm.solvers.exact import held_karp  # noqa: E402
from genomeasm.solvers.greedy import greedy  # noqa: E402
from genomeasm.solvers.lkh import LkhError  # noqa: E402

st.title("Solver Lab: watch a solver search")
ui.ensure_ready()
run_tab, step_tab = st.tabs(["Run a solver", "Step through a tiny instance"])

OUTCOME_TEXT = {
    "best": "a new best solution",
    "improved": "better than the current solution",
    "accepted": "worse, but accepted (simulated annealing allows it early in the run)",
    "lateral": "same fitness, a different order (accepted)",
    "unchanged": "the move rebuilt the same order",
    "rejected": "worse and rejected",
}


def _summary(
    result: Any, built: Any, params: lab.InstanceParams, extra: dict[str, Any]
) -> dict[str, Any]:
    order = result.order
    greedy_f = greedy(built.problem).fitness
    return {
        "params": params,
        "fitness": int(result.fitness),
        "merged_length": int(result.merged_length),
        "evaluations": int(result.evaluations),
        "seconds": float(result.seconds),
        "greedy_fitness": int(greedy_f),
        "true_fitness": int(built.instance.fitness(built.instance.true_order())),
        "trace": [(float(t), int(f)) for t, _, f in result.trace],
        "outcomes": result.info.get("outcomes"),
        "iterations": result.info.get("iterations"),
        "order": [int(i) for i in order],
        "quality": lab.quality(built, order),
        **extra,
    }


with run_tab:
    ui.live_badge()
    params = ui.instance_controls("lab")
    c1, c2, c3, c4 = st.columns([2.2, 2.2, 1.4, 1])
    solver = c1.selectbox("Algorithm", list(lab.SOLVERS), index=4, format_func=lab.SOLVERS.get)
    policy = "ucb1_rt"
    if solver == "alns":
        policy = c2.selectbox(
            "Operator selection", list(lab.POLICIES), index=5, format_func=lab.POLICIES.get
        )
    seconds = c3.slider("Time budget (s)", 1, int(lab.MAX_SECONDS), 6)
    solver_seed = int(c4.number_input("Solver seed", 0, 10**6, 1))
    tuned = st.toggle(
        "Use the parameters tuned in Phase 5", True,
        help="Off: the default parameters of the configuration file (trial 0 of the tuning)",
    )  # fmt: skip
    if solver == "lkh3":
        if lab.lkh_available():
            st.info(
                "LKH-3 runs in Ubuntu (WSL); its total time includes file writing and preprocessing."
            )
        else:
            st.warning(
                "LKH-3 is not available on this computer; it will be reported as unavailable, never replaced."
            )
    st.caption(
        f"Every live run is limited to {int(lab.MAX_SECONDS)} s of search and {lab.MAX_READS} reads, so "
        "that it finishes in under 30 s. Larger runs: page Experiments (precomputed)."
    )
    if st.button("Run", type="primary", icon=":material/play_arrow:", key="lab-run"):
        built = ui.build_instance(params)
        if built is not None:
            greedy_f = int(greedy(built.problem).fitness)
            true_f = int(built.instance.fitness(built.instance.true_order()))
            refs = {"greedy": float(greedy_f), "true read order": float(true_f)}
            bar = st.progress(0.0, text="Starting...")
            ph_live = st.empty()
            points: list[tuple[float, float]] = []
            history: list[tuple[float, list[float]]] = []
            next_draw = [0.0]  # the search clock at which the chart is redrawn next

            def on_progress(event: dict[str, Any]) -> None:
                t = float(event["seconds"])
                if event["best"] > 0:
                    points.append((t, float(event["best"])))
                if "probabilities" in event:
                    history.append((t, event["probabilities"]))
                if t >= next_draw[0] or event.get("final"):
                    drawn = time.perf_counter()
                    frame = figs.live_frame(
                        figs.live_dashboard(points, refs, history, float(seconds))
                    )
                    ph_live.image(frame, width="stretch")
                    # a redraw takes time from the search: keep it near 10% of the run
                    next_draw[0] = t + max(1.0, 8 * (time.perf_counter() - drawn))
                bar.progress(
                    min(t / seconds, 1.0),
                    text=f"{t:.1f} s of {seconds} s, best F = {int(event['best']):,}",
                )

            started = time.perf_counter()
            try:
                result = lab.run_solver(
                    built, solver, seconds=seconds, seed=solver_seed, policy=policy, tuned=tuned,
                    progress=on_progress,
                )  # fmt: skip
            except LkhError as error:
                bar.empty()
                st.error(f"LKH-3 could not run: {error}. It is reported as unavailable (rule R9).")
                st.session_state.pop("lab_run", None)
            else:
                bar.empty()
                st.session_state["lab_run"] = _summary(
                    result, built, params,
                    {"solver": solver, "policy": policy, "history": history, "refs": refs,
                     "wall": time.perf_counter() - started, "n": built.problem.n},
                )  # fmt: skip
                ph_live.empty()
    run = st.session_state.get("lab_run")
    if run is not None:
        st.divider()
        name = lab.SOLVERS[run["solver"]] + (
            f", {lab.POLICIES[run['policy']]}" if run["solver"] == "alns" else ""
        )
        st.subheader(f"Result: {name}")
        a, b, c, d = st.columns(4)
        a.metric(
            "Fitness F",
            f"{run['fitness']:,}",
            delta=run["fitness"] - run["greedy_fitness"],
            help="Delta: against greedy",
        )
        b.metric("Merged length L", f"{run['merged_length']:,}")
        c.metric("Evaluations", f"{run['evaluations']:,}")
        d.metric(
            "Time",
            f"{run['seconds']:.1f} s",
            help=f"Wall-clock of the whole run including the instance: {run['wall']:.1f} s",
        )
        st.caption(
            f"{run['n']} reads. Greedy F = {run['greedy_fitness']:,}; true read order F = {run['true_fitness']:,} "
            "(with read errors a solver can beat the true order: F is not the genome, see Learn)."
        )
        left, right = st.columns(2)
        footer = f"live run, seed {run['params'].seed}, solver seed"
        with left:
            points = [*run["trace"], (run["seconds"], run["fitness"])]
            ui.show_figure(
                figs.live_convergence(points, run["refs"], max(run["seconds"], 1.0)),
                "lab_convergence", caption="Best fitness found so far.", footer=footer,
            )  # fmt: skip
        with right:
            if run["solver"] == "alns" and len(run["history"]) >= 2:
                ui.show_figure(
                    figs.operator_probabilities(run["history"], max(run["seconds"], 1.0)),
                    "lab_operator_probabilities",
                    caption="Chance of choosing each destroy + repair pair (for UCB1: the pair it would pick next).",
                    footer=footer,
                )
            elif run["outcomes"] is None:
                st.info("This algorithm has no operators to show.")
        if run["outcomes"] is not None:
            ui.show_figure(
                figs.outcome_bars(run["outcomes"]), "lab_outcomes",
                caption=f"{run['iterations']:,} iterations. Green: new best; blue: improved; orange: accepted although worse.",
                footer=footer,
            )  # fmt: skip
        quality = run["quality"]
        st.markdown("**Quality against the real genome**")
        ui.table(
            pd.DataFrame(
                {
                    "genome fraction": [quality["genome_fraction"]],
                    "identity": [quality["identity"]],
                    "contigs": [quality["n_contigs"]],
                    "N50": [quality["n50"]],
                    "misjoins": [quality["misjoins"]],
                    "correct neighbours": [quality["adjacent_correct"]],
                }
            )
        )

with step_tab:
    ui.live_badge("Live, tiny instance")
    st.markdown(
        "Watch one algorithm make its decisions on a handful of reads. **Greedy** accepts or rejects "
        "overlaps one by one; **ALNS** breaks the current order, repairs it with a chosen tool and "
        "decides whether to keep the result."
    )
    s1, s2, s3 = st.columns(3)
    instance_kind = s1.radio(
        "Instance", ["Toy example (3 reads)", "Tiny simulated (7 to 8 reads)"], key="step_kind"
    )
    algorithm = s2.radio("Algorithm", ["Greedy", "ALNS"], key="step_algo")
    seed = int(s3.number_input("Seed", 0, 1000, 1, key="step_seed"))
    if instance_kind.startswith("Toy"):
        toy = content.make_toy(list(content.TOY_READS), 1)
        inst, labels = toy.instance, toy.labels
        problem = Problem.from_instance(inst)
    else:
        built_tiny = lab.tiny_instance(seed)
        inst, problem = built_tiny.instance, built_tiny.problem
        labels = [chr(ord("A") + i) for i in range(inst.n)]
        st.caption(f"{inst.n} reads from a window of phiX174, exact overlaps, seed {seed}.")
    best_f = int(held_karp(problem).fitness)
    if algorithm == "Greedy":
        steps, order = content.greedy_steps(inst.ov)
        k = st.slider("Decision", 0, len(steps), len(steps), key="greedy_step")
        taken = [s for s in steps[:k] if s["verdict"] == "accepted"]
        left, right = st.columns([1, 1])
        with left:
            ui.show_figure(
                figs.toy_graph(inst, labels, None, [(s["from"], s["to"]) for s in taken]),
                "greedy_steps_graph", caption="Red: overlaps accepted so far.",
            )  # fmt: skip
        with right:
            ui.table(
                [
                    {
                        "#": i + 1,
                        "overlap": f"{labels[s['from']]} -> {labels[s['to']]}",
                        "letters": s["overlap"],
                        "decision": s["verdict"],
                    }
                    for i, s in enumerate(steps[:k])
                ][::-1]
            )
        if k == len(steps):
            st.success(
                f"Greedy order: {''.join(labels[i] for i in order)} with F = {inst.fitness(order)}. "
                f"Best possible F = {best_f}."
            )
    elif instance_kind.startswith("Toy"):
        st.info("ALNS needs at least 4 reads: choose the tiny simulated instance.")
    else:
        log: list[dict[str, Any]] = []
        lab_budget = Budget(evaluations=6000)
        result = alns(
            problem, lab_budget, seed=seed, policy="ucb1_rt", step_log=log, calibration_moves=20
        )
        shown = log[:80]
        if not shown:
            st.warning("No ALNS iteration was recorded.")
            st.stop()
        k = st.slider("ALNS iteration", 1, len(shown), 1, key="alns_step")
        s = shown[k - 1]
        a, b, c, d = st.columns(4)
        a.metric("Chosen pair", s["arm_name"], help="destroy tool + repair tool")
        b.metric(
            "Fitness before -> after",
            f"{s['before_f']} -> {s['after_f']}",
            delta=s["after_f"] - s["before_f"],
        )
        c.metric("Best so far", s["best_f"], help=f"Best possible: {best_f}")
        d.metric(
            "Temperature",
            f"{s['temperature']:.2f}",
            help="Higher: worse moves are accepted more easily",
        )
        st.markdown(
            f"**Outcome:** {OUTCOME_TEXT[s['outcome']]}. "
            + ("Kept." if s["accepted"] else "Discarded.")
        )
        ui.show_figure(
            figs.tour_strip(s["before_order"], s["after_order"], labels, f"Iteration {k}: the order before and after"),
            "alns_step_orders", caption="Outlined places changed. The cycle is read from left to right.",
        )  # fmt: skip
        ui.table(
            pd.DataFrame({"pair": list(ARM_NAMES), "chance of being chosen": s["probabilities"]}),
            column_config={
                "chance of being chosen": st.column_config.ProgressColumn(
                    min_value=0.0, max_value=1.0, format="%.2f"
                )
            },
        )
        st.caption(
            f"The whole run ({result.info['iterations']} iterations) ended with F = {result.fitness}; the best possible is {best_f}. Only the first 80 iterations are shown."
        )
