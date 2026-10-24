"""Page 7, Assembly Viewer: dot plot of the assembled text against the reference, identity, contigs."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd  # noqa: E402
import streamlit as st  # noqa: E402
import ui  # noqa: E402

from genomeasm.appkit import figs, lab  # noqa: E402
from genomeasm.solvers.lkh import LkhError  # noqa: E402

st.title("Assembly Viewer: is the rebuilt genome right?")
ui.live_badge()
ui.ensure_ready()
params = ui.instance_controls("view", overlap_options=True)
built = ui.build_instance(params)
if built is None:
    st.stop()
inst = built.instance

c1, c2, c3 = st.columns([2.4, 2.4, 1.2])
choices = {"true": "The true order (the answer key)", **{k: v for k, v in lab.SOLVERS.items()}}
solver = c1.selectbox("Order of the reads from", list(choices), index=1, format_func=choices.get)
policy = "ucb1_rt"
if solver == "alns":
    policy = c2.selectbox(
        "Operator selection", list(lab.POLICIES), index=5, format_func=lab.POLICIES.get
    )
seconds = c3.slider(
    "Budget (s)", 1, 15, 5, disabled=solver in ("true", "greedy", "nearest_neighbour")
)

if st.button("Assemble", type="primary", icon=":material/construction:", key="view-assemble"):
    try:
        with st.spinner("Searching for an order and assembling..."):
            if solver == "true":
                order = inst.true_order()
            else:
                order = lab.run_solver(built, solver, seconds=seconds, policy=policy).order
            assembly = inst.assemble(order)
            st.session_state["view_run"] = {
                "params": params,
                "solver": solver,
                "order": [int(i) for i in order],
                "fitness": int(inst.fitness(order)),
                "quality": lab.quality(built, order),
                "dot": lab.dot_plot(built, assembly),
                "contigs": [
                    {
                        "contig": k + 1,
                        "reads": len(c.reads),
                        "letters": len(c.sequence),
                        "starts with": c.sequence[:30],
                    }
                    for k, c in enumerate(assembly.contigs)
                ],
                "fasta": "".join(
                    f">contig_{k + 1}\n{c.sequence}\n" for k, c in enumerate(assembly.contigs)
                ),
                "starts": inst.true_starts,
            }
    except LkhError as error:
        st.error(
            f"LKH-3 could not run: {error}. It is reported as unavailable, never replaced (rule R9)."
        )

run = st.session_state.get("view_run")
if run is None:
    st.info("Choose where the order comes from and press Assemble.")
    st.stop()

q = run["quality"]
st.subheader(f"{choices[run['solver']]}: F = {run['fitness']:,}")
a, b, c, d, e = st.columns(5)
a.metric(
    "Identity", f"{100 * q['identity']:.2f}%", help="Matching letters among the located windows"
)
b.metric(
    "Genome fraction",
    f"{100 * q['genome_fraction']:.1f}%",
    help="Share of the real window rebuilt with identity of at least 90%",
)
c.metric("Contigs", q["n_contigs"], help="Gap-free pieces of the assembly")
d.metric("N50", q["n50"], help="Half of the assembled letters lie in contigs at least this long")
e.metric("Misjoins", q["misjoins"], help="Jumps of a contig to a distant place of the genome")
left, right = st.columns(2)
footer = f"live run, seed {run['params'].seed}"
with left:
    ui.show_figure(
        figs.dot_plot(run["dot"]), "assembly_dot_plot",
        caption="Each dot is a 100-letter window of the assembled text located on the real genome. One straight rising line means a correct assembly; breaks are contig ends.",
        footer=footer,
    )  # fmt: skip
with right:
    ui.show_figure(
        figs.tour_blocks(run["order"], run["starts"]), "assembly_order_vs_truth",
        caption=f"Rank correlation of the order with the true order: {q['spearman']:.3f}; correct neighbours: {100 * q['adjacent_correct']:.1f}%.",
        footer=footer,
    )  # fmt: skip
st.markdown("**Contigs**")
ui.table(pd.DataFrame(run["contigs"]))
st.download_button(
    "Download the assembled contigs (FASTA)", run["fasta"], "assembled_contigs.fasta", "text/plain",
    icon=":material/download:", on_click="ignore",
)  # fmt: skip
st.caption(
    "Read errors break the chain into several contigs even for the true order (decision D-031); the metrics against the real genome are what counts, not F alone."
)
