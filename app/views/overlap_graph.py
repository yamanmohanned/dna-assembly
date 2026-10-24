"""Page 4, Overlap Graph: the overlap heatmap, the overlap-length and degree distributions."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
import streamlit as st  # noqa: E402
import ui  # noqa: E402

from genomeasm.appkit import figs  # noqa: E402

st.title("Overlap graph: which reads fit together?")
ui.live_badge()
params = ui.instance_controls("graph", overlap_options=True)
built = ui.build_instance(params)
if built is None:
    st.stop()
inst = built.instance

st.markdown(
    "Every read is a node; an arrow from read *i* to read *j* means *the end of i equals the beginning "
    "of j* for at least tau letters. With read errors we allow a few mismatches. The **overlap matrix** "
    "holds the arrow lengths (0 = no arrow)."
)
edges = int((inst.ov > 0).sum())
true_f = inst.fitness(inst.true_order())
m1, m2, m3, m4 = st.columns(4)
m1.metric(
    "Reads in the graph",
    inst.n,
    help=f"{inst.n_input} simulated; {inst.n_duplicates} duplicates and {inst.n_contained} contained reads removed",
)
m2.metric("Arrows (overlaps)", f"{edges:,}")
m3.metric("Arrows per read", f"{edges / max(inst.n, 1):.1f}")
m4.metric(
    "F of the true order",
    f"{true_f:,}",
    help="Letters saved if the reads were glued in their real order",
)
st.caption(
    f"Built in {built.build_seconds:.1f} s. Minimum overlap tau = {inst.min_overlap}, "
    f"mismatch tolerance {100 * inst.tolerance:.1f}%."
)

heat, lengths, degrees = st.tabs(["Heatmap", "Overlap lengths", "Degrees"])
with heat:
    max_n = st.slider(
        "Reads shown (first by true position)", 20, min(300, max(inst.n, 21)), min(150, inst.n)
    )
    ui.show_figure(
        figs.overlap_heatmap(inst, max_n),
        "overlap_heatmap",
        caption=(
            "Reads sorted by their true position on the genome: the bright band along the diagonal is "
            "the overlaps between true neighbours; scattered points elsewhere are chance overlaps."
        ),
        footer=f"live run, seed {params.seed}",
    )
with lengths:
    ui.show_figure(
        figs.overlap_lengths(inst),
        "overlap_lengths",
        caption="Log scale on the vertical axis. Many short overlaps arise by chance; long overlaps are real neighbours.",
        footer=f"live run, seed {params.seed}",
    )
with degrees:
    ui.show_figure(
        figs.degree_distribution(inst),
        "overlap_degrees",
        caption="A read in the middle of a well-covered genome has several successors and predecessors.",
        footer=f"live run, seed {params.seed}",
    )
    out_degree = (inst.ov > 0).sum(axis=1)
    st.write(
        f"Reads without any successor: {int((out_degree == 0).sum())}; "
        f"median arrows out of a read: {int(np.median(out_degree))}."
    )
