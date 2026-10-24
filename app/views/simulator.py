"""Page 3, Simulator: choose read length, coverage, error rate and seed; compare with theory."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd  # noqa: E402
import streamlit as st  # noqa: E402
import ui  # noqa: E402

from genomeasm.appkit import figs, lab  # noqa: E402

st.title("Simulator: cutting a genome into reads")
ui.live_badge()
params = ui.instance_controls("sim")
st.markdown(
    "A **read** starts at a random place of the genome window, has the chosen length and carries random "
    "typing mistakes (substitutions). The program knows where every read really came from, so it can "
    "later grade the assembly."
)
issues = params.problems()
if issues:
    for issue in issues:
        st.warning(issue)
    st.stop()

reads = lab.simulate(params)
fig, stats = figs.depth_vs_theory(reads)
errors = int(reads.n_substitutions.sum())
m1, m2, m3, m4 = st.columns(4)
m1.metric("Reads", f"{reads.n_reads:,}")
m2.metric(
    "Coverage (realised)",
    f"{reads.realised_coverage:.2f}x",
    help="Letters covered by reads / genome window",
)
m3.metric(
    "Uncovered letters",
    f"{100 * stats['uncovered_observed']:.2f}%",
    help="Theory: Lander-Waterman exp(-c) = "
    f"{100 * stats['uncovered_lander_waterman']:.2f}%, exact model of this simulator = "
    f"{100 * stats['uncovered_exact']:.2f}%",
)
m4.metric("Read errors", f"{errors:,}", help="Number of substituted letters over all reads")

tab_depth, tab_map, tab_reads = st.tabs(
    ["Coverage against theory", "Map of read origins", "The reads"]
)
with tab_depth:
    ui.show_figure(
        fig,
        "coverage_depth_vs_theory",
        caption=(
            f"One read set of {reads.n_reads} reads, window of {params.length} letters, seed {params.seed}. "
            "Lander-Waterman (Poisson) ignores the ends of the window; the exact model includes them, "
            "so it matches the simulation better on short windows (Phase 2)."
        ),
        footer=f"live run, seed {params.seed}",
    )
with tab_map:
    ui.show_figure(
        figs.read_map(reads),
        "read_origin_map",
        caption="Each bar is one read at its true place. Red reads contain at least one error; grey: no read covers the position.",
        footer=f"live run, seed {params.seed}",
    )
with tab_reads:
    frame = pd.DataFrame(
        {
            "read": reads.ids,
            "true start": reads.starts,
            "length": reads.lengths,
            "errors": reads.n_substitutions + reads.n_insertions + reads.n_deletions,
            "sequence": reads.sequences,
        }
    )
    ui.table(frame.head(200))
    st.download_button(
        "Download these reads (FASTA)", reads.to_fasta(), "simulated_reads.fasta", "text/plain",
        icon=":material/download:", on_click="ignore",
    )  # fmt: skip
