"""Shared widgets of the application: exportable charts, live/precomputed badges, instance controls.

Rules of the application (master prompt, Phase 6): it runs offline, every chart has SVG and PNG
export buttons, any live run ends in under 30 seconds, and larger runs are shown from precomputed
results with a label that says so.
"""

from __future__ import annotations

import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

import streamlit as st
from matplotlib.figure import Figure

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:  # run from a checkout without installing the package
    sys.path.insert(0, str(ROOT))

from genomeasm.appkit import figs, lab  # noqa: E402
from genomeasm.appkit.paths import GENOMES  # noqa: E402
from genomeasm.overlap import recommend_min_overlap  # noqa: E402


@st.cache_resource(show_spinner="Compiling the solver kernels (only the first time)...")
def _warm_up() -> bool:
    from genomeasm.solvers.warmup import warmup

    warmup()
    return True


def ensure_ready() -> None:
    """Compile the numba kernels once per server process, before any timed run."""
    _warm_up()


def show_figure(
    fig: Figure, name: str, caption: str | None = None, footer: str | None = None
) -> None:
    """Draw a chart with SVG and PNG download buttons (the files are made only on click)."""
    st.pyplot(fig, width="stretch", alt=caption or name.replace("_", " "))
    c1, c2, _ = st.columns([1, 1, 5])
    for column, fmt, mime in ((c1, "svg", "image/svg+xml"), (c2, "png", "image/png")):
        column.download_button(
            fmt.upper(),
            data=lambda f=fmt: figs.figure_bytes(fig, f, footer),  # type: ignore[misc]
            file_name=f"{name}.{fmt}",
            mime=mime,
            key=f"dl-{name}-{fmt}",
            icon=":material/download:",
            on_click="ignore",
            help=f"Export this chart as {fmt.upper()}",
        )
    if caption:
        st.caption(caption)


@st.cache_data(show_spinner=False)
def _helix_png(letters: int, rna: bool) -> bytes:
    sequence = lab.genome_sequence("NC_001422.1")[:letters]
    return figs.live_frame(figs.dna_helix(sequence, rna=rna, title=" "), dpi=110)


def banner() -> None:
    """A strip of real phiX174 DNA (the first 24 letters) as a double helix."""
    st.image(
        _helix_png(24, False), width="stretch",
        alt="Double helix drawn from the first 24 letters of the phiX174 genome",
    )  # fmt: skip


FORCE_LTR = """
<script>
(function () {
  // The sliders of Streamlit follow the language of the browser: with an Arabic browser the thumb runs
  // from right to left while the coloured bar runs from left to right. The text of this application is
  // English, so left-to-right is forced and the page is marked as not to be machine-translated.
  var root = document.documentElement;
  root.setAttribute("lang", "en");
  root.setAttribute("dir", "ltr");
  root.setAttribute("translate", "no");
  if (!document.querySelector('meta[name="google"]')) {
    var meta = document.createElement("meta");
    meta.name = "google";
    meta.content = "notranslate";
    document.head.appendChild(meta);
  }
  if (!String(navigator.language).toLowerCase().startsWith("en")) {
    Object.defineProperty(navigator, "language", { get: function () { return "en-US"; }, configurable: true });
    Object.defineProperty(navigator, "languages", { get: function () { return ["en-US"]; }, configurable: true });
    window.dispatchEvent(new Event("languagechange"));
  }
  window.__ltrForced = true;
})();
</script>
"""


def force_left_to_right() -> None:
    """Keep the sliders consistent on browsers set to a right-to-left language (D-054)."""
    st.html(FORCE_LTR, unsafe_allow_javascript=True)


def logo() -> None:
    st.logo(_helix_png(12, False), size="large")


def live_badge(text: str = "Live run on this computer") -> None:
    st.badge(text, icon=":material/bolt:", color="green")


def precomputed_badge(text: str = "Precomputed result (Phase 5)") -> None:
    st.badge(text, icon=":material/inventory_2:", color="orange")


def table(data: Any, **kwargs: Any) -> None:
    st.dataframe(data, hide_index=True, width="stretch", **kwargs)


def genome_length(accession: str) -> int:
    return len(lab.genome_sequence(accession))


def instance_controls(prefix: str, *, overlap_options: bool = False) -> lab.InstanceParams:
    """Sidebar controls of a read set. The values are remembered and offered on the other pages."""
    saved: lab.InstanceParams = st.session_state.get("instance", lab.InstanceParams())
    with st.sidebar:
        st.subheader("Read set")
        accessions = list(GENOMES)
        genome = st.selectbox(
            "Genome",
            accessions,
            index=accessions.index(saved.genome),
            format_func=lambda a: f"{a}: {GENOMES[a]}",
            key=f"{prefix}-genome",
        )
        total = genome_length(genome)
        max_window = min(lab.MAX_WINDOW, total)
        length = st.slider(
            "Window length (letters)", 500, max_window, min(saved.length, max_window), 100,
            key=f"{prefix}-length", format="%d",
            help="Reads are cut from this window of the real genome",
        )  # fmt: skip
        last_start = (total - length) // 100 * 100  # the largest start the steps of 100 can reach
        if st.session_state.get(f"{prefix}-start", 0) > last_start:
            st.session_state[f"{prefix}-start"] = (
                last_start  # the window got longer: keep it inside
            )
        start = st.slider(
            "Window start (letters)", 0, max(last_start, 100), min(saved.start, last_start), 100,
            key=f"{prefix}-start", format="%d",
        )  # fmt: skip
        read_length = st.slider(
            "Read length (letters)", 40, 200, saved.read_length, 10, key=f"{prefix}-rl", format="%d"
        )
        coverage = st.slider(
            "Coverage (x)", 1.0, 20.0, float(saved.coverage), 0.5, key=f"{prefix}-cov", format="%.1f",
            help="Average number of reads covering a letter of the genome",
        )  # fmt: skip
        error = (
            st.slider(
                "Substitution error rate (%)",
                0.0,
                3.0,
                round(100 * saved.error, 1),
                0.1,
                key=f"{prefix}-err",
                format="%.1f",
            )
            / 100
        )
        error = round(error, 4)
        seed = int(st.number_input("Seed", 0, 1_000_000, saved.seed, 1, key=f"{prefix}-seed"))
        tolerance, min_overlap = saved.tolerance, saved.min_overlap
        if overlap_options:
            tolerance = (
                st.slider(
                    "Mismatches allowed in an overlap (%)",
                    0.0,
                    10.0,
                    round(100 * saved.tolerance, 1),
                    0.5,
                    key=f"{prefix}-tol",
                    format="%.1f",
                    help="0 means exact overlaps; 5 lets one letter in 20 of an overlap differ",
                )  # fmt: skip
                / 100
            )
            min_overlap = st.slider(
                "Minimum overlap tau (letters)", 5, max(6, read_length - 1),
                min(max(saved.min_overlap, 5), max(read_length - 1, 6)), 1, key=f"{prefix}-tau", format="%d",
            )  # fmt: skip
        n_est = round(coverage * length / read_length)
        st.metric(
            "Reads in this set", f"about {n_est}", help=f"Live runs allow at most {lab.MAX_READS}"
        )
        st.caption(f"Statistical rule for tau: {recommend_min_overlap(max(n_est, 2))} letters")
    params = lab.InstanceParams(
        genome=genome, start=start, length=length, coverage=coverage, read_length=read_length,
        error=error, seed=seed, tolerance=tolerance, min_overlap=min_overlap,
    )  # fmt: skip
    st.session_state["instance"] = params
    return params


def build_instance(params: lab.InstanceParams) -> Any | None:
    """Build (or fetch from the cache) the instance; explain problems instead of crashing."""
    issues = params.problems()
    if issues:
        for issue in issues:
            st.warning(issue)
        return None
    with st.spinner("Simulating reads and computing the overlap graph..."):
        return lab.build(params)


ProgressFn = Callable[[dict[str, Any]], None]
