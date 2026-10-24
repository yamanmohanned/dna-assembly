"""Charts of the application, drawn with matplotlib (offline, exportable as SVG and PNG).

Standards (master prompt): colour-blind-safe palette, labelled axes, sample sizes stated, log
scales stated, an axis that does not start at zero says so. Figures are created with
``matplotlib.figure.Figure`` (no pyplot state), so a long-running server does not leak them.
"""

from __future__ import annotations

import io
import math
from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from matplotlib.figure import Figure
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle

from genomeasm.coverage import (
    coverage_depth,
    depth_histogram,
    expected_depth_histogram,
    lw_depth_pmf,
    uncovered_fraction,
)
from genomeasm.instance import AssemblyInstance
from genomeasm.plots import OKABE_ITO as C
from genomeasm.plots import apply_style
from genomeasm.simulate import ReadSet
from genomeasm.solvers.alns import ARM_NAMES, OUTCOMES

PALETTE = [C["blue"], C["sky"], C["green"], C["orange"], C["vermillion"], C["purple"]]


def _fig(width: float = 8.0, height: float = 4.6, **kwargs: Any) -> tuple[Figure, Any]:
    apply_style()
    fig = Figure(figsize=(width, height), layout="constrained")
    return fig, fig.subplots(**kwargs)


def figure_bytes(fig: Figure, fmt: str, footer: str | None = None) -> bytes:
    """The figure as SVG (identical bytes on every system) or as a 300-dpi PNG."""
    if fmt not in ("svg", "png"):
        raise ValueError("fmt must be 'svg' or 'png'")
    note = fig.text(0.005, 0.004, footer, fontsize=7.5, color="#555555") if footer else None
    buffer = io.BytesIO()
    if fmt == "svg":
        fig.savefig(buffer, format="svg", metadata={"Date": None})
    else:
        fig.savefig(buffer, format="png", dpi=300)
    if note is not None:
        note.remove()
    data = buffer.getvalue()
    return data.replace(b"\r\n", b"\n") if fmt == "svg" else data


def live_frame(fig: Figure, dpi: int = 100) -> bytes:
    """A PNG of a live chart (cheap: the application shows it with st.image while a solver runs)."""
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=dpi, bbox_inches="tight")
    return buffer.getvalue()


# ---- Learn: the toy example ---------------------------------------------------------------------


def toy_graph(
    inst: AssemblyInstance,
    labels: Sequence[str],
    order: Sequence[int] | None,
    highlight: Sequence[tuple[int, int]] | None = None,
) -> Figure:
    """Reads as nodes, overlaps as arrows; the chosen order (or the given arrows) in red."""
    fig, ax = _fig(5.8, 5.2)
    ax.axis("off")
    ax.set_aspect("equal")
    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-1.5, 1.5)
    n = inst.n
    pos = [
        (math.cos(2 * math.pi * k / n + math.pi / 2), math.sin(2 * math.pi * k / n + math.pi / 2))
        for k in range(n)
    ]
    chosen = set(zip(order[:-1], order[1:], strict=True)) if order is not None else set()
    chosen |= set(highlight or ())
    rad = 0.22
    for i in range(n):
        for j in range(n):
            if i == j or inst.ov[i, j] <= 0:
                continue
            (x1, y1), (x2, y2) = pos[i], pos[j]
            hot = (i, j) in chosen
            ax.add_patch(
                FancyArrowPatch(
                    (x1, y1),
                    (x2, y2),
                    connectionstyle=f"arc3,rad={rad}",
                    arrowstyle="-|>",
                    mutation_scale=16,
                    shrinkA=22,
                    shrinkB=22,
                    color=C["vermillion"] if hot else "#8a8a8a",
                    lw=2.6 if hot else 1.2,
                    zorder=1,
                )
            )
            mx, my = (x1 + x2) / 2, (y1 + y2) / 2  # apex of the curved arrow
            ax.text(
                0.5 * mx + 0.5 * (mx + rad * (y2 - y1)),
                0.5 * my + 0.5 * (my - rad * (x2 - x1)),
                str(int(inst.ov[i, j])),
                fontsize=11,
                ha="center",
                va="center",
                zorder=3,
                color=C["vermillion"] if hot else "#333333",
                bbox={"boxstyle": "round,pad=0.12", "fc": "white", "ec": "none", "alpha": 0.9},
            )
    for k, (x, y) in enumerate(pos):
        ax.add_patch(Circle((x, y), 0.17, fc="#ffffff", ec=C["blue"], lw=2, zorder=4))
        ax.text(x, y, labels[k], ha="center", va="center", fontsize=15, weight="bold", zorder=5)
        ax.text(
            x * 1.33,
            y * 1.33,
            (inst.reads[k] if len(inst.reads[k]) <= 12 else inst.reads[k][:11] + "..."),
            ha="center",
            va="center",
            fontsize=9.5,
            family="monospace",
        )
    ax.set_title("Arrow A to B: letters saved by gluing B after A")
    return fig


def matrix_figure(
    matrix: np.ndarray, labels: Sequence[str], title: str, colour_label: str
) -> Figure:
    """A small matrix with its numbers written in the cells (rows: first read, columns: second)."""
    fig, ax = _fig(5.0, 4.4)
    ax.grid(False)
    masked = np.ma.masked_array(matrix, mask=np.eye(matrix.shape[0], dtype=bool))
    image = ax.imshow(masked, cmap="Blues", vmin=0)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if i != j:
                ax.text(j, i, str(int(matrix[i, j])), ha="center", va="center", fontsize=13)
    ax.set_xticks(range(len(labels)), list(labels))
    ax.set_yticks(range(len(labels)), list(labels))
    ax.set_xlabel("second read")
    ax.set_ylabel("first read")
    ax.set_title(title)
    fig.colorbar(image, ax=ax, label=colour_label)
    return fig


def tour_strip(
    order_before: Sequence[int], order_after: Sequence[int], labels: Sequence[str], title: str
) -> Figure:
    """Two rows of boxes: the order before and after one ALNS step (changed places outlined)."""
    apply_style()
    fig = Figure(figsize=(max(5.5, 0.9 * len(order_before) + 2.0), 3.0))
    ax = fig.add_axes((0.02, 0.03, 0.96, 0.8))
    ax.axis("off")
    n = len(order_before)
    ax.set_xlim(-1.6, n + 0.2)
    ax.set_ylim(-0.2, 2.6)
    for row, (name, order) in enumerate((("after", order_after), ("before", order_before))):
        y = 0.2 + 1.2 * row
        ax.text(-0.2, y + 0.4, name, ha="right", va="center", fontsize=11)
        for p, read in enumerate(order):
            changed = order_before[p] != order_after[p]
            ax.add_patch(
                Rectangle(
                    (p, y),
                    0.85,
                    0.8,
                    fc=PALETTE[read % 6],
                    alpha=0.85,
                    ec=C["vermillion"] if changed else "#ffffff",
                    lw=3 if changed else 1,
                )
            )
            ax.text(
                p + 0.425,
                y + 0.4,
                labels[read],
                ha="center",
                va="center",
                fontsize=13,
                weight="bold",
            )
    ax.set_title(title)
    return fig


# ---- Data ---------------------------------------------------------------------------------------


def gc_profile(sequence: str, window: int = 200) -> Figure:
    """GC content (share of G and C letters) in a sliding window along the genome."""
    fig, ax = _fig(8.0, 3.6)
    codes = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)
    gc = ((codes == ord("G")) | (codes == ord("C"))).astype(float)
    kernel = np.ones(window) / window
    profile = np.convolve(gc, kernel, mode="valid")
    ax.plot(np.arange(profile.size) + window // 2, 100 * profile, color=C["blue"], lw=1.4)
    ax.axhline(
        100 * gc.mean(),
        color=C["vermillion"],
        ls="--",
        lw=1.2,
        label=f"whole genome {100 * gc.mean():.1f}%",
    )
    ax.set_xlabel("Position on the genome (letters)")
    ax.set_ylabel(f"GC content in a {window}-letter window (%)")
    ax.set_title("GC content along the genome (axis does not start at 0)")
    ax.legend(loc="best")
    return fig


def base_counts(card: dict[str, Any]) -> Figure:
    fig, ax = _fig(5.0, 3.6)
    counts = card["base_counts"]
    letters = [b for b in "ACGT"]
    ax.bar(
        letters,
        [counts[b] for b in letters],
        color=[C["green"], C["blue"], C["orange"], C["vermillion"]],
    )
    ax.set_ylabel("Number of letters")
    ax.set_xlabel("Base")
    ax.set_title(f"Base counts ({card['length_bp']:,} letters)")
    return fig


PIPELINE_ONLINE = [
    ("configs/genomes.yaml", "list of versioned accessions"),
    ("validate_accession()", "format check, e.g. NC_001422.1"),
    ("fetch_from_ncbi()", "Entrez esummary (facts) + efetch (FASTA)"),
    ("parse_single_record(), Genome", "one record, only A, C, G, T"),
    ("cross-checks with NCBI", "version, molecule type, length"),
    ("save (atomic)", "data/raw/<acc>.fasta, .meta.json, SHA256SUMS"),
]
PIPELINE_OFFLINE = [
    ("load_genome()", "SHA-256 of the file equals SHA256SUMS?"),
    ("parse + Genome.from_record", "compare with <acc>.meta.json"),
    ("LoadedGenome", "genome + NCBI metadata + fingerprint"),
    ("build_card()", "data/processed/<acc>.card.json / .md"),
    ("simulate_reads()", "reads from a window of the genome"),
]


def pipeline_diagram() -> Figure:
    """The data-loading pipeline: the online part runs once, the offline part on every run."""
    fig, ax = _fig(10.5, 5.6)
    ax.axis("off")
    ax.set_xlim(0, 10.5)
    ax.set_ylim(0, 7.2)
    for x0, title, colour, steps in (
        (0.1, "once, with internet:  genomeasm data fetch", C["orange"], PIPELINE_ONLINE),
        (5.5, "every run, offline", C["green"], PIPELINE_OFFLINE),
    ):
        ax.text(x0 + 2.4, 6.95, title, ha="center", va="center", fontsize=12, weight="bold")
        for k, (head, detail) in enumerate(steps):
            y = 5.7 - 1.05 * k
            ax.add_patch(
                FancyBboxPatch(
                    (x0, y - 0.38),
                    4.8,
                    0.76,
                    boxstyle="round,pad=0.03",
                    fc=colour,
                    alpha=0.22,
                    ec=colour,
                    lw=1.6,
                )
            )
            ax.text(
                x0 + 0.15,
                y + 0.12,
                head,
                fontsize=10.5,
                weight="bold",
                va="center",
                family="monospace",
            )
            ax.text(x0 + 0.15, y - 0.17, detail, fontsize=9.5, va="center")
            if k + 1 < len(steps):
                ax.add_patch(
                    FancyArrowPatch(
                        (x0 + 2.4, y - 0.4),
                        (x0 + 2.4, y - 0.67),
                        arrowstyle="-|>",
                        mutation_scale=12,
                        color="#444444",
                    )
                )
    ax.add_patch(
        FancyArrowPatch(
            (4.95, 0.75),
            (5.5, 4.9),
            connectionstyle="arc3,rad=-0.25",
            arrowstyle="-|>",
            mutation_scale=14,
            color="#444444",
            ls="--",
        )
    )
    ax.text(5.25, 2.6, "files on disk", rotation=90, fontsize=9, va="center", ha="center")
    ax.set_title("If any check fails nothing is written; if NCBI is unreachable the program stops")
    return fig


# ---- Simulator ----------------------------------------------------------------------------------


def depth_vs_theory(reads: ReadSet) -> tuple[Figure, dict[str, float]]:
    """Observed coverage depth against Lander-Waterman and the exact finite-genome model."""
    g = reads.genome_length
    depth = coverage_depth(g, reads.starts, reads.spans)
    support, pmf = reads.config.read_length.pmf()
    c = reads.config.nominal_coverage(g)
    n_bins = int(max(8, min(60, np.max(depth) + 3, 3 * c + 6)))
    observed = depth_histogram(depth, n_bins)
    exact = expected_depth_histogram(g, reads.n_reads, support, pmf, n_bins)
    lw = lw_depth_pmf(c, n_bins)
    fig, ax = _fig(8.0, 4.4)
    x = np.arange(n_bins)
    ax.bar(x, observed, width=0.8, color=C["sky"], label=f"observed ({reads.n_reads} reads)")
    ax.step(
        x,
        lw,
        where="mid",
        color=C["vermillion"],
        lw=2,
        label=f"Lander-Waterman: Poisson(c = {c:.1f})",
    )
    ax.step(
        x,
        exact,
        where="mid",
        color=C["green"],
        lw=2,
        ls="--",
        label="exact model of this simulator",
    )
    ax.set_xlabel("Depth: number of reads covering a position")
    ax.set_ylabel("Share of genome positions")
    ax.set_title("Coverage depth against theory (one simulated read set)")
    ax.legend(loc="best")
    stats = {
        "reads": float(reads.n_reads),
        "nominal_coverage": float(c),
        "uncovered_observed": uncovered_fraction(depth),
        "uncovered_lander_waterman": float(lw[0]),
        "uncovered_exact": float(exact[0]),
        "mean_depth": float(depth.mean()),
    }
    return fig, stats


def read_map(reads: ReadSet, max_reads: int = 150) -> Figure:
    """Where every read comes from: bars sorted by start, and the depth along the genome."""
    apply_style()
    fig = Figure(figsize=(9.0, 5.6), layout="constrained")
    top, bottom = fig.subplots(2, 1, sharex=True, gridspec_kw={"height_ratios": [3, 1.2]})
    order = np.argsort(reads.starts, kind="stable")
    shown = order if order.size <= max_reads else order[:: math.ceil(order.size / max_reads)]
    for row, r in enumerate(shown):
        has_error = int(reads.n_substitutions[r] + reads.n_insertions[r] + reads.n_deletions[r]) > 0
        top.plot(
            [reads.starts[r], reads.starts[r] + reads.spans[r]],
            [row, row],
            lw=2.2,
            color=C["vermillion"] if has_error else C["blue"],
            solid_capstyle="butt",
        )
    top.plot([], [], color=C["blue"], lw=3, label="no error")
    top.plot([], [], color=C["vermillion"], lw=3, label="at least one error")
    top.set_ylabel(f"Read (sorted by start; {len(shown)} of {reads.n_reads} drawn)")
    top.legend(loc="lower right")
    top.set_title("Map of read origins")
    depth = coverage_depth(reads.genome_length, reads.starts, reads.spans)
    bottom.fill_between(np.arange(depth.size), depth, step="mid", color=C["sky"])
    bottom.fill_between(
        np.arange(depth.size),
        0,
        np.max(depth) + 1,
        where=depth == 0,
        color="#bbbbbb",
        step="mid",
        label="not covered",
    )
    bottom.set_xlabel("Position on the genome window (letters)")
    bottom.set_ylabel("Depth")
    bottom.set_xlim(0, reads.genome_length)
    if (depth == 0).any():
        bottom.legend(loc="upper right")
    return fig


# ---- Overlap graph ------------------------------------------------------------------------------


def overlap_heatmap(inst: AssemblyInstance, max_n: int = 150) -> Figure:
    """Overlap lengths between reads sorted by their true start (a band along the diagonal)."""
    n = inst.n
    order = (
        np.argsort(inst.true_starts, kind="stable")
        if inst.true_starts is not None
        else np.arange(n)
    )
    shown = order[:max_n]
    fig, ax = _fig(6.4, 5.6)
    ax.grid(False)
    image = ax.imshow(
        inst.ov[np.ix_(shown, shown)], cmap="viridis", origin="upper", interpolation="nearest"
    )
    ax.set_xlabel("second read (sorted by true position)")
    ax.set_ylabel("first read (sorted by true position)")
    suffix = f" (first {len(shown)} of {n} reads)" if n > max_n else f" ({n} reads)"
    ax.set_title("Overlap matrix" + suffix)
    fig.colorbar(image, ax=ax, label="overlap length (letters)")
    return fig


def overlap_lengths(inst: AssemblyInstance) -> Figure:
    values = inst.ov[inst.ov > 0]
    fig, ax = _fig(6.4, 4.0)
    if values.size:
        bins = np.arange(inst.min_overlap, int(np.max(values)) + 2) - 0.5
        ax.hist(values, bins=bins, color=C["blue"])
    ax.set_yscale("log")
    ax.set_xlabel("Overlap length (letters)")
    ax.set_ylabel("Number of read pairs (log scale)")
    ax.set_title(
        f"Overlap lengths of {values.size:,} read pairs (minimum overlap {inst.min_overlap})"
    )
    return fig


def degree_distribution(inst: AssemblyInstance) -> Figure:
    """How many successors and predecessors each read has (positive overlaps)."""
    out_degree = (inst.ov > 0).sum(axis=1)
    in_degree = (inst.ov > 0).sum(axis=0)
    fig, ax = _fig(6.4, 4.0)
    top = int(max(np.max(out_degree, initial=0), np.max(in_degree, initial=0)))
    bins = np.arange(0, top + 2) - 0.5
    ax.hist(
        [out_degree, in_degree],
        bins=bins,
        color=[C["blue"], C["orange"]],
        label=["out-degree (reads that can follow)", "in-degree (reads that can precede)"],
    )
    ax.set_xlabel("Number of reads with a positive overlap")
    ax.set_ylabel(f"Number of reads (n = {inst.n})")
    ax.set_title("Degree distribution of the overlap graph")
    ax.legend(loc="best")
    return fig


# ---- Solver lab ---------------------------------------------------------------------------------


def live_convergence(
    points: Sequence[tuple[float, float]], references: dict[str, float], seconds: float
) -> Figure:
    """Best fitness found so far against time; dashed lines are reference solutions."""
    fig, ax = _fig(7.0, 4.2)
    if points:
        xs, ys = zip(*points, strict=True)
        ax.step(xs, ys, where="post", color=C["blue"], lw=2.2, label="best F so far")
    for k, (name, value) in enumerate(references.items()):
        ax.axhline(
            value,
            color=[C["vermillion"], C["green"], C["orange"]][k % 3],
            ls="--",
            lw=1.4,
            label=name,
        )
    ax.set_xlim(0, max(seconds, 0.5))
    ax.set_xlabel("Wall-clock time (s)")
    ax.set_ylabel("Fitness F (letters saved; axis does not start at 0)")
    ax.set_title("Convergence")
    ax.legend(loc="lower right")
    return fig


def operator_probabilities(
    history: Sequence[tuple[float, Sequence[float]]], seconds: float
) -> Figure:
    """The chance of choosing each destroy + repair pair over time (stacked to 1)."""
    fig, ax = _fig(7.0, 4.2)
    if len(history) >= 2:
        t = np.array([h[0] for h in history])
        p = np.array([list(h[1]) for h in history]).T
        ax.stackplot(t, p, colors=PALETTE, labels=list(ARM_NAMES), alpha=0.9)
        ax.legend(
            fontsize=8.5, loc="center left", bbox_to_anchor=(1.0, 0.5), title="destroy + repair"
        )
    ax.set_xlim(0, max(seconds, 0.5))
    ax.set_ylim(0, 1)
    ax.set_xlabel("Wall-clock time (s)")
    ax.set_ylabel("Probability of choosing the pair")
    ax.set_title("Operator selection over time")
    return fig


def live_dashboard(
    points: Sequence[tuple[float, float]],
    references: dict[str, float],
    history: Sequence[tuple[float, Sequence[float]]],
    seconds: float,
) -> Figure:
    """Convergence and (for ALNS) operator probabilities in one figure: cheap to redraw live."""
    apply_style()
    fig = Figure(figsize=(10.0, 3.6), layout="constrained")
    left, right = fig.subplots(1, 2)
    if points:
        xs, ys = zip(*points, strict=True)
        left.step(xs, ys, where="post", color=C["blue"], lw=2.2, label="best F so far")
    for k, (name, value) in enumerate(references.items()):
        left.axhline(value, color=[C["vermillion"], C["green"]][k % 2], ls="--", lw=1.4, label=name)
    left.set_xlim(0, max(seconds, 0.5))
    left.set_xlabel("Wall-clock time (s)")
    left.set_ylabel("Fitness F (axis does not start at 0)")
    left.set_title("Convergence")
    left.legend(loc="lower right", fontsize=9)
    if len(history) >= 2:
        t = np.array([h[0] for h in history])
        p = np.array([list(h[1]) for h in history]).T
        right.stackplot(t, p, colors=PALETTE, labels=list(ARM_NAMES), alpha=0.9)
        right.legend(fontsize=8, loc="center left", bbox_to_anchor=(1.0, 0.5))
    right.set_xlim(0, max(seconds, 0.5))
    right.set_ylim(0, 1)
    right.set_xlabel("Wall-clock time (s)")
    right.set_ylabel("Probability of choosing the pair")
    right.set_title("Operator selection")
    return fig


def outcome_bars(outcomes: Sequence[Sequence[int]]) -> Figure:
    """What happened after each use of an operator pair (share of its uses)."""
    table = np.array(outcomes, dtype=float)
    shares = table / np.maximum(table.sum(axis=1, keepdims=True), 1)
    fig, ax = _fig(7.0, 4.0)
    left = np.zeros(len(ARM_NAMES))
    colours = [C["green"], C["sky"], C["orange"], "#999999", "#cccccc", C["vermillion"]]
    for k, name in enumerate(OUTCOMES):
        ax.barh(range(len(ARM_NAMES)), shares[:, k], left=left, color=colours[k], label=name)
        left += shares[:, k]
    ax.set_yticks(
        range(len(ARM_NAMES)),
        [f"{a} ({int(u)} uses)" for a, u in zip(ARM_NAMES, table.sum(axis=1), strict=True)],
    )
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("Share of the uses of the pair")
    ax.set_title("Outcome of every move, by operator pair")
    ax.legend(ncol=3, fontsize=8.5, loc="lower center", bbox_to_anchor=(0.5, -0.38))
    return fig


# ---- Assembly viewer ----------------------------------------------------------------------------


def dot_plot(dp: dict[str, Any]) -> Figure:
    """Assembled text (x) against the reference (y); a perfect assembly is one 45 degree line."""
    fig, ax = _fig(6.6, 5.6)
    x, y, contig = dp["x"], dp["y"], dp["contig"]
    mapped = ~np.isnan(y)
    for k in np.unique(contig[mapped]) if mapped.any() else []:
        sel = mapped & (contig == k)
        ax.plot(x[sel], y[sel], ".", ms=7, color=PALETTE[int(k) % 6])
    if (~mapped).any():
        ax.plot(
            x[~mapped],
            np.zeros((~mapped).sum()),
            "|",
            color="#777777",
            ms=10,
            label="window not located",
        )
        ax.legend(loc="lower right")
    for b in dp["boundaries"][1:]:
        ax.axvline(b, color="#aaaaaa", lw=0.8, ls=":")
    ax.set_xlabel("Position in the assembled text (contigs written one after the other, letters)")
    ax.set_ylabel("Position of the window on the real genome (letters)")
    ax.set_title(
        f"Dot plot: assembly against reference ({len(dp['boundaries'])} contigs, dotted = ends)"
    )
    return fig


def tour_blocks(
    order: Sequence[int], true_starts: np.ndarray | None, max_reads: int = 200
) -> Figure:
    """The order of the reads: true start of each read at its place in the order."""
    fig, ax = _fig(7.0, 3.6)
    if true_starts is not None:
        ys = true_starts[np.asarray(order)][:max_reads]
        ax.plot(np.arange(ys.size), ys, ".", color=C["blue"], ms=5)
        ax.set_ylabel("True start of the read (letters)")
    ax.set_xlabel("Place of the read in the solver's order")
    ax.set_title(
        "Order of the reads against their true positions (a rising line is a correct order)"
    )
    return fig


# ---- Experiments --------------------------------------------------------------------------------


def box_by_algorithm(
    data: dict[str, np.ndarray], ylabel: str, title: str, log: bool = False
) -> Figure:
    """One horizontal box per algorithm (quartiles, 1.5 x IQR whiskers); n is in the label."""
    names = list(data)
    fig, ax = _fig(8.5, max(3.2, 0.62 * len(names) + 1.6))
    box = ax.boxplot(
        [data[a] for a in names],
        orientation="horizontal",
        patch_artist=True,
        widths=0.6,
        flierprops={"markersize": 3, "alpha": 0.4},
    )
    for patch, k in zip(box["boxes"], range(len(names)), strict=True):
        patch.set_facecolor(PALETTE[k % 6])
        patch.set_alpha(0.65)
    ax.set_yticks(range(1, len(names) + 1), [f"{a}  (n = {len(data[a])})" for a in names])
    ax.invert_yaxis()
    if log:
        ax.set_xscale("symlog", linthresh=0.01)
        ylabel += " (symlog scale)"
    ax.set_xlabel(ylabel)
    ax.set_title(title)
    return fig


def runs_frame_to_series(
    frame: pd.DataFrame, algorithms: Sequence[str], column: str
) -> dict[str, np.ndarray]:
    return {a: frame[frame["algorithm"] == a][column].to_numpy(dtype=float) for a in algorithms}


# ---- DNA and RNA pictures -----------------------------------------------------------------------

BASE_COLOUR = {
    "A": C["green"],
    "T": C["vermillion"],
    "G": C["orange"],
    "C": C["blue"],
    "U": C["purple"],
}
COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}


def dna_helix(sequence: str, *, rna: bool = False, title: str | None = None) -> Figure:
    """A short stretch of real sequence as a double helix (DNA) or as a single strand (RNA).

    DNA: the second strand is the complement (A with T, G with C), ten and a half letters per turn.
    RNA: one strand in which every T is written U.
    """
    seq = "".join(c for c in sequence.upper() if c in "ACGT")[:24]
    n = len(seq)
    apply_style()
    fig = Figure(figsize=(max(6.0, 0.42 * n + 1.5), 2.6), layout="constrained")
    ax = fig.subplots()
    ax.axis("off")
    ax.set_xlim(-1.0, n)
    if rna:
        ax.set_ylim(-1.9, 1.3)
        xs = np.linspace(-0.5, n - 0.5, 30 * n)
        ys = 0.3 * np.sin(2 * np.pi * xs / 4.0)
        ax.plot(xs, ys, color="#5b6b7b", lw=5, solid_capstyle="round", zorder=1)
        for i, base in enumerate(seq):
            letter = "U" if base == "T" else base
            y0 = 0.3 * np.sin(2 * np.pi * i / 4.0)
            ax.plot([i, i], [y0, y0 - 0.9], color="#5b6b7b", lw=2, zorder=2)
            ax.scatter(i, y0 - 1.0, s=330, color=BASE_COLOUR[letter], zorder=3, edgecolor="white")
            ax.text(
                i,
                y0 - 1.0,
                letter,
                ha="center",
                va="center",
                color="white",
                weight="bold",
                zorder=4,
            )
        ax.set_title(title or "RNA: a single strand, U instead of T", fontsize=11)
        return fig
    ax.set_ylim(-1.5, 1.5)
    xs = np.linspace(-0.5, n - 0.5, 30 * n)
    phase = 2 * np.pi * xs / 10.5
    for sign, colour in ((1, "#2f4b66"), (-1, "#8a5a1f")):
        y = sign * np.sin(phase)
        front = sign * np.cos(phase)  # > 0: this strand passes in front of the other one
        points = np.column_stack([xs, y]).reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        ax.add_collection(
            LineCollection(
                list(segments), colors=colour, linewidths=2.0 + 2.6 * (front[:-1] > 0), zorder=2
            )
        )
    for i, base in enumerate(seq):
        p = 2 * np.pi * i / 10.5
        y_top, other = np.sin(p), COMPLEMENT[base]
        z = 5 if np.cos(p) > 0 else 3
        ax.plot([i, i], [y_top, 0], color=BASE_COLOUR[base], lw=4, zorder=z - 1, alpha=0.85)
        ax.plot([i, i], [0, -y_top], color=BASE_COLOUR[other], lw=4, zorder=z - 1, alpha=0.85)
        for y, letter in ((y_top, base), (-y_top, other)):
            ax.scatter(i, y, s=300, color=BASE_COLOUR[letter], zorder=z, edgecolor="white")
            ax.text(
                i, y, letter, ha="center", va="center", color="white", weight="bold", zorder=z + 1
            )
    ax.set_title(title or "DNA: two strands, A pairs with T and G with C", fontsize=11)
    return fig
