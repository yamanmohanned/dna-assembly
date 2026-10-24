"""Figure style and the Phase 2 figures.

Chart standards (master prompt): vector export (SVG and PDF) plus a 300-dpi PNG, one colour-blind
safe palette (Okabe-Ito), large readable fonts, labelled axes with units, log scales stated,
uncertainty shown, sample sizes stated, and a caption with data provenance (config hash and
results file) written next to every figure and as a small footer inside it.
"""

from __future__ import annotations

import io
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from cycler import cycler  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

OKABE_ITO = {
    "blue": "#0072B2",
    "orange": "#E69F00",
    "green": "#009E73",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
    "sky": "#56B4E9",
    "yellow": "#F0E442",
    "black": "#000000",
}

CYCLE = ["blue", "vermillion", "green", "orange", "purple", "sky", "black"]  # default colours


def apply_style() -> None:
    """Global matplotlib style: readable at projector size, deterministic file output."""
    plt.rcParams.update(
        {
            "axes.prop_cycle": cycler(color=[OKABE_ITO[k] for k in CYCLE]),
            "font.family": "DejaVu Sans",
            "font.size": 12,
            "axes.titlesize": 12,
            "axes.labelsize": 12,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "legend.fontsize": 10.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.25,
            "grid.linewidth": 0.6,
            "figure.dpi": 100,
            "savefig.dpi": 300,
            "svg.hashsalt": "genomeasm",
            "svg.fonttype": "path",
            "pdf.fonttype": 42,
        }
    )


def save_figure(
    fig: Figure, stem: str, out_dir: Path, *, caption: str, source: str, config_hash: str
) -> list[Path]:
    """Save SVG, PDF and 300-dpi PNG plus a ``<stem>.caption.json`` with the provenance."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig.text(
        0.005,
        0.004,
        f"Source: {source}  |  config sha256 {config_hash[:12]}",
        fontsize=7.5,
        color="#555555",
        ha="left",
        va="bottom",
    )
    paths = [out / f"{stem}.svg", out / f"{stem}.pdf", out / f"{stem}.png"]
    # Matplotlib writes SVG as text, which on Windows turns "\n" into "\r\n". Normalise so the
    # file has identical bytes on every operating system.
    svg = io.BytesIO()
    fig.savefig(svg, format="svg", metadata={"Date": None})
    paths[0].write_bytes(svg.getvalue().replace(b"\r\n", b"\n"))
    fig.savefig(paths[1], format="pdf", metadata={"CreationDate": None, "Creator": "genomeasm"})
    fig.savefig(paths[2], format="png", dpi=300)
    sidecar = out / f"{stem}.caption.json"
    sidecar.write_text(
        json.dumps(
            {
                "figure": stem,
                "caption": caption,
                "source": source,
                "config_hash": config_hash,
                "files": [p.name for p in paths],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
        newline="\n",
    )
    plt.close(fig)
    return [*paths, sidecar]


def _sci(value: float) -> str:
    return f"{value:.2e}" if value < 0.01 else f"{value:.3f}"


def _labels(results: dict[str, Any]) -> dict[str, str]:
    return dict(results["config"].get("genome_labels", {}))


def figure_coverage_histogram(results: dict[str, Any]) -> tuple[Figure, str]:
    """Figure 1: simulated read depth versus Lander-Waterman (Poisson) and the exact model."""
    apply_style()
    cfg = results["config"]
    accessions = list(cfg["genomes"])
    cs = [float(c) for c in cfg["histogram_coverages"]]
    labels = _labels(results)
    fig, axes = plt.subplots(
        len(accessions),
        len(cs),
        figsize=(4.4 * len(cs), 3.5 * len(accessions) + 0.9),
        squeeze=False,
    )
    by_key = {(e["accession"], e["coverage_requested"]): e for e in results["coverage"]}
    worst = 0.0
    for row, accession in enumerate(accessions):
        for col, c in enumerate(cs):
            ax = axes[row][col]
            entry = by_key[(accession, c)]
            h = entry["histogram"]
            x = np.arange(h["bins"])
            mean = np.array(h["observed_mean"])
            lo = np.array(h["observed_q25"])
            hi = np.array(h["observed_q75"])
            c_nom = entry["nominal_coverage"]
            xmax = min(h["bins"] - 1, int(c_nom + 5 * np.sqrt(c_nom) + 4))
            ax.bar(
                x,
                mean,
                width=0.8,
                color=OKABE_ITO["sky"],
                alpha=0.75,
                label="Simulation (mean of runs)",
                yerr=[np.clip(mean - lo, 0, None), np.clip(hi - mean, 0, None)],
                error_kw={"elinewidth": 0.8, "capsize": 1.5, "ecolor": "#333333"},
            )
            ax.plot(
                x,
                h["poisson"],
                "o--",
                color=OKABE_ITO["vermillion"],
                ms=3.2,
                lw=1.2,
                label="Lander-Waterman: Poisson(c)",
            )
            ax.plot(
                x,
                h["exact"],
                "-",
                color=OKABE_ITO["black"],
                lw=1.6,
                label="Exact finite-genome model",
            )
            ax.set_xlim(-0.8, xmax + 0.8)
            ax.set_title(
                f"{labels.get(accession, accession)}, c = {c_nom:.1f} ({entry['n_reads']:,} reads)"
            )
            ax.set_xlabel("Reads covering a position (depth)")
            if col == 0:
                ax.set_ylabel("Fraction of genome positions")
            u = entry["uncovered"]
            ax.text(
                0.97,
                0.95,
                f"uncovered fraction\nsimulation {_sci(u['mean'])}\nexact {_sci(u['exact'])}\n"
                f"Lander-Waterman {_sci(u['lander_waterman'])}",
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=8.5,
                bbox={
                    "boxstyle": "round",
                    "facecolor": "white",
                    "alpha": 0.85,
                    "edgecolor": "#BBBBBB",
                },
            )
            worst = max(worst, float(np.max(np.abs(mean - np.array(h["exact"])))))
    handles, names = axes[0][0].get_legend_handles_labels()
    fig.legend(
        handles, names, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.995)
    )
    reps = results["coverage"][0]["replicates"]
    length = cfg["read_length"]
    fig.tight_layout(rect=(0, 0.025, 1, 0.93))
    caption = (
        f"Read depth along the genome for simulated single-strand reads of length {length}. Bars: "
        f"mean over {reps} simulated runs (whiskers: interquartile range across runs). Black line: "
        "exact distribution for this simulator (finite genome, reads fully inside it). Dashed red: "
        "Lander-Waterman Poisson approximation, which ignores the genome ends. Largest absolute "
        f"gap between simulation and the exact model in any bin: {worst:.1e}."
    )
    return fig, caption


def figure_validation_curves(results: dict[str, Any]) -> tuple[Figure, str]:
    """Figure S1: uncovered fraction, number of contigs and the Lander-Waterman check."""
    apply_style()
    cfg = results["config"]
    labels = _labels(results)
    colours = [OKABE_ITO["blue"], OKABE_ITO["green"], OKABE_ITO["purple"]]
    fig, (ax_u, ax_k, ax_l) = plt.subplots(1, 3, figsize=(14.4, 4.4))
    for index, accession in enumerate(cfg["genomes"]):
        entries = sorted(
            (e for e in results["coverage"] if e["accession"] == accession),
            key=lambda e: e["nominal_coverage"],
        )
        c = np.array([e["nominal_coverage"] for e in entries])
        colour = colours[index % len(colours)]
        name = labels.get(accession, accession)
        g = entries[0]["genome_length"]
        u = [e["uncovered"] for e in entries]
        ax_u.errorbar(
            c,
            [v["mean"] for v in u],
            yerr=[1.96 * v["se"] for v in u],
            fmt="o",
            color=colour,
            ms=5,
            capsize=2.5,
            label=f"{name} ({g:,} bp): simulation, 95% CI",
        )
        ax_u.plot(
            c, [v["exact"] for v in u], "-", color=colour, lw=1.6, label=f"{name}: exact model"
        )
        k = [e["contigs"] for e in entries]
        ax_k.errorbar(
            c,
            [v["mean"] for v in k],
            yerr=[1.96 * v["se"] for v in k],
            fmt="o",
            color=colour,
            ms=5,
            capsize=2.5,
            label=f"{name}: simulation, 95% CI",
        )
        ax_k.plot(
            c, [v["exact"] for v in k], "-", color=colour, lw=1.6, label=f"{name}: exact formula"
        )
        ax_k.plot(
            c,
            [v["lander_waterman"] for v in k],
            "--",
            color=colour,
            lw=1.2,
            label=f"{name}: Lander-Waterman",
        )
    grid = np.linspace(0.8, max(float(c) for c in cfg["coverages"]) + 0.5, 100)
    ax_u.plot(
        grid,
        np.exp(-grid),
        "--",
        color=OKABE_ITO["vermillion"],
        lw=1.4,
        label="Lander-Waterman exp(-c)",
    )
    ax_u.set_yscale("log")
    ax_u.set_xlabel("Nominal coverage c = N l / G (x)")
    ax_u.set_ylabel("Uncovered fraction of the genome (log scale)")
    ax_u.set_title("(a) Uncovered positions")
    ax_u.legend(fontsize=8.5, loc="lower left")
    ax_k.set_yscale("log")
    ax_k.set_xlabel("Nominal coverage c = N l / G (x)")
    ax_k.set_ylabel("Number of true contigs (log scale)")
    ax_k.set_title(f"(b) Contigs (overlap >= {cfg['min_overlap']} bases)")
    ax_k.legend(fontsize=8.5, loc="upper right")

    large = results["large_genome"]
    cl = np.array([e["nominal_coverage"] for e in large])
    gap = np.array([e["relative_gap"] for e in large]) * 100
    err = np.array([1.96 * e["uncovered_se"] / e["lander_waterman"] for e in large]) * 100
    tol = cfg["tolerances"]["lw_large_genome_rel"] * 100
    ax_l.axhspan(-tol, tol, color=OKABE_ITO["sky"], alpha=0.25, label=f"allowed band +-{tol:.0f}%")
    ax_l.axhline(0, color="#444444", lw=0.8)
    ax_l.errorbar(
        cl,
        gap,
        yerr=err,
        fmt="o",
        color=OKABE_ITO["blue"],
        ms=5,
        capsize=2.5,
        label="simulation, 95% CI",
    )
    ax_l.set_xlabel("Nominal coverage c = N l / G (x)")
    ax_l.set_ylabel("(simulated / exp(-c) - 1)  (%)")
    ax_l.set_title(f"(c) Lander-Waterman check on a\n{large[0]['genome_length']:,}-position genome")
    ax_l.legend(fontsize=8.5, loc="best")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    reps = results["coverage"][0]["replicates"]
    caption = (
        f"Simulator versus theory. (a, b) Mean over {reps} runs with 95% confidence intervals for "
        "two real genomes; solid lines: exact formulas derived for this simulator; dashed: "
        "Lander-Waterman. On short genomes the genome ends leave more positions uncovered than "
        "exp(-c) predicts, and the exact model explains this. (c) On a long genome, where the "
        f"approximation applies, exp(-c) holds ({large[0]['replicates']} runs per point)."
    )
    return fig, caption


def figure_overlap_length_distribution(results: dict[str, Any]) -> tuple[Figure, str]:
    """Figure 2: lengths of the overlaps found, split into true neighbours and other pairs."""
    apply_style()
    labels = _labels(results)
    panels = results["figure_data"]["overlap_length"]
    fig, axes = plt.subplots(1, len(panels), figsize=(6.2 * len(panels), 4.4), squeeze=False)
    summary: list[str] = []
    for ax, panel in zip(axes[0], panels, strict=True):
        runs = panel["runs"]
        true = np.array(panel["hist_true"], dtype=float) / runs
        other = np.array(panel["hist_other"], dtype=float) / runs
        k = np.arange(true.size)
        ax.bar(k, np.where(true > 0, true, np.nan), width=0.9, color=OKABE_ITO["blue"],
               label="true neighbours (positions overlap)")  # fmt: skip
        ax.bar(k, np.where(other > 0, other, np.nan), width=0.9, color=OKABE_ITO["vermillion"],
               label="other pairs (chance or repeat)")  # fmt: skip
        ax.axvline(
            panel["tau"],
            color="#444444",
            lw=1,
            ls="--",
            label=f"minimum overlap tau = {panel['tau']}",
        )
        ax.set_yscale("log")
        ax.set_ylim(top=float(max(true.max(), other.max())) * 25)  # headroom for the legend
        ax.set_xlabel("Overlap length k (bases)")
        ax.set_ylabel("Number of edges per run (log scale)")
        name = labels.get(panel["accession"], panel["accession"])
        ax.set_title(f"{name}, c = {panel['coverage']:g}, mean of {runs} runs")
        ax.legend(fontsize=9, loc="upper right")
        summary.append(f"{name}: {true.sum():.0f} true and {other.sum():.1f} other edges per run")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    caption = (
        "Overlap lengths between simulated error-free reads of length 100. Blue: pairs whose real "
        "positions overlap by that many bases. Red: all other edges, which are chance matches or "
        "repeats. Dashed line: the minimum overlap tau chosen by the random-sequence bound. "
        + "; ".join(summary)
        + "."
    )
    return fig, caption


def figure_overlap_heatmap(results: dict[str, Any]) -> tuple[Figure, str]:
    """Figure 3: overlap matrix of a small set of reads sorted by their true position."""
    apply_style()
    data = results["figure_data"]["heatmap"]
    labels = _labels(results)
    matrix = np.array(data["matrix"], dtype=float)
    shown = np.where(matrix > 0, matrix, np.nan)
    fig, ax = plt.subplots(figsize=(6.4, 5.6))
    image = ax.imshow(shown, cmap="viridis", origin="upper", interpolation="nearest")
    colorbar = fig.colorbar(image, ax=ax)
    colorbar.set_label("Overlap length ov(i, j) (bases)")
    ax.set_xlabel("Second read j (reads sorted by true position)")
    ax.set_ylabel("First read i (reads sorted by true position)")
    ax.grid(False)
    name = labels.get(data["accession"], data["accession"])
    ax.set_title(f"{name}: overlaps of {data['n_shown']} reads, c = {data['coverage']:g}")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    caption = (
        f"Overlap matrix of the {data['n_shown']} simulated error-free reads that start in the "
        f"first {results['config']['figure_heatmap']['max_start']} bases of {name}, sorted by "
        "their true start. Entry (i, j) is the longest suffix of read i that equals a prefix of "
        f"read j (minimum overlap {data['tau']}); white means no overlap. The band lies above the "
        "diagonal only: overlaps are directional, ov(i, j) is not ov(j, i)."
    )
    return fig, caption


def make_phase3_figures(results: dict[str, Any], out_dir: Path, source: str) -> list[Path]:
    """Create and save the Phase 3 figures from an overlap validation results dictionary."""
    config_hash = results["provenance"]["config_hash"]
    paths: list[Path] = []
    for stem, builder in (
        ("fig02_overlap_length_distribution", figure_overlap_length_distribution),
        ("fig03_overlap_heatmap", figure_overlap_heatmap),
    ):
        fig, caption = builder(results)
        paths += save_figure(
            fig, stem, out_dir, caption=caption, source=source, config_hash=config_hash
        )
    return paths


def make_phase2_figures(results: dict[str, Any], out_dir: Path, source: str) -> list[Path]:
    """Create and save the Phase 2 figures from a validation results dictionary."""
    config_hash = results["provenance"]["config_hash"]
    paths: list[Path] = []
    for stem, builder in (
        ("fig01a_coverage_histogram", figure_coverage_histogram),
        ("figS1_validation_curves", figure_validation_curves),
    ):
        fig, caption = builder(results)
        paths += save_figure(
            fig, stem, out_dir, caption=caption, source=source, config_hash=config_hash
        )
    return paths
