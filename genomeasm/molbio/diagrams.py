"""Scientific diagrams as SVG text, drawn from data (no pictures are stored in the project).

Used by the molecular-biology pages and by the lecture, so one drawing serves both. Everything
quantitative in a diagram comes from its arguments: the gene map from the NCBI annotation, the
read map from the simulator, the gel from a real sequence, the bars from the genomes and the
facts file. The purely explanatory drawings (nucleotide, packaging, central dogma) are schematic
and are labelled as such where they are shown.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from html import escape
from typing import Any

from genomeasm.molbio.library import Gene
from genomeasm.molbio.sequence import AMINO_ACIDS, CODON_TABLE

INK = "#10243b"
SOFT = "#4a5b70"
FAINT = "#8a97a8"
LINE = "#d9e1ea"
PANEL = "#f3f6f9"
ACCENT = "#0072b2"
SIGNAL = "#d55e00"
BASE_COLOUR = {"A": "#00a878", "T": "#e0602a", "G": "#e8a317", "C": "#2f86d6", "U": "#c46aa0"}
CLASS_COLOUR = {
    "nonpolar": "#eef1f5",
    "polar": "#dceefb",
    "acidic": "#fde3d5",
    "basic": "#d6f3ea",
    "stop": "#f7d4d4",
}
GENE_COLOURS = (
    "#0072b2", "#56b4e9", "#009e73", "#e69f00", "#d55e00", "#cc79a7",
    "#5b6f8a", "#8c6bb1", "#2a9d8f", "#b5651d", "#3d5a80",
)  # fmt: skip
FONT = "'Segoe UI', system-ui, -apple-system, Arial, sans-serif"
SMALL_TEXT_SCALE = 1.14  # annotations below 14 units are enlarged by this factor


def _svg(width: int, height: int, body: Sequence[str], label: str) -> str:
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
        f'role="img" aria-label="{escape(label)}" font-family="{FONT}" '
        f'style="width:100%;height:auto;max-height:100%;display:block">' + "".join(body) + "</svg>"
    )


def _text(x: float, y: float, value: str, *, size: float = 14, fill: str = INK, weight: int = 400,
          anchor: str = "middle", style: str = "") -> str:  # fmt: skip
    extra = f' font-style="{style}"' if style else ""
    size = round(size * SMALL_TEXT_SCALE, 1) if size < 14 else size  # legible on a projector
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" fill="{fill}" font-weight="{weight}" '
        f'text-anchor="{anchor}"{extra}>{escape(value)}</text>'
    )


def _line(x1: float, y1: float, x2: float, y2: float, stroke: str = INK, width: float = 2,
          dash: str = "") -> str:  # fmt: skip
    extra = f' stroke-dasharray="{dash}"' if dash else ""
    return (
        f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{stroke}" '
        f'stroke-width="{width}" stroke-linecap="round"{extra}/>'
    )


def _arrow_defs() -> str:
    def marker(name: str, colour: str) -> str:
        return (
            f'<marker id="{name}" viewBox="0 0 10 10" refX="8.5" refY="5" markerWidth="7" '
            f'markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" '
            f'fill="{colour}"/></marker>'
        )

    return (
        "<defs>"
        + marker("ah", INK) + marker("ah-accent", ACCENT) + marker("ah-soft", FAINT)
        + marker("ah-signal", SIGNAL)
        + "</defs>"
    )  # fmt: skip


# ---- the nucleotide ------------------------------------------------------------------------------


def nucleotide_svg() -> str:
    """One deoxyribonucleotide with its three parts and the numbering of the sugar carbons."""
    b = [_arrow_defs()]
    ring = {
        "O4": (360, 150),
        "C1": (428, 198),
        "C2": (402, 276),
        "C3": (318, 276),
        "C4": (292, 198),
    }
    points = " ".join(f"{x},{y}" for x, y in ring.values())
    b.append(f'<polygon points="{points}" fill="#ece8f8" stroke="#7a68c0" stroke-width="3"/>')
    b.append(_text(360, 232, "deoxyribose", size=15, fill="#5a4aa8", weight=700))
    b.append(_text(360, 251, "(sugar)", size=13, fill="#5a4aa8"))
    b.append('<circle cx="360" cy="150" r="15" fill="#fff" stroke="#7a68c0" stroke-width="2"/>')
    b.append(_text(360, 155, "O", size=14, fill="#5a4aa8", weight=700))
    for label, (x, y), dx, dy in (
        ("1′", ring["C1"], 20, -6), ("2′", ring["C2"], 22, 12),
        ("3′", ring["C3"], -22, 12), ("4′", ring["C4"], -20, -6),
    ):  # fmt: skip
        b.append(f'<circle cx="{x}" cy="{y}" r="5" fill="#7a68c0"/>')
        b.append(_text(x + dx, y + dy, label, size=15, fill=INK, weight=700))
    # 5' carbon and the phosphate
    b.append(_line(292, 198, 250, 140, "#7a68c0", 3))
    b.append('<circle cx="250" cy="140" r="5" fill="#7a68c0"/>')
    b.append(_text(262, 124, "5′", size=15, weight=700))
    b.append(_line(250, 140, 178, 140, "#b9770a", 3))
    b.append('<circle cx="126" cy="140" r="52" fill="#fdf0d9" stroke="#e69f00" stroke-width="3"/>')
    b.append(_text(126, 136, "phosphate", size=15, fill="#9a6400", weight=700))
    b.append(_text(126, 156, "PO₄", size=14, fill="#9a6400"))
    b.append(_text(126, 74, "negatively charged:", size=12.5, fill=SOFT))
    b.append(_text(126, 90, "DNA is an acid", size=12.5, fill=SOFT))
    # the base on carbon 1'
    b.append(_line(428, 198, 498, 150, "#2f86d6", 3))
    hexagon = " ".join(
        f"{566 + 62 * math.cos(math.radians(a)):.1f},{126 + 62 * math.sin(math.radians(a)):.1f}"
        for a in range(0, 360, 60)
    )
    b.append(f'<polygon points="{hexagon}" fill="#dceefb" stroke="#2f86d6" stroke-width="3"/>')
    b.append(_text(566, 122, "base", size=16, fill="#1b5fa3", weight=700))
    b.append(_text(566, 142, "A, T, G or C", size=13, fill="#1b5fa3"))
    b.append(_text(566, 44, "carries the information", size=12.5, fill=SOFT))
    # 3' hydroxyl and the 2' position
    b.append(_line(318, 276, 318, 330, "#d55e00", 3))
    b.append(_text(318, 352, "OH", size=17, fill=SIGNAL, weight=700))
    b.append(_text(318, 374, "3′ end: the next nucleotide", size=12.5, fill=SOFT))
    b.append(_text(318, 390, "is attached here", size=12.5, fill=SOFT))
    b.append(_line(402, 276, 402, 322, "#7a68c0", 3))
    b.append(_text(402, 344, "H", size=17, fill="#5a4aa8", weight=700))
    b.append(_text(486, 340, "no oxygen at 2′: “deoxy”", size=12.5, fill=SOFT, anchor="start"))
    b.append(_text(486, 356, "(RNA has OH here)", size=12.5, fill=SOFT, anchor="start"))
    b.append(_line(420, 338, 480, 338, FAINT, 1.2, "3 4"))
    b.append(_text(196, 214, "5′ end: the phosphate", size=12.5, fill=SOFT))
    b.append(_text(196, 230, "hangs on carbon 5′", size=12.5, fill=SOFT))
    return _svg(700, 410, b, "A nucleotide: phosphate, deoxyribose and base")


# ---- central dogma -------------------------------------------------------------------------------


def central_dogma_svg() -> str:
    """Flow of sequence information: the general transfers and the special ones (Crick 1970)."""
    b = [_arrow_defs()]
    boxes = (("DNA", 150, "#dceefb", "#2f86d6"), ("RNA", 450, "#f7e1ee", "#c46aa0"),
             ("protein", 750, "#d6f3ea", "#00a878"))  # fmt: skip
    for label, x, fill, stroke in boxes:
        b.append(f'<rect x="{x - 95}" y="110" width="190" height="84" rx="18" fill="{fill}" '
                 f'stroke="{stroke}" stroke-width="3"/>')  # fmt: skip
        b.append(_text(x, 162, label, size=30, weight=700))
    b.append(f'<path d="M250,152 L346,152" stroke="{INK}" stroke-width="4" marker-end="url(#ah)"/>')
    b.append(_text(298, 136, "transcription", size=16, weight=700))
    b.append(f'<path d="M550,152 L646,152" stroke="{INK}" stroke-width="4" marker-end="url(#ah)"/>')
    b.append(_text(598, 136, "translation", size=16, weight=700))
    b.append(f'<path d="M108,110 C80,34 220,34 192,104" fill="none" stroke="{INK}" '
             'stroke-width="4" marker-end="url(#ah)"/>')  # fmt: skip
    b.append(_text(150, 40, "replication", size=16, weight=700))
    b.append(f'<path d="M352,180 C320,228 280,228 248,186" fill="none" stroke="{FAINT}" '
             'stroke-width="2.5" stroke-dasharray="6 6" marker-end="url(#ah-soft)"/>')  # fmt: skip
    b.append(_text(300, 246, "reverse transcription", size=13, fill=SOFT))
    b.append(_text(300, 262, "(retroviruses)", size=12, fill=FAINT))
    b.append(f'<path d="M408,110 C388,50 512,50 492,104" fill="none" stroke="{FAINT}" '
             'stroke-width="2.5" stroke-dasharray="6 6" marker-end="url(#ah-soft)"/>')  # fmt: skip
    b.append(_text(450, 48, "RNA replication (RNA viruses)", size=13, fill=SOFT))
    b.append(f'<path d="M700,206 L600,246" stroke="{SIGNAL}" stroke-width="2.5" '
             'stroke-dasharray="3 6"/>')  # fmt: skip
    b.append(f'<path d="M638,214 l22,22 m0,-22 l-22,22" stroke="{SIGNAL}" stroke-width="3.5"/>')
    b.append(_text(750, 234, "never back from protein", size=13.5, fill=SIGNAL, weight=700))
    b.append(_text(750, 252, "to nucleic acid", size=13.5, fill=SIGNAL, weight=700))
    b.append(_line(60, 290, 110, 290, INK, 4))
    b.append(_text(120, 295, "general transfers (all cells)", size=13, fill=SOFT, anchor="start"))
    b.append(_line(400, 290, 450, 290, FAINT, 2.5, "6 6"))
    b.append(_text(460, 295, "special transfers (Crick 1970)", size=13, fill=SOFT, anchor="start"))
    return _svg(900, 312, b, "The central dogma: DNA to RNA to protein")


# ---- DNA packaging -------------------------------------------------------------------------------


def packaging_svg() -> str:
    """From the double helix to the mitotic chromosome of a eukaryote (schematic)."""
    b = [_arrow_defs()]
    titles = ("DNA double helix", "nucleosomes", "chromatin chain", "loops and domains",
              "mitotic chromosome")  # fmt: skip
    notes = (("2 nm wide", ""), ("146 bp around", "8 histone proteins"),
             ("in cells a disordered", "chain, 5 to 24 nm"), ("folded into loops", "by proteins"),
             ("most compact form,", "two sister chromatids"))  # fmt: skip
    centres = [100, 290, 480, 670, 870]
    for x, title, note in zip(centres, titles, notes, strict=True):
        b.append(_text(x, 30, title, size=15, weight=700))
        b.append(_text(x, 232, note[0], size=13, fill=SOFT))
        b.append(_text(x, 252, note[1], size=13, fill=SOFT))
    for x in centres[:-1]:
        b.append(f'<path d="M{x + 72},120 L{x + 112},120" stroke="{FAINT}" stroke-width="2.5" '
                 'marker-end="url(#ah-soft)"/>')  # fmt: skip
    # (1) double helix: two sine strands with rungs
    for phase, colour in ((0.0, "#2f86d6"), (math.pi, "#e69f00")):
        pts = " ".join(
            f"{40 + t * 2:.1f},{120 + 30 * math.sin(t * 2 * math.pi / 30 + phase):.1f}"
            for t in range(61)
        )
        b.append(f'<polyline points="{pts}" fill="none" stroke="{colour}" stroke-width="4"/>')
    for t in range(3, 60, 5):
        y = 30 * math.sin(t * 2 * math.pi / 30)
        b.append(_line(40 + t * 2, 120 - y, 40 + t * 2, 120 + y, "#9aa8b8", 2))
    # (2) beads on a string
    string = "M214,150 C232,150 236,110 252,110 S276,150 292,150 S316,110 330,110 S352,150 368,150"
    b.append(f'<path d="{string}" fill="none" stroke="#2f86d6" stroke-width="3"/>')
    for cx, cy in ((252, 112), (292, 148), (330, 112)):
        b.append(f'<circle cx="{cx}" cy="{cy}" r="19" fill="#f4d9a0" stroke="#b9770a" '
                 'stroke-width="2"/>')  # fmt: skip
        b.append(f'<ellipse cx="{cx}" cy="{cy}" rx="22" ry="8" fill="none" stroke="#2f86d6" '
                 f'stroke-width="3" transform="rotate(-18 {cx} {cy})"/>')  # fmt: skip
    # (3) irregular chromatin chain: nucleosomes along a wandering path
    path = [(420, 150), (436, 122), (452, 148), (470, 112), (488, 140), (508, 108), (522, 142),
            (540, 124), (530, 166), (504, 170), (476, 172), (448, 176)]  # fmt: skip
    pts = " ".join(f"{x},{y}" for x, y in path)
    b.append(f'<polyline points="{pts}" fill="none" stroke="#2f86d6" stroke-width="2.5"/>')
    for x, y in path:
        b.append(f'<circle cx="{x}" cy="{y}" r="9" fill="#f4d9a0" stroke="#b9770a" '
                 'stroke-width="1.5"/>')  # fmt: skip
    # (4) loops attached to a protein scaffold
    b.append(f'<rect x="604" y="150" width="132" height="12" rx="6" fill="{FAINT}"/>')
    for i in range(5):
        x = 618 + i * 26
        height = (62, 84, 70, 90, 60)[i]
        b.append(f'<path d="M{x - 9},150 C{x - 30},{150 - height} {x + 30},{150 - height} {x + 9},150" '
                 'fill="none" stroke="#2f86d6" stroke-width="5" stroke-linecap="round"/>')  # fmt: skip
    # (5) chromosome: two chromatids joined at the centromere
    for dx in (-17, 17):
        b.append(
            f'<rect x="{870 + dx - 14}" y="52" width="28" height="70" rx="14" fill="#2f86d6"/>'
        )
        b.append(
            f'<rect x="{870 + dx - 14}" y="118" width="28" height="84" rx="14" fill="#2f86d6"/>'
        )
    b.append('<circle cx="870" cy="120" r="10" fill="#e69f00" stroke="#fff" stroke-width="2"/>')
    b.append(_text(928, 124, "centromere", size=12, fill=SOFT, anchor="start"))
    b.append(_line(884, 120, 924, 120, FAINT, 1.2, "3 4"))
    return _svg(1000, 266, b, "Levels of DNA packaging in a eukaryotic cell")


# ---- the genes of phiX174 ------------------------------------------------------------------------


def gene_map_svg(
    genes: Sequence[Gene], length: int, transcripts: Sequence[dict[str, Any]] = ()
) -> str:
    """Linear map of a genome with one row per reading frame; overlapping genes sit on other rows."""
    width, left, right = 1000, 70, 30
    scale = (width - left - right) / length
    rows: dict[tuple[int, int], list[tuple[int, int]]] = {}
    placed = []
    for index, gene in enumerate(genes):
        frame = gene.frame()
        lane = 0
        while any(
            a <= end and start <= b
            for start, end in gene.parts
            for a, b in rows.get((frame, lane), [])
        ):
            lane += 1
        rows.setdefault((frame, lane), []).extend(gene.parts)
        placed.append((gene, frame, lane, GENE_COLOURS[index % len(GENE_COLOURS)]))
    lanes = {
        f: 1 + max((lane for (frame, lane) in rows if frame == f), default=0) for f in range(3)
    }
    row_height, top = 34, 56
    offsets, y = {}, top
    for frame in range(3):
        offsets[frame] = y
        y += lanes[frame] * row_height + 14
    axis_y = y + 6
    b = []
    for frame in range(3):
        height = lanes[frame] * row_height
        b.append(f'<rect x="{left}" y="{offsets[frame] - 4}" width="{width - left - right}" '
                 f'height="{height}" rx="6" fill="{PANEL}"/>')  # fmt: skip
        b.append(_text(left - 10, offsets[frame] + height / 2 + 1, f"frame {frame + 1}", size=12.5,
                       fill=SOFT, anchor="end"))  # fmt: skip
    for gene, frame, lane, colour in placed:
        y0 = offsets[frame] + lane * row_height
        for part_index, (start, end) in enumerate(gene.parts):
            x1, x2 = left + (start - 1) * scale, left + end * scale
            last = part_index == len(gene.parts) - 1
            tip = min(9.0, (x2 - x1) / 2) if last else 0.0
            b.append(
                f'<path d="M{x1:.1f},{y0} L{x2 - tip:.1f},{y0} L{x2:.1f},{y0 + 13} '
                f'L{x2 - tip:.1f},{y0 + 26} L{x1:.1f},{y0 + 26} z" fill="{colour}"/>'
            )
            wide_elsewhere = any((e - s + 1) * scale > 16 for s, e in gene.parts)
            if x2 - x1 > 16:
                b.append(_text((x1 + x2) / 2 - tip / 2, y0 + 18, gene.name, size=13.5, fill="#fff",
                               weight=700))  # fmt: skip
            elif not wide_elsewhere:
                b.append(_text(x1 - 4, y0 + 18, gene.name, size=12, fill=colour, weight=700,
                               anchor="end"))  # fmt: skip
        if len(gene.parts) > 1:  # the gene runs over the end of the record: mark the join
            b.append(_text(left + length * scale + 4, y0 + 18, "↩", size=14, fill=colour,
                           anchor="start"))  # fmt: skip
    b.append(_line(left, axis_y, left + length * scale, axis_y, INK, 2))
    for tick in range(0, length + 1, 500):
        x = left + tick * scale
        b.append(_line(x, axis_y, x, axis_y + 6, INK, 1.5))
        b.append(_text(x, axis_y + 21, f"{tick:,}", size=12, fill=SOFT))
    b.append(_text(left + length * scale, axis_y + 21, f"{length:,}", size=12, fill=SOFT))
    b.append(_text(left + length * scale / 2, axis_y + 40, "position in the genome (nucleotides)",
                   size=12.5, fill=SOFT))  # fmt: skip
    y_t = 16
    for t in transcripts:
        x1, x2 = left + (int(t["start"]) - 1) * scale, left + int(t["end"]) * scale
        b.append(_line(x1, y_t, x2, y_t, "#c46aa0", 3))
        b.append(_text(x2 + 6, y_t + 4, f"mRNA: {t['name']}", size=11.5, fill="#a04f82",
                       anchor="start"))  # fmt: skip
        y_t += 14
    return _svg(width, int(axis_y + 50), b, "Gene map of bacteriophage phiX174")


# ---- genome sizes --------------------------------------------------------------------------------


def genome_sizes_svg(rows: Sequence[tuple[str, int, str]]) -> str:
    """Horizontal bars on a logarithmic axis. ``rows`` = (name, base pairs, note)."""
    width, left, right, top, step = 1000, 250, 190, 16, 44
    low, high = 3.0, 10.0
    scale = (width - left - right) / (high - low)
    b = []
    height = top + len(rows) * step + 46
    for power in range(int(low), int(high) + 1):
        x = left + (power - low) * scale
        b.append(_line(x, top, x, top + len(rows) * step, LINE, 1))
        label = {3: "1 kb", 4: "10 kb", 5: "100 kb", 6: "1 Mb", 7: "10 Mb", 8: "100 Mb",
                 9: "1 Gb", 10: "10 Gb"}[power]  # fmt: skip
        b.append(_text(x, top + len(rows) * step + 18, label, size=12, fill=SOFT))
    for i, (name, size, note) in enumerate(rows):
        y = top + i * step + 8
        x2 = left + (math.log10(size) - low) * scale
        colour = ACCENT if "this project" in note else "#9db8cf"
        b.append(f'<rect x="{left}" y="{y}" width="{x2 - left:.1f}" height="24" rx="5" '
                 f'fill="{colour}"/>')  # fmt: skip
        b.append(_text(left - 12, y + 17, name, size=14, weight=600, anchor="end"))
        b.append(_text(x2 + 8, y + 12, f"{size:,} bp", size=13, weight=700, anchor="start"))
        b.append(_text(x2 + 8, y + 27, note, size=11, fill=FAINT, anchor="start"))
    b.append(_text(left + (high - low) * scale / 2, height - 6,
                   "genome size, logarithmic scale (each line is ten times the one before)",
                   size=12.5, fill=SOFT))  # fmt: skip
    return _svg(width, height, b, "Genome sizes on a logarithmic scale")


# ---- base composition ----------------------------------------------------------------------------


def composition_svg(rows: Sequence[tuple[str, dict[str, int]]]) -> str:
    """Base composition of one strand of each genome, in percent (axis starts at zero)."""
    width, height, left, bottom = 760, 300, 60, 250
    top_value = 35.0
    scale = (bottom - 46) / top_value
    b = []
    for value in range(0, 36, 5):
        y = bottom - value * scale
        b.append(_line(left, y, width - 20, y, LINE, 1))
        b.append(_text(left - 8, y + 4, f"{value}%", size=12, fill=SOFT, anchor="end"))
    group = (width - left - 40) / len(rows)
    for g, (name, counts) in enumerate(rows):
        total = sum(counts.values())
        x0 = left + g * group + 30
        for k, base in enumerate("ATGC"):
            share = 100.0 * counts[base] / total
            x, h = x0 + k * 62, share * scale
            b.append(f'<rect x="{x}" y="{bottom - h:.1f}" width="46" height="{h:.1f}" rx="4" '
                     f'fill="{BASE_COLOUR[base]}"/>')  # fmt: skip
            b.append(_text(x + 23, bottom - h - 6, f"{share:.1f}", size=12.5, weight=700))
            b.append(_text(x + 23, bottom + 18, base, size=14, weight=700))
        b.append(_text(x0 + 116, bottom + 40, name, size=13, fill=SOFT))
    b.append(_text(16, 18, "share of the letters of the sequenced strand", size=12.5, fill=SOFT,
                   anchor="start"))  # fmt: skip
    return _svg(width, height, b, "Base composition of the two genomes")


# ---- sequencing ----------------------------------------------------------------------------------


def sanger_gel_svg(sequence: str, first_position: int) -> str:
    """Chain-termination sequencing of a real stretch: four lanes, one band per stopped copy."""
    n = len(sequence)
    width, left, top, lane, row = 640, 150, 54, 78, 19
    height = top + n * row + 70
    b = [f'<rect x="{left - 16}" y="{top - 22}" width="{4 * lane + 32}" height="{n * row + 40}" '
         f'rx="10" fill="#101c2e"/>']  # fmt: skip
    for k, base in enumerate("ACGT"):
        x = left + k * lane + lane / 2
        b.append(_text(x, top - 30, f"dd{base}TP", size=13.5, fill=BASE_COLOUR[base], weight=700))
    for i, base in enumerate(sequence):
        y = top + (n - 1 - i) * row + 8  # the shortest copy runs furthest: bottom of the gel
        x = left + "ACGT".index(base) * lane + 12
        b.append(f'<rect x="{x}" y="{y - 5}" width="{lane - 24}" height="8" rx="3" '
                 f'fill="{BASE_COLOUR[base]}"/>')  # fmt: skip
        b.append(_text(left + 4 * lane + 28, y + 4, f"{base}", size=13, fill=BASE_COLOUR[base],
                       weight=700, anchor="start"))  # fmt: skip
        if i in (0, n - 1):
            b.append(_text(left - 26, y + 4, f"{i + 1} nt", size=12, fill=SOFT, anchor="end"))
    b.append(f'<path d="M{left + 4 * lane + 62},{top + n * row} L{left + 4 * lane + 62},{top}" '
             f'stroke="{INK}" stroke-width="2" marker-end="url(#ah)"/>')  # fmt: skip
    b.insert(0, _arrow_defs())
    b.append(_text(left + 4 * lane + 76, top + n * row / 2, "read", size=12.5, fill=SOFT,
                   anchor="start"))  # fmt: skip
    b.append(_text(left + 4 * lane + 76, top + n * row / 2 + 20, "5′ → 3′", size=12.5, fill=SOFT,
                   anchor="start"))  # fmt: skip
    b.append(_text(left - 26, top + n * row / 2, "longer", size=12, fill=FAINT, anchor="end"))
    b.append(
        _text(left - 26, top + n * row / 2 + 19, "copies ↑", size=12, fill=FAINT, anchor="end")
    )
    b.append(_text(width / 2, height - 30,
                   f"Read from the bottom: {sequence}", size=14, weight=700))  # fmt: skip
    b.append(_text(width / 2, height - 10,
                   f"positions {first_position:,} to {first_position + n - 1:,} of phiX174",
                   size=12, fill=SOFT))  # fmt: skip
    return _svg(width, height, b, "Chain-termination sequencing gel of a real sequence")


def shotgun_svg(starts: Sequence[int], lengths: Sequence[int], genome_length: int) -> str:
    """Reads of a shotgun experiment at their true positions, and the coverage they give."""
    width, left, right = 1000, 40, 40
    scale = (width - left - right) / genome_length
    order = sorted(range(len(starts)), key=lambda i: starts[i])
    lane_end: list[int] = []
    lane_of = {}
    for i in order:
        for lane, end in enumerate(lane_end):
            if starts[i] > end + 2:
                lane_end[lane] = starts[i] + lengths[i]
                lane_of[i] = lane
                break
        else:
            lane_of[i] = len(lane_end)
            lane_end.append(starts[i] + lengths[i])
    depth = [0] * genome_length
    for s, length in zip(starts, lengths, strict=True):
        for p in range(s, min(genome_length, s + length)):
            depth[p] += 1
    top, row = 78, 11
    lanes = len(lane_end)
    b = [
        _text(left, 18, "the genome (unknown to the sequencer)", size=13, fill=SOFT, anchor="start")
    ]
    b.append(f'<rect x="{left}" y="28" width="{genome_length * scale:.1f}" height="12" rx="6" '
             f'fill="{INK}"/>')  # fmt: skip
    b.append(_text(left, 66, f"{len(starts)} reads at random positions", size=13, fill=SOFT,
                   anchor="start"))  # fmt: skip
    for i in order:
        x, w = left + starts[i] * scale, lengths[i] * scale
        y = top + lane_of[i] * row
        b.append(f'<rect x="{x:.1f}" y="{y}" width="{w:.1f}" height="7" rx="3.5" '
                 f'fill="{ACCENT}" opacity="0.85"/>')  # fmt: skip
    base_y = top + lanes * row + 76
    peak = max(depth)
    b.append(_text(left, base_y - 56, "coverage: how many reads cover each letter", size=13,
                   fill=SOFT, anchor="start"))  # fmt: skip
    points = [f"{left},{base_y}"]
    for p, d in enumerate(depth):
        points.append(f"{left + p * scale:.1f},{base_y - 44 * d / peak:.1f}")
    points.append(f"{left + genome_length * scale:.1f},{base_y}")
    b.append(f'<polygon points="{" ".join(points)}" fill="#cfe5f3" stroke="{ACCENT}" '
             'stroke-width="1.5"/>')  # fmt: skip
    b.append(_text(left + genome_length * scale + 6, base_y - 40, f"{peak}×", size=12, fill=SOFT,
                   anchor="start"))  # fmt: skip
    b.append(_text(left + genome_length * scale + 6, base_y + 2, "0", size=12, fill=SOFT,
                   anchor="start"))  # fmt: skip
    uncovered = sum(1 for d in depth if d == 0)
    if uncovered:
        b.append(_text(width - right, base_y + 20, f"{uncovered} letters are covered by no read",
                       size=12.5, fill=SIGNAL, weight=700, anchor="end"))  # fmt: skip
    return _svg(width, base_y + 30, b, "Shotgun reads and their coverage of the genome")


# ---- the genetic code ----------------------------------------------------------------------------


def codon_table_html() -> str:
    """The standard genetic code as the classic 4 x 4 table (RNA letters), coloured by class."""
    bases = "UCAG"
    out = ['<table class="codon-table"><thead><tr><th>1st</th>']
    out += [f"<th>{b}</th>" for b in bases]
    out.append("<th>3rd</th></tr></thead><tbody>")
    for first in bases:
        for row, third in enumerate(bases):
            out.append("<tr>")
            if row == 0:
                out.append(f'<th rowspan="4" class="codon-first">{first}</th>')
            for second in bases:
                codon = first + second + third
                amino = CODON_TABLE[codon.replace("U", "T")]
                three, _, cls = AMINO_ACIDS[amino]
                mark = " codon-start" if codon == "AUG" else ""
                out.append(
                    f'<td class="codon{mark}" style="background:{CLASS_COLOUR[cls]}">'
                    f"<code>{codon}</code> <b>{three}</b></td>"
                )
            out.append(f"<th>{third}</th></tr>")
    out.append("</tbody></table>")
    return "".join(out)


CODON_TABLE_CSS = """
.codon-table { border-collapse: collapse; font-family: 'Segoe UI', system-ui, sans-serif;
  font-size: 14px; color: #10243b; width: 100%; }
.codon-table th { padding: 3px 8px; color: #4a5b70; font-weight: 700; }
.codon-table td { padding: 3px 10px; border: 1px solid #fff; white-space: nowrap; }
.codon-table code { font-family: Consolas, monospace; font-size: 13px; color: #4a5b70;
  background: none; }
.codon-table .codon-first { font-size: 18px; color: #10243b; }
.codon-table .codon-start { outline: 2px solid #00a878; outline-offset: -2px; }
.codon-table tr:nth-child(4n) td, .codon-table tr:nth-child(4n) th { border-bottom: 8px solid #fff; }
"""
