"""The slide deck: one content list, rendered to .pptx (python-pptx) and to a PDF copy.

Every slide title states its takeaway; speaker notes are in plain language; every number comes from
``Facts``. The PDF copy (headless Chrome) is the fallback if no PowerPoint is available (Q3, Q5).
"""

from __future__ import annotations

import base64
import html
import io
import math
from dataclasses import dataclass, field
from pathlib import Path

from matplotlib.figure import Figure
from PIL import Image
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Inches, Pt

from genomeasm.appkit import content, figs, lab
from genomeasm.appkit.paths import FIGURES
from genomeasm.docgen import md
from genomeasm.docgen.facts import Facts, name, p_text

W, H = 13.333, 7.5
TITLE_BOX = (0.6, 0.35, 12.1, 1.15)
BODY_TOP, BODY_H = 1.7, 5.15
INK, ACCENT, GREY = (
    RGBColor(0x1D, 0x27, 0x33),
    RGBColor(0x00, 0x72, 0xB2),
    RGBColor(0x66, 0x70, 0x7A),
)


@dataclass
class Slide:
    title: str
    bullets: list[str] = field(default_factory=list)
    image: bytes | None = None  # PNG
    equations: list[str] = field(default_factory=list)
    table: list[list[str]] | None = None  # the first row is the header
    source: str = ""
    notes: str = ""
    kind: str = "content"  # "title" or "content"


def png(fig: Figure, dpi: int = 200) -> bytes:
    return figs.live_frame(fig, dpi)


def figure_png(stem: str) -> bytes:
    return (FIGURES / f"{stem}.png").read_bytes()


def math_png(expression: str, size: float = 26) -> bytes:
    fig = Figure(figsize=(0.1, 0.1))
    fig.text(0, 0, f"${md.normalise(expression)}$", fontsize=size)
    buffer = io.BytesIO()
    fig.savefig(
        buffer, format="png", dpi=200, bbox_inches="tight", pad_inches=0.05, transparent=True
    )
    return buffer.getvalue()


def image_size(data: bytes) -> tuple[int, int]:
    with Image.open(io.BytesIO(data)) as image:
        return image.size


# ---- fitting text -----------------------------------------------------------------------------------


def text_height(bullets: list[str], width_in: float, font_pt: float) -> float:
    """Estimated height in inches of bulleted text (about 0.5 em per character, 1.2 line spacing)."""
    chars_per_line = max(8.0, width_in * 72.0 / (font_pt * 0.5))
    lines = sum(math.ceil((len(b) + 2) / chars_per_line) for b in bullets)
    return lines * font_pt * 1.2 / 72.0 + len(bullets) * 0.12


def fit_font(bullets: list[str], width_in: float, height_in: float) -> int:
    for size in (22, 20, 18, 16, 14, 12):
        if text_height(bullets, width_in, size) <= height_in:
            return size
    return 12


def layout(slide: Slide) -> dict[str, tuple[float, float, float, float]]:
    """Boxes (left, top, width, height in inches) of the parts of a content slide."""
    has_text = bool(slide.bullets)
    has_visual = bool(slide.image or slide.equations or slide.table)
    if not has_visual:
        return {"text": (0.7, BODY_TOP, W - 1.4, BODY_H)}
    if not has_text:
        return {"visual": (0.7, BODY_TOP, W - 1.4, BODY_H)}
    if slide.image:
        w_px, h_px = image_size(slide.image)
        if w_px / h_px >= 1.9:  # a wide chart: text above, chart below, both full width
            return {"text": (0.6, 1.5, 12.1, 1.95), "visual": (0.6, 3.5, 12.1, 3.4)}
    return {"text": (0.6, BODY_TOP, 5.7, BODY_H), "visual": (6.5, BODY_TOP, 6.3, BODY_H)}


# ---- pptx -------------------------------------------------------------------------------------------


def _add_text(slide_obj, box, bullets, font, bullet=True, bold=False, colour=INK):  # type: ignore[no-untyped-def]
    left, top, width, height = box
    frame = slide_obj.shapes.add_textbox(
        Inches(left), Inches(top), Inches(width), Inches(height)
    ).text_frame
    frame.word_wrap = True
    for k, text in enumerate(bullets):
        paragraph = frame.paragraphs[0] if k == 0 else frame.add_paragraph()
        paragraph.space_after = Pt(8)
        run = paragraph.add_run()
        run.text = ("• " if bullet else "") + text
        run.font.size = Pt(font)
        run.font.bold = bold
        run.font.color.rgb = colour
    return frame


def _place_image(slide_obj, data: bytes, box) -> None:  # type: ignore[no-untyped-def]
    left, top, width, height = box
    w_px, h_px = image_size(data)
    scale = min(width / w_px, height / h_px)
    w_in, h_in = w_px * scale, h_px * scale
    slide_obj.shapes.add_picture(
        io.BytesIO(data),
        Inches(left + (width - w_in) / 2),
        Inches(top + (height - h_in) / 2),
        Inches(w_in),
        Inches(h_in),
    )


def _add_table(slide_obj, rows, box) -> None:  # type: ignore[no-untyped-def]
    left, top, width, height = box
    font = 16 if len(rows) <= 7 else 13 if len(rows) <= 11 else 11
    row_h = min(0.5, height / len(rows))
    shape = slide_obj.shapes.add_table(
        len(rows), len(rows[0]), Inches(left), Inches(top), Inches(width), Inches(row_h * len(rows))
    )
    for r, row in enumerate(rows):
        for c, cell_text in enumerate(row):
            cell = shape.table.cell(r, c)
            cell.text = str(cell_text)
            for paragraph in cell.text_frame.paragraphs:
                for run in paragraph.runs:
                    run.font.size = Pt(font)
                    run.font.bold = r == 0


def _add_visual(slide_obj, slide: Slide, box) -> None:  # type: ignore[no-untyped-def]
    left, top, width, height = box
    if slide.table:
        _add_table(slide_obj, slide.table, box)
    elif slide.equations:
        gap = height / len(slide.equations)
        for k, expression in enumerate(slide.equations):
            _place_image(slide_obj, math_png(expression), (left, top + k * gap, width, gap * 0.9))
    elif slide.image:
        _place_image(slide_obj, slide.image, box)


def build_pptx(slides: list[Slide], path: Path) -> Path:
    deck = Presentation()
    deck.slide_width, deck.slide_height = Inches(W), Inches(H)
    for number, slide in enumerate(slides, start=1):
        obj = deck.slides.add_slide(deck.slide_layouts[6])
        if slide.kind == "title":
            _add_text(
                obj,
                (1.0, 2.2, W - 2.0, 1.8),
                [slide.title],
                40,
                bullet=False,
                bold=True,
                colour=ACCENT,
            )
            _add_text(obj, (1.0, 4.3, W - 2.0, 2.2), slide.bullets, 22, bullet=False)
            if slide.image:
                _place_image(obj, slide.image, (W - 5.0, 0.3, 4.4, 1.6))
        else:
            left, top, width, height = TITLE_BOX
            frame = obj.shapes.add_textbox(
                Inches(left), Inches(top), Inches(width), Inches(height)
            ).text_frame
            frame.word_wrap = True
            run = frame.paragraphs[0].add_run()
            run.text = slide.title
            run.font.size = Pt(28 if len(slide.title) <= 70 else 24)
            run.font.bold = True
            run.font.color.rgb = ACCENT
            boxes = layout(slide)
            if slide.bullets:
                box = boxes["text"]
                _add_text(obj, box, slide.bullets, fit_font(slide.bullets, box[2], box[3]))
            if "visual" in boxes:
                _add_visual(obj, slide, boxes["visual"])
            foot = obj.shapes.add_textbox(
                Inches(0.6), Inches(7.0), Inches(12.1), Inches(0.35)
            ).text_frame
            foot_run = foot.paragraphs[0].add_run()
            foot_run.text = f"{number}   {slide.source}".strip()
            foot_run.font.size = Pt(11)
            foot_run.font.color.rgb = GREY
            foot.paragraphs[0].alignment = PP_ALIGN.LEFT
        obj.notes_slide.notes_text_frame.text = slide.notes
    path.parent.mkdir(parents=True, exist_ok=True)
    deck.save(str(path))
    return path


# ---- PDF copy (HTML printed by Chrome) ------------------------------------------------------------------

SLIDE_CSS = """
@page { size: 1280px 720px; margin: 0; }
body { margin: 0; font-family: 'Segoe UI', Arial, sans-serif; color: #1d2733; }
.slide { width: 1280px; height: 720px; box-sizing: border-box; padding: 34px 58px; position: relative; break-after: page; overflow: hidden; background: white; }
.slide h2 { font-size: 34px; color: #0072B2; margin: 0 0 22px 0; line-height: 1.15; }
.row { display: flex; gap: 36px; align-items: center; height: 520px; }
.text { flex: 1 1 45%; } .vis { flex: 1 1 55%; text-align: center; max-height: 520px; }
.vis img { max-width: 100%; max-height: 500px; } .vis .eq { margin: 26px 0; } .vis .eq img { height: 78px; max-width: 100%; } .row.wide { flex-direction: column; align-items: stretch; height: 560px; } .row.wide .vis img { max-height: 340px; } ul { margin: 0; padding-left: 26px; } li { margin: 0 0 12px 0; line-height: 1.28; }
table { border-collapse: collapse; margin: auto; font-size: 20px; } th, td { border: 1px solid #c9d2db; padding: 6px 12px; text-align: left; } th { background: #e8eef4; }
.foot { position: absolute; bottom: 14px; left: 58px; font-size: 14px; color: #66707a; }
.title h1 { font-size: 52px; color: #0072B2; margin: 150px 0 20px 0; line-height: 1.1; } .title li { font-size: 26px; list-style: none; margin-bottom: 8px; } .title ul { padding: 0; }
.title img { position: absolute; right: 58px; top: 30px; width: 420px; }
"""


def _slide_html(slide: Slide, number: int, font_px: int, wide: bool = False) -> str:
    items = "".join(f"<li>{html.escape(b)}</li>" for b in slide.bullets)
    if slide.kind == "title":
        image = ""
        if slide.image:
            image = f"<img src='data:image/png;base64,{base64.b64encode(slide.image).decode()}'>"
        return f"<section class='slide title'>{image}<h1>{html.escape(slide.title)}</h1><ul>{items}</ul></section>"
    visual = ""
    if slide.table:
        head, *body = slide.table
        visual = "<table><tr>" + "".join(f"<th>{html.escape(str(c))}</th>" for c in head) + "</tr>"
        visual += (
            "".join(
                "<tr>" + "".join(f"<td>{html.escape(str(c))}</td>" for c in row) + "</tr>"
                for row in body
            )
            + "</table>"
        )
    elif slide.equations:
        visual = "".join(md.math_html(e, display=True) for e in slide.equations)
    elif slide.image:
        visual = f"<img src='data:image/png;base64,{base64.b64encode(slide.image).decode()}'>"
    text = (
        f"<div class='text'><ul style='font-size:{font_px}px'>{items}</ul></div>"
        if slide.bullets
        else ""
    )
    vis = f"<div class='vis'>{visual}</div>" if visual else ""
    return (
        f"<section class='slide'><h2>{html.escape(slide.title)}</h2><div class='row{' wide' if wide else ''}'>{text}{vis}</div>"
        f"<div class='foot'>{number}&nbsp;&nbsp;{html.escape(slide.source)}</div></section>"
    )


def build_pdf(slides: list[Slide], out_dir: Path) -> Path:
    sections = []
    for number, slide in enumerate(slides, start=1):
        boxes = layout(slide)
        box = boxes.get("text", (0, 0, 5.7, BODY_H))
        font_px = round(fit_font(slide.bullets, box[2], box[3]) * 96 / 72) if slide.bullets else 24
        sections.append(_slide_html(slide, number, font_px, wide=box[2] > 10 and "visual" in boxes))
    page = f"<!doctype html><html><head><meta charset='utf-8'><style>{SLIDE_CSS}</style></head><body>{''.join(sections)}<script>const s=new URLSearchParams(location.search).get('s');if(s){{document.querySelectorAll('.slide').forEach((e,i)=>{{if(i!=s-1)e.style.display='none'}})}}</script></body></html>"
    html_file = out_dir / "deck.html"
    html_file.write_text(page, encoding="utf-8", newline=chr(10))
    pdf = out_dir / "deck.pdf"
    md.print_pdf(html_file, pdf)
    return pdf


# ---- the content ------------------------------------------------------------------------------------


def build_slides(F: Facts, author: str = "") -> list[Slide]:
    ucb, eps = F.verdict_ucb, F.verdict_eps
    rt = F.test("alns_ucb1_rt", "alns_random")
    rt_roulette = F.test("alns_ucb1_rt", "alns_roulette")
    toy = content.make_toy(list(content.TOY_READS), 1)
    orders = content.all_orders(toy)
    helix = png(figs.dna_helix(lab.genome_sequence("NC_001422.1")[:24], title=" "), 120)
    cfg = F.config
    arm_rows = [
        ["Destroy (remove)", "Repair (re-insert)"],
        ["random reads", "greedy insertion"],
        ["worst-connected reads", "regret insertion"],
        ["a consecutive stretch", "(each destroy with each repair = 6 tools)"],
    ]
    bio_rows = [["Solver", "Genome fraction", "Identity", "Contigs", "N50", "Misjoins"]] + [
        [name(a), f"{F.bio(a, 'genome_fraction'):.3f}", f"{F.bio(a, 'identity'):.3f}", f"{F.bio(a, 'n_contigs'):.0f}", f"{F.bio(a, 'n50'):.0f}", f"{F.bio(a, 'misjoins'):.0f}"]
        for a in ("greedy", "tabu_cl", "alns_random", "alns_ucb1", "lkh3")
    ]  # fmt: skip
    n_rows = [["Reads n", "Orders n!"]] + [
        [
            str(n),
            f"{math.factorial(n):,}"
            if n <= 15
            else f"a number with {len(str(math.factorial(n)))} digits",
        ]
        for n in (5, 10, 15, 20, 50, 100)
    ]
    design_rows = [
        ["Item", "Choice"],
        ["Genomes", "phiX174 and lambda (real, NCBI)"],
        ["Read sets", f"{F.n_settings} settings x {F.seeds} seeds = {F.n_sets}"],
        ["Budget", f"{F.budget:g} s wall-clock for every solver"],
        ["Tuning", "24 random trials on separate data"],
        ["Statistics", "paired Wilcoxon, Holm, A12"],
    ]
    return [
        Slide(
            "Can a learning rule choose repair tools better than chance?",
            [
                "DNA fragment assembly as a travelling salesman problem",
                "Greedy, Tabu Search, LKH-3 and a bandit-driven ALNS",
                "Project presentation, 2026" + (" · " + author if author else ""),
            ],
            image=helix,
            kind="title",
            notes="Welcome. The project puts the pieces of a genome back in order, treats it as a travelling salesman problem, and tests one precise question: does a learning rule that picks repair tools beat picking them at random?",
        ),
        Slide(
            "A genome is read as thousands of shredded, overlapping pieces",
            [
                "DNA is a word over four letters A, C, G, T; a genome is the whole word",
                "A machine returns short reads, taken at random places, with typing mistakes",
                "Assembly means putting the reads back in the right order",
                "We use real virus genomes: phiX174 (5,386 letters) and lambda (48,502)",
            ],
            image=helix,
            source="Data: NCBI RefSeq NC_001422.1, NC_001416.1",
            notes="Think of a book shredded in many copies: each scrap has a few lines, and neighbouring scraps share words at their edges. The picture is a real stretch of the phiX174 virus, drawn as a double helix where A pairs with T and G with C.",
        ),
        Slide(
            "Two reads fit when the end of one equals the start of the other",
            [
                "The overlap is the number of letters they share at the join",
                "Longer overlap: more letters saved and a more trustworthy neighbour",
                "The overlap from A to B differs from B to A: the table is not symmetric",
                "Read errors force us to allow a few mismatches",
            ],
            image=png(figs.toy_graph(toy.instance, toy.labels, [0, 1, 2]), 130),
            source="Toy example: reads ATGCG, GCGTA, GTACC",
            notes="Reads A, B and C: the last three letters of A equal the first three of B, so gluing B after A saves three letters. The red arrows show the best order.",
        ),
        Slide(
            "Ordering the pieces to save the most letters is a travelling salesman problem",
            [
                "Fitness F: total overlap of neighbouring reads (maximise)",
                "L: length of the merged text, the same goal seen from the other side",
                "Add a dummy node 0: a tour through all nodes costs exactly L",
                "The costs are not symmetric: an asymmetric TSP (ATSP)",
            ],
            equations=[
                r"F(\pi)=\sum_{k=1}^{n-1}\mathrm{ov}(\pi_k,\pi_{k+1})",
                r"L(\pi)=\sum_{i}\ell_i-F(\pi)",
                r"c(i,j)=\ell_j-\mathrm{ov}(i,j)",
            ],
            source="Proof: docs/OBJECTIVE_FUNCTION.md, report section 3",
            notes="The key fact, proved in the report: a round trip through all reads costs exactly the length of the merged text. So the shortest round trip is the best assembly. That is why we can use travelling-salesman solvers.",
        ),
        Slide(
            "The toy example: three reads, six orders, one best",
            [
                f"Best order {orders[0]['order']}: F = {orders[0]['F']}, merged text of {orders[0]['L']} letters",
                f"Worst orders save only {orders[-1]['F']} letter",
                "Total length of the reads: 15; L = 15 - F",
            ],
            table=[["Order", "F", "L"]]
            + [[str(r["order"]), str(r["F"]), str(r["L"])] for r in orders],
            source="Computed live in the application, page Learn",
            notes="Six orders of three reads; you can type any order in the application and see F change. ABC glues all three reads into the nine-letter text ATGCGTACC.",
        ),
        Slide(
            "The number of orders explodes, so we need clever search",
            [
                "Real instances have hundreds of reads",
                f"A computer checking a billion orders per second needs about {math.factorial(20) / 1e9 / (365.25 * 86400):.0f} years for 20 reads",
                "Hence heuristics: greedy, tabu search, LKH-3, ALNS",
            ],
            table=n_rows,
            source="n! = number of orders of n reads",
            notes="Factorials grow faster than any computer can follow. For 20 reads there are more than two quintillion orders.",
        ),
        Slide(
            "A high F is not the genome, so we also grade against the real one",
            [
                "F counts letters saved, not correct letters",
                "With read errors, the best order we find can score higher than the true order",
                "We measure identity, genome fraction, contigs, N50 and misjoins against the answer key",
                "Reads come with known positions, so every answer can be graded",
            ],
            image=figure_png("fig03_overlap_heatmap"),
            source="Figure: overlap matrix sorted by true position",
            notes="The bright diagonal band are overlaps between true neighbours. Because errors break some overlaps, F alone can mislead; that is why the quality against the real genome is reported as well.",
        ),
        Slide(
            "Real genomes are loaded once, checked, and then everything runs offline",
            [
                "Download from NCBI once; check name, kind and length against NCBI's facts",
                "Save with a SHA-256 fingerprint; every later run verifies it first",
                "Nothing is written if a check fails; no data is ever made up",
                "The application and all experiments need no internet",
            ],
            image=png(figs.pipeline_diagram(), 130),
            source="docs/DATA_LOADING.md",
            notes="This is the data pipeline: the orange part runs once with internet, the green part on every run. A single changed letter in a stored genome would be noticed.",
        ),
        Slide(
            "Four solvers, and ALNS comes in six selection variants",
            [
                "Greedy: take the largest compatible overlaps first",
                "Tabu search: best small change, forbid undoing it for a while",
                "LKH-3: the reference TSP solver, verified, never replaced silently",
                "ALNS: destroy part of the order, repair it, accept sometimes worse",
            ],
            source="Methods, report section 4",
            notes="Greedy is fast and surprisingly good. Tabu search improves it locally. LKH-3 is the strong reference. ALNS is where the research question lives.",
        ),
        Slide(
            "ALNS breaks the order and repairs it with one of six tools",
            [
                "Three ways to remove reads times two ways to put them back",
                "Which tool is best changes during the search",
                "Worse results are sometimes accepted, early in the run more than late",
                "The question: how should the tool be chosen?",
            ],
            table=arm_rows,
            source="Operators: decision D-034",
            notes="Each step takes out some reads and puts them back. The six tools are the combinations of three removal methods and two insertion methods.",
        ),
        Slide(
            "A bandit learns which tool to use; we compare five rules",
            [
                "Random choice: the baseline of the claim",
                "Round robin and the classical roulette wheel",
                "Epsilon-greedy: usually the best tool, sometimes a random one",
                "UCB1: average reward plus a bonus for rarely tried tools",
            ],
            equations=[r"\mathrm{UCB1}(a)=\bar{x}_a+c\sqrt{\frac{2\ln N}{n_a}}"],
            source="Policies: decision D-038",
            notes="A bandit faces the exploration-exploitation trade-off: use the tool that looks best, or try another to learn more. The reward is large for a new best order, small for an improvement.",
        ),
        Slide(
            f"A fair test: {F.n_sets} read sets, equal {F.budget:g} s budgets, tuning on separate data",
            [
                "Same read sets and same wall-clock time for every solver",
                "Hyper-parameters tuned on different genome positions and seeds",
                "The claim and the two primary tests were fixed before the runs",
                "Differences judged on the same read set: paired tests",
            ],
            table=design_rows,
            source=f"Configuration: configs/phase5_experiments.yaml ({cfg['seeds']['evaluation']} seeds per setting)",
            notes="This is the experimental design. Every solver sees the same read sets, tuning never touches evaluation data, and the statistics are paired and corrected for multiple tests.",
        ),
        Slide(
            f"LKH-3 {'is never beaten' if F.lkh_never_beaten else 'is rarely beaten'}; ALNS clearly beats greedy and tabu search",
            [
                f"Mean gap to the best known: LKH-3 {F.mean_gap('lkh3'):.2f}%, ALNS random {F.mean_gap('alns_random'):.2f}%, tabu {F.mean_gap('tabu'):.2f}%, greedy {F.mean_gap('greedy'):.2f}%",
                f"Greedy already finds the best known solution in {100 * F.greedy_equals_best:.0f}% of the sets",
                "Differences appear on noisy, deep read sets",
            ],
            image=figure_png("fig05_final_fitness"),
            source="results/phase5, figure 5",
            notes="The gap is how far a solver ends below the best solution found by any solver on the same read set. LKH-3 is the best everywhere. ALNS comes next. Note that best known is not a proven optimum.",
        ),
        Slide(
            "The claim is not supported: bandits do not beat random choice",
            [
                f"UCB1 vs random: {p_text(ucb['p_holm'])}, A12 = {ucb['a12']:.2f} (0.5 = no difference)",
                f"Epsilon-greedy vs random: {p_text(eps['p_holm'])}, A12 = {eps['a12']:.2f}: slightly worse",
                "Round robin and roulette are slightly better than random",
                "All effects are small",
            ],
            image=figure_png("fig07_ablation"),
            source="results/phase5, figure 7",
            notes="This is the honest answer to the research question: no. The claim was fixed in advance, and the data do not support it. A negative result is still a result.",
        ),
        Slide(
            "We found a bug in our own UCB1; the fix is reported as exploratory",
            [
                f"The first operator pair received {100 * F.arm_zero_share:.0f}% of UCB1's uses",
                "Cause: equal indices, so the lowest index always wins",
                f"With random ties: A12 = {rt.a12:.2f} vs random, {p_text(rt_roulette.p_holm)} vs roulette",
                "Added after the results: exploratory, never presented as confirmation",
            ],
            image=figure_png("fig08_arm_probability"),
            source="results/phase5, figure 8; decision D-046",
            notes="In the right-hand picture the learning rule collapses onto the first tool. We kept the original result, added the corrected version afterwards, and labelled it as exploratory.",
        ),
        Slide(
            "All solvers rebuild the same genome, with no misjoins",
            [
                f"Genome fraction about {F.bio('lkh3', 'genome_fraction'):.3f} for every solver",
                "Stronger solvers: slightly higher identity, fewer contigs",
                "Fitness differences are small in absolute terms",
            ],
            table=bio_rows,
            source="Medians over all read sets; D-043",
            notes="Quality against the real genome: nothing is mis-joined, and the solvers differ mostly in how many pieces remain.",
        ),
        Slide(
            f"Time: LKH-3 needs preprocessing; live application runs end in under {lab.MAX_TOTAL_SECONDS:.0f} s",
            [
                "The overlap graph and LKH-3's setup grow fastest with the number of reads",
                f"Worst measured interactive run in the application: {F.worst_live_seconds:.0f} s (rule: 30 s)",
                "Larger runs are shown from precomputed results, clearly labelled",
            ],
            image=figure_png("fig10_runtime_vs_n"),
            source="results/phase5 figure 10; results/phase6",
            notes="The application keeps live runs under 30 seconds on this laptop; the largest experiments are shown from stored results.",
        ),
        Slide(
            "Limits: two viruses, substitution errors, and best known is not optimal",
            [
                "Real data have insertions, deletions, both strands and repeats",
                "Best known is the best solution we found; LKH-3 is a heuristic too",
                "Two post-hoc variants were added after the first results",
                "Tuning on 18 read sets did not always transfer",
            ],
            source="Threats to validity, report section 10",
            notes="We state the limits openly. The most important: the result concerns simulated reads from two small genomes.",
        ),
        Slide(
            "Conclusion: a rigorous pipeline and an honest negative result",
            [
                "The ATSP encoding is exact and proved",
                "ALNS beats greedy and tabu search; LKH-3 is the reference",
                "Learning the tool choice did not beat random mixing here",
                "Everything is reproducible: seeds, hashes, commit, one command",
            ],
            source="Report sections 9 to 11",
            notes="To summarise: a sound method, a fair comparison, and the answer to the question is no, with the reasons and limits stated.",
        ),
        Slide(
            "Live demo: the application shows every step",
            [
                "Guide and Learn: the story, a toy example, DNA and RNA pictures",
                "Solver Lab: a live run with convergence and tool choice",
                "Experiments: the verdict and all statistics",
                "Assembly Viewer and Reproducibility: grading and provenance",
            ],
            image=helix,
            source="Start: make app, then http://localhost:8501",
            notes="Switch to the application. Suggested path: Learn toy example, Solver Lab run, Experiments verdict, Assembly Viewer. The Guide page has the exact script.",
        ),
    ]


def build_deck(out_dir: Path, F: Facts, author: str = "") -> tuple[Path, Path]:
    slides = build_slides(F, author)
    out_dir.mkdir(parents=True, exist_ok=True)
    return build_pptx(slides, out_dir / "deck.pptx"), build_pdf(slides, out_dir)
