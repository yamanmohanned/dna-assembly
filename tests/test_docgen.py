"""Tests of the generators of the report, the slides and the oral-exam sheet (Phase 7)."""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from pptx import Presentation
from pptx.util import Emu

from genomeasm.docgen import deck, md, oral, report
from genomeasm.docgen.facts import Facts, p_text

ROOT = Path(__file__).resolve().parents[1]
needs_chrome = pytest.mark.skipif(md.chrome_path() is None, reason="Chrome is needed for the PDF")


@pytest.fixture(scope="module")
def facts() -> Facts:
    return Facts.load()


@pytest.fixture(scope="module")
def built(facts: Facts, tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    out = tmp_path_factory.mktemp("docs")
    slides = deck.build_slides(facts)
    return {"pptx": deck.build_pptx(slides, out / "deck.pptx"), "dir": out}


# ---- Markdown converter -------------------------------------------------------------------------


def test_markdown_converter_handles_the_constructs_the_documents_use() -> None:
    lines = ["## Title", "", "A **bold** and *italic* word, `a<b` code and [a link](http://x).", ""]
    lines += ["- one", "- two", "", "1. first", "2. second", ""]
    lines += ["| a | b |", "|---|---|", "| 1 | 2 |", "", "```", "x < y", "```", ""]
    lines += ["$$F=" + chr(92) + "sum_{k}x_k$$", "", "Inline $a" + chr(92) + "leq b$ here."]
    text = chr(10).join(lines)
    out = md.to_html(text)
    assert "<h2>Title</h2>" in out and "<strong>bold</strong>" in out and "<em>italic</em>" in out
    assert "<code>a&lt;b</code>" in out and '<a href="http://x">a link</a>' in out
    assert out.count("<li>") == 4 and "<ul>" in out and "<ol>" in out
    assert "<th>a</th>" in out and "<td>2</td>" in out and "x &lt; y" in out
    assert 'class="eq"' in out and 'class="im"' in out and "data:image/svg+xml;base64," in out


def test_formulas_with_latex_short_names_render() -> None:
    for expression in (r"a\le b", r"x\ge y", r"u\to v", r"\mathrm{ov}(i,j)\le\ell_j"):
        assert md.math_svg(expression, 11.0).startswith("data:image/svg+xml;base64,")
    assert md.normalise(r"a\le b\leq c") == r"a\leq b\leq c"


# ---- facts --------------------------------------------------------------------------------------


def test_facts_agree_with_the_committed_results(facts: Facts) -> None:
    assert facts.n_sets == 300 and facts.n_settings == 10 and facts.seeds == 30
    assert facts.lkh_never_beaten and 0 < facts.greedy_equals_best < 1
    assert facts.tabu_equals_greedy > 100 and 0.4 < facts.arm_zero_share < 0.8
    assert facts.verdict_ucb["p_holm"] > 0.05 and facts.worst_live_seconds < 30
    assert 0 < facts.lkh_over < 1 and facts.lkh_max > facts.budget
    assert (
        p_text(1e-6) == "p < 0.0001"
        and p_text(0.5) == "p = 0.50"
        and p_text(0.0123) == "p = 0.0123"
    )


# ---- the report ---------------------------------------------------------------------------------


def test_report_text_has_every_imrad_part_and_no_unfilled_placeholder(facts: Facts) -> None:
    text = report.report_markdown(facts)
    for heading in (
        "## Abstract",
        "## 1. Introduction",
        "## 2. Background",
        "## 3. Problem formalisation",
        "## 4. Methods",
        "## 5. Experimental design",
        "## 6. Results",
        "## 7. Ablation",
        "## 8. Biological quality",
        "## 9. Discussion",
        "## 10. Threats to validity",
        "## 11. Reproducibility statement",
        "## References",
    ):
        assert heading in text
    assert "Proposition 1" in text and "Proposition 2" in text and "not supported" in text
    assert text.count("TODO-VERIFY") >= 10  # every reference is unverified and says so
    for stem in (
        "fig01a_coverage_histogram",
        "fig05_final_fitness",
        "fig07_ablation",
        "fig08_arm_probability",
    ):
        assert (ROOT / "reports" / "figures" / f"{stem}.png").is_file()
    assert f"{facts.n_sets} read sets" in text and p_text(facts.verdict_ucb["p_holm"]) in text


@needs_chrome
def test_report_pdf_is_built_from_the_results(facts: Facts, tmp_path: Path) -> None:
    pdf = report.build_report(tmp_path, facts)
    data = pdf.read_bytes()
    assert data[:5] == b"%PDF-" and len(data) > 500_000
    pages = data.count(b"/Type /Page") - data.count(b"/Type /Pages")
    assert 9 <= pages <= 40
    assert "�" not in (tmp_path / "report.html").read_text(encoding="utf-8")


# ---- the slides ---------------------------------------------------------------------------------


def test_deck_has_at_least_fifteen_slides_with_takeaway_titles_and_notes(
    built: dict[str, Path],
) -> None:
    prs = Presentation(str(built["pptx"]))
    assert len(prs.slides) >= 15
    for number, slide in enumerate(prs.slides, start=1):
        title = slide.shapes[0].text_frame.text
        assert 15 <= len(title) <= 110, (number, title)
        assert len(slide.notes_slide.notes_text_frame.text) >= 80, number
    titles = [s.shapes[0].text_frame.text for s in prs.slides]
    assert len(set(titles)) == len(titles)


def test_nothing_in_the_deck_leaves_the_slide_or_overflows_its_box(
    facts: Facts, built: dict[str, Path]
) -> None:
    prs = Presentation(str(built["pptx"]))
    width, height = prs.slide_width or 0, prs.slide_height or 0
    for number, slide in enumerate(prs.slides, start=1):
        for shape in slide.shapes:
            left, top = shape.left or 0, shape.top or 0
            assert shape.left is not None and shape.width is not None, number
            assert left >= 0 and top >= 0, number
            assert left + (shape.width or 0) <= width + Emu(10), (number, shape.shape_type)
            assert top + (shape.height or 0) <= height + Emu(10), (number, shape.shape_type)
    for number, spec in enumerate(deck.build_slides(facts), start=1):
        if spec.bullets and spec.kind == "content":
            box = deck.layout(spec)["text"]
            size = deck.fit_font(spec.bullets, box[2], box[3])
            assert size >= 14, (number, "text too long for its box")
            assert deck.text_height(spec.bullets, box[2], size) <= box[3] + 1e-6


def test_deck_numbers_come_from_the_results(facts: Facts) -> None:
    text = " ".join(" ".join(s.bullets) + s.title for s in deck.build_slides(facts))
    assert p_text(facts.verdict_ucb["p_holm"]) in text and f"{facts.n_sets} read sets" in text
    assert f"{facts.mean_gap('lkh3'):.2f}%" in text and f"{facts.mean_gap('greedy'):.2f}%" in text


@needs_chrome
def test_deck_pdf_copy_has_one_page_per_slide(facts: Facts, tmp_path: Path) -> None:
    slides = deck.build_slides(facts)
    pdf = deck.build_pdf(slides, tmp_path)
    data = pdf.read_bytes()
    assert data.count(b"/Type /Page") - data.count(b"/Type /Pages") == len(slides)


# ---- the oral-exam sheet ------------------------------------------------------------------------


def test_oral_sheet_has_twenty_questions_and_every_tuned_parameter(facts: Facts) -> None:
    text = oral.oral_markdown(facts)
    assert len(re.findall(r"^\*\*\d+\. ", text, flags=re.M)) == 20
    for section in (
        "## 1. The objective function",
        "## 2. Solution representation",
        "## 3. The data-loading pipeline",
        "## 4. Every hyper-parameter",
        "## 5. Twenty questions",
    ):
        assert section in text
    for entry in facts.tuned.values():
        for key in oral._flat(entry["params"]):
            assert f"| {key} |" in text, key
    for file_name in ("data/raw/SHA256SUMS", "configs/genomes.yaml", "NC_001422.1"):
        assert file_name in text


def test_committed_oral_sheet_is_in_sync_with_the_results(facts: Facts) -> None:
    committed = (ROOT / "docs" / "ORAL_EXAM.md").read_text(encoding="utf-8")
    assert committed == oral.oral_markdown(facts)
