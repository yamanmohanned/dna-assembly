"""DNA Studio: the self-contained HTML documents, the lecture and the data the pages show.

The scripts in ``genomeasm/appkit/web`` draw in a browser and are checked there by hand (see
docs/phases/PHASE_9.md); here they are parsed with Node.js when it is installed, and everything
that is decided in Python (which slides, which numbers, which sources) is tested directly.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
from pathlib import Path

import pytest

from genomeasm.appkit import biology, lecture, studio
from genomeasm.docgen.facts import Facts
from genomeasm.molbio import library, sequence

EXTERNAL = re.compile(r"""(?:src|href)\s*=\s*["'](?:https?:)?//|@import|url\(\s*["']?https?:""")


def _data(html: str) -> dict:
    match = re.search(r'<script type="application/json" id="ds-data">(.*?)</script>', html, re.S)
    assert match, "the document carries its data as JSON"
    return json.loads(match.group(1).replace("<\\/", "</"))


@pytest.fixture(scope="module")
def deck() -> dict:
    return lecture.build()


@pytest.fixture(scope="module")
def deck_html() -> str:
    return lecture.html()


@pytest.mark.parametrize("name", studio.SCRIPTS)
def test_scripts_parse(name: str) -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node.js is not installed: the scripts are only checked in the browser")
    result = subprocess.run(
        [node, "--check", str(studio.WEB / name)], capture_output=True, text=True, check=False
    )
    assert result.returncode == 0, result.stderr


def test_scripts_load_nothing_from_the_network() -> None:
    for path in sorted(studio.WEB.iterdir()):
        text = path.read_text(encoding="utf-8")
        assert not re.search(r"https?://(?!www\.w3\.org)", text), path.name
        assert "fetch(" not in text and "XMLHttpRequest" not in text, path.name


def test_viewer_document_is_self_contained_and_carries_its_scene() -> None:
    scene = biology.scenes()["atomic"]
    html = studio.viewer_html(scene, height=500, colour="strand", spin=True)
    assert html.startswith("<!DOCTYPE html>") and 'lang="en" dir="ltr" translate="no"' in html
    assert not EXTERNAL.search(html) and "<link" not in html
    data = _data(html)
    assert data["scene"] == json.loads(json.dumps(scene))
    assert data["options"] == {"colour": "strand", "spin": True}
    assert "DNAStudio.MoleculeView" in html and "height:500px" in html


def test_json_cannot_close_the_script_element() -> None:
    text = studio.json_for_html({"note": "</script><script>alert(1)</script>", "c": "<!--"})
    assert "</script" not in text and "<!--" not in text
    assert json.loads(text.replace("<\\/", "</").replace("<\\!--", "<!--"))["c"] == "<!--"


def test_process_documents_exist_for_the_three_processes() -> None:
    for kind, scenario in biology.scenarios().items():
        html = studio.process_html(scenario)
        assert not EXTERNAL.search(html), kind
        assert _data(html)["scenario"]["kind"] == kind
        assert "DNAStudio.ProcessPlayer" in html


def test_lecture_covers_every_part_in_order(deck: dict) -> None:
    slides = deck["slides"]
    assert 24 <= len(slides) <= 32
    assert slides[0]["layout"] == "title" and slides[1]["layout"] == "outline"
    assert [s["layout"] for s in slides[-3:]] == ["summary", "references", "references"]
    sections = [s["section"] for s in slides]
    order = [name for name, _ in lecture.SECTIONS]
    assert [s for i, s in enumerate(sections) if s in order and s not in sections[:i]] == order
    for item in slides[1]["items"]:  # the outline jumps to the first slide of each part
        assert slides[item["slide"]]["section"] == item["name"]
        assert item["slide"] == sections.index(item["name"])
    for slide in slides:
        assert slide["title"] and slide["kicker"]


def test_every_figure_of_the_lecture_exists(deck: dict) -> None:
    kinds = set()
    for slide in deck["slides"]:
        figure = slide.get("figure")
        if not figure:
            continue
        kinds.add(figure["type"])
        if figure["type"] == "viewer":
            assert figure["scene"] in deck["scenes"]
        elif figure["type"] == "viewers":
            assert all(item["scene"] in deck["scenes"] for item in figure["items"])
        elif figure["type"] == "process":
            assert figure["scenario"] in deck["scenarios"]
        elif figure["type"] == "image":
            assert figure["src"].startswith("data:image/") and figure["alt"]
        else:
            assert figure["html"].lstrip().startswith(("<svg", "<table", "<style"))
    assert kinds == {"viewer", "viewers", "process", "image", "html"}
    used = {
        s["figure"]["scenario"]
        for s in deck["slides"]
        if s.get("figure", {}).get("type") == "process"
    }
    assert used == {"replication", "transcription", "translation"}


def test_lecture_numbers_come_from_the_data(deck: dict) -> None:
    text = json.dumps(deck["slides"], ensure_ascii=False)
    measured = library.measured_helix()
    assert f"{measured['rise_mean']:.2f} Å" in text and f"{measured['bp_per_turn']:.1f}" in text
    facts = Facts.load()
    assert f"<b>{facts.n_sets}</b> read sets" in text
    assert f"{facts.budget:.0f} seconds" in text
    # the honest verdict of Phase 5 stays visible, and the post-hoc result is labelled as such
    assert "not supported" in text and "Post hoc, exploratory only" in text
    assert f"A12 = {float(facts.verdict_ucb['a12']):.2f}" in text
    protein = sequence.translate(biology.gene_sequence("J"))
    assert f"{len(protein)} amino acids" in text
    for row in biology.mutation_rows():
        assert row["effect"] in text


def test_lecture_lists_every_reference_and_cites_only_known_sources(deck: dict) -> None:
    pages = [s["html"] for s in deck["slides"] if s["layout"] == "references"]
    listed = "".join(pages)
    for ref in library.references().values():
        assert ref["title"].replace("&", "&amp;").replace("'", "&#x27;") in listed, ref["key"]
    assert listed.count("TODO-VERIFY") == sum(
        1 for r in library.references().values() if r["check"] == "not-checked"
    )
    known = {library.cite(key) for key in library.references()}
    for slide in deck["slides"]:
        for source in slide.get("sources", []):
            assert source in known, source
        assert len(slide.get("sources", [])) == len(set(slide.get("sources", [])))


def test_lecture_document_is_one_file_that_needs_no_network(deck_html: str, tmp_path: Path) -> None:
    assert not EXTERNAL.search(deck_html) and "<link" not in deck_html
    assert "DNAStudio.Deck" in deck_html and len(deck_html) < 1_500_000
    data = _data(deck_html)
    assert data["options"] == {"embedded": False} and len(data["lecture"]["slides"]) >= 24
    path = lecture.write(tmp_path / "out" / "lecture.html", author="A. Student")
    written = path.read_text(encoding="utf-8")
    assert "A. Student" in written and "\r\n" not in written
    assert _data(lecture.html(embedded=True))["options"] == {"embedded": True}


def test_helix_table_and_genome_chart_rows_are_traceable() -> None:
    rows = biology.helix_rows()
    assert len(rows) == 6 and all(r["measured"] and r["published"] for r in rows)
    assert "Watson & Crick 1953" in rows[0]["published"] and "Wang 1979" in rows[2]["published"]
    sizes = biology.genome_rows()
    assert [size for _, size, _ in sizes] == sorted(size for _, size, _ in sizes)
    assert sizes[0][1] == len(biology.genome(biology.PHIX)) == 5386
    assert sizes[1][1] == len(biology.genome(biology.LAMBDA)) == 48502
    assert sizes[-1][1] == int(library.facts()["human"]["number"])


def test_mutation_examples_cover_the_four_textbook_outcomes() -> None:
    rows = biology.mutation_rows()
    effects = [r["effect"] for r in rows]
    assert {"silent", "missense", "nonsense", "frameshift"} <= set(effects)
    silent = rows[effects.index("silent")]
    assert silent["protein_after"] == silent["protein_before"]
    nonsense = rows[effects.index("nonsense")]
    assert nonsense["protein_before"].startswith(nonsense["protein_after"])
    assert len(nonsense["protein_after"]) < len(nonsense["protein_before"])


def test_shotgun_example_is_a_real_simulation_with_the_stated_coverage() -> None:
    shot = biology.shotgun_example()
    assert shot["n"] == round(shot["coverage"] * shot["length"] / shot["read_length"])
    assert len(shot["starts"]) == len(shot["lengths"]) == shot["n"]
    assert all(
        s >= 0 and s + n <= shot["length"]
        for s, n in zip(shot["starts"], shot["lengths"], strict=True)
    )
    summary = biology.overlap_summary()
    assert summary["n_genes"] == 11 and summary["sum_of_gene_lengths"] > summary["genome_length"]
    assert summary["coding_letters"] <= summary["genome_length"]


def test_assets_stamp_follows_the_newest_file() -> None:
    newest = max(p.stat().st_mtime_ns for p in studio.WEB.iterdir() if p.is_file())
    assert studio.assets_stamp() == str(newest)
    assert set(studio.SCRIPTS + studio.STYLES) <= {p.name for p in studio.WEB.iterdir()}
