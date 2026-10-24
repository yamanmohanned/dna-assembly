"""Phase 8: the audit of the documents and charts, make all, the demo script and the README."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from genomeasm.docgen import screens, verify
from genomeasm.docgen.facts import Facts

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results" / "phase8" / "verification.json"


@pytest.fixture(scope="module")
def facts() -> Facts:
    return Facts.load()


def test_number_tokens_are_found_and_identifiers_are_not_numbers() -> None:
    text = "LKH-3 solved A12 = 0.52 in 8 s, see D-046 and NC_001422.1; 5,386 letters, p = 0.0013."
    assert set(verify.numbers_in(text)) == {"8", "5,386", "0.0013", "0.52"}
    assert verify.numbers_in("phiX174 and UCB1 and L5-c10-e1.0") == []


def test_explained_numbers_follow_the_digits_shown(facts: Facts) -> None:
    universe = verify.number_universe(facts)
    mean_gap = round(facts.mean_gap("greedy"), 3)
    assert verify.explained(f"{mean_gap:.3f}", universe)
    assert verify.explained(f"{100 * facts.lkh_over:.0f}", universe)
    assert verify.explained("5,386", universe) and verify.explained("300", universe)
    assert not verify.explained("123.456", universe)  # a number that is nowhere in the results
    assert not verify.explained("0.7777", universe)


def test_every_number_of_the_slides_report_oral_sheet_and_guide_is_explained(facts: Facts) -> None:
    documents = verify.check_documents(facts)
    assert set(documents) == {"slides", "report", "oral", "guide"}
    for name, result in documents.items():
        assert result["numbers_checked"] > 50, name
        assert result["unexplained"] == [], (name, result["context"])


def test_the_statistics_tables_are_recomputed_from_the_raw_runs_without_a_difference() -> None:
    assert verify.recompute_statistics() == []


def test_every_figure_meets_the_chart_standards() -> None:
    assert verify.audit_figures() == []


def test_committed_verification_report_says_ok() -> None:
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    assert (
        report["ok"] is True
        and report["recomputation_problems"] == []
        and report["figure_problems"] == []
    )
    assert report["provenance"]["git"]["commit"]


def test_the_makefile_all_target_rebuilds_everything_then_checks() -> None:
    makefile = (ROOT / "Makefile").read_text(encoding="utf-8")
    assert "all: $(PY) reproduce check" in makefile
    for command in (
        "experiments analyze",
        "docs build",
        "docs verify",
        "simulate figures",
        "overlap figures",
        "data verify",
    ):
        assert command in makefile
    assert "experiments tune" in makefile and "all:" in makefile.split("experiments:")[0]


def test_the_demo_script_covers_every_page_and_the_screenshots_exist() -> None:
    script = (ROOT / "docs" / "DEMO_SCRIPT.md").read_text(encoding="utf-8")
    pages = ["Guide", "Learn", "Data", "Simulator", "Overlap Graph", "Solver Lab", "Experiments"]
    pages += ["Assembly Viewer", "Reproducibility"]
    for page in pages:
        assert page in script, page
    for name, _ in screens.PAGES:
        assert (ROOT / "reports" / "screenshots" / f"{name}.png").stat().st_size > 20_000, name
        assert name.split("_", 1)[1] in script.replace(" ", "_").lower() or f"{name}.png" in script
    for section in ("A. Before the talk", "B. The 15-minute run", "C. When something goes wrong"):
        assert section in script


def test_readme_gives_the_quick_start_and_the_honesty_notes() -> None:
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    assert "make app" in readme and "make all" in readme and "TODO-VERIFY" in readme
