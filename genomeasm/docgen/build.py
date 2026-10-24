"""One command for all three documents: ``genomeasm docs build``."""

from __future__ import annotations

from pathlib import Path

from genomeasm.appkit.paths import DOCS, ROOT
from genomeasm.docgen import deck, oral, report
from genomeasm.docgen.facts import Facts

TARGETS = ("report", "slides", "oral")


def build_all(
    only: str | None = None, root: Path = ROOT, author: str = ""
) -> dict[str, list[Path]]:
    """Build the chosen documents (all by default) from results/; returns the files written."""
    facts = Facts.load()
    written: dict[str, list[Path]] = {}
    if only in (None, "report"):
        written["report"] = [
            report.build_report(
                root / "reports" / "report", facts, author or "[name to be filled in]"
            )
        ]
    if only in (None, "slides"):
        written["slides"] = list(deck.build_deck(root / "reports" / "slides", facts, author))
    if only in (None, "oral"):
        written["oral"] = [oral.build_oral(DOCS / "ORAL_EXAM.md", facts)]
    return written
