"""DNA Studio: self-contained HTML for the 3D viewer, the process animations and the lecture.

Every document made here is one HTML string with its style sheet, scripts and data inside: no
file is fetched from anywhere, so it works offline, inside the application (``st.iframe``) and
as a file opened by double-click. The scripts live in ``genomeasm/appkit/web``; they only draw.
Sequences, geometry, numbers and texts are computed in Python and passed as JSON.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

WEB = Path(__file__).resolve().parent / "web"
SCRIPTS = ("engine3d.js", "molecule.js", "process.js", "deck.js")
STYLES = ("studio.css", "deck.css")


def asset(name: str) -> str:
    return (WEB / name).read_text(encoding="utf-8")


def json_for_html(data: Any) -> str:
    """JSON that is safe inside a <script> element (no ``</script>`` can appear in it)."""
    text = json.dumps(data, ensure_ascii=False, separators=(",", ":"))
    return text.replace("</", "<\\/").replace("<!--", "<\\!--")


def document(
    body: str,
    *,
    data: Any,
    boot: str,
    scripts: tuple[str, ...],
    styles: tuple[str, ...] = ("studio.css",),
    title: str = "DNA Studio",
    extra_css: str = "",
) -> str:
    """One self-contained HTML document: styles, body, data as JSON, scripts, start-up code."""
    css = "\n".join(asset(name) for name in styles) + extra_css
    js = "\n".join(asset(name) for name in scripts)
    return (
        '<!DOCTYPE html>\n<html lang="en" dir="ltr" translate="no">\n<head>\n'
        '<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        '<meta name="google" content="notranslate">\n'
        f"<title>{title}</title>\n<style>\n{css}\n</style>\n</head>\n<body>\n{body}\n"
        f'<script type="application/json" id="ds-data">{json_for_html(data)}</script>\n'
        f"<script>\n{js}\n</script>\n"
        "<script>\n(function () {\n"
        '  var DATA = JSON.parse(document.getElementById("ds-data").textContent);\n'
        f"{boot}\n}})();\n</script>\n</body>\n</html>\n"
    )


def viewer_html(scene: dict[str, Any], *, height: int = 560, **options: Any) -> str:
    """The interactive 3D viewer of one scene (``genomeasm.molbio.model``).

    Options are passed to ``DNAStudio.MoleculeView``: style, colour, hbonds, labels, spin,
    atomLabels, distanceLabels, note, rotation (quaternion), zoom, hint, selectable.
    """
    body = f'<div id="ds-root" style="height:{int(height)}px"></div>'
    boot = (
        '  window.dsView = new DNAStudio.MoleculeView(document.getElementById("ds-root"), '
        "DATA.scene, DATA.options);"
    )
    return document(
        body,
        data={"scene": scene, "options": options},
        boot=boot,
        scripts=("engine3d.js", "molecule.js"),
        title="DNA Studio: 3D model",
    )


def process_html(scenario: dict[str, Any], *, height: int = 620) -> str:
    """An animated, step-by-step biological process (``genomeasm.molbio.processes``)."""
    body = f'<div id="ds-root" style="height:{int(height)}px"></div>'
    boot = (
        '  window.dsProcess = new DNAStudio.ProcessPlayer(document.getElementById("ds-root"), '
        "DATA.scenario, DATA.options);"
    )
    return document(
        body,
        data={"scenario": scenario, "options": {}},
        boot=boot,
        scripts=("engine3d.js", "molecule.js", "process.js"),
        title="DNA Studio: " + str(scenario.get("title", "process")),
    )


def lecture_html(lecture: dict[str, Any], *, embedded: bool = False) -> str:
    """The whole lecture as one self-contained slide deck (``genomeasm.appkit.lecture``)."""
    body = '<div id="deck" class="deck"></div>'
    boot = (
        '  window.dsDeck = new DNAStudio.Deck(document.getElementById("deck"), DATA.lecture, '
        "DATA.options);"
    )
    return document(
        body,
        data={"lecture": lecture, "options": {"embedded": embedded}},
        boot=boot,
        scripts=SCRIPTS,
        styles=STYLES,
        title=str(lecture.get("title", "Lecture")),
    )


def fragment_html(body: str, css: str = "") -> str:
    """A static piece of HTML (a table, a drawing) as its own small document."""
    return document(body, data={}, boot="", scripts=(), extra_css=css, title="DNA Studio")


def assets_stamp() -> str:
    """Changes whenever a script or style sheet of DNA Studio is edited (a key for caches)."""
    newest = max(path.stat().st_mtime_ns for path in WEB.iterdir() if path.is_file())
    return str(newest)
