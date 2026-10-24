"""A small Markdown-to-HTML converter with formulas, and Chrome's print-to-PDF.

No pandoc or LaTeX is installed (D-007). Formulas are typeset by matplotlib's mathtext as SVG
images (``$...$`` inline, a line ``$$...$$`` as a display formula); the page is printed to PDF by
headless Chrome. Supported Markdown: headings, paragraphs, lists, tables, fenced code, raw HTML
lines, **bold**, *italic*, `code`, [links](url).
"""

from __future__ import annotations

import base64
import html
import io
import re
import subprocess
import tempfile
from functools import lru_cache
from pathlib import Path

from matplotlib.figure import Figure

CHROME_CANDIDATES = (
    r"C:\Program Files\Google\Chrome\Application\chrome.exe",
    r"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe",
    r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
)


def chrome_path() -> Path | None:
    for candidate in CHROME_CANDIDATES:
        if Path(candidate).is_file():
            return Path(candidate)
    return None


def normalise(expression: str) -> str:
    """LaTeX short names that matplotlib's mathtext spells differently."""
    for short, long in ((r"\le", r"\leq"), (r"\ge", r"\geq"), (r"\to", r"\rightarrow")):
        pattern = re.compile(re.escape(short) + r"(?![a-zA-Z])")
        expression = pattern.sub(long.replace(chr(92), chr(92) * 2), expression)
    return expression


@lru_cache(maxsize=512)
def math_svg(expression: str, size: float) -> str:
    """A formula as a base64 SVG data URI (text turned into paths, transparent background)."""
    expression = normalise(expression)
    fig = Figure(figsize=(0.1, 0.1))
    fig.text(0, 0, f"${expression}$", fontsize=size)
    buffer = io.BytesIO()
    fig.savefig(
        buffer, format="svg", bbox_inches="tight", pad_inches=0.03, transparent=True,
        metadata={"Date": None},
    )  # fmt: skip
    return "data:image/svg+xml;base64," + base64.b64encode(buffer.getvalue()).decode("ascii")


def math_html(expression: str, display: bool = False) -> str:
    size = 13.0 if display else 11.0
    label = html.escape(expression, quote=True)
    if display:
        return f'<div class="eq"><img alt="{label}" src="{math_svg(expression, size)}"></div>'
    return f'<img class="im" alt="{label}" src="{math_svg(expression, size)}">'


def inline(text: str) -> str:
    """Bold, italic, code, links and $formulas$ inside a line."""
    keep: list[str] = []

    def stash(value: str) -> str:
        keep.append(value)
        return f"\x00{len(keep) - 1}\x00"

    text = re.sub(r"`([^`]+)`", lambda m: stash(f"<code>{html.escape(m.group(1))}</code>"), text)
    text = re.sub(r"\$([^$]+)\$", lambda m: stash(math_html(m.group(1))), text)
    text = html.escape(text, quote=False)
    text = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", text)
    text = re.sub(r"(?<![\w*])\*(?!\s)(.+?)(?<!\s)\*(?![\w*])", r"<em>\1</em>", text)
    text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r'<a href="\2">\1</a>', text)
    return re.sub("\x00(\\d+)\x00", lambda m: keep[int(m.group(1))], text)


def to_html(markdown: str) -> str:
    out: list[str] = []
    lines = markdown.split(chr(10))
    i = 0
    paragraph: list[str] = []

    def flush() -> None:
        if paragraph:
            out.append("<p>" + inline(" ".join(paragraph)) + "</p>")
            paragraph.clear()

    while i < len(lines):
        line = lines[i]
        stripped = line.strip()
        if stripped.startswith("```"):
            flush()
            block = []
            i += 1
            while i < len(lines) and not lines[i].strip().startswith("```"):
                block.append(lines[i])
                i += 1
            out.append("<pre><code>" + html.escape(chr(10).join(block)) + "</code></pre>")
        elif not stripped:
            flush()
        elif re.fullmatch(r"\$\$.+\$\$", stripped):
            flush()
            out.append(math_html(stripped[2:-2], display=True))
        elif stripped.startswith("<"):
            flush()
            out.append(line)
        elif m := re.match(r"(#{1,4}) (.*)", stripped):
            flush()
            level = len(m.group(1))
            out.append(f"<h{level}>{inline(m.group(2))}</h{level}>")
        elif stripped.startswith("|"):
            flush()
            rows = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                rows.append([c.strip() for c in lines[i].strip().strip("|").split("|")])
                i += 1
            i -= 1
            head, body = rows[0], [r for r in rows[2:]]
            out.append(
                "<table><thead><tr>"
                + "".join(f"<th>{inline(c)}</th>" for c in head)
                + "</tr></thead><tbody>"
            )
            for r in body:
                out.append("<tr>" + "".join(f"<td>{inline(c)}</td>" for c in r) + "</tr>")
            out.append("</tbody></table>")
        elif re.match(r"([-*]|\d+\.) ", stripped):
            flush()
            ordered = stripped[0].isdigit()
            tag = "ol" if ordered else "ul"
            items: list[str] = []
            while i < len(lines) and re.match(r"\s*([-*]|\d+\.) ", lines[i]):
                items.append(re.sub(r"^\s*([-*]|\d+\.) ", "", lines[i]))
                i += 1
                while i < len(lines) and lines[i].startswith("  ") and lines[i].strip():
                    items[-1] += " " + lines[i].strip()
                    i += 1
            i -= 1
            out.append(f"<{tag}>" + "".join(f"<li>{inline(x)}</li>" for x in items) + f"</{tag}>")
        else:
            paragraph.append(stripped)
        i += 1
    flush()
    return chr(10).join(out)


def figure_html(path: Path, caption: str, width: str = "100%") -> str:
    """A figure as an embedded PNG with its caption."""
    data = base64.b64encode(path.read_bytes()).decode("ascii")
    return (
        f'<figure><img style="width:{width}" alt="{html.escape(caption[:90], quote=True)}" '
        f'src="data:image/png;base64,{data}"><figcaption>{inline(caption)}</figcaption></figure>'
    )


def print_pdf(html_file: Path, pdf_file: Path, timeout: int = 180) -> None:
    """Print an HTML file to PDF with headless Chrome (raises RuntimeError if that fails)."""
    chrome = chrome_path()
    if chrome is None:
        raise RuntimeError("Chrome (or Edge) is needed to make the PDF and was not found")
    pdf_file.parent.mkdir(parents=True, exist_ok=True)
    if pdf_file.exists():
        pdf_file.unlink()
    with (
        tempfile.TemporaryDirectory() as profile
    ):  # a private profile: the user's Chrome may be open
        command = [
            str(chrome), "--headless=new", "--disable-gpu", "--no-pdf-header-footer", "--no-sandbox",
            f"--user-data-dir={profile}", f"--print-to-pdf={pdf_file.resolve()}", html_file.resolve().as_uri(),
        ]  # fmt: skip
        done = subprocess.run(command, check=False, capture_output=True, text=True, timeout=timeout)
    if not pdf_file.is_file() or pdf_file.stat().st_size < 1000:
        raise RuntimeError(f"Chrome did not write {pdf_file}: {done.stderr[-300:]}")
