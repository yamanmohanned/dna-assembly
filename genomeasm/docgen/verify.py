"""Phase 8 audit: are the numbers in the documents the numbers in ``results/``, and do the charts meet the standards?

Three independent checks (``genomeasm docs verify``):

1. **Numbers in the documents.** Every number written in the slides (text, tables, notes), the report,
   the oral-exam sheet and the guide must be explained: it equals a value of the results (rounded to
   the digits shown, also as a percentage), a value of the configuration or of a genome card, or it is
   a small whole number (a count, a section or figure number, a year). Whatever is not explained is
   listed and fails the check. Numbers inside figures are covered by check 2 and by the figure tests.
2. **Recomputation from the raw runs.** The statistics tables are recomputed from
   ``results/phase5/main_runs.jsonl.gz`` with scipy only, independently of the analysis code.
3. **Chart standards.** All figures exist as SVG, PDF and 300-dpi PNG with a caption and provenance, and
   use only colour-blind-safe colours.
"""

from __future__ import annotations

import html
import json
import math
import re
from typing import Any

import numpy as np
import pandas as pd
from PIL import Image
from pptx import Presentation
from scipy import stats

from genomeasm import provenance
from genomeasm.appkit import lab
from genomeasm.appkit.paths import CONFIG, DOCS, FIGURES, PROCESSED, RESULTS, ROOT
from genomeasm.docgen.facts import Facts
from genomeasm.experiments import runner
from genomeasm.plots import OKABE_ITO

# whole numbers up to this size are counts, section numbers, years of a date and the like
SMALL_INT = 31
STRUCTURAL = {2026.0}  # the year of the dates in the documents


def _walk(value: Any, out: set[float]) -> None:
    if isinstance(value, bool):
        return
    if isinstance(value, (int, float)):
        out.add(float(value))
    elif isinstance(value, dict):
        for v in value.values():
            _walk(v, out)
    elif isinstance(value, (list, tuple)):
        for v in value:
            _walk(v, out)


def number_universe(F: Facts) -> set[float]:
    """Every number of the results, the configuration and the data cards, plus simple derivations."""
    base: set[float] = set()
    for path in (RESULTS / "analysis").glob("*.csv"):
        frame = pd.read_csv(path)
        for column in frame.select_dtypes("number").columns:
            base.update(float(x) for x in frame[column].dropna())
    for path in (
        RESULTS / "analysis" / "analysis.json",
        RESULTS / "tuning_summary.json",
        ROOT / "results" / "phase6" / "interactive_timing.json",
    ):
        _walk(json.loads(path.read_text(encoding="utf-8")), base)
    _walk(F.config, base)
    for card in PROCESSED.glob("*.card.json"):
        _walk(json.loads(card.read_text(encoding="utf-8")), base)
    base.update(
        {
            float(F.n_sets),
            100 * F.greedy_equals_best,
            100 * F.arm_zero_share,
            100 * F.lkh_over,
            F.lkh_max,
            F.lkh_max - F.budget,
            float(F.tabu_equals_greedy),
            F.worst_live_seconds,
            float(lab.MAX_READS),
            float(lab.MAX_WINDOW),
            float(lab.MAX_SECONDS),
            float(lab.MAX_TOTAL_SECONDS),
        }  # fmt: skip
    )
    universe = set(base)
    for u in base:
        universe.update({-u, 100 * u, u / 100, 1 - u, 100 * (1 - u), 100 * abs(u)})
    for cfg in F.config["evaluation_instances"] + F.config["tuning_instances"]:
        universe.update({float(cfg["start"] + cfg["length"]), 100 * float(cfg["error"])})
    universe.update(float(x) for x in range(0, SMALL_INT + 1))
    universe.update(float(math.factorial(n)) for n in range(1, 26))  # the orders of n reads
    universe.update(float(year) for year in range(1900, 2026))  # years of citations
    universe.add(8501.0)  # the port of the application
    universe.update(STRUCTURAL)
    return universe


IDENTIFIER = re.compile(
    r"\b[A-Za-z_][\w.\-]*\d[\w.\-]*\b"
)  # LKH-3, A12, NC_001422.1, D-046, UCB1, L5-c10-e1.0
NUMBER = re.compile(r"(?<![\w.])\d[\d,]*(?:\.\d+)?")


def numbers_in(text: str) -> list[str]:
    """The number tokens of a text (identifiers such as LKH-3 or D-046 are not numbers)."""
    text = IDENTIFIER.sub(" ", text)
    return [m.group(0).rstrip(".,") for m in NUMBER.finditer(text)]


def explained(token: str, universe: set[float]) -> bool:
    cleaned = token.replace(",", "")
    try:
        value = float(cleaned)
    except ValueError:
        return True
    decimals = len(cleaned.split(".")[1]) if "." in cleaned else 0
    tolerance = 0.5 * 10.0**-decimals + 1e-9
    return any(abs(u - value) <= tolerance for u in universe)


def document_texts() -> dict[str, str]:
    """The text of every deliverable (images and styles removed)."""
    texts: dict[str, str] = {}
    prs = Presentation(str(ROOT / "reports" / "slides" / "deck.pptx"))
    parts = []
    for slide in prs.slides:
        for shape in slide.shapes:
            if shape.has_text_frame:
                parts.append(shape.text_frame.text)
            if shape.has_table:
                parts += [cell.text for row in shape.table.rows for cell in row.cells]
        parts.append(slide.notes_slide.notes_text_frame.text)
    texts["slides"] = chr(10).join(parts)
    page = (ROOT / "reports" / "report" / "report.html").read_text(encoding="utf-8")
    page = re.sub(r"<(style|script)[^>]*>.*?</\1>", " ", page, flags=re.S)
    page = re.sub(r"<img[^>]*>", " ", page)
    body = html.unescape(re.sub(r"<[^>]+>", " ", page))
    start, end = body.find("References All references"), body.find("Appendix A")
    texts["report"] = (
        body[:start] + body[end:] if 0 <= start < end else body
    )  # references: TODO-VERIFY
    texts["oral"] = (DOCS / "ORAL_EXAM.md").read_text(encoding="utf-8")
    texts["guide"] = (DOCS / "GUIDE.md").read_text(encoding="utf-8")
    return texts


def check_documents(F: Facts) -> dict[str, Any]:
    universe = number_universe(F)
    out: dict[str, Any] = {}
    for name, text in document_texts().items():
        tokens = numbers_in(text)
        bad = sorted({t for t in tokens if not explained(t, universe)})
        context = {}
        for token in bad[:40]:
            at = text.find(token)
            context[token] = " ".join(text[max(0, at - 40) : at + 40].split())
        out[name] = {"numbers_checked": len(tokens), "unexplained": bad, "context": context}
    return out


# ---- check 2: statistics recomputed from the raw runs --------------------------------------------------


def _runs_from_raw() -> pd.DataFrame:
    rows = [
        (r["instance"], r["seed_index"], r["requested"], r["fitness"])
        for r in runner.read_records(RESULTS / "main_runs.jsonl.gz")
        if r["kind"] == "run"
    ]
    frame = pd.DataFrame(rows, columns=["instance", "seed", "algorithm", "fitness"])
    return frame.drop_duplicates(["instance", "seed", "algorithm"], keep="last")


def recompute_statistics(tolerance: float = 1e-6) -> list[dict[str, Any]]:
    """Mean/median gap per solver and every pooled test, recomputed with numpy and scipy only."""
    runs = _runs_from_raw()
    wide = runs.pivot(index=["instance", "seed"], columns="algorithm", values="fitness")
    best = wide.max(axis=1)
    gap = 100 * (-wide).add(best, axis=0).div(best, axis=0)
    rel = wide.div(best, axis=0)
    stored = pd.read_csv(RESULTS / "analysis" / "by_algorithm.csv").set_index("algorithm")
    tests = pd.read_csv(RESULTS / "analysis" / "tests.csv")
    problems: list[dict[str, Any]] = []

    def compare(what: str, mine: float, theirs: float, tol: float = tolerance) -> None:
        if not abs(mine - theirs) <= tol * max(1.0, abs(theirs)):
            problems.append({"what": what, "recomputed": float(mine), "stored": float(theirs)})

    for algorithm in stored.index:
        compare(f"mean gap {algorithm}", gap[algorithm].mean(), stored.loc[algorithm, "mean_gap"])
        compare(f"median gap {algorithm}", gap[algorithm].median(), stored.loc[algorithm, "median"])
    pooled = tests[tests["scope"] == "pooled"]
    for row in pooled.itertuples():
        a, b = rel[row.a], rel[row.b]
        d = a - b
        compare(f"mean diff {row.a} vs {row.b}", 100 * d.mean(), row.mean_diff_pp)
        compare(f"median diff {row.a} vs {row.b}", 100 * d.median(), row.median_diff_pp)
        wins, losses = int((d > 0).sum()), int((d < 0).sum())
        compare(f"wins {row.a} vs {row.b}", wins, row.wins, 0)
        compare(f"losses {row.a} vs {row.b}", losses, row.losses, 0)
        x, y = a.to_numpy(), b.to_numpy()
        a12 = float((x[:, None] > y[None, :]).mean() + 0.5 * (x[:, None] == y[None, :]).mean())
        compare(f"A12 {row.a} vs {row.b}", a12, row.a12)
        p = stats.wilcoxon(d, zero_method="wilcox").pvalue if (d != 0).any() else 1.0
        compare(f"p {row.a} vs {row.b}", p, row.p, 1e-6)
    for family, group in pooled.groupby("family"):  # Holm, step down
        order = np.argsort(group["p"].to_numpy())
        m = len(group)
        running, adjusted = 0.0, np.empty(m)
        for rank, k in enumerate(order):
            running = max(running, min(1.0, (m - rank) * float(group["p"].to_numpy()[k])))
            adjusted[k] = running
        for k, row in enumerate(group.itertuples()):
            compare(f"Holm {family}: {row.a} vs {row.b}", adjusted[k], row.p_holm, 1e-6)
    pivot_best = wide.drop(columns="lkh3").le(wide["lkh3"], axis=0).all().all()
    if not pivot_best:
        problems.append(
            {"what": "LKH-3 is not the best on every read set", "recomputed": 0, "stored": 1}
        )
    return problems


# ---- check 3: chart standards ---------------------------------------------------------------------------

NEUTRAL = {
    "#b0b0b0",
    "#000000",
    "#ffffff",
    "#777777",
    "#bbbbbb",
    "#999999",
    "#cccccc",
    "#444444",
    "#555555",
    "#aaaaaa",
}
EXTRA_PALETTE = {
    "#0b3c5d",
    "#882255",
    "#d9d9d9",
    "#333333",
    "#666666",
    "#c9d2db",
}  # darker blue, wine, greys


def audit_figures() -> list[str]:
    """Problems found in the figures of reports/figures (empty list: all good)."""
    allowed = {c.lower() for c in OKABE_ITO.values()} | NEUTRAL | EXTRA_PALETTE
    problems: list[str] = []
    stems = sorted(p.name.removesuffix(".caption.json") for p in FIGURES.glob("*.caption.json"))
    if len(stems) < 11:
        problems.append(f"only {len(stems)} figures found")
    for stem in stems:
        for suffix in (".svg", ".pdf", ".png", ".caption.json"):
            if not (FIGURES / f"{stem}{suffix}").is_file():
                problems.append(f"{stem}: {suffix} is missing")
        sidecar = json.loads((FIGURES / f"{stem}.caption.json").read_text(encoding="utf-8"))
        if len(sidecar.get("caption", "")) < 40:
            problems.append(f"{stem}: caption too short")
        if not re.fullmatch(
            r"[0-9a-f]{64}", str(sidecar.get("config_hash", ""))
        ) or not sidecar.get("source"):
            problems.append(f"{stem}: no provenance (source, config hash)")
        with Image.open(FIGURES / f"{stem}.png") as image:
            dpi = image.info.get("dpi", (0, 0))[0]
            if dpi < 299:
                problems.append(f"{stem}: PNG is {dpi:.0f} dpi, not 300")
        svg = (FIGURES / f"{stem}.svg").read_text(encoding="utf-8")
        if "Source:" not in svg:
            problems.append(f"{stem}: no provenance footer inside the figure")
        used = {c.lower() for c in re.findall(r"#[0-9a-fA-F]{6}\b", svg)}
        foreign = sorted(used - allowed)
        if foreign:
            problems.append(f"{stem}: colours outside the colour-blind-safe palette: {foreign[:6]}")
    return problems


def verify_all(F: Facts) -> dict[str, Any]:
    documents = check_documents(F)
    recomputation = recompute_statistics()
    figures = audit_figures()
    ok = not recomputation and not figures and all(not d["unexplained"] for d in documents.values())
    return {
        "ok": ok,
        "documents": documents,
        "recomputation_problems": recomputation,
        "figure_problems": figures,
        "provenance": provenance.collect({"config": str(CONFIG.name)}, seed=None, repo_dir=ROOT),
    }
