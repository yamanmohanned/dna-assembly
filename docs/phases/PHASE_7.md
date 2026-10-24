# Phase 7: report, slides, oral-exam sheet

Status: complete; waiting for your "continue" (Phase 8: QA, rehearsal, freeze). Date: 2026-10-24.

## 1. What was built

One command builds all three from `results/` (about 20 s): `genomeasm docs build` (or `make docs`; `--only report|slides|oral`).

| Document | File | Size |
|---|---|---|
| PDF report (IMRaD) | `reports/report/report.pdf` | 12 pages: abstract, introduction, background in plain language, problem formalisation with two proofs, methods, experimental design, results, ablation, biological quality, discussion, threats to validity, reproducibility statement, references, appendix |
| Slides | `reports/slides/deck.pptx` | 20 slides; every title states its takeaway; speaker notes in plain language |
| Slides as PDF (fallback) | `reports/slides/deck.pdf` | one page per slide, same content |
| Oral-exam sheet | `docs/ORAL_EXAM.md` | objective, representation, data pipeline with file names, every hyper-parameter with its tuned value and reason, 20 questions with model answers |

Every number in the three documents is read from `results/` through `genomeasm/docgen/facts.py`; the oral-exam sheet and the guide's glossary are checked against their sources by tests (`tests/test_docgen.py`, 11 tests).

## 2. How to use them

* **Report:** open the PDF. The author's name on the cover is a placeholder ("[name to be filled in]"): run `genomeasm docs build --author "Your Name"` (it goes on the cover and the title slide).
* **Slides:** open the pptx in PowerPoint; if no PowerPoint is available, show the PDF copy full screen. Speaker notes are in the pptx (and the Guide, section 8, has the 15-minute script that goes with the live demo).
* **Oral exam:** read `docs/GUIDE.md` first, then `docs/ORAL_EXAM.md`.

## 3. What is not verified (be honest about it)

* **References:** all ten are written from memory and are marked `TODO-VERIFY` in the report. Verify them (title, year, pages) before submitting.
* **The pptx has not been opened in PowerPoint** (none is installed here; Q3 is open). Its layout is checked by an estimate of text height per box, by a test that no shape leaves the slide, and by looking at the HTML/PDF copy of the same content (slides 4 and 14 were looked at in detail). Please open it once and tell me what looks wrong.
* **The report PDF was looked at as HTML screenshots** (sections 3, 4, 6, 7) rather than as the PDF pages (no PDF renderer here); page numbers are added by Chrome's page-margin box.
* The claims in the report text (hypotheses on why the bandit did not win) are marked as hypotheses.

## 4. Decisions

D-052 (toolchain; closes D-007). Questions still open: Q3 (program that opens the pptx), Q5 (is this the presenting laptop).
