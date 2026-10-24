# Phase 8: QA, rehearsal, freeze

Status: complete; release `v1.0` tagged. Date: 2026-10-24.

## 1. What was checked

| Check | Result |
|---|---|
| Fresh environment | a clean `git clone`, a new virtual environment from `requirements.lock` (internet needed only for this), `make all` |
| Lint, format, types, tests | `make check` (ruff, ruff format, mypy, pytest) inside the clean clone, see section 3 |
| Numbers in the documents | `genomeasm docs verify`: 118 numbers in the slides, 325 in the report, 310 in the oral-exam sheet, 94 in the guide, all explained by the results; none unexplained |
| Statistics recomputed from the raw runs | every pooled test and every per-solver gap recomputed with scipy only: 0 differences |
| Charts | all figures: SVG, PDF, 300-dpi PNG, caption with source and configuration hash, provenance footer, colour-blind-safe colours: 0 problems |
| Demo | `docs/DEMO_SCRIPT.md` (timed, a fallback per step), stored pictures of every page in `reports/screenshots/` |

Machine-readable result: `results/phase8/verification.json` (with the commit it was made at).

## 2. What QA found and fixed

* **Chart colours:** the audit found matplotlib's default colours in three figures (fig01a, figS1, fig05). The default colour cycle is now the Okabe-Ito palette, the boxplot median line is black, the figures were redrawn.
* **False alarm of the number audit:** the page numbers of the references; references are exempt (they are `TODO-VERIFY`).
* **Screenshots:** a first version with a virtual time budget captured half-drawn pages; the final version waits in real time through the browser's debugging protocol.
* **GNU make** was missing on this machine; reinstalled with winget (as chosen in Phase 0).

## 3. `make all` from a clean checkout

See the log summary in section 5 (filled from the last run).

## 4. What is not covered (honest list)

* `make all` does not rerun the Phase 4 validation or the Phase 5 experiments (hours); their raw outputs are committed and their analysis is rebuilt. `make experiments` reruns Phase 5.
* The LKH-3 tests need Ubuntu in WSL and the LKH-3 source in `tools/` (not redistributed); in a checkout without it those tests fail by design (LKH-3 is never replaced).
* The references of the report are unverified (`TODO-VERIFY`); the pptx was not opened in PowerPoint; no projector or other laptop was tested; the cover author is a placeholder (`genomeasm docs build --author`).
* Numbers inside figures are not read by the number audit; they come from the same results by construction and the figures were looked at.

## 5. Freeze

Release `v1.0`. From now on: fixes only. Rebuild with `make all` after any change; commit the rebuilt documents.
