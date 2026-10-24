# DNA fragment assembly as an asymmetric TSP

Greedy, Tabu Search, LKH-3 and an ALNS whose repair tool is chosen by a multi-armed bandit, compared on
real viral genomes under equal time budgets. **Result in one line:** bandit selection did *not* beat random
selection (pre-registered claim not supported); LKH-3 was never beaten; ALNS clearly beats greedy and tabu search.

## Start

```
make venv          # once: .venv with the pinned packages (needs internet only here)
make app           # the application, offline: http://localhost:8501  (start at the page "Guide")
make all           # rebuild every table, figure, report, deck, audit from results/, then lint, types, tests
```

Without `make`: `.venv\Scripts\python.exe -m genomeasm app run --browser`. LKH-3 needs Ubuntu in WSL and the
source in `tools/` (research-use licence, not redistributed, see `tools/README.md`); everything else runs without it.

## Where things are

| What | Where |
|---|---|
| Beginner guide, presentation script, glossary | `docs/GUIDE.md`, `docs/DEMO_SCRIPT.md`, `docs/GLOSSARY.md` |
| Report (PDF), slides (pptx and PDF), oral-exam sheet | `reports/report/report.pdf`, `reports/slides/deck.pptx`, `deck.pdf`, `docs/ORAL_EXAM.md` |
| Stored pictures of every page (demo fallback) | `reports/screenshots/` |
| Decisions, plan, phase reports | `docs/DECISIONS.md`, `docs/PLAN.md`, `docs/phases/` |
| Code | `genomeasm/` (library), `app/` (Streamlit), `tests/` |
| Data, results, figures | `data/raw` (NCBI genomes with SHA-256), `results/`, `reports/figures/` |

## Commands

`genomeasm data verify`, `genomeasm experiments tune|run|scale|analyze` (hours, resumable), `genomeasm docs build|verify|screenshots`,
`genomeasm app run|benchmark`. `make experiments` reruns the Phase 5 experiments (hours); `make all` does not.

## Honesty notes

The report's references are marked `TODO-VERIFY`; the pptx was not opened in PowerPoint; "best known" is the best solution
found by our solvers, not a proven optimum; two algorithms (`tabu_cl`, `ucb1_rt`) were added after the first results and are
labelled exploratory. See `docs/phases/PHASE_8.md`.
