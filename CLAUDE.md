# CLAUDE.md: project memory

Project: DNA fragment assembly as an asymmetric TSP, comparing Greedy, Tabu Search, LKH-3 and an ALNS whose operators are chosen by a multi-armed bandit.
Requirements: `docs/MASTER_PROMPT.md` (authoritative). Plan (approved): `docs/PLAN.md`. Decisions: `docs/DECISIONS.md`. Phase reports: `docs/phases/PHASE_N.md`. Data pipeline for the oral exam: `docs/DATA_LOADING.md`. Definitions: `docs/glossary.yaml`.

## Status
- **Phases 0, 1 and 2 complete; waiting for the user's "continue" to start Phase 3** (overlap graph and ATSP formulation).
- Deadline 2026-11-01 (confirmed). Language: English. Presentation: own laptop, offline, projector. Ahead of schedule by about half a day after Phase 2.
- Top blocker: LKH-3 source unreachable (official hosts time out, `tools/` empty). The user was asked to download it into `tools/` (Q6). Proposal D-025: final decision at the start of Phase 4 (about Oct 24). No source then means LKH-3 is reported as unavailable, never substituted (R9, D-011); the user decides before any extra baseline is added.
- Open questions: Q3 (program that opens the .pptx), Q5 (is this the presenting laptop), Q6 (LKH-3 tarball), `build-essential` in WSL (needs the user's sudo).

## Working protocol (from the master prompt)
- One phase at a time. End of phase: full tests, `docs/phases/PHASE_N.md`, update this file, commit, stop, wait for "continue".
- Real NCBI data only. No hard-coded numbers; everything comes from `results/` or the data cards. Mark unverified claims `TODO-VERIFY`. No silent substitution of LKH-3.
- Explain concepts first in plain language, then precisely. Code, comments and identifiers in English; report, slides and app text in English.
- The user's account email is never sent to a service (D-010); the user chose to send no email to NCBI.
- If a validation check fails, investigate before changing anything; record any change of method in `docs/DECISIONS.md` (see D-022).
- Generate results at a clean commit: commit the code first, then run the experiment, then commit `results/` and figures (provenance records commit and dirty flag; `results/` and `reports/` are ignored by the dirty flag).

## What exists
- `genomeasm/io.py`: FASTA parser, `Genome` (accession, sequence, base_counts, length, gc_fraction, sequence_sha256), NCBI download, offline cache. Alphabet rule D-016: only A/C/G/T in a `Genome`.
- `genomeasm/datacard.py`, `genomeasm/cli.py` (`genomeasm data fetch|verify|cards`, `genomeasm simulate run|validate|figures`; also `python -m genomeasm`).
- `genomeasm/simulate.py`: `ReadLengthModel`, `SimulationConfig`, `ReadSet` (sequences, starts, spans, error counts = ground truth), `simulate_reads`, `sample_geometry`, `streams`. Model D-023.
- `genomeasm/coverage.py`: depth, uncovered fraction, true contigs, Lander-Waterman, exact finite-genome model, exact expected contigs (D-021), bin grouping.
- `genomeasm/sim_validation.py`, `genomeasm/plots.py` (Okabe-Ito palette, `save_figure`), `genomeasm/provenance.py` (`collect`, `config_hash`, `git_info`).
- Real data (committed): `data/raw/NC_001422.1.*` (phiX174, circular ss) and `NC_001416.1.*` (lambda, linear ds) plus `SHA256SUMS` (LF endings, `-text` in `.gitattributes`); cards in `data/processed/`. Genomes are used as one linear string (D-017).
- Results (committed): `results/phase2/validation.json` (93 of 93 checks passed, run at commit 975eaa7); figures `reports/figures/fig01a_coverage_histogram.*` and `figS1_validation_curves.*`.
- Key fact for later phases: the true read order is `argsort(ReadSet.starts)`; the exact `exp(-c)` of Lander-Waterman does NOT hold on short linear genomes (edge effect); use the exact model in `coverage.py` for predictions.
- Next module planned: `genomeasm/overlap.py` (Phase 3) taking a `ReadSet`.

## Environment (Windows 10, i7-4800MQ 4C/8T, 15.6 GB RAM)
- Python 3.12.7 in `.venv`; exact packages in `requirements.lock` (direct list in `requirements.txt`); package installed editable (`pip install -e . --no-deps`).
- GNU make 4.4.1 via winget. Until a shell is restarted use `%LOCALAPPDATA%\Microsoft\WinGet\Links\make.exe`.
- MinGW GCC 6.3.0 (32-bit). WSL Ubuntu 24.04 has no compiler (`build-essential` needs the user's sudo).
- Plotly static export (SVG/PNG/PDF) works through kaleido with the installed Chrome; matplotlib figures use the Agg backend.
- No LaTeX/pandoc (PDF toolchain open, D-007). No PowerPoint/LibreOffice found.
- NCBI usage guidelines page unreadable (CAPTCHA); rate limit rule only verified via Biopython source (D-019).

## Commands (all through `.venv\Scripts\python.exe`)
- `make check`: ruff check, ruff format --check, mypy (`genomeasm tests`), pytest. Green at end of Phase 2: 202 tests.
- `make test`, `make lint`, `make format`, `make typecheck`, `make venv`, `make help`.
- `genomeasm data fetch` (needs internet), `genomeasm data verify`, `genomeasm data cards` (offline).
- `genomeasm simulate run ACC --coverage C ...` (reads + ground truth to `data/processed/reads/`, not committed), `genomeasm simulate validate` (about 30 s, rewrites `results/phase2/` and the figures), `genomeasm simulate figures` (redraw only).

## Conventions
- Small conventional commits. Generated results are never committed without config and seed. Commit message ends with the Co-Authored-By line given in the session.
- Layout: `genomeasm/`, `app/`, `configs/`, `data/raw`, `data/processed`, `results/`, `reports/{figures,report,slides}`, `tests/`, `docs/`. `tools/` holds local, uncommitted third-party programs (LKH-3).
- Ruff line length 100, target py311; mypy `disallow_untyped_defs`. With numpy typed arrays use `np.min(a)` / `np.max(a)` instead of `a.min()` (mypy overload false alarm). Tests that need the network use canned Entrez answers; the `no_network` fixture blocks sockets.
- Statistical checks: tolerances live in the config, fixed before the run; seeds come from `SeedSequence`; a smoke test with few runs cannot replace the full validation.
