# CLAUDE.md: project memory

Project: DNA fragment assembly as an asymmetric TSP, comparing Greedy, Tabu Search, LKH-3 and an ALNS whose operators are chosen by a multi-armed bandit.
Requirements: `docs/MASTER_PROMPT.md` (authoritative). Plan (approved): `docs/PLAN.md`. Decisions: `docs/DECISIONS.md`. Phase reports: `docs/phases/PHASE_N.md`. Data pipeline for the oral exam: `docs/DATA_LOADING.md`. Definitions: `docs/glossary.yaml`.

## Status
- **Phases 0 and 1 complete; waiting for the user's "continue" to start Phase 2** (sequencing simulator).
- Deadline 2026-11-01 (confirmed). Language: English. Presentation: own laptop, offline, projector.
- Top blocker: LKH-3 source unreachable (official hosts time out). The user was asked to download it into `tools/` (Q6). No source in hand means LKH-3 is reported as unavailable, never substituted (R9, D-011). Decision point: end of Phase 2 (about Oct 23).
- Open questions: Q3 (program that opens the .pptx), Q5 (is this the presenting laptop), Q6 (LKH-3 tarball), `build-essential` in WSL (needs the user's sudo).

## Working protocol (from the master prompt)
- One phase at a time. End of phase: full tests, `docs/phases/PHASE_N.md`, update this file, commit, stop, wait for "continue".
- Real NCBI data only. No hard-coded numbers; everything comes from `results/` or the data cards. Mark unverified claims `TODO-VERIFY`. No silent substitution of LKH-3.
- Explain concepts first in plain language, then precisely. Code, comments and identifiers in English; report, slides and app text in English.
- The user's account email is never sent to a service (D-010); the user chose to send no email to NCBI.

## What exists
- `genomeasm/io.py`: FASTA parser, `Genome` (accession, sequence, base_counts, length, gc_fraction, sequence_sha256), NCBI download, offline cache (`load_genome`, `download_genome`). Alphabet rule D-016: only A/C/G/T allowed in a `Genome`.
- `genomeasm/datacard.py`, `genomeasm/cli.py` (`genomeasm data fetch|verify|cards`, also `python -m genomeasm`).
- Real data (committed): `data/raw/NC_001422.1.*` (phiX174, circular ss) and `NC_001416.1.*` (lambda, linear ds) plus `SHA256SUMS`; cards in `data/processed/`. Genomes are used as one linear string (D-017).
- Next module planned: `genomeasm/simulate.py` (Phase 2) taking a `Genome` and nothing else.

## Environment (Windows 10, i7-4800MQ 4C/8T, 15.6 GB RAM)
- Python 3.12.7 in `.venv`; exact packages in `requirements.lock` (direct list in `requirements.txt`); package installed editable (`pip install -e . --no-deps`).
- GNU make 4.4.1 via winget. Until a shell is restarted use `%LOCALAPPDATA%\Microsoft\WinGet\Links\make.exe`.
- MinGW GCC 6.3.0 (32-bit). WSL Ubuntu 24.04 has no compiler (`build-essential` needs the user's sudo).
- Plotly static export (SVG/PNG/PDF) works through kaleido with the installed Chrome.
- No LaTeX/pandoc (PDF toolchain open, D-007). No PowerPoint/LibreOffice found.
- NCBI usage guidelines page unreadable (CAPTCHA); rate limit rule only verified via Biopython source (D-019).

## Commands (all through `.venv\Scripts\python.exe`)
- `make check`: ruff check, ruff format --check, mypy (`genomeasm tests`), pytest. Green at end of Phase 1: 124 tests, 97 % coverage of `genomeasm`.
- `make test`, `make lint`, `make format`, `make typecheck`, `make venv`, `make help`.
- `genomeasm data fetch` (needs internet), `genomeasm data verify`, `genomeasm data cards` (offline).

## Conventions
- Small conventional commits. Generated results are never committed without config and seed. Commit message ends with the Co-Authored-By line given in the session.
- Layout: `genomeasm/`, `app/`, `configs/`, `data/raw`, `data/processed`, `results/`, `reports/{figures,report,slides}`, `tests/`, `docs/`. `tools/` holds local, uncommitted third-party programs (LKH-3).
- Ruff line length 100, target py311; mypy `disallow_untyped_defs`. Tests that need the network are replaced by canned Entrez answers; `no_network` fixture blocks sockets.
