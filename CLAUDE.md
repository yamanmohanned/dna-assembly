# CLAUDE.md: project memory

Project: DNA fragment assembly as an asymmetric TSP, comparing Greedy, Tabu Search, LKH-3 and an ALNS whose operators are chosen by a multi-armed bandit.
Requirements: `docs/MASTER_PROMPT.md` (authoritative). Plan (approved): `docs/PLAN.md`. Decisions: `docs/DECISIONS.md`. Phase reports: `docs/phases/PHASE_N.md`.

## Status
- **Phase 0 complete; waiting for the user's "continue" to start Phase 1** (data acquisition and validation).
- Deadline 2026-11-01 (confirmed). Language: English. Presentation: own laptop, offline, projector.
- Top blocker: LKH-3 source unreachable (official hosts time out). The user was asked to download it into `tools/` (Q6). No source in hand means LKH-3 is reported as unavailable, never substituted (R9, D-011). Decision point: end of Phase 2.

## Working protocol (from the master prompt)
- One phase at a time. End of phase: full tests, `docs/phases/PHASE_N.md`, update this file, commit, stop, wait for "continue".
- Real NCBI data only. No hard-coded numbers; everything comes from `results/`. Mark unverified claims `TODO-VERIFY`. No silent substitution of LKH-3.
- Explain concepts first in plain language, then precisely. Code, comments and identifiers in English; report, slides and app text in English.
- Do not send the user's account email to any service unless told to (D-010).

## Environment (Windows 10, i7-4800MQ 4C/8T, 15.6 GB RAM)
- Python 3.12.7 in `.venv`; exact packages in `requirements.lock` (direct list in `requirements.txt`).
- GNU make 4.4.1 via winget. Until a shell is restarted use `%LOCALAPPDATA%\Microsoft\WinGet\Links\make.exe`.
- MinGW GCC 6.3.0 (32-bit). WSL Ubuntu 24.04 has no compiler (`build-essential` needs the user's sudo).
- Plotly static export (SVG/PNG/PDF) works through kaleido with the installed Chrome.
- No LaTeX/pandoc (PDF toolchain open, D-007). No PowerPoint/LibreOffice found.

## Commands (all through `.venv\Scripts\python.exe`)
- `make check`: ruff check, ruff format --check, mypy, pytest (all green at end of Phase 0, 27 tests).
- `make test`, `make lint`, `make format`, `make typecheck` (mypy covers `tests` only until Phase 1), `make venv`, `make help`.

## Conventions
- Small conventional commits. Generated results are never committed without config and seed.
- Layout: `genomeasm/`, `app/`, `configs/`, `data/raw`, `data/processed`, `results/`, `reports/{figures,report,slides}`, `tests/`, `docs/`. `tools/` holds local, uncommitted third-party programs (LKH-3).
- Ruff line length 100, target py311; mypy `disallow_untyped_defs`.
