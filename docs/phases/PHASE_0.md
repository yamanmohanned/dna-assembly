# PHASE 0: Plan and environment

Status: **complete, waiting for your "continue"** to start Phase 1. Dates: 2026-10-22. Plan approved and deadline 2026-11-01 confirmed by you the same day.

## What was done
1. Located the master prompt. `docs/MASTER_PROMPT.md` did not exist; the only file was `CLAUDE_PROMPT_genome_assembly.md` in the project root. I copied it to `docs/MASTER_PROMPT.md` with the three placeholders filled (original untouched). See D-001.
2. Collected the placeholders: deadline 2026-11-01, language English, own laptop offline with a projector.
3. Audited the environment (below) and wrote `docs/PLAN.md` (schedule, scope tiers, risk register, acceptance criteria per phase). You approved it.
4. After approval, built the scaffolding:
   - Git repository on branch `main`; layout from the master prompt (`genomeasm/`, `app/`, `configs/`, `data/raw`, `data/processed`, `results/`, `reports/{figures,report,slides}`, `tests/`, `docs/`; empty folders hold a `.gitkeep`).
   - Project virtual environment `.venv` (Python 3.12.7) with the stack from the master prompt; exact versions in `requirements.lock`, direct packages in `requirements.txt` (D-013).
   - `pyproject.toml` (pytest, ruff, mypy settings), `Makefile` (`help venv test lint format typecheck check`), `.gitignore`, `tools/README.md`.
   - `tests/test_environment.py`: 27 smoke tests (Python version, imports, a numba JIT check, layout, documents present).
   - Installed GNU make 4.4.1 with winget at your request (D-012).
5. Retried the LKH-3 download; still unreachable (below).

## Evidence

### Environment audit (2026-10-22, before approval)

| Check | Result |
|---|---|
| OS | Windows 10 Enterprise 10.0.19045, 64-bit |
| CPU | Intel Core i7-4800MQ @ 2.70 GHz, 4 cores, 8 logical processors |
| RAM / disk | 15.6 GB RAM; E: has 40.1 GB free |
| Python | 3.12.7, pip upgraded to 26.2.2 inside the venv |
| Git | 2.53.0.windows.1; `user.name` and `user.email` set |
| gcc / g++ | MinGW.org GCC 6.3.0-1, `gcc -dumpmachine` = `mingw32` (32-bit) |
| clang, cl, cmake, pandoc, LaTeX, tectonic | not found |
| WSL | Ubuntu-24.04 (WSL2, 8 CPUs, about 7 GB RAM): has python3, curl, wget, git; **no gcc, no make, no pip3** |
| Chrome, Edge | installed |
| PowerPoint, LibreOffice | not found at standard paths |
| HTTP | eutils.ncbi.nlm.nih.gov (GET einfo) 200; www.ncbi.nlm.nih.gov 200; pypi.org 200; github.com 200 |

Preliminary accession check (`esummary`, db `nuccore`; **not** the Phase 1 verification, no sequence downloaded, no SHA-256):

| Accession.version | Organism | slen | moltype |
|---|---|---|---|
| NC_001422.1 | Escherichia phage phiX174 | 5386 | dna |
| NC_001416.1 | Escherichia phage Lambda | 48502 | dna |

### After approval

- `make check` (GNU make 4.4.1, run by full path) exit code 0:
  - `ruff check .`: All checks passed
  - `ruff format --check .`: 7 files already formatted
  - `mypy tests`: no issues found in 1 source file
  - `pytest`: 27 passed (8.9 s on the second run; 12.6 s on the first, which included the numba compile)
- Key versions resolved in the venv: numpy 2.5.4, numba 0.68.0, pandas 3.1.0, scipy 1.18.1, biopython 1.88, plotly 7.1.0, kaleido 1.5.0, matplotlib 3.11.2, streamlit 1.66.1, pytest 9.1.1, hypothesis 6.170.5, ruff 0.17.1, mypy 2.4.0, python-pptx 1.0.2, PyYAML 6.0.3.
- numba 0.68.0 compiles and runs correctly with numpy 2.5.4 (the JIT test passes). This resolves the numba/numpy question from the audit for this machine.
- Plotly static export through kaleido works on this machine with the installed Chrome: SVG (7,174 bytes, 3.7 s), PNG (44,233 bytes, 8.1 s), PDF (9,754 bytes, 8.0 s). Chart export (R-15) is therefore feasible.
- LKH-3 retries on 2026-10-22: `webhotel4.ruc.dk` port 80 and 443 not reachable; `akira.ruc.dk` port 80 accepts the TCP connection but then times out or resets (curl exit 28/56); the built-in browser could not open the page either. Windows and WSL both failed.

## Not verified (honest gaps)
- **LKH-3**: no source in hand, so buildability on MinGW GCC 6.3.0 (32-bit) or in WSL is still `TODO-VERIFY`.
- `make venv` (rebuild `.venv` from `requirements.lock`) has not been run on a clean checkout. Phase 8 does the fresh-environment install.
- Other `make` targets besides `check`'s components (`lint`, `typecheck`, `test`) were run only through `make check`; `help` also ran. `format` was not run.
- The lock file is for Windows and Python 3.12.7 only.
- Whether PowerPoint or LibreOffice exists in a non-standard location is unknown (question Q3).
- NCBI Entrez current usage rules (rate limit, contact email): `TODO-VERIFY` in Phase 1.

## Deviations from the master prompt
- Before approval I created only documents and an empty git repository, because you said no code before approval (D-005). The venv, Makefile and scaffolding were created after approval.
- The "CI-style `make test`" has one test file with environment smoke tests; real tests start in Phase 1.
- `make` is installed per user through winget; new terminals pick it up on PATH, this session called it by full path.

## Open risks (full register in `docs/PLAN.md` Section 6)
- **R-02, LKH-3 unreachable: open and now the top blocker.** Action needed from you, Q6: download the LKH-3 tarball from the official page in your own browser and put it in `tools/`. Separately, to build it I need either WSL `build-essential` (needs your sudo password) or a native MinGW build to be tried first. Decision point from the plan: if no source is in hand by the end of Phase 2 (about Oct 23), LKH-3 is reported as unavailable (D-011).
- R-01, tight schedule: unchanged. Phase 0 finished within its 0.3-day window.
- R-04 (PDF toolchain) and R-05 (no slide renderer) remain open; neither blocks Phase 1.

## Test output
```
27 passed in 8.87s    (make check: ruff, ruff format --check, mypy, pytest all exit 0)
```
