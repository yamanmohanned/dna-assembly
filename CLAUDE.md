# CLAUDE.md: project memory

Project: DNA fragment assembly as an asymmetric TSP, comparing Greedy, Tabu Search, LKH-3 and an ALNS whose operators are chosen by a multi-armed bandit.
Requirements: `docs/MASTER_PROMPT.md` (authoritative). Plan (approved): `docs/PLAN.md`. Decisions: `docs/DECISIONS.md`. Phase reports: `docs/phases/PHASE_N.md`. Data pipeline for the oral exam: `docs/DATA_LOADING.md`. Objective function and ATSP proof for the oral exam: `docs/OBJECTIVE_FUNCTION.md`. Definitions: `docs/glossary.yaml`.

## Status
- **Phases 0 to 5 complete; waiting for the user's "continue" to start Phase 6** (Streamlit application, offline demo; figures come from `results/`). Phase 5 report: `docs/phases/PHASE_5.md` (generated from the results; tests `tests/test_phase5_results.py`).
- Deadline 2026-11-01 (confirmed). Language: English. Presentation: own laptop, offline, projector. About 2 days ahead of schedule after Phase 5.
- **LKH-3 is available** (D-033): version 3.0.6 supplied by the user in `tools/LKH3-master/` (research-use licence, never committed), built in WSL Ubuntu 24.04 (the user installed build-essential and prefers to run any sudo command themselves; keep token use low). Provenance caveat: identity with the official release not verifiable. LKH-3 is only used through `genomeasm/solvers/lkh.py` (verifies every tour) and is never replaced silently (R9).
- Open questions: Q3 (program that opens the .pptx), Q5 (is this the presenting laptop).
- **Phase 5 verdict (honest):** the pre-registered claim (UCB1 and epsilon-greedy beat random operator selection) is NOT supported (UCB1 p_holm 0.90, A12 0.52; epsilon-greedy slightly worse). LKH-3 was never beaten on 300 read sets; ALNS beats greedy and tabu. Post hoc (D-045, D-046, exploratory): `tabu_cl` (candidate-list tabu) and `alns_ucb1_rt` (UCB1 with random tie-breaking; the original UCB1 had a defect: it plays arm 0 on a plateau). `alns_ucb1_rt` is slightly better than random (A12 0.57) and equal to roulette. Keep both versions visible in any slide; never present post-hoc results as confirmation.
- **WSL is down again** (`Wsl/Service/E_UNEXPECTED`, 2026-10-24): LKH-3 cannot run and `tests/test_solvers.py::test_lkh_finds_the_optimum...` fails until the user restarts WSL (reboot, or admin `Restart-Service LxssManager`). `wsl --shutdown` does not help. First action next session: check `wsl -d Ubuntu-24.04 -- echo ok`, ask the user if it fails.
- Lessons: (1) never let the laptop sleep during timed runs (keep-awake); (2) a policy surprise (an arm used 60% of the time) is a bug until proven otherwise; (3) tuning results do not always transfer from the tuning set to the evaluation set; (4) in PowerShell never pass a single `@(@(a,b))` pair to a helper (it flattens): edit files with Python.
- Key findings so far: D-031 (errors fragment the true-order chain), D-036 (error-free instances are trivial: greedy = LKH-3 = true order; solver comparison needs noisy instances with tolerant overlaps; the true order has LOWER F than the heuristics' solutions there, so biological metrics are essential).

## Working protocol (from the master prompt)
- One phase at a time. End of phase: full tests, `docs/phases/PHASE_N.md`, update this file, commit, stop, wait for "continue".
- Real NCBI data only. No hard-coded numbers; everything comes from `results/` or the data cards. Mark unverified claims `TODO-VERIFY`. No silent substitution of LKH-3.
- Explain concepts first in plain language, then precisely. Code, comments and identifiers in English; report, slides and app text in English.
- The user's account email is never sent to a service (D-010); the user chose to send no email to NCBI.
- If a validation check fails or a test surprises you, investigate before changing anything; record any change of method in `docs/DECISIONS.md` (D-022, D-031 are examples).
- Generate results at a clean commit: commit the code first, then run the experiment, then commit `results/` and figures (provenance records commit and dirty flag; `results/` and `reports/` are ignored by the dirty flag).

## What exists
- `genomeasm/io.py`: FASTA parser, `Genome`, NCBI download, offline cache. Alphabet rule D-016: only A/C/G/T.
- `genomeasm/simulate.py`: `ReadLengthModel`, `SimulationConfig`, `ReadSet` (sequences, starts, spans, error counts = ground truth), `simulate_reads`. Model D-023.
- `genomeasm/coverage.py`: depth, uncovered fraction, true contigs, Lander-Waterman, exact finite-genome model. The exact `exp(-c)` does NOT hold on short linear genomes (edge effect).
- `genomeasm/overlap.py`: `make_substring_free`, `overlap_matrix` (fast), `overlap_matrix_naive` (reference), `overlap_matrix_tolerant` (numba, exhaustive), `recommend_min_overlap`, `truth_agreement`.
- `genomeasm/objective.py`: `fitness` (F, maximise), `merged_length` (L), `cost_matrix` ((n+1)x(n+1), dummy node 0, read i = node i+1), `tour_from_order`, `order_from_tour`, `cycle_cost`, `write_tsplib_atsp`, `read_tsplib_atsp`.
- `genomeasm/assembly.py`: `layout`, `consensus`, `assemble`, `merged_string`. `genomeasm/instance.py`: `AssemblyInstance.from_readset(reads, min_overlap="auto", tolerance=0.0)` with `.ov`, `.lengths`, `.cost()`, `.fitness(order)`, `.merged_length(order)`, `.assemble(order)`, `.true_order()`, `.write_tsplib(path, name)`, `.summary()`.
- `genomeasm/checks.py` (shared `add_check`, `mean_and_se`), `sim_validation.py`, `overlap_validation.py`, `plots.py`, `provenance.py`, `datacard.py`, `cli.py`.
- Real data (committed): `data/raw/NC_001422.1.*` (phiX174) and `NC_001416.1.*` (lambda) plus `SHA256SUMS` (LF endings, `-text`); cards in `data/processed/`. Genomes are one linear string (D-017).
- Results (committed): `results/phase2/validation.json` (93 checks, commit 975eaa7), `results/phase3/overlap_validation.json` (95 checks, commit 460ea19) and `results/phase3/skeleton_true|skeleton_trivial/`. Figures `reports/figures/`: `fig01a`, `figS1`, `fig02_overlap_length_distribution`, `fig03_overlap_heatmap` (SVG, PDF, 300-dpi PNG, caption JSON).
- Key facts for Phase 4/5: true order = `inst.true_order()`; its F is the reference line (46,235 for phiX174 at c = 10, seed 1; trivial order F = 435). Recommended tau by `recommend_min_overlap` (9 to 14 for the real instances). With read errors the true-order chain breaks into many contigs (D-031): measure contigs, N50, genome fraction, not only F. Delta evaluation (O(1) swap/insert) is a Phase 4 task (D-032).

- `genomeasm/solvers/`: `common.py` (`Problem` = cycle form, weights W (n+1)x(n+1) with dummy node 0; `Budget(seconds, evaluations)`; `SolveResult` with `order`, `fitness`, `merged_length`, `evaluations`, `seconds`, `trace`, `info`), `moves.py` (O(1) `swap_delta`, `insert_delta`), `exact.py` (`brute_force` n<=9, `held_karp` n<=15), `greedy.py` (`greedy`, `nearest_neighbour`), `tabu.py` (`tabu_search`), `alns.py` (`alns(problem, budget, seed, policy=..., fixed_arm=...)`, arms = 3 destroy x 2 repair, `ARM_NAMES`), `lkh.py` (`lkh`, `lkh_binary`, `verify_tour`), `warmup.py` (`warmup()`: ALWAYS call before timing). `genomeasm/bandit/policies.py` (random, roundrobin, roulette, epsilon_greedy, ucb1 with sliding window, fixed). `genomeasm/solver_validation.py`.
- `genomeasm/bandit/policies.py`: random, roundrobin, roulette, epsilon_greedy, ucb1 (lowest-index ties, kept for reproducibility), ucb1_rt (random ties, D-046), fixed. `genomeasm/experiments/`: `design.py` (instances, `run_algorithm`, bio metrics), `runner.py` (resumable pool), `plan.py` (tuning, tuned params, tasks, `missing_tasks`), `analysis.py` (statistics, verdicts, convergence, scaling), `figures.py` (figs 4 to 11), `bio.py`. Config `configs/phase5_experiments.yaml` (D-039 to D-046).
- Results (committed): `results/phase4/correctness.json` (12 checks, commit 07f63ea; includes a one-instance preview comparison, not a result).

## Environment (Windows 10, i7-4800MQ 4C/8T, 15.6 GB RAM)
- Python 3.12.7 in `.venv`; exact packages in `requirements.lock`; package installed editable (`pip install -e . --no-deps`). numba kernels use `cache=True` and `parallel=True`.
- GNU make 4.4.1 via winget. Until a shell is restarted use `%LOCALAPPDATA%\Microsoft\WinGet\Links\make.exe`.
- MinGW GCC 6.3.0 (32-bit). WSL Ubuntu 24.04 has no compiler (`build-essential` needs the user's sudo).
- Plotly static export works through kaleido with the installed Chrome; matplotlib figures use the Agg backend.
- No LaTeX/pandoc (PDF toolchain open, D-007). No PowerPoint/LibreOffice found.
- NCBI usage guidelines page unreadable (CAPTCHA); rate limit rule only verified via Biopython source (D-019).

## Commands (all through `.venv\Scripts\python.exe`)
- `make check`: ruff check, ruff format --check, mypy (`genomeasm tests`), pytest. Green at end of Phase 5 except the LKH-3 unit test while WSL is down (see Status).
- `make test`, `make lint`, `make format`, `make typecheck`, `make venv`, `make help`.
- `genomeasm solve ACC --coverage C --solver greedy|nn|tabu|alns|lkh [--policy ucb1] [--seconds S] [--sub-rate P --tolerance T --min-overlap N]`; `genomeasm solvers validate` (about 10 min, writes `results/phase4/correctness.json`). Run LKH-3 by hand: `wsl -d Ubuntu-24.04 -- bash -c 'cd "<wsl dir>" && "<wsl>/tools/LKH3-master/LKH" file.par'`.
- `genomeasm experiments tune|run|scale|analyze` (resumable; `rerun-lkh --target main|scaling` repeats LKH-3 only where missing; `add-runs --target main|scaling --algorithm A` adds one algorithm where it has no result). Work files in `data/processed/phase5_work/` (ignored); raw runs, tables and `stages.json` in `results/phase5/`. Never launch a timed run with a dirty code tree; request keep-awake first.
- `genomeasm data fetch|verify|cards`; `genomeasm simulate run|validate|figures`; `genomeasm overlap validate|figures` (about 170 s); `genomeasm assemble ACC --coverage C --order true|trivial|random [--min-overlap auto|N] [--tolerance T] [--sub-rate P] [--out-dir D]`.

## Conventions
- Small conventional commits. Generated results are never committed without config and seed. Commit message ends with the Co-Authored-By line given in the session.
- Layout: `genomeasm/`, `app/`, `configs/`, `data/raw`, `data/processed`, `results/`, `reports/{figures,report,slides}`, `tests/`, `docs/`. `tools/` holds local, uncommitted third-party programs (LKH-3). `data/processed/reads/` and `data/processed/assemblies/` are git-ignored.
- Ruff line length 100, target py311; mypy `disallow_untyped_defs`. With numpy typed arrays use `np.min(a)` / `np.max(a)` / `np.sum(a)` instead of `a.min()` (mypy overload false alarm). Comment banners and long strings are not wrapped by `ruff format`: keep them under 100 characters by hand. Tests that need the network use canned Entrez answers; the `no_network` fixture blocks sockets.
- Statistical checks: tolerances live in the config, fixed before the run; seeds come from `SeedSequence`; a smoke test with few runs cannot replace the full validation. Test expectations must be justified by theory, never copied from a first run.
- Windows file writing: text files written from Python use `newline="\n"`; SVG output is normalised to LF; never let git rewrite `data/raw`.
