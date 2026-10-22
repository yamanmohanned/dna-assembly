# PLAN: DNA Fragment Assembly as TSP with a Learning-Based ALNS

Status: **APPROVED by you on 2026-10-22**, with the deadline 2026-11-01 confirmed (Q1 answered). The text below is the approved snapshot; later changes are logged in `docs/DECISIONS.md` and `docs/phases/PHASE_0.md`.
Written: 2026-10-22 (Thursday). Source of requirements: `docs/MASTER_PROMPT.md`.
Evidence for the environment section: `docs/phases/PHASE_0.md`. Choices I made on your behalf: `docs/DECISIONS.md`.

| Input | Value |
|---|---|
| Deadline | **2026-11-01 (Sunday)**, read from your answer "2026-11-1". **Please confirm** (see D-002) |
| Language of report, slides and in-app text | English (code stays English) |
| Presentation setup | Live demo on your own laptop, no internet, projector. The app must run fully offline |
| Today | 2026-10-22 (Thursday) = day 0 |
| Days available | **10 working days** (2026-10-22 to 2026-10-31 inclusive). The deadline day itself is for presenting, not building |

---

## 1. Feasibility verdict (my judgment, not a measurement)

Ten days is very tight for this scope. The Must tier (core pipeline, three baselines plus ALNS-Bandit, ablation, app, report, 15+ slides) is reachable **only if all of the following hold**:

1. LKH-3 source can be fetched and built within about two days (it is currently unreachable, see Section 4 and risk R-02).
2. You answer each "continue" quickly. Every hour of waiting comes out of the schedule, and the only slack is Phase 8 (0.8 days).
3. The Should and Could tiers are dropped unless there is spare time after the Phase 5 gate.

What each deadline would realistically allow:

| If the deadline is | Calendar days | Scope that fits (my estimate) |
|---|---|---|
| 2026-11-01 | 10 | Must only, with 30 seeds on the main comparison and ablation. High schedule risk |
| 2026-11-10 | 19 | Must plus most of Should (sensitivity, scaling, error-tolerant overlaps) |
| 2026-11-19 | 28 | Must plus Should, with room for one Could item. Moderate risk |

If 2026-11-01 is wrong, tell me the right date and I will regenerate Section 2 before you approve.

---

## 2. Backward schedule from 2026-11-01

Method: the phase shares from the master prompt (3, 5, 10, 10, 22, 15, 12, 15, 8 %) applied to 10 days. Fractions are fractions of a working day, not clock hours. I assume every calendar day, weekends included, is a working day (D-003).

| # | Phase | Share | Length | Window (day 0 = Oct 22) | Calendar window |
|---|---|---|---|---|---|
| 0 | Plan and environment | 3% | 0.3 d | 0.0 to 0.3 | Oct 22 |
| 1 | Data acquisition and validation | 5% | 0.5 d | 0.3 to 0.8 | Oct 22 |
| 2 | Sequencing simulator | 10% | 1.0 d | 0.8 to 1.8 | Oct 22 to Oct 23 |
| 3 | Overlap graph and ATSP formulation | 10% | 1.0 d | 1.8 to 2.8 | Oct 23 to Oct 24 |
| 4 | Solvers | 22% | 2.2 d | 2.8 to 5.0 | Oct 24 to Oct 26 |
| 5 | Experiments and statistics | 15% | 1.5 d | 5.0 to 6.5 | Oct 27 to Oct 28 |
| 6 | Application | 12% | 1.2 d | 6.5 to 7.7 | Oct 28 to Oct 29 |
| 7 | Report, slides, oral-exam sheet | 15% | 1.5 d | 7.7 to 9.2 | Oct 29 to Oct 31 |
| 8 | QA, rehearsal, freeze | 8% | 0.8 d | 9.2 to 10.0 | Oct 31 |
| | Presentation | | | | **Nov 1** |

Checkpoints and rules tied to the schedule:

- **Walking skeleton** (end of Phase 3, about Oct 24): end-to-end CLI run on real phiX174 reads with a trivial ordering.
- **Feature freeze** at the end of Phase 6, about Oct 29. After that, only fixes, documents and rehearsal.
- **Final 20% reserve** (Oct 30 to Oct 31): shared by the tail of Phase 7 and all of Phase 8. The shares in the master prompt add to 100% and Phase 7 plus 8 are 23%, so the reserve and Phases 7 and 8 overlap (D-004). Report and slide *generators* are written early in Phase 7 so the last day is polish only.
- **Parallel work**: the long Phase 5 experiment runs go in the background while I build the Phase 6 app. This is the main way the 10 days become workable.
- **LKH-3 go/no-go**: the source must be in hand by the end of Phase 2 (Oct 23). If not, see the contingency under R-02.

---

## 3. Scope tiers and what is cut first

| Tier | Content | Rule |
|---|---|---|
| **Must** | Core pipeline; Greedy, Tabu, LKH-3, ALNS plus bandit policies; ablation; app; PDF report; 15+ slides; ORAL_EXAM.md | Never cut. Seeds can fall below 30 only with a written justification (R7) |
| **Should** | ALNS parameter sensitivity; runtime/gap scaling with n; error-tolerant overlaps | Cut in this order if Phase 4 ends late: scaling, sensitivity, error-tolerant |
| **Could** | Second genome family; two-strand reads | Not started unless Phase 8 finishes early. Not expected |

Two genomes are in scope by default: phiX174 (NC_001422.1) and lambda (NC_001416.1). With only two genomes, "held-out" (rule R6) will mean: tune on non-overlapping windows of one genome, evaluate on windows of the other genome and on unseen seeds. The exact split is fixed in Phase 5 and logged (D-009). A third genome family stays a Could item.

---

## 4. Environment audit (details and raw output in `docs/phases/PHASE_0.md`)

| Item | Found | Verdict / action |
|---|---|---|
| OS | Windows 10 Enterprise 10.0.19045, 64-bit | OK |
| CPU / RAM / disk | Core i7-4800MQ (2013), 4 cores / 8 threads, 15.6 GB RAM, about 40 GB free on E: | OK, but slow. Affects the 30-second live-demo limit and timing noise |
| Python | 3.12.7 with pip 25.3 and `venv` | OK (3.11+ needed) |
| Git | 2.53.0, user name and email configured | OK |
| **make** | **Missing** (not on PATH, not in `C:\MinGW\bin`) | **Needs your OK to install** (see Section 8). The Makefile and `make all` depend on it |
| C compiler | MinGW.org GCC 6.3.0, target `mingw32` (32-bit, old) | LKH-3 build **unverified** (TODO-VERIFY). Native build may need patching |
| WSL | Ubuntu 24.04 (WSL2), 8 CPUs; has python3, curl, wget, git; **no gcc, no make, no pip3** | Fallback build host for LKH-3 once `build-essential` is installed (needs sudo, your OK) |
| NCBI Entrez | Reachable (HTTP 200). `esummary` for both accessions worked | OK. Formal verification and SHA-256 are Phase 1 |
| PyPI, GitHub | Reachable (HTTP 200) | OK |
| **LKH-3 official site** | `webhotel4.ruc.dk` resolves to 130.225.220.230 but the TCP connection to port 80 **times out from both Windows and WSL** (3 attempts) | **Open blocker.** May be transient. See R-02 |
| Python packages (global) | numpy 2.3.5, pandas 2.3.3, scipy 1.16.3, matplotlib 3.10.8, streamlit 1.51.0, pyyaml, reportlab 4.4.9, jinja2, tqdm | Will not be used directly; a project venv will be created |
| Not installed | numba, biopython, plotly, pytest, hypothesis, ruff, mypy, python-pptx, kaleido, ortools | Install into the venv from PyPI after approval |
| PDF toolchain | No LaTeX, no pandoc. Chrome and Edge are installed | Report build tool is an open decision (D-007, risk R-04) |
| Slide software | PowerPoint and LibreOffice **not found** at standard paths | Cannot render the `.pptx` locally for visual QA. Question Q3 |

Genome accessions (preliminary check only, `esummary` on 2026-10-22; formal verification with SHA-256 in Phase 1):

| Accession | Organism (NCBI) | Length reported |
|---|---|---|
| NC_001422.1 | Escherichia phage phiX174 | 5,386 bp |
| NC_001416.1 | Escherichia phage Lambda (title: Enterobacteria phage lambda) | 48,502 bp |

---

## 5. Acceptance criteria per phase

Every phase ends the same way: full test suite green, `docs/phases/PHASE_N.md` written (done, evidence, test output, deviations, open risks), `CLAUDE.md` updated, small conventional commit, then I **stop and wait for your "continue"**.

**Phase 0: Plan and environment** (exit gate: your approval)
- This plan approved or amended by you.
- After approval: repo layout from the master prompt, project venv with the packages above, `Makefile` with `make test` and `make lint`, LKH-3 source obtained (or the blocker escalated), `CLAUDE.md` current.

**Phase 1: Data acquisition and validation**
- `genomeasm.io` downloads by accession via Entrez, parses FASTA, validates the alphabet, caches offline.
- Documented rule for `N` and IUPAC ambiguity codes, with tests for valid, empty, malformed, lowercase and ambiguity-code input.
- Data card for each genome (accession, version, organism, length, GC%, SHA-256) generated from the real files.
- Loading from cache works with the network off (tested).
- Glossary file started as the single source of definitions (reused by report, slides, app).

**Phase 2: Sequencing simulator**
- Parameters: read length (fixed or distribution), coverage, substitution rate, optional indel rate, seed. Output: reads plus the true origin of each.
- Statistical tests with tolerances fixed in advance: realised coverage vs `Nℓ/G`; uncovered fraction vs `e^{-c}`; empirical error rate vs configured rate (binomial interval); Monte-Carlo agreement.
- Figure 1 (coverage histogram vs Poisson) saved as SVG/PDF and 300-dpi PNG.
- Same seed gives byte-identical reads.

**Phase 3: Overlap graph and ATSP formulation**
- Mandatory tests (a) to (d) from the master prompt pass: toy example `ATGCG`, `GCGTA`, `GTACC` gives `[[–,3,1],[1,–,3],[0,0,–]]`, order A→B→C, merged `ATGCGTACC`; `L(π) = Σℓ − F(π)` equals real merged length; dummy-node cycle cost equals `L(π)`; fast overlap code equals a naive reference.
- Containment/duplicate removal with counts recorded; cost matrix; TSPLIB writer; layout → majority-vote consensus.
- **Walking skeleton**: CLI runs end to end on real phiX174 reads with a trivial ordering.

**Phase 4: Solvers**
- Greedy (global edge merging) and nearest-neighbour; Tabu (swap/insert, O(1) delta, tenure, aspiration, restarts); LKH-3 wrapper that re-verifies the returned tour on the cost matrix; ALNS with 3+ destroy and 2+ repair operators and calibrated simulated-annealing temperature; policies Random, RoundRobin, RouletteWheel, EpsilonGreedy, UCB1; each single operator alone; brute force (n ≤ 9) and Held–Karp (n ≤ 15).
- Mandatory tests: every heuristic reaches the exact optimum (or documents its gap) on random small instances; valid permutations; reproducible per seed; every delta equals a full recomputation.
- `PHASE_4.md` is the correctness report with evidence.

**Phase 5: Experiments and statistics**
- Tuning instances kept separate from evaluation instances; search space and results logged.
- Main comparison and ablation with equal budget (wall-clock and evaluation count), 30 seeds per configuration (fewer only with a written justification), paired by seed.
- Metrics from the master prompt (F, gap to LKH-3, runtime, convergence, biological quality, F of the true ordering).
- R7 statistics: median, IQR, bootstrap 95% CI, Wilcoxon signed-rank with Holm correction, Vargha–Delaney A12 or Cliff's delta. All as auto-generated tables naming the exact test.
- Every result file stores config hash, git commit, seed, library versions, machine info, timestamps. Conclusions drafted only from the numbers.

**Phase 6: Application**
- All eight Streamlit pages work fully offline; SVG and PNG export on every chart.
- Any interactive run finishes in under 30 s **on your laptop** (measured, not assumed); larger runs are shown from precomputed results and labelled as such.
- Feature freeze at the end of this phase.

**Phase 7: Report, slides, oral-exam sheet**
- PDF report (IMRaD) with proofs, injected numbers and figures, verified references only (anything unverified stays `TODO-VERIFY`).
- Slide deck of 15+ slides, each title stating its takeaway, speaker notes in plain language, generated from a script.
- `docs/ORAL_EXAM.md`: objective function, solution representation, full data-loading pipeline with file names, every hyper-parameter with its rationale, 20 examiner questions with model answers.
- One command builds all three.

**Phase 8: QA, rehearsal, freeze**
- `make all` reproduces every figure, table, report and deck from a clean checkout.
- Lint, type-check, tests and a fresh-environment install pass.
- A script checks every number in the slides against `results/`.
- Timed demo script with a fallback to precomputed results and screenshots. Release tagged.

---

## 6. Risk register

Likelihood (L) and impact (I): H high, M medium, L low. These ratings are my judgment.

| ID | Risk | L | I | Mitigation and trigger |
|---|---|---|---|---|
| R-01 | 10 days is too short for the full Must tier | H | H | Parallelise Phase 5 runs with app work; generators for report/slides built early; Should/Could cut first; confirm the deadline now (Section 1) |
| R-02 | LKH-3 source cannot be downloaded (official host unreachable) or cannot be built | H | H | Retry on a schedule; test other networks; ask you to download the tarball by hand and place it in `tools/`. If still missing by end of Phase 2: LKH-3 is reported as **unavailable** (R9). Clearly labelled extra baselines (e.g. OR-Tools) may be added, never called LKH-3; "gap" would then be measured against the best tour found by any solver, and said so. You decide before I switch |
| R-03 | Native Windows build of LKH-3 fails (old 32-bit MinGW, POSIX headers) | M | M | TODO-VERIFY. Fallback: build inside WSL after `build-essential` is installed. Any source patch is documented |
| R-04 | No LaTeX or pandoc for the PDF report; equations and proofs need good typesetting | M | M | Half-day spike early in Phase 3 to prove one candidate toolchain (D-007): HTML + KaTeX printed by headless Chrome, or Typst, or reportlab |
| R-05 | No PowerPoint/LibreOffice found, so the `.pptx` cannot be viewed locally | M | M | Structural checks with python-pptx plus a rendered image preview path; you open the deck once on the presenting laptop (Q3) |
| R-06 | Old laptop CPU: slow experiments, timing noise, 30 s demo limit missed | M | M | At most 6 parallel workers; interleave algorithms so load is equal; also report equal evaluations; measure demo runs early in Phase 3 |
| R-07 | Numba/numpy version conflict or long JIT compile times on Windows | M | L | Profile before using numba; pin versions in the venv; cache compiled functions |
| R-08 | Overlap matrix is n² (lambda at 20x with 100 bp reads is about 9,700 reads) | M | M | Use k-mer indexing; compare solvers on windowed instances of moderate n; document windows; int16/int32 storage |
| R-09 | Optimal tour differs from true genome order because of repeats | H | L | Expected. Report both F and biological quality, plus F of the true ordering (declared simplification) |
| R-10 | The bandit does not beat random selection | M | L | Acceptable. Report with effect size, CI and test whatever the outcome (R10: investigate surprises first) |
| R-11 | A result looks too good or contradicts theory | L | H | Debug before reporting (R10); exact-solver validation first (R5) |
| R-12 | Waiting time between phases (the "continue" gates) eats the schedule | M | H | Phase summaries kept short; you can reply "continue" as soon as you have skimmed them |
| R-13 | NCBI rate limits or Entrez contact-email requirement | L | L | TODO-VERIFY current rules in Phase 1; cache after one download; ask which email to use (Q4) |
| R-14 | Windows multiprocessing (spawn) plus numba caching problems | M | L | Test the experiment runner on a tiny config first; keep worker functions importable |
| R-15 | Plotly static export (kaleido) needs a working Chrome | L | M | Chrome exists on this machine; matplotlib is the fallback for vector figures |

---

## 7. Experiment compute budget (planning estimate only, not a result)

Assumption set, to be fixed properly in Phase 5: 12 evaluation instances, 10 time-limited configurations (Tabu, LKH-3, ALNS with 5 policies, 3 single-operator variants), 30 seeds, 10 s per run.

12 × 10 × 30 × 10 s = 36,000 CPU-seconds ≈ 10 CPU-hours. With 6 workers that is roughly 2 to 3 hours of wall time, plus about 25% for tuning. This fits the 1.5-day Phase 5 window if the runner is ready at the start and runs in the background. If LKH-3 or the runner slips, the seed count or instance count is what gets reduced, with the reason written down.

---

## 8. What I need from you

| # | Question | Why it matters | Blocking? |
|---|---|---|---|
| Q1 | Is the deadline **2026-11-01**? If you meant another date (for example 2026-11-10 or 2026-11-19), say so | Changes Section 1 and 2 | Yes, before approval |
| Q2 | May I install **GNU make** on Windows (via winget or chocolatey, chocolatey may need admin rights) and **`build-essential` in WSL Ubuntu** (needs your sudo password, which I will not type for you)? | `make all` and the LKH-3 fallback build | Needed right after approval |
| Q3 | Which program will open the final `.pptx` on the presenting laptop (PowerPoint, LibreOffice, other)? | Slide layout QA | No, by Phase 7 |
| Q4 | Which contact email may I send to NCBI Entrez in Phase 1? (I will not use your account email unless you say so) | Entrez etiquette | No, by Phase 1 |
| Q5 | Is the presenting laptop this machine (i7-4800MQ, 15.6 GB, Windows 10)? | The 30 s demo limit is measured on it | No, by Phase 3 |
| Q6 | Can you try to download the LKH-3 tarball in a browser from the official page if my retries keep failing? | R-02 | Only if retries fail |

---

## 9. Cross-cutting standards (apply to every phase)

- **No hard-coded numbers (R3)**: every number in the app, report and slides is generated from files in `results/`; every result file carries config hash, git commit, seed, library versions, machine info and timestamps.
- **Plain language first**: each concept is explained for a 13-year-old, then defined precisely; definitions live in one glossary file shared by app, report and slides.
- **Decisions log**: every ambiguity I resolve goes into `docs/DECISIONS.md`.
- **Unverified means marked**: formulas, citations and library behaviours I cannot verify stay `TODO-VERIFY` and are not presented as fact.
- **Commits**: small and conventional; generated results are never committed without their config and seed.

---

## 10. What happens when you approve

1. I take your answers to Q1 and Q2 (only these two block; Q3 to Q6 can wait for later phases).
2. I finish Phase 0: repo layout, venv and installs, `Makefile` with `make test`, `CLAUDE.md` current, LKH-3 retry, commit.
3. I write the Phase 0 summary and **stop until you say "continue"** to start Phase 1.

Until then I write no code. The only files created so far are documents (`docs/`) and an empty git repository.
