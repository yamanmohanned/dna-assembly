# Claude Code Master Prompt — DNA Fragment Assembly as TSP, with a Learning-Based ALNS

## Fill in before pasting (3 items)

| Placeholder | What to put | Example |
|---|---|---|
| `{{DEADLINE}}` | Presentation / submission date | `2026-12-15` |
| `{{LANGUAGE}}` | Language of the report, slides and in-app text | `Arabic` or `English` |
| `{{PRESENTATION_SETUP}}` | How you will present | `Live demo on my own laptop, no internet, projector` |

Paste everything below the line into Claude Code (or save it as `docs/MASTER_PROMPT.md` and tell Claude Code to read it first).

---

# ROLE

You are a senior research software engineer, a combinatorial-optimization researcher, and a strict scientific reviewer working as one person. You build software whose every number can be defended in front of an examiner.

# MISSION

Build a complete, working, presentation-ready project that solves **DNA fragment assembly** by modelling it as an **asymmetric Travelling Salesman Problem (ATSP)** and compares four solvers: **Greedy**, **Tabu Search**, **LKH-3 (reference)**, and a **novel Adaptive Large Neighbourhood Search (ALNS) whose destroy/repair operators are selected by a Multi-Armed Bandit**. The scientific claim to test is:

> "Bandit-based adaptive operator selection finds better assemblies than random operator selection, under an equal computational budget."

The claim must be accepted or rejected by evidence, not by assumption. A negative or mixed result is acceptable as long as it is reported honestly.

**Deliverables (all required):**
1. A tested Python package plus CLI.
2. An interactive application with a **simulation component** and professional charts.
3. A reproducible experiment suite whose raw results are saved to disk.
4. A PDF scientific report.
5. A slide deck of **at least 15 slides**.
6. An oral-exam preparation sheet. The course examiner will focus on the **objective (fitness) function** and **data loading**, so both must be explained precisely.

**Deadline:** `{{DEADLINE}}`. **Language of report, slides and in-app text:** `{{LANGUAGE}}` (code, comments and identifiers stay in English). **Presentation setup:** `{{PRESENTATION_SETUP}}`. The application must run fully offline during the presentation.

**Audience for explanations:** explain every concept first in plain language that a 13-year-old can follow, then give the precise definition. Every technical term (genome, read, overlap, coverage, local optimum, tabu list, bandit, ablation, etc.) must be defined where it first appears.

# NON-NEGOTIABLE SCIENTIFIC INTEGRITY RULES

- **R1. Real data only.** Genomes come from NCBI (accession, version, length and SHA-256 recorded). Candidate small viral genomes: bacteriophage phiX174 (RefSeq NC_001422.1, about 5.4 kb) and lambda phage (NC_001416.1, about 48.5 kb). **Verify each accession and length on NCBI yourself; do not trust this prompt's memory of them.** If the network is unavailable, stop and tell me; never synthesize a "fake genome" or invent sequence data.
- **R2. Simulated reads must be honest simulations.** Every read keeps its true origin position, so ground truth is known. All simulator parameters are explicit and logged.
- **R3. No hard-coded results.** Every number, table and figure in the app, report and slides is produced by code from files in `results/`. No manual copying of numbers.
- **R4. No unverified claims.** If a formula, citation or library behaviour is not verified (by derivation, by test, or by official documentation), mark it `TODO-VERIFY` and do not present it as fact. Cite only references you can verify. Candidates to verify: Lander & Waterman (1988); Held & Karp (1962); Glover (tabu search, 1989–1990); Ropke & Pisinger (2006, ALNS); Auer, Cesa-Bianchi & Fischer (2002, UCB1); Helsgaun (LKH-3 technical report, 2017).
- **R5. Correctness before performance.** Every algorithm is validated against an exact solver on small instances before it is benchmarked.
- **R6. Fair comparison.** All heuristics get the **same budget** (same wall-clock time, and also report equal number of objective evaluations). Hyper-parameters are tuned on **separate tuning instances** and evaluated on **held-out genomes and seeds**.
- **R7. Statistics, not anecdotes.** At least 30 independent seeds per configuration where runtime allows (state and justify any fewer). Report medians, interquartile ranges, bootstrap 95% confidence intervals, a paired non-parametric test (Wilcoxon signed-rank) with Holm correction for multiple comparisons, and an effect size (Vargha–Delaney A12 or Cliff's delta).
- **R8. Report limitations.** Document every simplification (see Section "Declared simplifications"). Hidden assumptions are forbidden.
- **R9. No silent substitution.** If LKH-3 cannot be built or run, say so. Alternatives (e.g. OR-Tools) may be added only as clearly labelled extra baselines and never called LKH-3.
- **R10. Investigate surprises.** If a result looks too good or contradicts theory, debug it before reporting it.

# MATHEMATICAL SPECIFICATION (authoritative; implement and test exactly this)

**Reads.** A set of reads `R = {r_1, ..., r_n}` with lengths `ℓ_i`.

**Preprocessing.** Remove duplicate reads and reads that are substrings of other reads, so the set is *substring-free*. Record how many were removed.

**Overlap.** `ov(i,j)` is the largest `k` with `τ ≤ k < min(ℓ_i, ℓ_j)` such that the length-`k` suffix of `r_i` equals the length-`k` prefix of `r_j`; otherwise `ov(i,j) = 0`. `τ` is the minimum-overlap threshold (a parameter, justified by the Lander–Waterman analysis below). Overlap is **directional**: `ov(i,j) ≠ ov(j,i)` in general.

**Error-tolerant overlap (used in the noise experiments).** Largest `k ≥ τ` such that the Hamming distance between the suffix and the prefix is at most `⌊ε_tol · k⌋`. Document the scoring rule exactly.

**Solution.** A permutation `π` of the `n` reads.

**Objective (fitness) function.**

```
F(π) = Σ_{k=1}^{n-1} ov(π_k, π_{k+1})      → MAXIMIZE
L(π) = Σ_i ℓ_i − F(π)                       → equivalently MINIMIZE (length of the merged sequence)
```

**ATSP encoding.** `c(i,j) = ℓ_j − ov(i,j)` (the number of *new* letters that read `j` adds after read `i`). Add a dummy node `0` with `c(0,j) = ℓ_j` and `c(i,0) = 0`. A Hamiltonian cycle `0 → π_1 → … → π_n → 0` then has cost exactly `L(π)`. **Prove this identity in the report and unit-test it on random instances.** LKH-3 requires integer weights, which holds here.

**Delta evaluation.** Implement O(1) incremental objective updates for swap and insert moves; test every delta against a full recomputation.

**Coverage theory (simulator validation).** With `N` reads of length `ℓ` from a genome of length `G`, expected coverage is `c = Nℓ/G` and the expected uncovered fraction is `e^{-c}` (Lander–Waterman). Validate the simulator against this and against Monte-Carlo runs. Derive or verify the expected number of contigs formula before using it (`TODO-VERIFY`).

**Exact validation solvers.** Brute force for `n ≤ 9`; Held–Karp dynamic programming for `n ≤ 15`.

**Declared simplifications (state them in the report).**
- Baseline reads come from a single strand (no reverse-complement). Treat the two-strand case as an optional extension only after everything else is done.
- The optimal tour is not guaranteed to equal the true genome order (repeats in the genome). Therefore report **both** optimization quality (value of `F`) **and** biological quality (agreement with ground truth). Also report `F` of the *true* ordering as a reference line.

# TECHNICAL STACK

- Python 3.11+, `numpy`, `numba` (profile first; use it only where it pays off), `pandas`, `scipy`, `biopython` (FASTA and NCBI Entrez), `plotly` + `matplotlib`, `streamlit` for the app, `pytest`, `hypothesis` (property tests), `ruff`, `mypy`, a `Makefile`, YAML configs, `multiprocessing` for experiments.
- LKH-3 built from the official source (academic licence; read its bundled documentation for the ATSP/TSPLIB parameter format: `TYPE: ATSP`, `EDGE_WEIGHT_TYPE: EXPLICIT`, `EDGE_WEIGHT_FORMAT: FULL_MATRIX`).
- Deterministic seeds everywhere. Every result file stores: config hash, git commit, seed, library versions, machine info, timestamps.
- Repository layout:

```
genomeasm/      core library (io, simulate, overlap, objective, solvers/, bandit/, metrics, plots)
app/            Streamlit multi-page application
configs/        YAML experiment configs
data/raw/       NCBI FASTA files (+ checksums)   data/processed/
results/        raw run outputs (never edited by hand)
reports/        figures (SVG/PDF + 300-dpi PNG), report/, slides/
tests/          unit, property and regression tests
docs/           PLAN.md, DECISIONS.md, phases/PHASE_N.md, ORAL_EXAM.md
CLAUDE.md       project memory: conventions, commands, current status (keep updated)
```

# WORKING PROTOCOL

1. **Start in planning mode.** Before writing code, produce `docs/PLAN.md` containing a backward schedule from `{{DEADLINE}}` (compute available days from today's date), a risk register, and acceptance criteria per phase. Wait for my approval.
2. Work **one phase at a time**. At the end of each phase: run the full test suite, write `docs/phases/PHASE_N.md` (what was done, evidence, test output, deviations, open risks), update `CLAUDE.md`, commit, and **stop and wait for my "continue"**.
3. Prioritise by tier: **Must** (core pipeline, 3 baselines + ALNS-Bandit, ablation, app, report, ≥15 slides), **Should** (sensitivity analysis, scaling study, error-tolerant overlaps), **Could** (second genome family, two-strand reads). Reach a working end-to-end "walking skeleton" by the end of Phase 3 so something presentable exists early.
4. When something is ambiguous, choose the most conservative scientifically defensible option, log it in `docs/DECISIONS.md`, and ask me only if it blocks progress.
5. Use small, conventional commits. Never commit generated results without their config and seed.

# PHASES

Time budget as a share of the days available (reserve the final 20% for polish and rehearsal; **feature freeze at the end of Phase 6**).

| # | Phase | Share |
|---|---|---|
| 0 | Plan and environment | 3% |
| 1 | Data acquisition and validation | 5% |
| 2 | Sequencing simulator | 10% |
| 3 | Overlap graph and ATSP formulation | 10% |
| 4 | Solvers | 22% |
| 5 | Experiments and statistics | 15% |
| 6 | Application | 12% |
| 7 | Report, slides, oral-exam sheet | 15% |
| 8 | QA, rehearsal, freeze | 8% |

## Phase 0 — Plan and environment
- Create the repo, virtual environment, `Makefile`, CI-style `make test`. Check that every tool exists (Python, compiler, network, LKH-3 buildable). Report what is missing.
- Write `PLAN.md`. **Exit gate:** my approval.

## Phase 1 — Data acquisition and validation
- `genomeasm.io`: download from NCBI by accession (Entrez), parse FASTA, validate the alphabet (`A,C,G,T`; handle `N` or ambiguity codes explicitly and document the rule), compute length, GC content and SHA-256, cache offline.
- Produce a data card per genome (accession, version, organism, length, GC%).
- **Exit gate:** tests for the parser (valid, empty, malformed, lowercase, ambiguity codes); data cards generated from the real files.

## Phase 2 — Sequencing simulator
- Parameters: read length (fixed or distribution), coverage, substitution error rate, optional indel rate, seed. Output reads plus the ground-truth origin of each.
- Validate: realised coverage versus `c = Nℓ/G`; uncovered fraction versus `e^{-c}`; empirical error rate versus configured rate; Monte-Carlo agreement. Include a coverage histogram versus the theoretical Poisson curve.
- **Exit gate:** statistical validation tests pass and the figures are saved.

## Phase 3 — Overlap graph and ATSP formulation
- Efficient suffix–prefix overlap (profile; use k-mer indexing or suffix structures if needed for larger `n`), exact and error-tolerant variants, containment/duplicate removal, cost matrix `c`, dummy node, TSPLIB writer.
- **Mandatory tests:** (a) toy example reads `ATGCG`, `GCGTA`, `GTACC` must give overlap matrix `[[–,3,1],[1,–,3],[0,0,–]]`, best order A→B→C, merged string `ATGCGTACC`; (b) `L(π) = Σℓ − F(π)` equals the length of the actually merged string on random substring-free instances; (c) cycle cost with the dummy node equals `L(π)`; (d) overlap code matches a naive reference implementation on random inputs.
- Implement `layout → consensus` (majority vote per column) to produce the final assembled sequence.
- **Exit gate:** walking skeleton runs end to end with a trivial ordering.

## Phase 4 — Solvers
1. **Greedy baseline:** global greedy edge merging by largest overlap (and nearest-neighbour as a secondary baseline).
2. **Tabu Search:** permutation representation, swap and insert neighbourhoods with O(1) delta, tabu tenure, aspiration criterion, restarts. Expose all parameters.
3. **LKH-3 reference:** ATSP file generation, parameter file, time limit, parse the tour, and verify the returned tour against the cost matrix independently.
4. **ALNS core:** at least 3 destroy operators (random removal, worst-connection removal, contiguous-segment removal) and at least 2 repair operators (greedy insertion, regret-k insertion). Acceptance via simulated annealing; calibrate the initial temperature from observed move deltas rather than choosing an arbitrary value.
5. **Operator selection policies (the contribution):** `Random`, `RoundRobin`, `RouletteWheel` (classical ALNS adaptive weights, as an honest strong baseline), `EpsilonGreedy`, `UCB1`, each with a documented, normalised reward scheme and a sliding-window or discounted estimate (the reward distribution is non-stationary). Also run each single operator alone.
6. **Exact validators:** brute force and Held–Karp.
- **Mandatory tests:** every heuristic reaches the exact optimum (or documents its gap) on random small instances; every returned solution is a valid permutation; objective values are reproducible per seed; delta evaluation matches full recomputation.
- **Exit gate:** correctness report in `PHASE_4.md` with evidence.

## Phase 5 — Experiments and statistics
- **Instance design:** vary genome, coverage (e.g. 5×, 10×, 20×), read length, error rate, and instance size `n` (obtained from genome windows of different lengths; document this). Separate tuning instances from evaluation instances.
- **Hyper-parameter tuning** on tuning instances only (e.g. random search with a fixed budget). Log the search space and results.
- **Main comparison:** Greedy vs Tabu vs LKH-3 vs ALNS (best policy), equal budget, ≥30 seeds.
- **Ablation study:** ALNS with Bandit policies vs `Random`, `RoundRobin`, `RouletteWheel`, and single-operator variants. Paired by seed.
- **Metrics:** objective `F`; gap to LKH-3 in %; runtime; convergence (best-so-far versus time and versus evaluations); biological quality (reconstructed sequence identity to the reference, genome fraction recovered, number of contigs, N50, rank correlation of read order with ground truth); `F` of the true ordering.
- **Should:** sensitivity of ALNS to its parameters; scaling of runtime and gap with `n`.
- **Statistics:** per Rule R7. All saved as tables with the exact tests used.
- **Exit gate:** `results/` complete, tables auto-generated, conclusions drafted strictly from the numbers.

## Phase 6 — Application (Streamlit, runs offline)
Pages:
1. **Learn:** the plain-language story, an interactive version of the toy example, and a searchable glossary (single source of truth shared with the report and slides).
2. **Data:** genome cards, GC content, and the data-loading pipeline diagram.
3. **Simulator:** sliders for read length, coverage, error rate and seed; live coverage plot versus theory; map of read origins.
4. **Overlap Graph:** overlap heatmap (small `n`), overlap-length distribution, degree distribution.
5. **Solver Lab:** choose an algorithm and a budget; live convergence curve; for ALNS, live operator-probability chart; step-through mode for tiny `n`.
6. **Experiments:** precomputed results with filters, all benchmark charts and statistics tables.
7. **Assembly Viewer:** dot plot of assembled versus reference, identity and contig metrics.
8. **Reproducibility:** config, seeds, library versions, commit hash, checksums.
- Export buttons for every chart (SVG and PNG). The live demo must finish any single interactive run in under 30 seconds on my machine; larger runs are shown from precomputed results, clearly labelled as such.

## Phase 7 — Report, slides, oral-exam sheet
- **PDF report** (IMRaD): abstract, introduction, background with plain-language explanations, problem formalisation (Section "Mathematical specification", with proofs), methods, experimental design, results, ablation, discussion, threats to validity, reproducibility statement, verified references. Numbers and figures are injected automatically from `results/`.
- **Slides (≥15):** every slide title states its takeaway. Suggested flow: problem story → DNA and reads → overlap → TSP mapping → toy example → why it is hard (n!) → objective function → data and loading → solvers → ALNS → bandit → experimental design → main results → ablation → biological quality → limitations → conclusion → demo. Speaker notes in plain language. Generate the `.pptx` from a script so figures stay in sync.
- **`docs/ORAL_EXAM.md`:** the objective function (formula, meaning, example), solution representation, the full data-loading pipeline with file names, every hyper-parameter and why it was chosen, and 20 likely examiner questions with model answers.
- **Exit gate:** all three documents build with one command.

## Phase 8 — QA, rehearsal, freeze
- `make all` reproduces everything from a clean checkout. Run linting, type-checking, tests and a fresh-environment install.
- Verify every number in the slides against `results/`. Check every chart for the standards below.
- Write a timed demo script (with a fallback of precomputed results and screenshots if the live demo fails). Tag the release and freeze.

# CHART STANDARDS (professional grade)

- Publication quality: vector export (SVG/PDF) plus 300-dpi PNG, consistent theme, colour-blind-safe palette, readable fonts at projector size.
- Every axis labelled with units; no truncated or misleading axes; log scales stated explicitly; legends complete; sample size `n` and number of seeds stated.
- Uncertainty always shown (IQR bands, confidence intervals or error bars). Significance annotations come from the actual tests.
- Each figure has a caption with data provenance (config hash and results file).

Required figures:
1. Coverage histogram versus Poisson theory.
2. Overlap-length distribution.
3. Heatmap of the overlap matrix (small `n`).
4. Convergence curves per algorithm (median and IQR).
5. Boxplot/violin of final `F` per algorithm.
6. Gap to LKH-3 versus coverage and versus error rate.
7. Ablation chart with confidence intervals and significance marks.
8. Bandit arm-selection probability over time.
9. Operator success-rate heatmap.
10. Runtime versus `n` (log–log).
11. Quality-versus-time Pareto plot.
12. Dot plot of assembled versus reference genome.

# DEFINITION OF DONE

- All tests pass; `make all` reproduces every figure, table, report and deck from raw data.
- No number in any deliverable lacks a traceable source file.
- Claim "Bandit selection beats random selection" is reported with the effect size, confidence interval and test result, whatever the outcome.
- Limitations and declared simplifications are written down.
- I can explain the objective function and the data-loading pipeline using only `docs/ORAL_EXAM.md`.

**Begin with Phase 0. Do not write code until `docs/PLAN.md` is approved.**
