# PHASE 4: Solvers

Status: **complete, waiting for your "continue"** to start Phase 5 (experiments and statistics). Date: 2026-10-23 (local time).
Plan window: 2.2 days (Oct 24 to Oct 26). Exit gate from the plan: the correctness report with evidence (this file). **Met**: 12 of 12 checks passed, run at clean commit `07f63ea`.

## What was done
1. **Solvers** (`genomeasm/solvers/`): `greedy` and `nearest_neighbour` ([greedy.py](../../genomeasm/solvers/greedy.py)); `tabu_search` ([tabu.py](../../genomeasm/solvers/tabu.py)); `alns` with 3 destroy and 2 repair operators, calibrated annealing and 5 policies ([alns.py](../../genomeasm/solvers/alns.py), [bandit/policies.py](../../genomeasm/bandit/policies.py)); `lkh` ([lkh.py](../../genomeasm/solvers/lkh.py)); `brute_force` and `held_karp` ([exact.py](../../genomeasm/solvers/exact.py)); O(1) `swap_delta` / `insert_delta` ([moves.py](../../genomeasm/solvers/moves.py)); shared `Problem`, `Budget`, `SolveResult`; `warmup()` (compiles every numba kernel before anything is timed).
2. **LKH-3 3.0.6 built from your tarball** in Ubuntu (D-033) and wrapped: it writes the TSPLIB file and the parameter file, runs `LKH`, and **verifies every tour** it gets back (Hamiltonian cycle, cost recomputed from the cost matrix, equality with LKH-3's own reported length).
3. **Correctness validation** ([solver_validation.py](../../genomeasm/solver_validation.py), [configs/phase4_correctness.yaml](../../configs/phase4_correctness.yaml)); results in [results/phase4/correctness.json](../../results/phase4/correctness.json).
4. Commands: `genomeasm solve ACC --coverage C --solver greedy|nn|tabu|alns|lkh` and `genomeasm solvers validate`.
5. Decisions D-033 to D-038 and 11 glossary terms.

## The four mandatory tests of the master prompt

| Requirement | Evidence | Result |
|---|---|---|
| Every heuristic reaches the exact optimum on small instances, or its gap is documented | 72 small instances with the optimum from Held-Karp (60 random overlap matrices with 6 to 14 reads, 12 windows of the real phiX174 genome with 8 to 15 reads); table below | tabu 71 of 72, LKH-3 12 of 12, ALNS 72, 72, 72, 71, 71 of 72 (limit set in advance: 95%); greedy and nearest neighbour gaps documented |
| Every returned solution is a valid permutation | 15 solvers x 72 instances, each reported fitness recomputed independently | pass |
| Objective values are reproducible per seed | tabu and the 5 ALNS policies run twice with the same seed and an evaluation budget: identical order, fitness, evaluations and operator choices | pass (6 of 6) |
| Delta evaluation equals a full recomputation | **160,000 random swap and insert moves** on 8 instances, plus an exhaustive test of every move on instances with 4 and 6 reads, plus a property test of random move sequences | 0 mismatches |

Also: brute force and Held-Karp agree on 60 random instances (0 disagreements).

## Evidence

### Small instances with a known optimum (72 instances; budget 500,000 evaluations for tabu and the five policies, 100,000 for single-operator runs, LKH-3: 1 s on 12 instances)

| Solver | runs | found the optimum | success rate | mean gap % | max gap % |
|---|---|---|---|---|---|
| greedy | 72 | 20 | 0.278 | 5.3385 | 19.485 |
| nearest_neighbour | 72 | 15 | 0.208 | 12.4441 | 36.612 |
| tabu | 72 | 71 | 0.986 | 0.0045 | 0.325 |
| alns_random | 72 | 72 | 1.000 | 0.0000 | 0.000 |
| alns_roundrobin | 72 | 72 | 1.000 | 0.0000 | 0.000 |
| alns_roulette | 72 | 72 | 1.000 | 0.0000 | 0.000 |
| alns_epsilon_greedy | 72 | 71 | 0.986 | 0.0102 | 0.735 |
| alns_ucb1 | 72 | 71 | 0.986 | 0.0045 | 0.325 |
| alns single operator: random+greedy | 72 | 63 | 0.875 | 0.2211 | 3.395 |
| alns single operator: random+regret | 72 | 61 | 0.847 | 0.2085 | 3.209 |
| alns single operator: worst+greedy | 72 | 60 | 0.833 | 0.2504 | 3.268 |
| alns single operator: worst+regret | 72 | 57 | 0.792 | 0.4605 | 5.263 |
| alns single operator: segment+greedy | 72 | 67 | 0.931 | 0.0748 | 1.948 |
| alns single operator: segment+regret | 72 | 67 | 0.931 | 0.0325 | 0.629 |
| lkh3 | 12 | 12 | 1.000 | 0.0000 | 0.000 |

(Gap = 100 x (optimum - F) / optimum.) The two constructive baselines are far from optimal on these instances, as expected. Every single operator used alone is clearly worse than the adaptive mixtures: this is the effect that operator selection is supposed to exploit.

### A look at a realistic instance (one instance, one seed, 10 s each: a smoke comparison, **not a result**)
phiX174, coverage 10, 1% substitution errors, tolerant overlaps (5%, minimum overlap 20): 536 reads, 3,961 overlap edges, F of the true order 43,311.

| Solver | F | L | evaluations | wall seconds |
|---|---|---|---|---|
| greedy | 44829 | 8771 | 3,961 | 0.01 |
| nearest_neighbour | 43912 | 9688 | 287,296 | 0.00 |
| tabu | 45154 | 8446 | 104,549,437 | 10.02 |
| alns_random | 45239 | 8361 | 210,238,008 | 10.00 |
| alns_roundrobin | 45247 | 8353 | 213,432,659 | 10.00 |
| alns_roulette | 45221 | 8379 | 204,908,478 | 10.00 |
| alns_epsilon_greedy | 45242 | 8358 | 185,831,836 | 10.00 |
| alns_ucb1 | 45227 | 8373 | 177,306,023 | 10.00 |
| lkh3 | 45203 | 8397 | n/a | 12.05 |

Read this carefully: with one seed the differences between the ALNS variants (45,221 to 45,247) are tiny and mean nothing yet; Phase 5 has the 30 seeds and the tests. What is solid: every solver found a valid order, all local searches are at least as good as their greedy start (by construction), and **the true order (43,311) scores below every heuristic**: maximising F does not mean recovering the genome (D-036). LKH-3 used 12.05 s of wall-clock for `TIME_LIMIT = 10` (D-035) and, with this limit, did not beat the ALNS runs; a longer LKH-3 run (and the best-known value across all solvers) will be the reference in Phase 5.

### Speed (this laptop)
Tabu about 10 million evaluations per second; ALNS about 18 to 21 million evaluations per second (an ALNS iteration costs many evaluations). Compilation is done once by `warmup()` (under a second once cached).

### Tests and checks
`make check` exit code 0: ruff, ruff format, mypy and **pytest 345 passed** (282 earlier + 63 new: 38 solver tests, 19 policy tests, 6 on the committed results). Coverage of `bandit` is 100%; the coverage tool cannot see inside numba-compiled kernels, so the reported 35 to 49% for `tabu`, `alns`, `greedy` and `moves` is an artefact, and their correctness rests on the exhaustive and property tests above.

## Problems found and how they were handled (R10)
1. **Error-free instances are trivial** (D-036). On the real phiX174 reads without errors, greedy, LKH-3 and the true order tie in all 5 settings tried. A benchmark there would show nothing. The comparison needs noisy instances.
2. **First timing run was dominated by the compiler** (8 to 11 s of numba compilation inside a 3 s budget). Fixed with `warmup()`; experiments must call it first.
3. **The ALNS rewarded the wrong thing in its first version.** Moves that keep the fitness were so frequent that UCB1 preferred the arm that produced the most of them (2,022 of 2,291 moves). Lateral moves now earn 0 (D-034).
4. **ALNS failed to reach the optimum on small instances.** About 90% of its small moves rebuilt the same tour; success on 40 small instances was 78% to 93% at 40,000 evaluations. Insertion jitter and a larger removal range helped; the remaining failures were all instances with 14 reads and disappeared at 1,000,000 evaluations, so the validation budget was set to 500,000 before the validation run. The diagnostics used other instance sets than the validation. Defaults are not tuned values (D-034).
5. **Two unit tests failed at first** (the tolerant-overlap and single-operator tests), because the calibration moves used up the whole small budget of the test. Calibration is now limited to 10% of the budget.
6. **The ALNS (every policy) is weaker than tabu on tiny random instances at small budgets**, and tabu missed one of 72 at 500,000 evaluations (gap 0.33%). Both are reported as measured, not hidden.

## Not verified (honest gaps)
- **Identity of the LKH-3 source with the official release** (D-033): the archive is version 3.0.6 and was supplied by you, but the official page was not reachable to compare.
- Citations: Ropke and Pisinger (2006), Auer, Cesa-Bianchi and Fischer (2002), Held and Karp (1962), Helsgaun (2017) and the NP-hardness of the ATSP are `TODO-VERIFY`.
- LKH-3 results depend on the machine and on its parameters; defaults plus `MAX_TRIALS` large were used, no tuning. It is not reproducible by evaluation count.
- The "random matrix" instances are synthetic and are used only to test algorithms. They support no biological claim.
- The preview is one instance and one seed.
- Tabu and ALNS parameters are untuned defaults.

## Deviations from the plan
- None in scope. Extras: single-operator runs and the preview table in the correctness results (the plan asks for the single-operator variants in Phase 4, item 5).
- The LKH-3 build used WSL as agreed with you (MinGW was not tried).
- Phase 4 took about one calendar day instead of the planned 2.2.

## Facts that matter for Phase 5
- API: `Problem.from_instance(inst)`; `greedy(p)`, `nearest_neighbour(p)`, `tabu_search(p, Budget(seconds=..), seed=..)`, `alns(p, Budget(..), seed=.., policy="ucb1", fixed_arm=None)`, `lkh(p, seconds=.., seed=..)`; every `SolveResult` has `order`, `fitness`, `merged_length`, `evaluations`, `seconds`, `trace` (best-so-far (seconds, evaluations, F)), `info` (ALNS: operator outcomes per arm, arm history, probability trace).
- Call `warmup()` once per process before timing. On Windows use `if __name__ == "__main__"` guards and spawn-safe workers for multiprocessing.
- Use noisy instances (D-036). Parameters to tune on tuning instances: tabu tenure, moves, restart; ALNS destroy fractions, noise, reward vector, window, epsilon, segment, reaction, final temperature ratio.
- LKH-3 wall-clock exceeds `TIME_LIMIT` by its preprocessing (about 2 s for 536 reads, more for larger n). Report total time.

## Open risks
- **R-02 (LKH-3)** is resolved for building and running; its provenance caveat (D-033) remains. LKH-3 run time grows with n and with each call starting Ubuntu (about 0.5 s overhead).
- R-06 (laptop speed): 30 seeds x many configurations on 8 threads; plan in Section 7 of the plan: at most 6 workers. A 10 s run with LKH-3 takes 12 s.
- R-01 (schedule): about 2 days ahead of the plan after Phase 4.
- Q3 and Q5 (slide software, presenting laptop) still unanswered; neither blocks Phase 5.

## Test output
```
345 passed in 62.15s     (make check: ruff, ruff format --check, mypy, pytest all exit 0)
12 checks, 0 failed      (genomeasm solvers validate, 611 s)
```
