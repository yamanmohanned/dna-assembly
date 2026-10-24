# Phase 5: experiments and statistics

Status: complete, waiting for your "continue" (Phase 6: Streamlit app). Dates: 2026-10-23 to
2026-10-24. All numbers below are read from `results/phase5/` by a script; the claims are guarded by
`tests/test_phase5_results.py`.

## 1. The short version

**In plain words.** We made 300 different "puzzles" out of real virus DNA (cut into short
overlapping reads with typing mistakes), and gave every method the same 8 seconds on each puzzle. The
method called ALNS keeps breaking the current answer and mending it with one of six repair tools.
The question was whether a *learning* chooser (a "bandit") picks the tool better than a dice roll.
The strongest competitor, LKH-3, is a famous solver for this kind of puzzle.

**What we found, honestly.**

1. **The pre-registered claim is not supported.** UCB1 against random tool choice: median paired
   difference +0.000 pp, Holm-corrected p = 0.90, A12 = 0.521
   (0.5 means no difference). Epsilon-greedy is even slightly *worse* than random (p < 0.0001,
   A12 = 0.462).
2. **We found and diagnosed a defect in our UCB1** (decision D-046): once nothing improves any more
   it always picked the first tool. With ties broken at random (a post-hoc fix, reported separately)
   UCB1 is slightly better than random selection (A12 = 0.567,
   p = 0.0002), but not different from the classical roulette
   wheel (p = 0.37). The effect is small and exploratory.
3. **ALNS clearly beats greedy and tabu search**, but **LKH-3 was never beaten**: in all 300 read
   sets no other solver reached a higher fitness than LKH-3 (checked by this script; LKH-3 also uses
   the whole 8 s, its median time is 7.97 s).
4. All solvers produce assemblies of the same *biological* quality (genome fraction about
   0.996, no misjoins); the fitness differences are small in
   absolute terms (the best-known gap of greedy is 0.88% on average).

## 2. What was run

| Item | Value |
|---|---|
| Read sets | 300 = 10 instance settings x 30 seeds (lambda windows of 5,000 and 10,000 bases and the rest of phiX174; coverage 5x, 10x, 20x; substitution errors 0.5%, 1%, 2%; read length 100 and 60) |
| Overlaps | error-tolerant (tolerance 0.05), minimum overlap 20 (15 for 60-base reads), decision D-036 |
| Budget | 8 s wall-clock per solver and read set (greedy and nearest neighbour take milliseconds); LKH-3 under the same total wall-clock (D-041) |
| Solvers | greedy, nearest neighbour, tabu, tabu + candidate lists, ALNS with 5 policies (random, round robin, roulette, epsilon-greedy, UCB1), 6 single-operator ALNS, LKH-3, plus the post-hoc UCB1 with random ties |
| Tuning | random search, 24 trials per algorithm on 6 tuning instances x 3 seeds (phiX174 positions 0 to 2500, seeds 1000 and up, no overlap with the evaluation), 5 s per run (D-040) |
| Scaling | lambda windows of 2,500 to 20,000 bases (n = 249 to 1986), 10 seeds |
| Statistics | gap to the best-known fitness of the read set; paired Wilcoxon, Holm within each family, bootstrap CI (10000 resamples), A12, Cliff's delta, win rate (D-042) |

"Best known" is the best fitness any of our solvers reached on that read set. It is **not** a proven
optimum (exact solvers only work for tiny instances, Phase 4), so every gap below is a gap to the best
solution *we found*.

## 3. Results

### 3.1 Final quality after the same budget

Gap = 100 x (best known - F) / best known on the same read set; lower is better. Figures:
`fig05_final_fitness`, `fig04_convergence`, `fig11_pareto`.


| Solver | Mean gap % | Median gap % | Q1 | Q3 | Median seconds |
|---|---|---|---|---|---|
| LKH-3 | 0.000 | 0.000 | 0.000 | 0.000 | 7.974 |
| ALNS roulette wheel | 0.174 | 0.087 | 0.000 | 0.238 | 8.001 |
| ALNS UCB1, random ties (post hoc) | 0.179 | 0.038 | 0.000 | 0.239 | 8.001 |
| ALNS round robin | 0.183 | 0.079 | 0.000 | 0.238 | 8.002 |
| ALNS random | 0.196 | 0.122 | 0.000 | 0.274 | 8.001 |
| ALNS UCB1 | 0.204 | 0.088 | 0.000 | 0.289 | 8.001 |
| single operator: segment+regret | 0.251 | 0.162 | 0.016 | 0.378 | 8.002 |
| single operator: segment+greedy | 0.253 | 0.163 | 0.008 | 0.406 | 8.001 |
| ALNS epsilon-greedy | 0.261 | 0.148 | 0.000 | 0.366 | 8.002 |
| single operator: worst+greedy | 0.267 | 0.105 | 0.000 | 0.390 | 8.001 |
| single operator: worst+regret | 0.269 | 0.105 | 0.000 | 0.362 | 8.001 |
| single operator: random+regret | 0.349 | 0.204 | 0.013 | 0.552 | 8.001 |
| single operator: random+greedy | 0.365 | 0.207 | 0.019 | 0.567 | 8.001 |
| Tabu + candidate lists (post hoc) | 0.441 | 0.154 | 0.002 | 0.578 | 8.004 |
| Tabu (swap + insert) | 0.593 | 0.387 | 0.095 | 0.803 | 8.004 |
| Greedy | 0.876 | 0.449 | 0.096 | 1.106 | 0.003 |
| Nearest neighbour | 4.942 | 3.701 | 2.176 | 8.544 | 0.002 |


The median is 0 for LKH-3 and for many ALNS runs because on many read sets (low error, low coverage)
the best solution is found by several solvers; the *mean* gap shows the hard read sets, so both are
shown. Per instance setting (mean gap %):


| Instance | n reads (median) | Greedy | Tabu | Tabu + candidate lists | ALNS random | ALNS roulette wheel | ALNS UCB1 | ALNS UCB1, random ties | LKH-3 |
|---|---|---|---|---|---|---|---|---|---|
| L5-c5-e1.0 | 249.000 | 0.096 | 0.086 | 0.016 | 0.013 | 0.006 | 0.021 | 0.002 | 0.000 |
| L5-c10-e1.0 | 496.500 | 0.388 | 0.388 | 0.127 | 0.077 | 0.059 | 0.042 | 0.015 | 0.000 |
| L5-c20-e1.0 | 987.000 | 0.774 | 0.763 | 0.543 | 0.228 | 0.213 | 0.273 | 0.215 | 0.000 |
| L5-c10-e0.5 | 491.000 | 0.039 | 0.039 | 0.008 | 0.010 | 0.002 | 0.002 | 0.000 | 0.000 |
| L5-c10-e2.0 | 500.000 | 2.907 | 1.256 | 1.292 | 0.280 | 0.283 | 0.409 | 0.305 | 0.000 |
| L5-c10-e1.0-l60 | 813.500 | 0.990 | 0.573 | 0.416 | 0.329 | 0.260 | 0.273 | 0.239 | 0.000 |
| L10-c5-e1.0 | 498.500 | 0.090 | 0.090 | 0.039 | 0.027 | 0.032 | 0.031 | 0.009 | 0.000 |
| L10-c10-e1.0 | 994.000 | 0.346 | 0.340 | 0.170 | 0.197 | 0.111 | 0.087 | 0.082 | 0.000 |
| L10-c10-e2.0 | 999.000 | 2.743 | 2.022 | 1.744 | 0.777 | 0.764 | 0.886 | 0.922 | 0.000 |
| P-c10-e1.0 | 287.000 | 0.383 | 0.378 | 0.056 | 0.019 | 0.012 | 0.014 | 0.007 | 0.000 |


### 3.2 The pre-registered claim: does bandit selection beat random selection?

Claim (master prompt): *bandit-based adaptive operator selection finds better assemblies than random
operator selection under an equal computational budget.* Test: paired by read set, two primary
comparisons, Holm over the two (alpha = 0.05). Positive differences mean that the first solver is
better. Figure: `fig07_ablation`.

| Comparison (a vs b) | Mean difference in pp [95% CI] | Median difference in pp | p (Holm) | A12 | wins / ties / losses of a |
|---|---|---|---|---|---|
| ALNS UCB1 vs ALNS random | -0.008 [-0.022, +0.006] | +0.000 | 0.90 | 0.521 | 118 / 87 / 95 |
| ALNS epsilon-greedy vs ALNS random | -0.065 [-0.082, -0.048] | +0.000 | <0.0001 | 0.462 | 78 / 75 / 147 |

**Verdict: not supported.** UCB1: not supported (no significant difference). Epsilon-greedy: rejected (random is significantly better); the median paired difference is 0, the effect is small.
In the single instance settings UCB1 was significantly different from random in
4 of 10
settings and not consistently in one direction (A12 above 0.5 in
6 settings);
epsilon-greedy was worse than random in
8 of 10.

Other comparisons (their own families, Holm inside each):

| Comparison (a vs b) | Mean difference in pp [95% CI] | Median difference in pp | p (Holm) | A12 | wins / ties / losses of a |
|---|---|---|---|---|---|
| ALNS round robin vs ALNS random | +0.013 [-0.001, +0.026] | +0.000 | 0.0013 | 0.543 | 130 / 87 / 83 |
| ALNS roulette wheel vs ALNS random | +0.022 [+0.009, +0.034] | +0.000 | <0.0001 | 0.537 | 134 / 93 / 73 |
| ALNS UCB1 vs ALNS roulette wheel | -0.030 [-0.045, -0.015] | +0.000 | 0.0014 | 0.488 | 91 / 97 / 112 |
| ALNS epsilon-greedy vs ALNS roulette wheel | -0.087 [-0.104, -0.070] | -0.015 | <0.0001 | 0.431 | 45 / 86 / 169 |

The two non-learning-but-structured choosers, round robin and the classical roulette wheel, are
slightly better than random (A12 = 0.543 and
0.537); the bandits are not better than the roulette wheel
(UCB1 vs roulette: A12 = 0.488, roulette better).

**Why might this be?** (hypotheses, `TODO-VERIFY`, not tested here) (a) With six similar operator pairs
and short runs there is little to learn: when every pair is used equally the share of moves that
improve the solution is 4.5% to 10.1% across pairs. (b) The reward signal is rare (a move
that improves is a few percent of the moves) and noisy. (c) The tuning set (phiX174, 18 read sets) may
not represent the evaluation set (see 3.3 and 5).

### 3.3 Post hoc: the UCB1 defect (D-046) and UCB1 with random ties

The arm-probability figure (`fig08_arm_probability`) showed that UCB1 moves onto the first operator
pair (random destroy + greedy repair) in the second half of the run. Pooled over all runs this pair got
59.5% of UCB1's uses (alns_random uses every pair 16.7%), although it is the pair with
the *lowest* share of improving moves when pairs are used equally. Cause, reproduced without ALNS
(`tests/test_bandit.py`): once every arm's window holds only zero rewards (a plateau), all UCB1 indices
are equal and `argmax` returns the lowest index. Epsilon-greedy breaks ties at random; UCB1 did not.

The original `alns_ucb1` and its results are kept; `alns_ucb1_rt` (ties at random) was tuned and run
afterwards, so its results are **exploratory**. Own family, Holm inside:

| Comparison (a vs b) | Mean difference in pp [95% CI] | Median difference in pp | p (Holm) | A12 | wins / ties / losses of a |
|---|---|---|---|---|---|
| ALNS UCB1, random ties (post hoc) vs ALNS random | +0.016 [+0.001, +0.031] | +0.000 | 0.0002 | 0.567 | 142 / 86 / 72 |
| ALNS UCB1, random ties (post hoc) vs ALNS UCB1 | +0.024 [+0.011, +0.037] | +0.000 | <0.0001 | 0.544 | 133 / 105 / 62 |
| ALNS UCB1, random ties (post hoc) vs ALNS roulette wheel | -0.005 [-0.019, +0.007] | +0.000 | 0.37 | 0.533 | 116 / 104 / 80 |

UCB1 with random ties is better than random selection in 8 of
10 instance settings and significantly different in 6; pooled it has the
smallest median gap of all ALNS variants (0.038% against 0.122% for random),
but the mean-gap difference is small (0.016 pp) and it is not better
than the roulette wheel. Reading: the original negative result was partly caused by our
tie-breaking, and a correctly implemented bandit is *at best slightly* better than random, about as good
as the classical adaptive method. The pre-registered verdict above does not change.

### 3.4 Comparison with the baselines

| Comparison (a vs b) | Mean difference in pp [95% CI] | Median difference in pp | p (Holm) | A12 | wins / ties / losses of a |
|---|---|---|---|---|---|
| ALNS random vs Tabu + candidate lists (post hoc) | +0.245 [+0.199, +0.293] | +0.028 | <0.0001 | 0.577 | 169 / 66 / 65 |
| ALNS random vs LKH-3 | -0.196 [-0.223, -0.169] | -0.122 | <0.0001 | 0.143 | 0 / 86 / 214 |
| ALNS random vs Greedy | +0.680 [+0.582, +0.781] | +0.311 | <0.0001 | 0.724 | 246 / 54 / 0 |
| Tabu + candidate lists (post hoc) vs LKH-3 | -0.441 [-0.509, -0.375] | -0.154 | <0.0001 | 0.123 | 0 / 74 / 226 |
| Tabu + candidate lists (post hoc) vs Greedy | +0.435 [+0.376, +0.495] | +0.223 | <0.0001 | 0.640 | 246 / 54 / 0 |
| LKH-3 vs Greedy | +0.876 [+0.758, +0.994] | +0.449 | <0.0001 | 0.948 | 269 / 31 / 0 |
| Tabu + candidate lists (post hoc) vs Tabu (swap + insert) | +0.152 [+0.130, +0.175] | +0.131 | <0.0001 | 0.612 | 221 / 54 / 25 |
| ALNS random vs Tabu (swap + insert) | +0.398 [+0.347, +0.448] | +0.205 | <0.0001 | 0.711 | 241 / 54 / 5 |

* **LKH-3 is never beaten** by any solver on any of the 300 read sets and has the smallest gap
  everywhere. Its total wall-clock time was within the budget in most runs, but in
  11% of the runs it exceeded the 8 s budget by up to 0.7 s (the overhead is calibrated, not exact),
  which slightly favours LKH-3 (limitation, D-041).
* **Tabu search (swap + insert)** was a weak baseline: even tuned, it ended with exactly the fitness of
  its greedy start in 193 of 300 read sets. On noisy instances it does help (mean gap
  1.26% against 2.91% for greedy at 2% errors,
  coverage 10x) and it is equal to greedy where greedy is already best. `tabu_cl` (candidate-list inserts,
  D-045, post hoc) is better than `tabu` in 221 read sets and worse in 25, but ALNS is
  still better than `tabu_cl` (ALNS random wins 169, loses
  65). The tuning data alone ranked `tabu_cl` *below* `tabu` (mean regret
  0.491 vs 0.344):
  the tuning set was not representative for this algorithm, a warning for all tuned results.
* **Greedy** (0.88% mean gap) is a remarkably strong start: it already finds the best known
  solution in 10% of the read sets.

### 3.5 Operator ablation

A single operator pair used alone is worse than the random mixture in all six cases
(`fig07_ablation`, orange; all Holm-corrected p below 0.001 in `tests.csv`, family
`random_vs_single_operator`). Mean paired difference of the single operator against the mixture: the two
random-destroy pairs lose 0.153 to 0.169 pp, the other four lose 0.055 to
0.073 pp. So the plain random mixture is already better than any single pair, and a bandit has
to beat that mixture, not a single tool.

### 3.6 Scaling (`fig10_runtime_vs_n`)


| instance | n | build s | greedy s | lkh3 s | tabu gap % | tabu_cl gap % | alns_random gap % | alns_ucb1 gap % |
|---|---|---|---|---|---|---|---|---|
| S2500 | 249.000 | 0.196 | 0.001 | 1.473 | 0.183 | 0.000 | 0.000 | 0.007 |
| S5000 | 497.000 | 0.674 | 0.003 | 1.996 | 0.358 | 0.131 | 0.117 | 0.026 |
| S10000 | 994.500 | 2.564 | 0.010 | 4.276 | 0.352 | 0.132 | 0.198 | 0.135 |
| S20000 | 1986.500 | 10.313 | 0.030 | 12.723 | 0.366 | 0.292 | 0.337 | 0.317 |

Building the overlap graph (error-tolerant overlaps are quadratic in the number of reads) grows from
0.2 s to 10.3 s and dominates the total time
for large instances; LKH-3's total wall-clock with a 1 s search grows from 1.5 s to
12.7 s because of its preprocessing, so at n about 1986 it no longer fits
the 8 s budget (we do not use such sizes in the main run). Times to within 0.1% of the best known are
censored at the budget for most solvers (see `scaling.csv`). Gap columns here are medians over
10 read sets with a 1 s LKH-3 search; they are not comparable to 3.1.

### 3.7 Biological quality

Metrics (D-043), medians over 300 read sets; the last row is the true read order of the same read sets
(the order in which the reads lie on the genome). It is fragmented because reads with errors do not
overlap (D-031) and it is **not** the fitness optimum: the best fitness we found is higher than the
fitness of the true order in 99.3% of the read sets (median 3.2% higher).


| Solver | Genome fraction | Identity | Contigs | N50 | Misjoins | Adjacent pairs correct |
|---|---|---|---|---|---|---|
| Greedy | 0.996 | 0.994 | 39.000 | 444.500 | 0.000 | 0.939 |
| Nearest neighbour | 0.996 | 0.992 | 61.000 | 303.000 | 0.000 | 0.908 |
| Tabu (swap + insert) | 0.996 | 0.994 | 39.000 | 448.500 | 0.000 | 0.939 |
| Tabu + candidate lists (post hoc) | 0.996 | 0.995 | 38.000 | 400.500 | 0.000 | 0.952 |
| ALNS random | 0.996 | 0.995 | 34.000 | 443.500 | 0.000 | 0.947 |
| ALNS round robin | 0.996 | 0.995 | 33.500 | 474.000 | 0.000 | 0.951 |
| ALNS roulette wheel | 0.996 | 0.995 | 34.000 | 461.000 | 0.000 | 0.947 |
| ALNS epsilon-greedy | 0.996 | 0.995 | 33.000 | 462.000 | 0.000 | 0.948 |
| ALNS UCB1 | 0.996 | 0.995 | 34.000 | 469.000 | 0.000 | 0.949 |
| ALNS UCB1, random ties (post hoc) | 0.996 | 0.995 | 33.000 | 477.500 | 0.000 | 0.949 |
| LKH-3 | 0.996 | 0.995 | 29.000 | 464.500 | 0.000 | 0.950 |
| True read order (reference) | 0.997 | nan | 61.750 | 241.000 | 0.000 | 1.000 |


All solvers have the same genome fraction and no misjoins; the stronger solvers have slightly higher
identity, fewer contigs and more correct neighbouring pairs. The fitness differences between the
solvers (a fraction of a percent, 3.1) therefore show only weakly in the recovered genome on these
instances; whether the objective is the right one is a question of D-036, not of the solvers.

### 3.8 Tuning outcome

Mean regret (% below the best fitness any trial reached) on the tuning set; "default" is trial 0.
Winner's curse: the best of 24 trials looks better than it is (glossary).


| Algorithm | Default regret % | Tuned regret % | Best trial | Tuned parameters |
|---|---|---|---|---|
| Tabu (swap + insert) | 0.5012 | 0.3441 | 1 | {"tabu":{"tenure_factor":0.5,"moves_factor":5,"restart_factor":10,"kick_fraction":0.05}} |
| ALNS random | 0.1601 | 0.0814 | 11 | {"alns":{"destroy_lo":0.05,"destroy_hi":0.25,"q_cap":15,"noise":0.3,"final_ratio":0.0001,"accept_probability0":0.8,"regret_k":2,"p_worst":2.0,"reward_preset":"balanced"},"policy":{"random":{}}} |
| ALNS round robin | 0.1328 | 0.0885 | 10 | {"alns":{"destroy_lo":0.02,"destroy_hi":0.6,"q_cap":30,"noise":0.3,"final_ratio":0.001,"accept_probability0":0.5,"regret_k":3,"p_worst":5.0,"reward_preset":"balanced"},"policy":{"roundrobin":{}}} |
| ALNS roulette wheel | 0.1032 | 0.0860 | 17 | {"alns":{"destroy_lo":0.1,"destroy_hi":0.15,"q_cap":15,"noise":0.3,"final_ratio":0.001,"accept_probability0":0.5,"regret_k":2,"p_worst":2.0,"reward_preset":"classical"},"policy":{"roulette":{"segment":100,"reaction":0.8}}} |
| ALNS epsilon-greedy | 0.1222 | 0.1105 | 1 | {"alns":{"destroy_lo":0.02,"destroy_hi":0.25,"q_cap":30,"noise":0.3,"final_ratio":0.01,"accept_probability0":0.8,"regret_k":2,"p_worst":2.0,"reward_preset":"progress"},"policy":{"epsilon_greedy":{"window":25,"epsilon":0.05}}} |
| ALNS UCB1 | 0.2323 | 0.1710 | 18 | {"alns":{"destroy_lo":0.05,"destroy_hi":0.4,"q_cap":15,"noise":0.3,"final_ratio":0.0001,"accept_probability0":0.2,"regret_k":2,"p_worst":2.0,"reward_preset":"progress"},"policy":{"ucb1":{"window":200,"c":0.3}}} |
| LKH-3 | 0.0232 | 0.0000 | 1 | {"lkh":{"CANDIDATE_SET_TYPE":"NEAREST-NEIGHBOR","MAX_CANDIDATES":"6"}} |
| Tabu + candidate lists (post hoc) | 0.7998 | 0.4909 | 23 | {"tabu":{"tenure_factor":2.0,"moves_factor":5,"restart_factor":10,"kick_fraction":0.2,"candidate_fraction":0.5,"candidates":6}} |
| ALNS UCB1, random ties (post hoc) | 0.1234 | 0.0967 | 6 | {"alns":{"destroy_lo":0.02,"destroy_hi":0.25,"q_cap":15,"noise":0.3,"final_ratio":0.0001,"accept_probability0":0.2,"regret_k":3,"p_worst":2.0,"reward_preset":"balanced"},"policy":{"ucb1_rt":{"window":50,"c":1.0}}} |


The best ALNS policy on the tuning data was **ALNS random** (it is `alns_best` in
the main comparisons). The tuned UCB1 had the *worst* tuning regret of all ALNS policies; UCB1 with
random ties improved it from 0.171 to
0.097, consistent with the defect.

## 4. What went wrong, and what was done about it

| Incident | What happened | Handling |
|---|---|---|
| Laptop suspend during tuning | 5 task blocks have a run with an impossible clock (one reported 10,513 s) | blocks dropped from the work file and re-run; keep-awake requested (D-040) |
| WSL crashed during the main run | 121 LKH-3 runs and the first scaling run ended in errors | re-run only where missing (stage `rerun-lkh`); error records kept as audit trail (D-044) |
| WSL down again (2026-10-24, after the LKH-3 runs) | `E_UNEXPECTED`; the LKH-3 unit test cannot run | not caused by code; **needs your action** (restart WSL), see section 7 |
| Weak tabu baseline | tabu equal to greedy in 193/300 read sets | tabu_cl added, disclosed, both reported (D-045) |
| UCB1 tie-breaking defect | UCB1 plays arm 0 on a plateau | original kept, ucb1_rt added post hoc (D-046) |
| Stage provenance | the raw run files carry no commit | `results/phase5/stages.json` records the commit of every stage (inferred where marked) |

## 5. Limitations (read before quoting a number)

* "Best known" is the best of our solvers, not the optimum. LKH-3 is a heuristic too (no proof).
* Two algorithms (tabu_cl, UCB1 random ties) were added after seeing the first results; their results
  are exploratory. Choosing which post-hoc variants to add is itself a forking path.
* The tuning set is small (18 read sets, one genome region) and tuned results do not always transfer
  (tabu_cl, UCB1 random ties improved on the evaluation set more than the tuning data suggested).
* All instances use substitution errors only (D-036) and one strand; real sequencing data also has
  insertions, deletions, both strands and repeats.
* Timing: five workers share a 4-core laptop (numba threads = 1); all solvers suffer alike, but the
  absolute times carry noise. LKH-3 exceeded the budget in some runs (3.4).
* Equal wall-clock does not mean equal effort: ALNS and tabu run compiled loops, LKH-3 runs C code
  through WSL; evaluation counts are recorded but not comparable to LKH-3.
* The explanations in 3.2 are hypotheses (`TODO-VERIFY`).
* Effect sizes of all policy comparisons are small (A12 between 0.43 and 0.57); with 300 paired read sets
  tests are very sensitive, so "significant" does not mean "large".

## 6. How to reproduce

```
.venv\Scripts\python.exe -m genomeasm experiments tune      # tuning (resumable)
.venv\Scripts\python.exe -m genomeasm experiments run       # 300 read sets (resumable)
.venv\Scripts\python.exe -m genomeasm experiments scale     # scaling study
.venv\Scripts\python.exe -m genomeasm experiments rerun-lkh --target main    # only where LKH-3 has no result
.venv\Scripts\python.exe -m genomeasm experiments add-runs --target main --algorithm tabu_cl
.venv\Scripts\python.exe -m genomeasm experiments analyze   # statistics, tables, figures 4 to 11
```

The raw runs are `results/phase5/{tuning,main,scaling}_runs.jsonl.gz`; the tables are in
`results/phase5/analysis/`; figures 4 to 11 are in `reports/figures/` (SVG, PDF, PNG, caption).

## 7. Open points and your decisions

* **WSL**: please restart it (reboot, or an administrator PowerShell: `Restart-Service LxssManager`),
  then run `make check`; the LKH-3 unit test needs it. Nothing else in Phase 5 needs WSL any more.
* Q3 (program that opens the .pptx) and Q5 (is this the presenting laptop) are still open.
* Question for the report: do you want the post-hoc variants shown in the main slides, or only in the
  appendix? (Recommendation: main results with the pre-registered algorithms, post-hoc in a clearly
  labelled box.)

