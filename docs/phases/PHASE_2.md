# PHASE 2: Sequencing simulator

Status: **complete, waiting for your "continue"** to start Phase 3 (overlap graph and ATSP formulation). Date: 2026-10-22.
Plan window: 1.0 day (Oct 22 to Oct 23); finished on Oct 22.
Exit gate from the plan: statistical validation tests pass and the figures are saved. **Both met**: 93 of 93 checks pass, figures are in `reports/figures/`.

## What was done
1. **Simulator** ([genomeasm/simulate.py](../../genomeasm/simulate.py)): reads from a real `Genome`, one strand, with a read-length model (fixed, uniform, truncated normal), substitution, insertion and deletion errors, two seeded random streams, and the ground truth of every read (true start, covered span, number of each error). Model definition: decision D-023.
2. **Theory** ([genomeasm/coverage.py](../../genomeasm/coverage.py)): coverage depth, uncovered fraction, true contigs, the Lander-Waterman (LW) formulas, an **exact finite-genome model** for this simulator, and a derived **exact expected-contig formula** (D-021).
3. **Validation suite** ([genomeasm/sim_validation.py](../../genomeasm/sim_validation.py)) driven by [configs/phase2_validation.yaml](../../configs/phase2_validation.yaml), where every tolerance is written before the run. It makes 93 checks and writes [results/phase2/validation.json](../../results/phase2/validation.json) with full provenance (config hash, git commit and dirty flag, seed, library versions, machine, timestamp).
4. **Figures** ([genomeasm/plots.py](../../genomeasm/plots.py)): `fig01a_coverage_histogram` (required figure 1) and `figS1_validation_curves`, each as SVG, PDF and 300-dpi PNG with a caption file.
5. **Provenance module** ([genomeasm/provenance.py](../../genomeasm/provenance.py)), reused by every later phase.
6. **CLI**: `genomeasm simulate run | validate | figures`.
7. **Glossary**: 14 new terms (read, coverage, Poisson, Lander-Waterman, edge effect, contig, ground truth, error types, seed, Monte Carlo, standard error, provenance).

## Evidence

Everything below is produced by `genomeasm simulate validate` and copied from `results/phase2/validation.json` (200 independent simulated read sets per genome and coverage, seed 20261022, run at clean commit `975eaa7`). Read length 100, contig overlap threshold theta = 30.

### Coverage and contigs on the two real genomes

| Genome | c | reads N | uncovered: simulation (+-SE) | exact model | Lander-Waterman | contigs: simulation (+-SE) | exact formula | Lander-Waterman |
|---|---|---|---|---|---|---|---|---|
| phiX174 | 1 | 54 | 3.663e-01 (+-2.3e-03) | 3.664e-01 | 3.669e-01 | 26.68 (+-0.18) | 26.67 | 26.77 |
| phiX174 | 5 | 269 | 1.228e-02 (+-5.3e-04) | 1.271e-02 | 6.776e-03 | 8.39 (+-0.14) | 8.24 | 8.15 |
| phiX174 | 10 | 539 | 3.626e-03 (+-1.8e-04) | 3.486e-03 | 4.506e-05 | 1.35 (+-0.04) | 1.39 | 0.49 |
| phiX174 | 20 | 1,077 | 1.651e-03 (+-8.8e-05) | 1.642e-03 | 2.069e-09 | 1.00 (+-0.00) | 1.00 | 0.00 |
| lambda | 1 | 485 | 3.680e-01 (+-7.4e-04) | 3.678e-01 | 3.679e-01 | 239.98 (+-0.51) | 239.69 | 240.85 |
| lambda | 5 | 2,425 | 7.171e-03 (+-1.5e-04) | 7.406e-03 | 6.739e-03 | 71.50 (+-0.48) | 71.71 | 73.24 |
| lambda | 10 | 4,850 | 4.340e-04 (+-2.1e-05) | 4.350e-04 | 4.542e-05 | 5.22 (+-0.13) | 5.13 | 4.42 |
| lambda | 20 | 9,700 | 1.898e-04 (+-1.0e-05) | 1.858e-04 | 2.063e-09 | 1.00 (+-0.01) | 1.01 | 0.01 |

(Eight coverages from 1 to 20 were checked; the table shows four. All eight are in the results file and in Figure S1.)

How to read it: the simulation agrees with the **exact model** at every coverage (all within 4 standard errors). The classic Lander-Waterman `exp(-c)` is **not** a good prediction for these short genomes at medium and high coverage: it says 2e-9 of phiX174 is uncovered at c = 20, the real answer is 1.6e-3. The reason is the **edge effect**: a read must fit inside the genome, so the first and last stretch of about one read length are covered less. The exact model captures this and the simulator reproduces it. This is a property of linear single-strand sampling (D-017), not a defect of the simulator.

### Lander-Waterman where it applies: a long genome (1,000,000 positions, 100 runs per point)

| c | uncovered, simulation (+-SE) | exp(-c) | gap |
|---|---|---|---|
| 1 | 3.6764e-01 (+-2.4e-04) | 3.6788e-01 | -0.06 % |
| 2 | 1.3564e-01 (+-2.2e-04) | 1.3534e-01 | +0.23 % |
| 5 | 6.7681e-03 (+-4.9e-05) | 6.7379e-03 | +0.45 % |

All within the pre-set limit (4 standard errors plus 2 % relative). So the simulator does follow Lander-Waterman when the genome is long compared with a read.

### Depth histogram (required figure 1)
For each genome, coverages 5, 10, 20 and also 1, 2, 3, 7, 15: the whole depth histogram agrees with the exact Binomial mixture. Largest |z| over all bin groups of the 16 histogram checks: **3.48** (limit 5). [fig01a_coverage_histogram.png](../../reports/figures/fig01a_coverage_histogram.png) shows simulation, Lander-Waterman Poisson and the exact model together, with the interquartile range across runs.

### Error model (phiX174 at c = 10, 100 runs per setting, 5.39 million reference bases each)

| Configured | Observed (pooled) | Standard error | Bases |
|---|---|---|---|
| substitution 0.001 | 0.00100 | 1.4e-05 | 5,390,000 |
| substitution 0.01 | 0.01004 | 4.3e-05 | 5,390,000 |
| substitution 0.05 | 0.04985 | 9.4e-05 | 5,390,000 |
| deletion 0.005 | 0.00503 | 3.0e-05 | 5,390,000 |
| insertion 0.005 | 0.00498 | 3.0e-05 | 5,390,000 |
| substitution among kept bases 0.01 (with indels) | 0.01003 | 4.3e-05 | 5,362,912 |

Further error checks, all passed:
- Per-run z-scores are calibrated (mean near 0, variance near 1): z mean +0.010, +0.094, -0.162 and z variance 0.90, 0.84, 1.25 for the three substitution rates (limits: mean within 0.40, variance in [0.6, 1.4]).
- **Independent recount**: for substitution-only reads, the Hamming distance between each sampled read and its true genome window equals the recorded substitution count.
- For every read, `length = span - deletions + insertions`.
- For sampled reads with indels, the Levenshtein edit distance to the true window never exceeds the recorded number of edits.

### Read lengths and start positions (200,000 reads on a 1,000,000-position genome)

| Length model | mean observed / expected | chi-square p | KS p (start positions) |
|---|---|---|---|
| fixed | 100.000 / 100.000 | 1.000 | 0.061 |
| uniform_80_120 | 100.029 / 100.000 | 0.744 | 0.893 |
| normal_100_15 | 100.024 / 100.000 | 0.319 | 0.156 |

(Limit: p >= 0.001 for both the length distribution and the uniformity of start positions.)

### Tests and checks
`make check` exit code 0 before the commit: ruff, ruff format, mypy (20 source files), **pytest 204 passed** (126 from earlier phases, 78 new), about 30 s. (The first Phase 2 commit had 202; two SVG line-ending tests were added afterwards.)

New test files and what they prove:

| File | What it proves |
|---|---|
| `test_coverage_theory.py` | hand examples for depth, uncovered fraction, contigs (including a contained read); the exact uncovered fraction, depth histogram and contig formula equal **full enumeration** on tiny cases (fixed and mixed read lengths); mean depth identity; flat interior and 1/S ends of the cover probability; bin grouping |
| `test_simulator.py` | error-free reads are exact substrings at the recorded position (real phiX174); variable lengths stay inside the genome; same seed gives identical reads, other seed does not; changing only an error rate does not move the reads; substitutions, insertions and deletions behave exactly as defined; reads that lose every base are reported; a property test (hypothesis) on random small genomes; CLI |
| `test_sim_validation.py` | edit distance; the histogram check accepts the true model and **rejects a model with 10 % more or fewer reads** (negative control); a reduced run of the whole validation passes; the **committed results**: all 93 checks passed, config equals the current YAML, provenance complete, every required validation present; figures exist in all formats, PNG is 300 dpi, SVG is redrawn byte-identical from the saved results; provenance helpers and git dirty flag |

## Problems found and how they were handled (R10: investigate surprises)
1. **First validation run: 16 of 93 checks failed.** All 16 were the depth-histogram check. I did not touch the tolerance. A 2000-run diagnostic against the exact model showed every depth bin within 1.5 standard errors, and the failing bins were sparse bins that never occurred in the runs (zero spread, so any difference counted as a failure). The check method was changed to merge sparse bins into groups (D-022), a negative-control test was added, and the result is 93 of 93. This is documented as a change made after seeing results.
2. **Floating-point noise in a zero-spread test** (found by my own unit test): a standard error of 1e-17 was treated as real. Spreads below 1e-12 are now treated as zero.
3. **A smoke test with too few runs** failed because the variance-range check needs about 100 runs. The test now keeps the full 100 runs for the error checks and says why.
4. **The textbook value `exp(-c)` disagrees with the simulator on short genomes.** This was investigated, not hidden: it is the edge effect, explained by the exact model and confirmed on a 1,000,000-position genome (gap 0.5 % or less).
5. **Line endings and provenance**: generated files are written with LF; the dirty flag ignores generated outputs, and the committed results were regenerated at the clean commit `975eaa7` (so `dirty: false`). After the results commit, `git ls-files --eol` showed that matplotlib had written the SVG files with Windows line endings. `save_figure` now normalises SVG output to LF, a test fails if a carriage return appears, and the figures were redrawn from the saved results (the results file is unchanged).

## Not verified (honest gaps)
- **Lander and Waterman (1988)** is cited as the origin of the model but the reference itself was not checked: `TODO-VERIFY` before the report. The formulas used in the code were derived or tested here and do not depend on that citation.
- **Reverse-complement reads and reads crossing the origin** of circular phiX174 are not simulated (D-017).
- Substitutions are uniform over the three other bases; real sequencers have position- and context-dependent errors and homopolymer indels. Stated as a limitation, not modelled.
- The validation checks the simulator against its own specification and against theory; it cannot show that real sequencing data look like this.
- Statistical tests are run at fixed seeds, so they are reproducible, but a pass is evidence, not proof. With about 90 checks at 4 standard errors, a rare chance failure on another seed would not be a bug by itself.

## Deviations from the plan
- None in scope. Extras: the exact contig formula and figure S1, which the plan did not require but which the report needs (rule R4 asks to verify the contig formula).
- The check rule change D-022 (above).
- Phase 2 took well under its 1.0-day window, so the schedule now has slack of roughly half a day.

## Facts that matter for Phase 3
- `simulate_reads(genome, config)` returns a `ReadSet` with `.sequences`, `.starts`, `.spans` and the error counts. The true order of the reads is `argsort(starts)`; this gives the **true ordering** whose objective value `F` is the reference line required in Phase 5.
- With about 100 base reads at c = 10 on phiX174 there are on average 1.4 true contigs at theta = 30, and 8 at c = 5. An assembler can only recover as many pieces as the data allow; Phase 5 must compare with this, not with a single genome.
- Reads from repeated parts of a genome look identical, so the best tour may differ from the true order (R-09).
- Duplicate reads and reads contained in other reads occur by chance (at c = 20 on phiX174, about 1,077 reads on 5,287 start positions); Phase 3 removes them and records the counts.

## Open risks
- **R-02, LKH-3 unreachable: still open.** Retried again at the end of this phase (both hosts time out; `tools/` is empty). The plan's decision point is now. Proposal D-025: continue with Phase 3 and take the final decision at the start of Phase 4 (about Oct 24), when LKH-3 is first needed. The options at that point are the ones in the plan: you place the tarball in `tools/` (and I build it), or LKH-3 is reported as unavailable and no other solver is called LKH-3.
- R-01 (schedule): ahead of plan by about half a day.
- Q3, Q5 unanswered, neither blocks Phase 3.

## Test output
```
204 passed in 30.37s     (make check: ruff, ruff format --check, mypy, pytest all exit 0)
93 checks, 0 failed      (genomeasm simulate validate, 31 s)
```
