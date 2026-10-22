# PHASE 3: Overlap graph and ATSP formulation

Status: **complete, waiting for your "continue"** to start Phase 4 (solvers). Date: 2026-10-23 (local time; timestamps inside result files are UTC).
Plan window: 1.0 day (Oct 23 to Oct 24). Exit gate from the plan: the four mandatory tests pass, and the walking skeleton runs end to end on real phiX174 reads with a trivial ordering. **Both met.**

## What was done
1. **Overlaps** ([genomeasm/overlap.py](../../genomeasm/overlap.py)):
   - `make_substring_free`: drops duplicates and contained reads and counts them (D-026).
   - `overlap_matrix`: exact overlaps with a prefix index (about `n * length` look-ups instead of `n^2`).
   - `overlap_matrix_naive`: the definition, written as simply as possible, as the test reference.
   - `overlap_matrix_tolerant`: error-tolerant overlaps, exhaustive search compiled with numba (D-028).
   - `recommend_min_overlap` and `expected_spurious_pairs`: the choice of tau (D-027).
   - `truth_agreement`: compares an overlap matrix with the true positions of simulated reads.
2. **Objective and ATSP** ([genomeasm/objective.py](../../genomeasm/objective.py)): `F`, `L`, the cost matrix with the dummy node, tour/order conversion, cycle cost, TSPLIB writer and reader (D-029).
3. **Layout and consensus** ([genomeasm/assembly.py](../../genomeasm/assembly.py)): chain layout, contigs split at zero overlap, majority-vote consensus (D-030).
4. **Instance** ([genomeasm/instance.py](../../genomeasm/instance.py)): one object joining all of the above, with the true ordering when the reads are simulated.
5. **Validation** ([genomeasm/overlap_validation.py](../../genomeasm/overlap_validation.py), [configs/phase3_overlap.yaml](../../configs/phase3_overlap.yaml)): 95 checks; results in [results/phase3/overlap_validation.json](../../results/phase3/overlap_validation.json), run at clean commit `460ea19`.
6. **Figures 2 and 3** of the required list: overlap-length distribution and overlap heatmap (SVG, PDF, 300-dpi PNG, with caption files).
7. **Walking skeleton**: `genomeasm assemble ACC --coverage C --order true|trivial|random` (simulate, preprocess, overlaps, order, layout, consensus, FASTA + summary with provenance). Two runs are committed in `results/phase3/`.
8. **Docs**: [docs/OBJECTIVE_FUNCTION.md](../OBJECTIVE_FUNCTION.md) (the objective, the proof of the ATSP identity, the worked toy example; for the oral exam), 11 glossary terms, decisions D-026 to D-032.
9. Small refactor: the check helpers shared by the validation suites moved to `genomeasm/checks.py` (Phase 2 tests re-run green).

## The four mandatory tests of the master prompt

| | Requirement | Test | Result |
|---|---|---|---|
| (a) | Toy reads `ATGCG`, `GCGTA`, `GTACC`: matrix `[[-,3,1],[1,-,3],[0,0,-]]`, best order A to B to C, merged `ATGCGTACC` | `test_toy_overlap_matrix_matches_the_master_prompt`, `test_toy_best_order_is_a_b_c_with_merged_string` (all 6 orders checked; fast, naive and tolerant code agree) | pass |
| (b) | `L(pi) = sum(len) - F(pi)` equals the length of the actually merged string, on random substring-free instances | `test_merged_length_equals_the_length_of_the_merged_string` (120 random instances x 4 orders, hypothesis) and `test_every_read_sits_in_the_merged_string_at_its_layout_offset` | pass |
| (c) | Cycle cost with the dummy node equals `L(pi)` | `test_cycle_cost_with_the_dummy_node_equals_merged_length`, `test_the_identity_holds_for_any_overlap_matrix_and_lengths`, brute force over all orders for small n | pass |
| (d) | Overlap code equals a naive reference on random inputs | `test_fast_overlap_equals_naive_on_random_reads` (150 random cases), `..._cut_from_one_genome` (60), and on the real data in the validation (68,382 ordered pairs of phiX174 reads compared) | pass |

## Evidence (all from `results/phase3/overlap_validation.json`, 5 runs per setting, read length 100)

### Every true overlap is found
Error-free reads, all 42 settings (2 genomes x 3 coverages x 7 values of tau): in total **7,060,551 true overlapping pairs checked, 7,060,551 found** (recall exactly 1 everywhere, as it must be: the true overlap is always a valid match). The code found the true overlap or a longer chance match for each of them.

### The recommended minimum overlap

| Genome | c | reads after preprocessing (duplicates removed) | tau | true overlapping pairs | edges | spurious edges (+-SE) | random-model bound | share of edges that are exactly the true overlap | F of true order | contigs of true order |
|---|---|---|---|---|---|---|---|---|---|---|
| phiX174 | 5 | 263 (6) | 9 | 1170 | 1171 | 1.0 (+-0.3) | 0.35 | 0.9992 | 21017 | 4.0 |
| phiX174 | 10 | 514 (25) | 10 | 4509 | 4509 | 0.8 (+-0.4) | 0.34 | 0.9998 | 46052 | 1.0 |
| phiX174 | 20 | 975 (102) | 11 | 15861 | 15862 | 0.4 (+-0.2) | 0.30 | 1.0000 | 92166 | 1.0 |
| lambda | 5 | 2364 (61) | 12 | 10104 | 10105 | 1.2 (+-0.4) | 0.44 | 0.9999 | 188260 | 32.2 |
| lambda | 10 | 4604 (246) | 13 | 38232 | 38233 | 1.6 (+-0.5) | 0.42 | 1.0000 | 411940 | 1.8 |
| lambda | 20 | 8794 (906) | 14 | 137207 | 137208 | 0.2 (+-0.2) | 0.38 | 1.0000 | 830866 | 1.0 |

(A "spurious" edge is an edge between two reads whose true positions do not overlap by at least tau. No contained reads were found to remove, because all reads have the same length; only exact duplicates occur.)

### Why tau matters (lambda, c = 10)

| tau | edges | spurious (+-SE) | random-model bound | true overlapping pairs | contigs of true order |
|---|---|---|---|---|---|
| 8 | 41071 | 683.4 (+-8.8) | 431.20 | 40387 | 1.4 |
| 10 | 39574 | 56.6 (+-1.1) | 26.95 | 39517 | 1.6 |
| 12 | 38675 | 5.2 (+-1.5) | 1.68 | 38670 | 1.6 |
| 16 | 36898 | 0.0 (+-0.0) | 0.01 | 36898 | 2.6 |
| 20 | 35155 | 0.0 (+-0.0) | 0.00 | 35155 | 2.8 |
| 30 | 30748 | 0.0 (+-0.0) | 0.00 | 30748 | 5.2 |

Too small a tau lets in hundreds of chance edges; too large a tau loses true overlaps and splits the genome into more contigs (1.6 contigs at tau = 12, 5.2 at tau = 30). The random-sequence bound is a **guide, not a guarantee**: on lambda at c = 10 the real count is 1.6 to 3.1 times the bound for tau = 8 to 12, and at the recommended tau it is between 0.5 and 3.8 times the bound over the six settings (still at most 1.6 chance edges per instance). (D-027)

### Objective identities and layout on real reads
- For the true order of every instance: `L = sum(len) - F` equals the length of the plain merge and equals the total length of the contigs (42 of 42 settings).
- The cycle cost with the dummy node equals `L` for the true order and a random order (35 of 42 settings; in the other 7 the dense cost matrix would not fit in memory, so only the layout and merge identities were checked).
- Gluing the error-free reads in their true order gives contigs that are **exact substrings of the real genome** in all 6 settings at the recommended tau (30 of 30 runs).

### Walking skeleton on real phiX174 (c = 10, seed 1, the same 516 reads and overlap graph)

| Order | F | L | contigs | contigs exactly inside the genome |
|---|---|---|---|---|
| true order (reads sorted by real position) | 46,235 | 5,365 | 1 | 1 of 1 |
| trivial order (as the sequencer produced them) | 435 | 51,165 | 508 | 508 of 508 |

The true order gives one contig of 5,365 letters out of 5,386 in the genome; a trivial order is almost useless, which is what a solver has to improve on.

### Why the error-tolerant overlap exists (phiX174, c = 10, 1% substitution errors, tau = 20, 5 runs)

| Overlap type | Share of true overlaps recovered | Edges | Spurious edges |
|---|---|---|---|
| exact | 0.316 (+-0.009) | 1357 | 0.0 |
| tolerant (5% of k) | 0.931 (+-0.003) | 3991 | 0.0 |

Pre-set pass rules (`configs/phase3_overlap.yaml`): exact at most 0.60 and tolerant at least 0.85; both met with margin.

### Figures
- [fig02_overlap_length_distribution.png](../../reports/figures/fig02_overlap_length_distribution.png): overlap lengths, true neighbours versus other pairs, for both genomes. The true overlaps are spread evenly over all lengths; the other pairs sit only just above tau.
- [fig03_overlap_heatmap.png](../../reports/figures/fig03_overlap_heatmap.png): overlap matrix of 69 phiX174 reads sorted by true position: a band above the diagonal only (overlaps are directional), with a visible gap where the chain breaks.

### Speed
Building the whole overlap graph (preprocessing plus prefix-index overlaps) takes 0.01 s for 263 reads, 0.16 s for 4,604 reads and 0.44 s for 8,794 reads on this laptop (mean of 5 runs). The tolerant search is quadratic: it handled the noise experiments (about 500 reads) in seconds.

### Tests and checks
`make check` exit code 0: ruff, ruff format, mypy and **pytest 282 passed**. That is 204 from earlier phases, 66 new tests of the Phase 3 code (committed with the code: 270) and 12 tests of the committed results (`tests/test_phase3_results.py`), which need the results and therefore join the suite in the second commit.

## Problems found and how they were handled (R10: investigate surprises)
1. **A test expectation of mine was wrong, and I diagnosed it instead of loosening it.** I had assumed that majority voting on noisy reads (1% errors) would remove at least 90% of the errors. It removed only about 30%. The diagnosis showed that voting itself was perfect where at least 5 reads vote (0 wrong columns in about 5,900), but with a 5% tolerance the true-order chain breaks into 61 to 71 contigs, so half of all columns have fewer than 5 votes. The test now asserts what is true and justified (at most 1 wrong column among well-covered columns, expected 0.006; overall never worse than the raw reads) and the finding is logged as D-031 because it changes how Phase 5 must measure quality.
2. **The random-sequence bound on chance overlaps is optimistic on real genomes** (1.6 to 3.1 times on lambda at tau = 8 to 12). Reported as such (D-027); the recommended tau still gives at most 1.6 chance edges per instance.
3. **Figure 2**: the legend first hid the bars; headroom was added and the figure redrawn.
4. **mypy and numpy**: `a.min()` on typed arrays triggers an overload error in the current stubs; `np.min(a)` is used (noted in CLAUDE.md).

## Not verified (honest gaps)
- **LKH-3 reading the TSPLIB file**: the file follows the TSPLIB format and is read back correctly by my own reader, but no LKH-3 exists on this machine, so the diagonal convention and parameter file are untested (`TODO-VERIFY`, D-029).
- **NP-hardness of the ATSP and the Held-Karp reference** are cited from general knowledge only (`TODO-VERIFY` before the report).
- The error-tolerant overlap handles substitutions only. Indel errors are not recovered by it (validated only with substitution errors).
- Overlap computations are tested for reads of one fixed length on real genomes. Variable-length reads are covered by unit tests (including contained reads) but not by the real-data validation.
- Memory: the dense overlap matrix of the largest instance (8,794 reads) takes about 0.3 GB; instances much larger than 10,000 reads would need a sparse representation. Not needed for the planned experiments (n up to about 1,000 to 2,000).

## Deviations from the plan
- None in scope. Extras: the truth comparison, the tolerance experiment, Figures 2 and 3 (planned for the report anyway), and the objective-function document, which the examiner focus calls for.
- Delta evaluation was moved to Phase 4 (D-032); the Phase 3 acceptance criteria do not include it.
- The plan's "walking skeleton" is a command line run; the Streamlit app stays in Phase 6.

## Facts that matter for Phase 4
- `AssemblyInstance.from_readset(reads)` gives `n`, `ov`, `lengths`, `cost()` (an `(n+1) x (n+1)` integer matrix with the dummy node 0), `fitness(order)`, `merged_length(order)`, `true_order()` and `write_tsplib(path, name)`.
- Orders are permutations of `0..n-1`; ATSP tours are `[0, pi_1 + 1, ..., pi_n + 1]`.
- Reference lines for every plot: `F` of the true order (e.g. 46,235 for phiX174 at c = 10), and the trivial order.
- Instance sizes in the plan (n up to about 1,000 to 2,000) leave Held-Karp (n <= 15) and brute force (n <= 9) for tiny sub-instances only. Tiny instances can be made from genome windows.
- Decision D-025 stands: the LKH-3 call is made at the start of Phase 4. The TSPLIB writer is ready.

## Open risks
- **R-02, LKH-3: still open** (not retried again today beyond the end of Phase 2; `tools/` is empty). Needed now at the start of Phase 4.
- R-01 (schedule): Phase 3 is done about 1.5 days before the end of its planned window (planned end: Oct 24, 80% of the day); that is the slack carried into Phase 4.
- Q3 and Q5 (slide software, presenting laptop) unanswered; neither blocks Phase 4.

## Test output
```
282 passed in 38.47s     (make check: ruff, ruff format --check, mypy, pytest all exit 0)
                         (126 from Phases 0 and 1, 78 from Phase 2, 66 new in Phase 3 code, 12 new on the committed Phase 3 results)
95 checks, 0 failed      (genomeasm overlap validate, 168 s)
```
