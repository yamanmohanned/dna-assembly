# Phase 6: the application

Status: complete, **feature freeze** (plan: end of Phase 6). From now on only fixes, documents and
rehearsal. Waiting for your "continue" (Phase 7: report, slides, oral-exam sheet). Date: 2026-10-24.

## 1. The short version

An offline Streamlit application with the eight pages of the master prompt. Start it with

```
make app                      # or: .venv\Scripts\python.exe -m genomeasm app run
```

and open http://localhost:8501. It reads the stored genomes (`data/raw`) and the precomputed results
(`results/phase5`); it never uses the network (`.streamlit/config.toml` turns off usage statistics).
Every chart has SVG and PNG download buttons; every page says whether what you see is a **live run
on this computer** or a **precomputed result**.

* **Measured, not assumed:** the worst interactive run I could construct (the largest instance the
  interface allows, the solver at its maximum budget, live charts) takes **21.2 s** on this laptop
  against the 30 s rule (section 3). LKH-3 live is included now that WSL answers again
  (18.4 s).
* **Honest content:** the Experiments page opens with the verdict *not supported* for the
  pre-registered claim and shows the post-hoc variants as exploratory (D-045, D-046, D-050).
* **A bug found by looking at the live page** (D-049): the first live run did about 10 ALNS iterations
  in 6 s because a slow chart redraw triggered the next redraw at once. Fixed with a test; after the fix
  the same 6 s run did about twelve thousand iterations (seen on the page).

## 2. The eight pages

| Page | What it shows | Live or precomputed | Main controls |
|---|---|---|---|
| Learn | the story in plain words; a toy example you can edit (overlap graph, matrix, your own order, all orders, greedy decision by decision); the objective and its TSP form; searchable glossary (from `docs/glossary.yaml`) | live | reads, tau, your order, search box |
| Data | genome cards, base counts, GC profile, the data-loading pipeline diagram, a fingerprint check button | stored files | genome, GC window |
| Simulator | coverage depth against Lander-Waterman and the exact model, map of read origins, the reads (FASTA download) | live | genome window, read length, coverage, error rate, seed |
| Overlap Graph | overlap heatmap, overlap-length distribution, degree distribution, counts | live | the above plus tolerance and tau |
| Solver Lab | any solver with a budget; live convergence curve and, for ALNS, the live operator probabilities; outcome chart; quality against the real genome; **step-through** of greedy and ALNS on a tiny instance | live | algorithm, policy, budget, seed, tuned or default parameters |
| Experiments | verdict, all Phase 5 figures with captions and downloads, every statistical test with filters, an interactive solver comparison, scaling, tuning | **precomputed** | family, scope, solver, instance, measure |
| Assembly Viewer | dot plot of the assembled text against the reference, identity, genome fraction, contigs, N50, misjoins, order against true positions, contigs as FASTA | live | order from: true order or any solver |
| Reproducibility | commit, library versions, configuration hash and seed scheme, data checksums, stage-by-stage origin of the precomputed results, result file fingerprints | both | none |

Live limits: at most 1000 reads, a window of at most 5000 letters, a solver search of at most
20 s; bigger runs are the Experiments page.

## 3. The 30-second rule, measured

`genomeasm app benchmark` (`genomeasm/appkit/timing.py`, result `results/phase6/interactive_timing.json`,
commit `4813d193a6`, clean tree) walks through what a page does after a click, with the worst values
of every slider (999 reads, 3% errors, tolerance and tau at their loosest for the overlap graph).
Seconds:

| Path | Instance | Search | Charts and quality | Total |
|---|---|---|---|---|
| Solver Lab, ALNS, 20 s | 0.5 | 20.3 | 0.3 | **21.2** |
| Solver Lab, tabu + candidate lists, 20 s | 0.5 | 20.0 | 0.3 | **20.8** |
| Solver Lab, LKH-3, 15 s (WSL) | 0.4 | 17.7 | 0.3 | **18.4** |
| Assembly Viewer, ALNS, 15 s | 0.0 | 15.0 | 0.9 | **15.9** |
| Overlap Graph, loosest settings, all charts | 0.7 | 0.0 | 1.8 | **2.6** |
| Simulator, largest set, both charts | 0.0 | - | 1.0 | **1.1** |

Notes. (1) The browser needs a little more time to draw than this benchmark measures; the margin is
about 9 s. (2) The overlap graph is built with all processor threads here (0.5 s for 999 reads); in
the Phase 5 experiments each of five workers used one thread (2.6 s for 994 reads), so build times differ, results
do not. (3) A live redraw costs 0.5 s on this laptop; the page spaces redraws so that
they use about a tenth of the run. A live run therefore searches for roughly 90% of its budget and its
numbers are not comparable with Phase 5 (the page says so). (4) Starting the server compiles the
numba kernels once (0.7 s of warm-up when the on-disk cache exists; a cold start without the cache is slower and was not
measured, `TODO-VERIFY`).

## 4. How it was checked

* `tests/test_appkit.py` (data layer: glossary, toy example, greedy steps against the solver, limits,
  tuned parameters, every solver, LKH-3 never substituted, dot plot on the true order, every chart as SVG
  and PNG, results and reproducibility information) and `tests/test_app.py` (every page runs headless with
  Streamlit's `AppTest`; interactions: edit the toy order, search the glossary, change genome, over-large
  instance refused, run a solver, step-through, assemble the true order, verdict page, navigation).
  `tests/test_phase6_results.py` guards the timing numbers above. `tests/test_solvers.py` gained the checks of
  the observation hooks. Whole suite: 460 tests.
* **By eye in the browser** (built-in browser, 1440 x 900): Learn (toy example), Data, Simulator, Overlap
  Graph, Solver Lab (a real 6 s run with live charts), Experiments (verdict), Assembly Viewer (greedy
  assembly), Reproducibility. Not looked at by eye: the other tabs of Experiments (their data is covered by
  the tests), a projector, other screen sizes, another laptop (open: Q5).

## 5. What was found while building (all handled)

| Finding | Handling |
|---|---|
| Live redraws starved the search (D-049) | throttle starts after the callback; unit test with a fake clock |
| Live redraws still cost about half the time at first (0.5 s each) | one combined frame, shown as a PNG, redraw spacing adapts to its cost |
| A click on a download button reruns a Streamlit page, which would rerun the solver | `on_click="ignore"` on every download; results kept in `st.session_state`; tested |
| `st.pyplot` with extra arguments is deprecated | live frames are PNG bytes shown with `st.image` |
| Default instance (1% errors) showed a flat convergence line: greedy is already optimal | default 3,000 letters, 10x, 2% errors (D-050) |
| Matplotlib warning: constrained layout failed with long read texts in the toy graph | texts cut to 12 letters |
| `app/` is not type-checked (`mypy` on Streamlit calls is noisy) | covered by the headless page tests instead |

## 6. Limits and open points

* Not measured: another computer, a projector, the browser's own drawing time; the demo script with a
  fallback of screenshots belongs to Phase 8.
* The application shows two post-hoc algorithms (`tabu_cl`, `ucb1_rt`) next to the pre-registered ones;
  they are labelled in the lists, but whoever presents must say it out loud (slides in Phase 7).
* Page texts are written for a 13-year-old first and precise second; they have not been read by anyone
  else yet.
* Questions still open: Q3 (program that opens the .pptx), Q5 (is this the presenting laptop).
* Decisions: D-047 to D-050 (50 decisions in total, D-001 to D-050).
