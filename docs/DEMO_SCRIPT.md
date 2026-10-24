# Demo script: timed, with a fallback for every step

Total 15 minutes live (a 5-minute version is at the end). Everything you can show live also exists as a
stored picture or figure, so a failed demo costs a few seconds, not the talk.

## A. Before the talk

| When | Do | Check |
|---|---|---|
| 30 min before | Plug in the laptop; close other programs; set the browser zoom to 110% and the screen to the projector's resolution | the battery is charging; no pop-ups |
| 15 min before | `make app` (or `.venv\Scripts\python.exe -m genomeasm app run --browser`) | the Guide page appears; the sidebar lists nine pages |
| 10 min before | Open **Solver Lab**, press **Run** once and wait 6 s | a curve is drawn and a result appears (this also warms up the code) |
| 10 min before | Optional, only if you want to show LKH-3 live: `wsl -d Ubuntu-24.04 -- echo ok` | prints `ok`; if not, skip LKH-3 (it stays "unavailable", never replaced) |
| 5 min before | Open the backup files in tabs: `reports/slides/deck.pdf`, `reports/report/report.pdf`, the folder `reports/screenshots/` | all open |
| 2 min before | Back to the Guide page; press F11 (full screen) | readable from the back row |

## B. The 15-minute run

| Min | Page | Click-by-click | What you should see | Say (one idea) | If it fails |
|---|---|---|---|---|---|
| 0:00-2:00 | Guide, then Learn tab "DNA and RNA" | Show the helix banner; open Learn; move the "Start" number; switch to "RNA single strand" | DNA with A-T and G-C pairs; RNA as one strand with U | "A genome is a long word over four letters; we get it as shredded pieces." | `reports/screenshots/01_guide.png`, `02_learn.png` |
| 2:00-4:00 | Learn, tab "Toy example" | Keep reads `ATGCG GCGTA GTACC`; order `ABC`, then `CBA` | F falls from 6 to 1; merged text 9 letters against 14 | "Each arrow says how many letters we save; the best order saves most." | slide 5 of the deck; Guide section 4 |
| 4:00-5:00 | Data | Choose the lambda genome; press "Check the fingerprints" | cards, GC curve, "exactly as downloaded" | "Real NCBI genomes, stored once and fingerprinted; the program runs offline." | `03_data.png` |
| 5:00-7:00 | Simulator, then Overlap Graph | Raise the error rate to 2%; open the Overlap Graph | red reads (errors); a bright diagonal band | "Reads carry errors; the band is the overlaps of true neighbours." | `04_simulator.png`, `05_overlap_graph.png` |
| 7:00-10:00 | Solver Lab | Keep ALNS and the defaults; press **Run** (6 s); then open "Step through a tiny instance", choose the tiny simulated instance, ALNS, move the slider | the curve rises above the greedy line; the operator probabilities move; one decision per step | "Same time for everyone; this panel shows the bandit choosing repair tools." | `06_solver_lab.png`; figures `fig04_convergence`, `fig08_arm_probability` |
| 10:00-13:00 | Experiments, tab "Verdict", then "Compare solvers" | Read the verdict; select Greedy, ALNS random, ALNS UCB1, LKH-3; the measure "gap" | the claim is **not supported**; LKH-3 at zero; ALNS close | "300 read sets, fair budgets, tuning on separate data: the bandit did not beat random choice." | `07_experiments.png`; slides 14 and 15 |
| 13:00-14:00 | Assembly Viewer | Order "The true order", **Assemble**; then "Greedy edge merging" | straight dot-plot line, no misjoin | "Fitness is not the genome, so we also grade against the real one." | `08_assembly_viewer.png`; slide 16 |
| 14:00-15:00 | Reproducibility | one scroll | commit, seeds, library versions, checksums | "Everything is repeatable; every number comes from stored results." | `09_reproducibility.png` |

**Say out loud (do not skip):** the result is negative and was fixed in advance; the corrected UCB1 and the
tabu variant were added afterwards and are exploratory; "best known" is not a proven optimum.
**Do not:** say "optimal" for ALNS; move sliders that are not in the table; start a second run while one is running.

## C. When something goes wrong

| Problem | Do |
|---|---|
| The application does not start | Say "the pages are stored"; show `reports/screenshots/` and the PDF deck; continue with slides 14 to 16 |
| A page shows a red error box | Press F5; if it persists, go to the next page and use its stored picture |
| Solver Lab seems stuck | A run never lasts more than 20 s; wait, then press F5 and use `06_solver_lab.png` plus `fig04_convergence.png` |
| LKH-3 "unavailable" | Correct and honest: the application never replaces it; show the stored results (Experiments) |
| Projector resolution is low | zoom 125% (Ctrl and +); the Experiments page also works in the PDF deck |
| The pptx does not open | use `reports/slides/deck.pdf` (same slides, same order) |

## D. Five-minute version

Guide helix (0:30) → Learn toy example (1:30) → Solver Lab run (1:30) → Experiments verdict (1:30). Say the three honest
sentences from section B.

## E. Questions

`docs/ORAL_EXAM.md` has 20 questions with model answers; `docs/GUIDE.md` section 9 has the short list; the glossary is on the Guide page.
