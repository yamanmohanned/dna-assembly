# Guide: from zero knowledge to presenting this project

This guide assumes you know nothing about genomes or optimisation. Read sections 1 to 6 once
(about 25 minutes), practise section 8 twice, keep section 9 next to you. Every number below comes
from the project's own results; the pages of the application show them live.

## 1. The big picture

**The problem in one sentence.** A machine that reads DNA cannot read a whole genome at once; it
reads millions of short, overlapping, slightly wrong pieces (*reads*), and a program must put them
back in the right order (*assembly*).

**The picture.** Take a book, shred many copies, and keep only random scraps of a few lines each.
Two scraps that end and begin with the same words probably sat next to each other. Rebuilding the
book means ordering the scraps so that neighbours overlap as much as possible.

**The mathematical idea of this project.** "Order the scraps so that the overlaps add up to as much
as possible" is exactly the **travelling salesman problem (TSP)**: visit every scrap once and travel
as little as possible, where travelling is cheap between scraps that overlap a lot. Because the
overlap from A to B differs from B to A, it is the *asymmetric* TSP (ATSP). The proof is in
`docs/OBJECTIVE_FUNCTION.md`.

**The research question.** Good TSP methods exist. One of them, *ALNS*, repeatedly breaks the
current order and mends it with one of six repair tools. Which tool is best changes during the search.
A *multi-armed bandit* is a learning rule that picks tools adaptively. **Claim to test: does a bandit
pick tools better than a coin toss, under the same time budget?** The answer, found by 300 experiments
and honest statistics, is in section 6.

## 2. DNA and RNA in five minutes

* **DNA** is the instruction book of a living thing. It is a very long word written with four
  letters, called *bases*: **A** (adenine), **C** (cytosine), **G** (guanine), **T** (thymine).
* In cells DNA is a **double helix**: two strands wound around each other. The strands are glued
  by *base pairs*: **A always pairs with T, and G always pairs with C**. So one strand determines the
  other. (Open the application, page Learn, tab "DNA and RNA", to see a real stretch drawn.)
* **RNA** is a close relative. It is usually a **single strand**, and it uses **U** (uracil) in
  place of T. Cells copy a stretch of DNA into RNA to carry instructions. Some viruses keep their
  whole genome as RNA instead of DNA.
* A **genome** is the complete DNA (or RNA) of an organism. Its **length** is the number of letters.
  A human genome has about three billion letters; the two virus genomes here have thousands.
* **Sequencing** reads letters. Machines return **reads**: pieces of 50 to 300 letters, taken at
  random places, with a few **errors** (a wrong letter). **Coverage** says how many reads cover an
  average letter: 10x means every letter is read about ten times.
* This project handles DNA letters as one string over A, C, G, T. RNA would be the same string with U
  instead of T, so the method carries over; we did not test RNA data.

## 3. The genomes used

Both are real sequences from NCBI, the public sequence library (RefSeq). They were downloaded once,
fingerprinted with SHA-256, and are stored in `data/raw/`; the application never needs the internet.

| | phiX174 | lambda |
|---|---|---|
| Accession | NC_001422.1 | NC_001416.1 |
| What it is | a virus that infects the bacterium *E. coli* | another *E. coli* virus |
| Length | 5,386 letters | 48,502 letters |
| Kind | circular, single-stranded DNA | linear, double-stranded DNA |
| GC content | 44.8% | 49.9% |
| Role in the project | tuning (positions 0 to 2,500) and some evaluation windows | evaluation windows (positions 10,000 and up) |

*Why viruses?* Their genomes are small, complete and public, so the answer key is known and every
experiment can be repeated on a laptop. *(phiX174 was the first DNA genome ever sequenced, by Sanger and
colleagues in 1977; citation `TODO-VERIFY`.)* We treat every genome as one linear string (decision D-017).

## 4. The examples used

1. **The toy example** (page Learn). Three reads: A = `ATGCG`, B = `GCGTA`, C = `GTACC`. The end of A
   (`GCG`) equals the start of B, so gluing B after A saves 3 letters; B and C overlap by 3 (`GTA`).
   Order ABC saves 6 letters: the merged text `ATGCGTACC` has 9 letters instead of 15. Order CBA saves only
   1. The *fitness* F is the number of letters saved; larger is better. You can type your own reads.
2. **Simulated reads** (pages Simulator, Overlap Graph, Solver Lab, Assembly Viewer). A window of a
   real genome is cut into random reads of a chosen length with a chosen error rate. The program knows
   where each read really came from, so it can grade any answer. Settings: window, read length,
   coverage, error rate, seed (the seed makes a run repeatable).
3. **The 300 experiments** (page Experiments). Ten settings times 30 seeds. A name such as
   `L5-c10-e1.0` reads: **L**ambda, **5**,000-letter window, **c**overage 10x, **e**rror 1.0%; `P` is a
   phiX174 window; `l60` means 60-letter reads. Every solver got the same 8 seconds on each of the
   300 read sets.
4. **The tiny instance** (Solver Lab, step-through). Seven or eight reads, small enough to watch
   every decision of greedy and of ALNS.

## 5. The techniques

**Overlap and its errors.** For two reads the *overlap* is the longest end of the first that equals
the start of the second (at least a minimum length tau, so that chance matches are ignored). With
sequencing errors we allow a few mismatches (for example 5%). The table of all overlaps is the
*overlap graph*; page Overlap Graph shows it as a heatmap.

**The objective.** Fitness F of an order = sum of the overlaps of neighbouring reads (maximise).
The merged text has length L = (total letters of all reads) - F. Both statements are the same goal.
Different orders: n reads have n! orders (20 reads: more than 2 quintillion), so we search cleverly.

**Greedy.** Sort all overlaps from largest to smallest; accept an overlap if the first read has no
successor yet, the second no predecessor yet, and it makes no loop. Fast, and surprisingly good.

**Tabu search.** Start from greedy; repeatedly make the best small change (swap two reads, or move one
read) even if it is worse, but forbid undoing recent changes for a while ("tabu list") so the search does
not walk in circles. We also tried a version with *candidate lists* (try only the most promising places).

**LKH-3.** A famous, very strong solver for the travelling salesman problem (Helsgaun's
Lin-Kernighan-Helsgaun heuristic, third version). We call it as an outside program and verify every tour
it returns. It is our reference: it was never beaten in the experiments. If it cannot run, the
application says so; nothing is silently put in its place.

**ALNS (adaptive large neighbourhood search).** Each step *destroys* part of the current order (removes
some reads: at random, the worst-connected ones, or a stretch) and *repairs* it (puts them back at the
best place, either greedily or by "regret": the read with most to lose goes first). That gives 3 x 2 = 6
tool pairs ("arms"). A worse result is sometimes accepted (*simulated annealing*), early in the run more
than late, to escape dead ends.

**Bandits (choosing the tool).** The policy picks one of the 6 arms each step and gets a *reward* (big
for a new best order, small for an improvement, none otherwise). Policies compared: **random** (coin toss,
the baseline), **round robin** (turns), **roulette wheel** (the classical adaptive rule: weights follow
recent rewards), **epsilon-greedy** (usually the best arm, sometimes a random one) and **UCB1** (best
average reward plus a bonus for arms tried rarely). The tension between *exploring* new arms and
*exploiting* the best one is the heart of a bandit.

**Fair comparison.** Same wall-clock time for everyone; settings of every method tuned on *separate*
read sets (phiX174 positions 0 to 2,500), then tested on lambda and other phiX174 positions; 30 seeds per
setting; the same read set for every method (paired); statistics: *Wilcoxon* test (is the difference
more than luck?), *Holm* correction (many tests need a stricter threshold), *A12* (chance that the first method
wins on a read set; 0.5 means no difference).

**Quality against the real genome.** Because F is not the genome, every result is also graded against the
answer key: *identity* (correct letters), *genome fraction* (share rebuilt), *contigs* (gap-free
pieces; fewer is better), *N50* (typical piece length), *misjoins* (a piece that jumps to a wrong place).

## 6. What we found, in plain words

1. **The claim is not supported.** UCB1 and epsilon-greedy did **not** beat random tool choice (UCB1: no
   difference; epsilon-greedy: slightly worse). This is a result, not a failure: the project reports it as it is.
2. **We found a bug in our own UCB1** after seeing the results: when nothing improved any more it always
   picked the first tool. A corrected version (ties broken at random) is slightly better than random and as
   good as the classical roulette wheel. It was added afterwards, so it is labelled **post hoc / exploratory**.
3. **LKH-3 was never beaten** in all 300 read sets. ALNS is clearly better than greedy and tabu search.
4. **All methods rebuild the same genome** (no misjoins; about 99.6% of the window); they differ in how
   many letters they save and in how many pieces remain.
5. **Limits:** only substitution errors and one strand; two viruses; one laptop; "best known" means the best
   solution *we* found, not a proven optimum.

## 7. Using the application

Start: `make app`, or `.venv\Scripts\python.exe -m genomeasm app run --browser`; address
http://localhost:8501. The first start compiles code for a second or two. Stop with Ctrl+C.

| Page | What to do | What you should see |
|---|---|---|
| Guide | read, or use the presentation script | this document |
| Learn | tab "DNA and RNA": move the slider; tab "Toy example": type the order `BCA` | a helix of real DNA; F falls from 6 to 3 |
| Data | choose the other genome; press "Check the fingerprints" | cards, GC curve, pipeline, "all as downloaded" |
| Simulator | raise the error rate; change the seed | read bars turn red; coverage matches theory |
| Overlap Graph | tick a higher tolerance | a bright diagonal band = true neighbours |
| Solver Lab | keep the defaults, press **Run** | curve climbs above greedy for 6 s |
| Solver Lab, step-through | tiny instance, ALNS, move the slider | one decision per step |
| Experiments | tab "Verdict", then "Compare solvers" | the honest verdict and boxplots |
| Assembly Viewer | order from "The true order", press Assemble | straight dot-plot line, 0 misjoins |
| Reproducibility | read | commit, seeds, checksums |

Tips: every chart has SVG and PNG buttons; sliders in the left bar are shared between pages; a green badge
means "computed now", an orange badge means "precomputed". If a page seems stuck, wait for the busy
indicator, then press F5. Restart the server after changing code.

## 8. Presenting to a professor (15 minutes)

**Before:** start the app 10 minutes early; open Solver Lab and press Run once (warms the code); close other
programs; plug in the power; browser zoom 110% and full screen (F11); check the projector shows the green
and orange badges legibly; have `docs/phases/PHASE_5.md` and this guide open as backup.

| Min | Page | Say | Do |
|---|---|---|---|
| 0-2 | Guide / Learn | "Genomes are read as shredded pieces; we must reorder them. I show that this is a travelling salesman problem and test whether a learning rule helps." | show the helix; the shredded-book picture |
| 2-4 | Learn, toy example | "Three pieces; each arrow says how many letters we save." | type `ABC`, then `CBA` |
| 4-5 | Data | "Real NCBI genomes, stored and fingerprinted; the program works offline." | pipeline diagram |
| 5-7 | Simulator, Overlap Graph | "Reads with errors; the bright band is the overlaps of true neighbours." | raise errors |
| 7-10 | Solver Lab | "Same time budget; watch the order improve; this panel is the bandit choosing tools." | Run (6 s), then step-through |
| 10-13 | Experiments, Verdict | "300 read sets, fair budgets, tuned on separate data. The claim is **not supported**; LKH-3 is never beaten." | Verdict, then boxplots |
| 13-14 | Assembly Viewer | "Fitness is not the genome, so we also grade against the real one." | true order, then greedy |
| 14-15 | Reproducibility | "Everything is repeatable: seeds, hashes, commit." | one scroll |

**Say out loud:** the main result is negative and was fixed in advance; the corrected UCB1 and the tabu
variant were added afterwards and are exploratory; "best known" is not a proven optimum.
**Do not:** claim the bandit wins; say "optimal" for ALNS; improvise sliders beyond the table during the talk.
**If the live demo fails:** open the Experiments page (it needs no computation), show `reports/figures/fig05`
and `fig07`, and tell the story from the toy example.

## 9. Questions to expect

* *Why is this a travelling salesman problem?* Visiting every read once while saving the most overlap equals a
  shortest tour with cost = read length minus overlap (proof in `docs/OBJECTIVE_FUNCTION.md`).
* *Why can't you try all orders?* n! grows faster than any computer; 20 reads already take years.
* *Why not just use LKH-3?* It is our reference and wins; the research question is about the learning rule
  inside ALNS, not about beating LKH-3.
* *Why did the bandit not win?* Six similar tools and short runs leave little to learn (a hypothesis, not tested);
  our first UCB1 also had a tie-breaking bug.
* *Isn't fixing the bug after seeing results cheating?* We kept the original, reported both, and labelled the
  new one exploratory.
* *Is 30 seeds enough?* It is the required minimum; effects are small (A12 about 0.5), so we report effect
  sizes, not only p-values.
* *How do you know the assembly is right?* Reads come with known positions; we grade identity, genome
  fraction, contigs, misjoins.
* *What about real data?* Real sequencers add insertions, deletions, both strands and repeats; not covered.
* *What is a hyper-parameter, and who chose them?* A setting of a method (for example the tabu tenure); chosen by
  random search on separate data (24 trials), the same effort for every method.
* *Does it work for RNA?* The method only sees letters; RNA with U in place of T is the same problem. Not tested.

The glossary (tab "Glossary" of the Guide page, or `docs/GLOSSARY.md`) explains every term.
