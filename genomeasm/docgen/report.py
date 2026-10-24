"""The PDF report (IMRaD). Text is Markdown; every number comes from ``Facts`` (results/)."""

from __future__ import annotations

from pathlib import Path

from genomeasm.appkit.paths import FIGURES, ROOT
from genomeasm.docgen import md
from genomeasm.docgen.facts import Facts, name, p_text

CSS = """
@page { size: A4; margin: 20mm 18mm 22mm 18mm; @bottom-center { content: counter(page); font-size: 9pt; color: #555; } }
body { font-family: Cambria, 'Times New Roman', serif; font-size: 10.5pt; line-height: 1.42; color: #1d2733; }
h1 { font-size: 22pt; margin: 0 0 4pt 0; } h2 { font-size: 15pt; margin: 20pt 0 6pt 0; border-bottom: 1px solid #c9d2db; padding-bottom: 2pt; }
h3 { font-size: 12pt; margin: 14pt 0 4pt 0; } h2, h3 { break-after: avoid; font-family: Arial, sans-serif; }
p { margin: 5pt 0; text-align: left; } code { font-family: Consolas, monospace; font-size: 9pt; background: #f2f5f8; padding: 0 2pt; }
pre { background: #f2f5f8; padding: 6pt; font-size: 8.5pt; overflow-wrap: anywhere; white-space: pre-wrap; }
table { border-collapse: collapse; width: 100%; margin: 8pt 0; font-family: Arial, sans-serif; font-size: 8.5pt; break-inside: avoid; }
th, td { border: 1px solid #c9d2db; padding: 3pt 5pt; text-align: left; vertical-align: top; } th { background: #e8eef4; }
figure { margin: 10pt 0; break-inside: avoid; text-align: center; } figcaption { font-family: Arial, sans-serif; font-size: 8.5pt; color: #444; text-align: left; margin-top: 3pt; }
img.im { vertical-align: -0.3em; } .eq { text-align: center; margin: 6pt 0; } .cover { margin-top: 70mm; }
.box { border-left: 4px solid #0072B2; background: #f2f7fb; padding: 4pt 10pt; margin: 8pt 0; } .small { font-size: 8.5pt; }
.toc li { margin: 2pt 0; } .cover p { text-align: left; }
"""


def table(header: list[str], rows: list[list[str]]) -> str:
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return chr(10).join(lines)


def eq(expression: str) -> str:
    return "$$" + expression + "$$"


def fig(stem: str, caption: str, width: str = "100%") -> str:
    return md.figure_html(FIGURES / f"{stem}.png", caption, width)


def test_row(F: Facts, a: str, b: str) -> list[str]:
    t = F.test(a, b)
    return [
        f"{name(a)} vs {name(b)}", f"{t.mean_diff_pp:+.3f}", f"{t.median_diff_pp:+.3f}",
        p_text(t.p_holm).replace("p = ", "").replace("p ", ""), f"{t.a12:.3f}", f"{int(t.wins)} / {int(t.ties)} / {int(t.losses)}",
    ]  # fmt: skip


TEST_HEADER = [
    "Comparison",
    "Mean diff. (pp)",
    "Median diff. (pp)",
    "p (Holm)",
    "A12",
    "wins / ties / losses",
]

REFERENCES = [
    "Sanger F. et al. (1977). Nucleotide sequence of bacteriophage phiX174 DNA. Nature 265, 687-695.",
    "Lander E. S., Waterman M. S. (1988). Genomic mapping by fingerprinting random clones: a mathematical analysis. Genomics 2, 231-239.",
    "Kececioglu J. D., Myers E. W. (1995). Combinatorial algorithms for DNA sequence assembly. Algorithmica 13, 7-51.",
    "Glover F. (1989, 1990). Tabu search, parts I and II. ORSA Journal on Computing 1 and 2.",
    "Ropke S., Pisinger D. (2006). An adaptive large neighborhood search heuristic for the pickup and delivery problem with time windows. Transportation Science 40, 455-472.",
    "Helsgaun K. (2017). An extension of the Lin-Kernighan-Helsgaun TSP solver for constrained traveling salesman and vehicle routing problems. Technical report, Roskilde University.",
    "Auer P., Cesa-Bianchi N., Fischer P. (2002). Finite-time analysis of the multiarmed bandit problem. Machine Learning 47, 235-256.",
    "Wilcoxon F. (1945). Individual comparisons by ranking methods. Biometrics Bulletin 1, 80-83.",
    "Holm S. (1979). A simple sequentially rejective multiple test procedure. Scandinavian Journal of Statistics 6, 65-70.",
    "Vargha A., Delaney H. D. (2000). A critique and improvement of the CL common language effect size statistics of McGraw and Wong. Journal of Educational and Behavioral Statistics 25, 101-132.",
]


def report_markdown(F: Facts, author: str = "[name to be filled in]") -> str:
    ucb, eps = F.verdict_ucb, F.verdict_eps
    rt = F.test("alns_ucb1_rt", "alns_random")
    rt_roulette = F.test("alns_ucb1_rt", "alns_roulette")
    cfg = F.config
    quality_rows = [
        [
            name(a),
            f"{F.mean_gap(a):.3f}",
            f"{F.median_gap(a):.3f}",
            f"{F.by_algorithm.loc[a, 'median_seconds']:.2f}",
        ]
        for a in sorted(
            (a for a in F.by_algorithm.index if not a.startswith("alns_fixed")), key=F.mean_gap
        )
    ]
    instance_rows = [
        [e["id"], "lambda" if e["genome"] == "NC_001416.1" else "phiX174", f"{e['start']}-{e['start'] + e['length']}",
         f"{e['coverage']:g}x", f"{100 * e['error']:.1f}%", str(e["read_length"])]
        for e in cfg["evaluation_instances"]
    ]  # fmt: skip
    tuned_rows = [
        [name(a), f"{v['default_mean_regret']:.3f}" if v["default_mean_regret"] is not None else "-", f"{v['best_mean_regret']:.3f}"]
        for a, v in F.tuned.items()
    ]  # fmt: skip
    bio_rows = [
        [name(a), f"{F.bio(a, 'genome_fraction'):.3f}", f"{F.bio(a, 'identity'):.3f}", f"{F.bio(a, 'n_contigs'):.0f}",
         f"{F.bio(a, 'n50'):.0f}", f"{F.bio(a, 'misjoins'):.0f}"]
        for a in ("greedy", "tabu", "tabu_cl", "alns_random", "alns_roulette", "alns_ucb1", "alns_ucb1_rt", "lkh3")
    ]  # fmt: skip
    main_pairs = [
        ("alns_random", "greedy"),
        ("alns_random", "tabu"),
        ("alns_random", "tabu_cl"),
        ("alns_random", "lkh3"),
        ("tabu_cl", "tabu"),
    ]
    sections: list[str] = []
    add = sections.append
    add(
        '<div class="cover"><h1>DNA fragment assembly as an asymmetric travelling salesman problem</h1>'
        "<p><strong>Does learning which repair tool to use help a large-neighbourhood search? "
        "An equal-budget comparison of Greedy, Tabu Search, LKH-3 and a bandit-driven ALNS</strong></p>"
        f"<p>Project report &middot; 2026-10-24 &middot; Author: {author}</p>"
        f"<p class='small'>All numbers in this report are generated from the committed results "
        f"(results/phase5, results/phase6) of {F.n_sets} read sets; the build script is "
        "<code>genomeasm docs build</code>.</p></div><div style='break-after:page'></div>"
    )
    add("## Abstract")
    add(
        f"Assembling a genome from short, overlapping, error-prone reads can be written as an asymmetric travelling "
        f"salesman problem (ATSP). We compare greedy edge merging, tabu search, the reference solver LKH-3 and an adaptive large "
        f"neighbourhood search (ALNS) whose destroy-and-repair tool is chosen by a multi-armed bandit, on {F.n_sets} read sets "
        f"simulated from two real viral genomes, with the same {F.budget:g} s wall-clock budget, hyper-parameters tuned on separate "
        f"instances, and paired non-parametric statistics. **The pre-registered claim, that bandit selection (UCB1, epsilon-greedy) "
        f"finds better assemblies than random selection, is not supported** (UCB1 vs random: {p_text(ucb['p_holm'])}, A12 = {ucb['a12']:.2f}; "
        f"epsilon-greedy is slightly worse). A tie-breaking defect in our UCB1 was found after the results; a corrected variant is "
        f"slightly better than random (A12 = {rt.a12:.2f}) and equal to the classical roulette wheel, an exploratory result. "
        f"ALNS clearly beats greedy and tabu search; LKH-3 was never beaten. All solvers rebuild the same genome fraction "
        f"({F.bio('lkh3', 'genome_fraction'):.3f}) without a single misjoin."
    )
    add("## 1. Introduction")
    add(
        "A sequencing machine cannot read a genome end to end. It returns many short *reads*, taken at random places, with a few "
        "wrong letters. *Assembly* puts the reads back in order. Two reads probably were neighbours if the end of one equals the "
        "beginning of the other; the longer this *overlap*, the safer the guess. Choosing an order that makes the overlaps as long as "
        "possible is a classical combinatorial problem: it is exactly a travelling salesman problem, and with the order of the "
        "overlap it is the *asymmetric* version (Kececioglu and Myers 1995)."
    )
    add(
        "Good heuristics exist. We study one of them, ALNS, which in every iteration destroys part of the current order and repairs it "
        "with one of six destroy-repair tools. Which tool works best changes during the search, and choosing among them is a "
        "*multi-armed bandit* problem. **Research question: under an equal computational budget, does bandit-based operator selection "
        "find better assemblies than random operator selection?** The claim was fixed in advance; we report the outcome whatever it is."
    )
    add(
        "**Contributions.** (i) A verified encoding of read ordering as an ATSP with proofs; (ii) a reproducible, offline pipeline from "
        "real NCBI genomes to simulated reads, overlap graphs, five solvers and ground-truth quality metrics; (iii) a fair, paired, "
        f"{F.n_sets}-read-set comparison with honest reporting of a negative result and two post-hoc variants; (iv) an interactive "
        "application for exploring all of it."
    )
    add("## 2. Background")
    add(
        "**DNA and RNA.** DNA is a long word over four letters (bases) A, C, G, T; in cells it forms a double helix in which A pairs "
        "with T and G with C. RNA is usually a single strand with U in place of T. A *genome* is the complete DNA (or RNA) of an organism; "
        "its length is the number of letters. We work with the order of the letters only, so RNA would be the same problem with U for T "
        "(not tested here)."
    )
    add(
        "**Reads and coverage.** A sequencer returns reads of 50 to 300 letters. *Coverage* c is the average number of reads covering a "
        "letter. If read starts are random, the share of letters covered by no read is about $e^{-c}$ (Lander and Waterman 1988); on a short "
        "linear window the ends matter, so we also use the exact finite-genome model, which matches the simulator better (Figure 1)."
    )
    add(
        fig(
            "fig01a_coverage_histogram",
            "Figure 1. Coverage depth of simulated reads against the Lander-Waterman (Poisson) prediction and the exact model (Phase 2 validation).",
        )
    )
    add(
        "**Overlap graph.** The overlap from read $i$ to read $j$ is a number; with $n$ reads these form an $n \\times n$ table "
        "that is not symmetric (Figure 2). With read errors the table must allow a few mismatches."
    )
    add(
        fig(
            "fig03_overlap_heatmap",
            "Figure 2. Overlap matrix of a simulated read set, reads sorted by true position: the diagonal band holds the overlaps of true neighbours.",
            "62%",
        )
    )
    add(
        "**Travelling salesman and large-neighbourhood search.** The TSP asks for the cheapest tour through all cities. ALNS "
        "(Ropke and Pisinger 2006) alternates *destroy* (remove some elements of the solution) and *repair* (re-insert them well), accepting "
        "a worse result sometimes (simulated annealing). **Bandits** (Auer et al. 2002) choose among several actions whose value is "
        "unknown and changes: *exploit* the best action so far or *explore* another. LKH-3 (Helsgaun 2017) is a state-of-the-art TSP solver "
        "based on Lin-Kernighan moves; tabu search (Glover 1989, 1990) escapes local optima by forbidding recent moves."
    )
    add("## 3. Problem formalisation")
    add(
        "Let $r_1,\\dots,r_n$ be substring-free reads (no duplicates, no read inside another) with lengths $\\ell_1,\\dots,\\ell_n$ and let $\\tau$ be the minimum overlap."
    )
    add(
        eq(
            r"\mathrm{ov}(i,j)=\max\{\,k:\ \tau\le k<\min(\ell_i,\ell_j),\ \mathrm{suffix}_k(r_i)=\mathrm{prefix}_k(r_j)\,\},\quad 0\ \mathrm{if\ none}"
        )
    )
    add(
        "With read errors the equality is replaced by *Hamming distance at most* $\\lfloor \\epsilon k\\rfloor$ ($\\epsilon$ = 0.05 here). "
        "A *solution* is a permutation $\\pi$ of the reads (the search space has $n!$ elements). Fitness and merged length are"
    )
    add(
        eq(
            r"F(\pi)=\sum_{k=1}^{n-1}\mathrm{ov}(\pi_k,\pi_{k+1})\ \ (\mathrm{maximise}),\qquad L(\pi)=\sum_{i=1}^{n}\ell_i-F(\pi)"
        )
    )
    add(
        "**Proposition 1** (length). Gluing the reads in the order $\\pi$ gives a text of length $L(\\pi)$ (counting pieces separately where "
        "an overlap is 0). *Proof.* By induction on $k$: the first $k$ reads give a text of length $\\sum_{m\\le k}\\ell_{\\pi_m}-\\sum_{m<k}\\mathrm{ov}(\\pi_m,\\pi_{m+1})$. "
        "Gluing read $\\pi_{k+1}$ appends its $\\ell_{\\pi_{k+1}}$ letters minus the $\\mathrm{ov}(\\pi_k,\\pi_{k+1})$ letters it shares with the end. "
        "Where an overlap is 0 the pieces are separate and their lengths still add up. □"
    )
    add("**ATSP encoding.** Add a dummy node 0 and define costs")
    add(eq(r"c(i,j)=\ell_j-\mathrm{ov}(i,j),\qquad c(0,j)=\ell_j,\qquad c(i,0)=0"))
    add(
        "**Proposition 2** (equivalence). Every Hamiltonian cycle through node 0 corresponds to exactly one order $\\pi$ (cut the cycle at 0), and its "
        "cost is $L(\\pi)$. Hence the cheapest tour gives the order of maximal fitness, and $\\min_\\pi L=\\sum\\ell_i-\\max_\\pi F$. *Proof.* The cycle "
        "$0\\to\\pi_1\\to\\dots\\to\\pi_n\\to0$ costs $c(0,\\pi_1)+\\sum_{k<n}c(\\pi_k,\\pi_{k+1})+c(\\pi_n,0)=\\ell_{\\pi_1}+\\sum_{k<n}(\\ell_{\\pi_{k+1}}-\\mathrm{ov}(\\pi_k,\\pi_{k+1}))+0=\\sum_i\\ell_i-F(\\pi)=L(\\pi)$. "
        "The costs are non-negative because $\\mathrm{ov}(i,j)<\\ell_j$. Since $\\mathrm{ov}(i,j)\\ne\\mathrm{ov}(j,i)$ in general, the matrix is not symmetric: an ATSP. □"
    )
    add(
        "**Delta evaluation.** Swapping two reads or moving one read changes at most six overlaps, so the change of $F$ is computed in constant "
        "time; all local-search moves use this. **What $F$ does not know:** with read errors the true order can have a *lower* $F$ than "
        "the best order found (it is not the optimum of $F$), which is why quality is also measured against the real genome (Section 8)."
    )
    add("## 4. Methods")
    add(
        "**Data.** Two real RefSeq genomes, phiX174 (NC_001422.1, 5,386 letters, circular single-stranded DNA) and lambda (NC_001416.1, "
        "48,502 letters, linear double-stranded DNA), downloaded once, checked against NCBI's metadata, fingerprinted with SHA-256 and "
        "stored offline; each genome is used as one linear string (D-016, D-017, D-018)."
    )
    add(
        "**Simulator.** Reads of fixed length start at uniformly random positions of a window; each letter is substituted with the error "
        "probability; a seed fixes everything (D-023). The simulator was validated against theory (Phase 2)."
    )
    add(
        "**Overlaps.** Substring-free preprocessing, then an exhaustive error-tolerant overlap search (compiled with numba), minimum overlap "
        f"$\\tau$ = {cfg['overlap']['min_overlap']} (15 for 60-letter reads), tolerance {cfg['overlap']['tolerance']:g} (D-026 to D-028, D-036)."
    )
    add(
        "**Solvers.** *Greedy*: accept overlaps from the largest down if no read gets two successors or predecessors and no cycle closes. "
        "*Tabu search*: swap and insert moves, tenure, aspiration, restarts; *tabu_cl* (post hoc, D-045) samples insert moves from candidate lists. "
        "*LKH-3*: the Helsgaun solver called as an external program, every returned tour verified, never replaced silently (D-033, D-041). "
        "*ALNS*: three destroy operators (random, worst-connected, consecutive segment) times two repair operators (greedy insertion, regret-2) = six arms; "
        "simulated-annealing acceptance with a calibrated start temperature; rewards for a new best, an improvement, an accepted worse move (D-034)."
    )
    add(
        "**Operator-selection policies.** *Random* (the baseline of the claim), *round robin*, *roulette wheel* (classical adaptive weights), "
        "*epsilon-greedy* and *UCB1* with a sliding window over the last rewards of each arm: "
    )
    add(eq(r"\mathrm{UCB1}(a)=\bar{x}_a+c\sqrt{\frac{2\ln N}{n_a}}"))
    add(
        "(D-038). Rewards are in $[0,1]$. Our UCB1 broke ties by the lowest arm index; when all rewards are equal it plays arm 0 for ever (D-046). "
        "The variant *UCB1 with random ties* was added after the results and is reported separately."
    )
    add(
        f"**Budgets and tuning.** Every timed solver gets {F.budget:g} s of wall-clock time per read set; LKH-3's search time is the budget minus its "
        "measured overhead (D-041). Hyper-parameters of every algorithm were chosen by random search (24 trials, trial 0 = defaults) on six "
        "tuning instances of phiX174 positions 0 to 2,500 with seeds from 1000, never used for evaluation (D-040)."
    )
    add("## 5. Experimental design")
    add(
        f"There are {F.n_settings} instance settings times {F.seeds} seeds = {F.n_sets} read sets. A name such as L5-c10-e1.0 means a lambda window of 5,000 "
        "letters, coverage 10x, error rate 1.0%."
    )
    add(table(["Setting", "Genome", "Window", "Coverage", "Errors", "Read length"], instance_rows))
    add(
        "**Fairness.** Same read sets and the same wall-clock budget for every solver; order of the timed solvers shuffled per read set; tuning data "
        "disjoint from evaluation data; seeds independent of the order of the runs (SeedSequence, D-039). **Comparison.** Fitness is turned into the "
        "gap to the best known fitness of the read set:"
    )
    add(eq(r"\mathrm{gap}=100\,\frac{F^{*}-F}{F^{*}}"))
    add(
        "where $F^{*}$ is the best fitness any solver reached (not a proven optimum). Solvers are compared pairwise on the same read sets with the "
        "two-sided Wilcoxon signed-rank test (Wilcoxon 1945), Holm correction within each family of tests (Holm 1979), a bootstrap interval "
        "(10,000 resamples) of the paired difference, and the Vargha-Delaney effect size "
    )
    add(eq(r"A_{12}=P(X>Y)+\frac{1}{2}P(X=Y)"))
    add(
        "(Vargha and Delaney 2000; 0.5 = no difference). The primary family has two tests, fixed before the runs: UCB1 vs random and epsilon-greedy vs random (D-042)."
    )
    add("## 6. Results")
    add(
        f"**Final quality.** After {F.budget:g} s the ranking by mean gap is LKH-3, the ALNS variants, tabu search and greedy (Table 1, Figure 3). "
        f"LKH-3 was {'never' if F.lkh_never_beaten else 'sometimes'} beaten on any of the {F.n_sets} read sets. Greedy already finds the best known solution in "
        f"{100 * F.greedy_equals_best:.0f}% of the read sets (the easy ones: low errors); the gaps appear on noisy, deep read sets."
    )
    add(table(["Solver", "Mean gap (%)", "Median gap (%)", "Median seconds"], quality_rows))
    add("*Table 1. Gap to the best known fitness over all read sets.*")
    add(
        fig(
            "fig05_final_fitness",
            "Figure 3. Gap to the best known fitness: (a) all read sets, (b) median per instance setting.",
        )
    )
    add(
        fig(
            "fig04_convergence",
            "Figure 4. Mean gap of the best solution so far against time and against evaluations (log axes). LKH-3 appears only after its preprocessing.",
        )
    )
    add(
        "**Against the baselines.** ALNS with random choice beats greedy and tabu search with large effect, and the post-hoc tabu variant improves on plain tabu "
        "but not on ALNS (Table 2, with the test of every pair against the first solver)."
    )
    add(table(TEST_HEADER, [test_row(F, a, b) for a, b in main_pairs]))
    add(
        "*Table 2. Paired comparisons over all read sets (positive differences favour the first solver).*"
    )
    add(
        fig(
            "fig06_gap_to_lkh3",
            "Figure 5. Gap to LKH-3 against coverage (a) and error rate (b), lambda 5,000-letter windows.",
        )
    )
    add(
        fig(
            "fig10_runtime_vs_n",
            "Figure 6. Runtime against the number of reads. The overlap-graph build and LKH-3's preprocessing grow fastest.",
            "70%",
        )
    )
    add("## 7. Ablation: does the learning rule matter?")
    add(table(TEST_HEADER, [test_row(F, "alns_ucb1", "alns_random"), test_row(F, "alns_epsilon_greedy", "alns_random"),
                           test_row(F, "alns_roundrobin", "alns_random"), test_row(F, "alns_roulette", "alns_random"),
                           test_row(F, "alns_ucb1", "alns_roulette")]))  # fmt: skip
    add(
        "*Table 3. Operator-selection policies against random selection (first two rows: the pre-registered family).*"
    )
    add(
        f"**The claim is not supported.** UCB1 is not different from random ({p_text(ucb['p_holm'])}, A12 = {ucb['a12']:.3f}); epsilon-greedy is significantly *worse* "
        f"({p_text(eps['p_holm'])}, A12 = {eps['a12']:.3f}); round robin and the roulette wheel are slightly better than random (A12 about 0.54). All effects are small; "
        "with 300 paired read sets the tests are sensitive, so significance is not size. A single operator pair used alone is worse than the random mixture in all six cases."
    )
    add(
        fig(
            "fig07_ablation",
            "Figure 7. Mean paired difference to ALNS random with 95% bootstrap intervals (left) and per instance setting (right).",
        )
    )
    add(
        f"**A defect found after the results (D-046).** Figure 8 shows UCB1 collapsing onto the first operator pair in the second half of a run: "
        f"pooled, that pair received {100 * F.arm_zero_share:.1f}% of its uses although it has the lowest success rate under uniform use. The cause is the "
        "tie-breaking: once the window of every arm holds zeros, all indices are equal and the lowest index wins. Reproduced without ALNS, "
        "the policy plays arm 0 in 100% of the plays. The original UCB1 and its results were kept; a variant with random ties was tuned and run afterwards. "
        f"It is better than random ({p_text(rt.p_holm)}, A12 = {rt.a12:.3f}) and not different from the roulette wheel ({p_text(rt_roulette.p_holm)}). "
        "It is an *exploratory* result: it was added after seeing the data and does not rescue the pre-registered claim."
    )
    add(
        fig(
            "fig08_arm_probability",
            "Figure 8. Share of the choices of each operator pair over a run (mean of 300 runs). UCB1 collapses onto the first pair; with random ties it does not.",
        )
    )
    add(
        fig(
            "fig09_operator_heatmap",
            "Figure 9. What each destroy-repair pair achieved (share of its uses ending in each outcome).",
            "85%",
        )
    )
    add("## 8. Biological quality")
    add(
        "Fitness is not the genome, so every assembly is also graded against the answer key: windows of 100 letters of every contig are located on the reference "
        "by shared 15-letter words (D-043). *Identity* is the share of matching letters, *genome fraction* the share of the window rebuilt, *N50* the contig length "
        "such that longer contigs hold half of the assembled letters, a *misjoin* a contig that jumps to a distant place."
    )
    add(table(["Solver", "Genome fraction", "Identity", "Contigs", "N50", "Misjoins"], bio_rows))
    add("*Table 4. Medians over all read sets.*")
    add(
        "All solvers rebuild practically the same genome fraction and none produces a misjoin; the stronger solvers have slightly higher identity and fewer "
        "contigs. The true read order itself is fragmented, because reads with errors do not overlap enough (D-031), and the best order found has a higher fitness than "
        "the true order in almost every read set (D-036)."
    )
    add("## 9. Discussion")
    add(
        "**The main result is negative.** Learning which repair tool to use did not help beyond random mixing. Hypotheses, not tested here (`TODO-VERIFY`): six similar "
        "tools leave little to learn; improving moves are rare (a few per cent), so rewards are sparse; the tuning set may not represent the evaluation set (two post-hoc "
        "variants improved on the evaluation set more than the tuning data suggested). A definitive answer needs more diverse operators, longer runs and rewards that are denser."
    )
    add(
        "**What does hold.** The ATSP encoding is exact; ALNS with a random mixture of tools is a strong heuristic; the solver choice matters for the objective but hardly for the "
        "recovered genome; and an honest pipeline (separate tuning data, equal budgets, paired tests, disclosure of post-hoc work) is the main safeguard against over-claiming. "
        "A bug can hide in a baseline: the defect in UCB1 was visible only because the policy's behaviour was plotted."
    )
    add("## 10. Threats to validity")
    for item in (
        f"**Internal.** 'Best known' is the best solution found by our solvers, not a proven optimum; LKH-3 is a heuristic. LKH-3's wall-clock exceeded the budget in {F.lkh_over:.0%} of the runs, by up to {F.lkh_max - F.budget:.1f} s, a slight advantage. Five workers shared four physical cores, so timings carry noise.",
        "**Post-hoc work.** `tabu_cl` and the UCB1 variant were added after the first results; their tests form separate families. Choosing which variants to add is a forking path.",
        "**Tuning.** 24 random trials on 18 tuning read sets give optimistic regrets (winner's curse) and did not always transfer.",
        "**External.** Two viral genomes, substitution errors only, one strand, windows up to 10,000 letters; real data have insertions, deletions, both strands and repeats.",
        "**Construct.** Fitness is not the genome (Section 8); the biological metrics rely on locating windows by 15-letter words.",
        "**Software.** The solver code was written for this project; correctness is covered by the test suite, exact validators on small instances and independent verification of every LKH-3 tour.",
    ):
        add("- " + item)
    add("## 11. Reproducibility statement")
    add(
        f"The repository is self-contained and offline. `make venv` installs the pinned packages; `make check` runs lint, types and all tests; `genomeasm experiments "
        f"tune|run|scale|analyze` regenerates the results (resumable); `genomeasm docs build` regenerates this report, the slides and the oral-exam sheet from "
        f"`results/`; `genomeasm app run` starts the application. Every result file stores the configuration hash, the commit, the seed and the library versions; "
        f"`results/phase5/stages.json` records the commit of each stage. Design decisions are in `docs/DECISIONS.md`. The worst interactive run of the application takes "
        f"{F.worst_live_seconds:.0f} s of the 30 s rule (`results/phase6`). LKH-3 (research-use licence) is not redistributed."
    )
    add("## References")
    add(
        '<p class="small">All references are given from memory and are marked <code>TODO-VERIFY</code> until they are checked against the source.</p>'
    )
    add(
        "<ol class='small'>"
        + "".join(f"<li>{r} <code>TODO-VERIFY</code></li>" for r in REFERENCES)
        + "</ol>"
    )
    add("## Appendix A. Tuned hyper-parameters")
    add(table(["Algorithm", "Default regret (%)", "Tuned regret (%)"], tuned_rows))
    add(
        "*Mean regret to the best fitness any trial reached on the tuning instances. The tuned values themselves are in `results/phase5/tuning_summary.json` and in `docs/ORAL_EXAM.md`.*"
    )
    return (chr(10) + chr(10)).join(sections)


def build_report(out_dir: Path, F: Facts, author: str = "[name to be filled in]") -> Path:
    """Write report.html and report.pdf into ``out_dir``; returns the PDF path."""
    out_dir.mkdir(parents=True, exist_ok=True)
    body = md.to_html(report_markdown(F, author))
    page = f"<!doctype html><html><head><meta charset='utf-8'><title>DNA fragment assembly as TSP</title><style>{CSS}</style></head><body>{body}</body></html>"
    html_file = out_dir / "report.html"
    html_file.write_text(page, encoding="utf-8", newline=chr(10))
    pdf = out_dir / "report.pdf"
    md.print_pdf(html_file, pdf)
    return pdf


__all__ = ["build_report", "report_markdown", "ROOT"]
